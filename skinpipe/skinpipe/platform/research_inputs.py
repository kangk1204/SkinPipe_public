"""Typed research inputs without implicit compound-to-protein substitution.

Identity is a graph or a provider-declared protein identity. It is separate
from preparation readiness, measured activity, and skin benefit. SMILES stay
local; only explicit name lookups contact public identity providers.
"""

from __future__ import annotations

from copy import deepcopy
import json
import re
import unicodedata
from urllib.parse import quote

from rdkit import Chem

from .artifacts import canonical_json
from .identity import chemical_identity, digest
from .targets import ACCESSION, SourceClient, UNIPROT, now, target_request

INPUT_SCHEMA = "research-input-v1"
DESCRIPTOR_SCHEMA = "research-input-descriptor-v1"
PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
MAX_NAME_MATCHES = 20
MAX_RESPONSE_BYTES = 8 * 1024 * 1024
MAX_INPUT_ATOMS = 512


def normalize_input(body):
    """Validate a small typed IR while retaining the exact submitted fields.

    A normalized value may cross the queue boundary, but its derived fields
    must still agree with its original input. Client verification assertions
    and arbitrary descriptor fields are never accepted here.
    """
    if not isinstance(body, dict):
        raise ValueError("input must be an object")
    if body.get("schema") == INPUT_SCHEMA:
        original = body.get("original_input")
        if not isinstance(original, dict) or original.get("schema") is not None:
            raise ValueError("Invalid original input")
        expected = normalize_input(original)
        if body != expected:
            raise ValueError("Normalized input differs from its original input")
        return expected
    allowed = {"kind", "format", "value", "taxon_id", "desired_action", "isoform", "pubchem_cid",
               "selection", "ligand_smiles", "ligand_sdf"}
    if set(body) - allowed:
        raise ValueError("Unsupported research input fields")
    kind, form, value = body.get("kind"), body.get("format", "name"), body.get("value")
    if not isinstance(kind, str) or kind not in {"protein", "compound", "bound_complex"}:
        raise ValueError("input.kind must be protein, compound or bound_complex")
    if kind == "bound_complex":
        from .research_structures import normalize_structure_input
        return normalize_structure_input(body)
    if set(body) & {"selection", "ligand_smiles", "ligand_sdf"}:
        raise ValueError("Structure selections only apply to bound_complex input")
    if not isinstance(form, str) or form not in {"name", "smiles"} or (kind == "protein" and form != "name"):
        raise ValueError("Proteins accept names; compounds accept names or SMILES")
    limit = 16384 if form == "smiles" else 256
    if (not isinstance(value, str) or not value.strip() or len(value) > limit
            or any(unicodedata.category(c).startswith("C") for c in value)):
        raise ValueError("Enter one bounded name or SMILES without control characters")
    result = {"schema": INPUT_SCHEMA, "kind": kind, "format": form,
              "value": value, "original_input": deepcopy(body)}
    if kind == "protein":
        if "pubchem_cid" in body:
            raise ValueError("pubchem_cid is only valid for compound names")
        request = target_request({"kind": "protein", "query": value.strip(),
            **{key: body[key] for key in ("taxon_id", "desired_action", "isoform") if key in body}})
        result.update(taxon_id=request["taxon_id"], desired_action=request["desired_action"],
                      species_defaulted=request["species_defaulted"],
                      desired_action_defaulted="desired_action" not in body)
        if "isoform" in body:
            result["isoform"] = body["isoform"]
    else:
        if "taxon_id" in body or "isoform" in body:
            raise ValueError("Compound identity has no protein species or isoform")
        if body.get("desired_action", "any") != "any":
            raise ValueError("A compound input alone does not define a protein action")
        cid = body.get("pubchem_cid")
        if "pubchem_cid" in body:
            if form != "name" or isinstance(cid, bool) or not isinstance(cid, int) or not 1 <= cid <= 2147483647:
                raise ValueError("pubchem_cid must be a positive integer for a compound name")
            result["pubchem_cid"] = cid
        result["desired_action"] = "any"
    return result


def _protein_request(spec, *, query=None):
    request = {"kind": "protein", "query": spec["value"].strip() if query is None else query,
               "desired_action": spec["desired_action"]}
    if not spec["species_defaulted"]:
        request["taxon_id"] = spec["taxon_id"]
    if "isoform" in spec:
        request["isoform"] = spec["isoform"]
    return request


class CompoundSourceClient:
    """Small official PubChem GET client with bounded immutable receipts."""

    def __init__(self, artifacts, *, get=None):
        self.artifacts, self.get, self.receipts = artifacts, get, []

    def request(self, url, params=None):
        import requests

        if not url.startswith(PUBCHEM + "/compound/"):
            raise ValueError("Only the official PubChem compound provider is allowed")
        params, attempted = params or {}, now()
        try:
            if self.get:
                data = self.get(url, params)
            else:
                with requests.get(url, params=params, timeout=(5, 20), allow_redirects=False,
                                  stream=True) as response:
                    response.raise_for_status()
                    if response.status_code != 200:
                        raise ValueError("Provider response must be a complete HTTP 200")
                    raw = bytearray()
                    for block in response.iter_content(65536):
                        raw.extend(block)
                        if len(raw) > MAX_RESPONSE_BYTES:
                            raise ValueError("PubChem response exceeds the identity lookup limit")
                    data = json.loads(raw, parse_constant=lambda _: (_ for _ in ()).throw(
                        ValueError("Nonfinite source JSON")))
            raw = canonical_json(data)
            if len(raw) > MAX_RESPONSE_BYTES or not isinstance(data, dict):
                raise ValueError("Invalid or oversized PubChem response")
            receipt = self.artifacts.put_bytes(raw, media_type="application/json")
            self.receipts.append({"source": "PubChem PUG REST", "url": url, "params": params,
                "retrieved_at": attempted, "raw_sha256": receipt["sha256"], "status": "complete"})
            return data
        except (requests.RequestException, ValueError, OSError) as error:
            not_found = isinstance(error, requests.HTTPError) and error.response is not None and error.response.status_code == 404
            self.receipts.append({"source": "PubChem PUG REST", "url": url, "params": params,
                "retrieved_at": attempted, "raw_sha256": None,
                "status": "not_found" if not_found else "source_unavailable"})
            if not_found:
                return {"IdentifierList": {"CID": []}}
            raise


def _base(spec):
    return {"schema": DESCRIPTOR_SCHEMA, "kind": spec["kind"], "input": deepcopy(spec),
            "status": "selection_required", "reason": "identity_pending", "target_id": None,
            "canonical_smiles": None, "desired_action": spec["desired_action"],
            "name": spec["value"].strip(), "choices": [], "source_receipts": [],
            "limitations": [], "proxy_target": None, "actual_protein_id": None,
            "compute_readiness": "identity_pending"}


def _class_name(name):
    normalized = unicodedata.normalize("NFKC", name).strip().casefold()
    normalized = re.sub(r"[\s_-]+", " ", normalized)
    return (normalized in {"ceramide", "ceramides", "세라마이드", "hyaluronan", "hyaluronic acid",
                           "sodium hyaluronate", "히알루론산", "collagen", "콜라겐", "elastin",
                           "sphingolipid", "sphingolipids", "lipid", "lipids", "essential oil"}
            or bool(re.fullmatch(r"ceramide\s+(?:np|ns|ap|as|eop|eos|eoh|nh|ah|ads)", normalized)))


def _structure_descriptor(spec, smiles, *, name=None, cid=None):
    result = _base(spec)
    parser = Chem.SmilesParserParams()
    parser.parseName = False
    mol = Chem.MolFromSmiles(smiles, parser)
    if mol is None or not mol.GetNumAtoms():
        return {**result, "reason": "invalid_smiles", "compute_readiness": "invalid_structure"}
    if mol.GetNumAtoms() > MAX_INPUT_ATOMS:
        return {**result, "reason": "structure_size_limit", "compute_readiness": "unsupported_structure_size",
                "limitations": [f"Identity input limit is {MAX_INPUT_ATOMS} atoms; provide one defined small molecule."]}
    if any(atom.GetAtomicNum() == 0 for atom in mol.GetAtoms()):
        return {**result, "reason": "concrete_structure_required", "compute_readiness": "generic_or_polymer_structure",
                "limitations": ["Wildcard atoms and unspecified repeat units do not identify a concrete molecule."]}
    # Exact identity retains the complete salt/charge/isotope/stereochemistry.
    # The selected largest fragment is never silently promoted to the input.
    identity = chemical_identity(Chem.MolToCXSmiles(mol))
    identity["original_smiles"] = smiles
    identity["original_structure_id"] = "raw:" + digest(smiles)
    result.update(identity=identity, canonical_smiles=identity["canonical_smiles"], exact_smiles=identity["identity_smiles"],
                  target_id="compound:" + identity["exact_structure_id"].split(":", 1)[1],
                  name=name or spec["value"].strip(), pubchem_cid=cid,
                  provenance={"identity_basis": "user_provided_structure" if spec["format"] == "smiles" else "provider_name_to_structure",
                              "identity_version": identity["identity_version"],
                              "activity_basis": "none", "protein_substitution": False})
    if identity["parent_ambiguous"]:
        return {**result, "reason": "multiple_organic_components_require_selection",
                "compute_readiness": "component_selection_required",
                "choices": [{"smiles": s} for s in identity["component_identities"]],
                "limitations": ["Multiple organic components are retained. Select a defined seed material explicitly."]}
    limitations = ["Structural similarity and skin-property predictions do not establish activity or skin benefit."]
    if identity["unassigned_stereo"] or identity["unassigned_bond_stereo"]:
        limitations.append("Unspecified stereochemistry is retained; no stereoisomer was assigned.")
    if identity["fragment_count"] > 1:
        limitations.append("Salt components are retained in exact identity; ligand preparation must declare its parent policy.")
    conventional = {1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 35, 53}
    unsupported_elements = sorted({atom.GetSymbol() for atom in mol.GetAtoms() if atom.GetAtomicNum() not in conventional})
    readiness = "ligand_only"
    if unsupported_elements:
        limitations.append("Elements outside the conventional organic prediction domain: " + ", ".join(unsupported_elements))
        readiness = "ligand_only_domain_review"
    if identity["heavy_atoms"] > 256:
        limitations.append("This graph exceeds the 256-heavy-atom conformer-comparison limit.")
        readiness = "ligand_only_domain_review"
    return {**result, "status": "verified", "reason": "exact_structure_verified",
            "compute_readiness": readiness, "limitations": limitations}


def _compound(spec, artifacts, source_get):
    if spec["format"] == "smiles":
        return _structure_descriptor(spec, spec["value"].strip())
    result = _base(spec)
    if _class_name(spec["value"]):
        return {**result, "reason": "compound_class_requires_exact_structure",
                "limitations": ["This name denotes a chemical class, polymer or variable composition. Provide an exact member name or SMILES; no protein was substituted."]}
    client = CompoundSourceClient(artifacts, get=source_get)
    try:
        response = client.request(PUBCHEM + "/compound/name/" + quote(spec["value"].strip(), safe="") + "/cids/JSON",
                                  {"name_type": "complete"})
        header = response.get("IdentifierList")
        ids = header.get("CID") if isinstance(header, dict) else None
        if not isinstance(ids, list) or any(isinstance(cid, bool) or not isinstance(cid, int) or cid <= 0 for cid in ids):
            raise ValueError("PubChem did not return a bounded exact identity list")
        ids = sorted(set(ids))
        if not ids:
            return {**result, "reason": "identity_not_found", "source_receipts": client.receipts}
        selected = spec.get("pubchem_cid")
        if len(ids) > MAX_NAME_MATCHES and selected is None:
            return {**result, "reason": "name_match_limit_requires_exact_structure",
                    "choices": [{"pubchem_cid": cid} for cid in ids[:MAX_NAME_MATCHES]],
                    "name_match_count": len(ids), "source_receipts": client.receipts}
        if selected is not None and selected not in ids:
            return {**result, "reason": "selected_cid_not_in_name_results",
                    "choices": [{"pubchem_cid": cid} for cid in ids[:MAX_NAME_MATCHES]], "source_receipts": client.receipts}
        # Check membership against the complete exact-name result, then fetch
        # only the user's chosen CID. Do not repeatedly offer an unusable list.
        if selected is not None and len(ids) > MAX_NAME_MATCHES:
            ids = [selected]
        response = client.request(PUBCHEM + "/compound/cid/" + ",".join(map(str, ids)) + "/property/Title,SMILES,InChIKey/JSON")
        table = response.get("PropertyTable")
        records = table.get("Properties") if isinstance(table, dict) else None
        if not isinstance(records, list):
            raise ValueError("PubChem property table missing")
        by_id = {}
        for record in records:
            cid = record.get("CID") if isinstance(record, dict) else None
            if isinstance(cid, bool) or cid not in ids or cid in by_id:
                raise ValueError("PubChem property identity conflict")
            by_id[cid] = record
        if set(by_id) != set(ids):
            raise ValueError("PubChem property page incomplete")
        choices = [{"pubchem_cid": cid, "name": by_id[cid].get("Title"),
                    "smiles": by_id[cid].get("SMILES") or by_id[cid].get("IsomericSMILES"),
                    "inchi_key": by_id[cid].get("InChIKey")} for cid in ids]
        if len(ids) != 1 and selected is None:
            return {**result, "reason": "ambiguous_compound_name", "choices": choices,
                    "source_receipts": client.receipts}
        cid = selected if selected is not None else ids[0]
        record = by_id[cid]
        smiles = record.get("SMILES") or record.get("IsomericSMILES")
        if not isinstance(smiles, str) or not smiles or len(smiles) > 16384:
            raise ValueError("PubChem exact stereochemical SMILES missing")
        descriptor = _structure_descriptor(spec, smiles, name=record.get("Title"), cid=cid)
        descriptor.update(source_receipts=client.receipts, name_choices=choices)
        if "provenance" in descriptor:
            descriptor["provenance"]["provider_inchi_key"] = record.get("InChIKey")
        return descriptor
    except (ValueError, OSError) as error:
        return {**result, "status": "source_unavailable", "reason": "compound_identity_source_unavailable",
                "source_receipts": client.receipts, "error_type": type(error).__name__}
    except Exception as error:
        # Requests is imported only by the provider client. Do not turn arbitrary
        # implementation faults into a scientific success or a chosen identity.
        import requests
        if not isinstance(error, requests.RequestException):
            raise
        return {**result, "status": "source_unavailable", "reason": "compound_identity_source_unavailable",
                "source_receipts": client.receipts, "error_type": type(error).__name__}


def _protein_names(record):
    description = record.get("proteinDescription", {})
    names = []
    if isinstance(description, dict):
        recommended = description.get("recommendedName", {})
        names.append(recommended)
        names.extend(description.get("alternativeNames", []))
        names.extend(description.get("submissionNames", []))
    result = set()
    for entry in names:
        if not isinstance(entry, dict):
            continue
        values = [entry.get("fullName", {})] + entry.get("shortNames", [])
        for value in values:
            if isinstance(value, dict) and isinstance(value.get("value"), str):
                result.add(value["value"].strip().casefold())
    return result


def _protein(spec, artifacts, resolver, source_get):
    request = _protein_request(spec)
    query = request["query"]
    if query.casefold() in {"collagen", "collagens", "keratin", "keratins", "integrin", "integrins", "interleukin", "ceramide", "ceramides", "hyaluronan"}:
        return {**_base(spec), "reason": "protein_family_requires_exact_member",
                "limitations": ["Select an exact protein family member; a compound or family label does not identify one protein."]}
    # Existing exact accession/gene/ChEMBL semantics remain authoritative. A
    # genuine full protein name gets an exact provider-declared name check.
    descriptor = None
    if ACCESSION.fullmatch(query.upper()) or re.fullmatch(r"(?:CHEMBL[0-9]+|[A-Za-z][A-Za-z0-9-]{0,29})", query):
        descriptor = resolver(request)
        if descriptor.get("status") == "verified" or descriptor.get("reason") not in {"identity_not_found"}:
            return {**descriptor, "input": spec, "source_receipts": descriptor.get("sources", []),
                    "proxy_target": None, "original_request": request}
    client = SourceClient(artifacts, get=source_get)
    try:
        term = query.replace("\\", "\\\\").replace('"', '\\"')
        response = client.request(UNIPROT + "/search", {"query": f'protein_name:"{term}" AND organism_id:{spec["taxon_id"]} AND reviewed:true',
            "format": "json", "size": 100})
        records = response.get("results", [])
        if not isinstance(records, list) or len(records) >= 100:
            return {**_base(spec), "reason": "protein_name_page_incomplete", "source_receipts": client.receipts}
        records = [r for r in records if isinstance(r, dict) and query.casefold() in _protein_names(r)
                   and r.get("organism", {}).get("taxonId") == spec["taxon_id"]]
        choices = [{"accession": r.get("primaryAccession"), "taxon_id": spec["taxon_id"],
                    "name": query, "source_record_sha256": artifacts.put_bytes(canonical_json(r))["sha256"]} for r in records]
        if len(choices) != 1:
            return {**_base(spec), "reason": "ambiguous_protein_name" if choices else "identity_not_found",
                    "choices": choices, "source_receipts": client.receipts}
        accession = choices[0]["accession"]
        if not isinstance(accession, str) or not ACCESSION.fullmatch(accession):
            raise ValueError("Provider protein accession invalid")
        descriptor = resolver(_protein_request(spec, query=accession))
        return {**descriptor, "input": spec, "original_request": request,
                "source_receipts": client.receipts + descriptor.get("sources", []),
                "protein_name_match": choices[0], "proxy_target": None}
    except (OSError, ValueError) as error:
        return {**_base(spec), "status": "source_unavailable", "reason": "protein_name_source_unavailable",
                "source_receipts": client.receipts, "error_type": type(error).__name__}
    except Exception as error:
        import requests
        if not isinstance(error, requests.RequestException):
            raise
        return {**_base(spec), "status": "source_unavailable", "reason": "protein_name_source_unavailable",
                "source_receipts": client.receipts, "error_type": type(error).__name__}


def resolve_input(input_spec, *, artifact_store, protein_resolver, source_get=None,
                  structure_resolver=None):
    """Resolve only the supplied entity. A compound never implies a protein."""
    spec = normalize_input(input_spec)
    if spec["kind"] == "bound_complex":
        if structure_resolver is None:
            return {**_base(spec), "reason": "structure_resolution_not_configured"}
        return structure_resolver(spec)
    if spec["kind"] == "protein":
        return _protein(spec, artifact_store, protein_resolver, source_get)
    return _compound(spec, artifact_store, source_get)
