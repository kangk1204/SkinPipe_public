"""Assemble the compound library that will be screened against a target.

Three sources, combined and de-duplicated on InChIKey:

1. The curated cosmetic actives set in `config/cosmetic_actives.yaml`, retained
   as annotated reference compounds. Skin or clinical endpoints do not establish
   affinity, action direction or a positive control against the selected node.
   Legacy potency/evidence annotations remain unreviewed claims until linked to
   actual target/assay observations.
2. ChEMBL ligands for the node being screened, provenance-filtered.
3. Any SDF or SMILES file the user supplies.
"""
from __future__ import annotations

import re

import csv
import copy
import hashlib
import json
import logging
from pathlib import Path
from typing import Iterable, Optional

import yaml
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, inchi
from rdkit.Chem.MolStandardize import rdMolStandardize

from .fetch import chembl, mechanism
from .types import LigandRecord
from . import measurements

#: Mirrors `pipeline.DEFAULT_POTENCY_FLOOR`. Defined here rather than
#: imported because `pipeline` imports this module; the test
#: `test_the_two_default_floors_are_the_same_number` pins them together.
DEFAULT_POTENCY_FLOOR = 5.0

RDLogger.DisableLog("rdApp.*")
log = logging.getLogger(__name__)

_ACTIVES = Path(__file__).resolve().parent.parent / "config" / "cosmetic_actives.yaml"


def _literature_context(entry: dict, definition: Path) -> dict | None:
    """Resolve inspectable external context; never supply an activity value.

    Identifiers must refer to this compound, so a moved or misspelled reference
    cannot borrow another compound's endpoint. The complete review is byte-bound
    and remains inside the quarantined annotation in JSON and CSV delivery.
    """
    ids = entry.get("literature_claim_ids")
    if not ids:
        return None
    if not isinstance(ids, list) or any(not isinstance(i, str) for i in ids) or len(set(ids)) != len(ids):
        raise ValueError("literature claim identifiers must be unique strings")
    path = definition.with_name("curated_literature.json")
    try:
        raw = path.read_bytes()
        data = json.loads(raw)
        if data["schema"] != "scoped-primary-literature-v1":
            raise ValueError("unsupported literature schema")
        claims = {c["id"]: c for c in data["claims"]}
        sources = {s["id"]: s for s in data["sources"]}
        if len(claims) != len(data["claims"]) or len(sources) != len(data["sources"]):
            raise ValueError("duplicate literature definition identifier")
        selected = []
        for key in ids:
            claim = copy.deepcopy(claims[key])
            if claim["compound"] != entry["name"] or claim["status"] != "SOURCE_VERIFIED_SCOPE_LIMITED":
                raise ValueError("literature compound or review status mismatch")
            if claim["project_validation"] != "NOT_EVALUATED":
                raise ValueError("literature context cannot establish project validation")
            claim["source"] = copy.deepcopy(sources[claim["source_id"]])
            selected.append(claim)
    except (OSError, KeyError, TypeError, json.JSONDecodeError) as exc:
        raise ValueError(f"literature context unavailable or inconsistent: {exc}") from exc
    return {"status": "SOURCE_REVIEWED_CONTEXT_ONLY", "claims": selected,
            "source_definition": path.name, "source_definition_sha256": hashlib.sha256(raw).hexdigest(),
            "interpretation": "External study context; target measurements, binding, direction, clinical score and project efficacy are not automatically established"}


def curated_literature_inventory() -> dict:
    """Expose reviewed external context that has no registered candidate."""
    path = _ACTIVES.with_name("curated_literature.json")
    try:
        raw = path.read_bytes()
        literature = json.loads(raw)
        definitions = yaml.safe_load(_ACTIVES.read_bytes())
        registered = {claim for row in definitions["actives"]
                      for claim in row.get("literature_claim_ids", [])}
        claims = {row["id"] for row in literature["claims"]}
    except (OSError, KeyError, TypeError, ValueError) as exc:
        return {"status": "NOT_EVALUATED", "reason": type(exc).__name__}
    return {"status": "SCOPED_CONTEXT_WITH_REGISTRATION_GAPS" if claims - registered else "SCOPED_CONTEXT_REGISTERED",
            "registered_claim_ids": sorted(claims & registered),
            "unregistered_claim_ids": sorted(claims - registered),
            "unknown_registered_claim_ids": sorted(registered - claims),
            "source_sha256": hashlib.sha256(raw).hexdigest(),
            "interpretation": "Unregistered claims are external context only; no candidate, activity or score is created."}


_UNCHARGER = rdMolStandardize.Uncharger()


def standardize(smiles: str) -> Optional[str]:
    """Neutralise an unambiguous organic parent; preserve ambiguous mixtures.

    Size cannot identify the active component of an organic salt. Tosylate is
    larger than methylamine. Multiple organic fragments retain their original
    graph and charges for downstream admission, instead of substituting the
    counterion. One organic fragment plus inorganic ions can be desalted.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fragments = Chem.GetMolFrags(mol, asMols=True)
    if len(fragments) > 1:
        organic = [part for part in fragments if any(a.GetAtomicNum() == 6 for a in part.GetAtoms())]
        if len(organic) != 1:
            log.warning("Ambiguous multi-fragment parent preserved: %s", smiles)
            return Chem.MolToSmiles(mol)
        mol = organic[0]
    try:
        mol = _UNCHARGER.uncharge(mol)
    except Exception:                     # permanently charged species (quats)
        pass
    if mol is None or mol.GetNumHeavyAtoms() == 0:
        return None
    return Chem.MolToSmiles(mol)


def structure_identity(original_smiles: str, standardized_smiles: str, *,
                       source_id: Optional[str] = None, hierarchy: Optional[dict] = None) -> dict:
    """Record exact source and working identities without guessing a salt parent.

    The existing standardizer may remove inorganic counterions around one
    organic component. Multiple organic components stay unresolved; neither a
    largest-fragment heuristic nor a database parent ID selects one silently.
    """
    hierarchy = hierarchy or {}
    original = Chem.MolFromSmiles(original_smiles)
    working = Chem.MolFromSmiles(standardized_smiles)
    raw_parts = Chem.GetMolFrags(original, asMols=True) if original is not None else ()
    parts = Chem.GetMolFrags(working, asMols=True) if working is not None else ()
    organic_parts = sum(any(atom.GetAtomicNum() == 6 for atom in part.GetAtoms())
                        for part in raw_parts)
    needs_review = len(parts) > 1
    return {"input_molecule_id": source_id,
            "parent_molecule_id": hierarchy.get("parent_chembl_id"),
            "active_molecule_id": hierarchy.get("active_chembl_id"),
            "original_smiles": original_smiles, "standardized_smiles": standardized_smiles,
            "input_structure_key": _key(original_smiles),
            "working_structure_key": _key(standardized_smiles),
            "input_fragment_count": len(raw_parts), "input_organic_fragment_count": organic_parts,
            "fragment_count": len(parts), "requires_parent_review": needs_review,
            "parent_selection": "NOT_SUBSTITUTED",
            "standardization": ("INORGANIC_COUNTERIONS_REMOVED_AND_CHARGE_NORMALIZED"
                                if len(raw_parts) > 1 and len(parts) == 1
                                else "RECORDED_GRAPH_WITH_CHARGE_NORMALIZATION"),
            "status": "AMBIGUOUS_MULTICOMPONENT" if needs_review else "SINGLE_COMPONENT",
            "interpretation": "Source identity and parent metadata are retained; an ambiguous multi-component graph is not a single active parent and requires explicit review before ranking or docking."}


def _descriptors(mol: Chem.Mol) -> dict:
    return {"mw": round(Descriptors.MolWt(mol), 2),
            "logp": round(Descriptors.MolLogP(mol), 2),
            "tpsa": round(Descriptors.TPSA(mol), 2),
            "hbd": Descriptors.NumHDonors(mol),
            "hba": Descriptors.NumHAcceptors(mol),
            "heavy_atoms": mol.GetNumHeavyAtoms(),
            "rotatable": Descriptors.NumRotatableBonds(mol)}


def _key(smiles: str) -> Optional[str]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        return inchi.MolToInchiKey(mol)
    except Exception:
        return Chem.MolToSmiles(mol)


def cosmetic_actives(marker: Optional[str] = None, node: Optional[str] = None,
                     path: Optional[Path] = None) -> list[LigandRecord]:
    """Curated topical actives, optionally narrowed to one marker or node."""
    definition = Path(path or _ACTIVES)
    raw = definition.read_bytes()
    data = yaml.safe_load(raw)
    definition_hash = hashlib.sha256(raw).hexdigest()
    out: list[LigandRecord] = []
    for index, a in enumerate(data["actives"]):
        if marker and a["marker"].lower() != marker.lower():
            continue
        if node and a["node"].upper() != node.upper():
            continue
        smi = standardize(a["smiles"]) or a["smiles"]
        mol = Chem.MolFromSmiles(smi)
        props = _descriptors(mol) if mol else {
            k: a.get(k) for k in ("mw", "logp", "tpsa", "hbd", "hba")}
        props["structure_identity"] = structure_identity(a["smiles"], smi,
                                                         source_id=str(a["pubchem_cid"]))
        annotation = copy.deepcopy(a) | {
            "claim_status": "UNREVIEWED_TARGET_ENDPOINT",
            "source_definition": definition.name,
            "source_definition_sha256": definition_hash,
            "source_entry_index": index,
            "interpretation": "Reference annotation; external literature context and legacy scalar potency are not automatically verified measurements against this node"}
        context = _literature_context(a, definition)
        if context is not None:
            annotation["literature_review"] = context
        out.append(LigandRecord(
            ligand_id=f"ACTIVE_{a['pubchem_cid']}", smiles=smi,
            name=a["name"], source="cosmetic_active",
            natural_product=a.get("natural_product"),
            natural_product_sources=[{"source": "curated", "value": a.get("natural_product"),
                                      "basis": a.get("natural_product_basis", "unreviewed"),
                                      "source_urls": a.get("natural_product_source_urls", [])}],
            original_structures=[{"source": "curated", "id": str(a["pubchem_cid"]),
                                  "smiles": a["smiles"]}],
            props=props,
            activity={"marker": a["marker"], "node": a["node"],
                      "evidence": "unknown", "note": a.get("note", ""),
                      "potency_note": a.get("potency_note", ""),
                      "pchembl": None,
                      "curated_claim_policy": "curated-target-claims-v1",
                      "curated_claim_status": "UNREVIEWED_TARGET_ENDPOINT",
                      "curated_claims": [annotation]}))
    return out


def beats(a: dict, cur: dict) -> bool:
    """Should activity `a` replace `cur` as the record kept for one compound?

    Module level, not a closure inside `from_chembl`, so the rule can be tested
    directly instead of only through a mocked ChEMBL fetch. It decides more than
    a number: the winning record also supplies the `provenance` and `assay`
    fields the report prints for that compound, so this is a labelling decision
    as much as a potency one.

    Strongest first, then confirmed-human, then the lowest assay id. The last
    step is what makes it deterministic: `a["pchembl"] > cur["pchembl"]` alone
    kept whichever row ChEMBL happened to serialise first, and pChEMBL is
    quantised to two decimals, so ties are common. The same compound came back
    labelled `human_protein` on one run and `surrogate` on the next from
    identical evidence. Under `trustworthy_only=False` that tie is between a
    human and a surrogate assay, which is the case that must not be settled by
    arrival order.
    """
    ka = (a["pchembl"], bool(a.get("provenance_trustworthy")))
    kc = (cur["pchembl"], bool(cur.get("provenance_trustworthy")))
    if ka != kc:
        return ka > kc
    # Everything that matters is equal; break the tie on the identifier so the
    # choice is a property of the data and not of the response.
    return str(a.get("assay_chembl_id") or "") < str(cur.get("assay_chembl_id") or "")


def from_chembl(node_gene: str, pchembl_min: float = DEFAULT_POTENCY_FLOOR,
                max_records: int = 3000,
                trustworthy_only: bool = True,
                target_identity: Optional[dict] = None,
                preserve_all_queried: bool = False,
                include_without_pchembl: bool = False) -> tuple[list[LigandRecord], dict]:
    """Ligands for a node, with the provenance audit returned alongside.

    The audit is returned rather than logged so the caller can put the warning
    in front of whoever reads the hit list. It also carries
    `target_chembl_id`, which is what lets `assemble` ask ChEMBL which
    direction each of these compounds pushes the node in without resolving the
    gene a second time.
    """
    from .fetch import provenance







    if target_identity is not None:
        if (target_identity.get("status") != "VERIFIED" or target_identity.get("taxon_id") != 9606
                or target_identity.get("reviewed") is not True
                or node_gene != (target_identity.get("gene") or target_identity.get("uniprot"))):
            raise ValueError("Unverified target identity cannot select a ChEMBL library")
        tid = target_identity.get("chembl_target")
    else:
        u = chembl.uniprot_for_gene(node_gene)
        if not u:
            return [], {"library_unavailable": f"no UniProt entry for {node_gene}"}
        tid = chembl.chembl_target_for_accession(u["accession"])
    if not tid:
        return [], {"library_unavailable": f"no human ChEMBL target verified for {node_gene}",
                    "target_identity": target_identity}

    fetch_args = {"include_without_pchembl": True} if include_without_pchembl else {}
    acts, fetch_meta = chembl.fetch_activities(
        tid, pchembl_min=pchembl_min, max_records=max_records,
        target_gene=node_gene, with_meta=True, **fetch_args)
    audit = provenance.audit(acts, node_gene)
    audit.pop("calls", None)
    audit["fetch"] = fetch_meta
    audit["target_chembl_id"] = tid
    if target_identity is not None:
        audit["target_identity"] = target_identity
    used = [a for a in acts if a["provenance_trustworthy"]] if trustworthy_only else acts
    if (trustworthy_only and not used and acts and not preserve_all_queried
            and fetch_meta.get("assay_context_complete") is not False):
        # Every record was thrown away, so the ChEMBL half really is empty --
        # the same state as the two give-ups above, reported under the same key.
        audit["library_unavailable"] = (
            f"every activity for {node_gene} failed provenance QC. "
            "Re-run with trustworthy_only=False only if you accept "
            "surrogate-species data, and label the result as such.")

    # Compound identity and measurements have different grains. Keep every
    # measurement, then choose a recorded-context group independently of potency.
    by_compound: dict[str, list[dict]] = {}
    eligible_ids = {a.get("molecule_chembl_id") for a in used if a.get("molecule_chembl_id")}
    audit["provenance_policy"] = {"trustworthy_only": trustworthy_only,
                                  "preserve_all_queried": preserve_all_queried,
                                  "n_eligible_compound_ids": len(eligible_ids)}
    audit["unassigned_measurements"] = []
    for activity_row in acts:
        cid = activity_row.get("molecule_chembl_id")
        if not cid:
            audit["unassigned_measurements"].append(activity_row)
            continue
        if not preserve_all_queried and cid not in eligible_ids:
            continue
        by_compound.setdefault(cid, []).append(activity_row)
    summaries = {cid: measurements.summarise(rows, require_human_target=trustworthy_only)
                 for cid, rows in by_compound.items()}
    best = {cid: measurements.representative(summary) for cid, summary in summaries.items()}
    audit["measurement_policy"] = measurements.POLICY_VERSION
    audit["n_preserved_measurements"] = sum(s["n_activities"] for s in summaries.values())
    audit["n_queried_compound_ids"] = len({a.get("molecule_chembl_id") for a in acts
                                         if a.get("molecule_chembl_id")})
    audit["unavailable_structures"] = []

    mols, mol_meta = chembl.fetch_molecules(best.keys(), with_meta=True)
    audit["molecules"] = {k: v for k, v in mol_meta.items()
                          if k not in ("missing_ids", "failed_batch_ids", "stale_ids")}
    audit["molecules"]["n_missing"] = len(mol_meta["missing_ids"])
    audit["molecules"]["n_stale"] = len(mol_meta.get("stale_ids", []))
    out: list[LigandRecord] = []
    n_local = 0
    for cid, act in best.items():
        m = mols.get(cid)
        if not m or not m.get("smiles"):
            audit["unavailable_structures"].append({"ligand_id": cid,
                "reason": "structure not returned or blank", **summaries[cid]})
            continue
        smi = standardize(m["smiles"])
        if smi is None:
            log.warning("%s: could not standardise %s", cid, m["smiles"][:60])
            audit["unavailable_structures"].append({"ligand_id": cid,
                "reason": "structure could not be standardised", "original_smiles": m["smiles"],
                **summaries[cid]})
            continue
        rd = Chem.MolFromSmiles(smi)
        props = _descriptors(rd) if rd else {}
        hierarchy = m.get("molecule_hierarchy") or {}
        props["structure_identity"] = structure_identity(m["smiles"], smi,
                                                         source_id=cid, hierarchy=hierarchy)
        # Finding N-20. A compound rescued from the local chemreps archive has
        # a structure and nothing else: `natural_product` arrives as None, and
        # `bool(None)` turned "the file cannot tell us" into "it is not a
        # natural product" -- a claim about the compound that nothing supports.
        # Keep None as None, and put the origin in `source`, which is the one
        # per-candidate field `report._row` carries into the CSV, the JSON and
        # the markdown table. Before this the audit's `local_ids` count was the
        # only trace of the rescue, and a count at the top of a report does not
        # qualify an individual row in the hit list.
        from_local = bool(m.get("partial")) or m.get("source") == "chembl_local"
        raw_np = m.get("natural_product")
        activity = {"pchembl": summaries[cid]["pchembl"], "type": act.get("standard_type"),
                    "provenance_policy_eligible": cid in eligible_ids,
                    "provenance_gate_deferred": preserve_all_queried,
                    "provenance": act["provenance"],









                    "provenance_basis": act.get("provenance_basis") or "unrecorded",
                    "assay": (act.get("assay_description") or "")[:200],
                    "node": node_gene,
                    "provenance_trustworthy": bool(act.get("provenance_trustworthy")),
                    **{k: act[k] for k in provenance.CONTEXT_FIELDS if k in act},
                    **{k: act[k] for k in ("assay_chembl_id", "document_chembl_id",
                         "target_organism", "assay_organism", "assay_description",
                         "assay_cell_type", "assay_tissue", "assay_type",
                         "assay_confidence_score", "target_chembl_id", "fetched_at",
                         "source_url") if k in act}}
        activity.update(summaries[cid])
        if activity["primary_measurement_status"] == "AMBIGUOUS_CONTEXTS":
            activity["assay"] = "multiple equally eligible contexts; inspect measurement groups"
        span = activity.get("pchembl_range")
        if span and span[1] - span[0] >= 1.0:
            activity["potency_disputed"] = (
                f"{activity['n_activities']} measurements span pChEMBL "
                f"{span[0]:.2f}-{span[1]:.2f}; separate contexts and exact-endpoint "
                "group medians are retained, not the strongest observation")
        if from_local:
            n_local += 1
            activity["structure_source"] = "chembl_local"
            # Whether the index that supplied it could be matched to its archive.
            # `chembl_local` serves an unproven index on purpose; a reader judging
            # a provenance line needs to know which kind they have.
            if "index_verified" in m:
                activity["structure_index_verified"] = bool(m["index_verified"])
                activity["structure_index_status"] = m.get("index_status")
            activity["annotation_gap"] = (
                "structure recovered from the local chemreps archive, which "
                "carries no pref_name, max_phase, natural_product flag or "
                "molecule_type; those four are unknown, not absent")
        out.append(LigandRecord(
            ligand_id=cid, smiles=smi, name=m.get("name"),
            source="chembl_local" if from_local else "chembl",
            natural_product=None if raw_np is None else bool(raw_np),
            max_phase=m.get("max_phase"), props=props, activity=activity,
            natural_product_sources=[{"source": "chembl", "value": raw_np,
                                      "basis": "database annotation; not an occurrence claim",
                                      "id": cid}],
            original_structures=[{"source": m.get("source", "chembl"),
                                  "id": cid, "smiles": m["smiles"],
                                  "molecule_hierarchy": hierarchy}]))
    audit["molecules"]["n_local_in_library"] = n_local
    return out, audit


#: Column names a user CSV may carry a registry number under.
_CAS_COLUMNS = ("cas", "CAS", "cas_no", "CAS_NO", "cas_number", "CAS Number",
                "casrn", "CASRN")


def _cas_of(row: dict) -> Optional[str]:
    """The registry number a user CSV supplied, if it supplied one."""
    row = {str(k).strip().casefold(): v for k, v in row.items()}
    for k in _CAS_COLUMNS:
        key = k.casefold()
        v = (row.get(key) or "").strip() if isinstance(row.get(key), str) else None
        if v:
            return v
    return None


def from_file(path: Path, source: str = "user", *, with_meta: bool = False):
    """Load an SDF, SMILES (.smi) or CSV with a `smiles` column."""
    path = Path(path)
    out: list[LigandRecord] = []
    audit = {"path": path.name, "source": source, "encoding": None,
             "n_rows": 0, "accepted": 0, "rejected_rows": [], "complete": True}
    def reject(index, reason):
        audit["rejected_rows"].append({"row_index": index, "reason": reason})
        audit["complete"] = False
        log.warning("%s row %d rejected: %s", path.name, index, reason)
    def result():
        audit["accepted"] = len(out)
        audit["n_rejected"] = len(audit["rejected_rows"])
        audit["status"] = "COMPLETE" if audit["complete"] else "PARTIAL_IMPORT"
        if not out:
            audit.update(complete=False, status="EMPTY_INPUT")
            log.warning("%s: no usable structures imported; inspect the input audit", path.name)
        return (out, audit) if with_meta else out
    if path.suffix.lower() == ".sdf":
        for i, mol in enumerate(Chem.SDMolSupplier(str(path))):
            audit["n_rows"] += 1
            if mol is None:
                reject(i, "invalid SDF record")
                continue
            name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"{path.stem}_{i}"
            smi = standardize(Chem.MolToSmiles(mol))
            if smi is None:
                reject(i, "structure could not be standardized")
                continue
            m2 = Chem.MolFromSmiles(smi)
            props = _descriptors(m2)
            props["structure_identity"] = structure_identity(Chem.MolToSmiles(mol), smi,
                                                             source_id=f"{source}_{i}")
            # An SDF carries its metadata as tagged properties, and only
            # `_Name` was ever read -- so a file that states a registry number
            # lost it at import, before anything could use or check it.
            cas = _cas_of({k: mol.GetProp(k) for k in mol.GetPropNames()})
            if cas:
                props["cas_supplied"] = cas
            out.append(LigandRecord(ligand_id=f"{source}_{i}", smiles=smi,
                                    name=name, source=source, props=props,
                                    original_structures=[{"source": source, "path": path.name,
                                                          "smiles": Chem.MolToSmiles(mol)}]))
        return result()

    rows: list[tuple[str, str, Optional[str], dict]] = []
    if path.suffix.lower() in (".csv", ".tsv"):
        delim = "\t" if path.suffix.lower() == ".tsv" else ","
        audit["encoding"] = "utf-8-sig"
        with path.open(newline="", encoding="utf-8-sig") as fh:
            reader = csv.DictReader(fh, delimiter=delim)
            headers = reader.fieldnames
            normalized = [str(h).strip().casefold() for h in headers or []]
            header_error = ("CSV/TSV requires a smiles column (case-insensitive)" if "smiles" not in normalized
                            else "duplicate CSV/TSV headers after normalization" if len(normalized) != len(set(normalized))
                            else None)
            if header_error:
                if not with_meta:
                    raise ValueError(f"{path.name}: {header_error}")
                audit.update(original_headers=headers, complete=False, status="INVALID_HEADER",
                             file_error=header_error, row_count_status="NOT_READ_INVALID_HEADER", n_rejected=0)
                return [], audit
            audit["original_headers"] = headers
            for i, raw_row in enumerate(reader):
                audit["n_rows"] += 1
                if None in raw_row:
                    reject(i, "row has more values than the header")
                    continue
                row = {str(k).strip().casefold(): v for k, v in raw_row.items()}
                smi = (row.get("smiles") or "").strip()
                if smi:
                    rows.append((smi, row.get("name") or row.get("id") or f"{path.stem}_{i}",
                                 _cas_of(row), dict(row, _input_row_index=i)))
                else:
                    reject(i, "blank smiles")
    else:                                            # .smi / .txt
        audit["encoding"] = "utf-8-sig"
        for i, line in enumerate(path.read_text(encoding="utf-8-sig").splitlines()):
            parts = line.split()
            if parts:
                audit["n_rows"] += 1
                rows.append((parts[0], parts[1] if len(parts) > 1 else f"{path.stem}_{i}",
                             None, {"_input_row_index": i}))

    for i, (smi, name, cas, metadata) in enumerate(rows):
        std = standardize(smi)
        if std is None:
            reject(metadata.get("_input_row_index", i), "invalid SMILES or empty structure")
            continue
        mol = Chem.MolFromSmiles(std)
        props = _descriptors(mol)
        props["structure_identity"] = structure_identity(smi, std, source_id=f"{source}_{i}")
        # A CAS the user supplied is evidence they already have; dropping it
        # made the pipeline look the same compound up again and, where the
        # lookup found nothing, ship a row with no registry number at all --
        # having been handed one. It is recorded as *theirs*, not as a checked
        # identifier: nothing here has verified that the number and the
        # structure describe the same substance.
        if cas:
            props["cas_supplied"] = cas
        out.append(LigandRecord(ligand_id=f"{source}_{i}", smiles=std,
                                name=name, source=source, props=props,
                                parent_id=metadata.get("parent_id") or None,
                                derivation=metadata.get("derivation") or None,
                                parent_smiles=metadata.get("parent_smiles") or None,
                                parent_metadata_source="input file; not independently verified"
                                         if metadata.get("parent_smiles") else None,
                                original_structures=[{"source": source, "path": path.name,
                                                      "smiles": smi}]))
    return result()


def assemble(node_gene: Optional[str] = None, marker: Optional[str] = None,
             extra_files: Optional[Iterable[Path]] = None,
             include_actives: bool = True,
             pchembl_min: float = DEFAULT_POTENCY_FLOOR,
             max_chembl: int = 3000, trustworthy_only: bool = True,
             include_offnode_references: bool = False,
             target_identity: Optional[dict] = None,
             preserve_all_queried: bool = False,
             include_without_pchembl: bool = False
             ) -> tuple[list[LigandRecord], dict]:
    'Build and de-duplicate the screening library.'
    records: list[LigandRecord] = []
    audit: dict = {}
    if include_actives:
        on_node = cosmetic_actives(node=node_gene) if node_gene else cosmetic_actives()
        records += on_node
        if include_offnode_references and node_gene:
            on_ids = {r.ligand_id for r in on_node}
            for r in cosmetic_actives():
                if r.ligand_id not in on_ids:
                    r.source = "reference_offnode"
                    records.append(r)
    if node_gene:
        identity_args = {"target_identity": target_identity} if target_identity is not None else {}
        if preserve_all_queried:
            identity_args["preserve_all_queried"] = True
        if include_without_pchembl:
            identity_args["include_without_pchembl"] = True
        chembl_recs, audit = from_chembl(node_gene, pchembl_min=pchembl_min,
                                         max_records=max_chembl,
                                         trustworthy_only=trustworthy_only, **identity_args)
        records += chembl_recs
    for i, f in enumerate(extra_files or []):
        # Namespaced per file: from_file numbers rows from 0, so two supplied
        # files both produce "user_0". Downstream stages key skin assessments
        # and ADMET profiles by ligand_id, and a collision there silently hands
        # one molecule another molecule's permeation and hazard numbers.
        imported, receipt = from_file(Path(f), source=f"user{i}" if i else "user", with_meta=True)
        records += imported
        audit.setdefault("input_files", []).append(receipt)

    seen: dict[str, LigandRecord] = {}
    merged_ids: dict[str, list[str]] = {}
    def precedence(record):
        # Explicit user uploads have a deliberate priority. An unknown source
        # must not acquire that priority just because it is unrecognised.
        source_rank = {"cosmetic_active": 5, "chembl": 3,
                       "reference_offnode": 2, "chembl_local": 1}.get(
                           record.source, 4 if re.fullmatch(r"user\d*", record.source) else 0)
        completeness = sum(value is not None for value in
                           (record.name, record.max_phase, record.natural_product))
        return (-source_rank, -completeness, record.ligand_id, record.smiles)
    # Choose from original records before merges can change their completeness.
    # This makes both representative identity and annotation consensus invariant
    # to ChEMBL's row order, including groups with more than two records.
    for r in sorted(records, key=precedence):
        k = _key(r.smiles)
        if k is None:
            continue
        prior = seen.get(k)
        if prior is None:
            seen[k] = r
            merged_ids.setdefault(k, [])
            continue
        # a curated on-node active beats a bare ChEMBL row for the same
        # structure -- including one rescued from the local archive, which is
        # the *least* annotated record of the three and must never displace the
        # curated one (finding N-20 renamed its source; this list has to follow)
        dropped = r






        keeper = seen[k]
        measurements.merge_activity(keeper.activity, dropped.activity)
        keeper.original_structures += dropped.original_structures
        for record in (keeper, dropped):
            if record.natural_product is not None and not record.natural_product_sources:
                record.natural_product_sources.append({"source": record.source,
                    "value": record.natural_product, "basis": "original_record_annotation"})
        keeper.natural_product_sources += dropped.natural_product_sources
        reviewed = [x for x in keeper.natural_product_sources
                    if x.get("source") == "curated" and x.get("value") is not None
                    and x.get("basis") != "unreviewed"]
        values = {x["value"] for x in reviewed}
        if not reviewed:
            values = {x["value"] for x in keeper.natural_product_sources
                      if x.get("value") is not None}
        keeper.natural_product = next(iter(values)) if len(values) == 1 else None
        for attribute in ("name", "max_phase"):
            if getattr(keeper, attribute) is None:
                setattr(keeper, attribute, getattr(dropped, attribute))
        merged_ids[k].append(dropped.ligand_id)
        merged_ids[k] += list(dropped.activity.get("duplicate_ids") or [])
    for k, keeper in seen.items():
        ids = list(keeper.activity.get("duplicate_ids") or []) + merged_ids.get(k, [])
        ids = sorted({i for i in ids if i and i != keeper.ligand_id})
        if ids:
            keeper.activity["duplicate_ids"] = ids
    deduped = _uniquify_ids(list(seen.values()))
    # Direction, last: which way each compound pushes the node. Done here, on
    # the finished library, so one mechanism fetch covers curated actives,
    # ChEMBL rows and user-supplied files alike -- a curated active carries an
    # `ACTIVE_<cid>` id that no ChEMBL query would ever return, and it is
    # exactly the kind of record (capsaicin) that gets mistaken for a positive
    # control. Potency says how hard a compound pulls; this says which way.
    audit["mechanism"] = mechanism.annotate(
        deduped, target_chembl_id=audit.get("target_chembl_id"), node_gene=node_gene)
    audit["curated_literature_registration"] = curated_literature_inventory()
    log.info("library: %d records -> %d unique structures, %d with a known direction",
             len(records), len(deduped), audit["mechanism"]["n_with_direction"])
    return deduped, audit


def _uniquify_ids(records: list[LigandRecord]) -> list[LigandRecord]:
    """Guarantee distinct `ligand_id`s, renaming duplicates rather than losing them.

    Every stage after this one builds `{ligand_id: ...}` maps. A duplicate id
    does not raise; it silently collapses two molecules into one entry and
    scores both on whichever assessment landed last.
    """
    # A generated suffix must not take a later record's original identifier.
    # A,A,A#2 previously became A,A#2,A#2 and collided downstream.
    reserved = {r.ligand_id for r in records}
    used: set[str] = set()
    next_suffix: dict[str, int] = {}
    for r in records:
        if r.ligand_id in used:
            original = r.ligand_id
            number = next_suffix.get(original, 2)
            new_id = f"{original}#{number}"
            while new_id in reserved or new_id in used:
                number += 1
                new_id = f"{original}#{number}"
            next_suffix[original] = number + 1
            log.warning("duplicate ligand_id %s; renamed to %s", r.ligand_id, new_id)
            r.ligand_id = new_id
        used.add(r.ligand_id)
    return records
