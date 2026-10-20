"""Owner-bound coordinate uploads and exact bound-ligand input identity.

Coordinate files supply observed atoms, never chemical bond orders. A seed
requires an explicit chemical graph from an official CCD or the user's exact
SMILES/SDF. Retained components, source coordinates and selection decisions are
content-addressed; unsupported preparation contexts are explicit holds.
"""
from __future__ import annotations

import base64
import binascii
from collections import defaultdict
from copy import deepcopy
from io import BytesIO, StringIO
import json
import math
from pathlib import Path
import re

from rdkit import Chem
from rdkit.Chem import AllChem

from .artifacts import ArtifactStore, atomic_bytes, canonical_json
from .identity import chemical_identity, digest
from .targets import now

MAX_UPLOAD_BYTES = 50 * 1024 * 1024
MAX_ATOMS = 300000
WATER = {"HOH", "WAT", "DOD"}
METALS = {"LI", "NA", "K", "RB", "CS", "MG", "CA", "SR", "BA", "ZN", "FE", "CU", "MN", "CO", "NI", "CD", "HG"}
STANDARD_AA = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "SEC", "PYL"}
COFACTORS = {"HEM", "HEC", "NAD", "NAP", "NDP", "FAD", "FMN", "PLP", "SAM", "SAH", "COA", "ATP", "ADP", "GTP", "GDP"}


def normalize_structure_input(body):
    from .research_inputs import INPUT_SCHEMA

    allowed = {"kind", "format", "value", "selection", "ligand_smiles", "ligand_sdf", "desired_action"}
    if not isinstance(body, dict) or set(body) - allowed or body.get("kind") != "bound_complex":
        raise ValueError("Unsupported bound complex input fields")
    form, value = body.get("format"), body.get("value")
    if not isinstance(form, str) or form not in {"pdb_id", "upload"} or not isinstance(value, str):
        raise ValueError("Bound complexes accept pdb_id or upload artifact identity")
    if form == "upload" and value.startswith("upload:"):
        value = value[len("upload:"):]
    if form == "pdb_id" and not re.fullmatch(r"[0-9][A-Za-z0-9]{3}", value):
        raise ValueError("Enter an exact four-character PDB identifier")
    if form == "upload" and not re.fullmatch(r"[a-f0-9]{64}", value):
        raise ValueError("An immutable uploaded structure SHA256 is required")
    if not isinstance(body.get("desired_action", "any"), str) or body.get("desired_action", "any") not in {"any", "inhibit", "activate"}:
        raise ValueError("Unsupported desired action")
    selection = body.get("selection", {})
    keys = {"model_id", "assembly_id", "protein_chains", "ligand_chain", "ligand_residue", "ligand_code", "altloc"}
    if not isinstance(selection, dict) or set(selection) - keys:
        raise ValueError("Unsupported structure selection fields")
    for key, item in selection.items():
        if key == "protein_chains":
            if (not isinstance(item, list) or not item or len(item) > 62
                    or any(not isinstance(c, str) or len(c) > 16 or re.search(r"[\x00-\x20]", c) for c in item)
                    or len(set(item)) != len(item)):
                raise ValueError("Select distinct protein chain identifiers")
        elif isinstance(item, bool) or not isinstance(item, (str, int)) or len(str(item)) > 64 or re.search(r"[\x00-\x1f]", str(item)):
            raise ValueError("Invalid bounded structure selection")
    result = {"schema": INPUT_SCHEMA, "kind": "bound_complex", "format": form, "value": value,
              "original_input": deepcopy(body), "selection": deepcopy(selection),
              "desired_action": body.get("desired_action", "any"), "species_defaulted": False,
              "taxon_id": None}
    for key in ("ligand_smiles", "ligand_sdf"):
        if key not in body:
            continue
        item = body[key]
        if key == "ligand_sdf" and isinstance(item, str) and item.startswith("upload:"):
            item = item[len("upload:"):]
        if not isinstance(item, str) or not item or len(item) > 16384 or any(c in item for c in "\r\n\x00"):
            raise ValueError("Supplemental ligand chemistry must be a bounded SMILES or uploaded SDF SHA256")
        if key == "ligand_sdf" and not re.fullmatch(r"[a-f0-9]{64}", item):
            raise ValueError("Supplemental SDF must be an immutable uploaded artifact")
        result[key] = item
    if "ligand_smiles" in body and "ligand_sdf" in body:
        raise ValueError("Supply one authoritative ligand chemical graph")
    return result


def _column(data, key, default=None, count=None):
    values = data.get(key)
    if values is None:
        return [default] * count if count is not None else []
    values = values if isinstance(values, list) else [values]
    if count is not None and len(values) != count:
        raise ValueError("Unequal mmCIF coordinate columns")
    return values


def _cif_data(text):
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    return MMCIF2Dict(StringIO(text))


def _clear(value):
    return "" if value in {".", "?", None} else str(value)


def parse_structure(text, form):
    """Read every model and alternate location without picking a first copy."""
    from Bio.PDB.Polypeptide import is_aa
    atoms, links, data = [], [], {}
    if form == "mmcif":
        data = _cif_data(text)
        groups = _column(data, "_atom_site.group_PDB")
        count = len(groups)
        if not count or count > MAX_ATOMS:
            raise ValueError("Coordinate atom count is missing or exceeds the input limit")
        def col(auth, label=None, default=None):
            key = auth if auth in data or label is None else label
            return _column(data, key, default, count)
        columns = [col("_atom_site.id"), col("_atom_site.auth_atom_id", "_atom_site.label_atom_id"),
            col("_atom_site.type_symbol"), col("_atom_site.auth_comp_id", "_atom_site.label_comp_id"),
            col("_atom_site.auth_asym_id", "_atom_site.label_asym_id"),
            col("_atom_site.auth_seq_id", "_atom_site.label_seq_id"), col("_atom_site.pdbx_PDB_ins_code", default=""),
            col("_atom_site.label_alt_id", default=""), col("_atom_site.pdbx_PDB_model_num", default="1"),
            col("_atom_site.Cartn_x"), col("_atom_site.Cartn_y"), col("_atom_site.Cartn_z"),
            col("_atom_site.occupancy", default="1")]
        for group, row in zip(groups, zip(*columns)):
            serial, name, element, code, chain, residue, insertion, alt, model, x, y, z, occupancy = row
            atoms.append({"serial": str(serial), "name": str(name), "element": str(element).upper(),
                "code": str(code), "chain": _clear(chain), "residue": str(residue) + _clear(insertion),
                "altloc": _clear(alt), "model": str(model), "xyz": [float(x), float(y), float(z)],
                "occupancy": float(occupancy), "record": group,
                "protein": is_aa(str(code).upper(), standard=False) or (group == "ATOM" and str(code).upper() == "UNK")})
        types = _column(data, "_struct_conn.conn_type_id")
        for index, kind in enumerate(types):
            item = {"type": kind}
            for n in (1, 2):
                for key in ("auth_asym_id", "auth_seq_id", "label_comp_id", "label_atom_id"):
                    values = _column(data, f"_struct_conn.ptnr{n}_{key}")
                    item[f"ptnr{n}_{key}"] = values[index] if index < len(values) else None
            links.append(item)
    elif form == "pdb":
        model = "1"
        for line in text.splitlines():
            if line.startswith("MODEL "):
                model = line[10:14].strip() or "1"
            if line.startswith(("LINK  ", "CONECT", "SSBOND")):
                links.append({"type": line[:6].strip(), "record": line})
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            if len(line) < 54:
                raise ValueError("Truncated PDB coordinate record")
            code, name = line[17:20].strip(), line[12:16].strip()
            element = line[76:78].strip().upper() if len(line) >= 78 else ""
            if not element:
                element = re.sub(r"[0-9]", "", name)[0:1].upper()
            atoms.append({"serial": line[6:11].strip(), "name": name, "element": element,
                "code": code, "chain": line[21:22].strip(), "residue": line[22:27].strip(),
                "altloc": line[16:17].strip(), "model": model,
                "xyz": [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                "occupancy": float(line[54:60].strip() or "1"), "record": line[:6].strip(),
                "protein": is_aa(code, standard=False) or (line.startswith("ATOM  ") and code == "UNK"), "source_line": line})
            if len(atoms) > MAX_ATOMS:
                raise ValueError("Coordinate atom count exceeds the input limit")
        if not atoms:
            raise ValueError("No PDB coordinate records")
    else:
        raise ValueError("Expected PDB or mmCIF coordinates")
    if any(not all(math.isfinite(v) and abs(v) < 1000000 for v in atom["xyz"]) for atom in atoms):
        raise ValueError("Nonfinite or excessive coordinates")
    keys = [(a["model"], a["chain"], a["residue"], a["code"], a["name"], a["altloc"]) for a in atoms]
    if len(set(keys)) != len(keys):
        raise ValueError("Duplicate atom identities in coordinates")
    proteins, ligands = defaultdict(list), defaultdict(list)
    for atom in atoms:
        if atom["protein"]:
            proteins[(atom["model"], atom["chain"])].append(atom)
        elif atom["code"].upper() not in WATER:
            ligands[(atom["model"], atom["chain"], atom["residue"], atom["code"])].append(atom)
    inventory = {"models": sorted({a["model"] for a in atoms}),
        "protein_chains": [{"model_id": model, "chain_id": chain,
            "residue_count": len({a["residue"] for a in rows}), "atom_count": len(rows)}
            for (model, chain), rows in sorted(proteins.items())],
        "ligands": [{"model_id": model, "chain_id": chain, "residue_id": residue, "code": code,
            "altlocs": sorted({a["altloc"] for a in rows if a["altloc"]}),
            "heavy_atom_count": sum(a["element"] not in {"H", "D"} for a in rows),
            "elements": sorted({a["element"] for a in rows}), "metal": any(a["element"] in METALS for a in rows),
            "likely_cofactor": code in COFACTORS}
            for (model, chain, residue, code), rows in sorted(ligands.items())],
        "assembly_ids": _column(data, "_pdbx_struct_assembly.id"),
        "water_atom_count": sum(a["code"].upper() in WATER for a in atoms),
        "atom_count": len(atoms), "connections": links}
    return atoms, inventory, data


def _pdb(atoms):
    lines = []
    for serial, atom in enumerate(atoms, 1):
        match = re.fullmatch(r"(-?[0-9]+)([A-Za-z]?)", atom["residue"])
        if (serial > 99999 or len(atom["chain"]) > 1 or len(atom["code"]) > 3 or not match
                or not -999 <= int(match[1]) <= 9999):
            raise ValueError("Selected coordinates cannot be represented without PDB identity loss")
        if any(v < -999.999 or v > 9999.999 for v in atom["xyz"]):
            raise ValueError("Coordinates cannot fit PDB fields")
        if atom.get("source_line"):
            lines.append(atom["source_line"])
            continue
        name = f" {atom['name']:<3}" if len(atom["element"]) == 1 and len(atom["name"]) < 4 else f"{atom['name']:<4}"
        lines.append(f"{atom['record']:<6}{serial:5d} {name}{atom['altloc'][:1]:1s}{atom['code']:>3s} {atom['chain']:1s}{int(match[1]):4d}{match[2]:1s}   "
                     + "".join(f"{v:8.3f}" for v in atom["xyz"]) + f"{atom['occupancy']:6.2f}{0:6.2f}          {atom['element']:>2s}")
    return "\n".join(lines) + "\nEND\n"


class ResearchStructures:
    def __init__(self, root, artifacts=None, *, source_get=None):
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)
        self.artifacts = artifacts or ArtifactStore(self.root / "artifacts")
        self.source_get = source_get

    def _owner_path(self, owner, sha):
        if not isinstance(owner, str) or not owner:
            raise ValueError("An authenticated upload owner is required")
        return self.root / "uploads" / digest(owner) / (sha + ".json")

    def upload(self, body, owner):
        self._owner_path(owner, "owner-check")
        allowed = {"format", "filename", "text", "content_base64"}
        if not isinstance(body, dict) or set(body) - allowed or not isinstance(body.get("format"), str) or body.get("format") not in {"pdb", "mmcif", "sdf"}:
            raise ValueError("Upload one PDB, mmCIF or supplemental SDF")
        if ("text" in body) == ("content_base64" in body):
            raise ValueError("Upload requires exactly one text or content_base64 field")
        filename = body.get("filename", "structure." + body["format"])
        if not isinstance(filename, str) or len(filename) > 255 or Path(filename).name != filename or any(c in filename for c in "\\\x00\r\n"):
            raise ValueError("Filename is a display basename, never a filesystem path")
        if "text" in body:
            if not isinstance(body["text"], str) or len(body["text"]) > MAX_UPLOAD_BYTES:
                raise ValueError("Upload exceeds 50 MiB decoded")
            raw = body["text"].encode("utf-8")
        else:
            value = body["content_base64"]
            if not isinstance(value, str) or len(value) > 4 * ((MAX_UPLOAD_BYTES + 2) // 3):
                raise ValueError("Upload exceeds 50 MiB decoded")
            try:
                raw = base64.b64decode(value, validate=True)
            except (ValueError, binascii.Error) as error:
                raise ValueError("Upload is not strict base64") from error
        if not raw or len(raw) > MAX_UPLOAD_BYTES:
            raise ValueError("Upload must contain at most 50 MiB decoded")
        try:
            text = raw.decode("utf-8")
        except UnicodeDecodeError as error:
            raise ValueError("Coordinate uploads must be UTF-8 text") from error
        if "\x00" in text:
            raise ValueError("Binary coordinate uploads are unsupported")
        if body["format"] == "sdf":
            mol = _read_sdf(text)
            inventory = {"format": "sdf", "heavy_atoms": mol.GetNumHeavyAtoms(), "canonical_smiles": Chem.MolToSmiles(mol)}
        else:
            _, inventory, _ = parse_structure(text, body["format"])
        receipt = self.artifacts.put_bytes(raw, media_type="chemical/x-mdl-sdfile" if body["format"] == "sdf" else "chemical/x-pdb" if body["format"] == "pdb" else "chemical/x-mmcif")
        record = {"schema": "research-structure-upload-v1", "owner_sha256": digest(owner),
            "artifact_sha256": receipt["sha256"], "size_bytes": len(raw), "format": body["format"],
            "filename": filename, "inventory": inventory, "uploaded_at": now()}
        path = self._owner_path(owner, receipt["sha256"])
        if path.exists():
            existing = json.loads(path.read_bytes())
            if existing["format"] != record["format"]:
                raise ValueError("Same uploaded bytes cannot be reinterpreted as another format")
            return {**existing, "upload_id": "upload:" + receipt["sha256"]}
        atomic_bytes(path, canonical_json(record))
        return {**record, "upload_id": "upload:" + receipt["sha256"]}

    def _uploaded(self, sha, owner):
        if not isinstance(sha, str) or not re.fullmatch(r"[0-9a-f]{64}", sha):
            raise ValueError("Invalid immutable upload identity")
        path = self._owner_path(owner, sha)
        if not path.is_file():
            raise ValueError("Upload is not available to this owner")
        record = json.loads(path.read_bytes())
        if record.get("owner_sha256") != digest(owner) or record.get("artifact_sha256") != sha:
            raise ValueError("Upload ownership receipt differs")
        return self.artifacts.verified_path(sha).read_text(), record

    def _request(self, url, params, receipts):
        import requests
        # Preserve a failed receipt and retry only the same official source.
        # Invalid content, redirects and HTTP errors still fail closed.
        for attempt in range(2):
            try:
                return self._request_once(url, params, receipts)
            except (requests.Timeout, requests.ConnectionError):
                if attempt:
                    raise

    def _request_once(self, url, params, receipts):
        import requests
        if not (url.startswith("https://files.rcsb.org/download/") or url.startswith("https://files.rcsb.org/ligands/download/")
                or url.startswith("https://models.rcsb.org/v1/")):
            raise ValueError("Only official RCSB coordinate and chemistry sources are allowed")
        receipt = {"url": url, "params": params or {}, "retrieved_at": now(), "status": "source_unavailable", "raw_sha256": None}
        try:
            if self.source_get:
                response = self.source_get(url, params or {})
                raw = response.encode() if isinstance(response, str) else bytes(response)
            else:
                with requests.get(url, params=params, timeout=(5, 25), allow_redirects=False, stream=True) as response:
                    response.raise_for_status()
                    if response.status_code != 200:
                        raise ValueError("RCSB requires a complete HTTP 200 response")
                    raw = bytearray()
                    for block in response.iter_content(65536):
                        raw.extend(block)
                        if len(raw) > MAX_UPLOAD_BYTES:
                            raise ValueError("Public structure exceeds 50 MiB")
                    raw = bytes(raw)
            if len(raw) > MAX_UPLOAD_BYTES or not raw or b"\x00" in raw:
                raise ValueError("Public structure is empty, binary or too large")
            text = raw.decode("utf-8")
            receipt.update(status="complete", raw_sha256=self.artifacts.put_bytes(raw)["sha256"])
            return text, receipt["raw_sha256"]
        finally:
            receipts.append(receipt)

    def resolve(self, input_spec, owner):
        from .research_inputs import normalize_input
        spec = normalize_input(input_spec)
        if spec["kind"] != "bound_complex":
            raise ValueError("Expected bound_complex input")
        receipts = []
        result = {"schema": "bound-complex-descriptor-v1", "kind": "bound_complex", "input": spec,
            "status": "selection_required", "target_id": None, "taxon_id": None, "actual_protein_id": None,
            "desired_action": spec["desired_action"], "proxy_target": None, "source_receipts": receipts,
            "limitations": ["A bound pose provides a structural hypothesis; it does not establish activity or skin benefit."],
            "docking_hold_reasons": [], "name": spec["value"], "compute_readiness": "identity_pending"}
        try:
            if spec["format"] == "upload":
                text, uploaded = self._uploaded(spec["value"], owner)
                form, structure_sha = uploaded["format"], spec["value"]
                if form not in {"pdb", "mmcif"}:
                    raise ValueError("Bound structure requires coordinate upload")
                result["upload_receipt"] = uploaded
            else:
                text, structure_sha = self._request("https://files.rcsb.org/download/" + spec["value"].upper() + ".cif", {}, receipts)
                form = "mmcif"
            atoms, inventory, data = parse_structure(text, form)
            result.update(inventory=inventory, choices=inventory["ligands"], original_structure_sha256=structure_sha,
                          original_structure_format=form)
            selection = deepcopy(spec["selection"])
            model = str(selection.get("model_id", inventory["models"][0] if len(inventory["models"]) == 1 else ""))
            if model not in inventory["models"]:
                return {**result, "reason": "model_selection_required"}
            chains = [r["chain_id"] for r in inventory["protein_chains"] if r["model_id"] == model]
            selected_chains = selection.get("protein_chains", chains if len(chains) == 1 else [])
            if not selected_chains or any(c not in chains for c in selected_chains):
                return {**result, "reason": "protein_chain_selection_required"}
            candidates = [r for r in inventory["ligands"] if r["model_id"] == model
                          and "C" in r["elements"] and r["heavy_atom_count"] >= 3]
            for field, key in (("ligand_chain", "chain_id"), ("ligand_residue", "residue_id"), ("ligand_code", "code")):
                if field in selection:
                    candidates = [r for r in candidates if str(r[key]) == str(selection[field])]
            if len(candidates) != 1:
                return {**result, "reason": "bound_ligand_selection_required", "choices": candidates}
            ligand = candidates[0]
            alternate = str(selection.get("altloc", ligand["altlocs"][0] if len(ligand["altlocs"]) == 1 else ""))
            if len(ligand["altlocs"]) > 1 and alternate not in ligand["altlocs"]:
                return {**result, "reason": "ligand_altloc_selection_required"}
            if "altloc" in selection and alternate and alternate not in ligand["altlocs"]:
                return {**result, "reason": "ligand_altloc_not_found"}
            selected_atoms = [a for a in atoms if a["model"] == model and a["chain"] == ligand["chain_id"]
                and a["residue"] == ligand["residue_id"] and a["code"] == ligand["code"] and a["altloc"] in {"", alternate}
                and a["element"] not in {"H", "D"}]
            if len({a["name"] for a in selected_atoms}) != len(selected_atoms):
                return {**result, "reason": "ligand_atom_altloc_conflict"}
            selected_protein = [a for a in atoms if a["model"] == model and a["chain"] in selected_chains and a["protein"]]
            source_selected_protein = deepcopy(selected_protein)
            protein_altlocs = {a["altloc"] for a in selected_protein if a["altloc"]}
            if protein_altlocs:
                result["docking_hold_reasons"].append("protein_alternate_locations_require_preparation_decision")
            assemblies = inventory["assembly_ids"]
            if "assembly_id" not in selection and len(assemblies) > 1:
                return {**result, "reason": "assembly_selection_required"}
            assembly = str(selection.get("assembly_id", assemblies[0] if len(assemblies) == 1 else "as_uploaded" if spec["format"] == "upload" else "asymmetric_unit"))
            if assembly not in {"as_uploaded", "asymmetric_unit"}:
                if assembly not in inventory["assembly_ids"]:
                    return {**result, "reason": "assembly_selection_not_found"}
                try:
                    from ..structure.assembly import build_context
                    work = self.root / "assembly-work" / digest([structure_sha, assembly, selected_chains])
                    work.mkdir(parents=True, exist_ok=True)
                    source = work / "source.cif"
                    atomic_bytes(source, text.encode())
                    context = build_context(source, _pdb([a for a in atoms if a["model"] == model]), selected_chains[0], assembly)
                    assembly_atoms, assembly_inventory, _ = parse_structure(context.pdb_text, "pdb")
                    selected_protein = [a for a in assembly_atoms if a["protein"]]
                    result["assembly_context"] = context.metadata
                    result["assembly_inventory"] = assembly_inventory
                    result["assembly_structure_sha256"] = self.artifacts.put_bytes(context.pdb_text.encode(), media_type="chemical/x-pdb")["sha256"]
                    result["assembly_metadata_sha256"] = self.artifacts.put_bytes(canonical_json(context.metadata))["sha256"]
                except (ValueError, RuntimeError, KeyError, IndexError) as error:
                    result["docking_hold_reasons"].append("deposited_assembly_expansion_required")
                    result["assembly_diagnostic"] = str(error)[:300]
            if any(a["code"] not in STANDARD_AA for a in selected_protein):
                result["docking_hold_reasons"].append("modified_polymer_requires_supported_preparation")
            selection.update(model_id=model, protein_chains=selected_chains, ligand_chain=ligand["chain_id"],
                ligand_residue=ligand["residue_id"], ligand_code=ligand["code"], altloc=alternate, assembly_id=assembly)
            result["selection"] = selection
            seed, chemistry_receipts = self._ligand(spec, selected_atoms, ligand, owner, receipts)
            result["chemistry_receipts"] = chemistry_receipts
            if seed is None:
                return {**result, "reason": "ligand_exact_chemistry_required", "compute_readiness": "chemistry_pending"}
            identity = chemical_identity(Chem.MolToCXSmiles(seed))
            if identity["parent_ambiguous"]:
                return {**result, "reason": "bound_ligand_multiple_organic_components"}
            writer = Chem.MolWriterParams()
            writer.forceV3000, writer.precision = True, 15
            observed = Chem.MolToMolBlock(seed, writer) + "\n$$$$\n"
            seed_sha = self.artifacts.put_bytes(observed.encode(), media_type="chemical/x-mdl-sdfile")["sha256"]
            protein_atoms = [{k: v for k, v in atom.items() if k != "source_line"} for atom in source_selected_protein if atom["element"] not in {"H", "D"}]
            protein_coords_sha = self.artifacts.put_bytes(canonical_json(protein_atoms), media_type="application/json")["sha256"]
            try:
                if any(abs(v - round(v, 3)) > 1e-9 for a in selected_protein for v in a["xyz"]):
                    raise ValueError("Source precision exceeds the native PDB coordinate fields")
                protein_text = _pdb(selected_protein)
                protein_sha = self.artifacts.put_bytes(protein_text.encode(), media_type="chemical/x-pdb")["sha256"]
            except ValueError:
                protein_sha = structure_sha
                result["docking_hold_reasons"].append("lossless_native_structure_format_support_required")
            retained = [r for r in inventory["ligands"] if r["model_id"] == model and r != ligand]
            result["retained_components"] = retained
            if any(r["metal"] or r["likely_cofactor"] for r in retained):
                result["docking_hold_reasons"].append("metal_or_cofactor_context_requires_supported_preparation")
            if any(r["chain_id"] in selected_chains and not r["metal"] and not r["likely_cofactor"] for r in retained):
                result["docking_hold_reasons"].append("retained_nonprotein_components_require_preparation_decision")
            if any(a["element"] in METALS for a in selected_atoms):
                result["docking_hold_reasons"].append("metal_containing_ligand_requires_supported_engine")
            if _covalent(inventory["connections"], selected_atoms, atoms):
                result["docking_hold_reasons"].append("covalent_or_coordination_link_requires_supported_engine")
            if identity["unassigned_stereo"] or identity["unassigned_bond_stereo"]:
                result["docking_hold_reasons"].append("bound_ligand_undefined_stereochemistry")
            result.update(status="verified", reason="bound_structure_and_ligand_verified",
                target_id="structure:" + digest([structure_sha, model, assembly, selected_chains, protein_coords_sha]),
                canonical_smiles=identity["canonical_smiles"], exact_smiles=identity["identity_smiles"], identity=identity,
                query_seed={"smiles": identity["identity_smiles"], "identity": identity,
                    "observed_sdf_sha256": seed_sha,
                    "atom_mapping": [{"template_atom_index": i, "source_atom_name": seed.GetAtomWithIdx(i).GetProp("source_atom_name"),
                        "source_serial": seed.GetAtomWithIdx(i).GetProp("source_serial"), "coordinates": list(seed.GetConformer().GetAtomPosition(i))}
                        for i in range(seed.GetNumAtoms())]},
                protein_structure_sha256=protein_sha, protein_structure_format="pdb" if protein_sha != structure_sha else form,
                original_protein_heavy_coordinates_sha256=protein_coords_sha,
                retained_components_sha256=self.artifacts.put_bytes(canonical_json(retained))["sha256"],
                pocket_box=_box(selected_atoms), compute_readiness="docking_held" if result["docking_hold_reasons"] else "structure_preparation_required",
                provenance={"protein_identity_basis": "selected_structure_coordinates", "ligand_identity_basis": chemistry_receipts["basis"],
                    "activity_basis": "observed_bound_pose_only", "protein_substitution": False,
                    "assembly_context": assembly})
            return result
        except (ValueError, OSError, KeyError, IndexError, RuntimeError) as error:
            import requests
            # RequestException inherits OSError: network failures reach this
            # handler before the generic exception handler below.
            if isinstance(error, requests.RequestException):
                return {**result, "status": "source_unavailable", "reason": "structure_source_unavailable",
                        "error_type": type(error).__name__}
            return {**result, "status": "selection_required", "reason": "structure_identity_unresolved",
                    "error_type": type(error).__name__, "diagnostic": str(error)[:300]}
        except Exception as error:
            import requests
            if not isinstance(error, requests.RequestException):
                raise
            return {**result, "status": "source_unavailable", "reason": "structure_source_unavailable", "error_type": type(error).__name__}

    def _ligand(self, spec, observed, ligand, owner, receipts):
        if "ligand_sdf" in spec:
            text, record = self._uploaded(spec["ligand_sdf"], owner)
            if record["format"] != "sdf":
                raise ValueError("Supplemental ligand must be an uploaded SDF")
            return _map_template(_read_sdf(text), observed), {"basis": "user_sdf", "sha256": spec["ligand_sdf"]}
        if "ligand_smiles" in spec:
            parser = Chem.SmilesParserParams()
            parser.parseName = False
            template = Chem.MolFromSmiles(spec["ligand_smiles"], parser)
            if template is None:
                raise ValueError("Supplemental ligand SMILES is invalid")
            return _map_template(template, observed), {"basis": "user_smiles", "original_smiles": spec["ligand_smiles"]}
        # Uploaded coordinate-only custom components never trigger a remote
        # lookup that might misidentify a user's private LIG/UNK residue code.
        if spec["format"] == "upload":
            return None, {"basis": "missing_exact_ligand_graph"}
        code = ligand["code"]
        if not re.fullmatch(r"[A-Za-z0-9]{1,8}", code):
            raise ValueError("Unsafe chemical component identifier")
        ccd, ccd_sha = self._request("https://files.rcsb.org/ligands/download/" + code + ".cif", {}, receipts)
        sdf, sdf_sha = self._request("https://files.rcsb.org/ligands/download/" + code + "_ideal.sdf", {}, receipts)
        seed = _map_ccd(_cif_data(ccd), _read_sdf(sdf), observed)
        # An instance SDF is a useful independent source receipt, but its atom
        # coordinates can never replace the selected original model/altloc.
        instance_sha = None
        try:
            _, instance_sha = self._request("https://models.rcsb.org/v1/" + spec["value"].lower() + "/ligand",
                {"auth_asym_id": ligand["chain_id"], "auth_seq_id": re.sub(r"[A-Za-z]+$", "", ligand["residue_id"]),
                 "auth_comp_id": code, "model_nums": ligand["model_id"], "encoding": "sdf"}, receipts)
        except Exception as error:
            import requests
            if not isinstance(error, (OSError, ValueError, requests.RequestException)):
                raise
        return seed, {"basis": "official_ccd_atom_names_and_sdf_exact_graph", "ccd_sha256": ccd_sha,
                      "template_sdf_sha256": sdf_sha, "instance_sdf_sha256": instance_sha}


def _read_sdf(text):
    from itertools import islice
    supplier = Chem.ForwardSDMolSupplier(BytesIO(text.encode()), sanitize=True, removeHs=True)
    mols = list(islice(supplier, 2))
    if len(mols) != 1 or mols[0] is None or not 1 <= mols[0].GetNumAtoms() <= 512:
        raise ValueError("Supplemental SDF must contain exactly one valid bounded molecule")
    return mols[0]


def _map_template(template, observed, *, names=None):
    template = Chem.RemoveHs(template)
    if template.GetNumHeavyAtoms() != len(observed) or len(observed) > 512:
        raise ValueError("Exact ligand graph and observed heavy atom counts differ")
    if names is None:
        # Proximity proposes only an atom correspondence; all output elements,
        # bond orders, charges and stereo come from the explicit exact graph.
        block = _pdb([{**a, "record": "HETATM", "code": "LIG", "chain": "A", "residue": "1", "altloc": "", "source_line": None} for a in observed])
        proposed = Chem.MolFromPDBBlock(block, sanitize=False, removeHs=True, proximityBonding=True)
        if proposed is None:
            raise ValueError("Cannot match exact chemistry to observed ligand atoms")
        assigned = AllChem.AssignBondOrdersFromTemplate(template, proposed)
        matches = assigned.GetSubstructMatches(template, uniquify=False, useChirality=False, maxMatches=1001)
        if not matches or len(matches) > 1000:
            raise ValueError("Exact ligand graph has no bounded atom correspondence")
        mapping = min(matches)
    else:
        mapping = names
    output = Chem.Mol(template)
    output.RemoveAllConformers()
    conformer = Chem.Conformer(output.GetNumAtoms())
    conformer.Set3D(True)
    for index, source_index in enumerate(mapping):
        atom, source = output.GetAtomWithIdx(index), observed[source_index]
        if atom.GetSymbol().upper() != source["element"]:
            raise ValueError("Ligand atom mapping element conflict")
        conformer.SetAtomPosition(index, source["xyz"])
        atom.SetProp("source_atom_name", source["name"])
        atom.SetProp("source_serial", source["serial"])
    output.AddConformer(conformer)
    spatial = Chem.Mol(output)
    Chem.AssignStereochemistryFrom3D(spatial, replaceExistingTags=True)
    Chem.AssignStereochemistry(output, cleanIt=True, force=True)
    Chem.AssignStereochemistry(spatial, cleanIt=True, force=True)
    for expected, measured in zip(output.GetAtoms(), spatial.GetAtoms()):
        if expected.HasProp("_CIPCode") and measured.HasProp("_CIPCode") and expected.GetProp("_CIPCode") != measured.GetProp("_CIPCode"):
            raise ValueError("Exact ligand stereochemistry conflicts with the selected observed coordinates")
    stereo = {Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ}
    for expected, measured in zip(output.GetBonds(), spatial.GetBonds()):
        if expected.GetStereo() in stereo and measured.GetStereo() in stereo and expected.GetStereo() != measured.GetStereo():
            raise ValueError("Exact ligand bond stereochemistry conflicts with the selected observed coordinates")
    return output


def _map_ccd(data, template, observed):
    names = _column(data, "_chem_comp_atom.atom_id")
    elements = _column(data, "_chem_comp_atom.type_symbol")
    charges = _column(data, "_chem_comp_atom.charge", "0", len(names))
    if len(names) != len(elements):
        raise ValueError("CCD atom identity table incomplete")
    known = {name: (element.upper(), charge) for name, element, charge in zip(names, elements, charges) if element.upper() not in {"H", "D"}}
    seen = {a["name"]: i for i, a in enumerate(observed)}
    if set(known) != set(seen):
        raise ValueError("CCD ligand is incomplete or observed atom names conflict")
    rw, index = Chem.RWMol(), {}
    for name, (element, charge) in known.items():
        atom = Chem.Atom(element.capitalize())
        atom.SetFormalCharge(int(charge))
        index[name] = rw.AddAtom(atom)
    first = _column(data, "_chem_comp_bond.atom_id_1")
    second = _column(data, "_chem_comp_bond.atom_id_2")
    order = _column(data, "_chem_comp_bond.value_order")
    if not first or len(first) != len(second) or len(first) != len(order):
        raise ValueError("CCD bond orders missing")
    kinds = {"SING": Chem.BondType.SINGLE, "DOUB": Chem.BondType.DOUBLE, "TRIP": Chem.BondType.TRIPLE, "AROM": Chem.BondType.AROMATIC}
    for a, b, kind in zip(first, second, order):
        if a not in index or b not in index:
            continue
        if kind not in kinds:
            raise ValueError("Unsupported CCD bond order")
        rw.AddBond(index[a], index[b], kinds[kind])
    graph = rw.GetMol()
    Chem.SanitizeMol(graph)
    template = Chem.RemoveHs(template)
    matches = graph.GetSubstructMatches(template, uniquify=False, useChirality=False, maxMatches=1001)
    if template.GetNumAtoms() != graph.GetNumAtoms() or not matches or len(matches) > 1000:
        raise ValueError("CCD bond graph and exact stereochemical template differ")
    graph_names = list(known)
    correspondence = [seen[graph_names[i]] for i in min(matches)]
    return _map_template(template, observed, names=correspondence)


def _covalent(links, selected, atoms):
    serials = {a["serial"] for a in selected}
    protein_serials = {a["serial"] for a in atoms if a["protein"]}
    ligand = selected[0]
    for link in links:
        if link["type"] == "CONECT":
            # PDB serials occupy adjacent five-column fields. Splitting digit
            # runs merges consecutive full-width serials (e.g. 1000110005)
            # and can silently miss a ligand-to-protein covalent bond.
            fields = [link["record"][start:start + 5].strip() for start in range(6, 31, 5)]
            if not fields[0] or any(value and not value.isdecimal() for value in fields):
                raise ValueError("Invalid fixed-width PDB CONECT atom serial")
            values = [value for value in fields if value]
            if values and values[0] in serials and any(v in protein_serials for v in values[1:]):
                return True
            if values and values[0] in protein_serials and any(v in serials for v in values[1:]):
                return True
        elif link["type"] == "LINK":
            line = link["record"]
            if ((line[21:22].strip(), line[22:27].strip()) == (ligand["chain"], ligand["residue"])
                    or (line[51:52].strip(), line[52:57].strip()) == (ligand["chain"], ligand["residue"])):
                return True
        elif str(link["type"]).lower().startswith(("covale", "metalc")):
            for n in (1, 2):
                if (str(link.get(f"ptnr{n}_auth_asym_id")), str(link.get(f"ptnr{n}_auth_seq_id"))) == (ligand["chain"], re.sub(r"[A-Za-z]+$", "", ligand["residue"])):
                    return True
    return False


def _box(atoms, padding=6.0):
    low = [min(a["xyz"][i] for a in atoms) for i in range(3)]
    high = [max(a["xyz"][i] for a in atoms) for i in range(3)]
    return {"center": [(a + b) / 2 for a, b in zip(low, high)],
            "size": [max(12.0, b - a + 2 * padding) for a, b in zip(low, high)],
            "padding_angstrom": padding, "source": "selected_observed_ligand", "coordinates_changed": False}
