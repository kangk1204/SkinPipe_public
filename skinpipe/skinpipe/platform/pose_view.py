"""Snapshot-bound MolStar data; PDBQT projections preserve saved coordinates."""
from __future__ import annotations

from collections import Counter
import hashlib
import io
import json
import math
import re

from .artifacts import ArtifactIntegrityError
from .engines import inspect_output
from .identity import candidate_identity

_ELEMENTS = {"A": "C", "C": "C", "N": "N", "NA": "N", "NS": "N", "OA": "O", "O": "O", "OS": "O",
             "S": "S", "SA": "S", "P": "P", "H": "H", "HD": "H", "HS": "H", "F": "F", "CL": "Cl",
             "BR": "Br", "I": "I", "ZN": "Zn", "FE": "Fe", "MG": "Mg", "MN": "Mn", "CA": "Ca",
             "CU": "Cu", "CO": "Co", "NI": "Ni", "CD": "Cd", "HG": "Hg", "SI": "Si",
             "B": "B", "SE": "Se", "CR": "Cr"}


def _hash(raw):
    return hashlib.sha256(raw).hexdigest()


def _read(artifacts, sha256):
    path = artifacts.verified_path(sha256)
    if path.stat().st_size > 16 * 1024 * 1024:
        raise ValueError("Viewer artifact exceeds its 16 MiB limit")
    return artifacts.get_bytes(sha256)


def _models(text):
    if not any(line.startswith("MODEL") for line in text.splitlines()):
        if any(line.startswith("ENDMDL") for line in text.splitlines()):
            raise ValueError("Unexpected model boundary")
        return [text]
    blocks, current, seen = [], None, set()
    for line in text.splitlines():
        if line.startswith("MODEL"):
            match = re.fullmatch(r"MODEL\s+([0-9]+)\s*", line)
            if current is not None or not match or int(match[1]) < 1 or match[1] in seen:
                raise ValueError("Ambiguous model boundaries")
            seen.add(match[1])
            current = []
        elif line.startswith("ENDMDL"):
            if current is None:
                raise ValueError("Unexpected model terminator")
            blocks.append("\n".join(current) + "\n")
            current = None
        elif current is not None:
            current.append(line)
        elif line.startswith(("ATOM  ", "HETATM", "REMARK VINA RESULT:")):
            raise ValueError("Coordinates or scores outside model boundaries")
    if current is not None or not blocks or len(blocks) > 1000:
        raise ValueError("Incomplete or excessive pose models")
    return blocks


def _project(text, *, ligand=False):
    # Same display contract as web/bridge/raw_docking.projection: no coordinate
    # inference or movement; Meeko macrocycle glue is explicitly omitted.
    records, points, omitted = [], [], 0
    for line in text.splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        if len(line) < 78:
            raise ValueError("Truncated PDBQT atom")
        xyz = tuple(float(line[start:start + 8]) for start in (30, 38, 46))
        if not all(math.isfinite(value) for value in xyz):
            raise ValueError("Nonfinite atom coordinate")
        atom_type = line[77:].strip().upper()
        if ligand and re.fullmatch(r"G[0-9]+", atom_type):
            omitted += 1
            continue
        element = "C" if ligand and re.fullmatch(r"CG[0-9]+", atom_type) else _ELEMENTS.get(atom_type)
        if not element:
            raise ValueError("Unsupported PDBQT display atom type")
        record = line[:66].ljust(76) + f"{element:>2}  "
        records.append("HETATM" + record[6:] if ligand else record)
        points.append((element, xyz))
    if not records or ligand and len(records) > 2000:
        raise ValueError("Empty or excessive ligand coordinates")
    data = "\n".join(records) + "\nTER\nEND\n"
    return {"format": "pdb", "data": data, "projection_sha256": _hash(data.encode()),
            "displayed_atoms": len(records), "omitted_glue_points": omitted,
            "coordinates_modified": False, "bond_orders_verified": False if ligand else None}, points


def _mapped_heavy_points(expected, prepared_text, points, pose_text):
    """Bind Meeko's one-based SMILES atom map to saved PDBQT serials.

    An element/coordinate inventory alone cannot distinguish a carbon swap
    that changes connectivity. Do not infer this mapping from bond lengths.
    """
    heavy = {a.GetIdx() for a in expected.GetAtoms() if a.GetAtomicNum() != 1}
    mapping = {}
    serials = set()
    for line in prepared_text.splitlines():
        if not line.startswith("REMARK SMILES IDX "):
            continue
        values = line[len("REMARK SMILES IDX "):].split()
        if len(values) % 2:
            raise ValueError("Invalid prepared atom mapping")
        try:
            pairs = [(int(values[i]) - 1, int(values[i + 1])) for i in range(0, len(values), 2)]
        except ValueError as exc:
            raise ValueError("Invalid prepared atom mapping") from exc
        for atom, serial in pairs:
            if atom < 0 or atom >= expected.GetNumAtoms() or serial < 1 or atom in mapping or serial in serials:
                raise ValueError("Duplicate or invalid prepared atom mapping")
            mapping[atom] = serial
            serials.add(serial)
    # A single heavy atom has exactly one possible correspondence, including
    # old methane fixtures. Multi-atom legacy files need explicit provenance.
    if not mapping and len(heavy) == 1:
        return {next(iter(heavy)): next(p for p in points if p[0] != "H")}
    if not heavy.issubset(mapping):
        raise ValueError("Prepared PDBQT lacks complete heavy atom mapping")
    selected = pose_text if pose_text is not None else prepared_text
    records = [line for line in selected.splitlines() if line.startswith(("ATOM  ", "HETATM"))
               and not re.fullmatch(r"G[0-9]+", line[77:].strip().upper())]
    if len(records) != len(points):
        raise ValueError("Pose atom mapping count differs from saved coordinates")
    by_serial = {}
    for line, point in zip(records, points):
        serial = int(line[6:11])
        if serial in by_serial:
            raise ValueError("Duplicate pose atom serial")
        by_serial[serial] = point
    if not {mapping[i] for i in heavy}.issubset(by_serial):
        raise ValueError("Pose lacks a mapped heavy atom")
    result = {i: by_serial[mapping[i]] for i in heavy}
    if any(expected.GetAtomWithIdx(i).GetSymbol() != point[0] for i, point in result.items()):
        raise ValueError("Mapped pose atom element differs from prepared graph")
    return result


def _checked_sdf(raw, prepared_text, points, *, pose_text=None):
    from rdkit import Chem
    molecules = list(Chem.ForwardSDMolSupplier(io.BytesIO(raw), removeHs=False, sanitize=True))
    if len(molecules) != 1 or molecules[0] is None:
        raise ValueError("Viewer SDF must contain one valid pose")
    molecule = molecules[0]
    if molecule.GetNumAtoms() > 2000 or molecule.GetNumConformers() != 1:
        raise ValueError("SDF pose size or conformer count is invalid")
    smiles = set(re.findall(r"^REMARK SMILES (?!IDX)([^\n]+)", prepared_text, re.M))
    if len(smiles) != 1:
        raise ValueError("Prepared PDBQT lacks unambiguous chemical graph metadata")
    expected = Chem.MolFromSmiles(next(iter(smiles)))
    if expected is None or Chem.MolToSmiles(Chem.RemoveHs(molecule), isomericSmiles=False) != Chem.MolToSmiles(
            Chem.RemoveHs(expected), isomericSmiles=False):
        raise ValueError("SDF bond-order graph differs from prepared chemical structure")
    expected_points = [(e, xyz) for e, xyz in points if e != "H"]
    actual_points = []
    conf = molecule.GetConformer()
    for atom in molecule.GetAtoms():
        position = conf.GetAtomPosition(atom.GetIdx())
        xyz = (position.x, position.y, position.z)
        if not all(math.isfinite(value) for value in xyz):
            raise ValueError("SDF has nonfinite coordinates")
        if atom.GetAtomicNum() != 1:
            actual_points.append((atom.GetSymbol(), xyz))
    if Counter(e for e, _ in actual_points) != Counter(e for e, _ in expected_points):
        raise ValueError("SDF heavy atom inventory differs from selected raw pose")
    mapped_points = _mapped_heavy_points(expected, prepared_text, points, pose_text)
    atom_map = {}
    for atom in molecule.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue
        xyz = tuple(conf.GetAtomPosition(atom.GetIdx()))
        matches = [i for i, (element, position) in mapped_points.items()
                   if element == atom.GetSymbol()
                   and all(abs(a - b) <= .002 for a, b in zip(xyz, position))]
        if not matches:
            raise ValueError("SDF coordinates differ from selected raw pose")
        if len(matches) != 1 or matches[0] in atom_map.values():
            raise ValueError("SDF atom-to-pose correspondence is ambiguous")
        atom_map[atom.GetIdx()] = matches[0]
        original = expected.GetAtomWithIdx(matches[0])
        if (atom.GetIsotope(), atom.GetFormalCharge()) != (original.GetIsotope(), original.GetFormalCharge()):
            raise ValueError("SDF mapped atom identity differs from prepared graph")
    if set(atom_map.values()) != set(mapped_points):
        raise ValueError("SDF heavy atom correspondence is incomplete")
    def edges(mol, remap):
        return {(*sorted((remap[b.GetBeginAtomIdx()], remap[b.GetEndAtomIdx()])),
                 str(b.GetBondType()), b.GetIsAromatic())
                for b in mol.GetBonds()
                if b.GetBeginAtomIdx() in remap and b.GetEndAtomIdx() in remap}
    if edges(molecule, atom_map) != edges(expected, {i: i for i in mapped_points}):
        raise ValueError("SDF atom-to-pose correspondence changes chemical connectivity")
    return {"format": "sdf", "data": raw.decode(), "coordinates_modified": False,
            "bond_orders_verified": True, "atom_correspondence_verified": True,
            "stereo_assignment_validation": "not_assessed_by_display_adapter"}


def _reference(artifacts, card_hash, receptor_hash):
    """Resolve only the native card's frozen co-crystal in its saved frame."""
    from rdkit import Chem

    raw_card = _read(artifacts, card_hash)
    if len(raw_card) > 1024 * 1024:
        raise ValueError("Receptor identity card exceeds viewer manifest limit")
    card = json.loads(raw_card)
    entries = card["artifacts"]
    if entries["pdbqt"].get("sha256") != receptor_hash or entries["pdbqt"].get("present") is not True:
        raise ArtifactIntegrityError("Reference receptor card differs from observed receptor")
    entry = entries.get("cocrystal_sdf")
    if not entry or entry.get("present") is not True:
        return None
    if entry.get("path") != card["receptor"].get("cocrystal_sdf"):
        raise ArtifactIntegrityError("Reference co-crystal differs from its receptor manifest")
    raw = _read(artifacts, entry["sha256"])
    molecules = list(Chem.ForwardSDMolSupplier(io.BytesIO(raw), removeHs=False, sanitize=True))
    if len(molecules) != 1 or molecules[0] is None:
        raise ValueError("Reference SDF must contain one valid molecule")
    molecule = molecules[0]
    if not 0 < molecule.GetNumAtoms() <= 2000 or molecule.GetNumConformers() != 1:
        raise ValueError("Reference SDF has invalid atom/conformer count")
    conformer = molecule.GetConformer()
    for atom in molecule.GetAtoms():
        point = conformer.GetAtomPosition(atom.GetIdx())
        if not all(math.isfinite(value) for value in (point.x, point.y, point.z)):
            raise ValueError("Reference SDF has nonfinite coordinates")
    return {"format": "sdf", "data": raw.decode(), "filename": "reference.sdf",
            "source_sha256": entry["sha256"], "identity_card_sha256": card_hash,
            "coordinates_modified": False, "bond_orders_verified": False,
            "validation_scope": "Hash-bound native co-crystal and finite saved coordinates; bond assignments not independently revalidated",
            "completeness": card["receptor"].get("cocrystal_completeness"),
            "is_fragment": card["receptor"].get("cocrystal_is_fragment")}


def candidate_docking_view(catalog, artifacts, candidate_id, *, snapshot_id, target_id, observation_index=0):
    if (not isinstance(target_id, str) or not target_id or isinstance(observation_index, bool)
            or not isinstance(observation_index, int) or not 0 <= observation_index <= 100):
        raise ValueError("Target context and bounded integer observation index required")
    candidate = catalog.get_candidate(snapshot_id, candidate_id, target_id=target_id)
    if candidate is None:
        raise KeyError("Candidate not found in this target snapshot")
    rows = candidate.get("docking") or []
    rows = [rows] if isinstance(rows, dict) else rows
    if not isinstance(rows, list) or observation_index >= len(rows):
        raise KeyError("Docking observation not found")
    row = rows[observation_index]
    if row.get("artifact_valid") is not True or row.get("status") not in {"complete", "review_required"}:
        raise ValueError("This docking observation has no admitted pose artifact")
    identity = candidate_identity(candidate, trust_stored=True)["exact_structure_id"]
    if (not identity or row.get("exact_structure_id") != identity or row.get("target_id") != target_id
            or row.get("candidate_id") != candidate_id):
        raise ArtifactIntegrityError("Docking observation differs from exact candidate identity")
    protocol = row["protocol"]
    receptor_hash = protocol["receptor_sha256"]
    if row.get("receptor_sha256", receptor_hash) != receptor_hash:
        raise ArtifactIntegrityError("Docking receptor identity mismatch")
    receptor_raw = _read(artifacts, receptor_hash)
    prepared = _read(artifacts, row["input_sha256"]).decode()
    pose_raw = _read(artifacts, row["pose_sha256"])
    checked = inspect_output(pose_raw.decode(), prepared)
    display_score = row.get("score") if row.get("score") is not None else row.get("affinity")
    if display_score != checked["affinity"] or row.get("status") != checked["status"]:
        raise ArtifactIntegrityError("Docking score/status differs from raw pose")
    models = _models(pose_raw.decode())
    if len(models) != checked["n_poses"]:
        raise ArtifactIntegrityError("Raw pose model count mismatch")
    pose_index = min(range(len(models)), key=lambda index: (checked["pose_scores"][index], index))
    receptor, _ = _project(receptor_raw.decode())
    ligand, points = _project(models[pose_index], ligand=True)
    warnings = []
    if row.get("sdf_artifact"):
        sdf_hash = row["sdf_artifact"]["sha256"]
        try:
            ligand = _checked_sdf(_read(artifacts, sdf_hash), prepared, points, pose_text=models[pose_index])
            ligand["source_sha256"] = sdf_hash
        except (ArtifactIntegrityError, OSError, ValueError, RuntimeError) as error:
            warnings.append({"component": "sdf", "status": "quarantined", "reason": str(error)})
    receptor.update(source_sha256=receptor_hash, filename="receptor.pdb")
    ligand.setdefault("source_sha256", row["pose_sha256"])
    ligand.update(filename="ligand." + ligand["format"], raw_pose_sha256=row["pose_sha256"],
                  selected_pose_index=pose_index, selected_pose_score=checked["affinity"])
    reference = None
    if row.get("receptor_identity_card_sha256"):
        try:
            reference = _reference(artifacts, row["receptor_identity_card_sha256"], receptor_hash)
        except (ArtifactIntegrityError, OSError, ValueError, RuntimeError, KeyError, TypeError) as error:
            warnings.append({"component": "reference", "status": "quarantined", "reason": str(error)})
    return {"candidate_id": candidate_id, "target_id": target_id, "snapshot_id": snapshot_id,
            "observation_index": observation_index, "receptor": receptor, "ligand": ligand,
            "original_pose": {"format": "pdbqt", "data": pose_raw.decode(), "source_sha256": row["pose_sha256"],
                              "filename": "docking-poses.pdbqt", "coordinates_modified": False},
            "reference": reference, "protocol": protocol, "warnings": warnings,
            "scope": "Saved coordinates and exact source binding; pose/affinity do not establish biological binding",
            "scientific_validation": "NOT_ESTABLISHED"}
