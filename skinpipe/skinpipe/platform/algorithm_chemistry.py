"""Optional, bounded chemistry transforms over registered molecular artifacts.

Admission, process ownership and aggregate resource reservations belong to the
worker. This module runs one thread per native operation, preserves every input
row, and never changes discovery weights or establishes biological performance.
Cancellation is cooperative between operations; CDPKit also receives its native
abort callback. A worker must enforce the enclosing process/lease timeout.
"""
from __future__ import annotations

from copy import deepcopy
import hashlib
import importlib
from importlib import metadata
import json
import math
from pathlib import Path
import resource
import sys
import time
import uuid

from .identity import chemical_identity, digest

CAPS = {"dimorphite": 32, "cdpl": 2000, "posebusters": 1000,
        "prolif": 192, "fpsim2": 10000}
DEPENDENCIES = {"dimorphite": ("dimorphite-dl",), "cdpl": ("CDPKit",),
                "posebusters": ("posebusters",), "prolif": ("prolif", "MDAnalysis"),
                "fpsim2": ("FPSim2",)}
CONTROL_ROLES = {"control", "positive_control", "negative_control", "counter_control", "reference_control"}
MAX_FILE_BYTES = 64 * 1024 * 1024
MAX_LIGAND_ATOMS = 256


class ChemistryHeld(ValueError):
    """An unavailable or inadmissible artifact holds only its candidate branch."""


class ChemistryStopped(RuntimeError):
    def __init__(self, status):
        self.status = status
        super().__init__(status)


def _load(name):
    return importlib.import_module(name)


def _sha(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def _software(provider):
    result = {"python": sys.version.split()[0], "packages": {}, "module_sha256": {}}
    for package in ("rdkit", *DEPENDENCIES[provider]):
        try:
            dist = metadata.distribution(package)
            result["packages"][package] = {"version": dist.version,
                "declared_license": dist.metadata.get("License-Expression") or dist.metadata.get("License")}
        except metadata.PackageNotFoundError:
            result["packages"][package] = {"version": None, "status": "missing"}
    return result


def _module(software, name):
    mod = _load(name)
    path = getattr(mod, "__file__", None)
    if path and Path(path).is_file():
        software["module_sha256"][name] = {"path": str(Path(path).resolve()), "sha256": _sha(path)}
    return mod


def _clean(value):
    """JSON finite scalars, preserving missing values rather than inventing zero."""
    if type(value).__name__ in {"NAType", "NaTType"}:
        return None
    if isinstance(value, dict):
        return {str(key): _clean(val) for key, val in value.items()}
    if isinstance(value, (tuple, list)):
        return [_clean(val) for val in value]
    if hasattr(value, "tolist") and getattr(value, "ndim", 0) > 0:
        return _clean(value.tolist())
    if hasattr(value, "item"):
        return _clean(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if value is None or isinstance(value, (str, bool, int, float)):
        return value
    return str(value)


def _artifact(path):
    path = Path(path).resolve()
    return {"path": str(path), "sha256": _sha(path), "size_bytes": path.stat().st_size}


def _write_json(path, payload):
    raw = json.dumps(_clean(payload), sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_bytes(raw)
    temporary.replace(path)
    return _artifact(path)


def _path(value, expected_sha256=None):
    if not isinstance(value, (str, Path)):
        raise ChemistryHeld("registered_artifact_path_required")
    path = Path(value).resolve(strict=True)
    if not path.is_file() or path.stat().st_size > MAX_FILE_BYTES:
        raise ChemistryHeld("artifact_missing_or_exceeds_64MiB")
    observed = _sha(path)
    if expected_sha256 is not None and observed != expected_sha256:
        raise ChemistryHeld("artifact_sha256_mismatch")
    return path, observed


def _count(context, key, default, cap):
    value = context.get(key, context.get("counts", {}).get(key, default))
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(key + " must be a positive integer")
    return min(value, cap)


def _check(cancel_event, deadline):
    if cancel_event is not None and cancel_event.is_set():
        raise ChemistryStopped("cancelled")
    if time.monotonic() >= deadline:
        raise ChemistryStopped("budget_deferred")


def _smiles_mol(row, Chem):
    value = row.get("smiles")
    if not isinstance(value, str) or not value or len(value) > 8192:
        raise ChemistryHeld("valid_registered_smiles_required")
    mol = Chem.MolFromSmiles(value)
    if mol is None or not mol.GetNumAtoms():
        raise ChemistryHeld("invalid_structure")
    if mol.GetNumAtoms() > MAX_LIGAND_ATOMS:
        raise ChemistryHeld("chemistry_ligand_atom_limit_256")
    return mol


def _exact(mol, Chem):
    return chemical_identity(Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=True))["exact_structure_id"]


def _mapped(mol, Chem):
    mol = Chem.Mol(mol)
    # Internal mapping refers to the input atom order, not untrusted SMILES labels.
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx() + 1)
    return mol


def _pose(row, Chem):
    path, file_hash = _path(row.get("pose_path"), row.get("pose_sha256"))
    if path.suffix.lower() in {".sdf", ".sd"}:
        supplier = Chem.SDMolSupplier(str(path), removeHs=False)
        mol = next(iter(supplier), None)  # One representative pose, never all poses.
    elif path.suffix.lower() == ".mol":
        mol = Chem.MolFromMolFile(str(path), removeHs=False)
    else:
        raise ChemistryHeld("representative_pose_requires_registered_SDF_or_MOL")
    if mol is None or not mol.GetNumConformers():
        raise ChemistryHeld("invalid_pose")
    if mol.GetNumHeavyAtoms() > MAX_LIGAND_ATOMS:
        raise ChemistryHeld("chemistry_ligand_atom_limit_256")
    if not all(math.isfinite(v) for xyz in mol.GetConformer().GetPositions() for v in xyz):
        raise ChemistryHeld("nonfinite_pose_coordinates")
    observed = _exact(mol, Chem)
    expected = row.get("prepared_exact_structure_id") or chemical_identity(row.get("smiles"))["exact_structure_id"]
    if observed != expected:
        raise ChemistryHeld("pose_exact_structure_identity_mismatch")
    return mol, {"pose_sha256": file_hash, "prepared_exact_structure_id": observed,
                 "pose_path": str(path), "representative_pose_index": 0}


def _dimorphite(row, context, directory, software, cancel_event, deadline):
    Chem = _module(software, "rdkit.Chem")
    dm = _module(software, "dimorphite_dl")
    std = _module(software, "rdkit.Chem.MolStandardize.rdMolStandardize")
    original = _smiles_mol(row, Chem)
    mapped = _mapped(original, Chem)
    limit = _count(context, "total_states", 2, 2)
    ph = context.get("pH", 7.4)
    ph_min, ph_max = context.get("ph_min", ph), context.get("ph_max", ph)
    if any(isinstance(v, bool) or not isinstance(v, (float, int)) or not math.isfinite(v)
           for v in (ph_min, ph_max)) or not 0 <= ph_min <= ph_max <= 14:
        raise ValueError("pH range must be finite and within 0..14")
    _check(cancel_event, deadline)
    variants = dm.protonate_smiles(Chem.MolToSmiles(mapped), ph_min=ph_min, ph_max=ph_max,
                                  max_variants=limit, precision=1.0, validate_output=True)
    enumerator = std.TautomerEnumerator()
    enumerator.SetMaxTautomers(32)
    enumerator.SetMaxTransforms(64)
    enumerator.SetRemoveSp3Stereo(False)
    states = {}
    expected_maps = {a.GetAtomMapNum(): (a.GetAtomicNum(), a.GetIsotope()) for a in mapped.GetAtoms()}
    for variant in variants:
        _check(cancel_event, deadline)
        protomer = Chem.MolFromSmiles(variant)
        if protomer is None:
            continue
        for mol in enumerator.Enumerate(protomer):
            observed_maps = {a.GetAtomMapNum(): (a.GetAtomicNum(), a.GetIsotope()) for a in mol.GetAtoms()}
            if observed_maps != expected_maps or mol.GetNumAtoms() != len(expected_maps):
                raise ChemistryHeld("enumerated_state_lost_input_atom_mapping")
            exact = _exact(mol, Chem)
            score = float(enumerator.ScoreTautomer(mol))
            smi = Chem.MolToSmiles(mol, isomericSmiles=True)
            states.setdefault(exact, (score, smi, mol))
    if not states:
        raise ChemistryHeld("no_enumerated_states")
    preparation = {"method": "Dimorphite-DL+RDKit-tautomer", "ph_min": ph_min, "ph_max": ph_max,
                   "total_state_limit": limit, "tautomer_enumeration_limit": 32,
                   "software_versions": software["packages"]}
    chosen = []
    for exact, (score, smi, mol) in sorted(states.items(), key=lambda item: (-item[1][0], item[0]))[:limit]:
        identity = chemical_identity(row["smiles"], prepared_smiles=smi, preparation=preparation)
        chosen.append({"smiles": smi, "exact_structure_id": exact,
                       "prepared_structure_id": identity["prepared_structure_id"], "tautomer_score": score,
                       "formal_charge": Chem.GetFormalCharge(mol),
                       "atom_mapping": [{"prepared_atom_index": a.GetIdx(),
                           "input_atom_index": a.GetAtomMapNum() - 1, "element": a.GetSymbol()}
                                        for a in mol.GetAtoms()]})
    metrics = {"method": preparation["method"], "states": chosen, "state_count": len(chosen),
               "enumerated_unique_states": len(states), "state_selection": "RDKit tautomer score then exact ID",
               "state_population_probabilities": None, "preparation": preparation}
    return metrics, [_write_json(directory / "states.json", {"candidate_id": row["candidate_id"], **metrics})]


def _cdpl_modules(software):
    result = tuple(_module(software, "CDPL." + name) for name in ("Chem", "ConfGen", "Pharm", "Shape", "Math"))
    # Bind native implementations as well as their Python loaders.
    for name in ("CDPL.Chem._chem", "CDPL.ConfGen._confgen", "CDPL.Pharm._pharm", "CDPL.Shape._shape"):
        _module(software, name)
    return result


def _cdpl_mol(row, directory, modules, software, context, cancel_event, deadline, *, reference=False):
    C, G, P, _, _ = modules
    Chem = _module(software, "rdkit.Chem")
    source_mol = _smiles_mol(row, Chem)
    mol = C.BasicMolecule()
    inputs = {"exact_structure_id": chemical_identity(row["smiles"])["exact_structure_id"],
              "atom_index_scope": "registered SMILES RDKit atom order"}
    pose_path = row.get("reference_path") if reference else row.get("pose_path")
    if pose_path:
        pose_row = {**row, "pose_path": pose_path,
                    "pose_sha256": row.get("reference_sha256") if reference else row.get("pose_sha256")}
        source_mol, pose_receipt = _pose(pose_row, Chem)
        inputs.update(pose_receipt)
        inputs.update(exact_structure_id=pose_receipt["prepared_exact_structure_id"],
                      atom_index_scope="registered representative pose RDKit atom order")
        mapped_path = directory / ("reference-input.sdf" if reference else "candidate-input.sdf")
        writer = Chem.SDWriter(str(mapped_path))
        writer.write(_mapped(source_mol, Chem))
        writer.close()
        reader = C.MoleculeReader(str(mapped_path))
        C.setMultiConfImportParameter(reader, False)
        if not reader.read(mol):
            raise ChemistryHeld("CDPKit_input_reader_failed")
        reader.close()
        C.calcBasicProperties(mol, False)
        num_conformers, generation_method = 1, "registered_representative_pose"
    else:
        if not C.parseSMILES(Chem.MolToSmiles(_mapped(source_mol, Chem)), mol):
            raise ChemistryHeld("CDPKit_SMILES_parse_failed")
        G.prepareForConformerGeneration(mol)
        generator = G.ConformerGenerator()
        count = 1 if reference else _count(context, "conformers", 10, 10)
        generator.settings.clearMaxNumOutputConformersRanges()
        generator.settings.maxNumOutputConformers = count
        generator.settings.maxNumSampledConformers = 64
        generator.settings.maxPoolSize = 64
        generator.settings.maxNumRefinementIterations = 200
        generator.settings.timeout = min(30000, max(1, int((deadline - time.monotonic()) * 1000)))
        generator.abortCallback = lambda: bool(cancel_event is not None and cancel_event.is_set()) or time.monotonic() >= deadline
        code = generator.generate(mol)
        _check(cancel_event, deadline)
        if code not in (G.ReturnCode.SUCCESS, G.ReturnCode.TOO_MUCH_SYMMETRY) or not generator.numConformers:
            raise ChemistryHeld("CONFORGE_generation_failed:" + str(code))
        generator.setConformers(mol)
        num_conformers, generation_method = min(generator.numConformers, count), "CDPKit-CONFORGE"
        for atom in mol.atoms:
            C.set3DCoordinates(atom, C.getConformer3DCoordinates(atom, 0))
    P.prepareForPharmacophoreGeneration(mol)
    inputs.update(conformer_count=num_conformers, conformer_method=generation_method,
                  seed=context.get("seed", 42),
                  seed_status="CONFORGE_default_internal_seed; requested seed not exposed by Python API")
    return mol, inputs


def _features(mol, P, C):
    ph4 = P.BasicPharmacophore()
    P.DefaultPharmacophoreGenerator().generate(mol, ph4)
    features = []
    for index, feature in enumerate(ph4):
        atoms = []
        if P.hasSubstructure(feature):
            for atom in P.getSubstructure(feature).atoms:
                map_id = C.getAtomMappingID(atom) if C.hasAtomMappingID(atom) else 0
                if map_id:
                    atoms.append(map_id - 1)
        features.append({"feature_index": index, "type": int(P.getType(feature)),
                         "position_angstrom": list(C.get3DCoordinates(feature).toArray()),
                         "input_atom_indices": sorted(atoms)})
    return ph4, features


def _references(row, context):
    if row.get("reference_path"):
        value = row.get("reference_smiles") or context.get("reference_smiles")
        if isinstance(value, list):
            value = value[0] if value else None
        return [{"smiles": value, "reference_path": row["reference_path"],
                 "reference_sha256": row.get("reference_sha256")}] if value else []
    values = context.get("reference_smiles", [])
    if isinstance(values, str):
        values = [values]
    if not isinstance(values, list):
        raise ValueError("reference_smiles must be a string or list")
    refs = {}
    for value in values:
        item = value if isinstance(value, dict) else {"smiles": value}
        identity = chemical_identity(item.get("smiles"))
        if identity["structure_status"] == "valid":
            refs.setdefault(identity["exact_structure_id"], item)
    return [refs[key] for key in sorted(refs)[:_count(context, "reference_limit", 5, 5)]]


def _cdpl(row, context, directory, software, cancel_event, deadline):
    refs = _references(row, context)
    if not refs:
        raise ChemistryHeld("valid_alignment_reference_required")
    modules = _cdpl_modules(software)
    C, _, P, S, M = modules
    mol, inputs = _cdpl_mol(row, directory, modules, software, context, cancel_event, deadline)
    gen = S.GaussianShapeGenerator()
    gen.generate(mol)
    scores, best, best_reference = [], None, None
    for ref_index, reference in enumerate(refs):
        _check(cancel_event, deadline)
        ref_dir = directory / ("reference-" + str(ref_index))
        ref_dir.mkdir()
        ref_mol, ref_inputs = _cdpl_mol(reference, ref_dir, modules, software, context,
                                      cancel_event, deadline, reference=True)
        ref_gen = S.GaussianShapeGenerator()
        ref_gen.generate(ref_mol)
        alignment = S.GaussianShapeAlignment()
        alignment.scoringFunction = S.TanimotoComboScore()
        alignment.maxNumOptIterations = 100
        alignment.addReferenceShapes(ref_gen.shapes)
        if not alignment.align(gen.shapes) or not alignment.numResults:
            scores.append({"reference": ref_inputs, "status": "held", "reason": "no_shape_alignment"})
            continue
        result = alignment.getResult(0)
        record = {"reference": ref_inputs, "status": "complete", "shape_tanimoto": float(S.calcShapeTanimotoScore(result)),
                  "pharmacophore_color_tanimoto": float(S.calcColorTanimotoScore(result)),
                  "tanimoto_combo": float(result.score), "candidate_conformer_index": int(result.alignedShapeIndex),
                  "transform": result.transform.toArray().tolist()}
        scores.append(record)
        if best is None or record["tanimoto_combo"] > best["tanimoto_combo"]:
            best, best_reference = record, C.BasicMolecule(ref_mol)
    if best is None:
        raise ChemistryHeld("no_shape_alignment")
    conformer = best["candidate_conformer_index"]
    if inputs["conformer_method"] == "CDPKit-CONFORGE":
        for atom in mol.atoms:
            C.set3DCoordinates(atom, C.getConformer3DCoordinates(atom, conformer))
    candidate_ph4, candidate_features = _features(mol, P, C)
    reference_ph4, reference_features = _features(best_reference, P, C)
    pharm = P.PharmacophoreAlignment(True)
    pharm.performExhaustiveSearch(False)
    pharm.setMinTopologicalMappingSize(min(3, len(candidate_features), len(reference_features)))
    pharm.addFeatures(reference_ph4, True)
    pharm.addFeatures(candidate_ph4, False)
    fit = P.PharmacophoreFitScore(True)
    feature_best, solutions = None, 0
    if candidate_features and reference_features:
        while solutions < 64:
            _check(cancel_event, deadline)
            if not pharm.nextAlignment():
                break
            solutions += 1
            score = float(fit(reference_ph4, candidate_ph4, pharm.transform))
            pairs = [{"reference_feature_index": int(pair.first), "candidate_feature_index": int(pair.second)}
                     for pair in pharm.topMapping]
            if feature_best is None or score > feature_best["fit_score"]:
                feature_best = {"fit_score": score, "feature_pairs": pairs,
                                "transform": pharm.transform.toArray().tolist()}
    # Save the molecular shape alignment. The separately reported pharmacophore
    # alignment may use a different transform; both are explicit, independent of MCS.
    coords = M.Vector3DArray()
    C.get3DCoordinates(mol, coords)
    transform = M.Matrix4D()
    transform.assign(_load("numpy").array(best["transform"]))
    M.transform(coords, transform)
    C.clearConformations(mol)
    C.set3DCoordinates(mol, coords)
    output = directory / "shape-aligned.sdf"
    writer = C.MolecularGraphWriter(str(output))
    if not writer.write(mol):
        raise ChemistryHeld("CDPKit_alignment_writer_failed")
    writer.close()
    Chem = _load("rdkit.Chem")
    reloaded = next(iter(Chem.SDMolSupplier(str(output), removeHs=False)), None)
    if reloaded is None or _exact(reloaded, Chem) != inputs["exact_structure_id"]:
        raise ChemistryHeld("aligned_output_exact_structure_identity_mismatch")
    metrics = {"method": "CDPKit-CONFORGE+GaussianShapeAlignment+PharmacophoreAlignment",
               "inputs": inputs, "references": scores, "best_shape_alignment": best,
               "pharmacophore_alignment": feature_best, "pharmacophore_solutions_examined": solutions,
               "pharmacophore_solution_limit": 64, "candidate_features": candidate_features,
               "reference_features": reference_features, "MCS_required": False,
               "coordinate_unit": "angstrom", "score_unit": "dimensionless"}
    return metrics, [_artifact(output), _write_json(directory / "correspondence.json", metrics)]


def _posebusters_binary_summary(config, report):
    """Read PoseBusters' normalized public columns without converting missing to failure."""
    selected = {}
    for module in config["modules"]:
        for name in module.get("chosen_binary_test_output", []):
            label = (module.get("rename_outputs") or {}).get(
                name, name + (module.get("rename_suffix") or ""))
            # PoseBusters normalizes both full and binary reports in _collect_in_table.
            column = label.lower().replace(" ", "_")
            value = report.get(column)
            selected[label] = value if type(value) is bool else None
    missing = [name for name, value in selected.items() if value is None]
    passed = None if not selected or missing else all(selected.values())
    return selected, missing, passed


class _PoseBustersResources:
    """One byte-verified conditioning molecule per chemistry request.

    Reuse parsing only, never tests, ligand energy ensembles or QC verdicts.
    Each call still hashes its registered receptor. A copy prevents a check
    from mutating the molecule seen by subsequent candidates.
    """

    def __init__(self):
        self.key, self.molecule = None, None

    def receptor(self, Chem, path, sha256, config, software):
        # This public Mol input path was checked against the pinned loader.
        # Other versions retain their native path loading until qualified.
        if software["packages"]["posebusters"].get("version") != "0.6.5":
            return path
        params = config.get("loading", {}).get("mol_cond", {})
        key = (str(path), sha256, digest(params), digest(software["module_sha256"]))
        if self.key != key:
            from posebusters.tools.loading import safe_load_mol

            molecule = safe_load_mol(path=path, **params)
            if molecule is None:
                return path
            self.key, self.molecule = key, molecule
        return Chem.Mol(self.molecule)


def _posebusters(row, context, directory, software, cancel_event, deadline, *, _resources=None):
    Chem = _module(software, "rdkit.Chem")
    pb = _module(software, "posebusters")
    mol, inputs = _pose(row, Chem)
    receptor, receptor_hash = _path(context.get("receptor_path"), context.get("receptor_sha256"))
    config = deepcopy(pb.PoseBusters(config="dock", top_n=1, max_workers=0).config)
    for module in config["modules"]:
        if module["function"] == "energy_ratio":
            module.setdefault("parameters", {}).update(ensemble_number_conformations=10, num_threads=1)
    checker = pb.PoseBusters(config=config, top_n=1, max_workers=0)
    _check(cancel_event, deadline)
    conditioned = (_resources.receptor(Chem, receptor, receptor_hash, config, software)
                   if _resources is not None else receptor)
    table = checker.bust(mol_pred=mol, mol_cond=conditioned, full_report=True)
    report = _clean(table.iloc[0].to_dict()) if len(table) else None
    if report is None:
        raise ChemistryHeld("PoseBusters_returned_no_report")
    # Select the same public binary columns from the single full report. Running
    # bust twice would repeat the expensive energy ensemble and alter its sample.
    selected, missing, passed = _posebusters_binary_summary(config, report)
    metrics = {"method": "PoseBusters-dock", "inputs": {**inputs, "receptor_sha256": receptor_hash},
               "checks": selected, "full_report": report,
               "pose_checks_pass": passed, "missing_binary_checks": missing,
               "pose_check_status": "unassessed" if passed is None else "pass" if passed else "failed",
               "configuration_sha256": digest(config), "energy_ensemble_conformers": 10,
               "biological_validation": "NOT_ESTABLISHED"}
    return metrics, [_write_json(directory / "posebusters.json", metrics)]


def _prolif_selection_guard(universe, receptor, context):
    """Compare selected atoms with source-bound observed polymer population.

    Native card residues and amino-acid MODRES declarations identify modified
    polymer, regardless of HETATM spelling. Undeclared ATOM residues with a
    peptide backbone are also protected. Loose ligand/ion HETATM records do not
    become polymer solely because of a residue-name or coordinate resemblance.
    This guard never changes selections or infers missing chemistry.
    """
    amino = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
             "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}
    native = set()
    card_hash, card_scope = None, "not_supplied"
    if context.get("receptor_identity_path"):
        expected_hash = context.get("receptor_identity_sha256")
        if not isinstance(expected_hash, str) or len(expected_hash) != 64:
            raise ChemistryHeld("ProLIF_bound_receptor_identity_checksum_required")
        path, card_hash = _path(context["receptor_identity_path"], expected_hash)
        if path.stat().st_size > 1024 * 1024:
            raise ChemistryHeld("ProLIF_receptor_identity_card_exceeds_1MiB")
        try:
            card = json.loads(path.read_text())
        except ValueError as error:
            raise ChemistryHeld("ProLIF_receptor_identity_card_JSON_invalid") from error
        if not isinstance(card, dict) or not isinstance(card.get("identity"), dict):
            raise ChemistryHeld("ProLIF_receptor_identity_card_schema_invalid")
        model = card["identity"].get("native_modified_repair_model")
        if card.get("schema_version") == 2:
            artifacts = card.get("artifacts")
            artifact = artifacts.get("clean_pdb", {}) if isinstance(artifacts, dict) else {}
            if not isinstance(artifact, dict):
                raise ChemistryHeld("ProLIF_identity_card_artifact_schema_invalid")
            if artifact.get("present") is not True or artifact.get("sha256") != _sha(receptor):
                raise ChemistryHeld("ProLIF_identity_card_receptor_checksum_mismatch")
            card_scope = "schema2_exact_prepared_receptor"
        elif model:
            raise ChemistryHeld("ProLIF_native_identity_card_requires_bound_schema2")
        else:
            card_scope = "legacy_card_without_native_assertions"
        if model:
            records = model.get("native_residues") if isinstance(model, dict) else None
            if not isinstance(records, list) or len(records) > 1000:
                raise ChemistryHeld("ProLIF_native_identity_residue_manifest_invalid")
            for record in records:
                if (not isinstance(record, dict) or any(not isinstance(record.get(key), str)
                                                       for key in ("chain", "residue", "name"))):
                    raise ChemistryHeld("ProLIF_native_identity_residue_entry_invalid")
                native.add((record["chain"], record["residue"], record["name"]))
    source_atoms = {}
    for line in Path(receptor).read_text().splitlines():
        if line.startswith("ENDMDL"):
            break
        if line.startswith("MODRES") and line[24:27].strip() in amino:
            native.add((line[16].strip(), line[18:23].strip(), line[12:15].strip()))
        if line.startswith("ATOM  "):
            key = (line[21].strip(), line[22:27].strip(), line[17:20].strip())
            source_atoms.setdefault(key, set()).add(line[12:16].strip())
    native.update(key for key, names in source_atoms.items() if {"N", "CA", "C"} <= names)
    by_residue = {}
    for atom in universe.atoms:
        key = (atom.chainID.strip(), str(atom.resid) + atom.icode.strip(), atom.resname)
        by_residue.setdefault(key, set()).add(int(atom.index))
    selected = set(map(int, universe.select_atoms("protein").indices))
    expected = set(selected)
    for key in native:
        expected.update(by_residue.get(key, set()))
    excluded, absent = expected-selected, native-by_residue.keys()
    excluded_residues = sorted(key for key in native if by_residue.get(key, set()) & excluded)
    names = [":".join(key) for key in excluded_residues]
    report = {"status": "held" if excluded or absent else "complete",
              "expected_polymer_atom_count": len(expected), "selected_atom_count": len(selected),
              "excluded_polymer_atom_count": len(excluded), "excluded_residues": names,
              "absent_identity_residues": [":".join(key) for key in sorted(absent)],
              "receptor_identity_sha256": card_hash, "identity_scope": card_scope,
              "predicate": "bound-native-card+MODRES-amino-parent+ATOM-peptide-backbone"}
    if excluded or absent:
        report["reason"] = ("ProLIF_protein_selection_excludes_observed_native_polymer:"
                            f"expected_atoms={len(expected)},selected_atoms={len(selected)},excluded_atoms={len(excluded)},"
                            f"absent_residues={len(absent)},excluded_residues={','.join(names[:10])}")
    return report


def _prolif(row, context, directory, software, cancel_event, deadline):
    Chem = _module(software, "rdkit.Chem")
    plf = _module(software, "prolif")
    mda = _module(software, "MDAnalysis")
    mol, inputs = _pose(row, Chem)
    receptor, receptor_hash = _path(context.get("receptor_path"), context.get("receptor_sha256"))
    _check(cancel_event, deadline)
    universe = mda.Universe(str(receptor))
    if universe.atoms.n_atoms > 100000:
        raise ChemistryHeld("ProLIF_receptor_atom_limit_100000")
    selection = _prolif_selection_guard(universe, receptor, context)
    if selection["status"] != "complete":
        raise ChemistryHeld(selection["reason"])
    try:
        # Explicit receptor hydrogens are needed for bond order/charge inference.
        # Do not silently convert a hydrogen-free receptor to an arbitrary model.
        protein = plf.Molecule.from_mda(universe, selection="protein", implicit_hydrogens=False)
    except Exception as error:
        raise ChemistryHeld("ProLIF_prepared_protein_required:" + str(error)[:300]) from error
    ligand = plf.Molecule.from_rdkit(Chem.AddHs(mol, addCoords=True))
    _check(cancel_event, deadline)
    fingerprint = plf.Fingerprint(count=True)
    ifp = fingerprint.generate(ligand, protein, metadata=True)
    interactions = []
    for (ligand_residue, protein_residue), types in ifp.items():
        for interaction, occurrences in types.items():
            for occurrence in occurrences:
                interactions.append({"ligand_residue": str(ligand_residue), "protein_residue": str(protein_residue),
                                     "interaction_type": interaction, **_clean(occurrence)})
    metrics = {"method": "ProLIF", "inputs": {**inputs, "receptor_sha256": receptor_hash},
               "interactions": interactions, "interaction_count": len(interactions),
               "fingerprint_keys": sorted({r["protein_residue"] + ":" + r["interaction_type"] for r in interactions}),
               "atom_index_scope": "parent_indices: ligand representative RDKit pose / inferred protein RDKit topology",
               "ligand_hydrogens": "RDKit AddHs addCoords; heavy atom coordinates retained",
               "protein_topology": "MDAnalysis coordinate bond inference; explicit receptor hydrogens required",
               "protein_selection": selection,
               "coordinate_unit": "angstrom", "biological_validation": "NOT_ESTABLISHED"}
    return metrics, [_write_json(directory / "prolif.json", metrics)]


def _fpsim2(rows, context, directory, software, cancel_event, deadline):
    Chem = _module(software, "rdkit.Chem")
    fps = _module(software, "FPSim2")
    io = _module(software, "FPSim2.io")
    refs = _references({}, context)
    if not refs:
        raise ChemistryHeld("valid_search_reference_required")
    valid, invalid = [], {}
    for row in rows:
        _check(cancel_event, deadline)
        try:
            mol = _smiles_mol(row, Chem)
            valid.append((row, Chem.MolToSmiles(mol, isomericSmiles=True)))
        except ChemistryHeld as error:
            invalid[row["candidate_id"]] = str(error)
    if not valid:
        return {}, invalid, []
    params = {"radius": 2, "fpSize": 2048}
    db = directory / "structures.h5"
    io.create_db_file([(smi, index + 1) for index, (_, smi) in enumerate(valid)], str(db),
                      "smiles", "Morgan", fp_params=params)
    mapping = _write_json(directory / "structure-index.json", {"fingerprint_type": "Morgan", "parameters": params,
        "records": [{"mol_id": index + 1, "candidate_id": row["candidate_id"], "smiles": smi,
                     "exact_structure_id": chemical_identity(smi)["exact_structure_id"]}
                    for index, (row, smi) in enumerate(valid)]})
    engine = fps.FPSim2Engine(str(db))
    threshold = context.get("similarity_threshold", 0.0)
    if isinstance(threshold, bool) or not isinstance(threshold, (int, float)) or not math.isfinite(threshold) or not 0 <= threshold <= 1:
        raise ValueError("similarity_threshold must be finite in 0..1")
    count = _count(context, "search_limit", len(valid), len(valid))
    output = {row["candidate_id"]: {"method": "FPSim2-Morgan", "fingerprint_parameters": params,
        "threshold": threshold, "best_tanimoto": None, "reference_scores": [], "search_limit": count,
        "structure_index_sha256": mapping["sha256"], "search_index_scope": "request-local registered structures"}
              for row, _ in valid}
    completed_queries, stop_status = 0, None
    for ref in refs:
        try:
            _check(cancel_event, deadline)
        except ChemistryStopped as error:
            stop_status = error.status
            break
        reference_id = chemical_identity(ref["smiles"])["exact_structure_id"]
        for hit in engine.top_k(ref["smiles"], count, threshold, n_workers=1):
            record = output[valid[int(hit["mol_id"]) - 1][0]["candidate_id"]]
            score = float(hit["coeff"])
            record["reference_scores"].append({"reference_exact_structure_id": reference_id, "tanimoto": score})
            record["best_tanimoto"] = max(record["best_tanimoto"] if record["best_tanimoto"] is not None else 0, score)
        completed_queries += 1
    for record in output.values():
        record.update(reference_queries_completed=completed_queries, reference_queries_requested=len(refs),
                      search_complete=completed_queries == len(refs), interruption=stop_status)
    return output, invalid, [_artifact(db), mapping]


def execute_chemistry(provider: str, rows: list[dict], context: dict, *, workdir: Path, cancel_event=None):
    """Execute real optional methods; report held rows without population removal.

    ``context`` supports candidate_limit, total_states<=2, conformers<=10,
    reference_limit<=5, budget_seconds<=28800, pH/ph_min/ph_max, reference_smiles,
    receptor_path/receptor_sha256, and FPSim2 similarity_threshold/search_limit.
    Counts may also be passed in ``context['counts']``. Larger requested counts
    are clamped to the immutable chemistry caps and reported in measurements.
    Purchase admission is verified by the registered-artifact worker upstream.
    """
    if provider not in CAPS:
        raise ValueError("unknown chemistry provider")
    if not isinstance(rows, list) or any(not isinstance(row, dict) for row in rows):
        raise ValueError("chemistry rows must be a list of objects")
    if not isinstance(context, dict) or not isinstance(context.get("counts", {}), dict):
        raise ValueError("chemistry context/counts must be objects")
    ids = [row.get("candidate_id") for row in rows]
    if any(not isinstance(cid, str) or not cid or len(cid) > 200 for cid in ids) or len(set(ids)) != len(ids):
        raise ValueError("distinct nonempty candidate IDs required")
    budget = context.get("budget_seconds", 28800)
    if isinstance(budget, bool) or not isinstance(budget, (int, float)) or not math.isfinite(budget) or not 0 < budget <= 28800:
        raise ValueError("chemistry budget_seconds must be positive and at most28800")
    started = time.monotonic()
    deadline = started + budget
    software = _software(provider)
    cap = CAPS[provider]
    if provider == "prolif" and context.get("on_demand") is True:
        cap = 1000
    limit = _count(context, "candidate_limit", cap, cap)
    result = {"provider": provider, "status": "held", "candidates": [], "software": software,
              "scientific_validation": "NOT_ESTABLISHED", "population_preserved": True}
    directory = Path(workdir).resolve() / (provider + "-" + uuid.uuid4().hex)
    directory.mkdir(parents=True)
    eligible = []
    receipts = {}
    controls = {chemical_identity(r.get("smiles"))["parent_id"] for r in rows
                if r.get("role") in CONTROL_ROLES or r.get("opposing_direction") is True}
    used = set()
    for row in rows:
        identity = chemical_identity(row.get("smiles"))
        receipt = {"candidate_id": row["candidate_id"], "status": "held", "metrics": {},
                   "input_exact_structure_id": identity["exact_structure_id"], "original_smiles": row.get("smiles")}
        receipts[row["candidate_id"]] = receipt
        parent = identity["parent_id"]
        if row.get("role") in CONTROL_ROLES or row.get("opposing_direction") is True or parent in controls:
            receipt["reason"] = "controls_separate"
        elif identity["structure_status"] != "valid" or identity["parent_ambiguous"]:
            receipt["reason"] = "invalid_or_ambiguous_structure"
        elif parent in used:
            receipt["reason"] = "duplicate_selection_parent"
        elif len(eligible) >= limit:
            receipt["reason"] = "candidate_compute_cap"
        else:
            used.add(parent)
            eligible.append(row)
    # Selection is frozen by the dispatcher. Respect its order; each result is
    # mapped back to the original request order, even when a row is held.
    try:
        _check(cancel_event, deadline)
        if provider == "fpsim2":
            scores, invalid, artifacts = _fpsim2(eligible, context, directory, software, cancel_event, deadline)
            result["artifacts"] = artifacts
            for row in eligible:
                receipt = receipts[row["candidate_id"]]
                if row["candidate_id"] in invalid:
                    receipt["reason"] = invalid[row["candidate_id"]]
                else:
                    metrics = scores[row["candidate_id"]]
                    receipt.update(status="complete" if metrics["search_complete"] else metrics["interruption"],
                                   metrics=metrics)
                    if metrics["interruption"]:
                        receipt["reason"] = metrics["interruption"]
        else:
            execute = {"dimorphite": _dimorphite, "cdpl": _cdpl,
                       "posebusters": _posebusters, "prolif": _prolif}[provider]
            resources = _PoseBustersResources() if provider == "posebusters" else None
            for index, row in enumerate(eligible):
                receipt = receipts[row["candidate_id"]]
                candidate_directory = directory / str(index)
                candidate_directory.mkdir()
                try:
                    _check(cancel_event, deadline)
                    if resources is not None:
                        metrics, artifacts = execute(row, context, candidate_directory, software, cancel_event,
                                                     deadline, _resources=resources)
                    else:
                        metrics, artifacts = execute(row, context, candidate_directory, software, cancel_event, deadline)
                    receipt.update(status="complete", metrics=_clean(metrics), artifacts=artifacts)
                except ChemistryStopped as error:
                    receipt.update(status=error.status, reason=str(error))
                except (ImportError, OSError, ChemistryHeld) as error:
                    receipt["reason"] = type(error).__name__ + ":" + str(error)[:500]
                except Exception as error:
                    receipt.update(status="failed", reason=type(error).__name__ + ":" + str(error)[:500])
    except ChemistryStopped as error:
        for row in eligible:
            receipts[row["candidate_id"]].update(status=error.status, reason=str(error))
    except (ImportError, OSError, ChemistryHeld) as error:
        for row in eligible:
            receipts[row["candidate_id"]]["reason"] = type(error).__name__ + ":" + str(error)[:500]
    except Exception as error:
        for row in eligible:
            receipts[row["candidate_id"]].update(status="failed", reason=type(error).__name__ + ":" + str(error)[:500])
    result["candidates"] = [receipts[cid] for cid in ids]
    statuses = [r["status"] for r in result["candidates"]]
    result["status"] = "complete" if statuses and all(s == "complete" for s in statuses) else "partial" if "complete" in statuses else "held"
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    result["measurements"] = {"wall_seconds": time.monotonic() - started,
        "process_lifetime_peak_rss_bytes": peak if sys.platform == "darwin" else peak * 1024,
        "memory_scope": "process lifetime high-water mark; not per-job incremental peak",
        "candidate_limit": limit, "selected_count": len(eligible), "input_count": len(rows),
        "completed_count": statuses.count("complete"), "requested_threads_per_native_operation": 1,
        "artifact_directory": str(directory)}
    result["receipt_artifact"] = _write_json(directory / "receipt.json", result)
    return result
