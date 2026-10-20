"""Fixed-coordinate, fixed-pocket research docking for a selected complex.

This path does not invent a UniProt identity. It preserves the selected
protein heavy atoms, qualifies the identical receptor/pocket with a query
ligand redock, and retains raw failures and coordinates as CAS artifacts.
"""

from copy import deepcopy
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import threading
import time

from .artifacts import atomic_bytes, canonical_json, sha256_file
from .engines import EngineRunner, ComputeBudgetExceeded, OwnedProcessCleanupError, terminate_owned
from .identity import chemical_identity, digest
from .native_prepare import PreparationUnavailable


def _docking_observation(source_row, execution, descriptor, receptor_sha256, reference_card_sha256=None):
    """Retain the native receipt and bind it to the exact snapshot molecule."""
    from .ranking import comparable_protocol
    protocol = deepcopy(execution["protocol"])
    return {
        **deepcopy(execution), "target_id": descriptor["target_id"],
        "candidate_id": source_row["candidate_id"],
        "exact_structure_id": chemical_identity(source_row["smiles"])["exact_structure_id"],
        "receptor_id": receptor_sha256, "receptor_sha256": receptor_sha256,
        "pocket_id": digest(descriptor["pocket_box"]), "protocol": protocol,
        "protocol_id": digest(comparable_protocol(protocol)), "unit": "kcal/mol",
        "score": execution.get("score", execution["affinity"]), "higher_is_better": False,
        "artifact_valid": True, "ranking_eligible": execution["status"] == "complete",
        "preparation_receipt_sha256": source_row.get("ligand_preparation_sha256"),
        "receptor_identity_card_sha256": reference_card_sha256,
        "biological_validation": "NOT_ESTABLISHED",
        "validation_scope": "Exact frozen source graph, preparation hash, raw pose atom inventory and finite saved coordinates",
    }


def freeze_docking_policy(config):
    settings = config.get("research_docking") or {}
    if not isinstance(settings, dict):
        raise ValueError("research_docking must be a server configuration object")
    engine = settings.get("engine", "vina")
    if engine not in {"vina", "unidock"}:
        raise ValueError("Only installed Vina or Uni-Dock may execute research docking")
    paths = settings.get("engine_paths", {})
    configured = paths.get(engine) if isinstance(paths, dict) else None
    binary = shutil.which(configured) if configured else None
    python = settings.get("preparation_python") or config.get("preparation_python") or sys.executable
    executable = Path(python).expanduser()
    return {
        "policy": "observed-complex-fixed-pocket-v1", "engine": engine,
        "status": "configured" if binary else "not_configured",
        "engine_binary_sha256": sha256_file(binary) if binary else None,
        "preparation_python_sha256": sha256_file(executable) if executable.is_file() else None,
        "qualification_sha256": settings.get("qualification_sha256"),
        "execution_source_sha256": sha256_file(__file__),
        "seed": 42, "n_poses": 9, "exhaustiveness": 16,
        "ligand_ph": None, "protein_heavy_atom_policy": "unchanged_observed_coordinates",
        "redocking_rmsd_threshold_angstrom": 2.0,
    }


def resolve_protein_context(service, descriptor, receptor_request, *, owner):
    """Freeze an explicitly chosen complex and prove its protein mapping.

    A deposited ligand is a redocking reference, not an activity observation.
    Missing/ambiguous structure choices hold docking while catalog research
    remains available under the user's original protein identity.
    """
    base = {"status": "held", "target_id": descriptor["target_id"],
            "requested_receptor": deepcopy(receptor_request),
            "hold_reasons": [], "bound_descriptor": None}
    if not receptor_request or not receptor_request.get("pdb_id") or not receptor_request.get("ligand_code"):
        return {**base, "hold_reasons": ["explicit_pdb_and_observed_ligand_selection_required"]}
    if receptor_request.get("explicit_pocket") is not None:
        return {**base, "hold_reasons": ["observed_ligand_required_for_fixed_pocket_redocking"]}
    if receptor_request.get("ph") not in (None, 7.4):
        return {**base, "hold_reasons": ["requested_protein_ph_not_supported_by_as_written_cpu_preparation"]}
    selection = {"ligand_code": receptor_request["ligand_code"].upper()}
    if receptor_request.get("chain"):
        selection["protein_chains"] = [receptor_request["chain"]]
        selection["ligand_chain"] = receptor_request["chain"]
    if receptor_request.get("assembly_id"):
        selection["assembly_id"] = receptor_request["assembly_id"]
    bound = service.research_structures.resolve(
        {"kind": "bound_complex", "format": "pdb_id", "value": receptor_request["pdb_id"],
         "selection": selection}, owner,
    )
    base["bound_descriptor"] = bound
    if bound.get("status") != "verified":
        return {**base, "hold_reasons": [bound.get("reason", "bound_complex_selection_required")]}
    try:
        from ..structure.identity import source_entities, assess
        if bound.get("original_structure_format") != "mmcif":
            raise ValueError("Deposited source mmCIF mapping is required")
        source_sha = bound["original_structure_sha256"]
        directory = service.root / "protein-context" / digest([descriptor["target_id"], source_sha, selection])
        directory.mkdir(parents=True, exist_ok=True)
        source_path = directory / "source.cif"
        atomic_bytes(source_path, service.artifacts.get_bytes(source_sha))
        entities = source_entities(source_path, [])
        accession = descriptor["components"][0]["accession"]
        mappings = []
        for chain in bound["selection"]["protein_chains"]:
            matching = [entity for entity in entities if chain in entity.get("chains", [])
                        and any(segment["accession"] == accession and segment["chain"] == chain
                                for segment in entity.get("reference_segments", []))]
            if len(matching) != 1:
                raise ValueError("Selected protein chain lacks an exact deposited observed mapping to the requested accession")
            mappings.append(assess(matching[0], chain, accession, "protein"))
        cpu_bound = deepcopy(bound)
        cpu_bound["structure_target_id"] = bound["target_id"]
        cpu_bound["target_id"] = descriptor["target_id"]
        cpu_bound["actual_protein_id"] = descriptor["target_id"]
        cpu_bound["taxon_id"] = descriptor.get("taxon_id")
        return {**base, "status": "verified", "bound_descriptor": cpu_bound,
                "observed_protein_mapping": mappings, "hold_reasons": deepcopy(bound.get("docking_hold_reasons", [])),
                "protein_preparation_ph": None, "requested_ph": receptor_request.get("ph"),
                "protonation_scope": "As-written deposited/template states; requested pH is not a validated protein protomer assignment"}
    except (ValueError, KeyError, OSError, ImportError) as error:
        return {**base, "hold_reasons": ["protein_structure_identity_not_verified"],
                "diagnostic": str(error)[:512]}


def _runtime_step(payload, work, python, *, stop, timeout=360):
    """Run only this installed module in the configured chemistry runtime."""
    stop()
    attempt = work / ("step-" + digest([payload, time.time_ns()])[:16])
    attempt.mkdir(parents=True)
    request_path, result_path = attempt / "request.json", attempt / "result.json"
    atomic_bytes(request_path, canonical_json(payload))
    environment = {
        **{key: os.environ[key] for key in ("PATH", "HOME", "TMPDIR", "LANG", "LC_ALL", "SYSTEMROOT", "WINDIR") if key in os.environ},
        "PYTHONPATH": str(Path(__file__).resolve().parents[2]),
        "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
    }
    started = time.monotonic()
    with (attempt / "execution.log").open("wb") as log:
        process = subprocess.Popen(
            [str(python), "-m", "skinpipe.platform.research_docking", str(request_path), str(result_path)],
            env=environment, stdout=log, stderr=log, start_new_session=True,
        )
        try:
            while process.poll() is None:
                stop()
                if time.monotonic() - started > timeout:
                    raise PreparationUnavailable("Research chemistry runtime exceeded its bounded step")
                time.sleep(0.05)
            stop()
        finally:
            terminate_owned(process)
    if process.returncode != 0 or not result_path.is_file():
        raise PreparationUnavailable("Research chemistry runtime did not produce a complete result")
    result = json.loads(result_path.read_bytes())
    result["runtime_compute_seconds"] = time.monotonic() - started
    if result.get("status") != "complete":
        raise PreparationUnavailable(result.get("reason", "Research chemistry step held"))
    return result


def _execute_engine_fenced(runner, request, artifacts, *, stop, progress, remaining_seconds):
    """Cancel owned native work while preserving its lease/budget stop reason."""
    stop()
    cancel, finished = threading.Event(), threading.Event()
    interruptions = []
    def checked_progress(value):
        try:
            stop()
        except BaseException:
            cancel.set()
            raise
        progress({"stage": "bound_complex_docking", **value})
    def watch():
        while not finished.wait(0.2):
            try:
                stop()
            except BaseException as error:
                interruptions.append(error)
                cancel.set()
                return
    monitor = threading.Thread(target=watch, daemon=True)
    monitor.start()
    try:
        try:
            result = runner.execute(request, fetch_artifact=artifacts.get_bytes,
                                    publish_artifact=artifacts.put_bytes, cancel_event=cancel,
                                    remaining_seconds=remaining_seconds, progress=checked_progress)
        except OwnedProcessCleanupError:
            # Cleanup uncertainty outranks the cancellation which triggered
            # it. Replacing it with Cancelled would falsely acknowledge free
            # resources even though the owned native group may remain live.
            raise
        except BaseException:
            # EngineCancelled describes the local kill, not the authoritative
            # coordinator fence that caused it. Preserve that original cause.
            if interruptions:
                raise interruptions[0]
            raise
    finally:
        finished.set()
        monitor.join()
    if interruptions:
        raise interruptions[0]
    return result


def _heavy_coordinates(text, *, pdbqt=False):
    positions = {}
    for line in text.splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        element = line.split()[-1] if pdbqt else (line[76:78].strip() or line[12:16].strip()[0])
        if element.upper() in {"H", "HD", "HS"}:
            continue
        key = (line[21:22], line[22:27], line[12:16].strip())
        if key in positions:
            raise ValueError("Duplicate selected protein atom identity")
        positions[key] = [float(line[start:start + 8]) for start in (30, 38, 46)]
    if not positions:
        raise ValueError("No protein heavy atoms available for receptor preparation")
    return positions


def _prepare_receptor(payload):
    from ..structure.clean import write_receptor_pdbqt
    source, output = Path(payload["source_path"]), Path(payload["output_path"])
    notes = write_receptor_pdbqt(source, [], output, observed_pdb=source)
    observed = _heavy_coordinates(source.read_text())
    prepared = _heavy_coordinates(output.read_text(), pdbqt=True)
    if set(observed) - set(prepared):
        raise ValueError("Receptor typing lost observed protein heavy atoms")
    deviations = [sum((a - b) ** 2 for a, b in zip(xyz, prepared[key])) ** 0.5 for key, xyz in observed.items()]
    maximum = max(deviations)
    if maximum > 0.001:
        raise ValueError("Receptor typing moved observed protein heavy atoms")
    return {"status": "complete", "path": str(output), "notes": notes,
            "observed_heavy_atom_count": len(observed), "max_observed_heavy_atom_displacement_angstrom": maximum,
            "source_sha256": sha256_file(source), "prepared_sha256": sha256_file(output)}


def _redock_geometry(payload):
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign
    from ..dock.poses import mols_from_docked_pdbqt, write_pose_sdf

    observed = Chem.SDMolSupplier(payload["observed_sdf"], removeHs=True)[0]
    poses = mols_from_docked_pdbqt(Path(payload["pose_path"]).read_text())
    if observed is None or poses is None or not poses.GetNumConformers():
        raise ValueError("Observed or redocked ligand geometry unavailable")
    poses = Chem.RemoveHs(poses)
    if (chemical_identity(Chem.MolToSmiles(observed))["exact_structure_id"]
            != chemical_identity(Chem.MolToSmiles(poses))["exact_structure_id"]):
        raise ValueError("Observed and redocked ligand exact graphs differ")
    values = [float(rdMolAlign.CalcRMS(poses, observed, prbId=cid, refId=0, maxMatches=10000))
              for cid in range(poses.GetNumConformers())]
    scores = payload["pose_scores"]
    if len(values) != len(scores):
        raise ValueError("Redocking pose-score correspondence differs")
    top = min(range(len(scores)), key=scores.__getitem__)
    output = Path(payload["output_sdf"])
    if not write_pose_sdf(poses, top, output, {"scope": "query redocking same receptor and fixed pocket"}):
        raise ValueError("Could not save the top redocking geometry")
    return {"status": "complete", "rmsd_angstrom": values[top], "all_pose_rmsd_angstrom": values,
            "best_pose_rmsd_angstrom": min(values), "top_pose_index": top + 1,
            "pass": values[top] <= 2.0, "threshold_angstrom": 2.0,
            "rmsd_policy": "symmetry-corrected heavy-atom RMSD in original receptor frame; no alignment",
            "top_pose_sdf": str(output), "atom_graph_mapping": "exact_graph_required",
            "interpretation": "Redocking checks this structure/pocket/protocol; it does not establish analogue binding or skin benefit"}


def _candidate_geometry(payload):
    from rdkit import Chem
    from ..dock.poses import mols_from_docked_pdbqt, write_pose_sdf
    poses = mols_from_docked_pdbqt(Path(payload["pose_path"]).read_text())
    if poses is None or not poses.GetNumConformers():
        raise ValueError("Docked candidate geometry unavailable")
    poses = Chem.RemoveHs(poses)
    if (chemical_identity(Chem.MolToSmiles(poses))["exact_structure_id"]
            != chemical_identity(payload["source_smiles"])["exact_structure_id"]):
        raise ValueError("Reconstructed candidate pose differs from its exact source graph")
    scores = payload["pose_scores"]
    if poses.GetNumConformers() != len(scores):
        raise ValueError("Candidate pose-score correspondence differs")
    top = min(range(len(scores)), key=scores.__getitem__)
    if not write_pose_sdf(poses, top, Path(payload["output_sdf"]), {"scope": "same fixed receptor and pocket docking"}):
        raise ValueError("Could not save candidate docking geometry")
    return {"status": "complete", "top_pose_sdf": payload["output_sdf"], "top_pose_index": top + 1,
            "atom_graph_mapping": "exact_graph_required"}


def _prepare_ligands(payload):
    """Use the native ligand preparer without inventing a pH transform."""
    from ..dock.ligprep import prepare_ligand_pdbqt, read_preparation_receipt
    ligands = payload.get("ligands")
    if not isinstance(ligands, list) or not 1 <= len(ligands) <= 32:
        raise ValueError("Research ligand preparation requires a bounded 1..32 batch")
    result = []
    for item in ligands:
        row = {"candidate_id": item["candidate_id"]}
        identity = chemical_identity(item["smiles"])
        if identity.get("unassigned_stereo") or identity.get("unassigned_bond_stereo") or identity.get("enhanced_stereo"):
            result.append({**row, "reason": "explicit_isomer_required_for_3d_preparation"})
            continue
        path = prepare_ligand_pdbqt(item["smiles"], item["candidate_id"], Path(payload["workdir"]), ph=None)
        if path is None or read_preparation_receipt(path, ligand_id=item["candidate_id"]) is None:
            result.append({**row, "reason": "ligand_preparation_failed"})
            continue
        result.append({**row, "path": str(path), "receipt_path": str(path.with_suffix(".preparation.json"))})
    return {"status": "complete", "ligands": result, "ph": None, "policy": "exact-as-written-research-ligands-v1"}


def evaluate_bound_complex(rows, descriptor, selected, service, config, request, *, stop, progress):
    if request.get('screening_policy') is not None:
        from .routine_docking import evaluate_routine_complex
        return evaluate_routine_complex(rows, descriptor, selected, service, config,
                                        request, stop=stop, progress=progress)
    started = time.monotonic()
    rows = deepcopy(rows)
    frozen = request.get("research_docking_policy") or freeze_docking_policy(config)
    receipt = {
        "policy": frozen, "status": "held", "receptor_valid": False,
        "protein_structure_sha256": descriptor.get("protein_structure_sha256"),
        "original_protein_heavy_coordinates_sha256": descriptor.get("original_protein_heavy_coordinates_sha256"),
        "pocket_box": deepcopy(descriptor.get("pocket_box")),
        "prepared_count": 0, "selected_count": 0, "completed_count": 0,
        "redocking": {"status": "not_run"}, "candidate_results": [],
    }
    if descriptor.get("docking_hold_reasons"):
        return rows, {**receipt, "hold_reasons": deepcopy(descriptor["docking_hold_reasons"])}
    if frozen != freeze_docking_policy(config):
        return rows, {**receipt, "hold_reasons": ["configured_docking_cohort_changed_since_submission"]}
    if frozen["status"] != "configured":
        return rows, {**receipt, "hold_reasons": ["docking_engine_not_configured"]}
    settings = config.get("research_docking") or {}
    python = settings.get("preparation_python") or config.get("preparation_python") or sys.executable
    work = service.root / "research-docking" / request["discovery_run_id"]
    work.mkdir(parents=True, exist_ok=True)
    source = work / "selected-protein.pdb"
    observed = work / "observed-query.sdf"
    try:
        stop()
        atomic_bytes(source, service.artifacts.get_bytes(descriptor["protein_structure_sha256"]))
        seed = next(r for r in rows if r.get("role") == "query_seed" or r.get("redock_reference") is True)
        atomic_bytes(observed, service.artifacts.get_bytes(seed["observed_sdf_sha256"]))
        receptor = _runtime_step({"operation": "receptor", "source_path": str(source), "output_path": str(work / "receptor.pdbqt")}, work, python, stop=stop)
        receptor_artifact = service.artifacts.put_file(receptor["path"], media_type="chemical/x-pdbqt")
        receipt.update(receptor_valid=True, receptor_preparation=receptor,
                       receptor_sha256=receptor_artifact["sha256"])
        reference_card = {"schema": "skinpipe-observed-complex-receptor-card-v1",
            "artifacts": {"pdbqt": {"present": True, "sha256": receptor_artifact["sha256"]},
                          "cocrystal_sdf": {"present": True, "sha256": seed["observed_sdf_sha256"], "path": "observed-query.sdf"}},
            "receptor": {"cocrystal_sdf": "observed-query.sdf"},
            "selected_structure_sha256": descriptor["protein_structure_sha256"],
            "pocket_box": descriptor["pocket_box"], "biological_validation": "NOT_ESTABLISHED"}
        reference_card_sha = service.artifacts.put_bytes(canonical_json(reference_card), media_type="application/json")["sha256"]
        receipt["receptor_identity_card_sha256"] = reference_card_sha
        maximum = request.get("generation", {}).get("max_3d", 50)
        chosen = [seed] + [r for cid in selected[:maximum] for r in rows if r["candidate_id"] == cid]
        prepared = {"ligands": []}
        for offset in range(0, len(chosen), 32):
            part = _runtime_step({"operation": "ligands", "workdir": str(work),
                                  "ligands": [{"candidate_id": r["candidate_id"], "smiles": r["smiles"]} for r in chosen[offset:offset + 32]]},
                                 work, python, timeout=360, stop=stop)
            prepared["ligands"].extend(part["ligands"])
        ligands = []
        from ..dock.ligprep import read_preparation_receipt
        chosen_by_id = {r["candidate_id"]: r for r in chosen}
        for item in prepared["ligands"]:
            if not item.get("path"):
                receipt.setdefault("preparation_holds", []).append(item)
                continue
            path = Path(item["path"]).resolve()
            if not path.is_relative_to(work.resolve()):
                raise ValueError("Prepared ligand path escaped the owned run directory")
            preparation = read_preparation_receipt(path, ligand_id=item["candidate_id"])
            source_row = chosen_by_id[item["candidate_id"]]
            if (preparation is None or preparation["request"]["ph"] is not None
                    or preparation["input_smiles"] != source_row["smiles"]
                    or chemical_identity(preparation["effective_smiles"])["exact_structure_id"]
                    != chemical_identity(source_row["smiles"])["exact_structure_id"]):
                raise ValueError("Prepared research ligand differs from its frozen exact source graph")
            artifact = service.artifacts.put_file(path, media_type="chemical/x-pdbqt")
            source_row["ligand_preparation_sha256"] = service.artifacts.put_bytes(canonical_json(preparation), media_type="application/json")["sha256"]
            ligands.append({"candidate_id": item["candidate_id"], "pdbqt_sha256": artifact["sha256"]})
        receipt["prepared_count"] = len(ligands)
        query = next((r for r in ligands if r["candidate_id"] == seed["candidate_id"]), None)
        if query is None:
            return rows, {**receipt, "hold_reasons": ["query_redocking_preparation_failed"]}
        engine_request = {"engine": frozen["engine"], "receptor_sha256": receptor_artifact["sha256"],
                          "center": descriptor["pocket_box"]["center"], "size": descriptor["pocket_box"]["size"],
                          "seed": frozen["seed"], "n_poses": frozen["n_poses"], "exhaustiveness": frozen["exhaustiveness"]}
        runner = EngineRunner(work / "engine", settings["engine_paths"], gpu_uuid=settings.get("gpu_uuid"))
        def dock(batch):
            remaining = max(0, request.get("compute_budget_seconds", 28800) - (time.monotonic() - started))
            return _execute_engine_fenced(runner, {**engine_request, "ligands": batch}, service.artifacts,
                                          stop=stop, progress=progress, remaining_seconds=remaining)
        redock = dock([query])
        receipt["redocking_engine_result"] = redock
        redocked = redock["candidates"][0]
        if redocked.get("status") not in {"complete", "review_required"} or not redocked.get("pose_sha256"):
            return rows, {**receipt, "hold_reasons": ["query_redocking_engine_failed"]}
        def observation(source_row, execution):
            return _docking_observation(source_row, execution, descriptor, receptor_artifact["sha256"], reference_card_sha)
        seed["docking_execution"] = deepcopy(redocked)
        seed["docking"] = [observation(seed, redocked)]
        seed["primary_docking"] = deepcopy(seed["docking"][0])
        pose = work / "redocked-query.pdbqt"
        atomic_bytes(pose, service.artifacts.get_bytes(redocked["pose_sha256"]))
        geometry = _runtime_step({"operation": "redock_geometry", "observed_sdf": str(observed), "pose_path": str(pose),
                                  "pose_scores": redocked["pose_scores"], "output_sdf": str(work / "redocked-query.sdf")}, work, python, stop=stop)
        geometry["top_pose_sdf_sha256"] = service.artifacts.put_file(geometry.pop("top_pose_sdf"), media_type="chemical/x-mdl-sdfile")["sha256"]
        receipt["redocking"] = geometry
        seed["query_redocking"] = deepcopy(geometry)
        seed["docking_execution"]["sdf_artifact"] = {"sha256": geometry["top_pose_sdf_sha256"]}
        seed["docking_execution"]["pose_geometry"] = deepcopy(geometry)
        seed["docking"] = [observation(seed, seed["docking_execution"])]
        seed["primary_docking"] = deepcopy(seed["docking"][0])
        if not geometry["pass"]:
            return rows, {**receipt, "hold_reasons": ["query_redocking_top_pose_rmsd_above_threshold"]}
        candidates = [r for r in ligands if r["candidate_id"] != seed["candidate_id"]]
        receipt["selected_count"] = len(candidates)
        if not candidates:
            return rows, {**receipt, "status": "complete", "hold_reasons": ["no_prepared_new_candidates"]}
        result = dock(candidates)
        receipt["candidate_results"] = result["candidates"]
        receipt["protocol"] = result["protocol"]
        by_id = {r["candidate_id"]: r for r in rows}
        for result_row in result["candidates"]:
            row = by_id[result_row["candidate_id"]]
            row["docking_execution"] = deepcopy(result_row)
            if result_row.get("status") in {"complete", "review_required"}:
                pose_path = work / ("candidate-" + digest(row["candidate_id"])[:24] + ".pdbqt")
                atomic_bytes(pose_path, service.artifacts.get_bytes(result_row["pose_sha256"]))
                try:
                    geometry = _runtime_step({"operation": "candidate_geometry", "pose_path": str(pose_path),
                                              "source_smiles": row["smiles"], "pose_scores": result_row["pose_scores"],
                                              "output_sdf": str(pose_path.with_suffix(".sdf"))}, work, python, stop=stop)
                    sdf = service.artifacts.put_file(geometry.pop("top_pose_sdf"), media_type="chemical/x-mdl-sdfile")
                    row["docking_execution"]["sdf_artifact"] = sdf
                    row["docking_execution"]["pose_geometry"] = geometry
                except PreparationUnavailable as error:
                    row["docking_execution"]["sdf_status"] = "unavailable: " + str(error)[:256]
                row["docking"] = [observation(row, row["docking_execution"])]
                row["primary_docking"] = row["docking"][0]
                receipt["completed_count"] += 1
        receipt["status"] = "complete"
        return rows, receipt
    except ComputeBudgetExceeded:
        raise TimeoutError("budget_deferred") from None
    except TimeoutError:
        raise
    except (PreparationUnavailable, ValueError, OSError, ImportError) as error:
        return rows, {**receipt, "hold_reasons": ["research_docking_execution_held"], "error_type": type(error).__name__, "diagnostic": str(error)[:512]}


def main():
    request_path, result_path = map(Path, sys.argv[1:])
    try:
        payload = json.loads(request_path.read_bytes())
        if payload["operation"] == "receptor":
            result = _prepare_receptor(payload)
        elif payload["operation"] == "redock_geometry":
            result = _redock_geometry(payload)
        elif payload["operation"] == "candidate_geometry":
            result = _candidate_geometry(payload)
        elif payload["operation"] == "ligands":
            result = _prepare_ligands(payload)
        elif payload['operation'] == 'routine_quality':
            from .routine_docking import pose_quality
            result = pose_quality(payload)
        else:
            raise ValueError("Unknown installed research chemistry operation")
    except Exception as error:
        result = {"status": "held", "reason": type(error).__name__ + ": " + str(error)[:512]}
    atomic_bytes(result_path, canonical_json(result))


if __name__ == "__main__":
    main()
