"""The Vina calls themselves: one ligand, and a resumable parallel screen.

Split out of :mod:`skinpipe.dock.vina` purely for file size; ``dock_one`` and
``dock_many`` are re-exported there.

Parallelism is by process, not thread: each worker builds its own ``Vina``
object and its own affinity maps, because the maps are C++ state that cannot be
shared across processes and Vina's own ``cpu`` threading scales far worse than
running independent ligands side by side. Each worker therefore runs with
``cpu=1``.

Resumption is per ligand and keyed on ``results/<ligand_id>.json``, which is
written by the worker as soon as that ligand finishes. A cached entry is only
trusted when its pose SDFs are still on disk *and* it was produced under the
same docking conditions -- same receptor bytes, chain, box and search settings
(see :mod:`skinpipe.dock.cachekey`). Anything else is re-docked, with the
difference named in a WARNING.
"""
from __future__ import annotations

import json
import logging
import math
import multiprocessing
import os
import tempfile
import uuid
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict
from pathlib import Path
from typing import TYPE_CHECKING, Any, Optional

from ..types import DockPose, LigandRecord, Receptor
from . import artifacts, cachekey, software
from . import poses as poseio
from .ligprep import prepare_ligand_pdbqt

if TYPE_CHECKING:  # pragma: no cover - import cost only paid by type checkers
    import pandas as pd

log = logging.getLogger(__name__)


class DockingEngineError(RuntimeError):
    'The engine failed to run. Not "this ligand has no pose".'


def _record_raw_output(workdir: Path, text: str, attempt: dict[str, Any]) -> None:
    """Preserve all returned poses, including finite scores outside our window."""
    from . import unidock
    directory = Path(workdir) / "raw_docking"
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / (uuid.uuid4().hex + ".pdbqt")
    with path.open("x", encoding="utf-8") as handle:
        handle.write(text)
    attempt["raw_output_artifact"] = artifacts.describe(workdir, path)
    attempt["score_review_policy"] = {
        "window": list(unidock.PLAUSIBLE_AFFINITY),
        "scope": "operational quality/routing window; not a universal physical affinity bound",
        "outside_window": "raw output retained; excluded from the current scored-pose subset"}
    scores = []
    try:
        pose_scores = poseio.scored_pose_remarks(text)
    except ValueError as exc:
        # Keep the raw bytes, but do not renumber a later valid score as the
        # failed MODEL's score in the execution receipt either.
        attempt["score_parse_error"] = str(exc)
        attempt["returned_pose_scores"] = []
        return
    for index, (value, _, _) in enumerate(pose_scores, start=1):
        finite = value is not None and math.isfinite(value)
        inside = finite and unidock.is_plausible(value)
        scores.append({"pose_index": index, "affinity": value if finite else None,
                       "raw_affinity": str(value), "finite": finite,
                       "within_operational_window": inside,
                       "flag": "within_window" if inside else
                               "finite_outside_operational_window" if finite else "nonfinite_affinity"})
    attempt["returned_pose_scores"] = scores


def _record_preparation_provenance(trace: dict[str, Any], workdir: Path,
                                   path: Path, ligand: LigandRecord) -> None:
    from . import ligprep
    receipt = ligprep.read_preparation_receipt(path, ligand_id=ligand.ligand_id)
    if receipt is None or receipt['request']['smiles'] != ligand.smiles:
        trace['prepared_ligand_provenance_status'] = 'unverified'
        return
    trace['prepared_ligand_provenance_status'] = 'verified_generation_record'
    trace['prepared_ligand_provenance_artifact'] = artifacts.describe(workdir, ligprep._receipt_path(path))
    trace['protonation'] = receipt['protonation']
    trace['preparation_reusable'] = receipt['reusable']


def dock_one(
    receptor: Receptor,
    ligand: LigandRecord,
    workdir: Path,
    exhaustiveness: int = 16,
    n_poses: int = 9,
    seed: int = 42,
    engine_name: str = "vina",
    unidock_search_mode: str = "balance",
    unidock_bin: Optional[str] = None,
    execution: Optional[dict[str, Any]] = None,
) -> list[DockPose]:
    'Dock one ligand and write every pose as a verified SDF.'
    from vina import Vina

    if engine_name not in ("vina", "unidock", "hybrid"):
        raise ValueError(f"engine must be vina, unidock or hybrid, got {engine_name!r}")
    if execution is not None and not isinstance(execution, dict):
        raise TypeError("execution must be a dictionary")
    trace = execution if execution is not None else {}
    trace.clear()
    trace.update(requested_engine=engine_name, actual_engine=None,
                 seed=int(seed), n_poses_requested=int(n_poses),
                 n_poses_returned=0,
                 vina_exhaustiveness=int(exhaustiveness),
                 unidock_search_mode=(unidock_search_mode if engine_name != "vina" else None),
                 unidock_binary_requested=unidock_bin, route_reason="selected_engine",
                 state="not_started", attempts=[])

    def completed(text: str, backend: str, attempt: dict[str, Any]) -> list[DockPose]:
        if "raw_output_artifact" not in attempt and execution is not None:
            _record_raw_output(workdir, text, attempt)
        try:
            docked = _poses_from_docked_text(text, ligand, workdir, receptor)
        except DockingEngineError as exc:
            attempt.update(state="error", error_type=type(exc).__name__, error=str(exc))
            trace.update(actual_engine=None, state="engine_error")
            raise
        attempt["state"] = "completed"
        trace.update(actual_engine=backend, state="completed",
                     engine_identity=software.engine_identity(backend, unidock_bin))
        trace['docking_conditions'] = cachekey.dock_conditions(
            receptor, exhaustiveness, n_poses, seed, engine=backend,
            search_mode=unidock_search_mode, unidock_bin=unidock_bin)
        trace["n_poses_returned"] = len(docked)
        return docked

    workdir = Path(workdir)
    if receptor.pdbqt is None or not Path(receptor.pdbqt).is_file():
        raise FileNotFoundError(f"receptor PDBQT missing for {receptor.pdb_id}")

    preparation_execution: dict[str, Any] = {}
    lig_pdbqt = prepare_ligand_pdbqt(ligand.smiles, ligand.ligand_id, workdir,
                                    execution=preparation_execution)
    if lig_pdbqt is None:
        if preparation_execution.get('state') == 'held':
            trace['preparation'] = preparation_execution
            trace['state'] = 'ligand_preparation_held'
            raise DockingEngineError(preparation_execution['reason'])
        trace["state"] = "ligand_preparation_rejected"
        return []
    if execution is not None:
        trace["prepared_ligand_sha256"] = artifacts.sha256(Path(lig_pdbqt))
        trace["receptor_pdbqt_sha256"] = artifacts.sha256(Path(receptor.pdbqt))
        trace["prepared_ligand_artifact"] = artifacts.describe(workdir, Path(lig_pdbqt))
        _record_preparation_provenance(trace, workdir, Path(lig_pdbqt), ligand)

    if engine_name == "hybrid":
        # The same two gates `_hybrid_gpu_pass` applies to a batch, for one
        # ligand. The redock control has to run through the routing the screen
        # actually uses, or it confirms a path the screen does not take -- and
        # a crystallographic macrocycle would otherwise be handed to an engine
        # that cannot score its glue atoms, silently, in the one job whose
        # whole purpose is to be checked.
        from . import unidock as ud

        text = Path(lig_pdbqt).read_text(encoding="utf-8")
        if ud.has_glue_atoms(text):
            trace["route_reason"] = "macrocycle_glue_cpu_route"
            log.info("%s: macrocycle glue atoms; the control runs on CPU Vina, "
                     "as the batch would route it", ligand.ligand_id)
        else:
            attempt = {"engine": "unidock", "state": "started",
                       "engine_identity": software.engine_identity("unidock", unidock_bin)}
            trace["attempts"].append(attempt)
            try:
                docked_text = ud.dock_one_text(
                    receptor, lig_pdbqt, search_mode=unidock_search_mode,
                    n_poses=int(n_poses), seed=int(seed), binary=unidock_bin)
                if execution is not None:
                    _record_raw_output(workdir, docked_text, attempt)
                affinity = ud.best_affinity(docked_text)
                attempt.update(state="completed", best_affinity=affinity)
                if affinity is not None and ud.is_plausible(affinity):
                    return completed(docked_text, "unidock", attempt)
                trace["route_reason"] = "unidock_implausible_score_cpu_fallback"
                log.warning("%s: unidock returned %r, outside %s; the control "
                            "falls back to CPU Vina", ligand.ligand_id,
                            affinity, ud.PLAUSIBLE_AFFINITY)
            except (ud.UnidockUnavailable, FileNotFoundError, ValueError) as exc:
                attempt.update(state="error", error_type=type(exc).__name__, error=str(exc))
                trace["actual_engine"] = None
                trace["route_reason"] = "unidock_error_cpu_fallback"
                log.warning("%s: the GPU control failed (%s); falling back to "
                            "CPU Vina", ligand.ligand_id, exc)
        # fall through to the Vina block below

    elif engine_name == "unidock":
        from . import unidock as ud
        attempt = {"engine": "unidock", "state": "started",
                   "engine_identity": software.engine_identity("unidock", unidock_bin)}
        trace["attempts"].append(attempt)
        try:
            docked_text = ud.dock_one_text(
                receptor, lig_pdbqt, search_mode=unidock_search_mode,
                n_poses=int(n_poses), seed=int(seed), binary=unidock_bin)
        except (ud.UnidockUnavailable, FileNotFoundError, ValueError) as exc:
            attempt.update(state="error", error_type=type(exc).__name__, error=str(exc))
            trace["state"] = "engine_error"




            log.error("%s: unidock failed against %s: %s",
                      ligand.ligand_id, receptor.pdb_id, exc)
            raise DockingEngineError(
                f"unidock failed for {ligand.ligand_id} against "
                f"{receptor.pdb_id}: {type(exc).__name__}: {exc}") from exc
        return completed(docked_text, "unidock", attempt)

    attempt = {"engine": "vina", "state": "started",
               "engine_identity": software.engine_identity("vina")}
    trace["attempts"].append(attempt)
    try:
        engine = Vina(sf_name="vina", cpu=1, seed=int(seed), verbosity=0)
        engine.set_receptor(str(receptor.pdbqt))
        engine.set_ligand_from_string(lig_pdbqt.read_text())
        engine.compute_vina_maps(center=list(receptor.center), box_size=list(receptor.size))
        engine.dock(exhaustiveness=int(exhaustiveness), n_poses=int(n_poses))
        docked_text = engine.poses(n_poses=int(n_poses))
    except (RuntimeError, ValueError, TypeError) as exc:
        attempt.update(state="error", error_type=type(exc).__name__, error=str(exc))
        trace["state"] = "engine_error"
        log.error("%s: vina failed against %s: %s", ligand.ligand_id, receptor.pdb_id, exc)
        raise DockingEngineError(
            f"vina failed for {ligand.ligand_id} against {receptor.pdb_id}: "
            f"{type(exc).__name__}: {exc}") from exc

    return completed(docked_text, "vina", attempt)


def _keep_plausible(poses: list[DockPose], ligand_id: str,
                    where: str) -> list[DockPose]:
    'Select the operational scored-pose subset and log excluded values.'
    from . import unidock as ud

    keep, dropped = [], []
    for pose in poses:
        (keep if (pose.affinity is not None and ud.is_plausible(pose.affinity))
         else dropped).append(pose)
    if dropped:
        log.warning("%s: %d of %d pose(s) from %s are outside the operational "
                    "review window %s (%s); excluded from scored poses", ligand_id, len(dropped),
                    len(poses), where, ud.PLAUSIBLE_AFFINITY,
                    ", ".join(f"{p.affinity}" for p in dropped[:4]))
    return keep


def _poses_from_docked_text(
    docked_text: str,
    ligand: LigandRecord,
    workdir: Path,
    receptor: Receptor,
) -> list[DockPose]:
    try:
        scores = poseio.scored_pose_remarks(docked_text)
    except ValueError as exc:
        raise DockingEngineError(f"{ligand.ligand_id}: invalid pose/score mapping: {exc}") from exc
    mol = poseio.mols_from_docked_pdbqt(docked_text)
    if mol is None or not scores:
        raise DockingEngineError(f"{ligand.ligand_id}: docking output could not be rebuilt as scored poses")
    conformers = list(mol.GetConformers())
    if len(conformers) != len(scores):
        raise DockingEngineError(
            f"{ligand.ligand_id}: docking output has {len(scores)} scored MODELs "
            f"but {len(conformers)} rebuilt conformers")

    heavy = mol.GetNumHeavyAtoms()
    pose_dir = workdir / "poses" / ligand.ligand_id
    out: list[DockPose] = []
    for index, conf in enumerate(conformers):
        affinity, rmsd_lb, rmsd_ub = scores[index]
        sdf_path = poseio.write_pose_sdf(
            mol,
            conf.GetId(),
            pose_dir / f"{ligand.ligand_id}_pose{index + 1}.sdf",
            props={
                "ligand_id": ligand.ligand_id,
                "name": ligand.name or "",
                "receptor": f"{receptor.pdb_id}_{receptor.chain}",
                "vina_affinity": f"{affinity:.3f}",
                "pose_index": index + 1,
            },
        )
        if sdf_path is None:
            raise DockingEngineError(
                f"{ligand.ligand_id} pose {index + 1}: SDF did not round-trip; "
                "refusing to promote a later pose")
        out.append(
            DockPose(
                ligand_id=ligand.ligand_id,
                pose_index=index + 1,
                affinity=float(affinity),
                sdf_path=sdf_path,
                rmsd_lb=float(rmsd_lb),
                rmsd_ub=float(rmsd_ub),
                ligand_efficiency=(-float(affinity) / heavy) if heavy else None,
            )
        )
    return _keep_plausible(out, ligand.ligand_id, "this run")


def _result_path(workdir: Path, ligand_id: str) -> Path:
    from .ligprep import valid_ligand_id
    if not valid_ligand_id(ligand_id):
        raise ValueError('ligand id must be a filename component')
    return Path(workdir) / "results" / f"{ligand_id}.json"


def _load_cached(workdir: Path, ligand: LigandRecord,
                 conditions: dict[str, Any],
                 quiet: bool = False,
                 execution: Optional[dict[str, Any]] = None) -> Optional[list[DockPose]]:
    """Cached poses for this exact ligand under these exact conditions.

    Returns ``None`` when the ligand must be re-docked, having said at WARNING
    which condition changed. Entries are invalidated by a corrupt file, a
    different SMILES (so correcting a bad structure in the library re-docks it
    instead of resurrecting the old result), an unbound, changed or missing pose SDF, and any
    difference in ``conditions`` -- receptor bytes, chain, box, exhaustiveness,
    n_poses or seed.

    Worked example: a screen run against 5M8L then re-run with ``--pdb 2Y9X``
    logs ``pdb_id '5M8L' -> '2Y9X'; receptor_sha256 ... -> ...`` and docks
    every ligand again, instead of reprinting the 5M8L affinities under the
    2Y9X heading.

    An empty pose list can be a completed operational answer under these
    settings. It does not establish that the ligand cannot bind.
    """
    path = _result_path(workdir, ligand.ligand_id)
    if not path.is_file():
        return None
    try:
        payload = json.loads(path.read_text())
    except (json.JSONDecodeError, OSError) as exc:
        log.warning("%s: unreadable cached result (%s); re-docking", ligand.ligand_id, exc)
        return None
    if not isinstance(payload, dict) or payload.get("smiles") != ligand.smiles:
        log.info("%s: cached result is for a different SMILES; re-docking", ligand.ligand_id)
        return None
    changed = cachekey.describe_mismatch(payload.get("conditions"), conditions)
    if changed:
        # `quiet` is the caller saying "I have another set of conditions to
        # try"; the miss is then expected and says nothing about the cache.
        (log.debug if quiet else log.warning)(
            "%s: cached poses were produced under different docking "
            "conditions (%s); re-docking", ligand.ligand_id, "; ".join(changed))
        return None

    rows = payload.get("poses")
    if not isinstance(rows, list):
        log.warning("%s: cached poses are not a list; re-docking", ligand.ligand_id)
        return None
    if payload.get("schema_version") not in (None, 2, 3):
        log.warning("%s: unsupported pose cache schema; re-docking", ligand.ligand_id)
        return None
    # Old absolute paths do not bind the bytes that produced the affinity.
    # Preserve the record on disk, but recompute its poses rather than infer a
    # hash from whatever happens to occupy that path now. An empty result has
    # no pose artifact; its existing conditions still determine reuse.
    if payload.get("schema_version") not in (2, 3) and rows:
        log.warning("%s: cached poses have no portable artifact identity; re-docking",
                    ligand.ligand_id)
        return None
    out: list[DockPose] = []
    try:
        for row in rows:
            if not isinstance(row, dict) or "pose_artifact" not in row:
                raise ValueError("missing pose artifact identity")
            sdf = artifacts.resolve(workdir, row.get("sdf_path"), row["pose_artifact"])
            fields = {k: v for k, v in row.items() if k != "pose_artifact"}
            out.append(DockPose(**{**fields, "sdf_path": sdf}))
    except (OSError, ValueError, TypeError, KeyError) as exc:
        log.warning("%s: cached pose artifact cannot be verified (%s); re-docking",
                    ligand.ligand_id, exc)
        return None
    if execution is not None:
        trace = payload.get("execution")
        if not isinstance(trace, dict) or trace.get("state") not in ("completed", "ligand_preparation_rejected"):
            log.warning("%s: cache has no verified execution receipt; re-docking", ligand.ligand_id)
            return None
        if trace["state"] == "completed":
            if (trace.get("actual_engine") != conditions.get("engine") or
                    trace.get("engine_identity") != conditions.get("software_identity") or
                    trace.get('docking_conditions', conditions) != conditions):
                log.warning("%s: cache engine receipt disagrees with conditions; re-docking", ligand.ligand_id)
                return None
            try:
                prepared = trace["prepared_ligand_artifact"]
                path = artifacts.resolve(workdir, prepared["path"], prepared)
                if artifacts.sha256(path) != trace.get("prepared_ligand_sha256"):
                    return None
                from . import ligprep
                generation = trace['prepared_ligand_provenance_artifact']
                generation_path = artifacts.resolve(workdir, generation['path'], generation)
                request = ligprep.preparation_request(ligand.smiles,
                    software_identity=conditions['software_identity']['preparation'])
                if (generation_path != ligprep._receipt_path(path).resolve() or
                        ligprep.read_preparation_receipt(path, request, ligand.ligand_id, for_reuse=True) is None):
                    return None
                attempts = trace.get("attempts", [])
                if not isinstance(attempts, list) or any(not isinstance(item, dict) for item in attempts):
                    raise ValueError("malformed execution attempts")
                for attempt in attempts:
                    raw = attempt.get("raw_output_artifact")
                    if raw is not None:
                        artifacts.resolve(workdir, raw["path"], raw)
            except (OSError, ValueError, KeyError, TypeError):
                log.warning("%s: prepared ligand identity cannot be verified; re-docking", ligand.ligand_id)
                return None
        elif rows or trace.get("actual_engine") is not None:
            return None
        execution.update(trace)
        execution['docking_conditions'] = payload['conditions']
        execution["cache_hit"] = True
    return out


def _save_cached(workdir: Path, ligand: LigandRecord, docked: list[DockPose],
                 conditions: dict[str, Any],
                 execution: Optional[dict[str, Any]] = None) -> None:
    """Write the poses *and* the conditions that produced them.

    SDFs use paths relative to the docking workdir and full byte hashes, so the
    package can move and the reader rejects changed files. The conditions
    travel with the record so a later run can say what changed,
    rather than a reader having to trust that the directory was only ever used
    for one receptor.
    """
    path = _result_path(workdir, ligand.ligand_id)
    path.parent.mkdir(parents=True, exist_ok=True)
    if execution is not None:
        execution['docking_conditions'] = conditions
    rows = _json_rows(docked)
    for row, pose in zip(rows, docked):
        entry = artifacts.describe(workdir, pose.sdf_path)
        row["pose_artifact"] = entry
        row["sdf_path"] = entry["path"] if entry else None
    payload = {
        "schema_version": 3 if execution is not None else 2,
        "ligand_id": ligand.ligand_id,
        "smiles": ligand.smiles,
        "conditions": conditions,
        "poses": rows,
        "execution": execution,
    }
    fd, temporary = tempfile.mkstemp(dir=path.parent, prefix=".pose-", suffix=".part")
    try:
        with os.fdopen(fd, "w") as handle:
            handle.write(json.dumps(payload, indent=1))
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def _json_rows(docked: list[DockPose]) -> list[dict[str, Any]]:
    rows = []
    for pose in docked:
        row = asdict(pose)
        row["sdf_path"] = str(pose.sdf_path) if pose.sdf_path else None
        rows.append(row)
    return rows


def _dock_task(payload: dict[str, Any]) -> tuple[str, list[dict[str, Any]], Optional[str]]:
    """Process-pool entry point. Returns JSON-safe rows so nothing exotic is pickled."""
    logging.basicConfig(level=logging.WARNING)
    receptor: Receptor = payload["receptor"]
    ligand: LigandRecord = payload["ligand"]
    workdir = Path(payload["workdir"])
    execution: dict[str, Any] = {}
    def record(error: Optional[str] = None) -> None:
        execution.setdefault('docking_conditions', payload['conditions'])
        if payload.get("requested_engine") == "hybrid":
            execution["requested_engine"] = "hybrid"
            execution["route_reason"] = payload.get("route_reason", "hybrid_cpu_route")
            execution["attempts"] = payload.get("prior_attempts", []) + execution.get("attempts", [])
            if payload.get("prior_gpu_batch_context") is not None:
                execution["prior_gpu_batch_context"] = payload["prior_gpu_batch_context"]
        if error:
            execution.update(state="engine_error", error=error)
        if payload.get("execution_receipt"):
            path = Path(payload["execution_receipt"])
            path.parent.mkdir(parents=True, exist_ok=True)
            content = {"ligand_id": ligand.ligand_id, "smiles": ligand.smiles,
                       "conditions": payload["conditions"], "execution": execution}
            fd, temporary = tempfile.mkstemp(dir=path.parent, prefix=".execution-", suffix=".part")
            try:
                with os.fdopen(fd, "w") as handle:
                    json.dump(content, handle, indent=1)
                os.replace(temporary, path)
            finally:
                if os.path.exists(temporary):
                    os.unlink(temporary)
    try:
        docked = dock_one(
            receptor,
            ligand,
            workdir,
            exhaustiveness=payload["exhaustiveness"],
            n_poses=payload["n_poses"],
            seed=payload["seed"],
            engine_name=payload.get("engine_name", "vina"),
            unidock_search_mode=payload.get("unidock_search_mode", "balance"),
            unidock_bin=payload.get("unidock_bin"),
            execution=execution,
        )
    except DockingEngineError as exc:



        log.warning("%s: %s; leaving it uncached so a later run retries",
                    ligand.ligand_id, exc)
        record(str(exc))
        return ligand.ligand_id, [], str(exc)
    except Exception as exc:  # a screen must survive any single-ligand blow-up
        log.exception("%s: docking raised", ligand.ligand_id)
        record(f"{type(exc).__name__}: {exc}")
        return ligand.ligand_id, [], f"{type(exc).__name__}: {exc}"
    record()
    _save_cached(workdir, ligand, docked, payload["conditions"], execution=execution)
    return ligand.ligand_id, _json_rows(docked), None


def _prep_one(args: "tuple[str, str, str]") -> "tuple[str, Any]":
    """Pool entry point: ``(ligand_id, smiles, workdir)`` -> path or exception.

    The exception travels back as a value rather than propagating, so one
    unpreparable ligand cannot end the pool and the caller keeps the same
    per-ligand fallback it had when this ran serially.
    """
    from .ligprep import prepare_ligand_pdbqt

    ligand_id, smiles, workdir = args
    try:
        return ligand_id, prepare_ligand_pdbqt(smiles, ligand_id, Path(workdir))
    except Exception as exc:                      # noqa: BLE001 - per ligand
        return ligand_id, exc


def _prepare_many(ligands: "list[LigandRecord]", workdir: Path,
                  workers: Optional[int] = None) -> "dict[str, Any]":
    """Prepare every ligand, in parallel when there are enough to be worth it.

    Falls back to the serial path for a handful of ligands, where the pool
    costs more than it saves, and whenever a pool cannot be started at all --
    inside an existing worker, say. A preparation failure is returned as the
    exception, never raised.
    """
    tasks = [(ligand_record.ligand_id, ligand_record.smiles, str(workdir)) for ligand_record in ligands]
    if len(tasks) < 8:
        return dict(_prep_one(t) for t in tasks)
    n = workers or max(1, min(len(tasks), (os.cpu_count() or 2) - 1))
    if multiprocessing.parent_process() is not None:
        n = 1
    if n == 1:
        return dict(_prep_one(t) for t in tasks)
    out: dict[str, Any] = {}
    try:
        with ProcessPoolExecutor(max_workers=n, mp_context=_mp_context()) as pool:
            for ligand_id, value in pool.map(_prep_one, tasks, chunksize=8):
                out[ligand_id] = value
    except Exception as exc:                      # noqa: BLE001 - pool itself
        log.warning("parallel ligand preparation unavailable (%s); preparing "
                    "%d ligand(s) serially", exc, len(tasks))
        return dict(_prep_one(t) for t in tasks)
    return out


def _hybrid_preparation_plan(
    library: list[LigandRecord], workdir: Path, conditions: dict[str, Any],
) -> dict[str, Any]:
    """Prepare and partition the full library independently of cache hits.

    A pending-only partition changes GPU co-members on partial resumption.
    The full eligible prepared order and exact bytes define fixed chunks;
    cached GPU receipts must agree with their current chunk before reuse.
    """
    from . import unidock as ud
    prepared: dict[str, Path] = {}
    eligible: list[LigandRecord] = []
    cpu: list[LigandRecord] = []
    reasons: dict[str, str] = {}
    outcomes = _prepare_many(library, workdir)
    for ligand in library:
        outcome = outcomes.get(ligand.ligand_id)
        if isinstance(outcome, BaseException):
            reasons[ligand.ligand_id] = "ligand_preparation_error_cpu_retry"
            cpu.append(ligand)
            continue
        if outcome is None:
            reasons[ligand.ligand_id] = "ligand_preparation_rejected_cpu_retry"
            cpu.append(ligand)
            continue
        try:
            path = Path(outcome)
            text = path.read_text(encoding="utf-8")
            artifacts.sha256(path)
        except (OSError, ValueError, TypeError) as exc:
            log.warning("%s: prepared ligand cannot be read (%s); retrying on CPU",
                        ligand.ligand_id, exc)
            reasons[ligand.ligand_id] = "ligand_preparation_error_cpu_retry"
            cpu.append(ligand)
            continue
        prepared[ligand.ligand_id] = path
        if ud.has_glue_atoms(text):
            reasons[ligand.ligand_id] = "macrocycle_glue_cpu_route"
            cpu.append(ligand)
        else:
            reasons[ligand.ligand_id] = "hybrid_cpu_route"
            eligible.append(ligand)
    contexts: dict[str, dict[str, Any]] = {}
    chunks: list[list[LigandRecord]] = []
    size = max(1, int(conditions.get("batch_size", ud.BATCH_SIZE)))
    for start in range(0, len(eligible), size):
        members = eligible[start:start + size]
        chunks.append(members)
        inputs = [{"ligand_id": item.ligand_id,
                   "prepared_ligand_sha256": artifacts.sha256(prepared[item.ligand_id])}
                  for item in members]
        context = {"batch_index": start // size, "members": inputs,
                   "batch_index_scope": "full_eligible_prepared_library",
                   "batch_policy": "fixed_full_eligible_library_v1",
                   "execution_mode": "gpu_batch", "batch_size": size,
                   "ordered_library_sha256": conditions.get("ordered_library_sha256")}
        for item in members:
            contexts[item.ligand_id] = context
    return {"prepared": prepared, "eligible": eligible, "cpu": cpu,
            "route_reasons": reasons, "batch_contexts": contexts, "chunks": chunks}


def _hybrid_gpu_pass(
    receptor: Receptor,
    pending: "list[LigandRecord]",
    workdir: Path,
    *,
    n_poses: int,
    seed: int,
    search_mode: str,
    binary: Optional[str],
    conditions: dict[str, Any],
    names: dict[str, str],
    executions: Optional[dict[str, dict[str, Any]]] = None,
    preparation_plan: Optional[dict[str, Any]] = None,
) -> "tuple[list[dict[str, Any]], list[LigandRecord]]":
    'Dock what the GPU is trusted for. Returns (rows, what is left for CPU).'
    from . import unidock as ud
    traces = executions if executions is not None else {}
    for ligand in pending:
        traces[ligand.ligand_id] = {"requested_engine": "hybrid", "actual_engine": None,
                                  "route_reason": "hybrid_cpu_route", "attempts": []}

    try:
        ud.executable(binary)
    except ud.UnidockUnavailable as exc:
        for ligand in pending:
            traces[ligand.ligand_id].update(route_reason="unidock_unavailable_cpu_fallback")
        log.info("hybrid docking: no GPU engine (%s); docking all %d ligand(s) "
                 "on CPU Vina", exc, len(pending))
        return [], pending

    plan = preparation_plan if preparation_plan is not None else _hybrid_preparation_plan(
        pending, workdir, conditions)
    prepared = plan["prepared"]
    pending_ids = {ligand.ligand_id for ligand in pending}
    small = [ligand for ligand in plan["eligible"] if ligand.ligand_id in pending_ids]
    large = [ligand for ligand in plan["cpu"] if ligand.ligand_id in pending_ids]
    for ligand in pending:
        traces[ligand.ligand_id]["route_reason"] = plan["route_reasons"][ligand.ligand_id]

    if not small:
        log.info("hybrid docking: %d pending ligand(s) have a CPU route "
                 "(macrocycle glue or preparation retry); none go to GPU", len(pending))
        return [], large

    log.info("hybrid docking: %d ligand(s) go to the GPU in batches of %d; %d "
             "carry macrocycle glue atoms the GPU engine does not implement and "
             "go to CPU Vina", len(small), ud.BATCH_SIZE, len(large))
    touched = [chunk for chunk in plan["chunks"]
               if any(item.ligand_id in pending_ids for item in chunk)]
    # Include completed co-members as execution context. Only pending results
    # are published; valid completed caches and their pose files are retained.
    gpu_inputs = [item for chunk in touched for item in chunk]
    log.info("hybrid docking: executing %d fixed GPU chunk(s), %d input(s), "
             "%d pending result(s)", len(touched), len(gpu_inputs), len(small))
    for ligand in small:
        traces[ligand.ligand_id]["batch_context"] = plan["batch_contexts"][ligand.ligand_id]
    try:
        docked = ud.dock_batch(
            receptor, [prepared[lig.ligand_id] for lig in gpu_inputs],
            search_mode=search_mode, n_poses=n_poses, seed=seed, binary=binary,
            batch_size=max(1, int(conditions.get("batch_size", ud.BATCH_SIZE))))
    except (ud.UnidockUnavailable, FileNotFoundError, ValueError) as exc:
        for ligand in small:
            traces[ligand.ligand_id].update(route_reason="unidock_batch_error_cpu_fallback",
                attempts=[{"engine": "unidock", "state": "error", "error_type": type(exc).__name__,
                           "error": str(exc), "engine_identity": conditions.get("software_identity"),
                           "batch_context": traces[ligand.ligand_id]["batch_context"]}])
        log.warning("hybrid docking: the GPU batch failed (%s); all %d ligand(s) "
                    "fall back to CPU Vina", exc, len(pending))
        return [], large + small

    rows: list[dict[str, Any]] = []
    rejected = 0
    left: list[LigandRecord] = list(large)
    for ligand in small:
        text = docked.get(Path(prepared[ligand.ligand_id]).stem)
        affinity = ud.best_affinity(text) if text else None
        trace = traces[ligand.ligand_id]
        trace["attempts"] = [{"engine": "unidock", "state": "completed" if text else "missing_output",
                              "best_affinity": affinity, "engine_identity": conditions.get("software_identity"),
                              "batch_context": trace["batch_context"]}]
        if text:
            _record_raw_output(workdir, text, trace["attempts"][-1])
        if affinity is None or not ud.is_plausible(affinity):
            trace["route_reason"] = "unidock_score_window_cpu_fallback" if affinity is not None else "unidock_missing_pose_cpu_fallback"
            if affinity is not None:
                rejected += 1
                log.warning("%s: unidock returned %.3f kcal/mol, outside the "
                            "operational score window; re-docking on CPU Vina",
                            ligand.ligand_id, affinity)
            left.append(ligand)
            continue
        try:
            poses = _poses_from_docked_text(text, ligand, workdir, receptor)
        except Exception as exc:                  # noqa: BLE001 - per ligand
            trace["route_reason"] = "unidock_pose_conversion_error_cpu_fallback"
            trace["attempts"][-1].update(state="pose_conversion_error", error=str(exc))
            log.warning("%s: could not read the GPU pose (%s); re-docking on "
                        "CPU Vina", ligand.ligand_id, exc)
            left.append(ligand)
            continue
        if not poses:
            trace["route_reason"] = "unidock_no_usable_pose_cpu_fallback"
            left.append(ligand)
            continue
        trace.update(actual_engine="unidock", state="completed", route_reason="hybrid_gpu_batch",
                     engine_identity=conditions.get("software_identity"), seed=int(seed),
                     n_poses_requested=int(n_poses), n_poses_returned=len(poses),
                     unidock_search_mode=search_mode,
                     prepared_ligand_sha256=artifacts.sha256(prepared[ligand.ligand_id]),
                     prepared_ligand_artifact=artifacts.describe(workdir, prepared[ligand.ligand_id]),
                     receptor_pdbqt_sha256=artifacts.sha256(Path(receptor.pdbqt)))
        _record_preparation_provenance(trace, workdir, prepared[ligand.ligand_id], ligand)
        _save_cached(workdir, ligand, poses, conditions, execution=trace)
        rows.extend(_pose_rows(poses, names))
    log.info("hybrid docking: %d ligand(s) scored on the GPU, %d rejected by the "
             "plausibility window, %d for CPU Vina",
             len(rows and {r["ligand_id"] for r in rows} or set()), rejected,
             len(left))
    return rows, left


def _unidock_cache_unsafe(workdir: Path, ligand: LigandRecord,
                          poses: "list[DockPose]") -> Optional[str]:
    'Why this cached Uni-Dock result must not be served to a hybrid run.'
    from . import unidock as ud
    from .ligprep import prepare_ligand_pdbqt

    # No pose at all is the case hybrid exists for. A fresh GPU run that comes
    # back empty falls through to CPU Vina; a *cached* empty result from an
    # earlier `--engine unidock` run was accepted as a finished ligand, so the
    # same library under the same hybrid settings routed a different set to the
    # CPU depending on what the workdir already held. Scoped to the alternate
    # GPU cache on purpose: an empty Vina result under the primary conditions
    # still means "Vina looked and found nothing", which is an answer.
    if not poses:
        return ("the cached GPU run produced no pose, and a hybrid run retries "
                "that on CPU Vina rather than recording it as docked")
    affinities = [p.affinity for p in poses
                  if isinstance(p.affinity, (int, float)) and not isinstance(p.affinity, bool)
                  and math.isfinite(p.affinity)]
    if len(affinities) != len(poses):
        return "the cached GPU pose lacks a finite numeric affinity; CPU retry is required"
    if affinities:
        best = min(affinities)
        if not ud.is_plausible(best):
            return (f"cached GPU affinity {best:.3f} kcal/mol is outside "
                    f"the operational score window {ud.PLAUSIBLE_AFFINITY}")
    try:
        path = prepare_ligand_pdbqt(ligand.smiles, ligand.ligand_id, workdir)
    except Exception as exc:                      # noqa: BLE001 - per ligand
        return f"could not re-derive the prepared ligand to check it ({exc})"
    if path is None:
        return "the prepared ligand is no longer available to check"
    if ud.has_glue_atoms(Path(path).read_text(encoding="utf-8")):
        return ("the ligand carries macrocycle glue atoms, which this GPU "
                "engine does not implement; it should never have been scored "
                "there")
    return None


def dock_many(
    receptor: Receptor,
    ligands: list[LigandRecord],
    workdir: Path,
    exhaustiveness: int = 16,
    n_poses: int = 9,
    seed: int = 42,
    n_workers: Optional[int] = None,
    progress: bool = True,
    engine: str = "vina",
    unidock_search_mode: str = "balance",
    unidock_bin: Optional[str] = None,
) -> "pd.DataFrame":
    """Dock a set of ligands in parallel, resumably. One row per pose.

    Resumption is per ligand: a ligand is skipped only when
    ``results/<id>.json`` exists, its pose SDFs are still on disk, and it was
    docked under *these* conditions -- this receptor file, chain and box, and
    this ``exhaustiveness``/``n_poses``/``seed``. Change the receptor or any
    search setting and every ligand is docked again, each one logging what
    changed. Delete the JSON to force a redock regardless.

    This matters because the CLI reuses one working directory per marker: a
    second run with ``--pdb`` pointing somewhere else lands in the same
    ``dock/`` directory as the first.
    """
    import pandas as pd

    from .ligprep import valid_ligand_id
    if any(not valid_ligand_id(ligand.ligand_id) for ligand in ligands):
        raise ValueError('ligand id must be a filename component')
    if len({ligand.ligand_id for ligand in ligands}) != len(ligands):
        raise ValueError('ligand ids must be unique within a docking library')
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    names = {lig.ligand_id: (lig.name or lig.ligand_id) for lig in ligands}
    # Hashed once for the whole screen, not once per ligand: the receptor PDBQT
    # is the same file for every one of them.
    if engine not in ("vina", "unidock", "hybrid"):
        raise ValueError(f"engine must be vina, unidock or hybrid, got {engine!r}")
    software.refresh()
    conditions = cachekey.dock_conditions(receptor, exhaustiveness, n_poses, seed,
                                          engine="vina" if engine == "hybrid"
                                          else engine,
                                          search_mode=unidock_search_mode,
                                          unidock_bin=unidock_bin)
    # A hybrid run legitimately holds both engines' records side by side, so a
    # cached ligand counts as done under either. Each record still names the
    # engine that produced it -- that is what the key is for -- so nothing is
    # served one engine's number under the other's name.
    alt_conditions = (cachekey.unidock_batch_conditions(receptor, ligands, exhaustiveness,
                                               n_poses, seed,
                                               search_mode=unidock_search_mode,
                                               unidock_bin=unidock_bin)
                      if engine == "hybrid" else None)
    hybrid_plan = (_hybrid_preparation_plan(ligands, workdir, alt_conditions)
                   if alt_conditions is not None else None)

    rows: list[dict[str, Any]] = []
    executions: dict[str, dict[str, Any]] = {}
    pending: list[LigandRecord] = []
    for ligand in ligands:







        trace: dict[str, Any] = {}
        cached = _load_cached(workdir, ligand, conditions,
                             quiet=alt_conditions is not None, execution=trace)
        if cached:
            # Cold and warm have to apply one policy or the cache decides what
            # the safety check was -- the lesson `_unidock_cache_unsafe`
            # already records, which guarded only the *alternate*-conditions
            # probe. RARG's +8.93 kcal/mol lives under the primary one.
            kept = _keep_plausible(cached, ligand.ligand_id, "the cache")
            if len(kept) != len(cached):
                log.warning("%s: the cached result is outside the operational "
                            "score window; re-docking", ligand.ligand_id)
                cached = None
        if cached is None and alt_conditions is not None:
            trace.clear()
            alt = _load_cached(workdir, ligand, alt_conditions, execution=trace)
            # Accepting it here would let a stored result skip the gates a
            # fresh GPU result has to clear. Re-apply them.
            if alt is not None:
                why = _unidock_cache_unsafe(workdir, ligand, alt)
                expected_batch = hybrid_plan["batch_contexts"].get(ligand.ligand_id)
                if not why and (expected_batch is None or trace.get("batch_context") != expected_batch):
                    why = "the cached GPU batch differs from the current full prepared library chunk"
                if not why and not any(
                    attempt.get("engine") == "unidock" and attempt.get("state") == "completed"
                    and attempt.get("batch_context") == expected_batch
                    and attempt.get("engine_identity") == alt_conditions.get("software_identity")
                    for attempt in trace.get("attempts", []) if isinstance(attempt, dict)
                ):
                    why = "the cached completed GPU attempt does not bind the current full chunk"
                if why:
                    log.warning("%s: discarding the cached GPU result -- %s; "
                                "re-docking", ligand.ligand_id, why)
                    alt = None
            cached = alt
        if cached is None:
            pending.append(ligand)
        else:
            trace["reuse_requested_engine"] = engine
            executions[ligand.ligand_id] = trace
            rows.extend(_pose_rows(cached, names))
    if rows:
        log.info("resuming: %d ligand(s) already done, %d to dock", len(ligands) - len(pending), len(pending))

    if engine == "hybrid" and pending:
        done, pending = _hybrid_gpu_pass(
            receptor, pending, workdir, n_poses=n_poses, seed=seed,
            search_mode=unidock_search_mode, binary=unidock_bin,
            conditions=alt_conditions or {}, names=names, executions=executions,
            preparation_plan=hybrid_plan)
        rows.extend(done)

    if pending:
        workers = n_workers if n_workers is not None else max(1, (os.cpu_count() or 2) - 1)
        workers = max(1, min(int(workers), len(pending)))
        payloads = [
            {
                "receptor": receptor,
                "ligand": ligand,
                "workdir": str(workdir),
                "exhaustiveness": exhaustiveness,
                "n_poses": n_poses,
                "seed": seed,
                "conditions": conditions,








                "engine_name": "vina" if engine == "hybrid" else engine,
                "unidock_search_mode": unidock_search_mode,
                "unidock_bin": unidock_bin,
                "requested_engine": engine,
                "route_reason": executions.get(ligand.ligand_id, {}).get("route_reason", "selected_engine"),
                "prior_attempts": executions.get(ligand.ligand_id, {}).get("attempts", []),
                "prior_gpu_batch_context": executions.get(ligand.ligand_id, {}).get("batch_context"),
                "execution_receipt": str(workdir / "executions" / (uuid.uuid4().hex + ".json")),
            }
            for ligand in pending
        ]
        done = 0
        if workers == 1:
            for payload in payloads:
                ligand_id, pose_rows, error = _dock_task(payload)
                done += 1
                rows.extend(_annotate(pose_rows, names))
                _report(progress, done, len(pending), ligand_id, pose_rows, error)
        else:
            with ProcessPoolExecutor(max_workers=workers, mp_context=_mp_context()) as pool:
                futures = {pool.submit(_dock_task, p): p["ligand"].ligand_id for p in payloads}
                for future in as_completed(futures):
                    ligand_id, pose_rows, error = future.result()
                    done += 1
                    rows.extend(_annotate(pose_rows, names))
                    _report(progress, done, len(pending), ligand_id, pose_rows, error)

        for payload in payloads:
            path = Path(payload["execution_receipt"])
            if path.is_file():
                receipt = json.loads(path.read_text())
                if (receipt.get("ligand_id") != payload["ligand"].ligand_id or
                        receipt.get("smiles") != payload["ligand"].smiles or
                        receipt.get("conditions") != payload["conditions"]):
                    raise ValueError("docking execution receipt identity mismatch")
                executions[receipt["ligand_id"]] = receipt["execution"]

    columns = ["ligand_id", "name", "pose_index", "affinity", "ligand_efficiency", "sdf_path"]
    frame = pd.DataFrame(rows, columns=columns)
    if rows:
        frame = frame.sort_values(["affinity", "ligand_id", "pose_index"]).reset_index(drop=True)
    frame.attrs["execution_by_ligand"] = executions
    return frame


def _mp_context() -> Optional["multiprocessing.context.BaseContext"]:
    """Prefer ``fork`` on POSIX.

    The default start method on macOS and Windows is ``spawn``, which re-imports
    the parent's ``__main__`` module in every worker. That breaks outright when
    the caller is a REPL, a notebook or a ``python - <<EOF`` heredoc (there is no
    file to import), and silently re-runs module-level code when it is a script
    without an ``if __name__ == "__main__"`` guard. ``fork`` inherits the parent
    instead, so a screen launched interactively just works. Workers run Vina with
    ``cpu=1`` and never touch parent threads, which is the condition under which
    forking a multithreaded process is safe.
    """
    try:
        return multiprocessing.get_context("fork")
    except ValueError:  # Windows: no fork, spawn is the only option
        log.debug("fork start method unavailable; using the platform default")
        return None


def _pose_rows(docked: list[DockPose], names: dict[str, str]) -> list[dict[str, Any]]:
    return _annotate(_json_rows(docked), names)


def _annotate(rows: list[dict[str, Any]], names: dict[str, str]) -> list[dict[str, Any]]:
    return [{**row, "name": names.get(row["ligand_id"], row["ligand_id"])} for row in rows]


def _report(
    progress: bool,
    done: int,
    total: int,
    ligand_id: str,
    rows: list[dict[str, Any]],
    error: Optional[str],
) -> None:
    if error:
        log.error("[%d/%d] %s FAILED: %s", done, total, ligand_id, error)
    elif not rows:
        log.warning("[%d/%d] %s produced no pose", done, total, ligand_id)
    elif progress:
        best = min(r["affinity"] for r in rows)
        log.info("[%d/%d] %s best %.2f kcal/mol", done, total, ligand_id, best)
