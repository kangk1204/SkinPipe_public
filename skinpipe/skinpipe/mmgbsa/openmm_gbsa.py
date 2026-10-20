'Single-trajectory MM/GBSA rescoring of docked poses.'
from __future__ import annotations

import logging
import multiprocessing
import os
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Optional

from ..types import GBSAResult
from . import cachekey, metal_guard
from . import cofactor_guard
from .cofactor_guard import CofactorSiteError  # noqa: F401 - re-exported
from .metal_guard import MetalSiteError  # noqa: F401 - re-exported
from .cachekey import cache_path as _cache_path
from .engine import (CONVERGENCE_RMS_FORCE, ENERGY_KEYS, Evaluator,
                     cpu_threads, execution_conditions, make_platform,
                     md_snapshots, mean_std, minimise)
from .ligand import LigandParameterisationError, parameterise_ligand  # noqa: F401
from .records import DATAFRAME_COLUMNS as _DATAFRAME_COLUMNS  # noqa: F401 - re-exported
from .records import _non_finite_energies, to_frame, to_result
from .system import (ACE_SURFACE_TENSION_KCAL_PER_A2, DEFAULT_PREPARATION_SEED,
                     ComplexSystem, SystemBuildError, build_complex_system,
                     prepare_receptor, validate_preparation_seed)

LOGGER = logging.getLogger(__name__)

__all__ = ["parameterise_ligand", "build_complex_system", "mmgbsa_single",
           "mmgbsa_batch", "ComplexSystem", "GBSAResult", "MetalSiteError",
           "CofactorSiteError"]


def _flag_non_finite_energies(record: dict[str, Any]) -> dict[str, Any]:
    """A converged record whose energies are not finite, demoted to a failure.

    Returns a new record; the argument is not mutated.

    ``mean_std`` averages whatever the frames produced, and a frame that blew
    up produces ``+/-Inf``: ``mean_std([-45.0, -44.0, inf])`` is ``(inf, nan)``.
    Reporting that as ``converged`` is a claim about a number that does not
    exist.  It is also expensive: :func:`skinpipe.score.score_screen` rejects a
    non-finite ``dG_bind`` on a converged result by raising, so the whole
    screen died *after* docking and MM/GBSA had been paid for (D-11).

    Done here, at the source, rather than only in
    :func:`skinpipe.mmgbsa.records.to_result`, so the batch DataFrame and the
    JSON written under ``workdir/cache/`` say "failed" too -- otherwise the
    cache file keeps ``"converged": true`` with ``"dG_bind": Infinity`` and
    resumes it on the next run.
    """
    corrupt = _non_finite_energies(record)
    if not corrupt:
        return record
    reported = ", ".join(f"{name}={record[name]!r}" for name in corrupt)
    LOGGER.error("MM/GBSA produced a non-finite energy for %s: %s",
                 record.get("ligand_id"), reported)
    note = (f"FAILED: non-finite MM/GBSA energy ({reported}); there is no "
            f"number to rank this pose by, so it is reported as a failed pose "
            f"rather than a converged one")
    return {**record, "converged": False,
            "notes": "; ".join(filter(None, [note, record.get("notes") or ""]))}


def _run_single(receptor_pdb: Path, ligand_sdf: Path, workdir: Path,
                minimise_steps: int = 500, implicit_solvent: str = "GBn2",
                ligand_forcefield: str = "gaff-2.11", n_snapshots: int = 1,
                md_ps: float = 0.0, salt_molar: float = 0.15,
                ligand_id: Optional[str] = None, pose_index: int = 0,
                charge_method: str = "auto", restraint_k: float = 5.0,
                temperature_k: float = 300.0, platform: Optional[str] = None,
                threads: Optional[int] = None, prepare_receptor_pdb: bool = True,
                nonbonded_cutoff_nm: Optional[float] = None,
                constrain_h_bonds: Optional[bool] = None, seed: int = 0,
                on_metal_site: str = metal_guard.REFUSE,
                preparation_seed: int = DEFAULT_PREPARATION_SEED,
                cpu_profile: str = "auto") -> dict[str, Any]:
    """Worker body: returns a plain dict so it can be JSON-cached and pickled."""
    cpu_threads(threads, cpu_profile)
    started = time.perf_counter()
    ligand_sdf = Path(ligand_sdf)
    record: dict[str, Any] = {
        "ligand_id": ligand_id or ligand_sdf.stem,
        "pose_index": int(pose_index),
        "converged": False,
        "notes": "",
        "sdf_path": str(ligand_sdf),
    }
    try:
        workdir = Path(workdir)
        workdir.mkdir(parents=True, exist_ok=True)
        cs = build_complex_system(
            receptor_pdb, ligand_sdf, implicit_solvent=implicit_solvent,
            ligand_forcefield=ligand_forcefield, charge_method=charge_method,
            salt_molar=salt_molar, temperature_k=temperature_k,
            restraint_k=restraint_k, prepare=prepare_receptor_pdb,
            nonbonded_cutoff_nm=nonbonded_cutoff_nm,
            # Off by default: HBonds constraints buy a 2 fs MD timestep but
            # make LocalEnergyMinimizer re-run L-BFGS until constraints are
            # satisfied, which costs ~19x on CPU.  Without them the MD
            # timestep drops to 1 fs, which is the cheaper trade here.
            constrain_h_bonds=bool(constrain_h_bonds),
            on_metal_site=on_metal_site, workdir=workdir,
            preparation_seed=preparation_seed)
        if cs.parameterization is not None:
            record["parameterization"] = cs.parameterization

        plat, properties = make_platform(platform, threads, cpu_profile)
        evaluator = Evaluator(cs, plat, properties)
        record["execution"] = execution_conditions(
            plat, properties, evaluator.complex_context)
        evaluator.complex_context.setPositions(cs.positions)
        diag = minimise(evaluator.complex_context, max(1, int(minimise_steps)))
        converged = diag["rms_force"] <= CONVERGENCE_RMS_FORCE

        if md_ps > 0 and n_snapshots >= 1:
            frames = md_snapshots(cs, evaluator, plat, properties, md_ps,
                                  max(1, int(n_snapshots)), temperature_k, seed)
            sampling = f"md {md_ps:g}ps/{len(frames)}frames"
        else:
            minimised = evaluator.complex_context.getState(
                getPositions=True).getPositions()
            frames = [evaluator.decompose(minimised)]
            sampling = "single-point on minimised pose"

        for key in ENERGY_KEYS:
            mean, std = mean_std([f[key] for f in frames])
            record[key] = mean
            if key == "dG_bind":
                record["dG_std"] = std
        record["n_frames"] = len(frames)
        record["converged"] = bool(converged)
        record["rms_force_kj_nm"] = diag["rms_force"]
        record["max_force_kj_nm"] = diag["max_force"]
        record["e_min_before_kcal"] = diag["e_before_kcal"]
        record["e_min_after_kcal"] = diag["e_after_kcal"]
        record["n_atoms_complex"] = cs.n_atoms
        record["n_atoms_ligand"] = cs.ligand.n_atoms

        notes = [sampling, f"minimise {minimise_steps} steps",
                 f"rms|F|={diag['rms_force']:.1f} kJ/mol/nm"
                 + ("" if converged else f" > {CONVERGENCE_RMS_FORCE:.0f} NOT CONVERGED"),
                 f"nonpolar={record['e_nonpolar']:.2f} kcal/mol "
                 f"(gamma={ACE_SURFACE_TENSION_KCAL_PER_A2} kcal/mol/A^2)"]
        if record.get("dG_std") is not None:
            notes.append(f"sd={record['dG_std']:.2f} kcal/mol over {len(frames)} frames")
        notes.extend(cs.provenance)
        execution = record["execution"]
        notes.append(f"OpenMM {execution['openmm_version']} "
                     f"platform={execution['resolved_platform']}")
        if execution["resolved_cpu_profile"] is not None:
            count = execution["resolved_platform_properties"]["Threads"]
            notes.append(f"CPU profile={execution['resolved_cpu_profile']} threads={count}")
            if int(count) > 1:
                notes.append("parallel CPU numerical repeats may differ")
        record["notes"] = "; ".join(notes)
    except (metal_guard.MetalSiteError, cofactor_guard.CofactorSiteError):
        raise          # not a failed pose: the receptor is the wrong input
    except (LigandParameterisationError, SystemBuildError) as exc:
        LOGGER.warning("MM/GBSA failed for %s: %s", record["ligand_id"], exc)
        record["notes"] = f"FAILED: {exc}"
    except Exception as exc:  # noqa: BLE001 - one bad pose must not kill the batch
        LOGGER.error("MM/GBSA crashed for %s: %s", record["ligand_id"], exc)
        LOGGER.debug("%s", traceback.format_exc())
        record["notes"] = f"FAILED: {type(exc).__name__}: {exc}"
    record["wall_s"] = time.perf_counter() - started
    # A pose whose energy came back +/-Inf computed nothing, whatever the
    # minimiser's rms|F| said about it; see `_flag_non_finite_energies`.
    return _flag_non_finite_energies(record)


def mmgbsa_single(receptor_pdb: Path, ligand_sdf: Path, workdir: Path,
                  minimise_steps: int = 500, implicit_solvent: str = "GBn2",
                  ligand_forcefield: str = "gaff-2.11", n_snapshots: int = 1,
                  md_ps: float = 0.0, salt_molar: float = 0.15,
                  **kwargs: Any) -> GBSAResult:
    """Rescore one docked pose by single-trajectory MM/GBSA.

    Read the module docstring before interpreting the number: it is a
    rank-ordering score for a congeneric series, not a binding free energy.

    Parameters
    ----------
    receptor_pdb:
        Protein structure.  Cleaned and protonated with PDBFixer unless
        ``prepare_receptor_pdb=False``.
    ligand_sdf:
        Single-molecule SDF holding the pose, with 3D coordinates.
    workdir:
        Scratch directory for the prepared receptor.
    minimise_steps:
        L-BFGS iterations of restrained minimisation.  500 relaxes docking
        clashes; it is not a converged minimum, and ``converged`` says so.
    implicit_solvent:
        ``GBn2`` (default), ``OBC2``, ``OBC1`` or ``HCT``.
    ligand_forcefield:
        ``gaff-2.11`` (needs AmberTools) or an OpenFF line such as
        ``openff-2.2.0``.  If GAFF is requested without AmberTools the run
        falls back to OpenFF and says so in ``notes``.
    n_snapshots, md_ps:
        With ``md_ps > 0``, run that much Langevin dynamics at
        ``temperature_k`` and average over ``n_snapshots`` evenly spaced
        frames; the standard deviation is reported in ``notes``.  Pass a
        nonzero ``seed`` for a reproducible trajectory (OpenMM reads seed 0
        as "choose randomly").
    salt_molar:
        Ionic strength for Debye-Huckel screening of the GB term.

    Other Parameters
    ----------------
    ligand_id, pose_index, charge_method, restraint_k, temperature_k,
    platform, threads, prepare_receptor_pdb, nonbonded_cutoff_nm,
    constrain_h_bonds, seed, on_metal_site.

    ``nonbonded_cutoff_nm=None`` (default) is exact and O(N^2); set ~1.8 to
    trade some accuracy for a large speed-up on CPU.

    Returns
    -------
    GBSAResult
        On any failure -- unparameterisable ligand, missing template, charge
        timeout -- returns ``converged=False`` with the reason in ``notes``
        instead of raising.

    Raises
    ------
    ~skinpipe.mmgbsa.metal_guard.MetalSiteError
        Not folded into ``notes``: the pose binds a catalytic metal this force
        field cannot represent, so there is no number to report rather than a
        bad one.  A caller who wants the apo number anyway passes
        ``on_metal_site="score_without_metal"``, and every result then carries
        ``SCORED WITHOUT METAL`` in its ``notes``.
    ~skinpipe.mmgbsa.cofactor_guard.CofactorSiteError
        Also not folded into ``notes``, and for the same reason one level
        along: a structural cofactor lines this pocket and is absent from the
        receptor, so the rescore would be against a site the pose was not
        docked into.
    """




    metal_guard.check_poses(Path(receptor_pdb), [Path(ligand_sdf)],
                            kwargs.get("on_metal_site", metal_guard.REFUSE))
    cofactor_guard.refuse_poses(Path(receptor_pdb), [Path(ligand_sdf)])
    record = _run_single(
        Path(receptor_pdb), Path(ligand_sdf), Path(workdir),
        minimise_steps=minimise_steps, implicit_solvent=implicit_solvent,
        ligand_forcefield=ligand_forcefield, n_snapshots=n_snapshots,
        md_ps=md_ps, salt_molar=salt_molar, **kwargs)
    return to_result(record)


@dataclass(frozen=True)
class _Job:
    "One pose's work order, handed to :func:`_worker` across a process."

    receptor: str
    ligand_id: str
    pose_index: int
    sdf: str
    workdir: str
    kwargs: dict[str, Any]
    conditions: dict[str, Any]
    #: Resolved before cache lookup; restarting cannot redistribute threads.
    threads: int = 1

    def __post_init__(self) -> None:
        if not isinstance(self.conditions, dict):
            raise TypeError(
                f"conditions must be the mapping from "
                f"cachekey.protocol_conditions, not "
                f"{type(self.conditions).__name__} ({self.conditions!r}); a "
                f"record cached under a non-mapping can never be verified and "
                f"every later run recomputes it")
        if isinstance(self.threads, bool) or not isinstance(self.threads, int):
            raise TypeError(
                f"threads must be an OpenMM thread count, not "
                f"{type(self.threads).__name__} ({self.threads!r})")
        if self.threads < 1:
            raise ValueError("threads must be a positive integer")


def _worker(job: _Job) -> dict[str, Any]:
    kwargs = dict(job.kwargs)
    if kwargs.get("threads") is None:
        kwargs["threads"] = job.threads
    record = _run_single(Path(job.receptor), Path(job.sdf), Path(job.workdir),
                         ligand_id=job.ligand_id, pose_index=job.pose_index,
                         **kwargs)
    # The parent computed these, so what is written is what it looked for.
    cachekey.save_record(
        _cache_path(Path(job.workdir), job.ligand_id, job.pose_index),
        record, job.conditions)
    return record


def _worker_failure(job: _Job, exc: Exception, *, resume: bool) -> dict[str, Any]:
    """Account for this exact pose without writing a fabricated cache result.

    A worker may have saved its completed result before its process/IPC died.
    Resume keeps that condition-verified record, with original cache bytes and
    protocol untouched. Explicit resume=False never falls back to old records.
    """
    if resume:
        saved = cachekey.load_record(
            _cache_path(Path(job.workdir), job.ligand_id, job.pose_index), job.conditions)
        if saved is not None and (saved.get("ligand_id"), saved.get("pose_index")) == (job.ligand_id, job.pose_index):
            return {**saved, "sdf_path": job.sdf,
                    "notes": "; ".join(filter(None, [str(saved.get("notes") or ""),
                        f"worker delivery failed ({type(exc).__name__}: {exc}); "
                        "recovered the condition-verified saved pose record"]))}
    return {"ligand_id": job.ligand_id, "pose_index": job.pose_index,
            "sdf_path": job.sdf, "converged": False,
            "notes": f"FAILED: worker died or did not deliver pose: {type(exc).__name__}: {exc}"}


def mmgbsa_batch(receptor_pdb: Path, poses: Iterable[tuple[str, int, Path]],
                 workdir: Path, n_workers: Optional[int] = None,
                 resume: bool = True, **kwargs: Any) -> "Any":
    """Rescore many poses in parallel across processes.  Returns a DataFrame.

    The receptor is prepared **once** in the parent process and the cleaned
    PDB is handed to the workers, so PDBFixer does not run N times and every
    pose is scored against byte-identical protein coordinates.

    Resumable: each pose's record is written to
    ``workdir/cache/<ligand>__pose<N>.json`` as it finishes.  Re-running with
    ``resume=True`` (default) reads those back and only computes what is
    missing, so an interrupted batch costs nothing to restart.

    A record is reused only when every condition that could move the number is
    the one asked for now: the receptor PDB's content, the pose SDF's content
    (so re-docking a pose invalidates its record) and the whole protocol --
    see :mod:`skinpipe.mmgbsa.cachekey`.  Re-run the same workdir at
    ``implicit_solvent="OBC2"`` and every pose is recomputed, each saying so.

    Raises :class:`~skinpipe.mmgbsa.metal_guard.MetalSiteError` before any pose
    is read or resumed when the receptor's site has lost a catalytic metal the
    force field cannot represent -- every row would be an apo number, so the
    whole table is refused rather than produced and then disbelieved.

    ``n_workers`` defaults to ``min(4, cpu_count)`` (one for a GPU). CPU poses
    default to one thread each. ``cpu_profile="parallel"`` opts into splitting
    cores between the planned workers, or pass ``threads=N`` explicitly.
    The allocation is fixed before cache lookup from the whole requested
    batch, not the number of missing poses. Thread changes invalidate reuse.

    .. warning::
       Workers are spawned, not forked, so Python re-imports the calling
       script in every worker.  **Your entry point must be guarded**::

           if __name__ == "__main__":
               mmgbsa_batch(...)

       Without it the script's top-level code re-runs in every worker, and if
       that code calls ``mmgbsa_batch`` again it escalates into a fork bomb.
       As a backstop this function drops to one worker when it detects it is
       already inside one, but the guard is your job.
    """
    kwargs = dict(kwargs)
    profile = kwargs.get("cpu_profile", "auto")
    count = cpu_threads(kwargs.get("threads"), profile)
    if n_workers is not None:
        if isinstance(n_workers, bool) or not isinstance(n_workers, int):
            raise TypeError("n_workers must be a positive integer or None")
        if n_workers < 1:
            raise ValueError("n_workers must be a positive integer")
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    pose_list = [(str(lid), int(idx), Path(p)) for lid, idx, p in poses]
    if not pose_list:
        return to_frame([])

    workers = min(n_workers or min(4, os.cpu_count() or 1), len(pose_list))
    if n_workers is None:
        try:
            actual, _props = make_platform(kwargs.get("platform"), count, profile)
            if actual.getName() in ("CUDA", "OpenCL", "HIP"):
                workers = 1
        except Exception as exc:              # noqa: BLE001 - sizing only
            LOGGER.debug("could not resolve platform for worker sizing: %s", exc)
    if multiprocessing.parent_process() is not None and workers > 1:
        LOGGER.warning("mmgbsa_batch called from inside a worker; using 1 "
                       "worker (guard your entry point with __main__)")
        workers = 1
    if profile == "parallel" and kwargs.get("threads") is None:
        count = max(1, (os.cpu_count() or 1) // workers)
    kwargs["threads"] = count
    preparation_seed = validate_preparation_seed(
        kwargs.get("preparation_seed", DEFAULT_PREPARATION_SEED))
    # Before the cache, before PDBFixer, before a single pose is scored: one
    # receptor decides every row, so a site that has lost its metal is refused
    # in the first second rather than N times over N hours.
    policy = kwargs.setdefault("on_metal_site", metal_guard.REFUSE)
    metal_guard.check_poses(receptor_pdb, [p for _, _, p in pose_list], policy)
    # Same contract one level along: a structural cofactor that lines this
    # pocket is absent from the receptor `prepare_receptor` is about to strip
    # further, so every row would be scored against a site the poses were not
    # docked into. `pipeline.run_screen` checks this too and refuses with its
    # own message; this is the guard for everyone who does not go through it.
    cofactor_guard.refuse_poses(receptor_pdb, [p for _, _, p in pose_list])
    if kwargs.pop("prepare_receptor_pdb", True):
        receptor_pdb = prepare_receptor(
            Path(receptor_pdb), workdir / f"{Path(receptor_pdb).stem}_prepared.pdb",
            preparation_seed=preparation_seed)
    kwargs["prepare_receptor_pdb"] = False

    # Hashed once: the prepared receptor is the same file for every pose.
    protocol = cachekey.protocol_conditions(_run_single, kwargs, Path(receptor_pdb))

    records: list[dict[str, Any]] = []
    pending: list[_Job] = []
    for ligand_id, pose_index, sdf in pose_list:
        conditions = cachekey.pose_conditions(protocol, ligand_id, pose_index, sdf)
        cached = (cachekey.load_record(_cache_path(workdir, ligand_id, pose_index),
                                       conditions) if resume else None)
        if cached is not None:
            # The conditions bind this supplied SDF's content, not its old
            # location. Return the verified current path after a package move;
            # preserve the original record and numerical result on disk.
            records.append({**cached, "sdf_path": str(sdf)})
            LOGGER.info("resumed %s pose %d from cache", ligand_id, pose_index)
            continue
        pending.append(_Job(receptor=str(receptor_pdb), ligand_id=ligand_id,
                            pose_index=pose_index, sdf=str(sdf),
                            workdir=str(workdir), kwargs=kwargs,
                            conditions=conditions, threads=count))

    if pending:
        workers = min(workers, len(pending))
        LOGGER.info("scoring %d poses on %d workers x %d OpenMM threads",
                    len(pending), workers, count)
        if workers == 1:
            for job in pending:
                try:
                    records.append(_worker(job))
                except Exception as exc:  # Preserve other requested poses on a Python/IPC failure.
                    LOGGER.error("worker for %s pose %d failed: %s", job.ligand_id, job.pose_index, exc)
                    records.append(_worker_failure(job, exc, resume=resume))
        else:
            ctx = multiprocessing.get_context("spawn")
            futures, by_job, accounted = {}, {}, set()
            try:
                with ProcessPoolExecutor(max_workers=workers, mp_context=ctx) as pool:
                    for job in pending:
                        try:
                            future = pool.submit(_worker, job)
                            futures[future] = job
                            by_job[id(job)] = future
                        except Exception as exc:  # A broken pool may reject later submissions.
                            records.append(_worker_failure(job, exc, resume=resume))
                            accounted.add(id(job))
                    for future in as_completed(futures):
                        job = futures[future]
                        try:
                            records.append(future.result())
                        except Exception as exc:  # noqa: BLE001
                            LOGGER.error("worker for %s pose %d died: %s", job.ligand_id, job.pose_index, exc)
                            records.append(_worker_failure(job, exc, resume=resume))
                        accounted.add(id(job))
            except Exception as exc:  # Startup/iteration/shutdown must not discard completed poses.
                LOGGER.error("MM/GBSA worker pool failed: %s", exc)
                for job in pending:
                    if id(job) in accounted:
                        continue
                    future = by_job.get(id(job))
                    try:
                        if future is None or not future.done():
                            raise exc
                        records.append(future.result())
                    except Exception as job_exc:  # noqa: BLE001
                        records.append(_worker_failure(job, job_exc, resume=resume))
                    accounted.add(id(job))

    return to_frame(records)
