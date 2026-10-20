"""Allowlisted docking subprocesses with owned cancellation and measured batches."""
from __future__ import annotations

from collections import Counter
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import re
import shutil
import signal
import subprocess
import threading
import time
import traceback

from .artifacts import ArtifactIntegrityError, atomic_bytes, canonical_json, sha256_file
from .resources import classify_failure, cpu_limit, gpu_admission


class EngineCancelled(RuntimeError):
    pass


class OwnedProcessCleanupError(RuntimeError):
    """Owned native resources could not be certified as stopped."""


class ResourceDeferred(RuntimeError):
    pass


class ComputeBudgetExceeded(RuntimeError):
    pass


def validate_engine_hashes(value):
    """A frozen execution cohort pins both primary and CPU fallback binaries."""
    if (not isinstance(value, dict) or set(value) != {"unidock", "vina"}
            or any(not isinstance(v, str) or not re.fullmatch(r"[0-9a-f]{64}", v)
                   for v in value.values())):
        raise ValueError("Engine cohort requires exact Uni-Dock and Vina binary hashes")


def validate_request(request):
    if "engine_hashes" in request:
        validate_engine_hashes(request["engine_hashes"])
    if request.get("engine", "unidock") not in {"unidock", "vina"}:
        raise ValueError("Only configured Uni-Dock and Vina engines execute docking tasks")
    for name in ("center", "size"):
        values = request.get(name)
        if not isinstance(values, list) or len(values) != 3 or any(
                isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) for value in values):
            raise ValueError(name + " must have three finite coordinates")
        if name == "size" and any(value <= 0 or value > 100 for value in values):
            raise ValueError("Box size must be greater than0 and at most100 Angstrom")
    if not re.fullmatch(r"[0-9a-f]{64}", str(request.get("receptor_sha256", ""))):
        raise ValueError("A receptor content hash is required")
    ligands = request.get("ligands")
    if not isinstance(ligands, list) or not 1 <= len(ligands) <= 1000:
        raise ValueError("Docking task needs1..1000 prepared ligands")
    ids = set()
    for ligand in ligands:
        identifier = ligand.get("candidate_id") if isinstance(ligand, dict) else None
        if not isinstance(identifier, str) or not identifier or len(identifier) > 200 or identifier in ids:
            raise ValueError("Ligand identifiers must be nonempty and unique")
        ids.add(identifier)
        if not re.fullmatch(r"[0-9a-f]{64}", str(ligand.get("pdbqt_sha256", ""))):
            raise ValueError("A prepared ligand content hash is required")
    for key, default, lower, upper in (("seed", 42, 0, 2147483647), ("n_poses", 9, 1, 100), ("exhaustiveness", 16, 1, 128)):
        value = request.get(key, default)
        if isinstance(value, bool) or not isinstance(value, int) or not lower <= value <= upper:
            raise ValueError("Invalid " + key)
    if request.get("search_mode", "balance") not in {"fast", "balance", "detail"}:
        raise ValueError("Unknown Uni-Dock search mode")
    if "canary" in request and not isinstance(request["canary"], bool):
        raise ValueError("canary must be boolean")
    if request.get("canary") and len(ligands) > 32:
        raise ValueError("A canary task has at most32 ligands")


def _gpu_memory(gpu_uuid, pid=None):
    binary = shutil.which("nvidia-smi")
    if not binary:
        return None
    try:
        if pid is None:
            query = [binary, "--id=" + gpu_uuid, "--query-gpu=memory.free", "--format=csv,noheader,nounits"]
        else:
            query = [binary, "--query-compute-apps=pid,gpu_uuid,used_gpu_memory", "--format=csv,noheader,nounits"]
        result = subprocess.run(query, capture_output=True, text=True, timeout=3)
        if result.returncode:
            return None
        if pid is None:
            return float(result.stdout.strip().splitlines()[0])
        for line in result.stdout.splitlines():
            fields = [part.strip() for part in line.split(",")]
            if len(fields) == 3 and fields[0] == str(pid) and fields[1] == gpu_uuid:
                return float(fields[2])
    except (OSError, ValueError, IndexError, subprocess.TimeoutExpired):
        return None
    return None


def _group_extinguished(group_id):
    """An independent POSIX check distinguishes extinct/zombie-only groups.

    Some macOS cleanup races report EPERM for signal-zero after the leader and
    child have already disappeared. Unknown or live membership never counts as
    successful cleanup; only the exact group created by this worker is inspected.
    """
    try:
        listing = subprocess.run(["ps", "-axo", "pgid=,stat="], capture_output=True, text=True, timeout=1)
        if listing.returncode:
            return False
        states = []
        for line in listing.stdout.splitlines():
            fields = line.split()
            if len(fields) != 2:
                return False
            if int(fields[0]) == group_id:
                states.append(fields[1])
        return not states or all(state.startswith("Z") for state in states)
    except (OSError, ValueError, subprocess.TimeoutExpired):
        return False


def terminate_owned(process):
    """Terminate only the session this worker created; never match names/PIDs globally."""
    def wait_owned():
        try:
            process.wait(timeout=3)
        except subprocess.TimeoutExpired as error:
            raise OwnedProcessCleanupError("Owned process did not exit after termination") from error

    if process.poll() is None:
        try:
            group = os.getpgid(process.pid)
        except ProcessLookupError:
            # The leader can exit after poll while its children remain. The
            # group we created still owns their lifetime; continue cleanup.
            group = process.pid
        if group != process.pid:
            raise OwnedProcessCleanupError("Refusing cancellation of an unowned process group")
    def owned_signal(sig):
        try:
            os.killpg(process.pid, sig)
            return True
        except ProcessLookupError:
            return False
        except PermissionError as error:
            if _group_extinguished(process.pid):
                return False
            raise OwnedProcessCleanupError("Owned process cleanup could not be verified") from error

    if not owned_signal(signal.SIGTERM):
        wait_owned()
        return
    deadline = time.monotonic() + 3
    while time.monotonic() < deadline:
        process.poll()
        if _group_extinguished(process.pid) or not owned_signal(0):
            break
        time.sleep(0.05)
    else:
        owned_signal(signal.SIGKILL)
    wait_owned()
    # wait() reaps only the leader. A TERM-ignoring child can survive its
    # parent; certify the entire exact group before acknowledging cleanup.
    deadline = time.monotonic() + 3
    while not _group_extinguished(process.pid):
        if time.monotonic() >= deadline:
            raise OwnedProcessCleanupError("Owned process group remains live after termination")
        time.sleep(0.05)


def cleanup_receipt(*, scope="owned_job_processes"):
    """Issued only after all synchronous owned-stage finally blocks returned.

    This certifies process lifetime, not global host free RAM (file caches and
    unrelated jobs deliberately remain outside this job's cancellation scope).
    """
    return {"owned_processes_stopped": True, "scope": scope,
            "memory_release": "owned_process_exit", "shared_service_preserved": True}


def release_stopped_exception_frames(error):
    """Drop completed chemistry stack locals before publishing its stop ACK."""
    while error is not None:
        following = error.__context__ or error.__cause__
        traceback.clear_frames(error.__traceback__)
        error.__traceback__ = error.__context__ = error.__cause__ = None
        error = following


def _process_rss(pid):
    try:
        path = Path('/proc') / str(pid) / 'statm'
        if path.is_file():
            return int(path.read_text().split()[1]) * os.sysconf('SC_PAGE_SIZE')
        result = subprocess.run(['ps','-o','rss=','-p',str(pid)],capture_output=True,text=True,timeout=1)
        return int(result.stdout.strip()) * 1024 if result.returncode == 0 and result.stdout.strip() else None
    except (OSError,ValueError,subprocess.TimeoutExpired):
        return None


def _group_measurement(group_id):
    """Sample the owned POSIX process group, including Java/provider children."""
    try:
        listing=subprocess.run(['ps','-axo','pid=,pgid=,rss='],capture_output=True,text=True,timeout=1)
        members=[]
        for line in listing.stdout.splitlines():
            pid,pgid,rss=(int(value) for value in line.split())
            if pgid==group_id:
                members.append((pid,rss*1024))
        return [pid for pid,_ in members],sum(rss for _,rss in members) or None
    except (OSError,ValueError,subprocess.TimeoutExpired):
        return [group_id],_process_rss(group_id)


def _gpu_group_memory(gpu_uuid,pids):
    try:
        query=subprocess.run(['nvidia-smi','--query-compute-apps=pid,gpu_uuid,used_memory',
                              '--format=csv,noheader,nounits'],capture_output=True,text=True,timeout=2)
        values=[]
        for line in query.stdout.splitlines():
            fields=[value.strip() for value in line.split(',')]
            if len(fields)==3 and int(fields[0]) in pids and fields[1]==gpu_uuid:
                values.append(float(fields[2]))
        return sum(values) if values else None
    except (OSError,ValueError,subprocess.TimeoutExpired):
        return None


def _wait_owned_gpu_release(gpu_uuid, pids, *, timeout=5):
    """Observe only the stopped job's PIDs, preserving every other GPU user."""
    deadline = time.monotonic() + timeout
    while True:
        try:
            query = subprocess.run(["nvidia-smi", "--query-compute-apps=pid,gpu_uuid,used_gpu_memory",
                "--format=csv,noheader,nounits"], capture_output=True, text=True, timeout=2)
            if query.returncode:
                raise ValueError("GPU process query unavailable")
            live = [fields for line in query.stdout.splitlines()
                if len(fields := [field.strip() for field in line.split(",")]) == 3
                and fields[0] in {str(pid) for pid in pids} and fields[1] == gpu_uuid]
            if not live:
                return
        except (OSError, ValueError, subprocess.TimeoutExpired) as error:
            raise OwnedProcessCleanupError("Owned GPU process release could not be verified") from error
        if time.monotonic() >= deadline:
            raise OwnedProcessCleanupError("Owned GPU context remains after process termination")
        time.sleep(0.1)


def run_owned(command, *, cwd, timeout, cancel_event, gpu_uuid=None, on_compute=None, env=None):
    if cancel_event.is_set():
        raise EngineCancelled("Cancellation requested before process launch")
    if timeout <= 0:
        raise ComputeBudgetExceeded("No actual compute budget remains")
    # Provider processes cannot inherit Python preload paths, dynamic loader
    # injection, cloud credentials or private API tokens from the API/worker.
    environment = {key: os.environ[key] for key in ("PATH", "HOME", "TMPDIR", "LANG", "LC_ALL", "SYSTEMROOT", "WINDIR") if key in os.environ}
    if env is not None:
        if not isinstance(env, dict) or any(not isinstance(k, str) or not isinstance(v, str)
                                            for k, v in env.items()):
            raise ValueError("Owned process environment overrides must be strings")
        environment.update(env)
    if gpu_uuid:
        environment["CUDA_VISIBLE_DEVICES"] = gpu_uuid
    started = time.monotonic()
    peak = None
    peak_rss = None
    log_path = Path(cwd) / "engine.log"
    with log_path.open("wb") as log:
        process = subprocess.Popen(command, cwd=cwd, stdout=log, stderr=subprocess.STDOUT,
                                   env=environment, start_new_session=True)
        try:
            while process.poll() is None:
                if cancel_event.is_set():
                    raise EngineCancelled("Lease lost or cancellation requested")
                elapsed = time.monotonic() - started
                if on_compute:
                    on_compute(elapsed)
                if elapsed >= timeout:
                    raise ComputeBudgetExceeded("Actual compute budget exhausted")
                pids,rss = _group_measurement(process.pid)
                if gpu_uuid:
                    current = _gpu_group_memory(gpu_uuid,pids)
                    if current is not None:
                        peak = max(peak or 0, current)
                if rss is not None:
                    peak_rss = max(peak_rss or 0,rss)
                cancel_event.wait(0.2)
        finally:
            owned_pids, _ = _group_measurement(process.pid)
            terminate_owned(process)
            if gpu_uuid:
                _wait_owned_gpu_release(gpu_uuid, {process.pid, *owned_pids})
            if on_compute:
                on_compute(time.monotonic() - started)
    with log_path.open("rb") as stream:
        stream.seek(max(0, log_path.stat().st_size - 1024 * 1024))
        tail = stream.read().decode(errors="replace")
    return {"returncode": process.returncode, "elapsed_seconds": time.monotonic() - started,
            "peak_mib": peak, "process_peak_rss_bytes": peak_rss,
            "process_tree_peak_rss_bytes": peak_rss,
            "memory_scope": "sampled_process_group_RSS_and_GPU_VRAM", "log": tail, "log_path": str(log_path)}


def _atom_signature(text):
    """Coordinate and atom-type audit without inferring bonding from coordinates."""
    models, current = [], []
    for line in text.splitlines():
        if line.startswith("MODEL"):
            current = []
        elif line.startswith(("ATOM  ", "HETATM")):
            values = [float(line[start:start + 8]) for start in (30, 38, 46)]
            if not all(math.isfinite(value) for value in values):
                raise ValueError("Nonfinite pose coordinates")
            current.append(line.split()[-1])
        elif line.startswith("ENDMDL"):
            models.append(Counter(current))
            current = []
    if current:
        models.append(Counter(current))
    if not models or any(not model for model in models):
        raise ValueError("PDBQT contains no atoms")
    return models


def inspect_output(text, input_text):
    from ..dock.poses import scored_pose_remarks
    from ..dock.unidock import is_plausible
    scores = scored_pose_remarks(text)
    expected = _atom_signature(input_text)[0]
    atoms = _atom_signature(text)
    if len(scores) != len(atoms) or any(signature != expected for signature in atoms):
        raise ValueError("Pose count/atom-type inventory differs from prepared ligand")
    values = [value[0] for value in scores]
    if not values or any(not math.isfinite(value) for value in values):
        raise ValueError("Missing or nonfinite pose score")
    # Keep every raw value; plausibility is an operational review annotation.
    best = min(values)
    return {"status": "complete" if is_plausible(best) else "review_required", "affinity": best,
            "n_poses": len(scores), "pose_scores": values, "within_operational_window": is_plausible(best),
            "pose_validation": "score_model_atom_inventory_and_finite_coordinates",
            "biological_validation": "NOT_ESTABLISHED"}


class EngineRunner:
    def __init__(self, workdir, engine_paths, *, gpu_uuid=None, reserve_mib=1536, memory_reader=None, runner=None,
                 runtime_context=None, required_engine_hashes=None, cpu_threads=None):
        if cpu_threads is not None and (type(cpu_threads) is not int or not 1 <= cpu_threads <= cpu_limit()):
            raise ValueError('Engine CPU threads exceed the host execution limit')
        self.cpu_threads = cpu_threads
        if required_engine_hashes is not None:
            validate_engine_hashes(required_engine_hashes)
        self.required_engine_hashes = dict(required_engine_hashes) if required_engine_hashes is not None else None
        self.workdir = Path(workdir)
        self.workdir.mkdir(parents=True, exist_ok=True)
        self.engine_paths = engine_paths
        self.gpu_uuid = gpu_uuid
        self.reserve_mib = reserve_mib
        self.memory_reader = memory_reader or _gpu_memory
        self.runner = runner or run_owned
        self.runtime_context = dict(runtime_context or {})

    def execute(self, request, *, fetch_artifact, publish_artifact, cancel_event=None, remaining_seconds=28800,
                progress=None, prior_result=None):
        validate_request(request)
        if (self.required_engine_hashes is not None
                and request.get("engine_hashes") != self.required_engine_hashes):
            raise ResourceDeferred("Worker requires its exact frozen engine cohort")
        cancel_event = cancel_event or threading.Event()
        requested_engine = request.get("engine", "unidock")
        configured = self.engine_paths.get(requested_engine)
        binary = shutil.which(configured) if configured else None
        if not binary:
            raise ResourceDeferred("Requested engine not installed on this worker")
        if requested_engine == "unidock" and not self.gpu_uuid:
            raise ResourceDeferred("No GPU UUID assigned to worker")
        binary_identity = {"path": str(Path(binary).resolve()), "sha256": sha256_file(binary)}
        pinned = request.get("engine_hashes")
        if pinned and binary_identity["sha256"] != pinned[requested_engine]:
            raise ResourceDeferred("Worker engine binary differs from frozen cohort")
        gpu = next((value for value in self.runtime_context.get("gpus", []) if value.get("uuid") == self.gpu_uuid), {})
        runtime = {"python": platform.python_version(), "platform": platform.platform(),
                   "hostname": self.runtime_context.get("hostname", platform.node()),
                   "worker_id": self.runtime_context.get("worker_id"),
                   "telemetry_observed_at": self.runtime_context.get("observed_at"),
                   "gpu_uuid": self.gpu_uuid, "gpu_name": gpu.get("name"), "gpu_driver": gpu.get("driver"),
                   "cuda_runtime_version": None, "cuda_runtime_version_status": "not_observed",
                   "engine_binary_sha256": binary_identity["sha256"],
                   "engine_version": self.runtime_context.get("engines", {}).get(requested_engine, {}).get("version"),
                   "execution_source_sha256": sha256_file(Path(__file__))}
        execution_policy = {}
        if requested_engine == "unidock":
            # Uni-Dock 1.x appends parsed ligands from an OpenMP critical
            # section. Parallel completion order changes ligand-index RNG
            # streams despite identical seed/batch bytes. Serial input parsing
            # leaves the GPU search and all scientific settings unchanged.
            execution_policy = {"gpu_input_parse_threads": 1,
                                "policy": "ordered-unidock-input-v1"}
            runtime["execution_policy"] = execution_policy
        protocol = {key: request.get(key, default) for key, default in (("seed", 42), ("n_poses", 9),
                    ("exhaustiveness", 16), ("search_mode", "balance"))}
        protocol.update(engine=requested_engine, center=request["center"], size=request["size"],
                        receptor_sha256=request["receptor_sha256"], engine_identity=binary_identity)
        fallback_binary = shutil.which(self.engine_paths.get("vina", "")) if self.engine_paths.get("vina") else None
        if requested_engine == "unidock":
            protocol["cpu_fallback_identity"] = ({"path": str(Path(fallback_binary).resolve()), "sha256": sha256_file(fallback_binary)}
                                                  if fallback_binary else None)
        if pinned and requested_engine == "unidock" and (
                (protocol.get("cpu_fallback_identity") or {}).get("sha256") != pinned["vina"]):
            raise ResourceDeferred("Worker CPU fallback differs from frozen cohort")
        execution_identity = {"request": request, "protocol": protocol}
        if execution_policy:
            execution_identity["execution_policy"] = execution_policy
        identity = hashlib.sha256(canonical_json(execution_identity)).hexdigest()
        request_content_sha256 = hashlib.sha256(canonical_json(request)).hexdigest()
        directory = self.workdir / identity
        directory.mkdir(exist_ok=True)
        receptor = directory / "receptor.pdbqt"
        self._fetch(request["receptor_sha256"], receptor, fetch_artifact)
        _atom_signature(receptor.read_text())
        inputs = []
        for ligand in request["ligands"]:
            path = directory / (hashlib.sha256(ligand["candidate_id"].encode()).hexdigest()[:24] + ".pdbqt")
            self._fetch(ligand["pdbqt_sha256"], path, fetch_artifact)
            inputs.append({**ligand, "path": path})
        receipt_path = directory / "checkpoint.json"
        complete = {}
        if receipt_path.exists():
            try:
                receipt = json.loads(receipt_path.read_text())
                checksum = receipt.pop("checkpoint_sha256")
                if checksum != hashlib.sha256(canonical_json(receipt)).hexdigest() or receipt.get("identity") != identity:
                    raise ArtifactIntegrityError("Docking checkpoint input/content mismatch")
                allowed = {ligand["candidate_id"] for ligand in inputs}
                for row in receipt.get("rows", []):
                    if row.get("candidate_id") not in allowed:
                        continue
                    artifacts = [row[key] for key in ("artifact", "sdf_artifact") if row.get(key)]
                    valid = True
                    for artifact in artifacts:
                        path = directory / artifact["local_name"]
                        if path.parent != directory or not path.is_file() or sha256_file(path) != artifact["sha256"]:
                            valid = False
                            break
                    if not valid:
                        continue  # Recompute only the corrupt/missing candidate.
                    complete[row["candidate_id"]] = row
            except (OSError, ValueError, KeyError, TypeError):
                # Preserve unusable receipt bytes for inspection before retry.
                quarantined = directory / ("checkpoint-invalid-" + str(time.time_ns()) + ".json")
                atomic_bytes(quarantined, receipt_path.read_bytes())
        checkpoint_errors = []
        if prior_result and prior_result.get("request_content_sha256") != request_content_sha256:
            checkpoint_errors.append({"candidate_id": None, "reason": "central checkpoint request identity mismatch"})
            atomic_bytes(directory / ("central-checkpoint-invalid-" + str(time.time_ns()) + ".json"), canonical_json(prior_result))
            prior_result = None
        if prior_result:
            by_id = {ligand["candidate_id"]: ligand for ligand in inputs}
            seen = set()
            for previous in prior_result.get("candidates", []):
                try:
                    identifier = previous.get("candidate_id")
                    if identifier not in by_id or identifier in seen:
                        raise ArtifactIntegrityError("Central checkpoint candidate identity mismatch")
                    seen.add(identifier)
                    ligand = by_id[identifier]
                    prior_protocol = previous.get("protocol", {})
                    if previous.get("input_sha256") != ligand["pdbqt_sha256"] or any(
                            prior_protocol.get(key) != protocol.get(key) for key in ("receptor_sha256", "center", "size", "seed", "n_poses", "exhaustiveness", "search_mode")):
                        raise ArtifactIntegrityError("Central checkpoint prepared input/protocol mismatch")
                    actual_engine = previous.get("engine")
                    if actual_engine != requested_engine and not (requested_engine == "unidock" and actual_engine == "vina" and
                                                                 previous.get("route_reason") == "macrocycle_glue_requires_cpu_vina"):
                        raise ArtifactIntegrityError("Central checkpoint engine mismatch")
                    expected_identity = (binary_identity if actual_engine == requested_engine
                                         else protocol.get("cpu_fallback_identity")) or {}
                    if (not expected_identity.get("sha256") or
                            prior_protocol.get("engine_identity", {}).get("sha256") != expected_identity["sha256"]):
                        raise ArtifactIntegrityError("Central checkpoint executing binary differs")
                    if (actual_engine == "unidock" and
                            previous.get("runtime", {}).get("execution_policy") != execution_policy):
                        raise ArtifactIntegrityError("Central checkpoint input parsing policy differs")
                    row = dict(previous)
                    if previous.get("pose_sha256"):
                        raw = fetch_artifact(previous["pose_sha256"])
                        if hashlib.sha256(raw).hexdigest() != previous["pose_sha256"]:
                            raise ArtifactIntegrityError("Central checkpoint pose checksum mismatch")
                        if row.get("status") in {"complete", "review_required"}:
                            checked = inspect_output(raw.decode(), ligand["path"].read_text())
                            if any(row.get(key) != checked[key] for key in ("status", "affinity", "pose_scores", "n_poses")):
                                raise ArtifactIntegrityError("Central checkpoint pose metadata differs from pose")
                        copy = directory / ("pose-" + ligand["path"].name)
                        atomic_bytes(copy, raw)
                        row["artifact"] = {"sha256": previous["pose_sha256"], "local_name": copy.name, "size_bytes": len(raw)}
                        if previous.get("sdf_artifact"):
                            sdf_hash = previous["sdf_artifact"]["sha256"]
                            sdf = fetch_artifact(sdf_hash)
                            if hashlib.sha256(sdf).hexdigest() != sdf_hash:
                                raise ArtifactIntegrityError("Central checkpoint SDF checksum mismatch")
                            sdf_copy = directory / ("pose-" + ligand["path"].stem + ".sdf")
                            atomic_bytes(sdf_copy, sdf)
                            row["sdf_artifact"] = {"sha256": sdf_hash, "local_name": sdf_copy.name, "size_bytes": len(sdf)}
                    elif row.get("status") not in {"failed", "unsupported_chemistry"}:
                        raise ArtifactIntegrityError("Central checkpoint lacks its completed pose")
                    complete[identifier] = row
                except (OSError, ValueError, KeyError, TypeError, AttributeError) as error:
                    rejected = {"candidate_id": previous.get("candidate_id") if isinstance(previous, dict) else None, "reason": str(error)}
                    checkpoint_errors.append(rejected)
                    atomic_bytes(directory / ("central-row-invalid-" + str(time.time_ns()) + ".json"), canonical_json({"row": previous, "error": rejected}))
        pending = [ligand for ligand in inputs if ligand["candidate_id"] not in complete]
        elapsed = 0.0
        measured_peak, measured_size, size = None, 1, 1
        measurements = []

        def checkpoint():
            receipt = {"identity": identity, "rows": list(complete.values()), "protocol": protocol}
            receipt["checkpoint_sha256"] = hashlib.sha256(canonical_json(receipt)).hexdigest()
            atomic_bytes(receipt_path, canonical_json(receipt))
            if progress:
                progress({"stage": "docking", "completed": len(complete), "total": len(inputs), "compute_seconds": elapsed,
                          "partial_result": {"request_content_sha256": request_content_sha256,
                                             "candidates": list(complete.values()), "protocol": protocol,
                                             "runtime": runtime}})

        while pending:
            if cancel_event.is_set():
                raise EngineCancelled("Lease lost or cancellation requested")
            if elapsed >= remaining_seconds:
                raise ComputeBudgetExceeded("Actual compute budget exhausted; completed poses checkpointed")
            batch = pending[:size if requested_engine == "unidock" else 1]
            from ..dock.unidock import has_glue_atoms
            unsupported = [ligand for ligand in batch if requested_engine == "unidock" and has_glue_atoms(ligand["path"].read_text())]
            for ligand in unsupported:
                if fallback_binary:
                    fallback_request = {**request, "engine": "vina", "ligands": [{key: ligand[key] for key in ("candidate_id", "pdbqt_sha256")}], "canary": False}
                    before_cpu = elapsed

                    def cpu_progress(value):
                        if progress:
                            update = {key: val for key, val in value.items() if key != "partial_result"}
                            update.update(compute_seconds=before_cpu + value.get("compute_seconds", 0), stage="cpu_fallback")
                            if "partial_result" in value:
                                merged = dict(complete)
                                for child_row in value["partial_result"].get("candidates", []):
                                    merged[child_row["candidate_id"]] = {**child_row, "route_reason": "macrocycle_glue_requires_cpu_vina"}
                                update["partial_result"] = {"request_content_sha256": request_content_sha256,
                                                            "candidates": list(merged.values()), "protocol": protocol,
                                                            "runtime": runtime}
                            progress(update)

                    fallback = EngineRunner(self.workdir / "cpu-fallback", self.engine_paths, runner=self.runner,
                                            runtime_context=self.runtime_context, cpu_threads=self.cpu_threads).execute(
                        fallback_request, fetch_artifact=fetch_artifact, publish_artifact=publish_artifact,
                        cancel_event=cancel_event, remaining_seconds=remaining_seconds - elapsed, progress=cpu_progress)
                    elapsed += fallback["compute_seconds"]
                    row = dict(fallback["candidates"][0])
                    row["route_reason"] = "macrocycle_glue_requires_cpu_vina"
                    row.setdefault("protocol", fallback["protocol"])
                    # Parent checkpoint references its own immutable local bytes.
                    if row.get("artifact"):
                        artifact_path = self.workdir / "cpu-fallback" / fallback["request_sha256"] / row["artifact"]["local_name"]
                        copy = directory / ("pose-" + ligand["path"].name)
                        atomic_bytes(copy, artifact_path.read_bytes())
                        row["artifact"] = {**row["artifact"], "local_name": copy.name}
                    if row.get("sdf_artifact"):
                        artifact_path = self.workdir / "cpu-fallback" / fallback["request_sha256"] / row["sdf_artifact"]["local_name"]
                        copy = directory / ("pose-" + ligand["path"].stem + ".sdf")
                        atomic_bytes(copy, artifact_path.read_bytes())
                        row["sdf_artifact"] = {**row["sdf_artifact"], "local_name": copy.name}
                    complete[ligand["candidate_id"]] = row
                else:
                    complete[ligand["candidate_id"]] = {"candidate_id": ligand["candidate_id"], "status": "unsupported_chemistry",
                                                        "affinity": None, "score": None, "reason": "macrocycle_glue_requires_cpu_vina",
                                                        "input_sha256": ligand["pdbqt_sha256"], "protocol": protocol,
                                                        "runtime": runtime, "engine": requested_engine}
                pending.remove(ligand)
                batch.remove(ligand)
            if unsupported:
                checkpoint()
            if not batch:
                continue
            if requested_engine == "unidock":
                free = self.memory_reader(self.gpu_uuid)
                if free is None:
                    raise ResourceDeferred("GPU free-memory telemetry unavailable")



                projected = 512 if measured_peak is None else measured_peak * size / measured_size
                admission = gpu_admission(free, projected, reserve_mib=self.reserve_mib)
                while not admission["admitted"] and size > 1:
                    size = max(1, size // 2)
                    projected = measured_peak * size / measured_size
                    admission = gpu_admission(free, projected, reserve_mib=self.reserve_mib)
                if not admission["admitted"]:
                    raise ResourceDeferred("Insufficient free VRAM for measured batch and reserved headroom")
                batch = batch[:size]
            attempt = directory / ("attempt-" + str(len(measurements)) + "-" + str(time.time_ns()))
            attempt.mkdir()
            command = [binary, "--receptor", str(receptor)]
            for axis, center, width in zip("xyz", request["center"], request["size"]):
                command += ["--center_" + axis, str(center), "--size_" + axis, str(width)]
            command += ["--seed", str(protocol["seed"]), "--num_modes", str(protocol["n_poses"])]
            if requested_engine == "unidock":
                command += ["--gpu_batch", *[str(ligand["path"]) for ligand in batch], "--search_mode", protocol["search_mode"], "--dir", str(attempt)]
            else:
                command += ["--ligand", str(batch[0]["path"]), "--out", str(attempt / (batch[0]["path"].stem + "_out.pdbqt")),
                            "--cpu", str(self.cpu_threads or cpu_limit()), "--exhaustiveness", str(protocol["exhaustiveness"])]
            before = elapsed

            def compute_update(seconds):
                if progress:
                    progress({"stage": "docking", "completed": len(complete), "total": len(inputs),
                              "batch_size": len(batch), "compute_seconds": before + seconds})

            outcome = self.runner(command, cwd=attempt, timeout=remaining_seconds - elapsed,
                                  cancel_event=cancel_event, gpu_uuid=self.gpu_uuid if requested_engine == "unidock" else None,
                                  on_compute=compute_update, env={'OMP_NUM_THREADS':'1' if requested_engine == 'unidock' else str(cpu_limit()),
                                      'MKL_NUM_THREADS':str(cpu_limit()),'OPENBLAS_NUM_THREADS':str(cpu_limit())})
            elapsed += outcome["elapsed_seconds"]
            measurements.append({key: outcome[key] for key in ("elapsed_seconds", "peak_mib", "returncode")})
            failure = classify_failure(outcome["log"]) if outcome["returncode"] else None
            if failure == "oom":
                if len(batch) == 1:
                    checkpoint()
                    raise ResourceDeferred("OOM on single-ligand canary; awaiting another worker or more VRAM")
                size = max(1, len(batch) // 2)
                continue
            if outcome["returncode"] and len(batch) > 1:
                size = max(1, len(batch) // 2)  # Isolate chemistry/engine failures to the affected ligand.
                continue
            if outcome["peak_mib"] is not None:
                measured_peak, measured_size = outcome["peak_mib"], len(batch)
            members = [{"candidate_id": ligand["candidate_id"], "pdbqt_sha256": ligand["pdbqt_sha256"]} for ligand in batch]
            actual_protocol = {**protocol, "batch": {"members": members, "size": len(batch),
                                                      "membership_sha256": hashlib.sha256(canonical_json(members)).hexdigest()}}
            for ligand in batch:
                identifier = ligand["candidate_id"]
                output = attempt / (ligand["path"].stem + "_out.pdbqt")
                row = {"candidate_id": identifier, "engine": requested_engine, "input_sha256": ligand["pdbqt_sha256"],
                       "affinity": None, "protocol": actual_protocol,
                       "runtime": runtime,
                       "execution_attempt": {"local_id": attempt.name, "batch_elapsed_seconds": outcome["elapsed_seconds"],
                                             "observed_peak_mib": outcome["peak_mib"]}}
                if outcome["returncode"] or not output.exists():
                    row.update(status="failed", reason=failure or "missing_pose", error=outcome["log"][-2000:])
                else:
                    raw = output.read_bytes()
                    try:
                        row.update(inspect_output(raw.decode(), ligand["path"].read_text()))
                    except (ValueError, UnicodeDecodeError) as error:
                        row.update(status="invalid_pose", reason=str(error))
                    artifact = publish_artifact(raw)
                    copy = directory / ("pose-" + ligand["path"].name)
                    atomic_bytes(copy, raw)
                    row["artifact"] = {**artifact, "local_name": copy.name}
                    row["pose_sha256"] = artifact["sha256"]
                    # The original Meeko reconstruction retains its atom mapping.
                    # If the optional stack is unavailable, raw PDBQT remains a
                    # traceable result and is not mislabeled as a validated SDF.
                    try:
                        from ..dock.poses import mols_from_docked_pdbqt, write_pose_sdf
                        mol = mols_from_docked_pdbqt(raw.decode())
                        if mol is not None and mol.GetNumConformers():
                            sdf = directory / ("pose-" + ligand["path"].stem + ".sdf")
                            scores = row.get("pose_scores", [])
                            chosen = min(range(len(scores)), key=scores.__getitem__) if scores else 0
                            if chosen < mol.GetNumConformers() and write_pose_sdf(mol, chosen, sdf, {"candidate_id": identifier, "pose_index": chosen + 1}):
                                row["sdf_artifact"] = {**publish_artifact(sdf.read_bytes()), "local_name": sdf.name}
                    except (ImportError, ValueError, RuntimeError, UnicodeDecodeError) as error:
                        row["sdf_status"] = "unavailable: " + type(error).__name__
                complete[identifier] = row
                row["score"] = row["affinity"]
                pending.remove(ligand)
                checkpoint()
            if requested_engine == "unidock":
                if measured_peak is None and pending:
                    if all(complete[ligand["candidate_id"]]["status"] in
                           {"failed", "unsupported_chemistry", "invalid_pose"} for ligand in batch):
                        # A rejected ligand can finish before any GPU process
                        # is measurable. Preserve its failure and use another
                        # single-ligand canary; it is not a resource shortage.
                        size = 1
                        continue
                    raise ResourceDeferred("Canary returned without process peak telemetry; batch expansion withheld")
                size = min(32, max(1, len(batch) * 2))
        return {"protocol": protocol, "request_sha256": identity, "request_content_sha256": request_content_sha256,
                "candidates": [complete[ligand["candidate_id"]] for ligand in inputs],
                "compute_seconds": elapsed, "measurements": measurements, "checkpoint_errors": checkpoint_errors,
                "runtime": runtime,
                "advanced_models": {"gnina": "not_configured", "boltz": "not_configured"}}

    @staticmethod
    def _fetch(digest, destination, fetch):
        if destination.exists() and sha256_file(destination) == digest:
            return
        data = fetch(digest)
        if hashlib.sha256(data).hexdigest() != digest:
            raise ArtifactIntegrityError("Downloaded docking input checksum mismatch")
        atomic_bytes(destination, data)
