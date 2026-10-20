"""Credentialed Tailscale HTTP worker. No shell commands are accepted from jobs.

Run: python -m skinpipe.platform.worker --config /private/worker.json
The config must have mode0600 and api_url, token, workdir and engine paths.
"""
from __future__ import annotations

import argparse
import base64
import hashlib
import ipaddress
import json
import os
from pathlib import Path
import signal
import stat
import threading
import time
from urllib import error, parse, request

from .artifacts import ArtifactIntegrityError, canonical_json
from .engines import ComputeBudgetExceeded, EngineCancelled, EngineRunner, OwnedProcessCleanupError, ResourceDeferred, validate_engine_hashes
from .resources import probe, gpu_admission
from .algorithm_profiles import PROVIDERS
from .algorithm_execution import AlgorithmRunner


class WorkerHTTPError(RuntimeError):
    def __init__(self, status, message):
        super().__init__(message)
        self.status = status


class _NoRedirect(request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        # Never forward a worker credential to a redirect destination.
        return None


class APIClient:
    def __init__(self, api_url, token, *, timeout=10):
        parsed = parse.urlsplit(api_url)
        if parsed.scheme not in {"http", "https"} or not parsed.hostname or parsed.username or parsed.password or parsed.query or parsed.fragment:
            raise ValueError("Invalid worker API URL")
        if parsed.scheme == "http":
            allowed = parsed.hostname in {"localhost", "127.0.0.1", "::1"} or parsed.hostname.endswith(".ts.net")
            try:
                address = ipaddress.ip_address(parsed.hostname)
                allowed = allowed or address in ipaddress.ip_network("100.64.0.0/10") or address.is_loopback
            except ValueError:
                pass
            if not allowed:
                raise ValueError("Plain HTTP worker API must use loopback or a Tailscale address")
        if not isinstance(token, str) or not token:
            raise ValueError("A per-node token is required")
        self.url = api_url.rstrip("/")
        self.token = token
        self.timeout = timeout
        self.opener = request.build_opener(request.ProxyHandler({}), _NoRedirect())

    def call(self, path, payload=None, *, method=None):
        if not path.startswith("/") or path.startswith("//"):
            raise ValueError("Only API-relative paths are allowed")
        req = request.Request(self.url + path, data=canonical_json(payload) if payload is not None else None,
                              headers={"Authorization": "Bearer " + self.token, "Content-Type": "application/json"},
                              method=method or ("POST" if payload is not None else "GET"))
        try:
            with self.opener.open(req, timeout=self.timeout) as response:
                raw = response.read(96 * 1024 * 1024 + 1)
                if len(raw) > 96 * 1024 * 1024:
                    raise ValueError("Worker response exceeds96MiB limit")
                return json.loads(raw)
        except error.HTTPError as exception:
            body = exception.read(4096).decode(errors="replace")
            raise WorkerHTTPError(exception.code, body) from exception

    def fetch_artifact(self, digest):
        response = self.call("/v2/worker/artifacts/" + digest)
        if response.get("sha256") != digest:
            raise ArtifactIntegrityError("Artifact response identity differs from request")
        data = base64.b64decode(response["base64"], validate=True)
        if hashlib.sha256(data).hexdigest() != digest:
            raise ArtifactIntegrityError("Artifact response checksum mismatch")
        return data

    def publish_artifact(self, data):
        digest = hashlib.sha256(data).hexdigest()
        response = self.call("/v2/worker/artifacts", {"base64": base64.b64encode(data).decode(), "sha256": digest})
        artifact = response.get("artifact", response)
        if artifact.get("sha256") != digest:
            raise ArtifactIntegrityError("Published artifact receipt checksum mismatch")
        return artifact


def load_config(path):
    path = Path(path)
    info = path.stat()
    if stat.S_IMODE(info.st_mode) != 0o600 or info.st_uid != os.getuid():
        raise ValueError("Worker credentials config must be owned by this user with mode0600")
    config = json.loads(path.read_text())
    for key in ("api_url", "token", "workdir"):
        if not isinstance(config.get(key), str) or not config[key]:
            raise ValueError("Missing worker config: " + key)
    paths = config.get("engines", {})
    if not isinstance(paths, dict) or set(paths) - {"unidock", "vina"}:
        raise ValueError("Worker engine paths allow only unidock/vina")
    for value in paths.values():
        if value is not None and (not isinstance(value, str) or not Path(value).is_absolute()):
            raise ValueError("Configured engine paths must be absolute")
    if "engine_hashes" in config:
        validate_engine_hashes(config["engine_hashes"])
    profiles = config.get("algorithm_profiles", {})
    if not isinstance(profiles, dict) or set(profiles) - (set(PROVIDERS)|{'evaluation'}) or any(not isinstance(p,dict) for p in profiles.values()):
        raise ValueError("Optional algorithm profiles must be operator-owned provider objects")
    for key, default, lower, upper in (("poll_seconds", 5, 0.1, 60), ("heartbeat_seconds", 2, 0.2, 20),
                                      ("lease_grace_seconds", 25, 2, 60)):
        value = config.get(key, default)
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not lower <= value <= upper:
            raise ValueError("Invalid worker setting: " + key)
    return config


class Worker:
    def __init__(self, config, *, client=None, engine_runner=None):
        if "engine_hashes" in config:
            validate_engine_hashes(config["engine_hashes"])
        self.config = config
        self.client = client or APIClient(config["api_url"], config["token"])
        self.workdir = Path(config["workdir"]).resolve()
        self.workdir.mkdir(parents=True, exist_ok=True)
        self.stop = threading.Event()
        self.active_cancel = None
        self.engine_runner = engine_runner

    def shutdown(self):
        self.stop.set()
        if self.active_cancel:
            self.active_cancel.set()

    def telemetry(self):
        data = probe(directory=self.workdir, engine_paths=self.config.get("engines", {}))
        data["capabilities"] = ["probe"] + (["dock_batch"] if "engine_hashes" in self.config
                and any(item.get("available") for item in data["engines"].values()) else [])
        data["algorithms"] = {name: {"status": "configured_qualification_required", "version": profile.get("version")}
                              for name,profile in self.config.get("algorithm_profiles", {}).items()}
        data["capabilities"].extend("algorithm_" + name for name in data["algorithms"])
        return data

    def run_once(self):
        telemetry = self.telemetry()
        if "dock_batch" in telemetry["capabilities"] and "engine_hashes" not in self.config:
            raise ResourceDeferred("Docking claims require a verified exact engine pair")
        registered = self.client.call("/v2/worker/register", {"telemetry": telemetry})
        claim = {"kinds": telemetry["capabilities"]}
        if "engine_hashes" in self.config:
            if "engine_hashes_v1" not in registered.get("claim_contracts", []):
                raise ResourceDeferred("Server lacks strict engine cohort claim admission")
            claim["engine_hashes"] = self.config["engine_hashes"]
        claimed = self.client.call("/v2/worker/claim", claim)
        job = claimed.get("job")
        if not job:
            return False
        ownership = {"job_id": job["id"], "lease_token": job["lease_token"]}
        cancel = threading.Event()
        finished = threading.Event()
        self.active_cancel = cancel
        base_compute = float(job.get("compute_seconds", 0))
        state = {"compute_seconds": base_compute, "progress": {**(job.get("progress") or {}), "stage": "preparing"},
                 "last_heartbeat": time.monotonic(), "lease_error": None}
        lock = threading.Lock()
        heartbeat_lock = threading.Lock()

        def update(progress):
            if "partial_result" in progress:
                progress = {**progress, "partial_result": {**progress["partial_result"], "worker_telemetry": telemetry}}
            with lock:
                current = base_compute + float(progress.get("compute_seconds", 0))
                state["compute_seconds"] = max(state["compute_seconds"], current)
                state["progress"].update({key: value for key, value in progress.items() if key != "compute_seconds"})
            if "partial_result" in progress:
                # Durably fence the receipt immediately after CAS publication;
                # do not leave completed poses only in a ten-second RAM window.
                heartbeat()

        def heartbeat():
            with heartbeat_lock:
                with lock:
                    payload = {**ownership, "progress": dict(state["progress"]), "compute_seconds": state["compute_seconds"]}
                response = self.client.call("/v2/worker/heartbeat", payload)
                state["last_heartbeat"] = time.monotonic()
                if response.get("status") != "running":
                    state["lease_error"] = response.get("status", "lease_lost")
                    cancel.set()

        def keep_lease():
            while not finished.wait(self.config.get("heartbeat_seconds", 2)):
                try:
                    heartbeat()
                except WorkerHTTPError as exception:
                    if exception.status in {400, 401, 403, 404, 409, 410}:
                        state["lease_error"] = "lease_or_authorization_lost"
                        cancel.set()
                        return
                except (OSError, ValueError):
                    pass
                except Exception as exception:
                    # Unexpected response/client failures must not silently
                    # kill the keeper while its native child keeps computing.
                    state["lease_error"] = "heartbeat_" + type(exception).__name__
                    cancel.set()
                    return
                if time.monotonic() - state["last_heartbeat"] >= self.config.get("lease_grace_seconds", 25):
                    state["lease_error"] = "heartbeat_deadline_exceeded"
                    cancel.set()
                    return

        keeper = threading.Thread(target=keep_lease, daemon=True, name="worker-lease")
        keeper.start()
        result = None
        cleanup_known, cleanup_ack_needed = True, False
        try:
            if (job["kind"] == "dock_batch" and "engine_hashes" in self.config
                    and job["request"].get("engine_hashes") != self.config["engine_hashes"]):
                raise ResourceDeferred("Worker requires its exact frozen engine cohort")
            if job["kind"] != "probe":
                cores = telemetry.get("cpu_limit",1) if job["kind"] == "dock_batch" else 1
                cpu = self.client.call("/v2/worker/reserve-cpu", {**ownership, "cores": cores})
                if not cpu.get("reserved"):
                    raise ResourceDeferred("Aggregate node CPU allocation is in use")
            if job["kind"] == "probe":
                result = {"telemetry": telemetry, "status": "complete", "scientific_validation": "NOT_ESTABLISHED"}
            elif job["kind"] == "dock_batch":
                payload = job["request"]
                gpu = None
                if payload.get("engine", "unidock") == "unidock":
                    selected = self.config.get("gpu_uuid")
                    available = [g for g in telemetry["gpus"] if selected is None or g["uuid"] == selected]
                    if not available:
                        raise ResourceDeferred("No configured GPU is present")
                    gpu = max(available, key=lambda row: row["free_mib"])["uuid"]
                    reservation = self.client.call("/v2/worker/reserve", {**ownership, "gpu_uuid": gpu})
                    if reservation.get("reserved") is False or reservation.get("status") in {"busy", "unavailable"}:
                        raise ResourceDeferred("GPU is reserved by another platform worker")
                engine = self.engine_runner or EngineRunner(self.workdir / "docking", self.config.get("engines", {}), gpu_uuid=gpu,
                                                            runtime_context={**telemetry, "worker_id": job.get("worker")},
                                                            required_engine_hashes=self.config.get("engine_hashes"))
                result = engine.execute(payload, fetch_artifact=self.client.fetch_artifact,
                                        publish_artifact=self.client.publish_artifact, cancel_event=cancel,
                                        remaining_seconds=float(job.get("budget_seconds", 28800)) - base_compute, progress=update,
                                        prior_result=(job.get("progress") or {}).get("partial_result"))
                state["compute_seconds"] = max(state["compute_seconds"], base_compute + result.get("compute_seconds", 0))
                result["worker_telemetry"] = telemetry
            elif job['kind']=='algorithm_evaluation':
                from .evaluation_execution import EvaluationRunner
                runner=EvaluationRunner(self.workdir/'evaluation',self.config.get('algorithm_profiles',{}).get('evaluation'))


                evaluation_request = {key: value for key, value in job['request'].items() if key != 'campaign_id'}
                result=runner.execute(evaluation_request,fetch_artifact=self.client.fetch_artifact,
                    publish_artifact=self.client.publish_artifact,cancel_event=cancel,
                    remaining_seconds=float(job.get('budget_seconds',3600))-base_compute,progress=update)
                state['compute_seconds']=max(state['compute_seconds'],base_compute+result.get('compute_seconds',0))
                result['worker_telemetry']=telemetry
            elif job["kind"].startswith("algorithm_"):
                provider = job["request"]["provider"]
                if job["kind"] != "algorithm_" + provider or provider not in self.config.get("algorithm_profiles", {}):
                    raise ValueError("Optional algorithm worker capability mismatch")
                gpu = None
                if PROVIDERS[provider][1] == "gpu":
                    available = [g for g in telemetry["gpus"] if self.config.get("gpu_uuid") in {None,g["uuid"]}]
                    if not available:
                        raise ResourceDeferred("No configured algorithm GPU is present")
                    selected = max(available, key=lambda g:g['free_mib'])
                    pinned = self.config['algorithm_profiles'][provider]
                    receipt = pinned.get('qualification_receipt')
                    if isinstance(receipt,dict):
                        from .artifacts import sha256_file
                        path = Path(receipt['path'])
                        if sha256_file(path) != receipt['sha256']:
                            raise ValueError('GPU qualification receipt changed')
                        measured = json.loads(path.read_text())
                        if measured.get('memory_scope') != 'gpu_process_peak' or measured.get('peak_memory_mib',0) <= 0:
                            raise ResourceDeferred('GPU process VRAM qualification is required')
                        if not gpu_admission(selected['free_mib'], measured['peak_memory_mib'])['admitted']:
                            raise ResourceDeferred('Measured algorithm VRAM plus reserve unavailable')
                    elif job['request'].get('qualification_mode') is not True:
                        raise ResourceDeferred('A bounded GPU qualification canary is required')
                    elif selected['free_mib'] <= 1536:
                        raise ResourceDeferred('Canary VRAM reserve unavailable')
                    gpu = selected['uuid']
                    reservation = self.client.call("/v2/worker/reserve", {**ownership, "gpu_uuid": gpu})
                    if not reservation.get("reserved"):
                        raise ResourceDeferred("Algorithm GPU is reserved by another platform worker")
                runner = AlgorithmRunner(self.workdir / "algorithms", self.config["algorithm_profiles"], gpu_uuid=gpu)
                def reserve_launch(key):
                    if not self.client.call('/v2/worker/reserve-model-launch',{**ownership,'work_key':key}).get('reserved'):
                        return False
                    return self.client.call('/v2/worker/start-model-launch',{**ownership,'work_key':key}).get('authorized',False)
                result = runner.execute(job["request"], fetch_artifact=self.client.fetch_artifact,
                    publish_artifact=self.client.publish_artifact, cancel_event=cancel,
                    remaining_seconds=float(job.get("budget_seconds",3600))-base_compute, progress=update,
                    prior_result=(job.get("progress") or {}).get("partial_result"),
                    reserve_launch=reserve_launch)
                state["compute_seconds"] = max(state["compute_seconds"],base_compute+result.get("compute_seconds",0))
                result["worker_telemetry"] = telemetry
            else:
                raise ValueError("Unimplemented task kind")
            if cancel.is_set() or self.stop.is_set():
                raise EngineCancelled("Job ownership ended before publication")
            # Serialize final accounting with heartbeat: stop the keeper before
            # the atomic publication to avoid a harmless post-finish heartbeat.
            finished.set()
            keeper.join(timeout=12)
            heartbeat()
            if cancel.is_set():
                raise EngineCancelled("Job lease ended before publication")
            from .engines import cleanup_receipt
            self.client.call("/v2/worker/finish", {**ownership, "result": result,
                "compute_seconds": state["compute_seconds"], "cleanup": cleanup_receipt()})
        except ResourceDeferred as exception:
            cleanup_ack_needed = True
            try:
                heartbeat()
                if not cancel.is_set():
                    self.client.call("/v2/worker/defer", {**ownership, "reason": str(exception)})
            except (WorkerHTTPError, OSError):
                pass  # Server lease expiration is the recovery boundary.
        except ComputeBudgetExceeded as exception:
            cleanup_ack_needed = True
            try:
                heartbeat()
                self.client.call("/v2/worker/fail", {**ownership, "error": str(exception), "status": "budget_deferred"})
            except (WorkerHTTPError, OSError):
                pass
        except OwnedProcessCleanupError:
            cleanup_known = False
            try:
                self.client.call("/v2/worker/fail", {**ownership,
                    "error": "owned_resource_cleanup_unverified", "status": "needs_data"})
            except (WorkerHTTPError, OSError):
                pass
        except (EngineCancelled, WorkerHTTPError) as error:
            cleanup_ack_needed = True
            # A disconnected/shutting-down worker is not a user cancellation.
            # The server accepts this acknowledgement only for an explicit
            # cancellation and the still-owned lease. All synchronous runner
            # cleanup has returned; cleanup failures use a distinct exception.
            from .engines import cleanup_receipt, release_stopped_exception_frames
            finished.set()
            keeper.join(timeout=12)
            release_stopped_exception_frames(error)
            result = None
            with lock:
                state["progress"].clear()
            try:
                self.client.call("/v2/worker/cancel-ack", {**ownership,
                    "compute_seconds": state["compute_seconds"], "cleanup": cleanup_receipt()})
            except (WorkerHTTPError, OSError, ValueError):
                pass  # An unacknowledged/expired lease never proves release.
        except OSError:
            # Lost transfer/network: retain local receipts and let the lease
            # expire rather than publishing a false terminal chemistry failure.
            cleanup_ack_needed = True
        except Exception as exception:
            cleanup_ack_needed = True
            try:
                self.client.call("/v2/worker/fail", {**ownership, "error": type(exception).__name__ + ": " + str(exception)[:1800], "status": "retryable"})
            except (WorkerHTTPError, OSError):
                pass
        finally:
            finished.set()
            keeper.join(timeout=12)
            if cleanup_known and (cleanup_ack_needed or state.get('lease_error') == 'budget_deferred'):
                from .engines import cleanup_receipt
                try:
                    self.client.call('/v2/worker/cleanup-ack', {**ownership, 'cleanup': cleanup_receipt()})
                except (WorkerHTTPError, OSError, ValueError):
                    pass
            self.active_cancel = None
        return True

    def run(self):
        while not self.stop.is_set():
            try:
                worked = self.run_once()
            except (WorkerHTTPError, OSError, ValueError):
                worked = False
            if not worked:
                self.stop.wait(self.config.get("poll_seconds", 5))


def main(argv=None):
    parser = argparse.ArgumentParser(description="SkinPipe internal leased worker")
    parser.add_argument("--config", required=True)
    parser.add_argument("--once", action="store_true")
    arguments = parser.parse_args(argv)
    worker = Worker(load_config(arguments.config))
    for signum in (signal.SIGTERM, signal.SIGINT):
        signal.signal(signum, lambda *_: worker.shutdown())
    if arguments.once:
        worker.run_once()
    else:
        worker.run()


if __name__ == "__main__":
    main()
