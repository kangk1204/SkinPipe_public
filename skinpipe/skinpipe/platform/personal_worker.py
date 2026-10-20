"""Credential-free in-process worker limited to a verified personal install.

The normal lease, budget, artifact and engine-cohort contracts stay intact.
An installed Vina executable alone cannot impersonate a pinned GPU/CPU pair.
"""

import json
import os
from pathlib import Path
import shutil
from urllib.parse import urlsplit

from .artifacts import sha256_file
from .engines import validate_engine_hashes
from .queue import Cancelled, LeaseLost
from .worker import APIClient, Worker, WorkerHTTPError


NODE_ID = "personal-cpu"


def verify_personal_install(service, root, config):
    root = Path(root).resolve()
    info = root.stat()
    if info.st_uid != os.getuid() or info.st_mode & 0o077:
        raise ValueError("Personal installation must be owned by this user and private")
    marker = json.loads((root / "personal-install.json").read_bytes())
    database = config.get("database", {})
    if (marker.get("schema") != "skinpipe-personal-install-v1"
            or marker.get("root") != str(root)
            or marker.get("database") != database
            or marker.get("legacy_recovery") is not False
            or config.get("legacy")
            or database.get("dbname") != "skinpipe_research"
            or not isinstance(database.get("host"), str)
            or not Path(database["host"]).is_absolute()
            or Path(config.get("root", "")).resolve() != root / "research-data"
            or Path(service.root).resolve() != root / "research-data"):
        raise ValueError("In-process worker requires its own fresh personal-install ledger and Unix socket")
    return root


class InProcessWorkerClient(APIClient):
    """Reuse artifact verification while avoiding a credential or HTTP server."""

    def __init__(self, service):
        self.service = service

    def call(self, path, payload=None, *, method=None):
        parsed = urlsplit(path)
        if (parsed.scheme or parsed.netloc or parsed.query or parsed.fragment
                or not parsed.path.startswith("/v2/worker/")
                or method not in {None, "GET", "POST"}):
            raise ValueError("Only installed worker API routes are available in process")
        try:
            if payload is None and method != "POST":
                return self.service.worker_get(path, NODE_ID)
            if not isinstance(payload, dict) or method == "GET":
                raise ValueError("Worker mutations require an object payload")
            return self.service.worker_post(path, payload, NODE_ID)
        except (LeaseLost, Cancelled) as error:
            raise WorkerHTTPError(409, type(error).__name__) from error
        except PermissionError as error:
            raise WorkerHTTPError(403, type(error).__name__) from error
        except KeyError as error:
            raise WorkerHTTPError(404, type(error).__name__) from error
        except ValueError as error:
            raise WorkerHTTPError(400, str(error)) from error


class PersonalWorker:
    def __init__(self, service, root, config):
        self.root = verify_personal_install(service, root, config)
        self.service = service
        discovery = config.get("discovery") or {}
        supplied = discovery.get("personal_worker") or {}
        if not isinstance(supplied, dict):
            raise ValueError("personal_worker must be an operator-owned object")
        engines = supplied.get("engines") or (discovery.get("research_docking") or {}).get("engine_paths", {})
        if not isinstance(engines, dict) or set(engines) - {"unidock", "vina"}:
            raise ValueError("Personal worker allows configured Uni-Dock and Vina binaries only")
        worker_config = {"workdir": str(self.root / "worker-work"), "engines": engines,
                         "heartbeat_seconds": 1, "lease_grace_seconds": 10}
        self.readiness = {"node_id": NODE_ID, "status": "held", "capabilities": ["probe"],
                          "reason": "exact_engine_pair_required_for_existing_campaigns",
                          "scope": "New personal installation only; no old-result recovery"}
        hashes = supplied.get("engine_hashes")
        if hashes is not None:
            validate_engine_hashes(hashes)
            for engine, expected in hashes.items():
                configured = engines.get(engine)
                executable = shutil.which(configured) if configured else None
                if not executable or sha256_file(executable) != expected:
                    raise ValueError("Personal worker binary differs from its frozen exact engine pair")
            worker_config["engine_hashes"] = dict(hashes)
            if "gpu_uuid" in supplied:
                worker_config["gpu_uuid"] = supplied["gpu_uuid"]
            self.readiness.update(status="configured", capabilities=["probe", "dock_batch"],
                                  reason="Actual job, GPU and frozen campaign cohort admission remains required",
                                  engine_hashes=dict(hashes))
        self.worker = Worker(worker_config, client=InProcessWorkerClient(service))

    def run_once(self):
        # Worker performs real telemetry, strict cohort claim admission and
        # owned subprocess/lease accounting, exactly as remote workers do.
        return self.worker.run_once()

    def shutdown(self):
        self.worker.shutdown()

    def status(self):
        return dict(self.readiness)
