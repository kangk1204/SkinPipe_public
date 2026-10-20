"""Owned short-lived CPU stages keep heavy chemistry out of the coordinator.

Only allowlisted pure chemistry operations run here. No API credentials, queue
connections, user commands, or source database writes are passed to the child.
"""
from pathlib import Path
import json
import os
import subprocess
import sys
import time
import uuid

from .artifacts import ArtifactStore, atomic_bytes, canonical_json
from .engines import terminate_owned


def run_compute(operation, payload, root, *, stop, progress=lambda _: None):
    if operation not in {"generation", "assessment", "hopping"}:
        raise ValueError("Unknown isolated chemistry stage")
    stop()
    work = Path(root) / "research-compute" / uuid.uuid4().hex
    work.mkdir(parents=True)
    request, result, status = work / "request.json", work / "result.json", work / "progress.json"
    atomic_bytes(request, canonical_json({"operation": operation, "payload": payload}))
    environment = {key: os.environ[key] for key in
        ("PATH", "HOME", "TMPDIR", "LANG", "LC_ALL", "SYSTEMROOT", "WINDIR") if key in os.environ}
    environment.update(PYTHONPATH=str(Path(__file__).resolve().parents[2]),
        PYTHONNOUSERSITE="1", CUDA_VISIBLE_DEVICES="", OMP_NUM_THREADS="1",
        OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1")
    last = None
    with (work / "execution.log").open("wb") as log:
        process = subprocess.Popen([sys.executable, "-m", "skinpipe.platform.research_compute",
            str(request), str(result), str(status)], env=environment, stdout=log, stderr=log,
            start_new_session=True)
        try:
            while process.poll() is None:
                stop()
                if status.is_file():
                    raw = status.read_bytes()
                    if raw != last:
                        progress(json.loads(raw))
                        last = raw
                time.sleep(0.05)
            stop()
        finally:
            terminate_owned(process)
    if process.returncode or not result.is_file():
        raise RuntimeError("Isolated chemistry stage failed; owned execution log retained")
    response = json.loads(result.read_bytes())
    if response.get("status") != "complete":
        raise RuntimeError("Isolated chemistry stage produced no complete result")
    stop()
    return response["result"]


def execute(operation, payload, progress):
    if operation == "generation":
        from .research_generation import generate_research
        return generate_research(payload["seed_smiles"], payload["policy"], payload["target_id"])
    if operation == "assessment":
        from .catalog_assessment import assess_catalog
        return assess_catalog(payload["rows"], skin_context=payload.get("skin_context"), progress=progress)
    if operation == "hopping":
        from .catalog_hopping import catalog_hopping
        return catalog_hopping(payload["rows"], payload["reference_set"], ArtifactStore(payload["artifact_root"]),
            progress=progress, **payload["settings"])
    raise ValueError("Unknown isolated chemistry stage")


def main():
    request, result, status = map(Path, sys.argv[1:])
    body = json.loads(request.read_bytes())
    def progress(value):
        atomic_bytes(status, canonical_json(value))
    response = execute(body["operation"], body["payload"], progress)
    atomic_bytes(result, canonical_json({"status": "complete", "result": response}))


if __name__ == "__main__":
    main()
