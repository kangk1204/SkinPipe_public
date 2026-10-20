"""Process-isolated preparation deadlines, preserving failed artifacts."""
from __future__ import annotations

from dataclasses import asdict
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time
from uuid import uuid4

from ..processes import run_calculation
from ..types import Receptor


class PreparationTimeoutError(RuntimeError):
    """A receptor attempt exceeded its explicit wall-time budget."""


def _write(path, value):
    temp = path.with_suffix(".tmp")
    temp.write_text(json.dumps(value, indent=2, default=str) + "\n")
    temp.replace(path)


def _close_interrupted(workdir, run_id, state):
    """Update only this child's receipts; invalidate its late cache card."""
    paths = list(workdir.glob("preparation_attempts/*.json"))
    paths += list(workdir.glob("receptors/*/*.preparation.*.json"))
    paths += list(workdir.glob("receptors/*/receptor_identity.json"))
    for path in paths:
        try:
            value = json.loads(path.read_text())
        except (OSError, ValueError):
            continue
        if value.get("execution_run_id") != run_id:
            continue
        if path.name == "receptor_identity.json":
            value["preparation_interrupted"] = True
        else:
            value.update(state=state, finished_at=time.time())
            for stage in value.get("stages", []):
                if stage.get("state") == "RUNNING":
                    stage.update(state=state, finished_at=time.time())
        _write(path, value)


def run_preparation(pdb_id, workdir, settings, timeout):
    if (isinstance(timeout, bool) or not isinstance(timeout, (int, float))
            or not math.isfinite(timeout) or not 0 < timeout <= 86400):
        raise ValueError("preparation_timeout must be finite, > 0 and <= 86400 seconds")
    workdir = Path(workdir).resolve()
    run_id = uuid4().hex
    directory = workdir / "bounded_preparations" / run_id
    directory.mkdir(parents=True)
    request = directory / "request.json"
    response = directory / "response.json"
    receipt_path = directory / "deadline.json"
    _write(request, {"pdb_id": pdb_id, "workdir": str(workdir), "settings": settings})
    receipt = {"schema_version": 1, "execution_run_id": run_id, "pdb_id": pdb_id,
               "state": "RUNNING", "started_at": time.time(), "timeout_seconds": timeout,
               "docking_validated": False, "functional_admissibility": "NOT_ESTABLISHED"}
    _write(receipt_path, receipt)
    environment = dict(os.environ, SKINPIPE_PREPARATION_RUN_ID=run_id)
    try:
        with (directory / "execution.log").open("w") as log:
            process = run_calculation([sys.executable, "-m", "skinpipe.structure.bounded",
                                       str(request), str(response)], timeout=timeout,
                                      env=environment, stdout=log, stderr=subprocess.STDOUT)
        if not response.is_file():
            raise RuntimeError(f"preparation worker exited {process.returncode} without a response; {directory}")
        value = json.loads(response.read_text())
        if process.returncode or value.get("state") != "PREPARED":
            from . import identity, prepare, rcsb
            errors = {c.__name__: c for c in (ValueError, RuntimeError, prepare.TargetIdentityError,
                                             identity.TargetStructureError, rcsb.RCSBError)}
            error_type = value.get("error_type", "RuntimeError")
            raise errors.get(error_type, RuntimeError)(f"{error_type}: {value.get('error', 'worker failed')}")
        fields = value["receptor"]
        for key in ("path", "pdbqt", "cocrystal_sdf"):
            if fields.get(key) is not None:
                fields[key] = Path(fields[key])
        for key in ("center", "size", "cocrystal_atoms"):
            if fields.get(key) is not None:
                fields[key] = tuple(fields[key])
        receptor = Receptor(**fields)
        receipt.update(state="PREPARED_REQUIRING_DOCKING_VALIDATION", finished_at=time.time())
        return receptor
    except subprocess.TimeoutExpired as error:
        _close_interrupted(workdir, run_id, "TIMED_OUT")
        receipt.update(state="TIMED_OUT", finished_at=time.time())
        raise PreparationTimeoutError(f"{pdb_id}: preparation timeout after {timeout:g}s; receipt {receipt_path}") from error
    except BaseException as error:
        _close_interrupted(workdir, run_id, "INTERRUPTED" if not isinstance(error, Exception) else "FAILED")
        receipt.update(state="FAILED", error_type=type(error).__name__, error=str(error), finished_at=time.time())
        raise
    finally:
        _write(receipt_path, receipt)


def main():
    request, response = map(Path, sys.argv[1:])
    value = json.loads(request.read_text())
    from .prepare import prepare_receptor
    try:
        receptor = prepare_receptor(value["pdb_id"], Path(value["workdir"]), **value["settings"])
    except Exception as error:
        _write(response, {"state": "FAILED", "error_type": type(error).__name__, "error": str(error)})
        raise
    _write(response, {"state": "PREPARED", "receptor": asdict(receptor)})


if __name__ == "__main__":
    main()
