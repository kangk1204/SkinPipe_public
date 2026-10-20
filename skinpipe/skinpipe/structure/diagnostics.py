"""Flushed stage receipts that survive a killed structure-preparation process."""
from __future__ import annotations

from contextlib import contextmanager
from functools import wraps
import hashlib
import inspect
import json
import os
from pathlib import Path
import time
from uuid import uuid4


class PreparationJournal:
    def __init__(self, source: Path, output: Path, settings: dict):
        self.path = Path(output).with_name(Path(output).name + ".preparation." + uuid4().hex + ".json")
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.value = {"schema_version": 1, "input_sha256": hashlib.sha256(Path(source).read_bytes()).hexdigest(),
                      "execution_run_id": os.environ.get("SKINPIPE_PREPARATION_RUN_ID"),
                      "settings": settings, "started_at": time.time(), "state": "RUNNING", "stages": []}
        self.flush()


    def flush(self):
        temp = self.path.with_suffix(".tmp")
        temp.write_text(json.dumps(self.value, indent=2) + "\n")
        temp.replace(self.path)

    @contextmanager
    def stage(self, name: str):
        row = {"name": name, "started_at": time.time(), "state": "RUNNING"}
        self.value["stages"].append(row)
        self.flush()
        start = time.monotonic()
        try:
            yield row
        except BaseException as exc:
            row.update(state="FAILED", error_type=type(exc).__name__, error=str(exc))
            self.value["state"] = "FAILED"
            raise
        else:
            row["state"] = "COMPLETED"
        finally:
            row.update(seconds=time.monotonic() - start, finished_at=time.time())
            self.flush()

    def finish(self, output: Path):
        self.value.update(state="COMPLETED", finished_at=time.time(),
                          output_sha256=hashlib.sha256(Path(output).read_bytes()).hexdigest())
        self.flush()


def record_attempt(function):
    """Keep rejected identities/chemistry and interrupted attempts inspectable."""
    signature = inspect.signature(function)

    @wraps(function)
    def wrapped(*args, **kwargs):
        arguments = signature.bind(*args, **kwargs)
        arguments.apply_defaults()
        values = arguments.arguments
        directory = Path(values["workdir"]) / "preparation_attempts"
        directory.mkdir(parents=True, exist_ok=True)
        path = directory / (uuid4().hex + ".json")
        receipt = {"schema_version": 1, "state": "RUNNING", "started_at": time.time(),
                   "execution_run_id": os.environ.get("SKINPIPE_PREPARATION_RUN_ID"),
                   "parameters": {k: v for k, v in values.items() if k != "workdir"},
                   "docking_validated": False, "functional_admissibility": "NOT_ESTABLISHED"}

        def write():
            temp = path.with_suffix(".tmp")
            temp.write_text(json.dumps(receipt, indent=2) + "\n")
            temp.replace(path)

        write()
        try:
            receptor = function(*args, **kwargs)
        except Exception as exc:
            error_type = type(exc).__name__
            message = str(exc)
            category = ("TARGET_IDENTITY_UNRESOLVED" if error_type == "TargetIdentityError"
                        else "TARGET_STRUCTURE_INADMISSIBLE" if error_type == "TargetStructureError"
                        else "NATIVE_CHEMISTRY_UNSUPPORTED" if "native modified" in message
                        else "OBSERVED_ATOM_PRESERVATION_FAILED" if "observed heavy atom" in message
                        else "PREPARATION_FAILED")
            receipt.update(state="FAILED", failure_category=category, error_type=error_type,
                           error=message, finished_at=time.time())
            write()
            raise
        receipt.update(state="PREPARED_REQUIRING_DOCKING_VALIDATION", finished_at=time.time(),
                       receptor_pdbqt_sha256=hashlib.sha256(Path(receptor.pdbqt).read_bytes()).hexdigest())
        write()
        return receptor

    return wrapped
