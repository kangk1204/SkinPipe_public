"""File boundary for the existing mapped-SMILES adaptive charge calculation.

Invoked by ``run_calculation`` so external producers share its cleanup scope.
The temporary JSON response also handles large vectors without a queue feeder.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
import sys


def main() -> None:
    from .ligand import _charge_values

    request, response = (Path(value) for value in sys.argv[1:])
    record = json.loads(request.read_text())
    if (not isinstance(record, dict) or not isinstance(record.get("smiles"), str)
            or not isinstance(record.get("method"), str)):
        raise ValueError("invalid adaptive charge request")
    status, payload = _charge_values(record["smiles"], record["method"])
    temporary = response.with_suffix(".part")
    temporary.write_text(json.dumps({"status": status, "payload": payload}, allow_nan=False))
    os.replace(temporary, response)


if __name__ == "__main__":
    main()
