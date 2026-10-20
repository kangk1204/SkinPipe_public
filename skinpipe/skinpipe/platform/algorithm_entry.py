"""Internal bounded chemistry entry point for isolated optional environments."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from .artifacts import atomic_bytes, canonical_json


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    source = Path(args.input)
    if source.is_symlink() or source.stat().st_size > 32 * 1024 * 1024:
        raise ValueError("Invalid owned chemistry request")
    payload = json.loads(source.read_text())
    from .algorithm_chemistry import execute_chemistry
    result = execute_chemistry(payload["provider"], payload["rows"], payload["context"], workdir=source.parent)
    atomic_bytes(Path(args.output), canonical_json(result))


if __name__ == "__main__":
    main()
