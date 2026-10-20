"""Portable pose references bound to the exact SDF bytes saved by a run."""
from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any, Optional


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def describe(root: Path, path: Optional[Path]) -> Optional[dict[str, Any]]:
    if path is None:
        return None
    target = Path(path).resolve()
    relative = target.relative_to(Path(root).resolve()).as_posix()
    if not target.is_file():
        raise ValueError("pose artifact is not a file")
    return {"path": relative, "sha256": sha256(target),
            "size_bytes": target.stat().st_size}


def resolve(root: Path, raw_path: Any, entry: Any) -> Optional[Path]:
    if raw_path is None and entry is None:
        return None
    if not isinstance(raw_path, str) or not isinstance(entry, dict):
        raise ValueError("pose has no artifact descriptor")
    relative = Path(raw_path)
    if relative.is_absolute() or ".." in relative.parts or entry.get("path") != raw_path:
        raise ValueError("pose artifact path is not within its package")
    target = (Path(root) / relative).resolve()
    if not target.is_relative_to(Path(root).resolve()) or not target.is_file():
        raise ValueError("pose artifact is missing or escapes its package")
    if (type(entry.get("size_bytes")) is not int
            or entry["size_bytes"] != target.stat().st_size
            or entry.get("sha256") != sha256(target)):
        raise ValueError("pose artifact bytes changed")
    return target
