"""Read-only acquisition checks for saved poses, separate from geometry checks.

A finite score is metadata. PASS here requires one readable 3D SDF with bytes
bound by its recorded portable descriptor. It certifies neither the geometry
of a binding pose nor a biological effect. Legacy coordinates remain inspectable
but cannot establish that the original scoring operation used those bytes.
"""
from __future__ import annotations

import hashlib
import io
import math
from pathlib import Path

from . import artifacts


def finite_number(value: object) -> bool:
    return (type(value) in (int, float) and math.isfinite(value))


def inspect_pose(dock: Path, pose: dict) -> dict:
    """Return an input receipt; absent, unbound or changed inputs are unknown."""
    receipt = {"status": "NOT_EVALUATED", "pose_index": pose.get("pose_index"),
               "binding": "unverified", "path": pose.get("sdf_path")}
    try:
        if type(pose.get("pose_index")) is not int or pose["pose_index"] < 1:
            raise ValueError("missing or invalid pose index")
        raw_path = pose.get("sdf_path")
        if "pose_artifact" in pose:
            path = artifacts.resolve(dock, raw_path, pose["pose_artifact"])
            binding = "portable_pose_bytes_verified"
        else:
            if not isinstance(raw_path, str) or not raw_path:
                raise ValueError("missing pose coordinates")
            path = Path(raw_path)
            if not path.is_absolute():
                path = dock / path
            path = path.resolve()
            if not path.is_relative_to(dock.resolve()) or not path.is_file():
                raise ValueError("legacy pose is missing or outside its docking package")
            binding = "legacy_current_bytes_only"
        if path is None:
            raise ValueError("missing pose coordinates")
        raw = path.read_bytes()
        digest = hashlib.sha256(raw).hexdigest()
        receipt.update(path=str(path), sha256=digest, size_bytes=len(raw), binding=binding)
        # Parse exactly the bytes that were hashed, rather than opening a second
        # path snapshot. A producer changing the file invalidates this receipt.
        entry = pose.get("pose_artifact")
        if binding == "portable_pose_bytes_verified" and (
                entry.get("sha256") != digest or entry.get("size_bytes") != len(raw)):
            raise ValueError("pose artifact bytes changed during acquisition")
        from rdkit import Chem

        supplier = iter(Chem.ForwardSDMolSupplier(io.BytesIO(raw), removeHs=False,
                                                  sanitize=True, strictParsing=True))
        mol = next(supplier, None)
        # Only a second record is needed to refuse an ambiguous multi-molecule
        # file; a large malformed input cannot make us retain all molecules.
        sentinel = object()
        if mol is None or next(supplier, sentinel) is not sentinel:
            raise ValueError("pose file must contain exactly one readable molecule")
        if not mol.GetNumAtoms() or mol.GetNumConformers() != 1:
            raise ValueError("pose has no unique coordinate conformer")
        conf = mol.GetConformer()
        if not conf.Is3D() or not all(math.isfinite(float(v))
                                    for xyz in conf.GetPositions() for v in xyz):
            raise ValueError("pose has no finite 3D coordinates")
        receipt.update(n_atoms=mol.GetNumAtoms(), coordinates_readable=True)
        if binding != "portable_pose_bytes_verified":
            raise ValueError("legacy coordinates have no recorded byte binding to the score")
        receipt["status"] = "PASS"
    except (OSError, ValueError, TypeError, RuntimeError, ImportError) as exc:
        receipt["reason"] = str(exc)
    return receipt


def inspect_record(dock: Path, record: dict) -> list[dict]:
    poses = record.get("poses") or []
    if record.get("_unreadable"):
        return []
    indices = [p.get("pose_index") for p in poses]
    paths = [p.get("sdf_path") for p in poses]
    ambiguous = (any(type(i) is not int or i < 1 for i in indices)
                 or len(set(indices)) != len(indices)
                 or any(not isinstance(p, str) or not p for p in paths)
                 or len(set(paths)) != len(paths))
    if ambiguous:
        return [{"status": "NOT_EVALUATED", "pose_index": p.get("pose_index"),
                 "reason": "pose indices or coordinate paths are missing or ambiguous"}
                for p in poses]
    return [inspect_pose(dock, p) for p in poses]
