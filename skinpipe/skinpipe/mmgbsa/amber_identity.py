"""Content identity of the AmberTools executables and antechamber tables.

Conda's package version does not detect an in-place binary or table replacement.
The digest excludes the installation directory, so relocating unchanged files
does not invalidate a record. Mixed tool prefixes, unavailable tables, and files
that change during inspection produce an unverified identity and cannot justify
a cache hit. This is a file identity, not a numerical determinism guarantee;
shared libraries, hardware and SQM's runtime choices are outside this manifest.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import shutil
from typing import Any

TOOLS = ("antechamber", "sqm", "am1bcc", "atomtype", "bondtype", "parmchk2")
_REQUIRED_TABLES = ("BCCPARM.DAT", "ATOMTYPE_BCC.DEF", "ATOMTYPE_SYBYL.DEF",
                    "APS.DAT", "PARMCHK.DAT")


def _signature(stat: os.stat_result) -> tuple[int, ...]:
    return (stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns,
            stat.st_ctime_ns, stat.st_mode)


def _file_identity(path: Path) -> dict[str, Any]:
    """Read once, rejecting an observed replacement or write during the read."""
    resolved = path.resolve(strict=True)
    with path.open("rb") as stream:
        before = os.fstat(stream.fileno())
        digest = hashlib.sha256()
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
        after = os.fstat(stream.fileno())
    if (path.resolve(strict=True) != resolved or
            _signature(before) != _signature(after) or
            _signature(after) != _signature(path.stat())):
        raise OSError("file changed during identity inspection")
    return {"bytes": after.st_size, "sha256": digest.hexdigest()}


def _paths() -> dict[str, Path]:
    tools: dict[str, Path] = {}
    for name in TOOLS:
        found = shutil.which(name)
        if found is None:
            raise OSError(f"missing-tool-{name}")
        tools[name] = Path(found).resolve(strict=True)
    prefixes = {path.parent.parent for path in tools.values()}
    if len(prefixes) != 1:
        raise OSError("mixed-tool-prefixes")
    root = next(iter(prefixes))
    amber_home = os.environ.get("AMBERHOME")
    if amber_home and Path(amber_home).resolve(strict=True) != root:
        # Executables can launch companions or read data under AMBERHOME.
        # An ambiguous combination must not be labelled as the PATH install.
        raise OSError("AMBERHOME-differs-from-tool-prefix")
    table_dir = root / "dat" / "antechamber"
    if not all((table_dir / name).is_file() for name in _REQUIRED_TABLES):
        raise OSError("missing-antechamber-tables")
    tables = {"dat/antechamber/" + p.relative_to(table_dir).as_posix(): p
              for p in table_dir.rglob("*") if p.is_file()}
    return {**{"bin/" + name: path for name, path in tools.items()}, **tables}


def ambertools_manifest() -> dict[str, Any]:
    """An inspectable file manifest and the portable value used by the cache.

    No memoization: a replacement under the same path must be noticed on the
    next cache lookup. Inspection failures contain no installation paths.
    """
    try:
        paths = _paths()
        stamps = {name: (path.resolve(strict=True), _signature(path.stat()))
                  for name, path in paths.items()}
        files = {name: _file_identity(path) for name, path in sorted(paths.items())}
        final_paths = _paths()
        if (paths != final_paths or
                stamps != {name: (path.resolve(strict=True), _signature(path.stat()))
                           for name, path in final_paths.items()}):
            raise OSError("installation-changed-during-inspection")
        contents = {"version": 1, "files": files}
        digest = hashlib.sha256(json.dumps(contents, sort_keys=True,
                                           separators=(",", ":")).encode()).hexdigest()
        return {**contents, "verified": True, "identity": "sha256:" + digest}
    except (OSError, RuntimeError, ValueError) as exc:
        # Error messages from pathlib can include absolute paths. Use a stable
        # generic status; explicit inventory failures are useful and path-free.
        message = str(exc)
        reasons = ("missing-tool-", "mixed-tool-prefixes",
                   "AMBERHOME-differs-from-tool-prefix", "missing-antechamber-tables",
                   "installation-changed-during-inspection",
                   "file changed during identity inspection")
        reason = message if any(message.startswith(x) for x in reasons) else "unreadable-files"
        return {"version": 1, "files": {}, "verified": False,
                "identity": "unverified:" + reason}


def ambertools_identity() -> str:
    return str(ambertools_manifest()["identity"])
