"Offline structure lookup from ChEMBL's bulk `chemreps` file."
from __future__ import annotations

import gzip
import hashlib
import json
import logging
import os
import sqlite3
import time
from pathlib import Path
from typing import Iterable, Optional

log = logging.getLogger(__name__)

#: Override for both the archive and, by extension, the index beside it.
ENV_VAR = "SKINPIPE_CHEMBL_CHEMREPS"
_DEFAULT = Path.home() / "chembl_local" / "chembl_37_chemreps.txt.gz"

#: Rows per transaction while indexing. Large enough that the commit overhead
#: disappears against 2.9M rows, small enough that a kill mid-build does not
#: leave the work of an hour uncommitted.
_CHUNK = 50_000

#: How much of the archive the fingerprint actually hashes. Hashing all 279 MB
#: costs real CPU and would be paid on every `lookup`; 8 MiB of the compressed
#: stream costs milliseconds and still changes whenever the file is rewritten,
#: because gzip's own header (mtime, name) sits in the first few bytes and its
#: deflate blocks diverge from the first differing row.
_PREFIX_BYTES = 8 << 20

#: Named in the stored metadata so a future scheme change is recognisable
#: rather than silently compared against digests it cannot match.
_SCHEME = "size+mtime_ns+blake2b16(first 8MiB)"

#: Where the fingerprint lives inside the index itself, so the two cannot be
#: separated by copying one file and not the other.
_META = "index_meta"


def archive_path() -> Path:
    """The bulk file to index, from the environment or the default location."""
    raw = os.environ.get(ENV_VAR)
    return Path(raw).expanduser() if raw else _DEFAULT


def index_path(archive: Optional[Path] = None) -> Path:
    """The SQLite index that `lookup` reads, kept beside its archive."""
    a = archive or archive_path()
    return a.with_name(a.name.split(".")[0] + ".sqlite")


def fingerprint(archive: Optional[Path] = None) -> Optional[dict]:
    """A cheap identity for the archive an index was built from.

    Deliberately *not* a full content hash. Rehashing 279 MB is about a minute
    of CPU, and the whole point of this module is to be the cheap path when the
    API is down; paying that on every lookup would make the fallback worse than
    the outage. What this does catch is every way the file realistically
    changes: a different release (size and prefix differ), a re-download (mtime
    and gzip header differ), a truncation (size differs), a rewrite of the same
    length at the same timestamp (prefix digest differs).

    What it does not catch is a rewrite that keeps the byte length, the
    nanosecond mtime *and* the first 8 MiB identical while changing a row
    further in. That is not an accident anyone has; it is a forgery, and this
    is not a tamper-proofing mechanism. Returns None when there is no archive.
    """
    a = Path(archive) if archive else archive_path()
    try:
        st = a.stat()
        h = hashlib.blake2b(digest_size=16)
        with a.open("rb") as fh:
            read = fh.read(_PREFIX_BYTES)
        h.update(read)
    except OSError:
        return None
    return {"scheme": _SCHEME, "name": a.name, "size": st.st_size,
            "mtime_ns": st.st_mtime_ns, "prefix_bytes": len(read),
            "digest": h.hexdigest()}


def _health(index: Path) -> str:
    """`missing` / `unreadable` / `empty` / `ok` for the index file itself.

    Deliberately `SELECT 1 ... LIMIT 1` and not `count(*)`: on 2.9M rows the
    count is a full table scan, and "does this hold anything at all" is the
    question being asked.
    """
    if not index.exists():
        return "missing"
    try:
        with sqlite3.connect(f"file:{index}?mode=ro", uri=True) as db:
            row = db.execute("SELECT 1 FROM structures LIMIT 1").fetchone()
    except sqlite3.Error:
        return "unreadable"
    return "ok" if row is not None else "empty"


def _read_meta(index: Path) -> dict:
    try:
        with sqlite3.connect(f"file:{index}?mode=ro", uri=True) as db:
            rows = db.execute(f"SELECT key, value FROM {_META}").fetchall()
    except sqlite3.Error:
        return {}
    return {k: v for k, v in rows}


def _write_meta(db: sqlite3.Connection, fp: Optional[dict],
                rows: Optional[int] = None) -> None:
    db.execute(f"CREATE TABLE IF NOT EXISTS {_META} (key TEXT PRIMARY KEY, "
               "value TEXT)")
    items = {"scheme": _SCHEME, "stamped_at": str(int(time.time()))}
    if fp is not None:
        items["archive_fingerprint"] = json.dumps(fp, sort_keys=True)
    if rows is not None:
        items["rows"] = str(rows)
    db.executemany(f"INSERT OR REPLACE INTO {_META} VALUES (?,?)",
                   sorted(items.items()))


def verify_index(archive: Optional[Path] = None) -> dict:
    """Is the index beside this archive readable, populated, and *its* index?

    The status is the useful part, because the three failures want three
    different responses:

    - `ok`          serve it.
    - `missing`     nothing to serve; build one.
    - `unreadable`  the file is not a database; rebuild over it.
    - `empty`       readable, schema-correct, zero structures; rebuild over it.
    - `mismatch`    built from a different archive than the one beside it now.
                    Refuse to serve, rebuild on request.
    - `unstamped`   built before fingerprints existed. Nothing is known about
                    it either way, so it is neither trusted nor rebuilt
                    unasked -- see `stamp_index`.
    - `no_archive`  the index is intact but its archive is gone, so there is
                    nothing to compare against.

    Never raises: this sits under a fallback, and a fallback that raises turns
    a partial screen into no screen.
    """
    a = Path(archive) if archive else archive_path()
    out = index_path(a)
    base = {"path": str(out), "archive": str(a), "scheme": _SCHEME}
    health = _health(out)
    if health != "ok":
        return {**base, "ok": False, "status": health, "rows": None,
                "recorded": None, "current": None,
                "reason": {"missing": f"no index at {out}",
                           "unreadable": f"the file at {out} is not a readable "
                                         "SQLite index",
                           "empty": f"the index at {out} holds no structures"}[health]}

    meta = _read_meta(out)
    raw = meta.get("archive_fingerprint")
    try:
        recorded = json.loads(raw) if raw else None
    except ValueError:
        recorded = None
    rows = int(meta["rows"]) if (meta.get("rows") or "").isdigit() else None
    current = fingerprint(a)

    if recorded is None:
        return {**base, "ok": False, "status": "unstamped", "rows": rows,
                "recorded": None, "current": current,
                "reason": (f"the index at {out} carries no record of the archive it "
                           "was built from, so it cannot be told apart from one "
                           "built from a different file. Re-validate it with "
                           "chembl_local.stamp_index() if you know it matches the "
                           "archive beside it, or rebuild with "
                           "build_index(force=True).")}
    if current is None:
        return {**base, "ok": False, "status": "no_archive", "rows": rows,
                "recorded": recorded, "current": None,
                "reason": f"no archive at {a} to check the index against"}

    differs = [k for k in ("size", "mtime_ns", "digest")
               if recorded.get(k) != current.get(k)]
    if differs:
        return {**base, "ok": False, "status": "mismatch", "rows": rows,
                "recorded": recorded, "current": current, "differs": differs,
                "reason": (f"the archive at {a} is not the file the index at {out} "
                           f"was built from ({', '.join(differs)} differ)")}
    return {**base, "ok": True, "status": "ok", "rows": rows,
            "recorded": recorded, "current": current, "differs": [],
            "reason": "the index is readable, populated and matches its archive"}


def stamp_index(archive: Optional[Path] = None) -> dict:
    """Record the current archive's fingerprint onto an existing index.

    This is the cheap way to re-validate an index built before fingerprinting
    existed, without spending minutes of CPU rebuilding one that is already
    correct. It writes three rows into the index; it does not read the archive
    beyond the first 8 MiB and does not touch the structures table.

    Be clear about what it means: this is an *operator's assertion* that the
    index beside this archive was built from it. Nothing here can prove that --
    if it could, the fingerprint would not be needed. Use it when you know the
    history of the two files; otherwise rebuild.
    """
    a = Path(archive) if archive else archive_path()
    out = index_path(a)
    health = _health(out)
    if health != "ok":
        return {"ok": False, "status": health, "path": str(out),
                "reason": f"no usable index at {out} to stamp ({health})"}
    fp = fingerprint(a)
    if fp is None:
        return {"ok": False, "status": "no_archive", "path": str(out),
                "reason": f"no archive at {a} to take a fingerprint of"}
    try:
        with sqlite3.connect(out) as db:
            _write_meta(db, fp)
    except sqlite3.Error as e:
        return {"ok": False, "status": "unwritable", "path": str(out),
                "reason": f"{type(e).__name__}: {e}"}
    log.info("chembl_local: stamped %s as built from %s", out, a)
    return {"ok": True, "status": "stamped", "path": str(out),
            "archive": str(a), "fingerprint": fp,
            "reason": "the index now records the archive beside it"}


#: Statuses `lookup` will still serve from: nothing is known to be wrong, only
#: unproven. Refusing these would silently disable the fallback on every index
#: built before this check existed, which trades a possible staleness for a
#: certain loss of compounds.
_SERVE_UNPROVEN = ("unstamped", "no_archive")


def available() -> bool:
    """True when a usable index exists. Never raises, never builds one.

    "Usable" now means readable, populated, and not known to disagree with the
    archive beside it. An index whose archive was swapped is no longer usable.
    """
    v = verify_index()
    if v["ok"]:
        return True
    if v["status"] in _SERVE_UNPROVEN:
        log.warning("chembl_local: %s", v["reason"])
        return True
    return False


def build_index(archive: Optional[Path] = None, *, force: bool = False,
                progress_every: int = 500_000) -> dict:
    """Index the archive into SQLite, returning a status dict rather than raising.

    The build writes to a temporary file and renames it into place, so an
    interrupted build never leaves a half-populated index that `lookup` would
    read as authoritative. Re-running without `force` is a no-op once a
    *verified* index exists, because this takes minutes and is not something to
    repeat on import.

    "Verified" is the part finding N-20 was about. The old test was
    `out.exists()`, so a swapped archive, a truncated index and an empty one
    were all reported ok. Now `verify_index` decides, and each answer gets the
    response it deserves: an index that is broken, empty or built from another
    archive is rebuilt over (the archive is right there); an index with no
    recorded fingerprint is refused rather than rebuilt, because rebuilding
    unasked is minutes of CPU on a file that is probably fine -- `stamp_index`
    is the cheap way out and the returned reason says so.
    """
    a = Path(archive) if archive else archive_path()
    out = index_path(a)
    if out.exists() and not force:
        v = verify_index(a)
        if v["ok"]:
            return {"ok": True, "built": False, "verified": True,
                    "status": "ok", "path": str(out), "rows": v["rows"],
                    "reason": "index already present and verified against its archive"}
        if v["status"] in ("unstamped", "no_archive"):
            log.warning("chembl_local: %s", v["reason"])
            return {"ok": False, "built": False, "verified": False,
                    "status": v["status"], "path": str(out), "rows": v["rows"],
                    "reason": v["reason"]}
        log.warning("chembl_local: %s; rebuilding", v["reason"])
    if not a.exists():
        return {"ok": False, "built": False, "verified": False,
                "status": "no_archive", "path": str(out),
                "reason": f"no archive at {a}"}

    # Taken before the read, so a file rewritten *during* the build fails the
    # next verification instead of being certified by a fingerprint of its
    # final state.
    fp = fingerprint(a)
    rebuilding = out.exists()
    tmp = out.with_suffix(".sqlite.partial")
    tmp.unlink(missing_ok=True)
    n = 0
    try:
        with sqlite3.connect(tmp) as db:
            db.execute("PRAGMA journal_mode=OFF")
            db.execute("PRAGMA synchronous=OFF")
            db.execute("CREATE TABLE structures ("
                       "chembl_id TEXT PRIMARY KEY, smiles TEXT, inchi_key TEXT)")
            with gzip.open(a, "rt", encoding="utf-8", errors="replace") as fh:
                header = fh.readline()
                if "chembl_id" not in header:
                    raise ValueError(f"unexpected header: {header[:80]!r}")
                batch: list[tuple[str, str, str]] = []
                for line in fh:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 4 or not parts[1]:
                        continue
                    batch.append((parts[0], parts[1], parts[3]))
                    if len(batch) >= _CHUNK:
                        db.executemany("INSERT OR REPLACE INTO structures VALUES (?,?,?)",
                                       batch)
                        db.commit()
                        n += len(batch)
                        batch.clear()
                        if progress_every and n % progress_every == 0:
                            log.info("chembl_local: indexed %d rows", n)
                if batch:
                    db.executemany("INSERT OR REPLACE INTO structures VALUES (?,?,?)",
                                   batch)
                    db.commit()
                    n += len(batch)
            _write_meta(db, fp, rows=n)
    except (OSError, sqlite3.Error, ValueError, EOFError) as e:
        tmp.unlink(missing_ok=True)
        log.error("chembl_local: index build failed: %s", e)
        return {"ok": False, "built": False, "verified": False,
                "status": "build_failed", "path": str(out),
                "reason": f"{type(e).__name__}: {e}", "rows": n}

    tmp.replace(out)
    log.info("chembl_local: indexed %d structures into %s", n, out)
    return {"ok": True, "built": True, "verified": True,
            "status": "rebuilt" if rebuilding else "built",
            "path": str(out), "rows": n, "fingerprint": fp}


def lookup(chembl_ids: Iterable[str]) -> dict[str, dict]:
    """SMILES and InChI key for the ids this index holds.

    Returns only what it found, so a caller can diff against what it asked for.
    An absent or unreadable index is an empty result and a log line, never an
    exception: this is a fallback path, and a fallback that raises turns a
    partial screen into no screen.

    An index that disagrees with the archive beside it is also an empty result.
    Serving a structure out of a superseded release under the current release's
    id is the one thing worse than losing the compound, because the loss is
    counted in `missing_ids` and the wrong structure is not counted anywhere.
    The verification costs one stat and 8 MiB of hashing per call.
    """
    ids = [i for i in dict.fromkeys(chembl_ids) if i]
    if not ids:
        return {}
    v = verify_index()
    if not v["ok"]:
        if v["status"] not in _SERVE_UNPROVEN:
            if v["status"] != "missing":          # "there is none" is not news
                log.warning("chembl_local: refusing this index: %s", v["reason"])
            return {}
        log.warning("chembl_local: serving an index that could not be matched to "
                    "its archive: %s", v["reason"])
    # Serving an unproven index is a deliberate policy -- an index with no archive
    # beside it is usually fine and rebuilding costs minutes -- but the record has
    # to say so. It carried `source="chembl_local"` and nothing else, so a
    # structure rescued from a verified index and one from an index nobody could
    # match to its archive reached the report as the same claim. The compound is
    # still usable; what changes is how much the provenance line is worth.
    _verif = {"index_verified": bool(v["ok"]), "index_status": v["status"]}
    p = index_path()
    out: dict[str, dict] = {}
    try:
        with sqlite3.connect(f"file:{p}?mode=ro", uri=True) as db:
            for i in range(0, len(ids), 900):          # under SQLite's var limit
                chunk = ids[i:i + 900]
                q = ("SELECT chembl_id, smiles, inchi_key FROM structures "
                     f"WHERE chembl_id IN ({','.join('?' * len(chunk))})")
                for cid, smi, key in db.execute(q, chunk):
                    out[cid] = {"smiles": smi, "inchi_key": key,
                                "source": "chembl_local", **_verif}
    except sqlite3.Error as e:
        log.warning("chembl_local: lookup failed, continuing without it: %s", e)
        return {}
    return out


def as_molecule_records(found: dict[str, dict]) -> dict[str, dict]:
    """Shape local hits like `/molecule` records, with the gaps left explicit.

    Every field the archive cannot supply is set to None rather than omitted,
    so a downstream `.get("max_phase")` sees "not known" instead of raising,
    and `partial` marks the record for any report that wants to say where the
    structure came from.
    """
    return {cid: {"name": None, "max_phase": None, "natural_product": None,
                  "molecule_type": None, "smiles": v["smiles"],
                  "inchi_key": v.get("inchi_key"),
                  "mw": None, "logp": None, "tpsa": None, "hbd": None,
                  "hba": None, "ro5_violations": None, "aromatic_rings": None,
                  "heavy_atoms": None,
                  "source": "chembl_local", "partial": True,
                  # Carried through the reshape. `lookup` records whether the
                  # index could be matched to its archive, and this function
                  # built a fresh dict that dropped it -- so the branch in
                  # `library.from_chembl` that reports it never ran, and no
                  # candidate, CSV or JSON could say which kind of index a
                  # rescued structure came from.
                  "index_verified": v.get("index_verified"),
                  "index_status": v.get("index_status")}
            for cid, v in found.items()}
