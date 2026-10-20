'Offline natural-product search from the LOTUS bulk export.'
from __future__ import annotations

import csv
import gzip
import hashlib
import math
import json
import logging
import os
import sqlite3
from pathlib import Path
from typing import Iterable, Optional

log = logging.getLogger(__name__)

ENV_VAR = "SKINPIPE_LOTUS_METADATA"
_DEFAULT = Path.home() / "np_local" / "lotus_metadata.csv.gz"

#: Matches `chembl_local`: hash a prefix rather than the whole file. Rehashing
#: 86 MB on every call would defeat the point of a fast local path.
_PREFIX_BYTES = 8 * 1024 * 1024
_SCHEME = "lotus-prefix-blake2b-16"
_FP_KEY = "archive_fingerprint"

#: Morgan radius 2 / 2048 bits, the same fingerprint `catalog.screen` ranks on,
#: so a natural-product candidate and a purchasable one are scored on one scale.
FP_RADIUS = 2
FP_BITS = 2048


def archive_path() -> Path:
    raw = os.environ.get(ENV_VAR)
    return Path(raw).expanduser() if raw else _DEFAULT


def index_path(archive: Optional[Path] = None) -> Path:
    a = archive or archive_path()
    return a.with_name(a.name.split(".")[0] + ".sqlite")


def fingerprint(archive: Optional[Path] = None) -> Optional[dict]:
    a = Path(archive) if archive else archive_path()
    if not a.exists():
        return None
    st = a.stat()
    h = hashlib.blake2b(digest_size=16)
    with a.open("rb") as fh:
        read = fh.read(_PREFIX_BYTES)
    h.update(read)
    return {"scheme": _SCHEME, "name": a.name, "size": st.st_size,
            "mtime_ns": st.st_mtime_ns, "prefix_bytes": len(read),
            "prefix_digest": h.hexdigest()}


def _meta(conn: sqlite3.Connection) -> dict:
    try:
        return {k: v for k, v in conn.execute("SELECT key, value FROM meta")}
    except sqlite3.Error:
        return {}


def verify_index(archive: Optional[Path] = None) -> dict:
    a = Path(archive) if archive else archive_path()
    out = index_path(a)
    if not out.exists():
        return {"ok": False, "status": "absent", "rows": 0,
                "reason": f"no index at {out}"}
    try:
        with sqlite3.connect(f"file:{out}?mode=ro", uri=True) as conn:
            meta = _meta(conn)
            rows = conn.execute("SELECT COUNT(*) FROM np").fetchone()[0]
    except sqlite3.Error as e:
        return {"ok": False, "status": "broken", "rows": 0, "reason": str(e)}
    if not rows:
        return {"ok": False, "status": "empty", "rows": 0,
                "reason": "index has no rows"}
    raw = meta.get(_FP_KEY)
    if not raw:
        return {"ok": False, "status": "unstamped", "rows": rows,
                "reason": "index carries no archive fingerprint; run stamp_index "
                          "if it is known good, or build_index(force=True)"}
    if not a.exists():
        return {"ok": False, "status": "no_archive", "rows": rows,
                "reason": f"no archive at {a} to verify against"}
    try:
        stored = json.loads(raw)
    except (ValueError, TypeError):
        return {"ok": False, "status": "broken", "rows": rows,
                "reason": "the stored fingerprint is not readable JSON; the "
                          "index cannot be matched to an archive and has to be "
                          "rebuilt"}
    if stored != fingerprint(a):
        return {"ok": False, "status": "stale", "rows": rows,
                "reason": "index was built from a different archive than the one present"}
    return {"ok": True, "status": "ok", "rows": rows,
            "reason": "index verified against its archive"}


def stamp_index(archive: Optional[Path] = None) -> dict:
    """Record the current archive's fingerprint on an existing index.

    The cheap way out of `unstamped`: rebuilding costs minutes of fingerprint
    CPU on a file that is probably fine. Only ever call this when the index is
    known to have come from the archive now present.
    """
    a = Path(archive) if archive else archive_path()
    out = index_path(a)
    fp = fingerprint(a)
    if fp is None:
        return {"ok": False, "reason": f"no archive at {a}"}
    if not out.exists():
        return {"ok": False, "reason": f"no index at {out}"}
    with sqlite3.connect(out) as conn:
        conn.execute("CREATE TABLE IF NOT EXISTS meta (key TEXT PRIMARY KEY, value TEXT)")
        conn.execute("INSERT OR REPLACE INTO meta VALUES (?, ?)",
                     (_FP_KEY, json.dumps(fp, sort_keys=True)))
    return {"ok": True, "reason": "stamped"}


def _row_stream(archive: Path):
    """Unique structures from the export, organisms aggregated per structure.

    The export is one row per (structure, organism, reference), so 709,741 rows
    collapse to 227,298 structures. Aggregating here rather than at query time
    keeps the index one row per molecule, which is what a similarity scan wants.
    """
    seen: dict[str, dict] = {}
    with gzip.open(archive, "rt", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh):
            key = (row.get("structure_inchikey") or "").strip()
            smi = (row.get("structure_smiles_2D") or row.get("structure_smiles") or "").strip()
            if not key or not smi:
                continue
            rec = seen.get(key)
            if rec is None:
                rec = seen[key] = {
                    "inchikey": key, "smiles": smi,
                    "mass": row.get("structure_exact_mass") or "",
                    "xlogp": row.get("structure_xlogp") or "",
                    "pathway": row.get("structure_taxonomy_npclassifier_01pathway") or "",
                    "superclass": row.get("structure_taxonomy_npclassifier_02superclass") or "",
                    "npclass": row.get("structure_taxonomy_npclassifier_03class") or "",
                    "name": row.get("structure_nameTraditional") or "",
                    "organisms": set(),
                }
            org = (row.get("organism_name") or "").strip()
            if org:
                rec["organisms"].add(org)
    return seen


def build_index(archive: Optional[Path] = None, *, force: bool = False,
                progress_every: int = 25_000) -> dict:
    'Index the export into SQLite with a Morgan fingerprint per structure.'
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdFingerprintGenerator
    RDLogger.DisableLog("rdApp.*")

    a = Path(archive) if archive else archive_path()
    out = index_path(a)
    if out.exists() and not force:
        v = verify_index(a)
        if v["ok"]:
            return {"ok": True, "built": False, "rows": v["rows"],
                    "path": str(out), "reason": v["reason"]}
        if v["status"] in ("unstamped", "no_archive"):
            return {"ok": False, "built": False, "rows": v["rows"],
                    "path": str(out), "reason": v["reason"]}
        log.warning("lotus_local: %s; rebuilding", v["reason"])
    if not a.exists():
        return {"ok": False, "built": False, "rows": 0, "path": str(out),
                "reason": f"no archive at {a}"}

    fp = fingerprint(a)          # before the read, as chembl_local does
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=FP_RADIUS, fpSize=FP_BITS)
    records = _row_stream(a)
    # Unique per process: every builder used the same path and
    # unlinked it first, so two concurrent builds of one archive would
    # delete each other's work-in-progress and one could rename a
    # half-written database over the real index.
    _sweep_stale_temps(out)
    tmp = out.with_suffix(f".sqlite.tmp.{os.getpid()}")
    tmp.unlink(missing_ok=True)

    n = skipped = 0
    with sqlite3.connect(tmp) as conn:
        conn.execute("""CREATE TABLE np (
            inchikey TEXT PRIMARY KEY, smiles TEXT NOT NULL, name TEXT,
            mass REAL, xlogp REAL, pathway TEXT, superclass TEXT, npclass TEXT,
            n_organisms INTEGER, organisms TEXT, fp BLOB NOT NULL)""")
        conn.execute("CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT)")
        for key, rec in records.items():
            mol = Chem.MolFromSmiles(rec["smiles"])
            if mol is None:
                skipped += 1
                continue
            bits = gen.GetFingerprint(mol).ToBitString().encode("ascii")
            orgs = sorted(rec["organisms"])
            conn.execute("INSERT OR REPLACE INTO np VALUES (?,?,?,?,?,?,?,?,?,?,?)",
                         (key, rec["smiles"], rec["name"] or None,
                          float(rec["mass"]) if rec["mass"] else None,
                          float(rec["xlogp"]) if rec["xlogp"] else None,
                          rec["pathway"] or None, rec["superclass"] or None,
                          rec["npclass"] or None, len(orgs),
                          json.dumps(orgs[:40], ensure_ascii=False), bits))
            n += 1
            if progress_every and n % progress_every == 0:
                log.info("lotus_local: indexed %d structures", n)
        conn.execute("CREATE INDEX np_pathway ON np(pathway)")
        conn.execute("INSERT OR REPLACE INTO meta VALUES (?, ?)",
                     (_FP_KEY, json.dumps(fp, sort_keys=True)))
        conn.execute("INSERT OR REPLACE INTO meta VALUES ('release', ?)", (a.name,))
        conn.execute("INSERT OR REPLACE INTO meta VALUES ('unparseable', ?)", (str(skipped),))
    tmp.replace(out)
    log.info("lotus_local: indexed %d structures (%d unparseable) into %s", n, skipped, out)
    return {"ok": True, "built": True, "rows": n, "unparseable": skipped,
            "path": str(out), "reason": "built"}


def _sweep_stale_temps(out: Path) -> None:
    """Remove build temporaries left by processes that are gone.

    The build writes to a per-process name so two concurrent builds cannot
    delete each other's work. The cost of that is debris: a crash used to leave
    one file at a shared path that the next build unlinked, and now it leaves
    one nobody claims. So the shared legacy name is removed unconditionally, and
    a `.tmp.<pid>` file is removed only when no such process is alive -- which is
    what makes this safe to call while another build is running.
    """
    legacy = out.with_suffix(".sqlite.tmp")
    legacy.unlink(missing_ok=True)
    for stale in out.parent.glob(out.with_suffix(".sqlite.tmp.*").name):
        suffix = stale.name.rsplit(".", 1)[-1]
        if not suffix.isdigit():
            continue
        try:
            os.kill(int(suffix), 0)
        except ProcessLookupError:
            stale.unlink(missing_ok=True)          # its builder is gone
        except (PermissionError, OSError):
            pass                                    # alive, or not ours to judge


def _validate_query(threshold: float, limit: int) -> tuple[float, int]:
    'Reject the inputs that silently disable the filter or the cap.'
    t = float(threshold)
    if math.isnan(t):
        raise ValueError("threshold is NaN, which disables the similarity floor "
                         "rather than setting it: every comparison passes")
    if not 0.0 <= t <= 1.0:
        raise ValueError(f"threshold must be in [0, 1]; got {t}")
    n = int(limit)
    if n < 1:
        raise ValueError(f"limit must be at least 1; got {n}. A negative limit "
                         "silently drops the highest-scoring rows")
    return t, n


class LotusIndexMissing(RuntimeError):
    """No usable index. Raised rather than returning an empty result set.

    An empty list and "the index was never built" are the same value to a
    caller ranking hits, and this package has been bitten by that confusion
    often enough to name it: a missing lookup must never read as a measured
    absence of natural products.
    """


def _require_index(archive: Optional[Path] = None) -> Path:
    v = verify_index(archive)
    if not v["ok"]:
        raise LotusIndexMissing(
            f"{v['reason']}. Build it with "
            f"`python -c \"from skinpipe.fetch import lotus_local as l; "
            f"print(l.build_index())\"` after downloading the LOTUS export.")
    return index_path(Path(archive) if archive else archive_path())


def similar_to(smiles: str, *, threshold: float = 0.35, limit: int = 200,
               pathway: Optional[str] = None, max_mass: Optional[float] = None,
               archive: Optional[Path] = None) -> list[dict]:
    """Natural products whose ECFP4 Tanimoto to `smiles` is at least `threshold`.

    The default threshold is deliberately low. A scaffold hop *wants* 2D
    distance from its parent -- `hop_score` gates on the Murcko scaffold having
    changed -- so screening natural products at the 0.7 a lead-optimisation
    search would use returns analogues, which is the opposite of the question.
    0.35 keeps enough shared pharmacophore to be worth embedding in 3D while
    leaving room for the scaffold to differ.

    `max_mass` is offered because natural-product space is full of glycosides
    and saponins that no skin barrier admits; filtering here saves embedding
    them. It is a pre-filter for cost, not a substitute for `filters.topical`.

    Returns rows sorted by similarity, each carrying its source organisms and
    NPClassifier pathway. Raises `LotusIndexMissing` if the index is not built.
    """
    from rdkit import Chem, DataStructs, RDLogger
    from rdkit.Chem import rdFingerprintGenerator
    from rdkit.DataStructs.cDataStructs import CreateFromBitString
    RDLogger.DisableLog("rdApp.*")

    threshold, limit = _validate_query(threshold, limit)
    path = _require_index(archive)
    query = Chem.MolFromSmiles(smiles)
    if query is None:
        raise ValueError(f"query SMILES could not be parsed: {smiles!r}")
    if query.GetNumHeavyAtoms() == 0:
        raise ValueError(f"query SMILES {smiles!r} has no heavy atoms, so every "
                         "Tanimoto against it is 0 and the scan ranks nothing")
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=FP_RADIUS, fpSize=FP_BITS)
    qfp = gen.GetFingerprint(query)

    sql = "SELECT inchikey, smiles, name, mass, xlogp, pathway, superclass, npclass, n_organisms, organisms, fp FROM np"
    clauses, params = [], []
    if pathway:
        clauses.append("pathway = ?")
        params.append(pathway)
    if max_mass is not None:
        clauses.append("(mass IS NULL OR mass <= ?)")
        params.append(float(max_mass))
    if clauses:
        sql += " WHERE " + " AND ".join(clauses)

    out: list[dict] = []
    with sqlite3.connect(f"file:{path}?mode=ro", uri=True) as conn:
        for row in conn.execute(sql, params):
            sim = DataStructs.TanimotoSimilarity(
                qfp, CreateFromBitString(row[10].decode("ascii")))
            if sim < threshold:
                continue
            out.append({
                "inchikey": row[0], "smiles": row[1], "name": row[2],
                "mass": row[3], "xlogp": row[4], "pathway": row[5],
                "superclass": row[6], "npclass": row[7],
                "n_organisms": row[8], "organisms": json.loads(row[9] or "[]"),
                "tanimoto_2d": round(sim, 4), "_sim": sim, "source": "lotus",
            })
    # Sorted on the full-precision similarity, not the rounded column.
    # `round(sim, 4)` makes scores that differ in the 5th decimal tie, and a
    # tie straddling `limit` drops whichever the sort happened to order
    # second -- losing a better candidate to a display convention.
    out.sort(key=lambda r: -r["_sim"])
    for r in out:
        r.pop("_sim", None)
    return out[:limit]


def lookup(inchikeys: Iterable[str], archive: Optional[Path] = None) -> dict[str, dict]:
    """Natural-product provenance for structures found elsewhere.

    The annotation path rather than the search path: a compound that arrived
    from ChEMBL or a vendor catalogue gets its organisms and NPClassifier
    pathway attached, which is what lets a report say "this is a plant phenol
    reported in 14 species" about a hit nobody went looking in nature for.
    """
    keys = [k for k in inchikeys if k]
    if not keys:
        return {}
    path = _require_index(archive)
    out: dict[str, dict] = {}
    with sqlite3.connect(f"file:{path}?mode=ro", uri=True) as conn:
        for i in range(0, len(keys), 900):        # SQLite variable limit
            chunk = keys[i:i + 900]
            q = ("SELECT inchikey, smiles, name, pathway, superclass, npclass, "
                 "n_organisms, organisms FROM np WHERE inchikey IN "
                 f"({','.join('?' * len(chunk))})")
            for r in conn.execute(q, chunk):
                out[r[0]] = {"inchikey": r[0], "smiles": r[1], "name": r[2],
                             "pathway": r[3], "superclass": r[4], "npclass": r[5],
                             "n_organisms": r[6],
                             "organisms": json.loads(r[7] or "[]"),
                             "natural_product": True, "np_source": "lotus"}
    return out


def lookup_skeleton(inchikeys: Iterable[str],
                    archive: Optional[Path] = None) -> dict[str, list[dict]]:
    """Same connectivity, any stereochemistry, keyed by the query's full key.

    `lookup` matches the whole InChIKey, which includes the stereochemistry
    layer, and that makes it brittle in the one direction that misleads.
    Staurosporine is in this index; a query built from a SMILES typed without
    its stereocentres hashes to `ZUURDKLOIVWVIB-UHFFFAOYSA-N` against the
    record's `HKSZLNNOFSGOKW-FYTWVXJKSA-N` and comes back "not found" -- a false
    negative that reads exactly like a real absence.

    The first block of an InChIKey is the connectivity skeleton, so matching on
    it finds the compound whatever a caller knew about its stereochemistry. The
    answer is a *list*, because several stereoisomers can share a skeleton and
    they are different compounds with different activities: the caller is told
    "this skeleton occurs in nature, in these forms", which is true, rather than
    being handed one of them as if it were the match.
    """
    keys = [k for k in inchikeys if k]
    if not keys:
        return {}
    path = _require_index(archive)
    out: dict[str, list[dict]] = {}
    with sqlite3.connect(f"file:{path}?mode=ro", uri=True) as conn:
        for key in keys:
            skeleton = str(key).split("-")[0]
            if len(skeleton) != 14:               # not an InChIKey; do not guess
                continue
            rows = conn.execute(
                "SELECT inchikey, smiles, name, pathway, superclass, npclass, "
                "n_organisms, organisms FROM np WHERE inchikey LIKE ?",
                (f"{skeleton}-%",)).fetchall()
            if rows:
                out[key] = [
                    {"inchikey": r[0], "smiles": r[1], "name": r[2],
                     "pathway": r[3], "superclass": r[4], "npclass": r[5],
                     "n_organisms": r[6], "organisms": json.loads(r[7] or "[]"),
                     "natural_product": True, "np_source": "lotus",
                     "exact": r[0] == key}
                    for r in rows]
    return out


def stats(archive: Optional[Path] = None) -> dict:
    """What the built index actually holds, for a report to quote."""
    path = _require_index(archive)
    with sqlite3.connect(f"file:{path}?mode=ro", uri=True) as conn:
        n = conn.execute("SELECT COUNT(*) FROM np").fetchone()[0]
        paths = conn.execute(
            "SELECT pathway, COUNT(*) FROM np WHERE pathway IS NOT NULL "
            "GROUP BY pathway ORDER BY 2 DESC LIMIT 10").fetchall()
        meta = _meta(conn)
    return {"structures": n, "pathways": dict(paths),
            "release": meta.get("release"), "unparseable": meta.get("unparseable")}
