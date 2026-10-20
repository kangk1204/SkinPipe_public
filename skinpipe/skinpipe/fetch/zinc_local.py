'Offline search over a purchasable ZINC subset, held as downloaded tranches.'
from __future__ import annotations

import hashlib
import math
import json
import logging
import os
import sqlite3
from pathlib import Path
from typing import Optional

log = logging.getLogger(__name__)

ENV_VAR = "SKINPIPE_ZINC_DIR"
_DEFAULT = Path.home() / "zinc_local"
_SCHEME = "zinc-tranche-set-blake2b-16"
_FP_KEY = "tranche_fingerprint"

FP_RADIUS = 2
FP_BITS = 2048          # same as lotus_local and catalog, so scores compare

#: ZINC20's fourth letter. Progressive: "in-stock" is also in "wait-ok" files,
#: so an exclusive selection means taking the A files only.
PURCHASABILITY = {"A": "in-stock", "B": "wait-ok", "C": "boutique",
                  "D": "annotated", "E": "agent", "F": "on-demand"}
#: The third letter, strictest first.
REACTIVITY = {"A": "anodyne", "B": "clean", "C": "mild", "E": "reactive",
              "G": "unstable", "I": "all"}


def tranche_dir() -> Path:
    raw = os.environ.get(ENV_VAR)
    return Path(raw).expanduser() if raw else _DEFAULT


def index_path(directory: Optional[Path] = None) -> Path:
    d = Path(directory) if directory else tranche_dir()
    return d / "zinc_tranches.sqlite"


def tranche_files(directory: Optional[Path] = None) -> list[Path]:
    d = Path(directory) if directory else tranche_dir()
    return sorted(d.glob("*.smi")) if d.is_dir() else []


def parse_tranche(name: str) -> dict:
    """What a filename says about everything inside it."""
    stem = Path(name).stem.upper()
    if len(stem) < 4:
        return {"tranche": stem, "mw_bin": None, "logp_bin": None,
                "reactivity": None, "purchasability": None}
    return {"tranche": stem,
            "mw_bin": stem[0], "logp_bin": stem[1],
            "reactivity": REACTIVITY.get(stem[2], stem[2]),
            "purchasability": PURCHASABILITY.get(stem[3], stem[3])}


def fingerprint(directory: Optional[Path] = None) -> Optional[dict]:
    """Identity of the tranche *set*, not of one file.

    Adding or removing a tranche changes which compounds a search can return,
    so it has to invalidate the index just as editing one would.
    """
    files = tranche_files(directory)
    if not files:
        return None
    h = hashlib.blake2b(digest_size=16)
    parts = []
    for f in files:
        st = f.stat()
        parts.append(f"{f.name}:{st.st_size}:{st.st_mtime_ns}")
    h.update("|".join(parts).encode())
    return {"scheme": _SCHEME, "n_files": len(files),
            "total_bytes": sum(f.stat().st_size for f in files),
            "digest": h.hexdigest()}


def _meta(conn: sqlite3.Connection) -> dict:
    try:
        return {k: v for k, v in conn.execute("SELECT key, value FROM meta")}
    except sqlite3.Error:
        return {}


def verify_index(directory: Optional[Path] = None) -> dict:
    d = Path(directory) if directory else tranche_dir()
    out = index_path(d)
    if not out.exists():
        return {"ok": False, "status": "absent", "rows": 0,
                "reason": f"no index at {out}"}
    try:
        with sqlite3.connect(f"file:{out}?mode=ro", uri=True) as conn:
            meta = _meta(conn)
            rows = conn.execute("SELECT COUNT(*) FROM zinc").fetchone()[0]
    except sqlite3.Error as e:
        return {"ok": False, "status": "broken", "rows": 0, "reason": str(e)}
    if not rows:
        return {"ok": False, "status": "empty", "rows": 0, "reason": "index has no rows"}
    raw = meta.get(_FP_KEY)
    if not raw:
        return {"ok": False, "status": "unstamped", "rows": rows,
                "reason": "index carries no tranche fingerprint; run stamp_index "
                          "if it is known good, or build_index(force=True)"}
    if not tranche_files(d):
        return {"ok": False, "status": "no_tranches", "rows": rows,
                "reason": f"no .smi tranches in {d} to verify against"}
    try:
        stored = json.loads(raw)
    except (ValueError, TypeError):
        return {"ok": False, "status": "broken", "rows": rows,
                "reason": "the stored fingerprint is not readable JSON; the "
                          "index cannot be matched to an archive and has to be "
                          "rebuilt"}
    if stored != fingerprint(d):
        return {"ok": False, "status": "stale", "rows": rows,
                "reason": "tranches were added, removed or re-downloaded since "
                          "this index was built"}
    return {"ok": True, "status": "ok", "rows": rows,
            "reason": "index verified against its tranche set"}


def stamp_index(directory: Optional[Path] = None) -> dict:
    d = Path(directory) if directory else tranche_dir()
    out, fp = index_path(d), fingerprint(d)
    if fp is None:
        return {"ok": False, "reason": f"no tranches in {d}"}
    if not out.exists():
        return {"ok": False, "reason": f"no index at {out}"}
    with sqlite3.connect(out) as conn:
        conn.execute("CREATE TABLE IF NOT EXISTS meta (key TEXT PRIMARY KEY, value TEXT)")
        conn.execute("INSERT OR REPLACE INTO meta VALUES (?, ?)",
                     (_FP_KEY, json.dumps(fp, sort_keys=True)))
    return {"ok": True, "reason": "stamped"}


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


class ZincIndexMissing(RuntimeError):
    'No usable index. Raised rather than returning an empty result set, so a'


def _require_index(directory: Optional[Path] = None) -> Path:
    v = verify_index(directory)
    if not v["ok"]:
        raise ZincIndexMissing(
            f"{v['reason']}. Download tranches from files.docking.org/2D/ into "
            f"{tranche_dir()} and run "
            "`python -c \"from skinpipe.fetch import zinc_local as z; "
            "print(z.build_index())\"`.")
    return index_path(Path(directory) if directory else tranche_dir())


def build_index(directory: Optional[Path] = None, *, force: bool = False,
                progress_every: int = 100_000) -> dict:
    """Fingerprint every tranche compound once, into SQLite beside the files.

    Duplicates across tranches are collapsed on ZINC id: a compound can appear
    in more than one file when reactivity levels overlap, and counting it twice
    would inflate every "how much purchasable space did we search" figure in a
    report.
    """
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdFingerprintGenerator
    RDLogger.DisableLog("rdApp.*")

    d = Path(directory) if directory else tranche_dir()
    out = index_path(d)
    if out.exists() and not force:
        v = verify_index(d)
        if v["ok"]:
            return {"ok": True, "built": False, "rows": v["rows"],
                    "path": str(out), "reason": v["reason"]}
        if v["status"] in ("unstamped", "no_tranches"):
            return {"ok": False, "built": False, "rows": v["rows"],
                    "path": str(out), "reason": v["reason"]}
        log.warning("zinc_local: %s; rebuilding", v["reason"])

    files = tranche_files(d)
    if not files:
        return {"ok": False, "built": False, "rows": 0, "path": str(out),
                "reason": f"no .smi tranches in {d}"}

    fp = fingerprint(d)
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=FP_RADIUS, fpSize=FP_BITS)
    # Unique per process: every builder used the same path and
    # unlinked it first, so two concurrent builds of one archive would
    # delete each other's work-in-progress and one could rename a
    # half-written database over the real index.
    _sweep_stale_temps(out)
    tmp = out.with_suffix(f".sqlite.tmp.{os.getpid()}")
    tmp.unlink(missing_ok=True)

    n = skipped = dup = 0
    seen: set[str] = set()
    with sqlite3.connect(tmp) as conn:
        conn.execute("""CREATE TABLE zinc (
            zinc_id TEXT PRIMARY KEY, smiles TEXT NOT NULL, tranche TEXT,
            mw_bin TEXT, logp_bin TEXT, reactivity TEXT, purchasability TEXT,
            fp BLOB NOT NULL)""")
        conn.execute("CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT)")
        for path in files:
            t = parse_tranche(path.name)
            with path.open("r", encoding="utf-8", errors="replace") as fh:
                for line in fh:
                    parts = line.split()
                    if len(parts) < 2 or parts[1].lower() == "zinc_id":
                        continue                       # header or malformed
                    smi, zid = parts[0], parts[1]
                    if zid in seen:
                        dup += 1
                        continue
                    mol = Chem.MolFromSmiles(smi)
                    if mol is None:
                        skipped += 1
                        continue
                    seen.add(zid)
                    conn.execute("INSERT OR REPLACE INTO zinc VALUES (?,?,?,?,?,?,?,?)",
                                 (zid, smi, t["tranche"], t["mw_bin"], t["logp_bin"],
                                  t["reactivity"], t["purchasability"],
                                  gen.GetFingerprint(mol).ToBitString().encode("ascii")))
                    n += 1
                    if progress_every and n % progress_every == 0:
                        log.info("zinc_local: indexed %d compounds", n)
        conn.execute("CREATE INDEX zinc_buy ON zinc(purchasability)")
        conn.execute("INSERT OR REPLACE INTO meta VALUES (?, ?)",
                     (_FP_KEY, json.dumps(fp, sort_keys=True)))
        conn.execute("INSERT OR REPLACE INTO meta VALUES ('tranches', ?)",
                     (json.dumps([f.name for f in files]),))
        conn.execute("INSERT OR REPLACE INTO meta VALUES ('unparseable', ?)", (str(skipped),))
        conn.execute("INSERT OR REPLACE INTO meta VALUES ('duplicates', ?)", (str(dup),))
    tmp.replace(out)
    log.info("zinc_local: indexed %d compounds from %d tranches (%d dup, %d unparseable)",
             n, len(files), dup, skipped)
    return {"ok": True, "built": True, "rows": n, "tranches": len(files),
            "duplicates": dup, "unparseable": skipped, "path": str(out),
            "reason": "built"}


def similar_to(smiles: str, *, threshold: float = 0.35, limit: int = 200,
               purchasability: Optional[str] = "in-stock",
               reactivity: Optional[str] = None,
               directory: Optional[Path] = None) -> list[dict]:
    """Purchasable compounds whose ECFP4 Tanimoto to `smiles` is >= `threshold`.

    `purchasability` defaults to "in-stock" because that is the reason to reach
    for ZINC over a generated library at all: a hop you can order this month is
    worth more than one you can only draw. Pass None to search every tranche
    present.

    The threshold default matches `lotus_local` and for the same reason -- a
    scaffold hop wants 2D distance from its parent, so the 0.7 of a
    lead-optimisation search returns analogues, which is the opposite question.

    Raises `ZincIndexMissing` if no verified index exists.
    """
    from rdkit import Chem, DataStructs, RDLogger
    from rdkit.Chem import rdFingerprintGenerator
    from rdkit.DataStructs.cDataStructs import CreateFromBitString
    RDLogger.DisableLog("rdApp.*")

    threshold, limit = _validate_query(threshold, limit)
    path = _require_index(directory)
    query = Chem.MolFromSmiles(smiles)
    if query is None:
        raise ValueError(f"query SMILES could not be parsed: {smiles!r}")
    if query.GetNumHeavyAtoms() == 0:
        raise ValueError(f"query SMILES {smiles!r} has no heavy atoms, so every "
                         "Tanimoto against it is 0 and the scan ranks nothing")
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=FP_RADIUS, fpSize=FP_BITS)
    qfp = gen.GetFingerprint(query)

    sql = ("SELECT zinc_id, smiles, tranche, mw_bin, logp_bin, reactivity, "
           "purchasability, fp FROM zinc")
    clauses, params = [], []
    if purchasability:
        clauses.append("purchasability = ?")
        params.append(purchasability)
    if reactivity:
        clauses.append("reactivity = ?")
        params.append(reactivity)
    if clauses:
        sql += " WHERE " + " AND ".join(clauses)

    out: list[dict] = []
    with sqlite3.connect(f"file:{path}?mode=ro", uri=True) as conn:
        for r in conn.execute(sql, params):
            sim = DataStructs.TanimotoSimilarity(
                qfp, CreateFromBitString(r[7].decode("ascii")))
            if sim < threshold:
                continue
            out.append({"zinc_id": f"ZINC{int(r[0]):012d}" if r[0].isdigit() else r[0],
                        "smiles": r[1], "tranche": r[2], "mw_bin": r[3],
                        "logp_bin": r[4], "reactivity": r[5],
                        "purchasability": r[6], "tanimoto_2d": round(sim, 4),
                        "_sim": sim, "source": "zinc"})
    # Sorted on the full-precision similarity, not the rounded column.
    # `round(sim, 4)` makes scores that differ in the 5th decimal tie, and a
    # tie straddling `limit` drops whichever the sort happened to order
    # second -- losing a better candidate to a display convention.
    out.sort(key=lambda x: -x["_sim"])
    for x in out:
        x.pop("_sim", None)
    return out[:limit]


def stats(directory: Optional[Path] = None) -> dict:
    """What was actually downloaded, for a report to quote instead of assume."""
    path = _require_index(directory)
    with sqlite3.connect(f"file:{path}?mode=ro", uri=True) as conn:
        n = conn.execute("SELECT COUNT(*) FROM zinc").fetchone()[0]
        buy = dict(conn.execute(
            "SELECT purchasability, COUNT(*) FROM zinc GROUP BY 1 ORDER BY 2 DESC"))
        tr = dict(conn.execute(
            "SELECT tranche, COUNT(*) FROM zinc GROUP BY 1 ORDER BY 2 DESC"))
        meta = _meta(conn)
    return {"compounds": n, "by_purchasability": buy, "by_tranche": tr,
            "tranches": json.loads(meta.get("tranches") or "[]"),
            "duplicates": meta.get("duplicates"),
            "unparseable": meta.get("unparseable")}
