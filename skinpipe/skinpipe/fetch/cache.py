"""Disk cache for the ChEMBL/UniProt/RCSB clients.

Why
---
ChEMBL answers 500 for hours at a time and its FTP mirror has been serving at
27 KB/s, so re-running a screen must not mean re-downloading it. Activity and
molecule pulls are by far the heaviest traffic in a run and previously cached
nothing at all.

The two rules that keep this from becoming its own silent-degradation bug
-----------------------------------------------------------------------
1.  **A partial result is never stored.** The caller decides what "complete"
    means and says so; anything else is fetched again next time. Caching a
    truncated activity pull would make a truncated library permanent, which is
    exactly the defect this pipeline was already carrying.
2.  **Stale data is never served silently.** Past its TTL an entry is only used
    when the network has actually failed, and then the result is flagged
    `stale=True` and logged at WARNING. A cache that quietly answers for a dead
    API turns an outage into a plausible-looking screen, which is the failure
    mode the whole review is about.

Layout: one gzipped JSON file per entry, named by a hash of the full query, so
a crash cannot corrupt more than the entry being written and two runs cannot
race over a shared file. Writes are temp-file-plus-rename.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import math
import logging
import os
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Optional

import requests

log = logging.getLogger(__name__)

#: ChEMBL releases quarterly, so a month-old answer is still the current one.
DEFAULT_TTL_DAYS = 30.0

#: The oldest a stale entry may be and still answer for a failed fetch.
#:
#: Absolute, not a multiple of the TTL, because the question it settles is
#: absolute: ChEMBL releases roughly quarterly, so past two releases the entry
#: is not "possibly out of date", it is a different database. Expressing it as
#: a TTL multiple also made it collapse to zero whenever the TTL was zero --
#: which is how callers ask for "always refetch" -- and that disabled the
#: outage fallback the TTL has nothing to do with.
MAX_STALE_DAYS = 180.0


class FetchError(RuntimeError):
    """A remote service did not answer. The base of `chembl.ChemblError`.

    Lives here, not in `chembl`, because `cached` is what has to recognise it
    and `chembl` imports this module rather than the other way round.
    """


#: What counts as "the network failed" for the stale fallback. Deliberately not
#: `Exception`: a KeyError or TypeError raised inside a fetch function is a
#: defect in this package, and serving last month's answer for it hides the
#: defect behind a plausible screen -- with the TypeError named in a log line
#: nobody reads and `stale=True` explaining it as an outage.
TRANSPORT_FAILURES: tuple[type[BaseException], ...] = (
    FetchError, requests.RequestException, OSError)

_ENV_DIR = "SKINPIPE_CACHE_DIR"
_ENV_TTL = "SKINPIPE_CACHE_TTL_DAYS"
_ENV_OFF = "SKINPIPE_CACHE_DISABLE"


def default_cache_dir() -> Path:
    """`$SKINPIPE_CACHE_DIR`, else a writable per-user cache directory."""
    env = os.environ.get(_ENV_DIR)
    if env:
        return Path(env).expanduser()
    from ..paths import user_cache_dir
    return user_cache_dir() / "requests"


def _default_ttl() -> float:
    try:
        return float(os.environ.get(_ENV_TTL, DEFAULT_TTL_DAYS))
    except ValueError:
        log.warning("%s is not a number; using %s days", _ENV_TTL, DEFAULT_TTL_DAYS)
        return DEFAULT_TTL_DAYS


@dataclass(frozen=True)
class CacheResult:
    """A value plus how it was obtained. `stale` is never silently True."""
    value: Any
    from_cache: bool = False
    age_days: Optional[float] = None
    stale: bool = False

    @property
    def note(self) -> Optional[str]:
        if not self.stale:
            return None
        return (f"served from a {self.age_days:.0f}-day-old cache entry because the "
                f"live fetch failed; the numbers may be out of date")


class DiskCache:
    def __init__(self, root: Optional[Path] = None,
                 ttl_days: Optional[float] = None,
                 enabled: Optional[bool] = None):
        self.root = Path(root) if root is not None else default_cache_dir()
        self.ttl_days = _default_ttl() if ttl_days is None else float(ttl_days)
        if enabled is None:
            enabled = os.environ.get(_ENV_OFF, "").lower() not in ("1", "true", "yes")
        self.enabled = enabled
        self.hits = self.misses = self.writes = 0

    # ------------------------------------------------------------------ keys
    @staticmethod
    def key(kind: str, params: Any) -> str:
        """Hash of the *full* query. pchembl cutoff, max_records and the field
        list all change the answer, so all of them are in the key."""
        blob = json.dumps({"kind": kind, "params": params},
                          sort_keys=True, default=str, separators=(",", ":"))
        return f"{kind}-{hashlib.sha256(blob.encode()).hexdigest()[:24]}"

    def _path(self, key: str) -> Path:
        return self.root / key[:2] / f"{key}.json.gz"

    # --------------------------------------------------------------- read
    def get(self, key: str) -> Optional[CacheResult]:
        if not self.enabled:
            return None
        path = self._path(key)
        try:
            with gzip.open(path, "rt", encoding="utf-8") as fh:
                entry = json.load(fh)
        except FileNotFoundError:
            self.misses += 1
            return None
        except (OSError, EOFError, ValueError, json.JSONDecodeError) as e:
            # A half-written or truncated entry is a miss, never a crash.
            log.warning("discarding unreadable cache entry %s: %s", path.name, e)
            path.unlink(missing_ok=True)
            self.misses += 1
            return None
        # Valid JSON of the wrong shape got past the decode above and broke
        # here: a top-level list raises AttributeError from `.get`, and a
        # `fetched_at` of the wrong type raises inside `float()`. Neither is in
        # the list caught above, so both surfaced as a crash in whatever was
        # doing the lookup rather than as the miss they are.
        stamp = entry.get("fetched_at", 0) if isinstance(entry, dict) else None
        try:
            fetched_at = float(stamp)                      # type: ignore[arg-type]
            if not math.isfinite(fetched_at):
                raise ValueError(f"fetched_at is {stamp!r}")
        except (TypeError, ValueError) as e:
            log.warning("discarding cache entry %s with an unusable shape: %s",
                        path.name, e)
            path.unlink(missing_ok=True)
            self.misses += 1
            return None
        age = max(0.0, (time.time() - fetched_at) / 86400.0)
        fresh = age <= self.ttl_days
        self.hits += 1 if fresh else 0
        self.misses += 0 if fresh else 1
        return CacheResult(entry.get("payload"), from_cache=True,
                           age_days=age, stale=not fresh)

    # -------------------------------------------------------------- write
    def put(self, key: str, payload: Any, kind: str = "") -> None:
        if not self.enabled:
            return
        path = self._path(key)
        path.parent.mkdir(parents=True, exist_ok=True)
        entry = {"kind": kind, "fetched_at": time.time(), "payload": payload}
        fd, tmp = tempfile.mkstemp(dir=str(path.parent), suffix=".part")
        os.close(fd)
        try:
            with gzip.open(tmp, "wt", encoding="utf-8") as fh:
                json.dump(entry, fh, default=str)
            os.replace(tmp, path)
            self.writes += 1
        except BaseException:
            Path(tmp).unlink(missing_ok=True)
            raise

    def stats(self) -> dict:
        return {"hits": self.hits, "misses": self.misses, "writes": self.writes,
                "root": str(self.root), "ttl_days": self.ttl_days,
                "enabled": self.enabled}


def cached(cache: Optional[DiskCache], kind: str, params: Any,
           fetch: Callable[[], Any],
           cacheable: Callable[[Any], bool] = lambda v: True) -> CacheResult:
    """Fresh cache hit, else fetch; fall back to a stale entry only on failure.

    `cacheable` is the caller's judgement on whether the value it just fetched
    is a complete answer. A False here means the value is returned but not
    stored, so the next run tries again instead of inheriting the gap.

    The stale fallback is for an *outage*, and the two limits below say so.
    `except Exception` caught programming errors too, so a `TypeError` in the
    fetch function was reported as "ChEMBL did not answer" and served from
    cache -- a bug disguised as weather, and the log line named the TypeError
    where nobody reads it. And there was no ceiling: a 400-day-old entry was
    served against a 30-day TTL, which is not "possibly out of date", it is a
    different ChEMBL release. Past `MAX_STALE_DAYS` the original
    failure is raised instead.
    """
    if cache is None or not cache.enabled:
        return CacheResult(fetch())
    key = cache.key(kind, params)
    hit = cache.get(key)
    if hit is not None and not hit.stale:
        log.debug("cache hit %s (%.1f days old)", key, hit.age_days or 0.0)
        return hit
    try:
        value = fetch()
    except TRANSPORT_FAILURES as e:
        if hit is None:
            raise
        if (hit.age_days or 0.0) > MAX_STALE_DAYS:
            log.error("%s fetch failed (%s) and the only cache entry is %.0f days "
                      "old, past the %.0f-day serving limit; refusing to answer "
                      "from it", kind, e, hit.age_days or 0.0, MAX_STALE_DAYS)
            raise
        log.warning("%s fetch failed (%s); falling back to a %.0f-day-old cache "
                    "entry. Results may be out of date.", kind, e, hit.age_days or 0.0)
        return CacheResult(hit.value, from_cache=True,
                           age_days=hit.age_days, stale=True)
    if cacheable(value):




        try:
            cache.put(key, value, kind=kind)
        except OSError as exc:
            log.warning("%s: fetched successfully but the result could not be "
                        "cached (%s); continuing without caching it", kind, exc)
    else:
        log.info("%s result is incomplete; not cached", kind)
    return CacheResult(value)
