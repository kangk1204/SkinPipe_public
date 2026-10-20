"""Can we buy it, and is it already someone else's?

Two questions decide whether a computational hit is worth anything, and neither
is about binding.

**Can it be bought.** A hit that exists only in a paper cannot be assayed this
month. ChEMBL is a bioactivity database, not a catalogue: a large share of its
potent compounds were made once, in one laboratory, for one publication. Ranking
them without asking who sells them produces a list that stalls the week it
lands on a chemist's desk.

**Is it already claimed.** For a company pursuing its own ingredient, a compound
that appears in the patent literature is a different proposition from one that
does not. This module reports presence in patents. It does **not** assess freedom
to operate: that needs a real patent search and legal review, and no structure
lookup can substitute for either.

Both answers come from UniChem (EMBL-EBI), which cross-references a structure
against 25 databases from one InChIKey. Two of those sources carry the answers:
MolPort for commercial availability, SureChEMBL for patent-literature presence.
It is free, needs no key, and covers the whole library in one pass.

Three more sources are worth reading for a cosmetic programme specifically,
because a 기능성화장품 dossier lives or dies on the safety package:

* **FooDB** -- the compound occurs in food, so there is human oral exposure
  history. That is the single most helpful thing a novel cosmetic active can
  have.
* **EPA CompTox** -- regulatory toxicology data already exists.
* **DrugCentral / ClinicalTrials** -- the compound is an approved or clinically
  studied drug. That does not bar it from cosmetics: niacinamide, adenosine and
  azelaic acid are all simultaneously drugs and cosmetic ingredients in Korea.
  What is regulated is the *claim*, not the molecule. The signal matters because
  a drug has a clinical safety record to read across from, and because the claim
  wording will need care.
"""
from __future__ import annotations

import contextlib
try:                                    # POSIX only
    import fcntl
except ImportError:                     # pragma: no cover
    fcntl = None                        # `_cache_file_lock` degrades
import json
import hashlib
import os
import logging
import re
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Callable, Iterable, Optional

import requests
from rdkit import Chem, RDLogger
from rdkit.Chem import inchi

RDLogger.DisableLog("rdApp.*")
log = logging.getLogger(__name__)

UNICHEM = "https://www.ebi.ac.uk/unichem/api/v1/compounds"
PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
_UA = {"User-Agent": "skinpipe/1.0 (cosmetic ingredient discovery)"}
_TIMEOUT = 60

#: UniChem short names we act on.
VENDOR_SOURCES = ("molport",)








VENDOR_NAME_PATTERNS: tuple[str, ...] = (
    "molport", "mcule", "sigma", "aldrich", "enamine", "chemspace", "ambeed",
    "key organics", "bionet", "apollo scientific", "tci ", "tokyo chemical",
    "alfa aesar", "thermo", "biosynth", "targetmol", "medchemexpress",
    "selleck", "cayman", "santa cruz", "toronto research", "combi-blocks",
    "matrix scientific", "fluorochem", "acros", "abcam", "vitas-m", "life chemicals",
    "chembridge", "asinex", "specs", "otava", "princeton biomolecular",
)
PATENT_SOURCES = ("surechembl",)
FOOD_SOURCES = ("foodb",)
TOX_SOURCES = ("comptox",)
#: Approved or clinically studied. Deliberately excludes FDA SRS, which registers
#: every regulated substance including food additives and cosmetic ingredients and
#: so says nothing about drug status, and DrugBank, which mixes approved drugs with
#: experimental and nutraceutical entries.
DRUG_SOURCES = ("drugcentral", "clinicaltrials")


def classification_policy() -> str:
    """Bind interpreted cache rows to the rules that gave them meaning."""
    return hashlib.sha256(json.dumps({
        "version": "sourcing-classification-v1", "vendors": VENDOR_SOURCES,
        "vendor_names": VENDOR_NAME_PATTERNS, "patents": PATENT_SOURCES,
        "food": FOOD_SOURCES, "toxicology": TOX_SOURCES, "drugs": DRUG_SOURCES,
    }, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


@dataclass
class SourcingRecord:
    inchikey: str
    classification_policy_sha256: str = field(default_factory=classification_policy)
    n_sources: int = 0
    purchasable: bool = False
    vendor_ids: dict[str, str] = field(default_factory=dict)
    n_vendors: int = 0                       # distinct suppliers naming this structure
    vendors: list[str] = field(default_factory=list)
    pubchem_cid: Optional[int] = None
    # None, not False: "no SureChEMBL record" and "the lookup never answered"
    # are different facts, and the buy list is where that difference costs a
    # freedom-to-operate opinion. The `purchasable` side of this same record
    # already refuses to cache a negative born of a 503 (see `lookup` below);
    # this field had no such guard, so a UniChem outage returned the dataclass
    # default and the Markdown printed it as though a lookup had cleared it.
    in_patents: Optional[bool] = None
    patent_ids: list[str] = field(default_factory=list)
    food_history: bool = False
    tox_data: bool = False
    is_drug: bool = False
    sources: list[str] = field(default_factory=list)
    #: Registry numbers PubChem lists as synonyms for this structure, check-digit
    #: verified. A list because a structure routinely has several and only one
    #: is the substance meant; `eligibility` checks them all rather than
    #: guessing. `cas_checked` separates "PubChem answered and has none" from
    #: "the lookup never got an answer", the same way `status` does for vendors.
    cas: list[str] = field(default_factory=list)
    cas_checked: bool = False
    status: str = "unknown"          # purchasable | known_not_sold | novel | unknown
    ip_note: str = ""
    dossier_note: str = ""
    cache_provenance: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        return asdict(self)


#: HTTP statuses that mean "ask again later", not "the answer is no". Treating a
#: 503 as an answer is how a throttled PubChem turned into a screen-wide report
#: that nothing was for sale; see ``_pubchem_vendors``.
RETRYABLE_STATUS = frozenset({408, 429, 500, 502, 503, 504})

#: PubChem PUG-REST asks for no more than five requests a second and enforces it
#: with HTTP 503. One compound costs up to three requests (UniChem, PubChem CID,
#: PubChem xrefs), so the budget is shared across every worker rather than
#: applied per thread -- N threads must not mean N times the rate.
DEFAULT_RATE = 5.0

#: Threads in the sourcing pool. Modest on purpose: the work is 100% network
#: bound, and the point of the pool is to overlap latency, not to flood a public
#: API that asks us not to.
DEFAULT_WORKERS = 4

#: Write the cache to disk every this many completed compounds. A 1493-compound
#: marker takes tens of minutes even parallelised; losing all of it to one
#: crash at minute forty is the failure this exists to prevent.
DEFAULT_SAVE_EVERY = 100


def _backoff(attempt: int) -> None:
    """Sleep between retries. Factored out so tests can run without waiting."""
    time.sleep(2 ** attempt)


class RateLimiter:
    """At most ``rate`` acquisitions per second, shared across threads.

    Each caller reserves the next free slot under the lock and then sleeps until
    that slot arrives, so the slots are handed out in request order and no
    caller can spin. ``clock`` and ``sleep`` are injectable because a test that
    asserts on a rate cannot afford to spend real seconds proving it.
    """

    def __init__(self, rate: float = DEFAULT_RATE,
                 clock: Callable[[], float] = time.monotonic,
                 sleep: Callable[[float], None] = time.sleep):
        self.rate = float(rate)
        self.min_interval = 1.0 / self.rate if self.rate > 0 else 0.0
        self._clock = clock
        self._sleep = sleep
        self._lock = threading.Lock()
        self._next_slot: Optional[float] = None

    def acquire(self) -> None:
        if self.min_interval <= 0:
            return
        with self._lock:
            now = self._clock()
            slot = now if self._next_slot is None else max(now, self._next_slot)
            self._next_slot = slot + self.min_interval
        delay = slot - self._clock()
        if delay > 0:
            self._sleep(delay)


def _acquire(limiter: Optional[RateLimiter]) -> None:
    if limiter is not None:
        limiter.acquire()


@contextlib.contextmanager
def _cache_file_lock(path: Path):
    """Exclusive lock for a cache file, held on a sidecar next to it.

    Not on the file itself: `save` replaces it by rename, and a lock held on
    the replaced inode protects nothing afterwards. Never raises -- a cache is
    an optimisation, and a platform or filesystem without `flock` should get
    the old unsynchronised behaviour rather than an exception in a `finally`.
    """
    lock_path = path.with_name(path.name + ".lock")
    fh = None
    try:
        lock_path.parent.mkdir(parents=True, exist_ok=True)
        if fcntl is None:
            raise OSError("fcntl is unavailable on this platform")
        fh = open(lock_path, "a+")
        fcntl.flock(fh.fileno(), fcntl.LOCK_EX)
    except (OSError, AttributeError) as e:
        log.debug("sourcing cache lock unavailable (%s); saving unsynchronised", e)
        fh = None
    try:
        yield
    finally:
        if fh is not None:
            try:
                fcntl.flock(fh.fileno(), fcntl.LOCK_UN)
            except OSError:
                pass
            fh.close()


_LEGACY_SOURCING_CACHE = (Path(__file__).resolve().parent.parent
                          / "config" / ".sourcing_cache.json")


class SourcingCache:
    """Disk cache keyed on InChIKey. UniChem is stable; re-querying is waste.

    Every accessor takes a lock: the lookup pool writes from several threads at
    once, and ``dict`` mutation plus a concurrent ``json.dumps`` of the same
    dict is how entries quietly go missing (or ``save`` dies on "dictionary
    changed size during iteration").
    """

    def __init__(self, path: Optional[Path] = None):
        from .paths import user_cache_dir

        override = os.environ.get("SKINPIPE_SOURCING_CACHE")
        self.path = Path(path or override or user_cache_dir() / "sourcing.json").expanduser()
        read_path = self.path
        # An installed package may be read only. Import its old cache only on
        # the default path's first use; an explicit override starts isolated.
        if path is None and not override and not self.path.exists():
            read_path = _LEGACY_SOURCING_CACHE
        self._legacy_data: dict[str, dict] = {}
        self._data: dict[str, dict] = {}
        #: Keys this process wrote. `save` merges only these over the file, so
        #: an entry another process refreshed is not reverted by our copy of
        #: what the file said when we loaded it. Before the CAS top-up below
        #: nothing ever rewrote an existing key, so "in-memory wins" could not
        #: lose anything; a top-up makes the same key writable from two runs
        #: at once and that assumption stops holding.
        self._dirty: set[str] = set()
        self._lock = threading.RLock()
        if read_path.is_file():
            try:
                loaded = json.loads(read_path.read_text())
            except (OSError, json.JSONDecodeError) as exc:
                log.warning("sourcing cache at %s could not be read (%s); "
                            "starting empty", read_path, exc)
            else:
                # Valid JSON of the wrong shape passed the decode and then
                # broke on first use: a top-level list made `self._data.get`
                # raise AttributeError from inside `get`, which reads as a bug
                # in the lookup rather than a bad file. The shape is part of
                # the format, so it is checked where the format is.
                if isinstance(loaded, dict):
                    self._data = loaded
                    if read_path != self.path:
                        self._legacy_data = dict(loaded)
                else:
                    log.warning("sourcing cache at %s holds %s, not an object; "
                                "starting empty", read_path,
                                type(loaded).__name__)

    def get(self, key: str) -> Optional[dict]:
        with self._lock:
            return self._data.get(key)

    def put(self, key: str, value: dict) -> None:
        with self._lock:
            self._data[key] = value
            self._dirty.add(key)

    def __len__(self) -> int:
        with self._lock:
            return len(self._data)

    def save(self) -> None:
        """Persist, merging whatever another process wrote since we loaded.

        The lock makes this safe between threads. It does nothing between
        *processes*, and two screens running at once is the normal case here:
        a coverage sweep over every marker alongside a docking run. Each holds
        its own copy of the whole cache, so a plain write means the second one
        to finish silently discards every lookup the first one made -- hours of
        network time, lost with no error. Re-reading and merging first costs one
        file read and makes the two runs cooperative instead of competing.

        Entries already in memory win over the ones on disk: they are from this
        process and are at least as fresh.
        """
        self.path.parent.mkdir(parents=True, exist_ok=True)
        # The merge below is a read-modify-write on one file. The instance lock
        # does nothing between processes, so two overlapping saves each read a
        # pre-merge file and the later rename wins whole: a reviewer measured
        # 50 of 100 entries lost that way. An exclusive flock on a sidecar --
        # not on the target, which is replaced by rename and would drop the
        # lock with the old inode -- closes the window. A platform without
        # flock degrades to the previous behaviour rather than failing to save.
        with _cache_file_lock(self.path):
            self._save_locked()

    def _save_locked(self) -> None:
        on_disk: dict[str, dict] = {}
        # Whether the file is still a readable record of what everyone else
        # has learned. When it is, our copy of a key we never wrote is stale
        # by definition and the file's version wins. When it is not, our copy
        # is the only surviving record and all of it goes back.
        disk_readable = False
        if self.path.is_file():
            try:
                loaded = json.loads(self.path.read_text())
                if isinstance(loaded, dict):
                    on_disk = loaded
                    disk_readable = True
            except (OSError, json.JSONDecodeError):
                # A corrupt or half-written file is a lost cache, not a reason
                # to lose what this process learned.
                log.warning("sourcing cache at %s unreadable while saving; "
                            "writing this process's entries only", self.path)
        with self._lock:
            # Only the keys this process wrote override the file. Everything
            # else we hold is just what the file said when we loaded it, and
            # the file may have moved on since -- another run's CAS top-up for
            # a key we merely read would be written back to its old value.
            written = set(self._dirty)
            mine = ({k: v for k, v in self._data.items() if k in written}
                    if disk_readable else dict(self._data))
            # Another process may have created the new user cache after we
            # read the package. Import absent old keys, never overwrite that
            # process's newer values with inherited package data.
            merged = {**self._legacy_data, **on_disk, **mine}
            self._data = merged
            if disk_readable and not written and merged == on_disk:
                return  # read-only compatibility must not rewrite original bytes
            payload = json.dumps(merged, indent=0, ensure_ascii=False)
        # Write-then-rename: a run killed mid-write leaves the previous cache
        # intact instead of a truncated file that the next run discards whole.
        # The temp name carries the pid as well as the thread id, because two
        # processes can hand out the same thread id.
        tmp = self.path.with_name(
            f"{self.path.name}.tmp.{os.getpid()}.{threading.get_ident()}")
        tmp.write_text(payload)
        tmp.replace(self.path)
        # Only now. Clearing before the durable replace meant a save that hit
        # a full disk left the entry in memory, absent from the file, and no
        # longer dirty -- so the obvious recovery, calling `save()` again,
        # wrote the file back without it and reported success. Measured by
        # injecting one OSError into `Path.replace`: memory new, disk old,
        # dirty empty, and the retry left disk old.
        #
        # `written` is the set as it stood when the payload was built, not
        # `self._dirty` as it stands now: another thread may have put a key
        # after that snapshot, and that key is not in the file yet.
        with self._lock:
            self._dirty -= written


def inchikey_for(smiles: str) -> Optional[str]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        return inchi.MolToInchiKey(mol)
    except Exception:                                    # pragma: no cover
        return None


#: A CAS Registry Number as it is written: 2-7 digits, 2 digits, 1 check digit.
#: Shape only -- `_valid_cas` is what decides whether it is a real one.
_CAS_SHAPE = re.compile(r"^(\d{2,7})-(\d{2})-(\d)$")


def _valid_cas(value: str) -> bool:
    """Does `value` pass the CAS check-digit test?

    The last digit is a checksum: read every other digit right to left,
    multiply by 1, 2, 3 ... and take the sum mod 10. Shape alone is not enough
    here. These numbers arrive as PubChem depositor synonyms, where a typo, a
    truncated registry string or another identifier that happens to look like a
    CAS all pass a regex -- and `eligibility` looks them up against a
    prohibition list, so a number that is merely well-formed matches a
    *different regulated substance*. `pipeline._CAS_PROP_KEYS` says the same
    thing about why its list is narrow: a wrong CAS is worse than none.
    """
    m = _CAS_SHAPE.match(value or "")
    if m is None:
        return False
    body, check = m.group(1) + m.group(2), int(m.group(3))
    return sum(i * int(d) for i, d in enumerate(reversed(body), start=1)) % 10 == check


def pubchem_cas(smiles: str, session: Optional[requests.Session] = None,
                limiter: Optional[RateLimiter] = None, retries: int = 3
                ) -> tuple[list[str], bool]:
    '``(cas_numbers, answered)`` for one structure, from PubChem synonyms.'
    session = session or requests.Session()
    for attempt in range(retries):
        _acquire(limiter)
        try:
            # POST, not GET: a SMILES carries "/", "=" and "#", and a URL-path
            # lookup of `Oc1ccc(/C=C/c2cc(O)cc(O)c2)cc1` comes back 400.
            r = session.post(f"{PUBCHEM}/compound/smiles/synonyms/JSON",
                             data={"smiles": smiles}, headers=_UA, timeout=_TIMEOUT)
        except requests.RequestException as exc:
            log.debug("PubChem synonyms attempt %d for %r failed: %s",
                      attempt + 1, smiles[:40], exc)
            _backoff(attempt)
            continue
        if r.status_code == 404:
            return [], True              # PubChem answered: it has no record
        if r.status_code in RETRYABLE_STATUS:
            log.debug("PubChem synonyms %r -> HTTP %s, retrying",
                      smiles[:40], r.status_code)
            _backoff(attempt)
            continue
        if r.status_code != 200:
            log.warning("PubChem synonyms %r -> HTTP %s", smiles[:40], r.status_code)
            return [], False
        try:
            payload = r.json()
        except ValueError:
            _backoff(attempt)
            continue
        # PubChem answered 200, so retrying will not change the body. Every
        # level is checked because each has been the wrong type in practice:
        # `InformationList` as a list makes `["Information"]` raise TypeError,
        # and a string `Information[0]` makes `.get` raise AttributeError.
        # Neither was caught, so one malformed compound ended the batch.
        node = payload.get("InformationList") if isinstance(payload, dict) else None
        info = node.get("Information") if isinstance(node, dict) else None
        if not isinstance(info, list):
            log.warning("PubChem synonyms %r: a 200 response with no usable "
                        "InformationList.Information; treating this compound as "
                        "unresolved", smiles[:40])
            return [], False
        first = info[0] if info else None
        synonyms = first.get("Synonym", []) if isinstance(first, dict) else []
        if not isinstance(synonyms, (list, tuple)):
            log.warning("PubChem synonyms %r: Synonym is %s, not a list; "
                        "treating this compound as unresolved",
                        smiles[:40], type(synonyms).__name__)
            return [], False
        seen: list[str] = []
        for s in synonyms:
            s = str(s).strip()
            if _valid_cas(s) and s not in seen:
                seen.append(s)
        return seen, True
    return [], False


def _pubchem_vendors(smiles: str, session: Optional[requests.Session] = None,
                     limiter: Optional[RateLimiter] = None, retries: int = 3
                     ) -> tuple[Optional[int], list[str], bool]:
    """``(cid, vendor_names, answered)``.

    ``answered`` is the part that matters. PubChem throttles with HTTP 503 and
    falls over with 500, and the old code read either as an empty vendor list --
    so a compound every supplier on earth stocks was filed as "known_not_sold"
    and cached, permanently, because the server was busy for one second. A
    retryable status is retried with backoff; if it never answers, the caller is
    told so rather than being handed a negative it cannot distinguish from a
    real one.
    """
    session = session or requests.Session()
    cid: Optional[int] = None

    for attempt in range(retries):
        _acquire(limiter)
        try:
            r = session.post(f"{PUBCHEM}/compound/smiles/cids/JSON",
                             data={"smiles": smiles}, headers=_UA, timeout=_TIMEOUT)
        except requests.RequestException as exc:
            log.debug("PubChem CID attempt %d for %r failed: %s",
                      attempt + 1, smiles[:40], exc)
            _backoff(attempt)
            continue
        if r.status_code == 404:
            return None, [], True                # PubChem genuinely has no CID
        if r.status_code in RETRYABLE_STATUS:
            log.debug("PubChem CID %r -> HTTP %s, retrying", smiles[:40], r.status_code)
            _backoff(attempt)
            continue
        if r.status_code != 200:
            log.warning("PubChem CID %r -> HTTP %s", smiles[:40], r.status_code)
            return None, [], False
        try:
            payload = r.json()
        except ValueError:
            _backoff(attempt)
            continue
        # PubChem answered 200, so the transport is fine; the body is not what
        # the contract says. `payload["IdentifierList"]` has been a list, and
        # `CID` a list of strings -- `.get` on the first raises AttributeError
        # and `int()` on the second ValueError, either of which ends the whole
        # batch from inside one compound's lookup. An unusable answer is this
        # compound's lookup failing, which the caller already handles.
        ident = payload.get("IdentifierList") if isinstance(payload, dict) else None
        cids = ident.get("CID") if isinstance(ident, dict) else None
        if not isinstance(cids, (list, tuple)) or not cids or not cids[0]:
            if cids is None:
                log.warning("PubChem CID %r: no IdentifierList.CID in a 200 "
                            "response; treating this compound as unresolved",
                            smiles[:40])
                return None, [], False
            return None, [], True                # answered, and the answer is none
        try:
            cid = int(cids[0])
        except (TypeError, ValueError):
            log.warning("PubChem CID %r: %r is not a CID; treating this "
                        "compound as unresolved", smiles[:40], cids[0])
            return None, [], False
        break
    if cid is None:
        return None, [], False

    for attempt in range(retries):
        _acquire(limiter)
        try:
            r2 = session.get(f"{PUBCHEM}/compound/cid/{cid}/xrefs/SourceName/JSON",
                             headers=_UA, timeout=_TIMEOUT)
        except requests.RequestException as exc:
            log.debug("PubChem xrefs attempt %d for CID %s failed: %s",
                      attempt + 1, cid, exc)
            _backoff(attempt)
            continue
        if r2.status_code == 404:
            return cid, [], True                 # no depositor cross-references
        if r2.status_code in RETRYABLE_STATUS:
            log.debug("PubChem xrefs CID %s -> HTTP %s, retrying", cid, r2.status_code)
            _backoff(attempt)
            continue
        if r2.status_code != 200:
            log.warning("PubChem xrefs CID %s -> HTTP %s", cid, r2.status_code)
            return cid, [], False
        try:
            payload = r2.json()
        except ValueError:
            _backoff(attempt)
            continue
        node = payload.get("InformationList") if isinstance(payload, dict) else None
        info = node.get("Information") if isinstance(node, dict) else None
        first = info[0] if isinstance(info, list) and info else None
        names = first.get("SourceName") if isinstance(first, dict) else None
        if (not isinstance(names, list)
                or any(not isinstance(name, str) for name in names)):
            log.warning("PubChem xrefs CID %s: unusable SourceName in a 200 "
                        "response; treating this compound as unresolved", cid)
            return cid, [], False
        seen: list[str] = []
        for n in names:
            low = n.lower()
            if any(p in low for p in VENDOR_NAME_PATTERNS) and n not in seen:
                seen.append(n)
        return cid, seen, True
    return cid, [], False


def pubchem_vendors(smiles: str, session: Optional[requests.Session] = None,
                    limiter: Optional[RateLimiter] = None
                    ) -> tuple[Optional[int], list[str]]:
    """Distinct chemical suppliers that have deposited this structure.

    Returns ``(cid, vendor_names)``. A supplier deposit is evidence the
    compound has been sold, not a live quote: catalogues go stale and a listing
    is not a price. It is still far better evidence than one aggregator's
    snapshot, and it is the difference between "nobody sells this" and "we did
    not look in the right place".

    Callers that need to tell "no vendors" from "PubChem did not answer" want
    :func:`_pubchem_vendors`, which returns that third value.
    """
    cid, vendors, _answered = _pubchem_vendors(smiles, session, limiter=limiter)
    return cid, vendors


def _query(inchikey: str, session: requests.Session, retries: int = 3,
           limiter: Optional[RateLimiter] = None) -> Optional[dict]:
    for attempt in range(retries):
        _acquire(limiter)
        try:
            r = session.post(UNICHEM, json={"type": "inchikey", "compound": inchikey},
                             headers=_UA, timeout=_TIMEOUT)
            if r.status_code == 200:
                return r.json()
            if r.status_code == 404:
                return {"compounds": []}                 # genuinely unknown structure
            if r.status_code in RETRYABLE_STATUS:
                _backoff(attempt)
                continue
            log.warning("UniChem %s -> HTTP %s", inchikey, r.status_code)
            return None
        except requests.RequestException as exc:
            log.debug("UniChem attempt %d for %s failed: %s", attempt + 1, inchikey, exc)
            _backoff(attempt)
    return None


#: Appended when the PubChem supplier check never answered. The record is not
#: cached in that case, so the next run asks again.
_PUBCHEM_SILENT = ("The PubChem supplier check did not answer, so commercial "
                   "availability was not established for this compound and nothing "
                   "here should be read as evidence that it cannot be bought.")


def _apply_notes(rec: SourcingRecord) -> None:
    'Derive ``ip_note`` and ``dossier_note`` from the record as it now stands.'
    #: UniChem resolved nothing at all for this structure.
    unichem_blind = rec.n_sources == 0
    #: ...but PubChem named suppliers, which means it is not novel, it is
    #: outside UniChem's cross-reference coverage.
    sold_per_pubchem = bool(rec.vendors)

    if rec.in_patents:
        rec.ip_note = ("appears in the patent literature (SureChEMBL). Presence is not "
                       "the same as an in-force claim covering your use; a freedom-to-"
                       "operate search is required before committing.")
    elif unichem_blind and sold_per_pubchem:
        rec.ip_note = (f"not found in any of UniChem's 25 databases, yet PubChem names "
                       f"{rec.n_vendors} supplier(s) for it. That makes the UniChem miss a "
                       "coverage gap rather than novelty: the structure is known and on "
                       "sale. No patent check was possible either, because SureChEMBL is "
                       "reached through UniChem and UniChem did not resolve it, so nothing "
                       "here speaks to freedom to operate.")
    elif unichem_blind:
        rec.ip_note = ("not found in any of UniChem's 25 databases, so it is structurally "
                       "novel and unclaimed as far as this check can see. Absence of "
                       "evidence is not freedom to operate.")
    else:
        rec.ip_note = "known compound with no patent-literature record found."

    notes = []
    if rec.food_history:
        notes.append("occurs in food, so human oral exposure history exists")
    if rec.tox_data:
        notes.append("regulatory toxicology data available (EPA CompTox)")
    if rec.is_drug:
        notes.append("approved or clinically studied drug, so a clinical safety record "
                     "exists to read across from. It does not bar cosmetic use -- what is "
                     "regulated is the claim, not the molecule -- but the claim wording "
                     "must stay clear of therapeutic language")
    if unichem_blind and not sold_per_pubchem:
        notes.append("a novel structure carries the full safety dossier burden: no read-"
                     "across, no exposure history, and a longer 기능성화장품 review")
    elif unichem_blind:
        notes.append("suppliers list it, but UniChem holds no record of the structure, so "
                     "there is no food, toxicology or clinical data to read across from "
                     "and the 기능성화장품 dossier burden is close to that of a novel "
                     "structure")
    rec.dossier_note = "; ".join(notes)

    if rec.status == "unknown":
        rec.ip_note = f"{rec.ip_note} {_PUBCHEM_SILENT}".strip()


def _interpret(inchikey: str, payload: dict) -> SourcingRecord:
    found: dict[str, list[str]] = {}
    for compound in payload.get("compounds", []):
        for src in compound.get("sources", []):
            name = (src.get("shortName") or str(src.get("id") or "")).lower()
            if name:
                found.setdefault(name, []).append(str(src.get("compoundId") or ""))

    rec = SourcingRecord(inchikey=inchikey, n_sources=len(found),
                         sources=sorted(found))
    rec.vendor_ids = {s: found[s][0] for s in VENDOR_SOURCES if s in found}
    rec.purchasable = bool(rec.vendor_ids)
    rec.patent_ids = [i for s in PATENT_SOURCES for i in found.get(s, [])][:5]
    rec.in_patents = bool(rec.patent_ids)
    rec.food_history = any(s in found for s in FOOD_SOURCES)
    rec.tox_data = any(s in found for s in TOX_SOURCES)
    rec.is_drug = any(s in found for s in DRUG_SOURCES)

    if rec.purchasable:
        rec.status = "purchasable"
    elif rec.n_sources:
        rec.status = "known_not_sold"
    else:
        rec.status = "novel"

    _apply_notes(rec)
    return rec


def _cas_unresolved(hit: dict) -> bool:
    """Was this entry's CAS never obtained, as opposed to obtained and empty?

    `cas_checked` distinguishes the two, and an entry predating the field
    reports the dataclass default, `False` -- the same value as a lookup that
    failed. Both mean the registry numbers were never established, which is a
    fact about a past run and not about the compound, so both are worth one
    retry. `cas_checked=True` with `cas=[]` is a real answer and is left alone.
    """
    return not hit.get("cas_checked")


def _cache_row_schema(hit: object, key: str) -> Optional[dict]:
    """A cached entry that `SourcingRecord(**hit)` will actually accept.

    The file's top-level shape is checked when it loads; this is the row. A
    record carrying a field the dataclass does not have, or missing one it
    requires, raised `TypeError` from the constructor -- which reads as a bug
    in the lookup rather than as the stale cache entry it is. Unknown fields
    are not silently dropped: an entry written by a newer version is not one
    this version can interpret, so it is refetched.
    """
    if hit is None:
        return None
    if not isinstance(hit, dict):
        log.warning("sourcing cache entry %s is %s, not an object; refetching",
                    key, type(hit).__name__)
        return None
    import dataclasses
    fields = {f.name for f in dataclasses.fields(SourcingRecord)}
    required = {f.name for f in dataclasses.fields(SourcingRecord)
                if f.default is dataclasses.MISSING
                and f.default_factory is dataclasses.MISSING}   # type: ignore[misc]
    unknown = set(hit) - fields
    missing = required - set(hit)
    if unknown or missing:
        log.warning("sourcing cache entry %s does not match this version's "
                    "record (unexpected %s, missing %s); refetching",
                    key, sorted(unknown) or "none", sorted(missing) or "none")
        return None
    if hit.get("inchikey") != key:
        log.warning("sourcing cache entry %s has a different structure identity; refetching", key)
        return None
    return hit


def _usable_hit(hit: object, key: str) -> Optional[dict]:
    hit = _cache_row_schema(hit, key)
    if hit is None:
        return None
    if hit.get("classification_policy_sha256") != classification_policy():
        log.info("sourcing cache entry %s predates the current classification policy; refetching", key)
        return None
    return hit


def _legacy_source_evidence(hit: object, key: str) -> Optional[SourcingRecord]:
    """Reinterpret retained positive evidence, without sealing a stale verdict.

    Old records have no policy identity. Their booleans/status cannot establish
    current negatives. Named source/vendor evidence can still be interpreted by
    the current rules if a refresh fails. The original row remains untouched;
    a content hash and explicit migration version travel with the fallback.
    """
    hit = _cache_row_schema(hit, key)
    if hit is None or hit.get("classification_policy_sha256"):
        return None
    sources, vendors, vendor_ids = hit.get("sources", []), hit.get("vendors", []), hit.get("vendor_ids", {})
    if not (isinstance(sources, list) and all(isinstance(s, str) for s in sources)
            and isinstance(vendors, list) and all(isinstance(v, str) for v in vendors)
            and isinstance(vendor_ids, dict)
            and all(isinstance(k, str) and isinstance(v, str) for k, v in vendor_ids.items())):
        return None
    sources = sorted(set(s.lower() for s in sources))
    rec = SourcingRecord(inchikey=key, sources=sources, n_sources=len(sources))
    rec.vendor_ids = {s: i for s, i in vendor_ids.items() if s in VENDOR_SOURCES and i}
    rec.vendors = sorted(set(v for v in vendors if any(p in v.lower() for p in VENDOR_NAME_PATTERNS)))
    rec.n_vendors = len(rec.vendors)
    rec.purchasable = bool(rec.vendor_ids or rec.vendors)
    rec.status = "purchasable" if rec.purchasable else "unknown"
    # Positive patent evidence survives; absence was never bound to a policy.
    patent_ids = hit.get("patent_ids", [])
    if isinstance(patent_ids, list) and all(isinstance(v, str) for v in patent_ids):
        rec.patent_ids = patent_ids[:5] if any(s in sources for s in PATENT_SOURCES) else []
    rec.in_patents = True if rec.patent_ids else None
    rec.food_history = any(s in sources for s in FOOD_SOURCES)
    rec.tox_data = any(s in sources for s in TOX_SOURCES)
    rec.is_drug = any(s in sources for s in DRUG_SOURCES)
    rec.pubchem_cid = hit.get("pubchem_cid") if isinstance(hit.get("pubchem_cid"), int) else None
    cas = hit.get("cas", [])
    if isinstance(cas, list) and all(isinstance(v, str) for v in cas):
        rec.cas = list(cas)
        rec.cas_checked = hit.get("cas_checked") is True
    try:
        original_sha = hashlib.sha256(json.dumps(hit, sort_keys=True,
            separators=(",", ":"), ensure_ascii=False, allow_nan=False).encode()).hexdigest()
    except (ValueError, TypeError):
        return None
    rec.cache_provenance = {
        "migration": "legacy-positive-source-evidence-v1",
        "original_record_sha256": original_sha,
        "original_policy_status": "unversioned", "refresh_status": "unavailable",
        "interpretation_policy_sha256": classification_policy(),
        "negative_verdicts_reused": False,
    }
    _apply_notes(rec)
    rec.ip_note = ("Historical source evidence was reinterpreted under the current policy; "
                   "the current lookup did not answer. Supplier records are not live stock confirmation. "
                   + (rec.ip_note if rec.in_patents else "Patent-literature absence is not established."))
    return rec


def lookup(smiles: str, session: Optional[requests.Session] = None,
           cache: Optional[SourcingCache] = None,
           limiter: Optional[RateLimiter] = None) -> Optional[SourcingRecord]:
    """Sourcing and IP status for one structure. ``None`` if it cannot be parsed."""
    key = inchikey_for(smiles)
    if key is None:
        log.warning("could not build an InChIKey for %r", smiles[:60])
        return None
    cache = cache if cache is not None else SourcingCache()
    raw_hit = cache.get(key)
    hit = _usable_hit(raw_hit, key)
    if hit is not None:
        rec = SourcingRecord(**hit)
        if not _cas_unresolved(hit):
            return rec
        # Refresh unresolved CAS metadata even for cached structures. Cache a
        # confirmed absence, but leave failed lookups retryable.
        session = session or requests.Session()
        cas, answered = pubchem_cas(smiles, session, limiter=limiter)
        if not answered:
            log.info("%s: legacy cache entry has no CAS and the top-up lookup "
                     "did not answer; leaving it for the next run", key)
            return rec
        rec.cas, rec.cas_checked = cas, True
        _apply_notes(rec)
        cache.put(key, rec.to_dict())
        return rec
    session = session or requests.Session()
    payload = _query(key, session, limiter=limiter)
    if payload is None:
        legacy = _legacy_source_evidence(raw_hit, key)
        if legacy is not None:
            return legacy
        return SourcingRecord(inchikey=key, status="unknown",
                              ip_note="UniChem lookup failed; status not established")
    rec = _interpret(key, payload)
    cid, vendors, answered = _pubchem_vendors(smiles, session, limiter=limiter)

    if not answered and not rec.purchasable:
        # PubChem is the broader vendor check and it stayed silent, so "nobody
        # sells it" is not something this run established. Say unknown, and do
        # not cache -- a cached negative born of a 503 outlives the outage.
        rec.status = "unknown"
        _apply_notes(rec)
        return rec

    rec.pubchem_cid = cid
    rec.vendors = vendors
    rec.n_vendors = len(vendors)
    # The registry numbers, from the same source and in the same round trip.
    # Without them `eligibility` has only the structure and the name to check
    # against 별표 1·2, and that list carries a structure for 15 of its 626 rows
    # against a CAS for 619 -- so "not on the list" meant "not looked up".
    rec.cas, rec.cas_checked = pubchem_cas(smiles, session, limiter=limiter)
    if vendors and not rec.purchasable:
        rec.purchasable = True
        rec.status = "purchasable"
    _apply_notes(rec)
    if not answered:
        # UniChem already found a vendor, so the status stands; only the vendor
        # names are incomplete. Not worth caching a partial answer.
        return rec
    if not rec.cas_checked:
        # The same rule as the vendor branch above, which this used to skip: a
        # negative born of an outage must not outlive it. `pubchem_cas` separates
        # "PubChem answered and has no registry number" -- `([], True)`, a real
        # answer and cacheable -- from "the lookup failed", `([], False)`. Only
        # the second lands here, and caching it freezes `cas_checked=False` for
        # the life of the entry, so no later run retries however healthy PubChem
        # has become. Since `pipeline.orderable` now holds a compound whose CAS
        # was never checked, a cached failure is not a silent gap any more: it
        # keeps a buyable compound off the buy list until the entry expires.
        log.info("%s: vendors answered but the CAS lookup did not; not caching "
                 "so the next run retries", key)
        return rec
    cache.put(key, rec.to_dict())
    return rec


def lookup_many(smiles_list: Iterable[str], pause: float = 0.0,
                cache: Optional[SourcingCache] = None,
                progress: bool = False,
                workers: int = DEFAULT_WORKERS,
                rate: float = DEFAULT_RATE,
                limiter: Optional[RateLimiter] = None,
                save_every: int = DEFAULT_SAVE_EVERY) -> dict[str, SourcingRecord]:
    'Sourcing status for a library, cached and gentle on the server.'
    cache = cache if cache is not None else SourcingCache()
    limiter = limiter if limiter is not None else RateLimiter(rate)

    # One lookup per distinct structure, first-seen order preserved. The serial
    # loop overwrote a repeated SMILES with an identical record, so the mapping
    # is the same and the duplicate round trip is not.
    unique: list[str] = []
    seen: set[str] = set()
    for smi in smiles_list:
        if smi not in seen:
            seen.add(smi)
            unique.append(smi)

    n_workers = max(1, min(int(workers), len(unique) or 1))
    local = threading.local()
    opened: list[requests.Session] = []
    opened_lock = threading.Lock()

    def _session() -> requests.Session:
        """A Session per thread. requests.Session is not thread-safe; sharing
        one across a pool corrupts connection state under load."""
        s = getattr(local, "session", None)
        if s is None:
            s = requests.Session()
            local.session = s
            with opened_lock:
                opened.append(s)
        return s

    def _one(smi: str) -> Optional[SourcingRecord]:
        rec = lookup(smi, session=_session(), cache=cache, limiter=limiter)
        if pause:
            time.sleep(pause)
        return rec

    results: dict[str, SourcingRecord] = {}
    results_lock = threading.Lock()
    done = 0

    def _record(smi: str, rec: Optional[SourcingRecord]) -> None:
        nonlocal done
        with results_lock:
            if rec is not None:
                results[smi] = rec
            done += 1
            n, total = done, len(unique)
        if progress and n % 25 == 0:
            log.info("sourcing %d/%d", n, total)
        # Periodic flush: a run that dies at minute forty keeps what it learnt.
        # A flush that cannot write is a fact about the disk. Raising here ended
        # the whole batch from inside one compound's bookkeeping, after the
        # lookups had already been paid for; `_dirty` keeps the unwritten keys
        # so the next successful save still carries them.
        if save_every and n % save_every == 0:
            try:
                cache.save()
            except OSError as exc:
                log.warning("could not flush the sourcing cache (%s); "
                            "continuing, and the entries stay queued", exc)

    try:
        if n_workers == 1:
            for smi in unique:
                _record(smi, _one(smi))
        else:
            with ThreadPoolExecutor(max_workers=n_workers,
                                    thread_name_prefix="sourcing") as pool:
                futures = {pool.submit(_one, smi): smi for smi in unique}
                for fut in as_completed(futures):
                    _record(futures[fut], fut.result())
    finally:
        # Sessions first, and the save in its own guard. The save used to come
        # first and unguarded, so an OSError from it skipped every `close()`
        # *and* replaced the return value -- losing lookups that had already
        # succeeded and leaking the connections that fetched them.
        for s in opened:
            try:
                s.close()
            except Exception:                            # pragma: no cover
                pass
        try:
            cache.save()
        except OSError as exc:
            log.warning("could not write the sourcing cache (%s); the results "
                        "below are still good, they are just not cached", exc)

    # Rebuilt in input order, so the mapping does not depend on completion order.
    return {smi: results[smi] for smi in unique if smi in results}


def annotate(ligands, cache: Optional[SourcingCache] = None, progress: bool = False,
             workers: int = DEFAULT_WORKERS):
    """Attach sourcing status to `LigandRecord` objects in place.

    Writes into ``ligand.props`` so the CSV and report layers pick it up with no
    further wiring, and appends a flag when the compound cannot be bought, since
    that is what stalls a screen after the fact.
    """
    recs = lookup_many([ligand.smiles for ligand in ligands], cache=cache, progress=progress,
                       workers=workers)
    for lig in ligands:
        rec = recs.get(lig.smiles)
        if rec is None:
            continue
        lig.props.update({
            "purchasable": rec.purchasable,
            "n_vendors": rec.n_vendors,
            "vendors": list(rec.vendors),
            "pubchem_cid": rec.pubchem_cid,
            "vendor_id": next(iter(rec.vendor_ids.values()), None),
            "in_patents": rec.in_patents,
            "sourcing_status": rec.status,
            "food_history": rec.food_history,
            "is_drug": rec.is_drug,
            # The freedom-to-operate caveat and the dossier-burden sentence were
            # computed for every compound and reached no output at all, because
            # they were never copied out of the record. They are the two
            # sentences a formulator actually has to read.
            "ip_note": rec.ip_note,
            "dossier_note": rec.dossier_note,
            "sourcing_cache_provenance": rec.cache_provenance,






            "cas": rec.cas[0] if rec.cas else None,
            "cas_all": list(rec.cas),
            "cas_checked": rec.cas_checked,
            "cas_source": "pubchem_synonyms" if rec.cas else None,
        })
        if not rec.purchasable:
            lig.flags.append(f"not purchasable ({rec.status})")
        if rec.in_patents:
            lig.flags.append("in patent literature")
    return ligands
