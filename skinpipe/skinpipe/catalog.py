'Screen what can actually be bought, not what has been published.'
from __future__ import annotations

import json
import logging
import time
from dataclasses import dataclass, field
from typing import Any, Iterable, Optional

import requests

from .types import SEVERITY_ORDER, normalise_severity_gate

log = logging.getLogger(__name__)

PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
_UA = {"User-Agent": "skinpipe/1.0 (cosmetic ingredient discovery)"}
_TIMEOUT = 150

#: PubChem PUG-REST publishes a limit of five requests per second and throttles
#: past it with 503 PUGREST.ServerBusy. Without a governor the funnel sails past
#: the limit and every later call returns an empty body -- a screen that reports
#: no hits because it was throttled looks exactly like a screen that found none,
#: which is the failure this pipeline exists to avoid elsewhere. The gap is
#: enforced globally, and a 503 is retried rather than swallowed.
_MIN_INTERVAL = 0.22
_last_call = 0.0


def _throttle() -> None:
    global _last_call
    wait = _MIN_INTERVAL - (time.monotonic() - _last_call)
    if wait > 0:
        time.sleep(wait)
    _last_call = time.monotonic()


class PubChemBusy(RuntimeError):
    """PubChem throttled us and kept throttling. Better raised than returned empty."""


class PubChemError(RuntimeError):
    """PubChem rejected the request (client error). Retrying will not help."""


class CIDSearchResult(list):
    """List-compatible search output with explicit partial-response evidence."""
    def __init__(self, values=(), *, unanswered=()):
        super().__init__(values)
        self.meta = {"complete": not unanswered, "unanswered": list(unanswered)}


def parse_cids(body: Any, what: str) -> CIDSearchResult:
    identifiers = body.get("IdentifierList") if isinstance(body, dict) else None
    values = identifiers.get("CID") if isinstance(identifiers, dict) else None
    if not isinstance(values, list):
        raise PubChemBusy(f"PubChem returned a malformed identifier list for {what}; search is unanswered")
    valid, invalid, seen = [], [], set()
    for index, value in enumerate(values):
        try:
            if isinstance(value, bool) or not (isinstance(value, int) or
                    isinstance(value, str) and value.strip().isdigit()):
                raise ValueError("not an integer CID")
            cid = int(value)
            if cid <= 0:
                raise ValueError("CID must be positive")
        except (TypeError, ValueError, OverflowError):
            invalid.append({"index": index, "value": value})
        else:
            if cid not in seen:
                valid.append(cid)
                seen.add(cid)
    if invalid:
        log.warning("%s: %d invalid CID row(s), %d valid retained; source response is partial",
                    what, len(invalid), len(valid))
        if not valid:
            raise PubChemBusy(f"PubChem returned no valid CIDs for {what}; response is malformed, not an empty search")
    return CIDSearchResult(valid, unanswered=invalid)


def empty_funnel(cids, *responses) -> list:
    metadata = [("identifiers", getattr(cids, "meta", {})), *responses]
    incomplete = [channel for channel, meta in metadata
                  if meta.get("unanswered") or meta.get("complete") is False]
    if incomplete:
        raise PubChemBusy("No validated search row survived, but source responses were partial "
                          f"({', '.join(incomplete)}); absence of purchasable candidates is not established")
    return []


def _request(session: requests.Session, method: str, url: str,
             retries: int = 5, **kw) -> Optional[requests.Response]:
    delay = 1.0
    for attempt in range(retries):
        _throttle()
        try:
            r = session.request(method, url, headers=_UA, timeout=_TIMEOUT, **kw)
        except requests.RequestException as exc:
            log.debug("pubchem %s attempt %d: %s", url, attempt + 1, exc)
            time.sleep(delay)
            delay *= 2
            continue
        if r.status_code == 200:
            return r
        if r.status_code == 404:
            return None                      # a genuine "no such compound"
        if r.status_code in (429, 500, 502, 503, 504):
            log.info("pubchem busy (HTTP %s), backing off %.1fs", r.status_code, delay)
            time.sleep(delay)
            delay *= 2
            continue
        log.warning("pubchem %s -> HTTP %s", url, r.status_code)
        raise PubChemError(
            f"PubChem returned HTTP {r.status_code} for {url}. "
            "The request was rejected; this will not succeed on retry.")
    raise PubChemBusy(
        f"PubChem kept returning a busy status for {url} after {retries} attempts. "
        "Re-run later rather than treating an empty result as no hits.")


def _json(response: requests.Response, what: str) -> Any:
    """The decoded body, or `PubChemBusy` if it is not JSON at all.

    A 200 whose body will not parse is a transport failure wearing a success
    code: an intercepting proxy, a captive portal or a gateway serving its own
    HTML error page. It is *not* an answer from PubChem, and the one thing it
    must never be read as is "no compound matched".

    That is exactly what happened, because `json.JSONDecodeError` subclasses
    `ValueError` and every caller here catches `ValueError` to mean "the body
    was shaped unexpectedly, treat it as empty". An HTTP 200 carrying a proxy
    error page returned `[]` in 0.4 seconds, the funnel printed "no purchasable
    compound survived" and exited 1 -- while a real outage, which is the *same
    class of event*, correctly raised `PubChemBusy` after backing off. The two
    were reported as opposite kinds of result.

    Decode failures are transport failures. Identifier parsing separately
    retains valid rows and marks malformed rows as partial acquisition;
    neither case is interpreted as a completed empty search.
    """
    try:
        return response.json()
    except json.JSONDecodeError as exc:
        body = (getattr(response, "text", "") or "")[:120].replace("\n", " ")
        raise PubChemBusy(
            f"PubChem returned HTTP 200 for {what} but the body is not JSON "
            f"({exc}); it begins {body!r}. That is a transport failure -- a "
            "proxy or gateway answering in PubChem's place -- not an empty "
            "result set. Re-run later rather than reading it as no hits."
        ) from exc


#: Supplier names that count as commercial availability. Same list the sourcing
#: module uses; kept in one place there.
from .sourcing import VENDOR_NAME_PATTERNS  # noqa: E402


@dataclass
class CatalogHit:
    cid: int
    smiles: str
    mw: Optional[float] = None
    logp: Optional[float] = None
    tpsa: Optional[float] = None
    hbd: Optional[int] = None
    vendors: list[str] = field(default_factory=list)
    n_vendors: int = 0
    shape_tanimoto: Optional[float] = None
    tanimoto_2d: Optional[float] = None


def _session(s: Optional[requests.Session] = None) -> requests.Session:
    return s or requests.Session()


def similarity_search(smiles: str, threshold: int = 80, max_records: int = 500,
                      session: Optional[requests.Session] = None) -> list[int]:
    """PubChem 2D similarity search. Returns CIDs, cheapest possible first pass."""
    s = _session(session)
    try:
        r = _request(s, "POST",
                     f"{PUBCHEM}/compound/fastsimilarity_2d/smiles/cids/JSON",
                     data={"smiles": smiles, "Threshold": threshold, "MaxRecords": max_records})
        if r is None:
            return []
        body = _json(r, "the 2D similarity search")
        return parse_cids(body, "the 2D similarity search")
    except (requests.RequestException, ValueError, KeyError) as exc:
        log.warning("similarity search failed: %s", exc)
        raise PubChemBusy(f"similarity search is unanswered: {exc}") from exc


def substructure_search(smiles: str, max_records: int = 500,
                        session: Optional[requests.Session] = None) -> list[int]:
    """Compounds containing the query as a substructure. Use for a pharmacophore core."""
    s = _session(session)
    try:
        r = _request(s, "POST",
                     f"{PUBCHEM}/compound/fastsubstructure/smiles/cids/JSON",
                     data={"smiles": smiles, "MaxRecords": max_records})
        if r is None:
            return []
        body = _json(r, "the substructure search")
        return parse_cids(body, "the substructure search")
    except (requests.RequestException, ValueError, KeyError) as exc:
        log.warning("substructure search failed: %s", exc)
        raise PubChemBusy(f"substructure search is unanswered: {exc}") from exc


def fetch_properties(cids: Iterable[int], batch: int = 150,
                     session: Optional[requests.Session] = None
                     ) -> tuple[dict[int, dict], dict[str, Any]]:
    """Structures and descriptors in bulk. One request per 150 compounds.

    Returns ``(properties, meta)``. ``meta["unanswered"]`` lists the CIDs this
    call did not get an answer for, which is not the same as the CIDs PubChem
    has nothing for -- see `vendor_lookup` below for why that distinction is
    the whole point.
    """
    s = _session(session)
    ids = [int(c) for c in cids]
    out: dict[int, dict] = {}
    unanswered: list[int] = []
    fields = "MolecularWeight,XLogP,TPSA,HBondDonorCount,HBondAcceptorCount,CanonicalSMILES,SMILES"
    for i in range(0, len(ids), batch):
        chunk = ids[i:i + batch]
        try:
            r = _request(s, "POST", f"{PUBCHEM}/compound/cid/property/{fields}/JSON",
                         data={"cid": ",".join(map(str, chunk))})
            if r is None:
                continue                     # 404: PubChem genuinely has none
            rows = _json(r, f"the property batch at {i}")["PropertyTable"]["Properties"]
        except (requests.RequestException, ValueError, KeyError) as exc:
            log.warning("property batch at %d failed: %s", i, exc)
            unanswered.extend(chunk)
            continue
        # CIDs this batch gave a usable answer about -- including "PubChem has
        # no structure on file", which is an answer. What is left over once a
        # row has failed is the set that cannot be accounted for either way.
        accounted: set[int] = set()
        row_failed = False
        for p in rows:
            # Per row. A malformed row used to abort the loop and take every
            # row after it with it -- see `vendor_lookup`.
            try:
                cid = int(p["CID"])
                smi = (p.get("SMILES") or p.get("CanonicalSMILES")
                       or p.get("ConnectivitySMILES"))
                if not smi:
                    accounted.add(cid)       # answered: no structure on file
                    continue
                out[cid] = {
                    "smiles": smi,
                    "mw": float(p["MolecularWeight"]) if p.get("MolecularWeight") else None,
                    "logp": p.get("XLogP"), "tpsa": p.get("TPSA"),
                    "hbd": p.get("HBondDonorCount"), "hba": p.get("HBondAcceptorCount")}
                accounted.add(cid)
            except (ValueError, KeyError, TypeError) as exc:
                row_failed = True
                log.warning("property row %r in batch at %d: %s",
                            p.get("CID"), i, exc)
        if row_failed:
            unanswered.extend(c for c in chunk if c not in accounted)
    answered = {c for c in ids if c in out}
    unanswered = [c for c in unanswered if c not in answered]
    return out, {"requested": len(ids), "answered": len(out),
                 "unanswered": unanswered}


def vendor_lookup(cids: Iterable[int], batch: int = 50,
                  session: Optional[requests.Session] = None
                  ) -> tuple[dict[int, list[str]], dict[str, Any]]:
    """Suppliers per compound. The expensive step, so run it on a cut list.

    A supplier deposit is evidence the compound has been sold, not a live quote.
    Catalogues go stale. It is still the difference between an orderable list
    and a wish list.

    Returns ``(vendors, meta)``. ``meta["unanswered"]`` is the part that matters:
    a CID absent from ``vendors`` because the lookup failed is not a CID nobody
    sells, and `screen` cuts on exactly that value. `sourcing.pubchem_cas`
    reports the same distinction as `answered`, and `chembl.fetch_molecules` as
    `failed_batch_ids`; this function was the one PubChem caller that did not.
    """
    s = _session(session)
    ids = [int(c) for c in cids]
    out: dict[int, list[str]] = {}
    unanswered: list[int] = []
    for i in range(0, len(ids), batch):
        chunk = ids[i:i + batch]
        try:
            r = _request(s, "POST", f"{PUBCHEM}/compound/cid/xrefs/SourceName/JSON",
                         data={"cid": ",".join(map(str, chunk))})
            if r is None:
                continue                     # 404: PubChem genuinely has none
            rows = _json(r, f"the vendor batch at {i}")["InformationList"]["Information"]
        except (requests.RequestException, ValueError, KeyError) as exc:
            # The batch produced no rows at all. Its CIDs are *unchecked*, and
            # the caller has to be able to tell that from "checked, no vendor":
            # `screen` drops on `len(vendors.get(cid, [])) >= min_vendors`, so a
            # missing entry and a genuine zero are the same compound-killing
            # value. Most PubChem failures never reach here -- a 503 or a
            # non-JSON 200 raises `PubChemBusy` past this handler by design --
            # so what lands here is a body that parsed and was not what was
            # asked for.
            log.warning("vendor batch at %d failed: %s", i, exc)
            unanswered.extend(chunk)
            continue
        # PubChem omits a CID it has no supplier xref for, so in a clean batch
        # "absent" genuinely means "nobody sells it" and must stay that way.
        # Once a row has failed to parse that inference is no longer available:
        # the failed row's own CID is unreadable, so an absent CID is either a
        # real zero or the casualty. Both then have to travel as *unchecked*.
        accounted: set[int] = set()
        row_failed = False
        for info in rows:
            # Per row, not per batch. Measured before this change: one row whose
            # CID PubChem rendered as "n/a" raised ValueError out of the loop
            # and cost 26 of 50 compounds, every one of which then read as
            # "no vendor found" downstream.
            try:
                cid = int(info["CID"])
                names = info.get("SourceName", []) or []
                seen: list[str] = []
                for n in names:
                    low = n.lower()
                    if any(p in low for p in VENDOR_NAME_PATTERNS) and n not in seen:
                        seen.append(n)
                out[cid] = seen
                accounted.add(cid)
            except (ValueError, KeyError, TypeError) as exc:
                row_failed = True
                log.warning("vendor row %r in batch at %d: %s",
                            info.get("CID"), i, exc)
        if row_failed:
            unanswered.extend(c for c in chunk if c not in accounted)
    answered = set(out)
    unanswered = [c for c in unanswered if c not in answered]
    return out, {"requested": len(ids), "answered": len(out),
                 "unanswered": unanswered}


def screen(query_smiles: str, *, threshold: int = 80, max_records: int = 500,
           min_topical: float = 6.0, min_vendors: int = 1,
           shape_top_n: int = 60, n_confs: int = 10,
           admet_gate: Optional[str] = "critical",
           target_gene: Optional[str] = None,
           progress: bool = True) -> "list[dict]":
    'The whole funnel: query in, orderable and shape-matched compounds out.'
    from rdkit import Chem, DataStructs, RDLogger
    from rdkit.Chem import rdFingerprintGenerator

    from . import hopping
    from .filters import admet as admet_f
    from .filters import topical as topical_f
    from .types import LigandRecord

    # Before any network call: a funnel that has already spent its requests is
    # a bad place to discover the gate was never going to exclude anything.
    admet_gate = validate_gate(admet_gate)

    RDLogger.DisableLog("rdApp.*")
    s = requests.Session()
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    qmol = Chem.MolFromSmiles(query_smiles)
    if qmol is None:
        raise ValueError(f"query SMILES could not be parsed: {query_smiles!r}")
    if qmol.GetNumHeavyAtoms() == 0:
        # RDKit parses "" into a perfectly valid molecule with no atoms. Its
        # fingerprint is empty, so every Tanimoto against it is 0.0 and the
        # whole funnel runs its network budget to rank nothing by nothing.
        raise ValueError(
            f"query SMILES {query_smiles!r} has no heavy atoms, so there is "
            "nothing to search for. RDKit parses it, but an empty molecule "
            "has an empty fingerprint and matches nothing meaningfully.")
    qfp = gen.GetFingerprint(qmol)

    cids = similarity_search(query_smiles, threshold=threshold,
                             max_records=max_records, session=s)
    if progress:
        log.info("similarity: %d compounds", len(cids))
    if not cids:
        return empty_funnel(cids)

    props, prop_meta = fetch_properties(cids, session=s)
    if progress:
        log.info("properties: %d retrieved", len(props))
    if prop_meta["unanswered"]:
        log.warning("%d of %d compounds got no property answer and are "
                    "absent from this screen -- not filtered out, unchecked",
                    len(prop_meta["unanswered"]), prop_meta["requested"])
    if not props:
        raise PubChemBusy(
            "the similarity search returned compounds but no properties came back. "
            "That is throttling, not an empty result set; re-run later.")

    # Skin filter first: local, free, and it removes most of the list.
    kept: dict[int, dict] = {}
    # A compound PubChem has no XLogP for was not judged and failed nothing.
    # Counted separately because the progress line below read "skin filter:
    # N/M pass", which filed "the descriptor never arrived" under "this
    # compound does not permeate skin" -- a property of PubChem's coverage
    # reported as a property of the molecule. Measured live at 1.3% of a real
    # similarity result (2 of 150), so it is small, not zero.
    unscored: list[int] = []
    for cid, p in props.items():
        if p["mw"] is None or p["logp"] is None:
            unscored.append(cid)
            continue
        lig = LigandRecord(ligand_id=str(cid), smiles=p["smiles"], props=p)
        a = topical_f.assess(lig)
        if a.topical_score >= min_topical:
            kept[cid] = p | {"topical_score": a.topical_score,
                             "log_kp": a.log_kp,
                             "permeation": a.permeation_class}
    if progress:
        log.info("skin filter: %d/%d pass topical >= %.1f",
                 len(kept), len(props) - len(unscored), min_topical)
    if unscored:
        log.warning("%d of %d compounds carry no PubChem MW/XLogP and were not "
                    "scored against the skin filter at all; they are absent from "
                    "the funnel, not rejected by it", len(unscored), len(props))
    if not kept:
        return empty_funnel(cids, ("properties", prop_meta),
                            ("skin descriptors", {"unanswered": unscored}))

    vendors, vendor_meta = vendor_lookup(kept, session=s)
    # `unchecked` is carried onto the row rather than folded into
    # `n_vendors`, because a compound the lookup never reached must not
    # read as one nobody sells. It is still cut here -- an unchecked
    # compound cannot go on a buy list -- but it leaves saying which it is.
    unchecked = set(vendor_meta["unanswered"])
    if unchecked:
        log.warning("%d of %d vendor lookups did not answer; for those "
                    "compounds 'not purchasable' means 'not checked'",
                    len(unchecked), vendor_meta["requested"])
    buyable = {cid: d | {"vendors": vendors.get(cid, []),
                         "n_vendors": len(vendors.get(cid, [])),
                         "vendor_status": ("unchecked" if cid in unchecked
                                          else "checked")}
               for cid, d in kept.items()
               if cid not in unchecked
               and len(vendors.get(cid, [])) >= min_vendors}
    if progress:
        log.info("vendors: %d/%d purchasable", len(buyable), len(kept))
    if not buyable:
        return empty_funnel(cids, ("properties", prop_meta), ("vendors", vendor_meta),
                            ("skin descriptors", {"unanswered": unscored}))

    # 2D similarity is free and orders the list for the expensive 3D step.
    for cid, d in buyable.items():
        m = Chem.MolFromSmiles(d["smiles"])
        d["tanimoto_2d"] = (DataStructs.TanimotoSimilarity(qfp, gen.GetFingerprint(m))
                            if m is not None else 0.0)
    ranked = sorted(buyable.items(), key=lambda kv: -kv[1]["tanimoto_2d"])[:shape_top_n]

    out: list[dict] = []
    for cid, d in ranked:
        shape = hopping.best_shape_from_smiles(query_smiles, d["smiles"], n_confs=n_confs)
        if shape is None:
            continue
        # Pass the target through: a hazard the target itself activates cannot
        # be judged without knowing the target. See admet.MECHANISM_ESCALATION.
        prof = admet_f.summarise(d["smiles"], ligand_id=str(cid),
                                 target_gene=target_gene)
        hazard = _hazard(prof)
        if admet_gate and _at_least(hazard, admet_gate):
            continue
        out.append({"cid": cid, "smiles": d["smiles"], "mw": d["mw"], "logp": d["logp"],
                    "tpsa": d.get("tpsa"), "log_kp": d.get("log_kp"),
                    "permeation": d.get("permeation"),
                    "topical_score": d.get("topical_score"),
                    "tanimoto_2d": round(d["tanimoto_2d"], 3),
                    "shape_tanimoto": round(shape, 3),
                    "n_vendors": d["n_vendors"], "vendors": d["vendors"][:4],
                    # A hazard severity without the target it was judged against
                    # is not readable: the same phenol is 'context' off a
                    # pigmentation programme and 'high' on one. Both travel.
                    "target_gene": target_gene or "",
                    # `admet_gate` is the *normalised* gate, so "none" here
                    # means what `validate_gate("none")` means -- no gate --
                    # and a row read back out of the CSV configures the same
                    # screen that produced it. Before the validator was shared
                    # this column recorded "none" for a gate-off run and that
                    # spelling, fed back in, dropped every compound.
                    "admet_gate": admet_gate or "none",
                    "admet_hazard": hazard, "admet_overall": prof.overall})
    out.sort(key=lambda r: -(r["shape_tanimoto"] or 0))
    if progress:
        log.info("shape + ADMET: %d survivors", len(out))
    return out if out else empty_funnel(cids, ("properties", prop_meta), ("vendors", vendor_meta),
                                       ("skin descriptors", {"unanswered": unscored}))


#: The severity ladder, taken from `types` rather than restated. A second copy
#: is a second thing to forget to update, and this one is indexed against a
#: gate the pipeline validates with its own list.
_SEV = tuple(SEVERITY_ORDER)


def _hazard(profile) -> str:
    """Worst hazard severity, ignoring systemic exposure (a dose, not a liability)."""
    worst = "none"
    for c in getattr(profile, "concerns", []):
        if c.domain == "systemic" or c.severity not in _SEV:
            continue
        if _SEV.index(c.severity) < _SEV.index(worst):
            worst = c.severity
    return worst


def validate_gate(gate: Optional[str]) -> Optional[str]:
    'The gate to apply, or `None` for "no gate". Raises on anything else.'
    return normalise_severity_gate(gate, field="ADMET gate")


def _at_least(severity: str, gate: Optional[str]) -> bool:
    """Is `severity` at least as bad as `gate`?

    An unknown *severity* is not a hazard this module ranks, so it does not trip
    the gate. An unknown *gate* is an input error and is refused: silently
    answering False for it is what left the gate off for a whole run. `"none"`
    and `None` are the off switch, not a threshold, and gate nothing.
    """
    gate = validate_gate(gate)
    if gate is None:
        return False
    return severity in _SEV and _SEV.index(severity) <= _SEV.index(gate)
