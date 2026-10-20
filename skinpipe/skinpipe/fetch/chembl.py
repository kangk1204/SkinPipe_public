"""ChEMBL client for target mapping and bioactivity retrieval.

Every activity that comes out of here carries a provenance call, because for
skin targets the raw database is not safe to rank on directly. See
`skinpipe.fetch.provenance` for why.
"""
from __future__ import annotations

import logging
import math
import time
import hashlib
import json
from typing import Any, Iterable, Optional

import requests

from . import provenance as prov
from . import cache as cache_mod
from .cache import DiskCache, FetchError, cached
from . import chembl_local

log = logging.getLogger(__name__)

#: Shared disk cache. Pass `cache=` to override, `cache=False` to bypass.
_CACHE: Optional[DiskCache] = None


def cache() -> DiskCache:
    global _CACHE
    if _CACHE is None:
        _CACHE = DiskCache()
    return _CACHE


def _c(arg):
    """None -> the shared cache; False -> no caching; a DiskCache -> itself."""
    return None if arg is False else (arg if arg is not None else cache())

CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"
UNIPROT = "https://rest.uniprot.org/uniprotkb"
_UA = {"User-Agent": "skinpipe/1.0 (cosmetic ingredient discovery)"}
_TIMEOUT = 90


class ChemblError(FetchError):
    """A ChEMBL/UniProt/RCSB call did not answer.

    Inherits `cache.FetchError` so `cache.cached` can tell "the service is
    down" -- which a stale entry may legitimately answer for -- from a defect
    in this package, which it must not.
    """


def _get(url: str, params: dict[str, Any], retries: int = 3,
         *, timeout: float | tuple[float, float] = _TIMEOUT) -> dict:
    """GET and decode JSON, retrying the transient failures.

    The status code is checked before `.json()` is touched: ChEMBL answers 500
    with an HTML error page, and decoding that as JSON raises somewhere far from
    the cause. A 200 carrying something that is not JSON is treated as a flake
    too -- `requests.JSONDecodeError` is a `RequestException` -- so it retries
    rather than surfacing as a decode error in the middle of a screen.
    """
    last: Optional[str] = None
    for attempt in range(retries):
        try:
            r = requests.get(url, params=params, headers=_UA, timeout=timeout)
            if r.status_code == 200:
                return r.json()
            if r.status_code in (429, 500, 502, 503, 504):
                last = f"HTTP {r.status_code}"
                if attempt + 1 < retries:
                    time.sleep(2 ** attempt)
                continue
            raise ChemblError(f"{url} -> HTTP {r.status_code}")
        except requests.RequestException as e:      # network flake, or bad JSON
            last = f"{type(e).__name__}: {e}"
            if attempt + 1 < retries:
                time.sleep(2 ** attempt)
    raise ChemblError(f"{url} failed after {retries} attempts; last failure: {last}")


def uniprot_for_gene(gene: str, organism_id: int = 9606, cache=None) -> Optional[dict]:
    """One reviewed entry with this exact primary gene, never an alias hit."""
    from .target_identity import validate_query, _ACCESSION, _GENE
    gene, kind = validate_query(gene)
    if not _GENE.fullmatch(gene) or kind == "chembl_target":
        raise ChemblError("A primary gene symbol is required")
    q = f"gene_exact:{gene} AND organism_id:{organism_id} AND reviewed:true"
    params = {"query": q, "format": "json", "size": 100}
    data = cached(_c(cache), "uniprot", params,
                  lambda: _get(f"{UNIPROT}/search", params)).value
    rows = data.get("results")
    if not isinstance(rows, list) or len(rows) >= 100:
        raise ChemblError("UniProt result list is invalid or incomplete")
    entries = {}
    for row in rows:
        if not isinstance(row, dict):
            raise ChemblError("UniProt returned an invalid protein record")
        names = {g.get("geneName", {}).get("value", "").upper() for g in row.get("genes", [])}
        accession = row.get("primaryAccession", "")
        if (names == {gene} and _ACCESSION.fullmatch(accession) and
                row.get("entryType") == "UniProtKB reviewed (Swiss-Prot)" and
                str((row.get("organism") or {}).get("taxonId")) == str(organism_id)):
            entries[accession] = row
    if len(entries) > 1:
        raise ChemblError(f"Multiple reviewed primary-gene entries for {gene}; specify an accession")
    if not entries:
        log.warning("no reviewed UniProt entry for gene %s", gene)
        return None
    e = next(iter(entries.values()))
    name = (e.get("proteinDescription", {}).get("recommendedName", {})
             .get("fullName", {}).get("value"))
    return {"accession": e["primaryAccession"], "id": e.get("uniProtkbId"),
            "name": name, "length": e.get("sequence", {}).get("length")}


def chembl_target_for_accession(accession: str, organism: str = "Homo sapiens",
                                cache=None) -> Optional[str]:
    """One verified single-protein component; ambiguous mappings need a choice."""
    params = {"target_components__accession": accession,
              "target_type": "SINGLE PROTEIN", "limit": 100}
    data = cached(_c(cache), "target", params,
                  lambda: _get(f"{CHEMBL}/target.json", params)).value
    from .target_identity import _ACCESSION, _TARGET
    rows = data.get("targets")
    page = data.get("page_meta") or {}
    if (not isinstance(rows, list) or page.get("next") or
            int(page.get("total_count", len(rows))) != len(rows)):
        raise ChemblError("ChEMBL target mapping is invalid or incomplete")
    ids = set()
    for row in rows:
        if not isinstance(row, dict):
            raise ChemblError("ChEMBL returned an invalid target record")
        components = row.get("target_components") or []
        tid = row.get("target_chembl_id", "")
        if (row.get("target_type") != "SINGLE PROTEIN" or len(components) != 1 or
                components[0].get("accession") != accession or not _ACCESSION.fullmatch(accession) or
                not _TARGET.fullmatch(tid)):
            raise ChemblError("ChEMBL returned an unexpected target component or type")
        if row.get("organism") != organism:
            raise ChemblError("ChEMBL returned an unexpected organism")
        if organism == "Homo sapiens" and str(row.get("tax_id")) != "9606":
            raise ChemblError("Human ChEMBL target taxonomy was not verified")
        ids.add(tid)
    if not ids:
        return None
    if len(ids) > 1:
        raise ChemblError(f"Accession {accession} maps to multiple targets; specify an exact target ID")
    return next(iter(ids))


def _activity_filter(target_chembl_id, pchembl_min, include_without_pchembl=False):
    params = {"target_chembl_id": target_chembl_id}
    if not include_without_pchembl:
        params.update(pchembl_value__isnull="false", pchembl_value__gte=pchembl_min)
    return params


def activity_count(target_chembl_id: str, pchembl_min: float = 6.0, cache=None,
                   include_without_pchembl: bool = False) -> int:
    params = _activity_filter(target_chembl_id, pchembl_min, include_without_pchembl)
    params["limit"] = 1
    data = cached(_c(cache), "activity_count", params,
                  lambda: _get(f"{CHEMBL}/activity.json", params)).value
    return int(data["page_meta"]["total_count"])


def fetch_activities(target_chembl_id: str, pchembl_min: float = 6.0,
                     max_records: int = 2000, target_gene: Optional[str] = None,
                     page: int = 200, with_meta: bool = False, cache=None,
                     include_without_pchembl: bool = False):
    """Potent activities for a target, each tagged with an assay provenance call.

    A positive `max_records` caps the pull; zero reads the full live query.
    Pages are ordered by activity identity, never by potency. A finite slice
    therefore cannot be interpreted as the most potent compounds.

    With `with_meta=True` returns `(records, meta)`. The meta carries the
    server's own `total_count` against what was actually retrieved, because a
    short page under load ends the pagination loop and leaves a truncated set
    that no downstream ratio can detect: every statistic computed from these
    records shrinks with them, so the result stays self-consistent and wrong.

    Cached on the full query -- cutoff, cap and field list all change the
    answer. Complete raw pages are checkpointed with their original retrieval
    times even if a later page fails. A partial whole query is never cached.
    """
    if type(max_records) is not int or max_records < 0:
        raise ValueError("max_records must be a nonnegative integer; zero means all")
    if type(page) is not int or page < 1:
        raise ValueError("page must be a positive integer")
    fields = ("activity_id,molecule_chembl_id,pchembl_value,standard_type,standard_value,"
              "standard_units,standard_relation,assay_type,assay_description,"
              "assay_chembl_id,target_organism,document_chembl_id,"
              "activity_comment,data_validity_comment")
    # Broad acquisition includes inactive/inconclusive records and endpoint
    # properties. Retain the full API record so its comments and validity flags
    # cannot disappear behind a missing pChEMBL cell.
    if include_without_pchembl:
        fields = None
    ckey = {"target": target_chembl_id, "pchembl_min": pchembl_min,
            "max_records": max_records, "gene": target_gene, "page": page,
            "fields": fields, "activity_schema_rev": 4,
            "include_without_pchembl": include_without_pchembl,
            # Bumped when the *meaning* of a stored record changes, not only
            # its fields. v2 classifies on the organism the assay was run in
            # rather than on the target record's organism, so every cached
            # entry written before it carries a provenance call that was made
            # on the wrong evidence. Without this, a cache hit serves those
            # calls unchanged and the fix reaches nobody who has run the
            # screen before.
            "provenance_rev": prov.CLASSIFIER_VERSION}
    res = cached(_c(cache), "activities", ckey,
                 lambda: _activities_live(target_chembl_id, pchembl_min,
                                          max_records, target_gene, page, fields,
                                          cache, include_without_pchembl),
                 cacheable=lambda v: bool(v[1]["complete"])
                 and v[1].get("assay_context_complete") is True
                 and not v[1].get("stale"))
    out, meta = res.value
    pages = meta.get("activity_pages") or []
    meta = dict(meta) | {"from_cache": res.from_cache or bool(pages and all(p.get("from_cache") for p in pages)),
                         "stale": res.stale or bool(meta.get("stale")),
                         "cache_age_days": (round(res.age_days, 1) if res.age_days is not None
                                            else max((p.get("cache_age_days") or 0.0
                                                      for p in meta.get("activity_pages", []) if p.get("stale")), default=None))}
    if meta["stale"]:
        meta["stale_note"] = res.note or (
            "Some raw activity pages were served from stale checkpoints because "
            "the live fetch failed; the records may be out of date.")
    return (out, meta) if with_meta else out


#: Assay ids per `assay_chembl_id__in` request. ChEMBL accepts a comma-joined
#: list; 50 keeps the URL well inside any gateway's limit.
_ASSAY_BATCH = 50


def assay_organisms(assay_ids: "Iterable[str]", cache=None) -> "dict[str, Optional[str]]":
    'The organism each assay was actually run in, keyed by assay id.'
    return {aid: row.get("assay_organism")
            for aid, row in assay_metadata(assay_ids, cache).items()}


def assay_metadata(assay_ids: Iterable[str], cache=None) -> dict[str, dict]:
    """Join context once per assay, preserving raw species and cell metadata.

    A failed or omitted assay stays absent; target organism is never used as a
    substitute. The separate cache kind leaves older raw joins untouched.
    """
    ids = [i for i in dict.fromkeys(assay_ids) if i]
    out: dict[str, dict] = {}
    for start in range(0, len(ids), _ASSAY_BATCH):
        chunk = ids[start:start + _ASSAY_BATCH]
        ckey = {"assays": chunk, "schema_version": 2}
        try:
            res = cached(_c(cache), "assay_context", ckey,
                         lambda chunk=chunk: _assay_metadata_live(chunk),
                         cacheable=lambda value, chunk=chunk: set(value) == set(chunk))
            # A stale context can still be displayed as evidence, but cannot
            # certify a freshly completed provenance join.
            out.update({aid: dict(row, _context_stale=True) if res.stale is True else row
                        for aid, row in res.value.items()})
        except ChemblError as exc:
            log.warning("assay organism lookup failed for %d assay(s) (%s); "
                        "those rows will be classified without it, which reads "
                        "as unknown rather than as human", len(chunk), exc)
    return out


def _assay_organisms_live(assay_ids: "list[str]") -> "dict[str, Optional[str]]":
    return {aid: row.get("assay_organism")
            for aid, row in _assay_metadata_live(assay_ids).items()}


def _assay_metadata_live(assay_ids: list[str]) -> dict[str, dict]:
    data = _get(f"{CHEMBL}/assay.json",
                {"assay_chembl_id__in": ",".join(assay_ids),
                 "limit": len(assay_ids),
                 "only": "assay_chembl_id,assay_organism,assay_cell_type,assay_tissue,assay_type,confidence_score,bao_format"})
    rows = data.get("assays")
    if (not isinstance(rows, list) or any(not isinstance(a, dict) for a in rows)
            or (data.get("page_meta") or {}).get("next")):
        raise ChemblError("Assay context response is invalid or paginated")
    identities = [a.get("assay_chembl_id") for a in rows]
    if (any(aid not in assay_ids for aid in identities)
            or len(set(identities)) != len(identities)):
        raise ChemblError("Assay context response has unexpected or duplicate identities")
    return {a["assay_chembl_id"]: dict(a) for a in rows}


def _validate_activity_page(data: Any, params: dict) -> int:
    """Validate before checkpointing; missing totals never establish coverage."""
    if not isinstance(data, dict) or not isinstance(data.get("page_meta"), dict):
        raise ChemblError("Activity response has no page metadata")
    metadata = data["page_meta"]
    total = metadata.get("total_count")
    if isinstance(total, str) and total.isascii() and total.isdigit():
        total = int(total)
    if type(total) is not int or total < 0:
        raise ChemblError("Activity response has no verified nonnegative total_count")
    rows = data.get("activities")
    if (not isinstance(rows, list) or any(not isinstance(row, dict) for row in rows)
            or len(rows) > params["limit"] or params["offset"] + len(rows) > total):
        raise ChemblError("Activity response violates the requested page bounds")
    for field in ("offset", "limit"):
        if field in metadata and (type(metadata[field]) is not int or metadata[field] != params[field]):
            raise ChemblError(f"Activity response has an unexpected {field}")
    return total


def _activity_page(params: dict, cache=None):
    """Retain one complete raw response without claiming a complete query."""
    def fetch():
        response = _get(f"{CHEMBL}/activity.json", params)
        _validate_activity_page(response, params)
        return {"response": response, "fetched_at": time.time()}

    def complete(value):
        total = _validate_activity_page(value["response"], params)
        return len(value["response"]["activities"]) == min(params["limit"], max(0, total-params["offset"]))

    try:
        result = cached(_c(cache), "activity_page", {"schema": 1, "params": params}, fetch,
                        cacheable=complete)
        value = result.value
        _validate_activity_page(value["response"], params)
        if (type(value.get("fetched_at")) not in (int, float)
                or not math.isfinite(value["fetched_at"]) or value["fetched_at"] <= 0):
            raise ChemblError("Activity checkpoint has no valid original retrieval time")
        return result
    except (ChemblError, KeyError, TypeError) as error:
        raise ChemblError(f"Activity {params['target_chembl_id']} offset={params['offset']} "
                          f"limit={params['limit']}: {error}") from error


def _activities_live(target_chembl_id: str, pchembl_min: float, max_records: int,
                     target_gene: Optional[str], page: int,
                     fields: Optional[str], cache=None,
                     include_without_pchembl: bool = False) -> tuple[list[dict], dict]:
    out: list[dict] = []
    raw: list[dict] = []
    raw_times: list[float] = []
    page_receipts: list[dict] = []
    pagination_consistent = True
    offset = 0
    total: Optional[int] = None
    while max_records == 0 or offset < max_records:
        params = _activity_filter(target_chembl_id, pchembl_min, include_without_pchembl)
        params.update(limit=page if max_records == 0 else min(page, max_records-offset),
                      offset=offset, order_by="activity_id")
        if fields is not None:
            params["only"] = fields
        fetched = _activity_page(params, cache)
        data = fetched.value["response"]
        page_total = _validate_activity_page(data, params)
        if total is None:
            total = page_total
        pagination_consistent = pagination_consistent and page_total == total
        batch = data["activities"]
        page_receipts.append({"params": dict(params), "retrieved": len(batch), "total_count": page_total,
            "fetched_at": fetched.value["fetched_at"], "from_cache": fetched.from_cache,
            "stale": fetched.stale, "cache_age_days": fetched.age_days,
            "raw_response_sha256": hashlib.sha256(json.dumps(data, sort_keys=True, separators=(",", ":"),
                                                            ensure_ascii=False).encode()).hexdigest()})
        if not batch:
            break
        raw.extend(batch)
        raw_times.extend([fetched.value["fetched_at"]] * len(batch))
        offset += len(batch)
        if len(batch) < params["limit"] or offset >= page_total:
            break

    # The organism is resolved for the whole pull at once, before anything is
    # classified, because it lives on the assay and not on the activity row.
    contexts = assay_metadata((a.get("assay_chembl_id") for a in raw), cache)
    requested_assays = {a.get("assay_chembl_id") for a in raw if a.get("assay_chembl_id")}
    unresolved_assays = sorted(requested_assays - contexts.keys())
    stale_assays = sorted(aid for aid in requested_assays
                          if contexts.get(aid, {}).get("_context_stale"))
    n_missing_assay_ids = sum(not a.get("assay_chembl_id") for a in raw)
    assay_context_complete = not (unresolved_assays or stale_assays or n_missing_assay_ids)
    fetched_at = time.time()
    for row_index, a in enumerate(raw):
        context = contexts.get(a.get("assay_chembl_id")) or {}
        organism = context.get("assay_organism")
        call = prov.classify(a.get("assay_description"),
                             organism, target_gene,
                             target_organism=a.get("target_organism"),
                             assay_type=a.get("assay_type") or context.get("assay_type"),
                             assay_cell_type=context.get("assay_cell_type"),
                             assay_tissue=context.get("assay_tissue"))
        try:
            pchembl = float(a.get("pchembl_value"))
            if not math.isfinite(pchembl) or isinstance(a.get("pchembl_value"), bool):
                pchembl = None
        except (ValueError, TypeError):
            pchembl = None
        out.append({
            "activity_id": a.get("activity_id"),
            "target_chembl_id": target_chembl_id,
            "target_gene": target_gene,
            "molecule_chembl_id": a.get("molecule_chembl_id"),
            "pchembl": pchembl,
            "source_row_index": row_index,
            "standard_type": a.get("standard_type"),
            "standard_value": a.get("standard_value"),
            "standard_units": a.get("standard_units"),
            "activity_comment": a.get("activity_comment"),
            "data_validity_comment": a.get("data_validity_comment"),
            # Asked for in `fields` and then dropped on the floor. It is
            # the qualifier: a ">" or "<" means the number is a bound, not
            # a measurement, and treating one as exact is a scientific
            # error nothing downstream could have caught because the field
            # never arrived. ChEMBL only computes `pchembl_value` for
            # `standard_relation = '='` -- measured live, 295/295 rows for
            # CHEMBL1973 -- so this should always be "=", and `_relation_
            # census` below says so out loud instead of assuming it.
            "standard_relation": a.get("standard_relation"),
            "assay_type": a.get("assay_type"),
            "assay_description": a.get("assay_description"),
            "assay_chembl_id": a.get("assay_chembl_id"),
            "document_chembl_id": a.get("document_chembl_id"),






















            "assay_organism": organism,
            # Kept beside it rather than replaced: the two disagreeing is
            # itself evidence, and a reader comparing them can see which
            # call was made and why.
            "target_organism": a.get("target_organism"),
            "provenance": call.label,
            "provenance_detail": call.detail,
            "provenance_trustworthy": call.trustworthy,




            "provenance_basis": call.basis,
            "assay_cell_type": context.get("assay_cell_type"),
            "assay_tissue": context.get("assay_tissue"),
            "assay_confidence_score": context.get("confidence_score"),
            "assay_target_confidence": context.get("confidence_score"),
            "assay_context_status": ("UNAVAILABLE" if a.get("assay_chembl_id") not in contexts
                                     else "STALE" if context.get("_context_stale") else "RESOLVED"),
            "bao_format": context.get("bao_format"),
            "raw_activity": dict(a),
            "raw_assay_context": {k: v for k, v in context.items() if k != "_context_stale"},
            "fetched_at": raw_times[row_index],
            "source_url": f"{CHEMBL}/activity.json?target_chembl_id={target_chembl_id}",
            **call.record_fields(),
        })
    log.info("target %s: %d activities at pChEMBL>=%.1f",
             target_chembl_id, len(out), pchembl_min)
    expected = total if max_records == 0 else min(total, max_records) if total is not None else None
    ids = [row.get("activity_id") for row in raw]
    valid_ids = [str(value) for value in ids if type(value) in (int, str) and str(value)]
    activity_id_unique = (len(valid_ids) == len(set(valid_ids))) if len(valid_ids) == len(ids) else None
    complete = (expected is not None and len(out) == expected and pagination_consistent
                and activity_id_unique is not False and (max_records != 0 or activity_id_unique is True))
    if not complete:
        log.warning("target %s: retrieved %d of %s expected activities; pagination/identity complete=%s",
                    target_chembl_id, len(out), expected, complete)
    censored = sorted({str(r["standard_relation"]) for r in out
                       if r["standard_relation"] not in (None, "", "=")})
    if censored:
        log.warning("target %s: %d of %d activities carry a censored relation "
                    "(%s); their pChEMBL is a bound, not a measurement",
                    target_chembl_id,
                    sum(1 for r in out if r["standard_relation"] not in (None, "", "=")),
                    len(out), ", ".join(censored))
    return out, {"total_count": total, "retrieved": len(out),
                 "activity_scope": "all_recorded" if include_without_pchembl else "pchembl_floor",
                 "include_without_pchembl": include_without_pchembl,
                 "pchembl_min": None if include_without_pchembl else pchembl_min,
                 "query_filter": _activity_filter(target_chembl_id, pchembl_min, include_without_pchembl),
                 "selected_activity_fields": fields,
                 "requested_max": max_records,
                 "provenance_classifier_version": prov.CLASSIFIER_VERSION,
                 "fetched_at": fetched_at,
                 "expected": expected, "complete": complete,
                 "assay_context_complete": assay_context_complete,
                 "n_unresolved_assays": len(unresolved_assays),
                 "unresolved_assay_ids": unresolved_assays,
                 "stale_assay_ids": stale_assays,
                 "n_activities_without_assay_id": n_missing_assay_ids,
                 "assay_context_interpretation": (
                     "Missing assay context is an acquisition limitation, not evidence of absent human activity."),
                 "capped_at_max_records": bool(max_records > 0 and total is not None and total > max_records),
                 "pagination_consistent": pagination_consistent, "activity_id_unique": activity_id_unique,
                 "activity_pages": page_receipts, "stale": any(p["stale"] for p in page_receipts),
                 # Expected to be empty. ChEMBL computes pchembl_value only for
                 # an exact relation, so anything here means that assumption
                 # has stopped holding and the potencies are not comparable.
                 "censored_relations": censored,
                 "n_censored": sum(1 for r in out
                                   if r["standard_relation"] not in (None, "", "="))}


def fetch_molecules(chembl_ids: Iterable[str], batch: int = 50,
                    with_meta: bool = False, cache=None):
    """Structures and calculated properties for a set of ChEMBL molecules.

    Cached per molecule, not per request. A compound's structure does not
    depend on which screen asked for it, so the entry is reused across every
    target that touches it -- and these fetches are the slowest part of a run.

    A batch that fails is logged and skipped, because one bad chunk should not
    lose a whole screen. With `with_meta=True` the ids that never arrived are
    returned too: without them the caller sees a smaller library and no reason,
    which is how a screen ends up ranking a tenth of its compounds and saying
    nothing about the other nine tenths. Failed ids are never cached.
    """
    ids = sorted(set(chembl_ids))
    fields = ("molecule_chembl_id,pref_name,max_phase,natural_product,"
              "molecule_type,molecule_properties,molecule_structures,molecule_hierarchy")
    ch = _c(cache)
    out: dict[str, dict] = {}
    todo: list[str] = []
    n_cached = 0
    for cid in ids:
        hit = ch.get(ch.key("molecule", {"id": cid, "fields": fields})) if ch else None
        if hit is not None and not hit.stale and hit.value is not None:
            out[cid] = hit.value
            n_cached += 1
        else:
            todo.append(cid)

    failed: list[str] = []
    for i in range(0, len(todo), batch):
        chunk = todo[i:i + batch]
        try:
            data = _get(f"{CHEMBL}/molecule.json",
                        {"molecule_chembl_id__in": ",".join(chunk),
                         "limit": batch, "only": fields})
        except ChemblError as e:
            log.error("molecule batch at %d failed: %s", i, e)
            failed.extend(chunk)
            continue
        for m in data.get("molecules", []):
            p = m.get("molecule_properties") or {}
            st = m.get("molecule_structures") or {}
            rec = {
                "name": m.get("pref_name"),
                "max_phase": m.get("max_phase"),
                # ChEMBL answers 0 or 1 here, and either is a fact. A field
                # that is absent is not: `bool()` over it reported every
                # unannotated molecule as "not a natural product" (N-20, the
                # same coercion that swallowed the local fallback's None).
                "natural_product": (None if m.get("natural_product") is None
                                    else bool(m.get("natural_product"))),
                "molecule_type": m.get("molecule_type"),
                "molecule_hierarchy": dict(m.get("molecule_hierarchy") or {}),
                "smiles": st.get("canonical_smiles"),
                "mw": p.get("full_mwt"), "logp": p.get("alogp"),
                "tpsa": p.get("psa"), "hbd": p.get("hbd"), "hba": p.get("hba"),
                "ro5_violations": p.get("num_ro5_violations"),
                "aromatic_rings": p.get("aromatic_rings"),
                "heavy_atoms": p.get("heavy_atoms"),
            }
            cid = m["molecule_chembl_id"]
            out[cid] = rec
            if ch:
                try:
                    ch.put(ch.key("molecule", {"id": cid, "fields": fields}),
                           rec, kind="molecule")
                except OSError as exc:
                    log.warning("molecule %s fetched successfully but could not be "
                                "cached (%s); continuing without caching it", cid, exc)

    # A batch that failed is an outage, not an answer. Serve a stale entry for
    # those ids rather than dropping the compound, and say that is what happened.
    stale_ids: list[str] = []
    too_old: list[str] = []
    for cid in failed:
        hit = ch.get(ch.key("molecule", {"id": cid, "fields": fields})) if ch else None
        legacy_fields = fields.removesuffix(",molecule_hierarchy")
        legacy_hit = False
        if hit is None and ch:
            # Existing exact structures remain useful during a service outage.
            # Do not overwrite/re-date the prior cache, and do not invent the
            # newly requested parent hierarchy when that old record lacks it.
            hit = ch.get(ch.key("molecule", {"id": cid, "fields": legacy_fields}))
            legacy_hit = hit is not None
        if hit is None or hit.value is None:
            continue
        # The same ceiling `cache.cached` enforces, which this path bypassed by
        # reading the store directly. ChEMBL releases roughly quarterly, so past
        # `MAX_STALE_DAYS` the entry is not "possibly out of date", it is a
        # different database -- and the reason the ceiling exists is that an
        # outage must not be answered with one. Measured before this check: a
        # 210-day-old entry was served against a 180-day limit and counted as a
        # merely stale one.
        if (hit.age_days or 0.0) > cache_mod.MAX_STALE_DAYS:
            too_old.append(cid)
            continue
        out[cid] = (dict(hit.value, molecule_hierarchy={},
                         hierarchy_status="NOT_EVALUATED_LEGACY_CACHE")
                    if legacy_hit else hit.value)
        stale_ids.append(cid)
    if stale_ids:
        log.warning("%d molecules served from stale cache entries after a failed "
                    "fetch; structures may be out of date", len(stale_ids))
    if too_old:
        log.error("%d molecules had only cache entries older than the %.0f-day "
                  "serving limit after a failed fetch; they are reported missing "
                  "rather than answered from a different ChEMBL release",
                  len(too_old), cache_mod.MAX_STALE_DAYS)

    # Last resort: the bulk chemreps index, if one has been built. It carries
    # structures only -- no name, phase or property block -- so a compound
    # recovered here is scorable (RDKit recomputes every descriptor the ranking
    # uses from SMILES) but under-annotated, and the ids are reported so a
    # screen can say which rows came from a local file instead of the API.
    local_ids: list[str] = []
    still_missing = [i for i in ids if i not in out]
    if still_missing:
        try:
            found = chembl_local.lookup(still_missing)
        except Exception as e:                       # never let a fallback throw
            log.warning("local chemreps lookup failed, skipping it: %s", e)
            found = {}
        if found:
            out.update(chembl_local.as_molecule_records(found))
            local_ids = sorted(found)
            log.warning("%d molecules recovered from the local chemreps index; "
                        "these carry structures only, with no name or phase",
                        len(local_ids))

    missing = [i for i in ids if i not in out]
    if missing:
        log.warning("%d of %d molecules were not retrieved (%d in failed batches)",
                    len(missing), len(ids), len(failed))
    if not with_meta:
        return out
    return out, {"requested": len(ids), "retrieved": len(out),
                 "missing_ids": missing, "failed_batch_ids": failed,
                 "from_cache": n_cached, "stale_ids": stale_ids,
                 # Ids whose only cache entry was past the serving limit. Kept
                 # apart from `missing_ids`: nobody answered for them *and* the
                 # fallback refused, which is a different thing to report.
                 "too_stale_ids": too_old,
                 "stale": bool(stale_ids), "complete": not missing,
                 "local_ids": local_ids, "from_local": len(local_ids)}


def pdb_count(accession: str, cache=None) -> int:
    """Number of RCSB entries mapped to a UniProt accession.

    Raises on a transport failure rather than returning 0. A silent 0 here is
    indistinguishable from "this protein has no structures", and the caller
    grades a target on exactly that difference.
    """
    return cached(_c(cache), "pdb_count", {"accession": accession},
                  lambda: _pdb_count_live(accession)).value


def _pdb_count_live(accession: str) -> int:
    q = {"query": {"type": "terminal", "service": "text", "parameters": {
            "attribute": ("rcsb_polymer_entity_container_identifiers"
                          ".reference_sequence_identifiers.database_accession"),
            "operator": "exact_match", "value": accession}},
         "return_type": "entry", "request_options": {"return_counts": True}}





    from ..structure.rcsb import RCSBError, _post_json
    try:
        data = _post_json("https://search.rcsb.org/rcsbsearch/v2/query", q)
        return int(data.get("total_count", 0))
    except (RCSBError, ValueError) as e:
        raise ChemblError(f"RCSB search for {accession} failed: {e}") from e
