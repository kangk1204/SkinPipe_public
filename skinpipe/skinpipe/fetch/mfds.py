"""식약처 화장품 원료성분정보: the registered cosmetic raw materials, with CAS.

What this is, and what it is not
--------------------------------
This is an **inventory**, not a permission. A substance appearing here means
the regulator recognises it as a cosmetic raw material and has recorded a
standard Korean name, an English name and a CAS number for it. It does not mean
the substance may be used, at any concentration, in any product -- 별표 1 still
prohibits what it prohibits, and a 고시원료 permission is a separate instrument
(`regulatory.find_notified`). Reading inventory membership as a clearance is the
mistake this module's callers must not make, so `InventoryHit` carries no
verdict and `eligibility` records it as a note.

Why it is worth fetching at all: it is the only free, officially licensed
source that maps a **Korean cosmetic ingredient name to a CAS number**. The
pipeline screens ChEMBL compounds, which arrive with neither, and 별표 1 carries
a CAS for 619 of its 626 rows against a structure for 15. PubChem synonyms
(`sourcing.pubchem_cas`) close the structure->CAS half; this closes the
name->CAS half and adds the positive signal "this chemotype is already used in
cosmetics", which is exactly what a discovery screen wants to know.

Access
------
`https://apis.data.go.kr/1471000/CsmtcsIngdCpntInfoService01`, operation
`getCsmtcsIngdCpntInfoService01`. Licence: 이용허락범위 제한 없음 (no usage
restriction), free. It needs a `serviceKey`, issued per account at
data.go.kr after 활용신청 on dataset 15111774 -- so the key cannot be shipped
and is read from `MFDS_SERVICE_KEY` (or passed in).

Response fields, from the published Swagger:

    INGR_KOR_NAME          표준명
    INGR_ENG_NAME          영문명
    CAS_NO                 CASNo
    ORIGIN_MAJOR_KOR_NAME  기원및정의
    INGR_SYNONYM           이명
"""
from __future__ import annotations

import logging
import os
import time
import urllib.parse
from typing import Any, Iterator, Optional

import requests

log = logging.getLogger(__name__)

HOST = "https://apis.data.go.kr/1471000"

#: The MFDS cosmetics services this module knows how to page, keyed by the name
#: used on the command line. Each is `(service, operation, data.go.kr pk)`.
#: They are listed rather than passed as free strings so a typo fails here
#: instead of arriving as an empty dataset that looks like a real zero.
SERVICES: dict[str, tuple[str, str, str]] = {
    # 화장품 원료성분정보: the inventory of registered raw materials.
    "ingredients": ("CsmtcsIngdCpntInfoService01",
                    "getCsmtcsIngdCpntInfoService01", "15111774"),
    # 화장품 사용제한 원료정보: prohibited and limited substances, per country,
    # with the proviso and limit text. This is the authoritative form of what
    # `config/kr_prohibited_ingredients.yaml` currently holds as a PDF parse --
    # the one whose own notes record 47 of 626 Korean names arriving truncated.
    "restrictions": ("CsmtcsUseRstrcInfoService",
                     "getCsmtcsUseRstrcInfoService", "15111772"),
    # The per-country detail behind the same list.
    "restrictions_by_country": ("CsmtcsUseRstrcInfoService",
                                "getCsmtcsUseRstrcNatnInfoService", "15111772"),
    # 화장품 규제정보: standard name, English name, and which countries prohibit
    # or restrict it.
    "regulations": ("CsmtcsReglMaterialInfoService",
                    "getCsmtcsReglMaterialInfoService", "15111773"),
    # 기능성화장품 보고품목정보: reported functional-cosmetic products. Product
    # level, not ingredient level -- market context rather than a gate input.
    "functional_products": ("FtnltCosmRptPrdlstInfoService01",
                            "getRptPrdlstInq01", "15095680"),
    "functional_firms": ("FtnltCosmRptPrdlstInfoService01",
                         "getRptEntpInq01", "15095680"),
}

BASE = f"{HOST}/{SERVICES['ingredients'][0]}/{SERVICES['ingredients'][1]}"

#: data.go.kr caps a development key at 10,000 calls a day, so pages are taken
#: large. 100 is the portal's usual ceiling for `numOfRows`.
PAGE_SIZE = 100

#: The gateway answers 401 with `SERVICE_KEY_IS_NULL` / `SERVICE_KEY_IS_NOT_
#: REGISTERED_ERROR` rather than a transport error, so those are read from the
#: body, not the status line.
_KEY_ERRORS = ("SERVICE_KEY_IS_NULL", "SERVICE_KEY_IS_NOT_REGISTERED_ERROR",
               "SERVICE_ACCESS_DENIED_ERROR")

_RETRYABLE = frozenset({408, 429, 500, 502, 503, 504})
_TIMEOUT = 60


class MfdsKeyError(RuntimeError):
    """No usable `serviceKey`. Distinct from an outage: retrying will not help."""


class MfdsError(RuntimeError):
    """The service answered, and the answer was not data."""


def service_key(explicit: Optional[str] = None) -> str:
    """The key, in the one form `requests` can send.

    data.go.kr issues every key twice -- "Encoding" (percent-escaped) and
    "Decoding" (raw) -- and the portal shows them side by side, so either is
    what someone has in hand. Only the raw one can go through `params=`:
    `requests` escapes what it is given, so an already-escaped key arrives as
    `%252B` and the gateway answers SERVICE_KEY_IS_NOT_REGISTERED_ERROR, which
    reads as a bad key rather than a double-encoded one. Unquoting an
    already-raw key is a no-op (`+`, `/` and `=` carry no percent escapes), so
    normalising here is safe in both directions.
    """
    key = (explicit or os.environ.get("MFDS_SERVICE_KEY", "")).strip()
    if not key:
        raise MfdsKeyError(
            "no MFDS service key. Register at data.go.kr, 활용신청 on dataset "
            "15111774 (식품의약품안전처_화장품 원료성분정보), then export the "
            "issued key as MFDS_SERVICE_KEY. The key is per account and cannot "
            "be shipped with this package.")
    return urllib.parse.unquote(key)


def _page(session: requests.Session, key: str, page_no: int,
          num_rows: int, retries: int, url: str = BASE) -> dict[str, Any]:
    params = {"serviceKey": key, "pageNo": page_no,
              "numOfRows": num_rows, "type": "json"}
    for attempt in range(retries):
        try:
            r = session.get(url, params=params, timeout=_TIMEOUT)
        except requests.RequestException as exc:
            log.debug("MFDS page %d attempt %d failed: %s", page_no, attempt + 1, exc)
            time.sleep(min(2 ** attempt, 8))
            continue
        # A key problem is terminal and an outage is not; the gateway reports
        # both in the body, so the body is what decides.
        body_text = r.text or ""
        for marker in _KEY_ERRORS:
            if marker in body_text:
                raise MfdsKeyError(
                    f"the service rejected the key ({marker}). Check the key is "
                    "the decoded one and that 활용신청 for dataset 15111774 has "
                    "been approved.")
        if r.status_code in _RETRYABLE:
            log.debug("MFDS page %d -> HTTP %s, retrying", page_no, r.status_code)
            time.sleep(min(2 ** attempt, 8))
            continue
        if r.status_code != 200:
            raise MfdsError(f"MFDS page {page_no} -> HTTP {r.status_code}: "
                            f"{body_text[:200]}")
        try:
            return r.json()
        except ValueError:
            # XML fault dressed as a 200. Not retryable on its own, but one
            # more attempt costs little and the gateway does this transiently.
            log.debug("MFDS page %d returned a non-JSON body", page_no)
            time.sleep(min(2 ** attempt, 8))
            continue
    raise MfdsError(f"MFDS page {page_no} never returned a usable body after "
                    f"{retries} attempts")


def _rows(payload: dict[str, Any]) -> list[dict[str, Any]]:
    """The item list, whatever shape the gateway chose this time.

    `items` is a list on a multi-row page, a bare object on a single-row one,
    and `{"item": [...]}` when the gateway keeps the XML nesting. All three
    have been observed across data.go.kr services; a reader that assumes one
    of them silently returns nothing for the other two.
    """
    body = (payload.get("body") or
            (payload.get("response") or {}).get("body") or {})
    items = body.get("items")
    if items is None:
        return []
    if isinstance(items, dict):
        inner = items.get("item", items)
        return [inner] if isinstance(inner, dict) else list(inner or [])
    if isinstance(items, list):
        out: list[dict[str, Any]] = []
        for entry in items:
            if isinstance(entry, dict) and set(entry) == {"item"}:
                inner = entry["item"]
                out += inner if isinstance(inner, list) else [inner]
            elif isinstance(entry, dict):
                out.append(entry)
        return out
    return []


def _total(payload: dict[str, Any]) -> Optional[int]:
    body = (payload.get("body") or
            (payload.get("response") or {}).get("body") or {})
    try:
        return int(body.get("totalCount"))
    except (TypeError, ValueError):
        return None


def iter_rows(dataset: str = "ingredients", key: Optional[str] = None, *,
                     session: Optional[requests.Session] = None,
                     page_size: int = PAGE_SIZE, retries: int = 3,
                     max_pages: Optional[int] = None,
                     start_page: int = 1,
                     pause: float = 0.0) -> Iterator[dict[str, Any]]:
    """Yield every row of one MFDS dataset, a page at a time.

    `dataset` is a key of `SERVICES`.

    Termination is on the page that comes back short or empty, cross-checked
    against `totalCount`: a loop that trusts `totalCount` alone runs forever
    when the gateway reports a count it will not serve, and one that trusts
    only the short page stops early when a page is dropped. Both conditions
    are logged when they disagree, because a silent early stop here is a
    truncated inventory that looks complete.
    """
    if dataset not in SERVICES:
        raise MfdsError(f"unknown dataset {dataset!r}; "
                        f"known: {', '.join(sorted(SERVICES))}")
    service, operation, _pk = SERVICES[dataset]
    url = f"{HOST}/{service}/{operation}"
    key = service_key(key)
    session = session or requests.Session()
    seen = 0
    total: Optional[int] = None
    page_no = max(1, start_page)
    last_page = None if max_pages is None else page_no + max_pages - 1
    while last_page is None or page_no <= last_page:
        payload = _page(session, key, page_no, page_size, retries, url)
        if total is None:
            total = _total(payload)
        rows = _rows(payload)
        if not rows:
            break
        for row in rows:
            yield row
        seen += len(rows)
        if len(rows) < page_size:
            break
        if total is not None and seen >= total:
            break
        page_no += 1
        if pause:
            time.sleep(pause)
    # `seen` counts what *this* call yielded, so a resumed run has to add back
    # the rows the earlier run already took, or every resume would report the
    # copy as short and the warning would stop meaning anything.
    already = (max(1, start_page) - 1) * page_size
    if total is not None and max_pages is None and seen + already < total:
        # Collectors only publish complete=true after this iterator finishes.
        # Raise after yielding the retained rows so their normal checkpoint
        # path records complete=false rather than blessing a truncated copy.
        raise MfdsError(f"MFDS {dataset} stopped at {seen + already} rows but "
                        f"the service reports {total}; the copy is incomplete")


def iter_ingredients(key: Optional[str] = None, **kw) -> Iterator[dict[str, Any]]:
    """The raw-material inventory. `iter_rows("ingredients", ...)` by another name."""
    return iter_rows("ingredients", key, **kw)
