"""Bounded Europe PMC receipts; source retrieval never certifies a target claim.

Official API: https://europepmc.org/RestfulWebService
`search` with resultType=core includes abstracts and indexed metadata;
`/{PMCID}/fullTextXML` is for the open-access subset. This implementation
fetches one bounded relevance page, not an exhaustive or negative search.
"""
from __future__ import annotations

import copy
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re
import uuid
from urllib.parse import urlencode
import xml.etree.ElementTree as ET

import requests

from .review import canonical_bytes, validate_review

SCHEMA = "skinpipe-literature-evidence-v1"
BASE = "https://www.ebi.ac.uk/europepmc/webservices/rest"
MAX_PAPERS = 8
MAX_SEARCH_BYTES = 2_000_000
MAX_XML_BYTES = 4_000_000
MAX_REVIEW_TEXT_BYTES = 8_000
_GENE = re.compile(r"[A-Z][A-Z0-9-]{0,29}\Z")
_QUERY = re.compile(r"(?:[A-Z][A-Z0-9-]{0,29}|CHEMBL[0-9]{1,12})\Z")
_PMID = re.compile(r"[0-9]{1,12}\Z")
_PMCID = re.compile(r"PMC[0-9]{1,12}\Z")
_ACCESSION = re.compile(r"(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-[0-9]+)?\Z")
_ACTIONS = {"inhibit", "activate", "increase", "decrease", "suppress", "unknown"}
_HEADERS = {"User-Agent": "SkinPipe-literature/1.0 (bounded public literature evidence)"}
PUBLIC_COMPOUND_SCHEMA = "skinpipe-public-compound-literature-v1"
_USAGE_LIMIT_CODES = {"usage_exhausted", "usage_limit_exceeded", "insufficient_quota", "usage_limit_reached",
                      "subscription_sharing_usage_limit_exceeded", "billing_hard_limit_reached",
                      "credit_balance_exhausted", "organization_spend_limit_exceeded",
                      "project_spend_limit_exceeded", "organization_usage_limit_exceeded"}
_RATE_LIMIT_CODES = {"rate_limit_exceeded", "rate_limited", "too_many_requests", "slow_down"}
_REVIEW_DIAGNOSTIC_CODES = {
    "network_error", "request_timed_out", "review_timed_out", "interrupted_response_stream",
    "response_not_completed", "streaming_response_required", "selected_model_unavailable",
    "model_selection_required", "inference_http_error", "openai_http_error",
    "invalid_provider_response", "provider_response_too_large", "invalid_model_catalog",
    "review_response_too_large", "invalid_response_stream", "missing_response_identity",
    "inconsistent_response_text", "review_json_required", "review_claims_required",
    "invalid_reasoning_effort", "sign_in_required", "insufficient_scope", "invalid_token",
    "access_denied", "temporarily_unavailable", "subscription_sharing_usage_unavailable",
    "subscription_sharing_unsupported_capability", "subscription_sharing_user_not_eligible",
    "subscription_sharing_route_not_supported", "subscription_sharing_invalid_user",
    "subscription_sharing_user_unavailable", "chatpass_v2_scope_not_authorized",
    "chatpass_v2_invalid_authorization_context",
}


def _hash(raw):
    return hashlib.sha256(raw).hexdigest()


def _now():
    return datetime.now(timezone.utc).isoformat()


def _review_failure(exc):
    """Persist known error attributes only, never exception/body/ID text.

    HTTP 429 alone establishes a call limit, not exhausted account usage.
    Provider timing hints are recorded only; this function never waits/retries.
    """
    def field(name):
        try:
            return getattr(exc, name, None)
        except Exception:
            return None
    code, provider_code = field("code"), field("provider_code")
    usage = (field("usage_exhausted") is True
             or isinstance(code, str) and code in _USAGE_LIMIT_CODES
             or isinstance(provider_code, str) and provider_code in _USAGE_LIMIT_CODES)
    http_status = field("http_status")
    http_status = http_status if type(http_status) is int and 100 <= http_status <= 599 else None
    rate = (isinstance(code, str) and code in _RATE_LIMIT_CODES
            or isinstance(provider_code, str) and provider_code in _RATE_LIMIT_CODES
            or http_status == 429)
    safe = {"code": "usage_exhausted" if usage else "rate_limit_exceeded" if rate else "review_unavailable"}
    if isinstance(code, str) and code in _REVIEW_DIAGNOSTIC_CODES:
        safe["diagnostic_code"] = code
    body_shape = field("body_shape")
    if isinstance(body_shape, str) and body_shape in {"error_object", "error_string", "detail_object", "other_object"}:
        safe["body_shape"] = body_shape
    if isinstance(provider_code, str) and provider_code in _USAGE_LIMIT_CODES | _RATE_LIMIT_CODES:
        safe["provider_code"] = provider_code
    if http_status is not None:
        safe["http_status"] = http_status
    if usage or rate:
        retry_after = field("retry_after_seconds")
        if type(retry_after) in (int, float) and 0 <= retry_after <= 31_536_000:
            safe["retry_after_seconds"] = retry_after
        reset_at = field("reset_at")
        try:
            if isinstance(reset_at, str) and len(reset_at) <= 64:
                reset = datetime.fromisoformat(reset_at.replace("Z", "+00:00"))
            elif type(reset_at) in (int, float) and 0 <= reset_at <= 4_133_980_800:
                reset = datetime.fromtimestamp(reset_at, timezone.utc)
            else:
                reset = None
            if reset is not None and reset.tzinfo is not None and 2000 <= reset.year <= 2100:
                safe["reset_at"] = reset.astimezone(timezone.utc).isoformat().replace("+00:00", "Z")
        except (ValueError, OverflowError, OSError):
            pass
    return {"status": "HOLD_USAGE_LIMIT" if usage else "HOLD_RATE_LIMIT" if rate else "HOLD_REVIEW_UNAVAILABLE",
            "claims": [], "error": safe,
            "errors": [{"code": safe["code"],
                        "reason": "Reviewer call/output unavailable; no scientific conclusion recorded."}],
            "grounding_check": "HOLD", "human_review_required": True,
            "experimental_validation": "NOT_EVALUATED", "scientific_status": "NOT_ESTABLISHED"}


def _strict_json(raw):
    def pairs(items):
        out = {}
        for key, value in items:
            if key in out:
                raise ValueError("duplicate source JSON key")
            out[key] = value
        return out
    def invalid_constant(_):
        raise ValueError("nonfinite source JSON")
    return json.loads(raw, object_pairs_hook=pairs, parse_constant=invalid_constant)


def public_compound_name(value):
    if (not isinstance(value, str) or not 1 <= len(value.strip()) <= 512
            or any(ord(char) < 32 for char in value) or re.match(r"[A-Za-z][A-Za-z0-9+.-]*://", value.strip())):
        raise ValueError("Public literature requires a bounded individual compound name")
    return value.strip()


def public_skin_topic(value=""):
    if not isinstance(value, str):
        raise ValueError("Public skin topic must be text")
    value = value.strip()
    if value and (not re.fullmatch(r"[A-Za-z][A-Za-z0-9 -]{0,127}", value)
                  or len(value.split()) > 12):
        raise ValueError("Public skin topic requires at most 128 English letters, digits, spaces or hyphens and 12 words")
    return value


def _compound_descriptor(raw):
    allowed = {"schema", "kind", "query", "pubchem_cid", "public_topic", "taxon_id", "action",
               "identity_status", "identity_source_sha256", "include_full_text"}
    if (set(raw) - allowed or raw.get("schema") != PUBLIC_COMPOUND_SCHEMA
            or raw.get("identity_status") != "NAME_RESOLVED_BY_PUBCHEM"
            or type(raw.get("pubchem_cid")) is not int or raw["pubchem_cid"] < 1
            or not isinstance(raw.get("identity_source_sha256"), list) or not 1 <= len(raw["identity_source_sha256"]) <= 8
            or any(not isinstance(value, str) or not re.fullmatch(r"[0-9a-f]{64}", value) for value in raw["identity_source_sha256"])):
        raise ValueError("Compound literature requires an officially resolved individual compound")
    name, topic = public_compound_name(raw.get("query")), public_skin_topic(raw.get("public_topic", ""))
    taxon = raw.get("taxon_id")
    if type(taxon) is not int or not 1 <= taxon <= 100_000_000:
        raise ValueError("an explicit positive integer taxon_id is required")
    action = raw.get("action", "unknown")
    action = "unknown" if action == "any" else action
    if not isinstance(action, str) or action.lower() not in _ACTIONS:
        raise ValueError("Unsupported literature action")
    full_text = raw.get("include_full_text", True)
    if type(full_text) is not bool:
        raise ValueError("include_full_text must be a boolean")
    return {"kind": "compound", "query": name, "compound_name": name, "pubchem_cid": raw["pubchem_cid"],
            "public_topic": topic, "taxon_id": taxon, "action": action.lower(), "include_full_text": full_text,
            "descriptor_basis": "PubChem-resolved public compound name and CID; no protein substitution",
            "identity_status": "NAME_RESOLVED_BY_PUBCHEM; literature relevance not established"}


def _search_query(target):
    def literal(value):
        return value.replace("\\", "\\\\").replace('"', '\\"')
    if target.get("kind") == "compound":
        terms = [target["compound_name"], *target["public_topic"].split()]
        subject = " AND ".join('TITLE_ABS:"' + literal(term) + '"' for term in terms)
    else:
        subject = f'TITLE_ABS:"{target["gene"]}"'
    return subject + ' AND SRC:MED NOT (PUB_TYPE:"Review" OR PUB_TYPE:"Editorial" OR PUB_TYPE:"Comment")'


def _descriptor(raw):
    if not isinstance(raw, dict):
        raise ValueError("target descriptor must be an object")
    request = raw.get("request") if isinstance(raw.get("request"), dict) else raw
    kind = raw.get("kind", request.get("kind", "protein"))
    if kind == "compound":
        return _compound_descriptor(raw)
    if raw.get("mode") == "endpoint" or kind == "endpoint":
        raise ValueError("endpoint biomarkers need an explicitly selected molecular target")
    if kind != "protein":
        raise ValueError("a complex/proxy needs an explicitly selected single molecular target")
    gene = raw.get("gene")
    query = request.get("query", raw.get("query", gene))
    component = None
    if raw.get("schema") == "target-descriptor-v2":
        components = raw.get("components")
        if (raw.get("status") != "verified" or not isinstance(components, list)
                or len(components) != 1 or not isinstance(components[0], dict)):
            raise ValueError("an unresolved or multi-component descriptor needs explicit target selection")
        component = components[0]
        names = component.get("gene_names", [])
        component_gene = component.get("gene")
        if component_gene is None and isinstance(names, list) and len(names) == 1:
            component_gene = names[0]
        if component_gene is None:
            raise ValueError("the selected component has no unique recorded gene symbol")
        if gene is not None and str(gene).upper() != str(component_gene).upper():
            raise ValueError("top-level and selected component gene identities conflict")
        gene = component_gene
    elif gene is None:
        # The standalone CLI accepts an explicitly typed gene without claiming
        # identity verification. Accessions/target IDs need resolver metadata;
        # searching their text as a gene would silently choose another namespace.
        if not isinstance(query, str) or _ACCESSION.fullmatch(query.upper()) or query.upper().startswith("CHEMBL"):
            raise ValueError("an accession or target ID needs a resolved single-component gene")
        gene = query
    if not isinstance(gene, str) or not _GENE.fullmatch(gene.strip().upper()):
        raise ValueError("a single gene symbol is required; URLs and search expressions are not accepted")
    gene = gene.strip().upper()
    query = query or gene
    if not isinstance(query, str) or not _QUERY.fullmatch(query.strip().upper()):
        raise ValueError("target query must be one identifier, not a URL or search expression")
    taxon = raw.get("taxon_id", raw.get("taxon", request.get("taxon_id")))
    if type(taxon) is not int or not 1 <= taxon <= 100_000_000:
        raise ValueError("an explicit positive integer taxon_id is required")
    if component is not None and component.get("taxon_id") != taxon:
        raise ValueError("selected component taxon conflicts with requested taxon")
    action = raw.get("action", raw.get("desired_action", raw.get("node_action", request.get("desired_action", "unknown"))))
    action = "unknown" if action in (None, "any") else action
    if not isinstance(action, str) or action.lower() not in _ACTIONS:
        raise ValueError("action must be inhibit, activate, increase, decrease, suppress or unknown")
    include_full_text = raw.get("include_full_text", True)
    if type(include_full_text) is not bool:
        raise ValueError("include_full_text must be a boolean")
    return {"query": query.strip().upper(), "gene": gene, "taxon_id": taxon,
            "action": action.lower(), "include_full_text": include_full_text,
            "descriptor_basis": "selected single-component descriptor" if component else "typed target/gene; identity not independently verified",
            "identity_status": "NOT_VERIFIED_BY_LITERATURE_SEARCH"}


class _Artifacts:
    def __init__(self, directory):
        self.directory = directory
        self.records = []

    def write(self, name, raw):
        path = self.directory / name
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("xb") as handle:
            handle.write(raw)
        path.chmod(0o444)
        record = {"path": name, "sha256": _hash(raw), "size_bytes": len(raw)}
        self.records.append(record)
        return record

    def json(self, name, data):
        return self.write(name, canonical_bytes(data))


def _request(session, url, params, timeout, limit, artifacts, prefix):
    """One nonredirecting bounded request; retain even unsuccessful raw bytes."""
    requested_url = url + ("?" + urlencode(params) if params else "")
    receipt = {"requested_url": requested_url, "method": "GET", "retrieved_at": _now(),
               "timeout_seconds": timeout, "byte_limit": limit, "retries": 0,
               "status": "UNAVAILABLE", "http_status": None}
    raw = bytearray()
    response = None
    try:
        response = session.get(url, params=params or None, timeout=timeout,
                               headers=_HEADERS, allow_redirects=False, stream=True)
        receipt["http_status"] = response.status_code
        receipt["response_url"] = str(getattr(response, "url", "") or requested_url)
        receipt["headers"] = {key: value for key, value in getattr(response, "headers", {}).items()
                              if key.lower() in {"content-type", "etag", "last-modified", "date"}}
        for chunk in response.iter_content(chunk_size=65536):
            if not chunk:
                continue
            remaining = limit - len(raw)
            raw.extend(chunk[:remaining])
            if len(chunk) > remaining:
                receipt.update(status="PARTIAL", error="RESPONSE_BYTE_LIMIT")
                break
        else:
            receipt["status"] = "COMPLETE" if response.status_code == 200 else "UNAVAILABLE"
            if response.status_code != 200:
                receipt["error"] = "HTTP_STATUS_NOT_200"
    except (requests.RequestException, OSError, TimeoutError) as exc:
        receipt.update(status="PARTIAL" if raw else "UNAVAILABLE", error=type(exc).__name__)
    finally:
        if response is not None:
            response.close()
    receipt["body"] = artifacts.write(prefix + ".body", bytes(raw))
    receipt["receipt"] = {"path": prefix + ".receipt.json"}
    receipt_record = artifacts.json(prefix + ".receipt.json", receipt)
    receipt["receipt"] = receipt_record
    return bytes(raw), receipt


def _plain_markup(value):
    if not isinstance(value, str):
        return ""
    try:
        root = ET.fromstring("<text>" + value + "</text>")
        text = " ".join(root.itertext())
    except ET.ParseError:
        # Failed markup is retained literally, not silently discarded.
        text = value
    return re.sub(r"\s+", " ", text).strip()


def _strings(value, key):
    if isinstance(value, dict):
        value = value.get(key, [])
    if isinstance(value, str):
        value = [value]
    return sorted({item.strip() for item in value if isinstance(item, str) and item.strip()}) if isinstance(value, list) else []


def _model_scope(row):
    container = row.get("meshHeadingList")
    headings = container.get("meshHeading", []) if isinstance(container, dict) else []
    if not isinstance(headings, list):
        headings = []
    terms = sorted({str(item.get("descriptorName", "")).strip() for item in headings if isinstance(item, dict)} - {""})
    lowered = {term.casefold() for term in terms}
    animals = bool(lowered & {"animals", "mice", "rats", "rabbits", "dogs", "swine", "murinae"})
    humans = "humans" in lowered
    label = ("HUMAN_AND_ANIMAL_INDEX_TERMS_PRESENT" if humans and animals else
             "ANIMAL_INDEX_TERMS_PRESENT" if animals else
             "HUMAN_INDEX_TERM_PRESENT" if humans else "UNKNOWN")
    return label + "; actual experimental organism and assay require human review", terms


def _article_class(types):
    joined = " ".join(types).casefold()
    if "retract" in joined:
        return "RETRACTED_PUBLICATION_OR_NOTICE"
    if "review" in joined:
        return "REVIEW"
    if any(term in joined for term in ("editorial", "comment", "letter", "news")):
        return "NON_PRIMARY_OR_UNCERTAIN"
    if "preprint" in joined:
        return "PREPRINT"
    return "PRIMARY_RESEARCH_CANDIDATE" if any(term in joined for term in ("journal article", "research-article", "clinical trial")) else "UNKNOWN"


def _xml_text(raw, pmcid, pmid):
    decoded = raw.decode("utf-8")
    if "<!ENTITY" in decoded.upper() or "\x00" in decoded:
        raise ValueError("unsafe or unsupported XML entity declaration")
    root = ET.fromstring(decoded)
    if root.tag.split("}")[-1] != "article":
        raise ValueError("full text is not one JATS article")
    found_pmc, found_pmid = set(), set()
    for element in root.iter():
        if element.tag.split("}")[-1] != "article-id":
            continue
        value = "".join(element.itertext()).strip()
        kind = element.attrib.get("pub-id-type", "").lower()
        if kind in {"pmc", "pmcid", "pmcaid"}:
            found_pmc.add(value if value.startswith("PMC") else "PMC" + value)
        elif kind == "pmid":
            found_pmid.add(value)
    if (found_pmc and found_pmc != {pmcid}) or (found_pmid and pmid and found_pmid != {pmid}):
        raise ValueError("full text article identifiers conflict with search metadata")
    if not ((pmcid in found_pmc) or (pmid and pmid in found_pmid)):
        raise ValueError("full text article identity was not verified")
    bodies = [element for element in root if element.tag.split("}")[-1] == "body"]
    if len(bodies) != 1:
        raise ValueError("full text has no unique article body")
    text = re.sub(r"\s+", " ", " ".join(bodies[0].itertext())).strip()
    if not text:
        raise ValueError("full text body is empty")
    licenses = []
    for element in root.iter():
        if element.tag.split("}")[-1] == "license":
            licenses.append({"type": element.attrib.get("license-type"),
                             "url": element.attrib.get("{http://www.w3.org/1999/xlink}href"),
                             "text": " ".join(element.itertext()).strip()[:4000]})
    return text, licenses, root.attrib.get("article-type", "")


def _excerpt(text):
    raw = text.encode("utf-8")
    excerpt = raw[:MAX_REVIEW_TEXT_BYTES].decode("utf-8", errors="ignore")
    return excerpt, len(raw) > MAX_REVIEW_TEXT_BYTES, _hash(raw), len(raw)


def _source(row, number, target, session, timeout, artifacts):
    pmid = str(row.get("pmid") or row.get("id") or "")
    if row.get("source") != "MED" or not _PMID.fullmatch(pmid):
        raise ValueError("search record is not a verifiable MED publication identifier")
    pmcid = str(row.get("pmcid") or "")
    pmcid = pmcid if _PMCID.fullmatch(pmcid) else None
    doi = row.get("doi")
    doi = doi if isinstance(doi, str) and re.fullmatch(r"10\.[0-9]{4,9}/[^\s]{1,240}", doi) else None
    canonical = canonical_bytes(row)
    source_id = "EPMC_MED_" + pmid + "_" + _hash(canonical)[:12]
    prefix = f"sources/{number:02d}_{source_id}"
    metadata = artifacts.write(prefix + ".metadata.json", canonical)
    types = _strings(row.get("pubTypeList"), "pubType")
    retracted = (True if any("retract" in value.casefold() for value in types) or row.get("isRetracted") in (True, "Y") else
                 False if row.get("isRetracted") in (False, "N") else None)
    model_scope, model_terms = _model_scope(row)
    abstract = _plain_markup(row.get("abstractText"))
    text = abstract
    coverage = "ABSTRACT_ONLY" if abstract else "METADATA_ONLY"
    retrieval = "NOT_REQUESTED"
    full_receipt = None
    licenses = []
    if isinstance(row.get("license"), str):
        licenses.append({"type": "Europe PMC index field", "text": row["license"]})
    if target["include_full_text"] and pmcid and row.get("isOpenAccess") == "Y":
        raw, full_receipt = _request(session, f"{BASE}/{pmcid}/fullTextXML", {}, timeout,
                                     MAX_XML_BYTES, artifacts, prefix + ".fulltext")
        retrieval = full_receipt["status"]
        if retrieval == "COMPLETE":
            try:
                text, xml_licenses, xml_article_type = _xml_text(raw, pmcid, pmid)
                licenses += xml_licenses
                if "retract" in xml_article_type.casefold():
                    retracted = True
                    types = sorted(set(types) | {"Retraction publication or notice (JATS XML)"})
                coverage = "OA_FULL_TEXT"
            except (ValueError, ET.ParseError, UnicodeError):
                retrieval = "HOLD_UNUSABLE_OR_IDENTITY_UNVERIFIED"
    elif target["include_full_text"]:
        retrieval = "UNAVAILABLE_NOT_INDEXED_AS_OA_PMC" if not (pmcid and row.get("isOpenAccess") == "Y") else "NOT_REQUESTED"
    text, truncated, full_text_hash, original_text_bytes = _excerpt(text)
    if truncated:
        coverage = "OA_FULL_TEXT_EXCERPT" if coverage == "OA_FULL_TEXT" else "ABSTRACT_EXCERPT"
    text_file = artifacts.write(prefix + ".text.txt", text.encode("utf-8"))
    record = {"source_id": source_id, "title": _plain_markup(row.get("title"))[:1000],
        "url": f"https://europepmc.org/article/MED/{pmid}", "pmid": pmid, "pmcid": pmcid, "doi": doi,
        "author_string": str(row.get("authorString") or "")[:1000], "year": str(row.get("pubYear") or ""),
        "coverage": coverage, "text_sha256": text_file["sha256"],
        "metadata_sha256": _hash(canonical), "publication_types": types[:30],
        "retracted": retracted,
        "retraction_status": "RETRACTED_OR_NOTICE_FLAGGED" if retracted else "NO_FLAG_IN_RETRIEVED_METADATA" if retracted is False else "UNKNOWN",
        "article_class": _article_class(types), "model_scope": model_scope, "indexed_model_terms": model_terms,
        "full_text_retrieval": retrieval, "review_text_truncated": truncated,
        "extracted_full_text_sha256": full_text_hash, "extracted_full_text_bytes": original_text_bytes,
        "extraction": "Whitespace-normalized abstract markup or JATS body only; references excluded; review text is a bounded prefix.",
        "source_identity_status": "INDEX_IDENTIFIERS_RECORDED; target/organism/assay match not established",
        "open_access_index_flag": row.get("isOpenAccess") if row.get("isOpenAccess") in ("Y", "N") else "UNKNOWN",
        "licenses": licenses, "license_status": "RECORDED_NOT_INDEPENDENTLY_VERIFIED" if licenses else "UNKNOWN",
        "artifacts": {"metadata": metadata, "text": text_file},
        "limitations": "Single-provider indexed metadata; no comprehensive retraction search, semantic support, species identity, direct binding, direction or efficacy validation."}
    if full_receipt:
        record["artifacts"]["fulltext_body"] = full_receipt["body"]
        record["artifacts"]["fulltext_receipt"] = full_receipt["receipt"]
    return record, text


def collect_evidence(descriptor: dict, output_dir: Path, *, reviewer=None,
                     max_papers=8, timeout=20, session=None) -> dict:
    """Create one exclusive evidence run, returning a dossier and path envelope.

    `output_dir` must not exist. Persisted paths are relative. The returned
    `artifact_dir` and `dossier_artifact` are envelope fields outside dossier.json.
    A reviewer is a callable accepting a public-data packet and proposing JSON;
    errors/forgeries are explicit HOLD states and never scientific acceptance.
    """
    if type(max_papers) is not int or not 1 <= max_papers <= MAX_PAPERS:
        raise ValueError(f"max_papers must be an integer in 1..{MAX_PAPERS}")
    if isinstance(timeout, bool) or not isinstance(timeout, (int, float)) or not math.isfinite(timeout) or not 0 < timeout <= 120:
        raise ValueError("timeout must be finite and in (0, 120] seconds")
    directory = Path(output_dir)
    directory.mkdir(parents=True, exist_ok=False)
    artifacts = _Artifacts(directory)
    dossier = {"schema": SCHEMA, "run_id": str(uuid.uuid4()), "created_at": _now(),
        "target": {}, "sources": [], "claims": [], "checks": {},
        "human_review_required": True, "experimental_validation": "NOT_EVALUATED",
        "scientific_status": "NOT_ESTABLISHED", "status": "HOLD_UNAVAILABLE",
        "search_scope": {"provider": "Europe PMC MED index", "max_papers": max_papers,
            "exhaustive": False, "pages_requested": 1, "aliases_expanded": False,
            "interpretation": "Bounded relevance search; missing hits never establish biological absence. Taxon and desired action describe the request, not the species/direction demonstrated by a keyword hit."},
        "human_review_checklist": ["Verify exact protein/isoform, organism, tissue and assay in each primary paper.",
            "Verify direction, concentration, endpoint and controls; distinguish direct binding from downstream association.",
            "Check retraction/correction status with an independent current source and inspect complete primary papers.",
            "Do not promote literature or quote checks to project efficacy, clinical readiness, screening policy or structure admission."],
        "coverage": "UNAVAILABLE"}
    own_session = session is None
    client = session if session is not None else requests.Session()
    try:
        try:
            target = _descriptor(descriptor)
        except ValueError as exc:
            dossier.update(status="HOLD_UNSUPPORTED_TARGET", target={"descriptor_status": "UNSUPPORTED"})
            dossier["checks"]["descriptor"] = {"status": "HOLD", "reason": str(exc)}
            return _finish(dossier, artifacts)
        dossier["target"] = target
        dossier["checks"]["descriptor"] = {"status": "COMPLETE", "scope": "syntax and explicit fields only; molecular/source identity not established"}
        artifacts.json("target.json", target)
        query = _search_query(target)
        params = {"query": query, "resultType": "core", "format": "json", "pageSize": max_papers, "cursorMark": "*"}
        raw, receipt = _request(client, BASE + "/search", params, timeout,
                                MAX_SEARCH_BYTES, artifacts, "search")
        dossier["search_scope"].update(query=query, request_url=receipt["requested_url"])
        if receipt["status"] != "COMPLETE":
            dossier["checks"]["bounded_search"] = {"status": "HOLD", "reason": receipt.get("error", receipt["status"])}
            if receipt["status"] == "PARTIAL":
                dossier.update(status="HOLD_PARTIAL_RETRIEVAL", coverage="PARTIAL_RETRIEVAL")
            return _finish(dossier, artifacts)
        try:
            payload = _strict_json(raw)
            rows = payload["resultList"]["result"]
            total = payload["hitCount"]
            if not isinstance(rows, list) or any(not isinstance(row, dict) for row in rows) or type(total) is not int or total < 0:
                raise ValueError("invalid search record shape or hit count")
        except (ValueError, KeyError, TypeError):
            dossier["checks"]["bounded_search"] = {"status": "HOLD", "reason": "INVALID_SEARCH_RESPONSE"}
            return _finish(dossier, artifacts)
        expected = min(total, max_papers)
        partial = len(rows) != expected
        dossier["checks"]["bounded_search"] = {"status": "PARTIAL" if partial else "COMPLETE",
            "hit_count": total, "retrieved_page_rows": len(rows), "expected_page_rows": expected,
            "scope": "one bounded relevance page only; not complete literature coverage"}
        dossier["search_scope"]["provider_hit_count"] = total
        texts, errors, seen = {}, [], set()
        for index, row in enumerate(rows[:max_papers], start=1):
            try:
                source, text = _source(row, index, target, client, timeout, artifacts)
            except (ValueError, TypeError, KeyError) as exc:
                artifacts.json(f"quarantine/search_record_{index:02d}.json", row)
                errors.append({"index": index, "status": "HOLD", "reason": str(exc)})
                continue
            identity = (source["pmid"], source["pmcid"])
            if identity in seen:
                errors.append({"index": index, "status": "HOLD", "reason": "DUPLICATE_PUBLICATION_ID"})
                continue
            seen.add(identity)
            dossier["sources"].append(source)
            texts[source["source_id"]] = text
        retrieval_partial = any(source["full_text_retrieval"] in {"PARTIAL", "UNAVAILABLE", "HOLD_UNUSABLE_OR_IDENTITY_UNVERIFIED"} for source in dossier["sources"])
        partial = partial or bool(errors) or retrieval_partial
        dossier["checks"]["source_snapshots"] = {"status": "PARTIAL" if errors or retrieval_partial else "COMPLETE", "errors": errors, "n_sources": len(dossier["sources"])}
        dossier["coverage"] = "PARTIAL_RETRIEVAL" if partial else "BOUNDED_SOURCE_PACKET" if dossier["sources"] else "NO_SOURCES_IN_BOUNDED_SEARCH"
        dossier["status"] = "HOLD_PARTIAL_RETRIEVAL" if partial else "COLLECTED_HUMAN_REVIEW_REQUIRED" if dossier["sources"] else "HOLD_NO_SOURCES"
        dossier["checks"]["claim_grounding"] = {"status": "NOT_REQUESTED"}
        if reviewer is not None and any(texts.values()):
            packet = {"schema": "skinpipe-literature-review-packet-v1", "public_literature_only": True,
                      "target": copy.deepcopy(target), "search_scope": copy.deepcopy(dossier["search_scope"]),
                      "sources": [{key: copy.deepcopy(source[key]) for key in (
                          "source_id", "title", "url", "pmid", "doi", "coverage", "text_sha256",
                          "publication_types", "retracted", "retraction_status", "article_class", "model_scope")}
                          | {"text": texts[source["source_id"]]} for source in dossier["sources"]
                          if texts[source["source_id"]]]}
            for packet_source in packet["sources"]:
                packet_source["title"] = packet_source["title"][:500]
                packet_source["publication_types"] = [value[:80] for value in packet_source["publication_types"][:8]]
            dossier["checks"]["review_packet"] = {"status": "COMPLETE", "metadata_bounded": True,
                "text_sources": len(packet["sources"]),
                "metadata_only_sources_omitted": len(dossier["sources"]) - len(packet["sources"])}
            artifacts.json("review/packet.json", packet)
            try:
                if len(canonical_bytes(packet)) > 128_000:
                    raise ValueError("review packet exceeds the 128 KB bound")
                proposal = reviewer(copy.deepcopy(packet))
                proposal_bytes = canonical_bytes(proposal)
                artifacts.write("review/proposal.json", proposal_bytes)
                reviewed = validate_review(proposal, packet)
                if isinstance(proposal, dict) and isinstance(proposal.get("llm_receipt"), dict):
                    reviewed["llm_receipt"] = copy.deepcopy(proposal["llm_receipt"])
                if reviewed["status"] == "HOLD_INVALID_REVIEW":
                    reviewed["quarantined_proposal"] = {"path": "review/proposal.json", "sha256": _hash(proposal_bytes)}
            except Exception as exc:
                reviewed = _review_failure(exc)
            artifacts.json("review/checks.json", reviewed)
            dossier["claims"] = reviewed["claims"]
            dossier["checks"]["claim_grounding"] = {"status": reviewed["grounding_check"], "errors": reviewed["errors"]}
            if "error" in reviewed:
                dossier["checks"]["claim_grounding"]["error"] = copy.deepcopy(reviewed["error"])
            dossier["review"] = reviewed
            if reviewed["status"].startswith("HOLD"):
                dossier["status"] = reviewed["status"]
            elif not partial:
                dossier["status"] = reviewed["status"]
        elif reviewer is not None and dossier["sources"]:
            dossier["status"] = "HOLD_REVIEW_UNAVAILABLE"
            dossier["checks"]["claim_grounding"] = {"status": "HOLD", "reason": "NO_REVIEWABLE_SOURCE_TEXT"}
        return _finish(dossier, artifacts)
    finally:
        if own_session:
            client.close()


def _finish(dossier, artifacts):
    dossier["checks"]["human_scientific_review"] = {"status": "REQUIRED"}
    dossier["checks"]["real_experiment"] = {"status": "NOT_EVALUATED"}
    dossier["artifact_paths"] = copy.deepcopy(artifacts.records)
    record = artifacts.json("dossier.json", dossier)
    return copy.deepcopy(dossier) | {"artifact_dir": str(artifacts.directory.resolve()), "dossier_artifact": record}
