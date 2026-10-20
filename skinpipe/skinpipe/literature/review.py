"""Deterministic quote grounding of untrusted literature-review proposals."""
from __future__ import annotations

import copy
import hashlib
import json
import re

MAX_CLAIMS = 20
MAX_QUOTED_WORDS_PER_SOURCE = 25
MAX_REVIEW_BYTES = 128_000


def canonical_bytes(value):
    return json.dumps(value, sort_keys=True, ensure_ascii=False, allow_nan=False,
                      separators=(",", ":")).encode("utf-8")


def validate_review(proposal: object, packet: dict) -> dict:
    """Check source/quote identity; semantic support always needs human review.

    Invalid output is quarantined as a whole. Exact quotations can still be
    cherry-picked, negated, animal-only, abstract-only or irrelevant to the
    proposed target. A completed grounding check is not a scientific PASS.
    """
    errors = []
    checked = []
    ignored = []
    try:
        encoded = canonical_bytes(proposal)
    except (TypeError, ValueError, OverflowError):
        encoded = b""
        errors.append({"code": "NON_JSON_REVIEW"})
    if len(encoded) > MAX_REVIEW_BYTES:
        errors.append({"code": "REVIEW_SIZE_LIMIT"})
    if not isinstance(proposal, dict) or not isinstance(proposal.get("claims"), list):
        errors.append({"code": "INVALID_REVIEW_SCHEMA"})
        claims = []
    else:
        claims = proposal["claims"]
        if len(claims) > MAX_CLAIMS:
            errors.append({"code": "CLAIM_COUNT_LIMIT"})
            claims = []
        ignored = [key for key in ("status", "scientific_status", "experimental_validation",
                                   "human_review_required", "checks") if key in proposal]
    sources = {source["source_id"]: source for source in packet.get("sources", [])}
    used_words = {}
    for index, claim in enumerate(claims):
        reason = None
        if not isinstance(claim, dict):
            errors.append({"claim_index": index, "code": "INVALID_CLAIM_SCHEMA"})
            continue
        source_id = claim.get("source_id")
        source = sources.get(source_id) if isinstance(source_id, str) else None
        text = source.get("text") if source else None
        quote = claim.get("quote")
        statement = claim.get("claim")
        if source is None:
            reason = "UNKNOWN_SOURCE_ID"
        elif not isinstance(text, str) or hashlib.sha256(text.encode("utf-8")).hexdigest() != source.get("text_sha256"):
            reason = "SOURCE_TEXT_HASH_MISMATCH"
        elif not isinstance(statement, str) or not statement.strip() or len(statement) > 1500:
            reason = "INVALID_CLAIM_TEXT"
        elif not isinstance(quote, str) or not quote.strip() or len(quote) > 500 or quote not in text:
            reason = "QUOTE_NOT_EXACT_SOURCE_SUBSTRING"
        if reason:
            errors.append({"claim_index": index, "code": reason})
            continue
        start = text.find(quote)
        end = start + len(quote)
        if "quote_start" in claim or "quote_end" in claim:
            supplied_start, supplied_end = claim.get("quote_start"), claim.get("quote_end")
            if (type(supplied_start) is not int or type(supplied_end) is not int
                    or supplied_start < 0 or supplied_end > len(text)
                    or text[supplied_start:supplied_end] != quote):
                errors.append({"claim_index": index, "code": "QUOTE_OFFSET_MISMATCH"})
                continue
            start, end = supplied_start, supplied_end
        # Count CJK characters conservatively as separate words rather than
        # allowing one unspaced paragraph to evade the short-quotation budget.
        quote_words = len(re.findall(r"[\u3400-\u4dbf\u4e00-\u9fff]|[^\s\u3400-\u4dbf\u4e00-\u9fff]+", quote))
        used_words[source_id] = used_words.get(source_id, 0) + quote_words
        if used_words[source_id] > MAX_QUOTED_WORDS_PER_SOURCE:
            errors.append({"claim_index": index, "code": "PER_SOURCE_QUOTE_BUDGET_EXCEEDED"})
            continue
        normalized = {"source_id": source_id, "claim": statement.strip(), "quote": quote,
                      "quote_start": start, "quote_end": end,
                      "status": "HUMAN_REVIEW_REQUIRED", "scientific_status": "NOT_ESTABLISHED",
                      "grounding": "EXACT_QUOTE_PRESENT",
                      "source_coverage": source.get("coverage"),
                      "source_retracted": source.get("retracted"),
                      "source_article_class": source.get("article_class"),
                      "source_model_scope": copy.deepcopy(source.get("model_scope"))}
        for field in ("organism", "evidence_type", "relation", "reasoning", "limitations"):
            value = claim.get(field)
            if value is not None and (not isinstance(value, str) or len(value) > 3000):
                errors.append({"claim_index": index, "code": "INVALID_CLAIM_FIELD", "field": field})
                reason = "INVALID_CLAIM_FIELD"
                break
            normalized[field] = value.strip() if isinstance(value, str) and value.strip() else "UNKNOWN; human review required"
        if not reason:
            checked.append(normalized)
    summary = proposal.get("summary", "") if isinstance(proposal, dict) else ""
    if not isinstance(summary, str) or len(summary) > 6000:
        errors.append({"code": "INVALID_SUMMARY"})
        summary = ""
    return {"status": "HOLD_INVALID_REVIEW" if errors else "REVIEW_PROPOSED_HUMAN_REVIEW_REQUIRED",
            "claims": [] if errors else checked, "errors": errors,
            "summary": "" if errors else summary,
            "grounding_check": "HOLD" if errors else "COMPLETE",
            "human_review_required": True, "experimental_validation": "NOT_EVALUATED",
            "scientific_status": "NOT_ESTABLISHED", "model_status_fields_ignored": ignored,
            "quote_word_counts": used_words,
            "interpretation": "Exact quote presence is checked, not scientific entailment, organism identity, biological direction, efficacy or completeness of the literature."}
