"""Lossless measurements and a deterministic, condition-specific summary.

A compound identity does not erase assays. Only exact concentration endpoints
within the same recorded biological/experimental context share a median. Missing
conditions remain explicitly unknown; this summary is not a pooled meta-analysis.
"""
from __future__ import annotations

import copy
import hashlib
import json
import math
import statistics
from collections import Counter, defaultdict

POLICY_VERSION = "measurements-v4-human-scoped-range"
_UNEVALUATED = {"not evaluated", "not determined", "not tested", "not assayed",
                "not measured", "not available"}
_QUALITATIVE = {"active", "inactive", "not active"}
_MOLAR = {"M": 1.0, "mM": 1e-3, "uM": 1e-6, "µM": 1e-6,
          "μM": 1e-6, "nM": 1e-9, "pM": 1e-12}
_ENDPOINT_ORDER = {"Ki": 0, "Kd": 1, "IC50": 2, "EC50": 3, "AC50": 4,
                   "XC50": 5, "Potency": 6}
_CONTEXT = ("target_chembl_id", "target_organism", "target_species",
            "assay_type", "assay_cell_type", "assay_tissue", "assay_organism",
            "assay_host_species", "assay_description", "assay_conditions",
            "assay_strain", "assay_tax_id", "assay_target_confidence")


def _finite(value):
    if isinstance(value, bool):
        return None
    try:
        val = float(value)
    except (ValueError, TypeError, OverflowError):
        return None
    return val if math.isfinite(val) else None


def _canonical(row):
    return json.dumps(row, sort_keys=True, ensure_ascii=False, default=str,
                      separators=(",", ":"))


def unique_measurements(rows):
    """Remove exact duplicate rows only; preserve conflicting repeat records."""
    unique = {_canonical(row): copy.deepcopy(row) for row in rows}
    return [unique[key] for key in sorted(unique)]


def observation_status(row):
    """An activity record is not necessarily a completed measurement.

    Explicit unavailable/unevaluated comments override numeric placeholders.
    Censored numeric observations remain measurements, but not exact potencies.
    """
    raw = row.get("raw_activity") or {}
    comment = str(row.get("activity_comment") or "").strip().casefold()
    if not comment:
        comment = str(raw.get("activity_comment") or "").strip().casefold()
    if comment in _UNEVALUATED:
        return "NOT_EVALUATED"
    if any(_finite(value) is not None for value in
           (row.get("standard_value"), row.get("pchembl"),
            raw.get("standard_value"), raw.get("pchembl_value"))):
        return "MEASURED"
    if comment in _QUALITATIVE:
        return "MEASURED_QUALITATIVE"
    return "ACTIVITY_RECORD_ONLY"


def activity_status(activity):
    """Aggregate observation state without mistaking a missing median for no data."""
    rows = activity.get("measurements") or []
    if rows:
        states = [observation_status(row) for row in rows]
        if "MEASURED" in states:
            return "MEASURED"
        if "MEASURED_QUALITATIVE" in states:
            return "MEASURED_QUALITATIVE"
        return "NOT_EVALUATED" if all(s == "NOT_EVALUATED" for s in states) else "ACTIVITY_RECORD_ONLY"
    if activity.get("curated_claim_status") == "UNREVIEWED_TARGET_ENDPOINT":
        return "CURATED_CLAIM_ONLY"
    return "ACTIVITY_VALUE_RECORDED" if _finite(activity.get("pchembl")) is not None else "NOT_MEASURED"


def _potency_quality_excluded(row):
    # Retain the observation, but do not promote a database quality warning to
    # an exact ranking potency. A manual validation explicitly clears it.
    validity = str(row.get("data_validity_comment") or "").strip()
    if not validity:
        validity = str((row.get("raw_activity") or {}).get("data_validity_comment") or "").strip()
    return validity.casefold() not in ("", "manually validated")


def _pvalue(row):
    if observation_status(row) != "MEASURED" or _potency_quality_excluded(row):
        return None
    relation = row.get("standard_relation", "=")
    if relation != "=":
        return None
    endpoint = row.get("standard_type") or row.get("type")
    if endpoint not in _ENDPOINT_ORDER:
        return None
    value = _finite(row.get("standard_value"))
    unit = row.get("standard_units")
    if value is not None and value > 0 and unit in _MOLAR:
        return -math.log10(value) - math.log10(_MOLAR[unit])
    # pChEMBL is already a molar logarithm, retained when raw units are absent.
    return _finite(row.get("pchembl"))


def summarise(rows, *, require_human_target=False):
    measurements = unique_measurements(rows)
    bins = defaultdict(list)
    for row in measurements:
        key = {"endpoint": row.get("standard_type") or row.get("type"),
               "human_target_eligible": bool(row.get(
                   "human_target_activity_eligible", row.get("provenance_trustworthy"))),
               "provenance": row.get("provenance"),
               **{field: row.get(field) for field in _CONTEXT}}
        bins[_canonical(key)].append(row)
    groups = []
    for key, members in sorted(bins.items()):
        # Republished identical assay/document measurements count once in the
        # median; all original rows are still retained above.
        values = {}
        for row in members:
            pv = _pvalue(row)
            if pv is not None:
                ident = (row.get("assay_chembl_id"), row.get("document_chembl_id"),
                         row.get("standard_value"), row.get("standard_units"), pv)
                values[ident] = pv
        pvals = sorted(values.values())
        group = {"group_id": hashlib.sha256(key.encode()).hexdigest()[:16],
                 "context": json.loads(key), "n_rows": len(members),
                 "n_exact_measurements": len(pvals),
                 "pchembl_median": statistics.median(pvals) if pvals else None,
                 "pchembl_range": [min(pvals), max(pvals)] if pvals else None,
                 "assay_ids": sorted({str(r["assay_chembl_id"]) for r in members
                                      if r.get("assay_chembl_id")}),
                 "document_ids": sorted({str(r["document_chembl_id"]) for r in members
                                         if r.get("document_chembl_id")})}
        groups.append(group)
    eligible = [g for g in groups if g["pchembl_median"] is not None
                and (not require_human_target or g["context"]["human_target_eligible"])]
    def priority(group):
        ctx = group["context"]
        return (not ctx["human_target_eligible"], ctx.get("assay_type") != "B",
                _ENDPOINT_ORDER.get(ctx["endpoint"], 99),
                -len(group["document_ids"]), -group["n_exact_measurements"])
    preferred = [g for g in eligible if priority(g) == min(map(priority, eligible))] if eligible else []
    # Equal-quality but different contexts do not have a scientifically
    # justified winning condition. An identifier must not choose the potency.
    primary = preferred[0] if len(preferred) == 1 else None
    all_values = [v for r in measurements if (v := _pvalue(r)) is not None]
    scoped_values = [v for r in measurements if (v := _pvalue(r)) is not None
                     and (not require_human_target or bool(r.get(
                         "human_target_activity_eligible", r.get("provenance_trustworthy"))))]
    return {"measurement_policy": POLICY_VERSION, "measurements": measurements,
            "exact_potency_scope": "human_target_only" if require_human_target else "all_recorded_contexts",
            "activity_status": activity_status({"measurements": measurements}),
            "observation_status_counts": dict(sorted(Counter(
                observation_status(row) for row in measurements).items())),
            "measurement_groups": groups, "n_activities": len(measurements),
            "primary_measurement_group": primary["group_id"] if primary else None,
            "primary_measurement_status": "EVALUATED" if primary else "AMBIGUOUS_CONTEXTS" if preferred else "NOT_EVALUATED",
            "primary_measurement_limitation": (
                "Exact measurements exist in distinct contexts; no pooled potency was selected. Inspect measurement_groups."
                if len(preferred) > 1 else
                "No eligible exact potency; measured bounds, qualitative observations and excluded quality records remain preserved."
                if not preferred else None),
            "n_exact_measurements": sum(g["n_exact_measurements"] for g in groups),
            "n_potency_quality_excluded": sum(
                observation_status(r) == "MEASURED" and _potency_quality_excluded(r)
                for r in measurements),
            "pchembl": primary["pchembl_median"] if primary else None,
            "pchembl_range": [min(scoped_values), max(scoped_values)] if scoped_values else None,
            "pchembl_range_scope": "human_target_only" if require_human_target else "all_recorded_contexts",
            "pchembl_range_all_contexts": [min(all_values), max(all_values)] if all_values else None,
            "potency_aggregation": "condition-specific exact-endpoint median",
            "conditions_comparability": "recorded context only; unrecorded conditions unknown"}


def representative(summary):
    gid = summary["primary_measurement_group"]
    if gid is None:
        order = {"MEASURED": 0, "MEASURED_QUALITATIVE": 1,
                 "ACTIVITY_RECORD_ONLY": 2, "NOT_EVALUATED": 3}
        rows = summary["measurements"]
        if summary.get("exact_potency_scope") == "human_target_only":
            rows = [r for r in rows if bool(r.get(
                "human_target_activity_eligible", r.get("provenance_trustworthy")))] or rows
        return min(rows, key=lambda r: order[observation_status(r)]) if rows else {}
    group = next(g for g in summary["measurement_groups"] if g["group_id"] == gid)
    ctx = group["context"]
    return next(row for row in summary["measurements"]
                if _pvalue(row) is not None
                and bool(row.get("human_target_activity_eligible", row.get("provenance_trustworthy"))) == ctx["human_target_eligible"]
                and (row.get("standard_type") or row.get("type")) == ctx["endpoint"]
                and row.get("provenance") == ctx["provenance"]
                and all(row.get(k) == ctx[k] for k in _CONTEXT))


def merge_activity(keeper, other):
    """Keep curated claims separate while combining every measured observation."""
    claims = keeper.get("curated_claims", []) + other.get("curated_claims", [])
    for activity in (keeper, other):
        if not activity.get("curated_claims") and (activity.get("evidence") or activity.get("potency_note")):
            claims.append({k: copy.deepcopy(v) for k, v in activity.items()
                           if k not in ("measurements", "measurement_groups", "curated_claims")})
    rows = keeper.get("measurements", []) + other.get("measurements", [])
    if rows:
        summary = summarise(rows, require_human_target=any(
            activity.get("exact_potency_scope") == "human_target_only"
            for activity in (keeper, other)))
        keeper.update(summary)
        if summary['primary_measurement_group']:
            primary = representative(summary)
            from .fetch.provenance import CONTEXT_FIELDS
            for field in (*_CONTEXT, *CONTEXT_FIELDS, 'provenance', 'provenance_basis',
                          'provenance_trustworthy', 'assay_chembl_id', 'document_chembl_id',
                          'source_url', 'fetched_at'):
                if field in primary:
                    keeper[field] = copy.deepcopy(primary[field])
    if claims:
        keeper["curated_claims"] = unique_measurements(claims)
    return keeper
