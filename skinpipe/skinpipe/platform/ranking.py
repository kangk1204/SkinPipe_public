"""Cohort-bound, missing-aware priority heuristic; never a binding probability.

Each candidate contributes one median value per compatible cohort. Midrank
percentiles are O(n log n) overall, with O(n + observations) storage. Repeated
measurements cannot increase the cohort denominator or receive extra votes.
"""
from __future__ import annotations

from collections import defaultdict
from copy import deepcopy
import math
from statistics import median
from typing import Any

from .identity import digest

WEIGHTS = {"measured": 0.4, "structural": 0.3, "discovery": 0.3}
PRIORITY_POLICY_VERSION = "priority-v2.3"
COMPARABILITY_VERSION = "scientific-cohort-v2"
MEASUREMENT_CONTEXT = ("target_id", "species", "assay_type", "endpoint", "unit", "assay_context")
DOCKING_CONTEXT = ("target_id", "receptor_id", "pocket_id", "engine", "protocol_id", "unit")


def validate_requested_action(value):
    """A research hypothesis, never a verified biological effect."""
    if not isinstance(value, str) or value not in {'any', 'activate', 'inhibit'}:
        raise ValueError('desired_action must be any, activate or inhibit')
    return value


def bind_target_goals(candidates, *, desired_action=None, marker=None):
    'Validate the selected structure and retain explicit calculation provenance.'
    from .planning import target_relations
    relations = target_relations()
    output = deepcopy(candidates)
    for row in output:
        binding = row.get("target_goal_binding")
        binding = binding if isinstance(binding, dict) else {}
        scope = marker or binding.get("marker") or row.get("marker")
        goal = desired_action
        origin = "campaign" if isinstance(goal, str) and goal in {"activate", "inhibit"} else None
        explicit_unknown = desired_action is not None or binding.get("explicit_unknown") is True
        if origin is None and desired_action is None:
            goal = row.get("desired_action")
            if isinstance(goal, str) and goal in {"activate", "inhibit"}:
                origin = (binding.get("origin") if binding.get("status") == "resolved" else None) or "candidate"
        support = row.get("goal_support")
        support = support if isinstance(support, dict) else {}
        requested = row.get("requested_action", support.get("requested_action"))
        if desired_action is None and requested is not None and origin not in {"campaign", "candidate"}:
            if isinstance(requested, str) and requested in {"activate", "inhibit"}:
                goal, origin = requested, "source_requested_action"
            else:
                goal, origin, explicit_unknown = None, None, True
        if origin is None and not explicit_unknown:
            goal = support.get("required_action")
            explicit_unknown = goal is not None and (not isinstance(goal, str) or goal not in {"activate", "inhibit"})
            if isinstance(goal, str) and goal in {"activate", "inhibit"}:
                origin = "source_goal_support"
        if origin is None and not explicit_unknown:
            target = str(row.get("target_id", ""))
            gene = target.split(":")[1] if target.startswith("gene:") else target
            matches = [r for r in relations if r["gene"] == gene
                       and (not scope or r["marker_id"] == scope)]
            actions = {r["desired_action"] for r in matches}
            if len(actions) == 1 and actions <= {"activate", "inhibit"}:
                goal, origin = actions.pop(), "configured_marker" if scope else "unambiguous_configured_gene"
        row["desired_action"] = goal if origin else None
        row["target_goal_binding"] = {"status": "resolved" if origin else "unresolved",
                                      "origin": origin, "marker": scope, "explicit_unknown": explicit_unknown and not origin}
        observations = row.get("measurements") or []
        observations = [observations] if isinstance(observations, dict) else observations
        observations = observations if isinstance(observations, (list, tuple)) else []
        actions = {m.get("action") for m in observations if isinstance(m, dict)
                   and m.get("is_measured") is True and isinstance(m.get("action"), str)
                   and m.get("target_id") == row.get("target_id")
                   and m.get("status") not in {"failed", "invalid", "not_evaluated", "unmeasured",
                       "measurement_failed", "measurement_unavailable"}
                   and finite(m.get("value")) is not None
                   and m["action"] in {"activate", "inhibit"}}
        if origin and actions and goal not in actions and row.get("role", "candidate") == "candidate":
            row["opposing_direction"] = True
            row["role"] = "counter_control"
    return output


def comparable_protocol(protocol: dict[str, Any]) -> dict[str, Any]:
    """Project execution provenance onto the scientific comparison contract.

    Binary hashes, seeds, receptor, pocket and search settings remain cohort
    variables. Batch membership/size and installation paths govern execution
    and exact replay, not the population in which an affinity is ranked. The
    full original protocol remains in receipts and cache identities.
    """
    result = deepcopy(protocol)
    result.pop("batch", None)
    result.pop("batch_membership_sha256", None)
    # An installed fallback is not the engine that produced this observation.
    # Actual CPU fallback rows carry their own Vina protocol and binary hash.
    result.pop("cpu_fallback_identity", None)
    for key in ("engine_identity", "cpu_fallback_identity"):
        if isinstance(result.get(key), dict):
            result[key].pop("path", None)
    return result


def finite(value: Any) -> float | None:
    if isinstance(value, bool) or not isinstance(value, (float, int)):
        return None
    try:
        number = float(value)
    except OverflowError:
        return None
    return number if math.isfinite(number) else None


def _observations(candidate: dict[str, Any], channel: str):
    raw = candidate.get("measurements" if channel == "measured" else "docking") or []
    records = [raw] if isinstance(raw, dict) else raw
    if not isinstance(records, (list, tuple)):
        return
    fields = MEASUREMENT_CONTEXT if channel == "measured" else DOCKING_CONTEXT
    for row in records:
        if not isinstance(row, dict):
            continue
        # Do not impute target IDs or protocol metadata from the candidate row.
        context = tuple(row.get(key) for key in fields)
        if any(not isinstance(value, str) or not value.strip() for value in context):
            continue
        if row["target_id"] != candidate.get("target_id"):
            continue
        if row.get("ranking_eligible") is False or row.get("artifact_valid") is False or row.get("status") in {
            "failed", "invalid", "unsupported", "not_evaluated", "unmeasured",
            "measurement_failed", "measurement_unavailable"}:
            continue
        if channel == "measured":
            if row.get("relation", "=") != "=" or row.get("is_measured") is not True:
                continue  # Censored values remain in source, not exact ranks.
            if row.get("opposing_direction") is True:
                continue
            value = finite(row.get("value"))
        else:
            # New discovery runs freeze one primary protocol on the whole
            # population. Confirmation, alternate receptors and fallback
            # engines remain visible without adding votes to this channel.
            if ("priority_primary_protocol_id" in candidate and
                    row.get("protocol_id") != candidate["priority_primary_protocol_id"]):
                continue
            if "priority_primary_protocol_id" in candidate and candidate["priority_primary_protocol_id"] is None:
                continue
            # Score only admitted pose artifacts; missing artifact checks are not passes.
            if row.get("artifact_valid") is not True:
                continue
            value = finite(row.get("score"))
        higher = row.get("higher_is_better")
        if value is None or not isinstance(higher, bool):
            continue
        direction = row.get("action", "unspecified")
        if not isinstance(direction, str) or not direction:
            continue
        desired = candidate.get("desired_action")
        if channel == "measured" and direction in {"activate", "inhibit"} and desired not in {"activate", "inhibit"}:
            continue
        if channel == "measured" and desired and direction not in {desired, "unspecified"}:
            continue
        # Lower/higher semantics and action direction are part of comparability.
        key = (channel, *context, higher, direction)
        yield key, value if higher else -value


def _percentiles(values: dict[str, float]) -> dict[str, float]:
    ordered = sorted(values.items(), key=lambda item: (item[1], item[0]))
    n = len(ordered)
    if n == 1:
        return {ordered[0][0]: 50.0}
    result: dict[str, float] = {}
    start = 0
    while start < n:
        end = start + 1
        while end < n and ordered[end][1] == ordered[start][1]:
            end += 1
        percentile = 100.0 * (start + end - 1) / (2 * (n - 1))
        for candidate_id, _ in ordered[start:end]:
            result[candidate_id] = percentile
        start = end
    return result


def rank_candidates(candidates: list[dict[str, Any]], *,
                    policy_version: str = PRIORITY_POLICY_VERSION) -> list[dict[str, Any]]:
    """Augment copied rows, preserving all original observations including nulls.

    Sort descending priority then candidate ID within target. Purchase preference
    belongs to selection/UI ordering, not a hidden numerical score adjustment.
    Cohort n=1 and equal-score ties contribute neutral 50. Coverage measures
    channel availability, while uncertainty additionally penalizes tiny cohorts.
    """
    out = bind_target_goals(candidates)
    identities: set[tuple[str, str]] = set()
    groups: dict[tuple, dict[str, list[float]]] = defaultdict(lambda: defaultdict(list))
    for candidate in out:
        cid, target = candidate.get("candidate_id"), candidate.get("target_id")
        if not isinstance(cid, str) or not cid or not isinstance(target, str) or not target:
            raise ValueError("candidate_id and target_id must be nonempty strings")
        if (target, cid) in identities:
            raise ValueError("duplicate candidate_id within target")
        identities.add((target, cid))
        for channel in ("measured", "structural"):
            for key, value in _observations(candidate, channel):
                groups[key][cid].append(value)
    ranked: dict[tuple[str, str], dict[str, list[dict]]] = defaultdict(lambda: defaultdict(list))
    for key, group in groups.items():
        percentiles = _percentiles({cid: median(values) for cid, values in group.items()})
        population_id = digest(sorted((cid, median(values)) for cid, values in group.items()))
        for cid, percentile in percentiles.items():
            ranked[(key[1], cid)][key[0]].append({
                "cohort_id": digest(key), "cohort_context": list(key[1:]),
                "cohort_size": len(group), "percentile": percentile,
                "cohort_population_sha256": population_id,
                "candidate_observation_count": len(group[cid]),
                "aggregation": "median within candidate then midrank percentile"})
    populations = defaultdict(list)
    for row in out:
        populations[row["target_id"]].append(row["candidate_id"])


    population_hashes = {target: digest(sorted(ids)) for target, ids in populations.items()}
    for candidate in out:
        channels: dict[str, dict[str, Any]] = {}
        terms: dict[str, float] = {}
        coverage = reliability = 0.0
        for name, weight in WEIGHTS.items():
            evidence = ranked[(candidate["target_id"], candidate["candidate_id"])].get(name, [])
            if name == "discovery":
                score = finite(candidate.get("discovery_score"))
                if score is not None and not 0 <= score <= 100:
                    score = None
                recorded_coverage = candidate.get("discovery_coverage", 1.0)
                strength = finite(recorded_coverage)
                strength = strength if strength is not None and 0 <= strength <= 1 and score is not None else 0.0
            else:
                score = median([item["percentile"] for item in evidence]) if evidence else None
                # A singleton establishes observation availability, not relative rank.
                strength = (median([min(1.0, (item["cohort_size"] - 1) / 19)
                                    for item in evidence]) if evidence else 0.0)
            channels[name] = {"score": score, "weight": weight, "cohorts": evidence,
                              "status": "observed" if score is not None else "not_evaluated"}
            terms[name] = 50.0 if score is None else score
            coverage += weight * (strength if name == "discovery" else score is not None)
            reliability += weight * strength
        candidate.update(priority_score_v2=round(sum(WEIGHTS[k] * v for k, v in terms.items()), 6),
            priority_policy="priority-v2.4-primary" if "priority_primary_protocol_id" in candidate else policy_version, priority_weights=dict(WEIGHTS),
            priority_population_sha256=population_hashes[candidate["target_id"]],
            priority_channels=channels, priority_coverage=round(coverage, 6),
            priority_missing_channels=[key for key, value in channels.items() if value["score"] is None],
            priority_imputed_weight=round(sum(WEIGHTS[key] for key, value in channels.items()
                                              if value["score"] is None), 6),
            priority_missing_policy="neutral-50-placeholder-v1; not a measured percentile",
            priority_score_observed_channels=(round(sum(WEIGHTS[key] * value["score"]
                for key, value in channels.items() if value["score"] is not None)
                / sum(WEIGHTS[key] for key, value in channels.items() if value["score"] is not None), 6)
                if any(value["score"] is not None for value in channels.values()) else None),
            priority_scientific_evidence_status=("OBSERVED" if any(
                channels[key]["score"] is not None for key in ("measured", "structural")) else "NOT_EVALUATED"),
            priority_uncertainty=round(1.0 - reliability, 6),
            priority_uncertainty_method="heuristic missingness and cohort-size; not statistical interval",
            priority_interpretation="experimental priority heuristic; missing channels use an unmeasured neutral placeholder that can outrank low observations; compare channels and coverage, not efficacy or binding probability")
    out.sort(key=lambda row: (row["target_id"], -row["priority_score_v2"], row["candidate_id"]))
    ranks: dict[str, int] = defaultdict(int)
    for row in out:
        ranks[row["target_id"]] += 1
        row["priority_rank_v2"] = ranks[row["target_id"]]
    return out
