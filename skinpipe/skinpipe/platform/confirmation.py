"""Frozen, score-ordered CPU checks after verified broad docking."""
from __future__ import annotations

from copy import deepcopy

from .identity import digest, candidate_identity
from .ranking import comparable_protocol, finite

CONFIRMATION_VERSION = "score-ordered-cpu-confirmation-v1"


def confirmation_policy(value):
    """New plans default to96 candidates and one seed; frozen plans stay unchanged."""
    if not isinstance(value, dict) or set(value) - {"candidate_limit", "ranking", "seeds"}:
        raise ValueError("confirmation_policy accepts candidate_limit, ranking and seeds")
    limit = value.get("candidate_limit", 96)
    if isinstance(limit, bool) or not isinstance(limit, int) or not 1 <= limit <= 1000:
        raise ValueError("confirmation candidate_limit must be an integer from1 to1000")
    ranking = value.get("ranking", "priority_score_v2")
    if not isinstance(ranking, str) or ranking not in {"priority_score_v2", "gpu_docking_score"}:
        raise ValueError("confirmation ranking must be priority_score_v2 or gpu_docking_score")
    seeds = value.get("seeds", [137])
    if (not isinstance(seeds, list) or not 1 <= len(seeds) <= 2
            or any(isinstance(seed, bool) or not isinstance(seed, int) or not 0 < seed < 2**31 for seed in seeds)
            or len(set(seeds)) != len(seeds)):
        raise ValueError("confirmation seeds must contain one or two distinct positive32-bit integers")
    return {"version": CONFIRMATION_VERSION, "candidate_limit": limit, "ranking": ranking,
            "seeds": deepcopy(seeds), "engine": "vina", "selection": "score_order_unique_parent",
            "higher_is_better": ranking == "priority_score_v2", "controls_separate": True}


def select_confirmation(ranked, results, policy):
    """Select existing parent-unique rows; never turn failed or absent scores into zero.

    Raw GPU affinities require one scientific cohort. CPU fallbacks remain in
    the source population, but cannot be mixed into a GPU-only score ordering.
    """
    successful = {row["candidate_id"] for row in results if row.get("status") == "complete"}
    eligible = [row for row in ranked if row["candidate_id"] in successful]
    if policy["ranking"] == "priority_score_v2":
        eligible = [row for row in eligible if finite(row.get("priority_score_v2")) is not None]
        eligible.sort(key=lambda row: (-row["priority_score_v2"], row["candidate_id"]))
    else:
        scores, cohorts = {}, set()
        for row in results:
            if row.get("status") != "complete" or row.get("engine") != "unidock":
                continue
            score = finite(row.get("affinity"))
            protocol = row.get("protocol") or {}
            if score is None:
                continue
            if (protocol.get("engine") != "unidock" or not protocol.get("receptor_sha256")
                    or not (protocol.get("engine_identity") or {}).get("sha256")):
                raise ValueError("GPU score ordering requires verified scientific protocol identity")
            cohorts.add(digest(comparable_protocol(protocol)))
            scores[row["candidate_id"]] = score
        if len(cohorts) > 1:
            raise ValueError("Raw GPU scores span incompatible cohorts; use priority_score_v2")
        eligible = [row for row in eligible if row["candidate_id"] in scores]
        eligible.sort(key=lambda row: (scores[row["candidate_id"]], row["candidate_id"]))
    controls = set()
    def parent(row):
        identity = candidate_identity(row, trust_stored=True)
        return identity.get("parent_id") or row["candidate_id"]
    for row in ranked:
        if row.get("role") in {"control", "positive_control", "negative_control", "counter_control", "reference_control"} or row.get("opposing_direction"):
            controls.add(parent(row))
    seen, selected = set(controls), []
    for row in eligible:
        key = parent(row)
        if key in seen:
            continue
        seen.add(key)
        selected.append(row["candidate_id"])
        if len(selected) == policy["candidate_limit"]:
            break
    return selected
