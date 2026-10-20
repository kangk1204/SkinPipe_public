"""Deterministic same-target experiment selection with distinct control accounting."""
from __future__ import annotations

from collections import Counter
from copy import deepcopy
from typing import Any

from .identity import chemical_identity, candidate_identity, digest
from .ranking import finite, rank_candidates
from .procurement import has_purchase_record, documented_offers


def _purchase(row: dict[str, Any]) -> bool:
    return has_purchase_record(row)


def _order(row: dict[str, Any]):
    return (not _purchase(row), -(finite(row.get("priority_score_v2")) or 0), row["candidate_id"])


def select_experiment_batch(candidates: list[dict[str, Any]], target_id: str, *,
                            references=(), counts: tuple[int, int, int] = (60, 24, 12),
                            trust_stored_identity: bool = False,
                            require_documented_procurement: bool = False,
                            preserve_frozen_ranking: bool = False) -> dict[str, Any]:
    """Select novel/evidence/exploration quotas, then disclosed same-target backfill.

    Invariant: every selected parent is unique and absent from the control set;
    every selected candidate belongs to target_id. An empty reference set cannot
    establish novelty. Acyclic compounds occupy a separate group. Complexity is
    O(n log n) plus RDKit parsing/scaffold construction, with O(n) storage.
    """
    if not isinstance(target_id, str) or not target_id:
        raise ValueError("target_id is required")
    if len(counts) != 3 or any(isinstance(x, bool) or not isinstance(x, int) or x < 0 for x in counts):
        raise ValueError("counts must be three nonnegative integers")
    rows = [deepcopy(row) for row in candidates if row.get("target_id") == target_id]
    # Stored scores belong to their original cohorts. Planning and experiment
    # selection must rank the identical pool, with the same hypothesis.
    if not preserve_frozen_ranking:
        rows = rank_candidates(rows)
    if len({row["candidate_id"] for row in rows}) != len(rows):
        raise ValueError("duplicate candidate IDs in target")
    reference_cores: set[str] = set()
    reference_exact: list[str] = []
    for reference in references:
        if isinstance(reference, dict):
            if reference.get("target_id", target_id) != target_id:
                continue
            smiles = reference.get("smiles")
        else:
            smiles = reference
        identity = chemical_identity(smiles)
        if identity["structure_status"] == "valid":
            reference_exact.append(identity["exact_structure_id"])
            reference_cores.add(identity["scaffold_id"])
    excluded: list[dict] = []
    controls: list[dict] = []
    eligible: list[dict] = []
    for row in rows:
        identity = candidate_identity(row, trust_stored=trust_stored_identity)
        row["identity"] = identity
        row["parent_id"] = identity["parent_id"]
        row.update({key: identity.get(key) for key in ("scaffold_id", "scaffold_smiles", "scaffold_kind")})
        row["novel_scaffold"] = ((row["scaffold_kind"] == "murcko"
                                  and row["scaffold_id"] not in reference_cores)
                                 if reference_exact else None)
        row["novelty_scope"] = "frozen target reference Murcko set; not literature or patent novelty"
        if identity["structure_status"] != "valid" or identity["parent_ambiguous"] or row.get("identity_conflict"):
            excluded.append({"candidate_id": row["candidate_id"], "reason": "unresolved_material_identity"})
            continue
        if require_documented_procurement and row.get('discovery_filter_pass') is not True:
            excluded.append({"candidate_id": row["candidate_id"], "reason": "frozen_discovery_filter"})
            continue
        if row.get("role") in {"control", "positive_control", "negative_control", "counter_control", "reference_control"} or row.get("opposing_direction") is True:
            controls.append(row)
        elif require_documented_procurement and not documented_offers(row):
            excluded.append({"candidate_id": row["candidate_id"], "reason": "procurement_unconfirmed"})
        elif not _purchase(row):
            excluded.append({"candidate_id": row["candidate_id"], "reason": "purchase_availability_unconfirmed"})
        else:
            eligible.append(row)
    controls.sort(key=_order)
    control_parents: set[str] = set()
    unique_controls = []
    for row in controls:
        if row["parent_id"] not in control_parents:
            control_parents.add(row["parent_id"])
            unique_controls.append(row)
    eligible.sort(key=_order)
    chosen: list[dict] = []
    used_parents = set(control_parents)
    used_cores: Counter = Counter()
    cohort_counts: Counter = Counter()

    def take(row, cohort, reason):
        if row["parent_id"] in used_parents:
            return False
        used_parents.add(row["parent_id"])
        used_cores[row["scaffold_id"]] += 1
        row.update(selection_group=cohort, selection_reason=reason, selection_rank=len(chosen) + 1)
        chosen.append(row)
        cohort_counts[cohort] += 1
        return True

    for row in eligible:
        if cohort_counts["novel_scaffold"] >= counts[0]:
            break
        if row["novel_scaffold"] and not used_cores[row["scaffold_id"]]:
            take(row, "novel_scaffold", "distinct Murcko core absent from frozen reference set")
    def evidence_score(row):
        channel = row.get("priority_channels", {}).get("measured", {})
        return finite(channel.get("score"))
    evidence_rows = sorted(eligible, key=lambda row: (
        not _purchase(row), -(evidence_score(row) if evidence_score(row) is not None else -1),
        *_order(row)[1:]))
    for row in evidence_rows:
        if cohort_counts["evidence"] >= counts[1]:
            break
        if evidence_score(row) is not None:
            take(row, "evidence", "condition-compatible measured target evidence")
    exploration_rows = sorted(eligible, key=lambda row: (
        not _purchase(row), -(finite(row.get("priority_uncertainty")) or 0), *_order(row)[1:]))
    for row in exploration_rows:
        if cohort_counts["exploration"] >= counts[2]:
            break
        associated = row.get("target_association") is True or evidence_score(row) is not None or bool(row.get("target_evidence"))
        uncertain = (finite(row.get("priority_uncertainty")) or 0) > 0 or row.get("conflicting_evidence") is True
        if associated and uncertain:
            take(row, "exploration", "target-related candidate with missing or conflicting evidence")
    deficits = {name: max(0, requested - cohort_counts[name]) for name, requested in zip(
        ("novel_scaffold", "evidence", "exploration"), counts)}
    # Diversity-first passes avoid repeatedly rescanning/sorting the whole pool.
    remaining = [row for row in eligible if row["parent_id"] not in used_parents]
    buckets: dict[str, list[dict]] = {}
    for row in remaining:
        buckets.setdefault(row["scaffold_id"], []).append(row)
    rounds: list[tuple[int, tuple, dict]] = []
    for core, bucket in buckets.items():
        for offset, row in enumerate(bucket):
            rounds.append((used_cores[core] + offset, _order(row), row))
    rounds.sort(key=lambda item: (item[1][0], item[0], item[1][1:]))
    for _, _, row in rounds:
        if len(chosen) >= sum(counts):
            break
        take(row, "backfill", "same-target diversity backfill; requested category unavailable")
    selected_parents = {row["parent_id"] for row in chosen}
    if len(selected_parents) != len(chosen) or not selected_parents.isdisjoint(control_parents):
        raise ValueError("Selected candidates and controls must have distinct parents")
    policy = {"version": "experimental-96-v2.3", "target_id": target_id,
        "procurement_requirement": "dated-exact-form-v1" if require_documented_procurement else "catalog-record-v1",
        "counts": dict(zip(("novel_scaffold", "evidence", "exploration"), counts)),
        "purchase_preference": True, "purchase_required": True, "reference_exact_ids": sorted(set(reference_exact)),
        "reference_scaffold_ids": sorted(reference_cores),
        "reference_status": "frozen" if reference_exact else "not_established"}
    return {"target_id": target_id, "candidates": chosen, "controls": unique_controls,
        "candidate_count": len(chosen), "control_count": len(unique_controls),
        "requested_candidate_count": sum(counts), "category_deficits": deficits,
        "backfill_count": cohort_counts["backfill"], "unfilled_count": sum(counts) - len(chosen),
        "selected_counts": dict(cohort_counts), "excluded_from_batch": excluded,
        "population_preserved": True, "policy": policy, "policy_id": digest(policy)}
