"""Deterministic research stage plans, with explicit dependencies and no execution claims."""
from __future__ import annotations

from collections import Counter, defaultdict
from copy import deepcopy
import hashlib
from pathlib import Path
from typing import Any

import yaml

from .identity import chemical_identity, digest
from .confirmation import confirmation_policy as normalize_confirmation_policy
from .policy import COMPUTE_DEFAULTS, POLICY_VERSION, assess_candidate
from .algorithm_profiles import compute_profile as resolve_compute_profile
from .ranking import PRIORITY_POLICY_VERSION, bind_target_goals, finite, rank_candidates
from .procurement import has_purchase_record


def target_relations(config_path=None) -> list[dict[str, Any]]:
    """Read configured marker/node relations without network identity inference.

    Configuration provides a gene-level hypothesis, not verified UniProt mapping.
    Each marker/node action is retained separately even if the node is reused.
    Endpoint-only markers have no protein and cannot accidentally become receptors.
    """
    path = Path(config_path) if config_path else Path(__file__).resolve().parents[2] / "config" / "skin_targets.yaml"
    raw = path.read_bytes()
    config = yaml.safe_load(raw)
    if not isinstance(config, dict) or not isinstance(config.get("markers"), list):
        raise ValueError("invalid target registry")
    source_hash = hashlib.sha256(raw).hexdigest()
    relations = []
    for marker in config["markers"]:
        marker_id = marker["marker"]
        mode = marker.get("mode")
        nodes = marker.get("nodes") or []
        if mode not in {"direct", "proxy", "endpoint"}:
            raise ValueError("unknown target relation mode")
        if mode == "endpoint" and nodes:
            raise ValueError("endpoint marker cannot silently declare protein nodes")
        for node in nodes or [None]:
            gene = node.get("gene") if node else None
            action = node.get("node_action", "unknown") if node else None
            if node and (not isinstance(gene, str) or not gene or action not in {"activate", "inhibit", "unknown"}):
                raise ValueError("invalid gene/action relation")
            relation = {"marker_id": marker_id, "marker_direction": marker.get("direction"),
                "mode": mode, "gene": gene, "actual_target_id": "gene:" + gene + ":human" if gene else None,
                "desired_action": action, "species": "human" if gene else None,
                "identity_status": "configured_gene_unresolved_protein" if gene else "nonprotein_endpoint",
                "is_proxy": mode == "proxy", "preferred": bool(node and (node.get("preferred") or node.get("primary"))),
                "goal_support": deepcopy(node.get("goal_support", {"status": "not_established"})) if node else None,
                "rationale": node.get("rationale") if node else marker.get("note"),
                "action_rationale": node.get("action_why") if node else None,
                "structure_usable": node.get("structure_usable") if node else False,
                "config_sha256": source_hash,
                "interpretation": "node binding does not establish marker endpoint effect"}
            relation["relation_id"] = digest(relation)
            relations.append(relation)
    return relations


def _preference(row):
    return (not has_purchase_record(row), -row["priority_score_v2"], row["candidate_id"])


def _diverse(rows, limit):
    """One row per parent; interleave scaffold groups without quadratic rescans."""
    buckets = defaultdict(list)
    seen = set()
    for row in sorted(rows, key=_preference):
        parent = row["identity"]["parent_id"]
        if parent in seen:
            continue
        seen.add(parent)
        buckets[row["identity"]["scaffold_id"]].append(row)
    rounds = [(depth, _preference(row), row) for bucket in buckets.values()
              for depth, row in enumerate(bucket)]
    rounds.sort(key=lambda item: (item[1][0], item[0], item[1][1:]))
    return [item[2] for item in rounds[:limit]]


def plan_target(candidates: list[dict[str, Any]], target_id: str, *, snapshot_id: str,
                references=(), target_context: dict[str, Any] | None = None,
                compute_budget_seconds: int = 28_800,
                consumed_compute_seconds: float = 0, trust_stored_identity: bool = False,
                confirmation_policy: dict[str, Any] | None = None,
                compute_profile: str | None = "balanced-research") -> dict[str, Any]:
    """Prepare a frozen plan; returning it never submits or completes a job.

    Downstream candidate lists are provisional and must be reselected against
    verified upstream receipts before dispatch. The stable budget_scope must be
    enforced transactionally by the queue across all stage/chunk requests.
    """
    if not isinstance(target_id, str) or not target_id or not isinstance(snapshot_id, str) or not snapshot_id:
        raise ValueError("target_id and immutable snapshot_id are required")
    budget = finite(compute_budget_seconds)
    consumed = finite(consumed_compute_seconds)
    if budget is None or not 0 < budget <= 604800 or consumed is None or consumed < 0:
        raise ValueError("invalid target computation budget")
    context = deepcopy(target_context or {})
    profile = resolve_compute_profile(compute_profile)
    if profile and budget > profile["budget"]["target_compute_seconds"]:
        raise ValueError("Balanced target budget cannot exceed8h actual compute")
    limits = {**dict(COMPUTE_DEFAULTS), **(profile["counts"] if profile else {})}
    confirmation_config = normalize_confirmation_policy(confirmation_policy if confirmation_policy is not None else
        {'ranking':'gpu_docking_score'} if profile and profile['version']=='generalized-research-v1' else {})
    if profile and (confirmation_config["candidate_limit"] > limits["confirmation"]
                    or confirmation_config["seeds"] != profile["multiplicities"]["confirmation_seeds"]):
        raise ValueError("Balanced confirmation permits at most96 candidates and seed137 once")
    if context.get("target_id", target_id) != target_id:
        raise ValueError("target context does not match candidate target")
    rows = bind_target_goals([r for r in candidates if r.get("target_id") == target_id],
        desired_action=context.get("desired_action"), marker=context.get("marker_id"))
    if not (profile and profile['version']=='generalized-research-v1' and all(r.get('discovery_run_id') for r in rows)):
        rows=rank_candidates(rows)
    all_ids = [r["candidate_id"] for r in rows]
    valid_identity = context.get("target_identity_valid") is True
    receptor = context.get("receptor_id")
    pocket = context.get("pocket_id")
    actual_protein = context.get("actual_protein_id")
    structure_context = bool(valid_identity and actual_protein and receptor and pocket)
    additional = context.get("additional_receptors", [])
    if not isinstance(additional, list):
        raise ValueError("additional_receptors must be a list of verified contexts")
    extra_contexts = [r for r in additional if isinstance(r, dict)
                      and r.get("target_identity_valid") is True
                      and r.get("actual_protein_id") == actual_protein
                      and isinstance(r.get("receptor_id"), str) and r["receptor_id"]
                      and isinstance(r.get("pocket_id"), str) and r["pocket_id"]]
    unique_extra = {}
    for extra in extra_contexts:
        key = (extra["receptor_id"], extra["pocket_id"])
        if key != (receptor, pocket):
            unique_extra.setdefault(key, extra)
    extra_contexts = [unique_extra[key] for key in sorted(unique_extra)][:limits["additional_receptors"]]
    assessed = []
    for row in rows:
        item = assess_candidate(row, target_identity_valid=valid_identity, trust_stored_identity=trust_stored_identity)
        if not structure_context:
            item["states"]["compute"]["allowed"] = False
            item["states"]["compute"]["status"] = "held"
            item["states"]["compute"]["reasons"].append("verified_protein_receptor_pocket_required")
        assessed.append(item)
    control_roles = {"control", "positive_control", "negative_control", "counter_control", "reference_control"}
    controls = [r for r in assessed if r.get("role") in control_roles]
    candidate_rows = [r for r in assessed if r.get("role") not in control_roles]
    control_parents = {r["identity"]["parent_id"] for r in controls if r["identity"]["parent_id"]}
    purchase_eligible = [r for r in candidate_rows if has_purchase_record(r)]
    generalized = bool(profile and profile["version"] == "generalized-research-v1")
    if generalized:
        purchase_eligible = [r for r in candidate_rows if r.get("discovery_filter_pass") is True]
    purchase_held = [{"candidate_id": r["candidate_id"], "reason": "purchase_availability_unconfirmed",
                      "searchable": True} for r in candidate_rows if not has_purchase_record(r)]
    cheap = sorted(purchase_eligible, key=_preference)[:limits["cheap_rerank"]]
    broad = _diverse([r for r in cheap if r["states"]["compute"]["allowed"]
                      and r["identity"]["parent_id"] not in control_parents], limits["broad_docking"])
    if generalized:
        from .catalog_hopping import select_docking
        ids = select_docking(
            [r for r in cheap if r["states"]["compute"]["allowed"]] + controls,
            limits["broad_docking"],
        )
        by_id = {r["candidate_id"]: r for r in cheap}
        broad = [by_id[cid] for cid in ids]
    confirmation = _diverse(broad, limits["confirmation"])
    if confirmation_config:


        confirmation = sorted(broad, key=lambda row: (-row["priority_score_v2"], row["candidate_id"]))[
            :confirmation_config["candidate_limit"]]
    sensitivity = _diverse(confirmation, limits["receptor_sensitivity"])
    control_pool = _diverse([r for r in controls if r["states"]["compute"]["allowed"]], len(controls))
    reference_records = []
    for reference in references:
        if isinstance(reference, dict) and reference.get("target_id", target_id) != target_id:
            continue
        smiles = reference.get("smiles") if isinstance(reference, dict) else reference
        identity = chemical_identity(smiles)
        reference_records.append({"smiles": smiles, "exact_structure_id": identity["exact_structure_id"],
                                  "scaffold_id": identity["scaffold_id"], "structure_status": identity["structure_status"]})
    reference_records.sort(key=lambda row: (row["exact_structure_id"] or "", str(row["smiles"])))
    reference_snapshot = {"target_id": target_id, "records": reference_records,
                          "status": "frozen" if any(r["exact_structure_id"] for r in reference_records) else "not_established"}
    reference_snapshot["sha256"] = digest(reference_snapshot)
    policy = {"version": POLICY_VERSION, "defaults": limits,
              "selection": "purchase required, " + PRIORITY_POLICY_VERSION + ", scaffold round-robin, unique parent",
              "purchase_required": True,
              "reference_sha256": reference_snapshot["sha256"]}
    if confirmation_config:
        policy["confirmation"] = confirmation_config
    if profile:
        policy["compute_profile"] = profile
    policy_id = digest(policy)
    budget_scope = digest({"snapshot_id": snapshot_id, "target_id": target_id, "policy_id": policy_id})
    input_signature = digest([{"candidate_id": r["candidate_id"],
        "exact_structure_id": r["identity"]["exact_structure_id"],
        "priority_score_v2": r["priority_score_v2"], "purchasable": r.get("purchasable"),
        "compute_allowed": r["states"]["compute"]["allowed"], "role": r.get("role")}
        for r in sorted(assessed, key=lambda item: item["candidate_id"])])
    plan_id = digest({"budget_scope": budget_scope, "input_signature": input_signature,
                      "target_context": context, "compute_budget_seconds": budget})
    pools = {"cheap_rerank": cheap, "broad_docking": broad,
             "confirmation": confirmation, "receptor_sensitivity": sensitivity, "automatic_deep": []}
    dependencies = {"cheap_rerank": [], "broad_docking": [], "confirmation": ["broad_docking"],
                    "receptor_sensitivity": ["confirmation"], "automatic_deep": []}
    job_ids = {stage: digest([plan_id, stage]) for stage in pools}
    stages, jobs = {}, []
    for stage, pool in pools.items():
        ids = [r["candidate_id"] for r in pool]
        control_ids = [r["candidate_id"] for r in control_pool] if stage == "confirmation" else []
        omissions = []
        if stage != "cheap_rerank" and not structure_context:
            omissions.append("verified_protein_receptor_pocket_required")
        if stage == "automatic_deep":
            omissions.append("automatic_deep_disabled_by_policy")
        if stage == "receptor_sensitivity" and not extra_contexts:
            omissions.append("additional_receptor_contexts_unavailable")
        if consumed >= budget and stage != "cheap_rerank":
            omissions.append("budget_deferred")
        status = "held" if omissions else "awaiting_upstream" if dependencies[stage] else "planned"
        if not ids and not control_ids and not omissions:
            status = "empty_pool"
        request = {"plan_id": plan_id, "snapshot_id": snapshot_id, "target_id": target_id,
            "target_context": context, "policy_id": policy_id, "reference_sha256": reference_snapshot["sha256"],
            "budget_scope": budget_scope, "target_budget_seconds": budget,
            "stage": stage, "candidate_ids": ids, "control_ids": control_ids,
            "candidate_limit": (confirmation_config["candidate_limit"]
                                if stage == "confirmation" and confirmation_config else limits[stage]),
            "seed": (confirmation_config["seeds"][0] if confirmation_config else 137)
                    if stage == "confirmation" else COMPUTE_DEFAULTS["seed"],
            "seeds": (confirmation_config["seeds"] if confirmation_config else [137, 271])
                     if stage == "confirmation" else [COMPUTE_DEFAULTS["seed"]],
            "requires_stage_receipts": [job_ids[s] for s in dependencies[stage]],
            "selection_provisional": bool(dependencies[stage]),
            "engine_preference": ([] if stage == "cheap_rerank" else
                                  ["unidock", "vina_cpu"] if stage == "broad_docking" else ["vina_cpu"]),
            "additional_receptors": extra_contexts if stage == "receptor_sensitivity" else []}
        if profile:
            request["compute_profile_sha256"] = profile["sha256"]
        if stage == "confirmation" and confirmation_config:
            request["confirmation_policy"] = deepcopy(confirmation_config)
        stages[stage] = {"status": status, "candidate_count": len(ids), "control_count": len(control_ids),
            "candidate_ids": ids, "control_ids": control_ids, "omissions": omissions,
            "requires_stage_receipts": request["requires_stage_receipts"]}
        if ids or control_ids:
            jobs.append({"job_key": job_ids[stage], "kind": "research_stage", "status": status,
                         "runnable": status == "planned", "submitted": False, "request": request})
    broad_ids = {r["candidate_id"] for r in broad}
    cheap_ids = {r["candidate_id"] for r in cheap}
    deferred = []
    for row in candidate_rows:
        cid = row["candidate_id"]
        if cid in broad_ids:
            continue
        reason = ("compute_held" if not row["states"]["compute"]["allowed"] else
                  "purchase_availability_unconfirmed" if not has_purchase_record(row) else
                  "cheap_rerank_cap" if cid not in cheap_ids else "broad_docking_cap_or_parent_duplicate")
        deferred.append({"candidate_id": cid, "reason": reason,
                         "details": row["states"]["compute"]["reasons"], "searchable": True})
    return {"plan_id": plan_id, "target_id": target_id, "snapshot_id": snapshot_id,
        "status": "planned_not_submitted", "population_count": len(rows), "population_ids": all_ids,
        "stages": stages, "jobs": jobs, "deferred": deferred, "purchase_held": purchase_held,
        "deferred_counts": dict(Counter(row["reason"] for row in deferred)),
        "policy": policy, "policy_id": policy_id, "reference_snapshot": reference_snapshot,
        "budget_scope": budget_scope, "budget_seconds": budget,
        "consumed_compute_seconds": consumed, "remaining_compute_seconds": max(0, budget - consumed),
        "queue_eta_seconds": None, "compute_eta_seconds": None,
        "additional_receptor_contexts_omitted": len(additional) - len(extra_contexts),
        "optional_stages": {"mmgbsa": limits["optional_mmgbsa"],
                            "advanced_models": profile["counts"]["complex_models"] if profile else 12},
        "compute_profile": profile,
        "candidate_population_preserved": True, "actual_completed_count": 0,
        "warnings": ["Downstream pools are provisional until verified upstream results are ranked",
                     "Reference-relative scaffold difference is not literature or patent novelty"]}
