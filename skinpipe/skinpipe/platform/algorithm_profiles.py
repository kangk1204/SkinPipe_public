"""One frozen count/workload contract; optional failures do not remove candidates."""
from __future__ import annotations

from copy import deepcopy
import math

from .identity import candidate_identity, digest
from .procurement import has_purchase_record
from .ranking import finite

PROFILE_VERSION = "balanced-research-v1"
_PROFILE = {
    "version": PROFILE_VERSION,
    "counts": {"cheap_rerank": 10000, "broad_docking": 1000, "confirmation": 96,
               "receptor_sensitivity": 16, "additional_receptors": 1,
               "shape_pharmacophore": 2000, "pose_quality": 1000, "admet_uncached": 1000,
               "gnina": 192, "interaction_explanation": 192, "state_sensitivity": 32,
               "boltzina": 48, "complex_models": 4, "optional_mmgbsa": 4,
               "automatic_deep": 0, "experiment": 96},
    "multiplicities": {"conformers": 10, "alignment_references": 5,
                       "prepared_states_total": 2, "representative_poses": 1,
                       "confirmation_seeds": [137]},
    "budget": {"target_compute_seconds": 28800, "optional_inference_seconds": 3600,
               "optional_mmgbsa_seconds": 3600, "cost_reserve_multiplier": 1.25},
    "resources": {"canary_candidates": 16, "canary_max": 32,
                  "gpu_batch_ladder": [8, 16, 32], "reserve_vram_mib": 1536,
                  "peak_memory_multiplier": 1.25, "cpu_fraction": 0.25, "cpu_max": 8},
    "bands": {"cdpl": [1000, 600, 400], "gnina": [128, 32, 32],
              "prolif": [128, 32, 32], "boltzina": [32, 8, 8]},
    "preview_candidates": 300,
    "score_changes": False, "candidate_population_preserved": True,
    "new_model_execution": "requires_exact_asset_real_input_qualification",
}
PROVIDERS = {
    "fpsim2": ("cheap_rerank", "cpu"), "cdpl": ("shape_pharmacophore", "cpu"),
    "dimorphite": ("state_sensitivity", "cpu"), "posebusters": ("pose_quality", "cpu"),
    "prolif": ("interaction_explanation", "cpu"), "gnina": ("gnina", "gpu"),
    "boltzina": ("boltzina", "gpu"), "boltz2": ("complex_models", "gpu"),
    "chai1": ("complex_models", "gpu"), "openadmet": ("admet_uncached", "cpu"),
    "chemeleon": ("admet_uncached", "cpu"), "unidock2": ("broad_docking", "gpu"),
    "openfe": ("automatic_deep", "gpu"), "fpocket": ("receptor_sensitivity", "cpu"),
    "p2rank": ("receptor_sensitivity", "cpu"),
}


def compute_profile(name="balanced-research"):
    if name is None or name == "historical":
        return None
    if not isinstance(name, str) or name not in {"balanced-research", PROFILE_VERSION, "generalized-research", "generalized-research-v1"}:
        raise ValueError("Unknown compute profile")
    value = deepcopy(_PROFILE)
    if name in {"generalized-research", "generalized-research-v1"}:
        value.update(version="generalized-research-v1", selection="reference hopping then same-filter catalog backfill",
                     ranking="priority-v2.4-primary", procurement="dated-exact-form-v1")
        value["multiplicities"].update(primary_seed=42, automatic_additional_states=0,
                                       automatic_complex_models=0)
    value["sha256"] = digest(value)
    return value


def validate_profile(value):
    expected = compute_profile(value.get("version") if isinstance(value, dict) else None)
    if value != expected:
        raise ValueError("Compute profile is not the pinned balanced workload contract")
    return deepcopy(expected)


def provider_limit(provider, *, qualification=False):
    if provider not in PROVIDERS:
        raise ValueError("Unsupported research algorithm")
    stage = PROVIDERS[provider][0]
    limit = _PROFILE["counts"][stage]
    if provider in {"fpocket", "p2rank"}:
        limit = 1  # one receptor search, not 16 ligand invocations
    if qualification:
        limit = min(limit, 16)
    return limit


def select_bands(rows, target_id, provider, *, trust_stored_identity=False):
    'O(n log n) deterministic score/diversity/uncertainty union, unique parents.'
    cap = provider_limit(provider)
    rows = [r for r in rows if r.get("target_id") == target_id]
    if any(not isinstance(r.get("candidate_id"), str) or not r["candidate_id"] for r in rows):
        raise ValueError("Candidate IDs are required")
    if len({r.get("candidate_id") for r in rows}) != len(rows):
        raise ValueError("Duplicate candidate IDs in branch selection")
    controls, eligible, omitted = set(), [], []
    for row in rows:
        if row.get("target_id") != target_id:
            continue
        identity = candidate_identity(row, trust_stored=trust_stored_identity)
        if row.get("role") in {"control", "positive_control", "negative_control", "counter_control"} or row.get("opposing_direction"):
            if identity.get("parent_id"):
                controls.add(identity["parent_id"])
            continue
        if identity.get("structure_status") != "valid" or identity.get("parent_ambiguous") or row.get("identity_conflict"):
            omitted.append({"candidate_id": row["candidate_id"], "reason": "unresolved_material_identity"})
        elif not has_purchase_record(row):
            omitted.append({"candidate_id": row["candidate_id"], "reason": "purchase_availability_unconfirmed"})
        else:
            eligible.append({**deepcopy(row), "identity": identity})
    def order(r):
        return (-(finite(r.get("priority_score_v2")) or 0), r["candidate_id"])
    eligible.sort(key=order)
    selected, used, cores, counts = [], set(controls), set(), {"score": 0, "diversity": 0, "uncertainty": 0, "backfill": 0}
    bands = _PROFILE["bands"].get(provider, [cap, 0, 0])

    def take(row, group):
        parent = row["identity"]["parent_id"]
        if parent in used or len(selected) >= cap:
            return False
        used.add(parent)
        cores.add(row["identity"].get("scaffold_id"))
        row = deepcopy(row)
        row["algorithm_selection_group"] = group
        selected.append(row)
        counts[group] += 1
        return True

    for row in eligible:
        if counts["score"] >= bands[0]:
            break
        take(row, "score")
    for row in eligible:
        if counts["diversity"] >= bands[1]:
            break
        if row["identity"].get("scaffold_id") not in cores:
            take(row, "diversity")
    uncertain = sorted(eligible, key=lambda r: (-(finite(r.get("priority_uncertainty")) or 0), *order(r)))
    for row in uncertain:
        if counts["uncertainty"] >= bands[2]:
            break
        associated = row.get("target_association") is True or bool(row.get("target_evidence")) or bool(row.get("measurements"))
        if associated and ((finite(row.get("priority_uncertainty")) or 0) > 0 or row.get("conflicting_evidence")):
            take(row, "uncertainty")
    deficits = {name: max(0, n-counts[name]) for name,n in zip(("score", "diversity", "uncertainty"), bands)}
    for row in eligible:
        take(row, "backfill")
    return {"provider": provider, "candidate_ids": [r["candidate_id"] for r in selected],
            "candidates": selected, "cap": cap, "category_deficits": deficits,
            "selected_counts": counts, "omitted": omitted, "population_preserved": True}


def cost_envelope(seconds_per_unit, *, units, remaining_seconds, measured_context_matches,
                  reserve=1.25):
    """Unknown or unmatched costs stay null; budget excludes queue wait."""
    if isinstance(units, bool) or not isinstance(units, int) or units < 0:
        raise ValueError("Work units must be a nonnegative integer")
    if finite(remaining_seconds) is None or remaining_seconds < 0 or finite(reserve) is None or reserve < 1:
        raise ValueError("Invalid resource budget")
    rate = finite(seconds_per_unit)
    if not measured_context_matches or rate is None or rate <= 0:
        return {"status": "benchmark_required", "compute_eta_seconds": None,
                "queue_eta_seconds": None, "admitted_units": 0, "deferred_units": units}
    allowed = min(units, math.floor(remaining_seconds / (rate * reserve)))
    return {"status": "ready" if allowed == units else "budget_deferred",
            "compute_eta_seconds": rate * units * reserve, "queue_eta_seconds": None,
            "admitted_units": allowed, "deferred_units": units-allowed}
