"""Research visibility, compute, interpretation, and experiment states are independent."""
from __future__ import annotations

from copy import deepcopy
from types import MappingProxyType
from typing import Any

from .identity import candidate_identity
from .procurement import has_purchase_record
from .ranking import finite
from .algorithm_profiles import compute_profile

POLICY_VERSION = "research-platform-v2.5"
_DEFAULT_PROFILE = compute_profile()
COMPUTE_DEFAULTS = MappingProxyType({
    **_DEFAULT_PROFILE["counts"],
    "optional_advanced": _DEFAULT_PROFILE["counts"]["complex_models"],
    "target_compute_seconds": 28_800, "canary_max_candidates": 32,
    "reserve_vram_mib": 1536, "peak_memory_multiplier": 1.25,
    "cpu_fraction": 0.25, "cpu_max": 8, "seed": 42,
})


def assess_candidate(candidate: dict[str, Any], *, target_identity_valid: bool = True,
                     engine_supported: bool = True, trust_stored_identity: bool = False) -> dict[str, Any]:
    """Return a new candidate; absent evidence never erases a source row."""
    out = deepcopy(candidate)
    identity = candidate_identity(candidate, trust_stored=trust_stored_identity)
    warnings: list[str] = []
    reasons: list[str] = []
    if identity["structure_status"] != "valid":
        reasons.append(identity["structure_status"])
    if candidate.get("identity_conflict"):
        reasons.append("conflicting_structure_identity")
    if identity.get("parent_ambiguous"):
        reasons.append("multicomponent_material_requires_preparation")
    if target_identity_valid is not True or candidate.get("target_identity_valid") is False:
        reasons.append("unresolved_target_identity")
    if engine_supported is not True or candidate.get("engine_supported") is False:
        reasons.append("engine_unsupported_chemistry")
    if candidate.get("receptor_available") is False:
        reasons.append("receptor_unavailable")
    if not candidate.get("admet"):
        warnings.append("admet_not_evaluated")
    purchase_available = has_purchase_record(candidate)
    if not purchase_available:
        warnings.append("purchase_availability_unconfirmed")
    skin = candidate.get("skin_score")
    if skin is None:
        skin = candidate.get("topical_score")
    if finite(skin) is None:
        warnings.append("skin_score_not_evaluated")
    elif skin < 4:
        warnings.append("low_skin_score")
    if identity.get("molecular_weight", 0) > 800:
        warnings.append("high_molecular_weight")
    if identity.get("unassigned_stereo") or identity.get("unassigned_bond_stereo"):
        warnings.append("stereochemistry_unassigned")
    if identity.get("fragment_count", 1) > 1:
        warnings.append("salt_or_multicomponent_material")
    if candidate.get("redock_validated") is not True:
        warnings.append("redocking_not_validated")
    if candidate.get("critical_warning") or candidate.get("hazard_severity") == "critical":
        warnings.append("critical_hazard_review")
    role = candidate.get("role", "candidate")
    if candidate.get("opposing_direction") is True or role == "counter_control":
        role = "counter_control"
        warnings.append("opposing_direction_control")
    limitations = [name for name in ("covalent", "metal_dependent", "membrane_context")
                   if candidate.get(name)]
    warnings.extend("special_interpretation:" + x for x in limitations)
    artifact = candidate.get("artifact_valid")
    interpretation = "invalid_artifact" if artifact is False else (
        "context_limited" if limitations else "provisional")
    if artifact is True and candidate.get("protocol_validated") is True and not limitations:
        interpretation = "validated_protocol_scope"
    out.update(identity=identity, role=role, platform_policy=POLICY_VERSION,
        states={
            "discovery": {"visible": True, "status": "control" if role != "candidate" else "candidate"},
            "compute": {"allowed": not reasons, "status": "held" if reasons else "eligible",
                        "reasons": reasons, "scope": "docking",
                        "evidence_route_allowed": True,
                        "ligand_route_allowed": identity["structure_status"] == "valid"},
            "interpretation": {"status": interpretation, "limitations": limitations,
                               "efficacy_established": False},
            "experiment": {"status": "needs_review" if reasons or warnings else "ready_to_plan",
                           "purchasable": True if purchase_available else candidate.get("purchasable"),
                           "controls_ready": candidate.get("controls_ready", False)}},
        research_warnings=warnings)
    return out


def compute_budget(stage: str, requested: int | None = None) -> int:
    'Validate the selected structure and retain explicit calculation provenance.'
    stages = {"cheap_rerank", "broad_docking", "confirmation", "receptor_sensitivity",
              "automatic_deep", "optional_mmgbsa", "optional_advanced"}
    if stage not in stages:
        raise ValueError("unknown compute stage")
    limit = COMPUTE_DEFAULTS[stage]
    if requested is None:
        return limit
    if isinstance(requested, bool) or not isinstance(requested, int) or requested < 0:
        raise ValueError("requested count must be a nonnegative integer")
    return min(requested, limit)
