"""Bounded, optional graph proposals with a complete exclusion census.

The source graph is never uncharged, desalted or replaced by a parent before
generation. A new graph is a proposal, not a binding or synthesis result.
"""

from collections import Counter

from .identity import chemical_identity, digest
from .ranking import finite


DEFAULT_GENERATION = {
    "enabled": False, "brics_products": 200, "max_3d": 50,
    "conformers": 10, "max_sa": 6.0,
}


def freeze_generation(body=None):
    body = {} if body is None else body
    if not isinstance(body, dict) or set(body) - set(DEFAULT_GENERATION):
        raise ValueError("Unknown generation setting")
    result = {**DEFAULT_GENERATION, **body}
    if not isinstance(result["enabled"], bool):
        raise ValueError("generation.enabled must be boolean")
    for key, lo, hi in (("brics_products", 0, 2000), ("max_3d", 0, 200), ("conformers", 1, 30)):
        value = result[key]
        if isinstance(value, bool) or not isinstance(value, int) or not lo <= value <= hi:
            raise ValueError(f"generation.{key} must be an integer within {lo}..{hi}")
    if finite(result["max_sa"]) is None or not 1 <= result["max_sa"] <= 10:
        raise ValueError("generation.max_sa must be within 1..10")
    return result


def generate_research(seed_smiles, policy, target_id, *, stop=lambda: None):
    from ..hop_generate import _bioisostere_candidates, _brics_candidates, BRICS_SEED
    from ..hop_sanity import is_sane, sa_assessment

    seed = chemical_identity(seed_smiles)
    receipt = {
        "policy": "research-exact-seed-generation-v1", "enabled": policy["enabled"],
        "original_seed_smiles": seed_smiles,
        "original_seed_exact_structure_id": seed["exact_structure_id"],
        "seed_standardization": "none", "settings": dict(policy),
        "rng_seed": BRICS_SEED, "proposals": [],
        "interpretation": "Graph proposals and exploratory scaffold comparisons; no activity inheritance, synthesis route, patent novelty or freedom-to-operate determination",
    }
    if not policy["enabled"]:
        return [], {**receipt, "raw_proposal_count": 0, "unique_structure_count": 0, "exclusion_counts": {}}
    if seed.get("enhanced_stereo"):
        return [], {**receipt, "status": "held", "reason": "enhanced_stereochemistry_requires_explicit_isomer",
                    "raw_proposal_count": 0, "unique_structure_count": 0, "exclusion_counts": {}}
    if seed["structure_status"] != "valid" or seed["parent_ambiguous"]:
        raise ValueError("Generation requires one valid exact source structure")
    stop()
    raw = _bioisostere_candidates(seed_smiles)
    stop()
    raw += _brics_candidates(seed_smiles, policy["brics_products"], None, 2, BRICS_SEED)
    rows, seen, exclusions = [], {}, Counter()
    for ordinal, (smiles, derivation) in enumerate(raw, 1):
        stop()
        identity = chemical_identity(smiles)
        exact = identity["exact_structure_id"]
        reasons = []
        sane, sanity_reason = is_sane(smiles)
        if not sane:
            reasons.append("generation_policy_" + sanity_reason)
        assessment = sa_assessment(smiles)
        if assessment["score"] is None:
            reasons.append("generation_sa_not_assessed")
        elif assessment["score"] > policy["max_sa"]:
            reasons.append("generation_sa_above_threshold")
        if exact == seed["exact_structure_id"]:
            reasons.append("regenerated_query_structure")
        if exact is not None and exact in seen:
            reasons.append("duplicate_generated_structure")
        proposal = {
            "ordinal": ordinal, "raw_smiles": smiles,
            "canonical_smiles": identity["identity_smiles"],
            "identity": identity, "derivation": derivation,
            "sa_assessment": assessment, "exclusion_reasons": reasons,
            "accepted_by_generation_policy": not reasons,
        }
        receipt["proposals"].append(proposal)
        exclusions.update(reasons)
        if exact is None or "regenerated_query_structure" in reasons or "duplicate_generated_structure" in reasons:
            continue
        cid = "gen-" + exact.split(":", 1)[1][:32]
        seen[exact] = cid
        proposal["candidate_id"] = cid
        rows.append({
            "candidate_id": cid, "target_id": target_id, "smiles": identity["identity_smiles"],
            "identity": identity, "name": cid, "role": "candidate",
            "derivation": {"method": derivation, "seed_exact_structure_id": seed["exact_structure_id"], "generation_policy": receipt["policy"]},
            "generation_filter_reasons": list(reasons),
            "generation_proposal_ordinals": [ordinal], "purchasable": False,
            "admet": [], "admet_status": "unassessed_model_not_configured",
            "admet_recalculation_required": True,
            "measurements": [], "docking": [], "target_association": False,
            "generation_origin_sha256": digest([seed["exact_structure_id"], derivation, ordinal, exact]),
        })
    return rows, {**receipt, "raw_proposal_count": len(raw),
                  "unique_structure_count": len(rows),
                  "generation_policy_pass": sum(not r["generation_filter_reasons"] for r in rows),
                  "exclusion_counts": dict(sorted(exclusions.items()))}
