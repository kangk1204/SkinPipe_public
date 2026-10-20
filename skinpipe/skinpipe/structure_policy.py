"""Keep geometric controls, independent ranking validation and user weights separate."""
from __future__ import annotations

import math

VERSION = "structure-weight-gate-v1"


def _rmsd(value):
    if isinstance(value, (int, float)) and not isinstance(value, bool) \
            and math.isfinite(value) and value >= 0:
        return float(value)
    return None


def weights_for_run(weights, *, explicit, docked, redock, threshold):
    'No automatic structure bonus before an independent ranking evaluation.'
    applied = dict(weights or {})
    best, top = _rmsd(redock.get("best_pose_rmsd")), _rmsd(redock.get("rmsd"))
    def status(value):
        return "NOT_APPLICABLE" if not docked else "NOT_EVALUATED" if value is None else "PASS" if value <= threshold else "FAIL"
    search, pose_order = status(best), status(top)
    requested = applied.get("structure") if explicit else None
    valid_weight = isinstance(requested, (int, float)) and not isinstance(requested, bool) \
        and math.isfinite(requested) and requested >= 0
    if explicit and not valid_weight:
        raise ValueError("structure weight must be a finite non-negative number")
    policy = {"policy": VERSION, "legacy_overrides": "RETIRED_IN_SAMPLE_DIAGNOSTICS",
              "search_geometry_status": search, "top_pose_geometry_status": pose_order,
              "best_pose_rmsd": best, "top_pose_rmsd": top, "rmsd_threshold": threshold,
              "independent_ranking_status": "NOT_EVALUATED" if docked else "NOT_APPLICABLE",
              "requested_structure_weight": requested if valid_weight else repr(requested) if explicit else None,
              "explicit_weight": bool(explicit)}
    if not docked:
        applied["structure"] = 0.0
        reason = "Docking did not run; no structure contribution."
        application = "NOT_APPLICABLE"
    elif search == "FAIL":
        applied["structure"] = 0.0
        reason = f"Redock search failed: best sampled pose {best:.2f} A exceeds {threshold:g} A. This does not identify the cause of the search failure."
        application = "DISABLED_GEOMETRY"
    elif pose_order == "FAIL":
        applied["structure"] = 0.0
        reason = f"Redock pose ordering failed: top-ranked pose {top:.2f} A exceeds {threshold:g} A. A close sampled pose does not validate pose ranking."
        application = "DISABLED_GEOMETRY"
    elif explicit and valid_weight and requested == 0:
        applied["structure"] = 0.0
        reason = "Structure contribution explicitly disabled."
        application = "EXPLICITLY_DISABLED"
    elif explicit and valid_weight and search == pose_order == "PASS":
        reason = "Explicit experimental structure weight; independent ranking validation is missing. Geometric control success is not an affinity-ranking validation."
        application = "USER_DEFINED_UNVALIDATED"
    else:
        applied["structure"] = 0.0
        reason = "Structure contribution disabled: independent ranking validation is missing. Docking/MMGBSA values remain available as diagnostics; candidates are retained."
        application = "DISABLED_UNVALIDATED"
    policy.update(application=application, applied_structure_weight=applied["structure"], reason=reason)
    return applied, policy
