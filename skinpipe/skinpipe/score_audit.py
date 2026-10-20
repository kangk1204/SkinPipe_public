"""Independent reconstruction of exported scoring inputs; never writes results.

This module intentionally does not call scorer functions. Missing historical
inputs produce NOT_EVALUATED instead of substituting zeros and declaring PASS.
"""
from __future__ import annotations
import math
from fractions import Fraction

SEVERITY = {"high": 3.0, "moderate": 1.5, "low": 0.5, "context": 0.5}
COLOUR = {"low": 0.0, "moderate (may tint)": 1.0,
          "strong (expect visible colour)": 2.5, "unknown": 0.0}
# Independent copy of the exported protocol, not an import from the scorer.
# The evidence class itself still requires a separate raw-biological audit.
EVIDENCE = {"clinical": 1.0, "ex_vivo_human": .95, "human_cell": .9,
            "biochemical": .8, "nonhuman_cell": .65, "surrogate": .6,
            "animal": .65, "predicted": .55, "unknown": .5}


def finite(v):
    if not isinstance(v, (int, float)) or isinstance(v, bool):
        return False
    try:
        return math.isfinite(v)
    except OverflowError:
        return False


def alerts(v):
    return [s for s in v.split("|") if s] if isinstance(v, str) else list(v or [])


def structure_terms(rows, resolution=2.0):
    energies = [r.get("dock_affinity") for r in rows]
    present = [v for v in energies if finite(v)]
    out = []
    for val in energies:
        if not finite(val):
            out.append(0.0)
            continue
        better = sum(v < val for v in present)
        tied = sum(v == val for v in present)
        out.append(round(10 * (len(present) - better - (tied - 1) / 2) / len(present), 2))
    subset = [(i, r["mmgbsa_dG"]) for i, r in enumerate(rows) if finite(r.get("mmgbsa_dG"))]
    slots = sorted((out[i] for i, _ in subset), reverse=True)
    ordered = sorted(subset, key=lambda item: item[1])
    position = 0
    while position < len(ordered):
        end = position + 1
        # Exported protocol: adjacent unresolved differences form a transitive
        # tier. A chain can span more than resolution. This independent loop
        # deliberately does not import the scorer's tier implementation.
        while end < len(ordered) and ordered[end][1] - ordered[end - 1][1] <= resolution:
            end += 1
        # The frozen term protocol sums descending slots with binary float
        # arithmetic before rounding. statistics.mean changes that arithmetic
        # at half-cent boundaries (10, 8.33, 5, 3.33 gives 6.67 instead of the
        # recorded 6.66). Composite v4's exact arithmetic begins after these
        # already-rounded terms; it does not redefine their generation.
        binary_total = 0.0
        for slot in slots[position:end]:
            binary_total += slot
        share = round(binary_total / (end - position), 2)
        for i, _ in ordered[position:end]:
            out[i] = share
        position = end
    return out


def audit(rows):
    struct = structure_terms(rows)
    output = []
    for i, row in enumerate(rows):
        required = ("skin_assessed", "n_pains", "weights_used", "weighted_terms",
                    "weight_denominator", "evidence_factor", "score_terms", "scoring_policy", "evidence")
        current = row.get("scoring_policy") == "screen-score-v5"
        if current:
            required = tuple(k for k in required if k != "n_pains") + (
                "sensitisation_status", "stability_status")
        missing = [k for k in required if k not in row or row[k] is None]
        if current:
            for field in ("sensitisation_status", "stability_status"):
                if row.get(field) not in ("EVALUATED", "NOT_EVALUATED"):
                    missing.append("supported " + field)
            if row.get("sensitisation_status") == "EVALUATED" and not finite(row.get("n_pains")):
                missing.append("finite n_pains")
        if row.get("evidence") not in EVIDENCE:
            missing.append("supported evidence class")
        if row.get("structure_tier_policy", "adjacent-gap-transitive-v1") != "adjacent-gap-transitive-v1":
            missing.append("supported structure tier policy")
        if row.get("scoring_policy") not in ("screen-score-v2", "screen-score-v3", "screen-score-v4", "screen-score-v5"):
            missing.append("supported scoring policy")
        if row.get("scoring_policy") in ("screen-score-v4", "screen-score-v5") and row.get("composite_rounding_policy") != "decimal-half-even-2-v1":
            missing.append("supported composite rounding policy")
        if missing:
            output.append({"ligand_id": row.get("ligand_id"), "status": "NOT_EVALUATED",
                           "missing_inputs": missing})
            continue
        skin = row["skin_assessed"]
        pain = row.get("n_pains")
        if skin and not finite(pain):
            output.append({"ligand_id": row.get("ligand_id"), "status": "NOT_EVALUATED",
                           "missing_inputs": ["finite n_pains"]})
            continue
        p = row.get("pchembl")
        burden = sum(SEVERITY.get(a.split(":")[0], 0) for a in alerts(row.get("sensitiser_alerts"))
                     if not a.split(":", 1)[-1].startswith("PAINS("))
        rebuilt = {"potency": round(max(0, min(10, (p - 4) * 2)), 2) if finite(p) else 0.0,
                   "structure": struct[i], "topical": row.get("topical_score") if skin else 0.0,
                   "safety": round(max(0, 10 - 2 * burden - 1.5 * pain), 2) if skin else 5.0,
                   "stability": round(max(0, 10 - 1.5 * len(alerts(row.get("stability_alerts")))
                                          - COLOUR.get(row.get("colour") if skin else "unknown", 2.5)), 2)}
        weights, included = row["weights_used"], row["weighted_terms"]
        factor = EVIDENCE[row["evidence"]]
        if current:
            # v5 records each filter's missingness independently. A neutral
            # placeholder is reconstructed as such, never inferred from an
            # empty alert list. Topical assessment may exist on a partial row.
            rebuilt["topical"] = row.get("topical_score") if finite(row.get("topical_score")) else 0.0
            rebuilt["safety"] = (round(max(0, 10 - 2 * burden - 1.5 * pain), 2)
                                 if row["sensitisation_status"] == "EVALUATED" else 5.0)
            rebuilt["stability"] = (
                5.0 if row["stability_status"] != "EVALUATED" or row.get("colour") == "unknown"
                else round(max(0, 10 - 1.5 * len(alerts(row.get("stability_alerts")))
                               - COLOUR.get(row.get("colour"), 2.5)), 2))
        if row["scoring_policy"] in ("screen-score-v4", "screen-score-v5"):
            rational_weights = {k: Fraction(str(weights[k])) for k in included}
            exact_denominator = sum(rational_weights.values(), Fraction()) or Fraction(1)
            exact_mean = sum((rational_weights[k] * Fraction(str(rebuilt[k])) for k in included), Fraction())
            scaled = exact_mean * Fraction(str(factor)) / exact_denominator * 100
            nearest = scaled.numerator // scaled.denominator
            distance = scaled - nearest
            if distance > Fraction(1, 2) or (distance == Fraction(1, 2) and nearest % 2 != 0):
                nearest += 1
            value, denominator = nearest / 100.0, float(exact_denominator)
        else:
            # Historical v2/v3 used binary arithmetic in this insertion order.
            # Their sorted weighted_terms field did not preserve calculation
            # order. Reproduce the documented version, never widen tolerance.
            legacy_order = [k for k in ("topical", "safety", "stability", "potency", "structure") if k in included]
            denominator, numerator = 0.0, 0.0
            for k in legacy_order:
                denominator += weights[k]
                numerator += weights[k] * rebuilt[k]
            denominator = denominator or 1.0
            value = round(numerator / denominator * factor, 2)
        mismatch = {k: {"rebuilt": v, "reported": row["score_terms"].get(k)}
                    for k, v in rebuilt.items() if v != row["score_terms"].get(k)}
        if value != row.get("composite"):
            mismatch["composite"] = {"rebuilt": value, "reported": row.get("composite")}
        if abs(denominator - row["weight_denominator"]) > 1e-12:
            mismatch["weight_denominator"] = {"rebuilt": denominator,
                                             "reported": row["weight_denominator"]}
        if factor != row["evidence_factor"]:
            mismatch["evidence_factor"] = {"rebuilt": factor, "reported": row["evidence_factor"]}
        if "composite_without_evidence_factor" in row:
            if row.get("evidence_diagnostic_policy") != "pre-evidence-weighted-mean-v1":
                output.append({"ligand_id": row.get("ligand_id"), "status": "NOT_EVALUATED",
                               "missing_inputs": ["supported evidence diagnostic policy"]})
                continue
            diagnostic_weights = {key: Fraction(str(weights[key])) for key in included}
            diagnostic_mean = sum((diagnostic_weights[key] * Fraction(str(rebuilt[key]))
                                   for key in included), Fraction()) / (
                                       sum(diagnostic_weights.values(), Fraction()) or Fraction(1))
            # Fraction.__round__ implements exact half-even; no scorer call.
            unmultiplied = round(diagnostic_mean * 100) / 100.0
            for key, expected in (("composite_without_evidence_factor", unmultiplied),
                                  ("evidence_factor_contribution", round(value - unmultiplied, 2))):
                if row.get(key) != expected:
                    mismatch[key] = {"rebuilt": expected, "reported": row.get(key)}
        output.append({"ligand_id": row.get("ligand_id"), "status": "FAIL" if mismatch else "PASS",
                       "mismatches": mismatch})
    return {"status": "FAIL" if any(r["status"] == "FAIL" for r in output)
            else "NOT_EVALUATED" if any(r["status"] == "NOT_EVALUATED" for r in output)
            else "PASS" if output else "NOT_EVALUATED", "rows": output,
            "scope": "Arithmetic and evidence multiplier reconstructed from exported class; raw biological classification is not independently validated"}
