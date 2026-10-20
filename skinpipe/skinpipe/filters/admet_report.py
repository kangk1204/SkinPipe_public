"""Ranking and plain-language reporting for `filters.admet`.

Turning six domain dictionaries into one call a formulator can act on is its
own problem, and it is where a tox module usually goes wrong: it either buries
the finding in a wall of alerts, or it collapses everything into a single score
that hides which experiment to run next. Neither is useful at a bench.

So the contract here is: every concern keeps its mechanism and its next
experiment, concerns are ranked by severity and then by how much the domain
matters for a sun-exposed leave-on product, and the one-line overall verdict
names the binding constraint rather than averaging over it.
"""
from __future__ import annotations

from typing import Any

from ..types import ADMETConcern, ADMETProfile

# Reporting priority for screening concerns; a CMR alert is not a legal
# classification. Actual compound identity and exposure require separate review.
DOMAIN_ORDER = ("cmr", "phototoxicity", "genotoxicity", "irritation",
                "systemic", "cyp")


def collect(domain: str, block: dict[str, Any], concerns: list[ADMETConcern],
             fallback: str = "") -> None:
    """Lift a domain's alerts into ranked concerns, never losing the verdict.

    The fallback matters: a domain can reach a concern level with no named
    alert firing. Tretinoin is a phototoxicity concern because of its polyene,
    not a chemotype match, and a strong acid is an irritation concern because
    of its pKa class. Collecting only alerts would drop exactly the findings
    that come from reasoning rather than pattern-matching.
    """
    key = "chemotypes" if "chemotypes" in block else "alerts"
    before = len(concerns)
    for a in block.get(key, ()):
        concerns.append(ADMETConcern(domain=domain, name=a["name"],
                                     severity=a["severity"],
                                     mechanism=a.get("mechanism", ""),
                                     action=a.get("action", "")))
    level = block.get("concern", "none")
    if level in ("none", "unknown", "low"):
        return
    rank = ADMETProfile.SEVERITY_ORDER
    collected = [c.severity for c in concerns[before:]]
    if collected and min(rank.index(c) for c in collected) <= rank.index(level):
        return
    concerns.append(ADMETConcern(domain=domain, name=fallback or domain,
                                 severity=level, mechanism=block["verdict"],
                                 action=block.get("basis", "")))



def rank(concerns: list[ADMETConcern]) -> None:
    """Sort in place: severity first, then topical relevance of the domain."""
    order = {s: i for i, s in enumerate(ADMETProfile.SEVERITY_ORDER)}
    concerns.sort(key=lambda c: (order.get(c.severity, 99),
                                 DOMAIN_ORDER.index(c.domain)
                                 if c.domain in DOMAIN_ORDER else 99))


def overall(concerns: list[ADMETConcern], photo_block: dict,
             exposure: dict, pct: float, product: str) -> str:
    """Describe computed alerts and an exposure scenario without a use approval."""
    use = f"{pct:g}% in a {product.replace('_', ' ')}"
    scope = f"Scenario: {use}. Predictions only; safety and a suitable use level are not established. "
    if exposure.get("upper_bound_only"):
        sed_line = (f"Modeled SED upper bound {exposure['sed_mg_kg_day']:.3g} mg/kg bw/day "
                    "within the capped scenario; a human upper bound is not validated "
                    "and actual dermal absorption has not been measured")
    else:
        sed_line = (f"Modeled SED {exposure['sed_mg_kg_day']:.3g} mg/kg bw/day "
                    f"({exposure['band']}); actual dermal absorption has not been measured")
    concerns = [c for c in concerns if c.domain != "systemic"]
    blocking = [c for c in concerns if c.severity == "critical"]
    if blocking:
        return (scope + "Critical screening alerts: "
                + "; ".join(f"{c.name} ({c.domain})" for c in blocking)
                + ". Investigation: " + blocking[0].action + ". " + sed_line + ".")
    if photo_block["concern"] == "high":
        return (scope + "Predicted photosafety concern from structural features; "
                "experimental photosafety assessment is needed. " + sed_line + ".")
    high = [c for c in concerns if c.severity == "high"]
    if high:
        return (scope + "High screening alerts: "
                + "; ".join(f"{c.name} ({c.domain})" for c in high)
                + ". Investigation: " + high[0].action + ". " + sed_line + ".")
    moderate = [c for c in concerns if c.severity == "moderate"]
    if moderate:
        return (scope + "Moderate screening alerts to investigate: "
                + ", ".join(sorted({c.name for c in moderate})) + ". " + sed_line + ".")
    if exposure["band"] in ("material", "high"):
        return (scope + "No high-priority structural alert matched the evaluated rules. "
                + sed_line + ". A hypothetical margin-of-safety target of 100 implies "
                f"a NOAEL of {exposure['noael_required_for_mos_100']:.3g} mg/kg bw/day "
                "under this exposure model; this arithmetic does not establish a safety margin.")
    return (scope + "No high-priority structural alert matched the evaluated rules. "
            + sed_line + ". Absence of an alert is not an experimental safety finding.")
