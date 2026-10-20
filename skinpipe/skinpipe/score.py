"""Heuristic composite scoring for topical candidates.

These weights have not been calibrated against independent efficacy data.
The pipeline disables the structure contribution by default until independent
ranking validation exists; standalone weights below support explicit diagnostic
comparisons and reconstruction of historical results.

The values below are an exploratory preference profile, not measured relative
efficacy, safety, permeability or formulation-success probabilities. Their
ordering and size are not an experimentally calibrated optimum.

    topical    0.30   can it get through the stratum corneum at all
    potency    0.25   exact primary potency; 0 when no exact value is available
    safety     0.20   sensitisation and frequent-hitter alerts
    structure  0.15   docking score and MM/GBSA rescoring
    stability  0.10   survives formulation, packaging and shelf life

Potency and structure are separate terms and neither substitutes for the other.
A missing exact primary potency scores 0 even if a censored or qualitative
observation exists. Its observation status is exported separately. A pose does
not establish a measured potency.

Missing terms
-------------
A term is dropped from the weighted mean only when the *stage* never ran for
anyone in the screen -- no docking, no MM/GBSA. When the stage did run and this
one ligand has no result, the term scores 0 and stays in the denominator.
Excusing it instead rewards failure: a ligand Vina could not prepare would
outrank the same ligand docked badly, and ligand-prep failures are silent.

Change them with a `weights` dict when the product format warrants it: a
rinse-off cleanser can tolerate a poor permeant, an anti-ageing serum cannot.

Non-finite numbers
------------------
NaN and +/-Inf are missing, never scored. Every number here arrives from
outside the module -- pChEMBL from the ChEMBL API, affinity from Vina, dG_bind
from an MM/GBSA worker, weights from config -- and CPython lets a corrupt one
through both checks this module used to rely on: ``min(10.0, nan)`` is ``10.0``,
so the clamp hands NaN the best score on the 0-10 scale, and ``nan is not None``
is ``True``, so the "was anything measured" test passes for a value that
measured nothing. A NaN pChEMBL therefore scored a perfect potency 10, outranked
a real 7.0, and printed no caveat. `_finite` is the single gate; a rejected value
is named in a caveat so the CSV says which record was unusable and why.

On the structure term
---------------------
Vina scores and MM/GBSA energies are not comparable across targets, and
MM/GBSA single-trajectory numbers are not binding free energies at all. Both
are therefore normalised *within a screen* by percentile rank rather than
against absolute anchors. Percentiles are a descriptive convention; they do
not validate affinity ordering or transfer between candidate populations.

MM/GBSA is also affordable only on a subset -- `pipeline.run_screen` scores with
`use_gbsa=False`, rescores the leaders of that pass, then scores again -- so one
ranked column ends up fed by two kinds of evidence. An earlier build averaged
them, putting two reference populations of different size inside one number:
with six of twelve rescored, c02 (Vina -10.9, rescored) came out at structure
4.33 and c03 (Vina -10.5, never rescored) at 7.50, because inside a six-member
pool somebody comes last: selection was itself a demotion. Instead, docking
sets the scale for the whole screen and MM/GBSA reorders the rescored group
within the scores that group already holds, redistributing tied slots by their
rounded mean, so
the column carries one scale. `_rescored_structure` works an example through,
and `ScoreBreakdown.structure_evidence` names the method behind each row.
"""
from __future__ import annotations

import json
import logging
from pathlib import Path
import math
import random
from fractions import Fraction
from dataclasses import dataclass, field
from typing import Optional, Sequence

from .filters.safety import SEVERITY_WEIGHT
from .types import Candidate, Evidence

log = logging.getLogger(__name__)

DEFAULT_WEIGHTS: dict[str, float] = {
    "topical": 0.30, "potency": 0.25, "safety": 0.20,
    "structure": 0.15, "stability": 0.10,
}

#: Exploratory tie threshold, informed by a historical kinase-complex spread
#: (+/- 2.0 kcal/mol over 5 snapshots). This is a ranking convention, not an
#: established error bound or target-independent MM/GBSA resolution.
GBSA_RESOLUTION_KCAL: float = 2.0

#: What set a row's structure term: (rationale tag, report label). On every
#: ScoreBreakdown, so a part-rescored screen names its rescored rows.
STRUCTURE_EVIDENCE: dict[str, tuple[str, str]] = {
    "mmgbsa": ("MM/GBSA", "docked, then MM/GBSA rescored"),
    "dock": ("dock", "docked only"),
    "none": ("no pose", "docking ran, no pose for this compound"),
    "budget_skipped": ("not docked", "outside the --dock-top budget; never "
                       "offered to the engine"),
    "not_run": ("not run", "docking did not run for this screen"),
}

#: Structure-evidence keys that mean the stage answered for this compound.
#: "none" belongs here: docking ran, returned nothing, and that is a result.
EVALUATED: frozenset[str] = frozenset({"mmgbsa", "dock", "none"})

#: Every source label that means "this row came out of the ChEMBL activity
#: table", and is therefore a library comparator rather than a control.
#: `library.from_chembl` mints two: ``chembl`` for a row fetched whole from the
#: API, and ``chembl_local`` (finding N-20) for one whose *structure* had to be
#: rescued from the local chemreps archive because the API was unreachable. Only
#: the structure was ever missing -- the pChEMBL comes from the activity table,
#: which is a separate call from ``/molecule`` -- so a rescued row carries just
#: as real a measured potency as any other and belongs in the same pool.
#:
#: This is a named constant and not a literal at the point of use because the
#: literal is exactly what went wrong. N-20 added the second label and widened
#: `library.assemble`'s de-duplication precedence for it; `control_check`'s
#: comparator pool kept its bare ``== "chembl"`` and nothing tied the two
#: together, so the pool silently emptied during precisely the outages the local
#: fallback exists to survive. Read this set; do not spell the labels out again.
#:
#: DEPENDENCY: `library.assemble` uses a wider source set, and it is NOT the same set: it
#: adds `reference_offnode`, so it is a strict superset of this one. Importing
#: this constant there would silently narrow that comparison. If the two ever
#: need to agree, the wider set needs its own name -- do not merge them.
#: The wider source priority remains local to `library.assemble`; the sets
#: serve different comparisons and are tested independently.
CHEMBL_LIBRARY_SOURCES: frozenset[str] = frozenset({"chembl", "chembl_local"})

#: Evidence multiplier. A measured human number outranks a docking pose.
EVIDENCE_FACTOR: dict[Evidence, float] = {
    Evidence.CLINICAL: 1.00,
    Evidence.EX_VIVO_HUMAN: 0.95,
    Evidence.HUMAN_CELL: 0.90,
    Evidence.BIOCHEMICAL: 0.80,
    Evidence.NONHUMAN_CELL: 0.65,
    # Below a non-human cell line because the protein itself is a different
    # one, above `PREDICTED` because a measurement was still made.
    Evidence.SURROGATE: 0.60,
    Evidence.ANIMAL: 0.65,
    Evidence.PREDICTED: 0.55,
    Evidence.UNKNOWN: 0.50,
}


#: `pipeline.assess_skin` fills `SkinAssessment.sensitiser_alerts` with
#: ``f"{severity}:{name}"`` entries from `filters.safety.sensitisation_alerts`,
#: and then appends one synthetic entry, ``f"low:PAINS({n})"``, carrying the
#: frequent-hitter count. This is the name that marks the second kind.
PAINS_ALERT_NAME: str = "PAINS"


def _electrophilic_alerts(alerts: Optional[Sequence[str]]) -> list[str]:
    'The genuine sensitisation alerts in a `sensitiser_alerts` list.'
    out = []
    for a in alerts or []:
        name = a.split(":", 1)[1] if ":" in a else a
        if not name.startswith(PAINS_ALERT_NAME):
            out.append(a)
    return out


def _finite(value: object) -> Optional[float]:
    """``value`` as a finite float, or ``None`` when it is not a usable number.

    ``None``, NaN, +/-Inf, booleans (a YAML ``pchembl: yes`` is not a potency of
    1.0) and anything ``float()`` refuses all come back as ``None``, so a caller
    that asks ``is not None`` gets the answer it meant to ask for.
    """
    if value is None or isinstance(value, bool):
        return None
    try:
        f = float(value)
    except (TypeError, ValueError, OverflowError):
        return None
    return f if math.isfinite(f) else None


def _nonfinite_note(field_name: str, value: object) -> str:
    'The caveat for a rejected number. It names the value; a reader of the'
    return (f"recorded {field_name} {value!r} is not a finite number; it is treated "
            "as missing and scores nothing, because a non-finite value defeats the "
            "0-10 clamp and would take the top of the scale")


def _nonnegative_count(value: object) -> Optional[float]:
    """A finite cardinal count, or None when it cannot describe alert matches."""
    count = _finite(value)
    return count if count is not None and count >= 0 and count.is_integer() else None


@dataclass
class ScoreBreakdown:
    ligand_id: str
    potency: float = 0.0
    structure: float = 0.0
    topical: float = 0.0
    safety: float = 0.0
    stability: float = 0.0
    evidence_factor: float = 1.0
    composite: float = 0.0
    #: Diagnostic on the same terms/denominator before the historical global
    #: evidence multiplier. Neither field changes the versioned ranking.
    composite_without_evidence_factor: float = 0.0
    evidence_factor_contribution: float = 0.0
    evidence_diagnostic_policy: str = "pre-evidence-weighted-mean-v1"
    structure_tier_policy: str = "adjacent-gap-transitive-v1"
    docking_tie_tolerance_kcal: float = 0.0
    structure_missing_policy: str = "screen-v1-no-pose-zero; never-run-excluded"
    rationale: str = ""
    caveats: list[str] = field(default_factory=list)
    structure_evidence: str = "not_run"      # key of STRUCTURE_EVIDENCE
    #: Which block of the ranking this row sits in. "evaluated" when the
    #: structure stage answered for it -- including answering "no pose" --
    #: and "not_evaluated" when it was never offered to the engine. A screen
    #: with only one kind calls everything "evaluated", because then the
    #: column really is one scale.
    rank_group: str = "evaluated"
    rank_within_group: int = 0
    rank_group_size: int = 0
    rank_interpretation: str = "composite order within shared scored terms"
    weights_used: dict[str, float] = field(default_factory=dict)
    weighted_terms: list[str] = field(default_factory=list)
    weight_denominator: float = 0.0
    scoring_policy: str = "screen-score-v5"
    composite_rounding_policy: str = "decimal-half-even-2-v1"


def _tolerance_tiers(ordered: Sequence[tuple[int, float]],
                     tolerance: float) -> list[list[int]]:
    """Group an ordered (index, value) list into tiers the method cannot split.

    ``tolerance`` is the smallest difference the method can actually resolve.
    Grid-rounding does not achieve it -- two values 0.4 apart can straddle a bin
    edge and land in different bins. Every adjacent pair within tolerance is
    connected into the same tier. This is a transitive uncertainty grouping:
    a chain of unresolved differences can span more than the tolerance. It
    deliberately avoids assigning different ranks to any unresolved pair.
    """
    tiers: list[list[int]] = []
    previous: Optional[float] = None
    for idx, val in ordered:
        if previous is None or abs(val - previous) > tolerance:
            tiers.append([])
        tiers[-1].append(idx)
        previous = val
    return tiers


def _percentile_scores(values: Sequence[Optional[float]], lower_is_better: bool,
                       tolerance: float = 0.0) -> list[float]:
    """Map values to 0-10 by rank within the screen. ``None`` scores 0.

    Values within ``tolerance`` share a score; `_tolerance_tiers` groups them.

    A non-finite value scores 0 like a missing one, and is dropped *before* the
    sort rather than ranked: NaN compares False against everything, so `sorted`
    returns an arbitrary order and `_tolerance_tiers` -- whose test is
    ``abs(val - anchor) > tolerance`` -- opens no new tier, collapsing the whole
    screen into one tie. One corrupt Vina score used to flatten the structure
    column for every other compound in the run.
    """
    present = [(i, f) for i, f in ((j, _finite(v)) for j, v in enumerate(values))
               if f is not None]
    out = [0.0] * len(values)
    if not present:
        return out
    ordered = sorted(present, key=lambda t: t[1], reverse=not lower_is_better)
    n = len(ordered)
    placed = 0
    for tier in _tolerance_tiers(ordered, tolerance):
        # every member of a tier gets the score of the tier's midpoint rank
        mid = placed + (len(tier) - 1) / 2.0
        score = 10.0 * (n - mid) / n if n > 1 else 10.0
        for idx in tier:
            out[idx] = round(score, 2)
        placed += len(tier)
    return out


def _rescored_structure(dock_s: Sequence[float],
                        gbsa_vals: Sequence[Optional[float]],
                        tolerance: float) -> list[float]:
    """Let MM/GBSA reorder the rescored candidates inside their own slots.

    The rescored group is redistributed across the structure scores docking
    gave that group and nothing else, so MM/GBSA decides who holds which --
    the only comparison it may make, having run here.

    What remains on the docking scale
    ---------------------------------
    Tied candidates share the mean of their group's slots, rounded per row to
    two decimals. Their rounded sum can therefore differ slightly: slots
    10.00, 8.33, 5.00 and 3.33 sum to 26.66, but four tied scores of 6.66 sum
    to 26.64. Every candidate outside the rescored group keeps its own docking
    score untouched, and no score from outside that group is handed in.

    The multiset survives only when MM/GBSA separates every member: then each
    one takes a whole slot and the same values come back permuted. As soon as
    two candidates sit within `tolerance` they tie, and a tie is *averaged* --
    they share the mean of the slots they jointly hold, which spends both slots
    and returns neither. ``_rescored_structure([10.0, 5.0], [-10.0, -9.0], 2)``
    is ``[7.5, 7.5]``: 15.0 in, 15.0 out, and the 10.0 and the 5.0 are gone.
    That is the intended exploratory tie policy: the configured threshold
    treats that 1 kcal/mol difference as unresolved. It is not multiset
    preservation and this docstring used to say it was.

    Worked example, the six rescored of a twelve-compound screen. Docking gives
    that six the slots 10.00, 8.33, 6.67, 5.83, 4.17, 1.67; MM/GBSA reports c07
    -55.0, c00 -52.0, c04 -49.5, c05 -49.0, c10 -47.0, c02 -41.0. c04/c05 and
    c05/c10 each sit within 2 kcal/mol, so all three form one connected tier.
    They consume the third through fifth slots: each gets their mean,
    (6.67 + 5.83 + 4.17) / 3 = 5.56; the rest take 10.00, 8.33, 1.67 in MM/GBSA
    order, so c07 rises from 4.17 on docking
    alone to the best slot the group owns and c02 falls from 8.33 to the worst.
    Had it agreed with docking, no rank would have moved.
    """
    finite = [(i, f) for i, f in ((j, _finite(v)) for j, v in enumerate(gbsa_vals))
              if f is not None]
    tier = [i for i, _ in finite]
    out = list(dock_s)
    if not tier:
        return out
    slots = sorted((dock_s[i] for i in tier), reverse=True)
    ordered = sorted(finite, key=lambda t: t[1])
    placed = 0
    for group in _tolerance_tiers(ordered, tolerance):
        # This term's frozen protocol predates Python 3.12's compensated sum.
        # Retain ordered binary additions across runtimes; composite v4 applies
        # its separately versioned exact arithmetic after term generation.
        total = 0.0
        for slot in slots[placed:placed + len(group)]:
            total += slot
        share = round(total / len(group), 2)
        for i in group:
            out[i] = share
        placed += len(group)
    return out


#: Where the measured structure-weight overrides live, relative to the package.
STRUCTURE_OVERRIDES_PATH = (Path(__file__).resolve().parent.parent / "config"
                            / "structure_weight_overrides.json")


def structure_override(node: Optional[str]) -> Optional[dict]:
    """No project-specific structure weights are included in this distribution."""
    return None


def potency_score(pchembl: Optional[float]) -> float:
    """pChEMBL 6 -> 4, 7 -> 6, 8 -> 8, >=9 -> 10. Unmeasured -> 0.

    Non-finite is unmeasured. ``min(10.0, nan)`` is ``10.0``, so before this was
    guarded a NaN pChEMBL scored a perfect 10 -- the best potency in the screen,
    from a record that measured nothing.
    """
    p = _finite(pchembl)
    if p is None:
        return 0.0
    return round(max(0.0, min(10.0, (p - 4.0) * 2.0)), 2)


def safety_score(sensitisation_burden: float, n_pains: int = 0) -> float:
    """10 is clean. Burden and PAINS both subtract.

    An unreadable burden or PAINS count is an *unknown* penalty, not an absent
    one, so it scores 0 rather than leaving the compound looking clean.
    """
    burden, pains = _finite(sensitisation_burden), _nonnegative_count(n_pains)
    if burden is None or burden < 0 or pains is None:
        return 0.0
    return round(max(0.0, 10.0 - 2.0 * burden - 1.5 * pains), 2)


#: Colour penalties apply only to evaluated structures. An unknown judgement
#: uses the neutral placeholder below and is explicitly marked unmeasured.
_COLOUR_PENALTY = {
    "low": 0.0,
    "moderate (may tint)": 1.0,
    "strong (expect visible colour)": 2.5,
}


def skin_filter_status(candidate, domain: str) -> str:
    """Require recorded execution evidence and an assessable current structure.

    A parseable SMILES proves only that a filter *could* run. Legacy empty
    alerts/statuses cannot distinguish a clean result from a skipped filter.
    """
    if domain not in ("sensitisation", "stability"):
        raise ValueError(f"unknown skin filter domain: {domain}")
    if candidate.skin is None:
        return "NOT_EVALUATED"
    recorded = getattr(candidate.skin, domain + "_status", None)
    if recorded != "EVALUATED":
        return "NOT_EVALUATED"
    from .filters.safety import assessable_molecule
    return ("EVALUATED" if assessable_molecule(candidate.ligand.smiles) is not None
            else "NOT_EVALUATED")


def colour_of(candidate) -> str:
    """The colour judgement that scores this candidate's stability.

    Prefers the call `pipeline.assess_skin` already made and recorded on the
    assessment, so the sentence in the report and the number in the score are
    the same judgement. Falls back to `stability.colour_risk` on the SMILES for
    an assessment built by some other path -- a hand-built one in a test, an
    older cached screen -- because a colour term that only works when one
    particular caller populated one particular attribute is a term that quietly
    does nothing, which is how this repository keeps shipping inert guards.
    """
    if candidate.skin is None:
        return "unknown"
    recorded = getattr(candidate.skin, "colour", None)
    if isinstance(recorded, str) and recorded:
        return recorded
    smiles = getattr(candidate.ligand, "smiles", "") or ""
    if not smiles:
        return "unknown"
    try:
        from .filters import stability as _stability
        return _stability.colour_risk(smiles)
    except Exception:                      # never let a colour call stop a screen
        return "unknown"


def stability_score(n_alerts: int, colour: str = "low") -> float:
    """10 survives the jar. An unreadable alert count scores 0, not 10.

    `colour` is a `stability.colour_risk` verdict; see `_COLOUR_PENALTY` for
    what each evaluated call costs. "unknown" is a neutral placeholder of 5,
    not a clean 10 or a measured adverse finding.
    """
    alerts = _nonnegative_count(n_alerts)
    if alerts is None:
        return 0.0
    if colour == "unknown":
        return 5.0
    penalty = 1.5 * alerts + _COLOUR_PENALTY.get(colour, 2.5)
    return round(max(0.0, 10.0 - penalty), 2)


def validate_weights(weights: Optional[dict[str, float]] = None) -> dict[str, float]:
    """Validate the preference profile independently of acquisition or scoring."""
    if weights is not None and not isinstance(weights, dict):
        raise ValueError("scoring weights must be a mapping of score terms to numbers")
    unknown = [key for key in (weights or {}) if key not in DEFAULT_WEIGHTS]
    if unknown:
        raise ValueError("unknown scoring weight terms: " + ", ".join(sorted(map(str, unknown))))
    w = DEFAULT_WEIGHTS | (weights or {})
    bad_w = sorted(k for k, v in w.items()
                   if not isinstance(v, (int, float)) or isinstance(v, bool)
                   or _finite(v) is None or _finite(v) < 0.0)
    if bad_w:
        raise ValueError(
            "scoring weights must be finite and non-negative: "
            + ", ".join(f"{k}={w[k]!r}" for k in bad_w)
            + ". A non-finite weight makes every composite in the screen NaN, which "
              "sorts arbitrarily and would be delivered as a ranking.")
    # Match score_screen's exact decimal denominator before acquisition or
    # output: individually finite weights can still overflow their total.
    try:
        denominator = float(sum((Fraction(str(value)) for value in w.values()), Fraction(0)))
    except OverflowError as exc:
        raise ValueError("scoring weights exceed the finite denominator range") from exc
    if not math.isfinite(denominator):
        raise ValueError("scoring weights exceed the finite denominator range")
    return w


def score_screen(candidates: Sequence[Candidate],
                 weights: Optional[dict[str, float]] = None,
                 use_gbsa: bool = True,
                 dock_attempted: Optional[bool] = None) -> list[ScoreBreakdown]:
    """Score a whole screen at once, so structure terms can be rank-normalised.

    When MM/GBSA covers only part of the screen, docking still sets the scale for
    every row and MM/GBSA only permutes the rescored group; see above.

    `dock_attempted` is the caller's statement of whether docking was *run*, as
    opposed to whether it *worked*. Inferring it from the poses -- `any(c.pose)`
    -- is right whenever docking partly succeeded and wrong in the one case that
    matters most: when it ran and failed for every ligand, no candidate has a
    pose, so the whole screen is scored `not_run` and the structure term is
    excused. That is the same hole the comment below was written to close, which
    it does for a corrupt affinity on one row and not for a total failure. Pass
    it explicitly; None keeps the old inference for callers that cannot know.
    """
    # A weight is configuration, not evidence: there is no defensible way to
    # guess what a caller meant by NaN. One non-finite weight makes every
    # composite in the screen NaN, and NaN sorts arbitrarily, so the run would
    # deliver a ranked table built from nothing. Refuse it instead.
    w = validate_weights(weights)

    # Vina's number is parsed from its stdout, so a malformed run can hand us a
    # non-finite affinity. It is no pose at all; see `_percentile_scores` for
    # what one used to do to everybody else's structure score.
    dock_vals = [_finite(c.pose.affinity) if c.pose else None for c in candidates]
    # A converged MM/GBSA result with no dG_bind, or no pose for it to have
    # refined, breaks that stage's contract rather than being a data condition
    # like a silent Vina prep failure: falling through to "ranked on docking
    # alone" would caveat a converged run as unconverged. Fail loudly instead.
    gbsa_vals: list[Optional[float]] = []
    for i, c in enumerate(candidates):
        g = c.gbsa if (use_gbsa and c.gbsa and c.gbsa.converged) else None
        # A non-finite dG_bind breaks the same contract a missing one does, and
        # is the more dangerous of the two: `None` already raised here, while
        # `nan` was handed the best structure slot in the screen.
        dg = _finite(g.dG_bind) if g is not None else None
        if g is not None and (dg is None or dock_vals[i] is None):
            pose_state = ("present" if dock_vals[i] is not None
                          else ("unusable" if c.pose else "missing"))
            raise ValueError(
                f"{c.ligand.ligand_id}: MM/GBSA reports converged but cannot be scored "
                f"(dG_bind={g.dG_bind!r}, pose={pose_state}); it "
                "refines a pose, and without both there is no scale to place it on.")
        gbsa_vals.append(dg)

    dock_s = _percentile_scores(dock_vals, lower_is_better=True)
    # One scale: docking ranks the screen, and MM/GBSA redistributes the
    # rescored subset across its docking slots. Unresolved tiers share a rounded
    # mean; rounding can slightly change the group sum. See `_rescored_structure`.
    structure_s = _rescored_structure(dock_s, gbsa_vals, GBSA_RESOLUTION_KCAL)
    n, n_rescored = len(candidates), sum(1 for v in gbsa_vals if v is not None)
    # Docking ran if anything came back with a pose, even one whose affinity is
    # unreadable. Keying this off a usable number instead would excuse the
    # structure term for the whole screen precisely when the docking output was
    # corrupt, which is the failure the "Missing terms" rule exists to prevent.
    dock_ran = (any(c.pose is not None for c in candidates)
                if dock_attempted is None else bool(dock_attempted))

    structure_enabled = w.get("structure", 0) > 0
    out: list[ScoreBreakdown] = []
    for i, c in enumerate(candidates):
        lig = c.ligand
        caveats: list[str] = []

        raw_pchembl = lig.activity.get("pchembl")
        pchembl = _finite(raw_pchembl)
        pot = potency_score(pchembl)
        measured = pchembl is not None
        if raw_pchembl is not None and pchembl is None:
            caveats.append(_nonfinite_note("pChEMBL", raw_pchembl))
        recorded_status = lig.activity.get("activity_status")
        if not measured and (recorded_status in {"MEASURED", "MEASURED_QUALITATIVE"}
                             or lig.activity.get("primary_measurement_status") == "AMBIGUOUS_CONTEXTS"):
            caveats.append("recorded activity measurements exist; no eligible exact primary potency "
                           "was selected. Inspect measurement_groups and primary_measurement_limitation; "
                           "potency 0 is a missing ranking term, not measured inactivity")
        elif not measured and dock_vals[i] is None:
            caveats.append("no affinity evidence against this node, measured or docked; "
                           "it is a reference compound here, not a hit")
        elif not measured:
            caveats.append("no measured affinity; ranked on docking alone" if structure_enabled
                           else "no exact primary potency; docking is diagnostic and has zero ranking weight")

        if gbsa_vals[i] is not None:
            structure, has_structure, ev = structure_s[i], True, "mmgbsa"
            if n_rescored < n:
                caveats.append(f"one of the {n_rescored} of {n} candidates MM/GBSA "
                               "rescored; structure is MM/GBSA order within that group, "
                               "on the screen's docking scale")
        elif dock_vals[i] is not None:
            structure, has_structure, ev = structure_s[i], True, "dock"
            if use_gbsa and c.gbsa is not None:
                caveats.append("MM/GBSA ran on this compound and did not converge; "
                               + ("ranked on docking alone" if structure_enabled
                                  else "the docking value is diagnostic only"))
            elif n_rescored:
                caveats.append(f"not among the {n_rescored} of {n} candidates MM/GBSA "
                               + ("rescored; ranked on docking alone" if structure_enabled
                                  else "rescored; docking is diagnostic only"))
        else:
            # No pose, and three different reasons for it.
            #
            # `dock_ran` is a screen-level fact and was used as a per-compound
            # one, so every candidate outside a `--dock-top` budget was told
            # "docking returned no pose for this compound ... it is not
            # excused" -- a failure report about a run that never happened.
            # `Candidate.dock_attempted` answers for the compound; None means
            # nothing set it and the screen-level flag is all there is.
            attempted = getattr(c, "dock_attempted", None)
            if attempted is None:
                attempted = dock_ran
            structure, has_structure = 0.0, bool(attempted)
            if c.pose is not None:
                ev = "none"
                caveats.append(
                    f"docking returned an unusable affinity ({c.pose.affinity!r}) for "
                    "this compound; it counts as no pose and the structure term "
                    "scores 0, it is not excused")
            elif attempted:
                ev = "none"
                caveats.append("docking returned no pose for this compound; "
                               "the structure term scores 0, it is not excused")
            elif dock_ran:
                # Docking ran for the screen; this compound was not in the
                # slice. The structure term is 0 here the same way it is for a
                # screen that never docked -- `has_structure` is False, so the
                # weights renormalise without it and the 0 does not count
                # against the compound. Changing that is a scoring decision,
                # not a labelling one, and is deliberately not taken here.
                ev = "budget_skipped"
                caveats.append(
                    "outside the --dock-top budget, so this compound was never "
                    "offered to the docking engine; it has no structure "
                    "evidence and this is not a docking failure. Raise "
                    "--dock-top, or pass 0, to place it on the same footing as "
                    "the rows above it")
            else:
                ev = "not_run"

        raw_topical = c.skin.topical_score if c.skin else 0.0
        topical = _finite(raw_topical)
        if topical is None:
            # Left unguarded this is the worst of the lot: it does not merely
            # mis-score one term, it makes the composite itself NaN, and a NaN
            # composite takes whatever rank `sorted` happens to give it.
            topical = 0.0
            caveats.append(_nonfinite_note("topical score", raw_topical))
        elif not 0 <= topical <= 10:
            caveats.append(f"recorded topical score {raw_topical!r} is outside the 0-10 domain; "
                           "it is treated as missing and the topical term scores 0")
            topical = 0.0





        burden = sum(SEVERITY_WEIGHT.get(a, 0.0)
                     for a in [x.split(":")[0] for x in _electrophilic_alerts(
                         c.skin.sensitiser_alerts if c.skin else [])])
        n_pains = c.skin.n_pains if c.skin else 0
        safety_evaluated = skin_filter_status(c, "sensitisation") == "EVALUATED"
        stability_evaluated = skin_filter_status(c, "stability") == "EVALUATED"
        safety = safety_score(burden, n_pains=n_pains) if safety_evaluated else 5.0
        if safety_evaluated:
            from .filters.safety import coverage_limitations
            caveats.extend(coverage_limitations(lig.smiles))
        if c.skin and not safety_evaluated:
            caveats.append("sensitisation and PAINS risk not assessed; the safety "
                           "score is a neutral placeholder, not a clean result")
        # The colour call is a stability judgement that is deliberately not a
        # stability *alert*, so it never reached the alert count and the score
        # was blind to it: anthracene came back "moderate (may tint)" in the
        # verdict sentence and 10.0 -- a clean sweep -- in the column a reader
        # ranks on.
        colour = colour_of(c)
        stab = (stability_score(len(c.skin.stability_alerts), colour)
                if stability_evaluated else 5.0)
        if c.skin is None:
            # The caveat below already says this whole row is a placeholder;
            # adding "colour unknown" on top would read as a colour finding.
            caveats.append("no skin assessment; safety and stability are placeholders, "
                           "not measurements")
        elif not stability_evaluated or colour == "unknown":
            caveats.append("stability and colour could not be assessed for this structure; "
                           "the stability score is a neutral placeholder, not a clean result")
        elif colour != "low":
            caveats.append(f"predicted colour risk {colour}; it may tint a white base; "
                           "verify colour in the intended formulation")
        if safety_evaluated and _nonnegative_count(n_pains) is None:


            if _finite(n_pains) is None:
                caveats.append(_nonfinite_note("PAINS count", n_pains))
            else:
                caveats.append(f"recorded PAINS count {n_pains!r} is not a finite non-negative "
                               "integer; it is treated as missing and the safety term scores 0")
        elif safety_evaluated and n_pains:
            caveats.append(f"{n_pains} PAINS substructure(s): a frequent-hitter "
                           "chemotype, so treat any measured potency as unconfirmed")

        ef = EVIDENCE_FACTOR.get(c.evidence, 0.6)
        terms = {"topical": topical, "safety": safety,
                 "stability": stab, "potency": pot}
        if has_structure:
            terms["structure"] = structure
        # Exact rational arithmetic on decimal representations defines the result.
        # Binary summation order changed a real F2RL1 score from 4.36 to 4.37
        # at an exact halfway point; a tolerance would conceal that defect.
        denominator = sum((Fraction(str(w[k])) for k in terms), Fraction(0)) or Fraction(1)
        numerator = sum((Fraction(str(w[k])) * Fraction(str(v))
                         for k, v in terms.items()), Fraction(0))
        def rounded_mean(factor):
            cents = numerator * factor / denominator * 100
            whole, remainder = divmod(cents.numerator, cents.denominator)
            twice_remainder = 2 * remainder
            if twice_remainder > cents.denominator or (twice_remainder == cents.denominator and whole % 2):
                whole += 1
            return whole / 100.0

        composite = rounded_mean(Fraction(str(ef)))
        unmultiplied = rounded_mean(Fraction(1))
        denom = float(denominator)
        if not math.isfinite(denom):
            raise ValueError("included scoring weights exceed the finite denominator range")

        if c.skin and c.skin.permeation_class == "blocked":
            caveats.append("passive permeation model predicts limited penetration through intact skin; "
                           "delivery performance requires formulation testing")
        # `sensitiser_alerts` being non-empty was the whole test, and
        # `assess_skin` appends `low:PAINS(n)` into that same list. So a
        # compound that matched no electrophilic alert at all, and whose only
        # entry was a frequent-hitter flag, went out with a caveat commissioning
        # three sensitisation assays -- DPRA, KeratinoSens and h-CLAT -- for an
        # assay-interference finding that is already caveated separately above.
        # The two endpoints answer different questions and cost different money.
        if c.skin and _electrophilic_alerts(c.skin.sensitiser_alerts):
            caveats.append("electrophilic alert; run DPRA/KeratinoSens/h-CLAT before use")

        if not structure_enabled:
            caveats.append("structure weight 0; raw docking/MMGBSA values do not affect the composite or rank group")

        bits = [f"topical {topical:.1f}", f"potency {pot:.1f}", f"safety {safety:.1f}"]
        if has_structure:
            bits.append(f"structure {structure:.1f} ({STRUCTURE_EVIDENCE[ev][0]}"
                        + (", diagnostic, weight 0)" if not structure_enabled else ")"))
        bits.append(f"stability {stab:.1f}")
        rationale = ", ".join(bits) + f" (evidence x{ef:.2f})"

        out.append(ScoreBreakdown(ligand_id=lig.ligand_id, potency=pot, structure=structure,
                                  topical=topical, safety=safety, stability=stab,
                                  evidence_factor=ef, composite=composite,
                                  composite_without_evidence_factor=unmultiplied,
                                  evidence_factor_contribution=round(composite - unmultiplied, 2),
                                  rationale=rationale, caveats=caveats,
                                  structure_evidence=ev, weights_used=dict(w),
                                  weighted_terms=sorted(terms), weight_denominator=denom))

    # Only a positive structure weight makes the term sets incomparable.
    # At zero weight, applying an evaluated-first partition would still give
    # docking an unsupported rank bonus despite its zero numeric contribution.
    evaluated = [i for i in range(len(out)) if out[i].structure_evidence in EVALUATED]
    backlog = [i for i in range(len(out)) if out[i].structure_evidence not in EVALUATED]
    if w.get("structure", 0) == 0 or not evaluated or not backlog:
        groups = [list(range(len(out)))]
    else:
        groups = [evaluated, backlog]
    rank = 0
    for group in groups:
        for local_rank, i in enumerate(sorted(group, key=lambda j: out[j].composite, reverse=True), 1):
            rank += 1
            candidates[i].composite = out[i].composite
            candidates[i].rank = rank
            out[i].rank_group = ("shared_terms" if w.get("structure", 0) == 0 else "evaluated" if len(groups) == 1
                                 else ("evaluated" if i in set(evaluated)
                                       else "not_evaluated"))
            out[i].rank_within_group = local_rank
            out[i].rank_group_size = len(group)
            if len(groups) > 1:
                out[i].rank_interpretation = ("grouped display ordinal: evaluated group precedes never-docked "
                    "group; composites use different denominators and global rank is not a global score sort")
    return out


def structure_evidence_summary(breakdowns: Sequence[ScoreBreakdown]) -> dict:
    """Counts and a header line for the structure column, for `report.py`.

    With MM/GBSA on 6 of 12 this returns ``mixed=True`` and a ``note`` to print
    above the ranked table, so a reader is told outright that half the rows carry
    a rescored term. Per row: ``ScoreBreakdown.structure_evidence``.
    """
    counts = {k: sum(b.structure_evidence == k for b in breakdowns)
              for k in STRUCTURE_EVIDENCE}
    n, k = len(breakdowns), counts["mmgbsa"]
    note = f"MM/GBSA rescored {k} of {n} candidates"
    if 0 < k < n:
        note += ("; docking sets the structure scale for the whole screen and MM/GBSA "
                 "reorders that group within the scores it already holds, so the column "
                 "is one scale throughout and structure_evidence names the method row "
                 "by row")












    n_unevaluated = counts["budget_skipped"] + counts["not_run"]
    n_structure = counts["mmgbsa"] + counts["dock"] + counts["none"]
    active = any(b.weights_used.get("structure", DEFAULT_WEIGHTS["structure"]) > 0 for b in breakdowns)
    incomparable = bool(active and n_unevaluated and n_structure)
    if not active:
        note += ". Structure weight is zero; these values are diagnostics and do not contribute to the composite or partition the ranking."
    if incomparable:
        note += (f". {n_unevaluated} of {n} candidates were never docked, so "
                 "their composite is a mean over fewer terms than the docked "
                 "rows and the two are not on one scale. Compare within a "
                 "`structure_evidence` group, not down the column")
    return {"counts": counts, "n": n, "n_rescored": k, "mixed": 0 < k < n,
            "n_unevaluated": n_unevaluated, "incomparable": incomparable,
            "structure_weight_active": active,
            "note": note}




















MIN_CONTROLS_FOR_POWER = 3

#: Below this, a rank correlation over so few pairs is noise with a decimal
#: point, and printing one invites the reader to act on it.
MIN_PAIRS_FOR_RHO = 10
RHO_BOOTSTRAP_DRAWS = 200

#: |rho| under this is reported as "not a potency ordering". Not a
#: significance test -- it is the threshold at which the composite has stopped
#: tracking the one measured number in the screen, which is the thing a reader
#: has to know before treating rank 1 as the most potent compound.
RHO_UNCORRELATED = 0.2


def _spearman(pairs: Sequence[tuple[float, float]]) -> Optional[float]:
    """Spearman rho over (x, y) pairs, ties averaged. None when undefined.

    Pure Python and O(n log n) on a few hundred rows: cheap enough to run
    inside `control_check`, which is the only reason it is computed there
    rather than quoted from a review document. A constant column has no rank
    variance and returns None rather than a ZeroDivisionError or a 0.0 that
    reads as a measurement.
    """
    if len(pairs) < 2:
        return None

    def ranks(values: list[float]) -> list[float]:
        order = sorted(range(len(values)), key=lambda i: values[i])
        out = [0.0] * len(values)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
                j += 1
            average = (i + j) / 2 + 1
            for k in range(i, j + 1):
                out[order[k]] = average
            i = j + 1
        return out

    xs, ys = ranks([x for x, _ in pairs]), ranks([y for _, y in pairs])
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    num = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    den = math.sqrt(sum((a - mx) ** 2 for a in xs)
                    * sum((b - my) ** 2 for b in ys))
    return num / den if den else None


def _potency_composite_pairs(candidates: Sequence[Candidate]
                             ) -> list[tuple[float, float]]:
    """Ranked candidates carrying both a finite pChEMBL and a finite composite."""
    out = []
    for c in candidates:
        if not c.rank:
            continue
        p = _finite(c.ligand.activity.get("pchembl"))
        s = _finite(c.composite)
        if p is not None and s is not None:
            out.append((p, s))
    return out


def _unusable_count(candidates: Sequence[Candidate]) -> int:
    """How many candidates carry a pChEMBL that is recorded but not a number."""
    return sum(1 for c in candidates
               if c.ligand.activity.get("pchembl") is not None
               and _finite(c.ligand.activity.get("pchembl")) is None)


def _descriptive_rho_interval(pairs):
    """Seeded row bootstrap, conditional on exchangeable rows, not validation.

    Study/scaffold clusters are not available here. This interval therefore
    cannot establish independence, efficacy or prospective ranking accuracy.
    """
    if len(pairs) < MIN_PAIRS_FOR_RHO:
        return None
    ordered = sorted(pairs)
    rng = random.Random(20261017)
    draws = [_spearman(rng.choices(ordered, k=len(ordered))) for _ in range(RHO_BOOTSTRAP_DRAWS)]
    draws = sorted(value for value in draws if value is not None)
    if len(draws) < RHO_BOOTSTRAP_DRAWS * .9:
        return None

    def quantile(fraction):
        position = fraction * (len(draws) - 1)
        low = int(position)
        high = min(low + 1, len(draws) - 1)
        return round(draws[low] + (draws[high] - draws[low]) * (position - low), 3)

    return [quantile(.025), quantile(.975)]


def control_check(candidates: Sequence[Candidate]) -> dict:
    """Did the screen rank the known actives in the better half?

    A screen that buries its positive controls is not reporting hits, it is
    reporting noise. This is the first thing to look at in any result.

    Every branch below has to be defensible on its own, because this verdict is
    the line a reader acts on. In particular it must never claim the controls
    lost to "comparably potent" library members when no library member carries a
    measured potency at all, and it must not report a missing positive control
    when the real problem is that nothing was ranked.

    That rule is symmetric and the branch order enforces it: with no comparator
    carrying a measured potency the check says so *before* it looks at where the
    controls landed, because a good rank against an unmeasured library is worth
    no more than a bad one. It is also the case that matters most -- a screen
    containing nothing but its own positive controls ranks them first by
    construction, and the old order certified that as "the ranking is behaving".

    The same rule applies to the *number* of controls, which is the other way
    this check certified a screen on nothing. The median split is a coin flip
    per control, so a pass below `MIN_CONTROLS_FOR_POWER` says only that the
    screen did not bury its controls -- never that the ranking is sound -- and
    the verdict names K so the reader can see how thin the evidence is. The
    check is not skipped; it states what it can and cannot conclude.

    Every verdict also carries the rank correlation between measured potency
    and the composite, because "the controls are near the top" and "the
    ranking orders compounds by potency" are different claims and only the
    first one has ever been reported here.
    """
    node = next((c.target for c in candidates), None)

    def on_node(c) -> bool:
        return (c.ligand.source == "cosmetic_active"
                and (node is None
                     or c.ligand.activity.get("node", "").upper() == str(node).upper()))

    references = [c for c in candidates if on_node(c)]
    controls = [c for c in references if c.ligand.activity.get('curated_claim_status')
                != 'UNREVIEWED_TARGET_ENDPOINT']
    n_unreviewed = len(references) - len(controls)
    if not any(c.rank for c in candidates):
        return {"n_controls": 0, "n_screened": len(candidates),
                "n_unreviewed_references": n_unreviewed,
                "n_unusable_potency": _unusable_count(candidates),
                "verdict": ("the screen was never scored, so there are no ranks to "
                            "check. This says nothing about the controls"
                            + (f"; {len(controls)} curated active(s) for {node} are "
                               "present and unranked." if controls else "."))}
    ranks = [(c.ligand.name, c.rank) for c in controls if c.rank]
    if not ranks:
        return {"n_controls": 0, "n_screened": len(candidates),
                "status": "NOT_EVALUATED", "n_unreviewed_references": n_unreviewed,
                "n_unusable_potency": _unusable_count(candidates),
                "verdict": (f"{n_unreviewed} annotated reference(s) assigned to {node} are unreviewed "
                            "for the target endpoint; no validated positive control available"
                            if n_unreviewed else
                            f"no curated active targets {node}; no positive control available")}

    n = len(candidates)
    # "Better half", not "strictly above n/2": in a one-compound screen the only
    # rank is 1 and n/2 is 0.5, which would call the top-ranked control a failure.
    median_rank = math.ceil(n / 2)
    above = [r for _, r in ranks if r <= median_rank]
    frac = len(above) / len(ranks)

    # `is not None` was the whole test here, so a NaN pChEMBL counted as a
    # measurement: `max(measured) < min(others)` is False against NaN, and the
    # verdict fell through to "controls rank below comparably potent library
    # members" -- a screen failure declared on the strength of a corrupt record.
    def _potencies(rows) -> tuple[list[float], int]:
        vals = [(c, _finite(c.ligand.activity.get("pchembl"))) for c in rows]
        unusable = sum(1 for c, f in vals
                       if f is None and c.ligand.activity.get("pchembl") is not None)
        return [f for _, f in vals if f is not None], unusable

    measured, bad_controls = _potencies([c for c in controls if c.rank])
    # CHEMBL_LIBRARY_SOURCES, not a literal: a compound rescued from the local
    # chemreps archive is labelled `chembl_local` and was silently dropped from
    # this pool, which emptied it during exactly the API outages that rescue
    # exists for. See the constant for why it is a constant.
    others, bad_others = _potencies([c for c in candidates
                                     if c.ligand.source in CHEMBL_LIBRARY_SOURCES])
    n_unusable = bad_controls + bad_others

    if not measured:
        verdict = ("controls carry no measured potency for this node, so their rank "
                   "says nothing about the screen. Add literature values to the "
                   "actives file to make this check meaningful.")
    elif not others:
        # Above the `frac` test, not below it. Ranking well against a library
        # that measured nothing is no more informative than ranking badly
        # against it, and when the library *is* the controls -- a one-compound
        # screen, or an actives file run on its own -- the controls ranked
        # against themselves and the old order printed "the ranking is
        # behaving" over `1 / 1`. The docstring commits to the mirror rule for
        # the failure direction; it holds in both directions or it holds in
        # neither.
        verdict = ("no library compound carries a measured potency, so whichever "
                   "way the controls ranked there is nothing to rank them "
                   "against. This check cannot pass or fail on this screen.")
    elif len(ranks) < MIN_CONTROLS_FOR_POWER:




        verdict = (
            f"underpowered: this screen has {len(ranks)} control(s), and a "
            f"median split over K={len(ranks)} cannot tell a working ranking "
            "from a random one -- a random ranking passes this test about "
            "half the time at K=1 and one time in four at K=2. What it does "
            f"show: {len(above)} of {len(ranks)} control(s) ranked at or above "
            f"{median_rank} of {n}. What it "
            "cannot show is that the ranking is sound. At least "
            f"{MIN_CONTROLS_FOR_POWER} controls on this node are needed before "
            "a pass here is evidence.")
    elif frac >= 0.6:
        verdict = "controls rank as expected; the ranking is behaving"
    elif max(measured) < min(others):
        verdict = (f"controls have lower recorded potency "
                   f"(best control pChEMBL {max(measured):.1f} vs library minimum "
                   f"{min(others):.1f}); library inclusion thresholds and assay "
                   "context may explain this comparison. It does not establish "
                   "that the ranking works or excuse a failed control check.")
    else:
        verdict = "controls rank below comparably potent library members; investigate before acting"

    if n_unusable:
        # The count is part of the verdict, not a side note: a reader has to be
        # able to tell "no library compound carries a measured potency" from
        # "the ones that did were unreadable and were thrown away".
        verdict += (f" ({n_unusable} pChEMBL record(s) were not finite numbers and "
                    "were discarded before this check; they are neither controls "
                    "nor comparators here.)")
    if len(ranks) < MIN_CONTROLS_FOR_POWER and "underpowered" not in verdict:
        verdict += (f" Underpowered: only {len(ranks)} control(s); neither a high "
                    "nor a low rank validates this screen.")







    pairs = _potency_composite_pairs(candidates)
    rho = _spearman(pairs) if len(pairs) >= MIN_PAIRS_FOR_RHO else None
    rho = None if rho is None else round(rho, 3)
    rho_interval = _descriptive_rho_interval(pairs) if rho is not None else None
    if rho is not None:
        verdict += (f" Separately, Spearman rho between measured pChEMBL and the "
                    f"composite is {rho:+.3f} over the {len(pairs)} ranked "
                    "compounds carrying both"
                    + ("; the ranking is not a potency ordering."
                       if abs(rho) < RHO_UNCORRELATED else "."))
    else:
        verdict += (f" Potency/composite correlation not evaluated: at least "
                    f"{MIN_PAIRS_FOR_RHO} usable pairs and nonconstant columns are required.")
    verdict += (" This is a descriptive same-data association: composite already "
                "includes potency, so its correlation with pChEMBL is partly "
                "circular and is not independent ranking validation.")

    return {"n_controls": len(ranks), "n_above_median": len(above),
            "fraction_above_median": round(frac, 2),
            "n_with_measured_potency": len(measured),
            "n_comparators_with_potency": len(others),
            "n_unusable_potency": n_unusable, "n_screened": n,
            "median_rank": median_rank,
            "underpowered": len(ranks) < MIN_CONTROLS_FOR_POWER,
            "min_controls_for_power": MIN_CONTROLS_FOR_POWER,
            "spearman_potency_composite": rho,
            "spearman_potency_composite_ci_95": rho_interval,
            "spearman_status": "DESCRIPTIVE_ONLY" if rho is not None else "NOT_EVALUATED",
            "spearman_policy": "same-data-row-bootstrap-v1",
            "spearman_bootstrap_draws": RHO_BOOTSTRAP_DRAWS if rho is not None else 0,
            "spearman_interval_assumption": "exchangeable rows; study/scaffold dependence unmodelled",
            "spearman_independent_validation": False,
            "min_pairs_for_rho": MIN_PAIRS_FOR_RHO,
            "n_potency_composite_pairs": len(pairs),
            "control_ranks": sorted(ranks, key=lambda t: t[1]),
            "best": min(ranks, key=lambda t: t[1]), "verdict": verdict}
