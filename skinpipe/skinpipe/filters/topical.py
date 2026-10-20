'Stratum-corneum permeation and topical feasibility.'
from __future__ import annotations

import logging
import math
from typing import Any, Optional

from ..types import LigandRecord, SkinAssessment

log = logging.getLogger(__name__)

# Potts-Guy coefficients
_PG_INTERCEPT, _PG_LOGP, _PG_MW = -2.72, 0.71, -0.0061

# Permeation classes, keyed on log Kp (cm/h)
_KP_BANDS = (
    (-2.0, "good"),      # >= -2.0
    (-3.5, "moderate"),  # -3.5 .. -2.0
    (-5.0, "poor"),      # -5.0 .. -3.5
)

# The 500 Da rule, with a soft shoulder rather than a cliff
_MW_IDEAL, _MW_SOFT, _MW_HARD = 350.0, 500.0, 700.0
_LOGP_LO, _LOGP_HI = 1.0, 3.5          # window where both partition steps work

# Unit conversion for `estimate_flux`, which reports ug/h. Kp*C*A comes out in
# mg/h (cm/h * mg/cm3 * cm2); this is the only factor that should be applied.
_UG_PER_MG = 1000.0


def potts_guy_log_kp(logp: float, mw: float) -> float:
    """log10 permeability coefficient in cm/h through human stratum corneum."""
    return _PG_INTERCEPT + _PG_LOGP * float(logp) + _PG_MW * float(mw)


def permeation_class(log_kp: float, mw: float) -> str:
    """Coarse band, with an explicit override for molecules too large to pass."""
    if mw > _MW_HARD:
        return "blocked"
    for threshold, label in _KP_BANDS:
        if log_kp >= threshold:
            return label
    return "poor"


def estimate_flux(log_kp: float, solubility_mg_ml: float, area_cm2: float = 100.0) -> float:
    'Validate the selected structure and retain explicit calculation provenance.'
    kp = 10.0 ** log_kp                                  # cm/h
    flux_mg_per_h = kp * solubility_mg_ml * area_cm2     # cm/h * mg/cm3 * cm2
    return flux_mg_per_h * _UG_PER_MG                    # -> ug/h


def _score_window(value: float, lo: float, hi: float, soft: float) -> float:
    """10 inside [lo, hi], decaying to 0 at `soft` distance outside."""
    if lo <= value <= hi:
        return 10.0
    gap = lo - value if value < lo else value - hi
    return max(0.0, 10.0 * (1.0 - gap / soft))


def _score_mw(mw: float) -> float:
    if mw <= _MW_IDEAL:
        return 10.0
    if mw <= _MW_SOFT:
        return 10.0 - 6.0 * (mw - _MW_IDEAL) / (_MW_SOFT - _MW_IDEAL)   # 10 -> 4
    if mw <= _MW_HARD:
        return max(0.0, 4.0 - 4.0 * (mw - _MW_SOFT) / (_MW_HARD - _MW_SOFT))
    return 0.0


#: The value each polarity term takes when the descriptor is absent: the
#: worst end of its own window in `topical_score`, so a missing number can
#: never help a compound through a gate.
_WORST_HBD = 99
_WORST_TPSA = 1000.0


def topical_score(mw: float, logp: float, hbd: int, tpsa: float,
                  weights: Optional[dict[str, float]] = None) -> float:
    """0-10 composite of the four descriptors that govern SC permeation.

    Weights are exposed because the right balance is formulation-dependent: a
    leave-on serum tolerates a slower permeant than a rinse-off product.
    """
    w = {"mw": 0.40, "logp": 0.30, "hbd": 0.15, "tpsa": 0.15} | (weights or {})
    s_mw = _score_mw(float(mw))
    s_logp = _score_window(float(logp), _LOGP_LO, _LOGP_HI, soft=3.0)
    s_hbd = 10.0 if hbd <= 3 else (6.0 if hbd <= 5 else 2.0)
    s_tpsa = 10.0 if tpsa <= 90 else (6.0 if tpsa <= 130 else 2.0)
    total = w["mw"] * s_mw + w["logp"] * s_logp + w["hbd"] * s_hbd + w["tpsa"] * s_tpsa
    return round(total / sum(w.values()) * 1.0, 2) if sum(w.values()) else 0.0


def _finite_or_none(value: Any) -> Optional[float]:
    """A real number, or None. NaN, +/-inf and non-numerics are all absence."""
    if value is None or isinstance(value, bool):
        return None
    try:
        out = float(value)
    except (TypeError, ValueError, OverflowError):
        return None
    return out if math.isfinite(out) else None


def assess(ligand: LigandRecord, props: Optional[dict[str, Any]] = None) -> SkinAssessment:
    """Build the permeation half of a SkinAssessment from cached descriptors.

    Safety and stability alerts are filled in by the sibling modules; this
    function never overwrites them.
    """
    p = dict(ligand.props) | (props or {})
    # `bool` is a subclass of `int`, so `float(True)` is 1.0 and sails through
    # every numeric check below -- measured: mw=True scored 7.6 with
    # permeation "good", evaluated as a 1 Da molecule, and logp=False scored
    # 6.6. A flag that reached a descriptor slot is a plumbing error, not a
    # measurement, and it must not be ranked as one.
    for _name in ("mw", "logp"):
        if isinstance(p.get(_name), bool):
            log.warning("ligand %s has a boolean in props[%r]; permeation not "
                        "assessed", ligand.ligand_id, _name)
            return SkinAssessment(
                ligand_id=ligand.ligand_id,
                verdict=f"descriptors not usable ({_name} is a boolean, not a number)")
    try:
        mw = float(p["mw"])
        logp = float(p["logp"])
    except (KeyError, TypeError, ValueError, OverflowError):
        log.warning("ligand %s lacks mw/logp; cannot assess permeation", ligand.ligand_id)
        return SkinAssessment(ligand_id=ligand.ligand_id, verdict="descriptors missing")
    # `float()` accepts "nan" and "inf", and Potts & Guy is linear in both terms,
    # so a non-finite descriptor propagates into `log_kp` and then scores. The
    # measured outcome, before this check: logp=inf gave log_kp=inf, a finite
    # topical_score of 4.6 and `permeation_class="good"` -- a garbage descriptor
    # producing a *favourable* verdict. mw=nan gave log_kp=nan with a score of
    # 3.6. A number that is not a number is not a measurement, so it is treated
    # as the absence it is rather than carried into the model.
    # A mass or a surface area cannot be negative, and Potts & Guy is linear in
    # log MW, so a negative one does not fail -- it scores *well*. Measured:
    # mw=-100 came back `permeation_class="good"`, `topical_score=10.0`,
    # verdict "permeates unaided", outranking a real 350 Da molecule. That is
    # the same class of defect as the non-finite check below: a descriptor that
    # is not a measurement must not produce a favourable verdict.
    if mw <= 0.0:
        log.warning("ligand %s has a non-physical mass (mw=%r); permeation not "
                    "assessed", ligand.ligand_id, mw)
        return SkinAssessment(
            ligand_id=ligand.ligand_id,
            verdict=f"descriptors not usable (mw={mw!r}: mass must be positive)")
    if not (math.isfinite(mw) and math.isfinite(logp)):
        log.warning("ligand %s has non-finite descriptors (mw=%r, logp=%r); "
                    "permeation not assessed", ligand.ligand_id, mw, logp)
        return SkinAssessment(
            ligand_id=ligand.ligand_id,
            verdict=f"descriptors not usable (mw={mw!r}, logp={logp!r})")

    # Absence must not score better than a measurement. These defaulted to 0,
    # which is the *best* value on both windows -- `topical_score` awards 10.0
    # for hbd <= 3 and 10.0 for tpsa <= 90 -- so a PubChem record missing both
    # gained over a point on a gate set at 6.0. `catalog.fetch_properties`
    # returns None for any of the four, and `screen` guards only mw and logp,
    # so this is the normal shape of a catalog row, not an edge case.
    #
    # Missing terms now take their worst value and are named in the verdict:
    # the compound stays in the screen, which is right, but it cannot pass a
    # threshold on a number nobody measured.
    missing: list[str] = []
    # Same rule as mw/logp above, and for the same reason: `int(True)` is 1 and
    # `float(False)` is 0.0, both of which sit at the *best* end of these two
    # windows -- measured, `hbd=True` and `tpsa=False` each scored 10.0. A flag
    # in a descriptor slot is treated as the absence it is, which takes the worst
    # value and is named in the verdict.
    raw_hbd = None if isinstance(p.get("hbd"), bool) else p.get("hbd")
    raw_tpsa = None if isinstance(p.get("tpsa"), bool) else p.get("tpsa")
    # Finite, not merely convertible. `int(float("inf"))` raises
    # `OverflowError`, which was not in the caught set, so one contaminated
    # descriptor aborted the whole screen; and `float("nan")` and
    # `float("inf")` convert without complaint, so `tpsa=nan` scored 8.8 and
    # "permeates with a suitable vehicle" with nothing saying the number was
    # missing. Both are the absence they represent: worst value, named in the
    # verdict, same rule as `bool` above and as mw/logp.
    hbd = _finite_or_none(raw_hbd)
    # HBD counts atoms; truncating 3.9 to 3 promoted invalid data into the
    # best scoring band. Negative counts/areas are missing descriptors too.
    hbd = int(hbd) if hbd is not None and hbd >= 0 and hbd.is_integer() else None
    tpsa = _finite_or_none(raw_tpsa)
    if tpsa is not None and tpsa < 0:
        tpsa = None
    if hbd is None:
        missing.append("HBD")
        hbd = _WORST_HBD
    if tpsa is None:
        missing.append("TPSA")
        tpsa = _WORST_TPSA

    kp = potts_guy_log_kp(logp, mw)
    cls = permeation_class(kp, mw)
    score = topical_score(mw, logp, hbd, tpsa)

    if cls == "blocked":
        verdict = (f"MW {mw:.0f} exceeds the model's {_MW_HARD:g} Da screening threshold; "
                   "passive penetration is predicted to be limited")
    elif cls == "good":
        verdict = "Model predicts relatively favorable passive penetration; skin uptake requires measurement"
    elif cls == "moderate":
        verdict = "Model predicts moderate passive penetration; vehicle effects require measurement"
    else:
        verdict = "Model predicts low passive penetration; delivery and surface effects require measurement"

    if missing:
        verdict += (f" ({' and '.join(missing)} not available; scored at the "
                    "worst value, so this score is a lower bound)")

    return SkinAssessment(
        ligand_id=ligand.ligand_id,
        log_kp=round(kp, 2),
        permeation_class=cls,
        topical_score=score,
        verdict=verdict,
    )


def needs_delivery_system(assessment: SkinAssessment) -> bool:
    """True when the molecule is worth pursuing only with encapsulation."""
    return assessment.permeation_class in {"blocked", "poor"}
