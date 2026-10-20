"""Systemic exposure after topical use, and the liabilities that depend on it.

The question this module answers
--------------------------------
For an oral drug the relevant quantity is bioavailability. For a cosmetic
applied to skin it is the Systemic Exposure Dosage: how many milligrams per
kilogram of body weight per day actually reach the circulation. The SCCS Notes
of Guidance (SCCS/1647/22, section 3-4.5) defines it as

    SED [mg/kg bw/day] = A [mg product/kg bw/day]
                         * C [% active] / 100
                         * DAp [% dermal absorption] / 100

and the safety case is then the Margin of Safety,

    MoS = NOAEL / SED,      which a dossier is expected to clear at 100.

A is fixed by the product format, and this module tabulates the SCCS daily-use
amounts. C is what the formulator chooses, which is why `applied_percent` is
exposed as an argument: the point is to be able to ask "what use level can I
defend?" rather than to receive one number.

DAp is the interesting term. SCCS lets a dossier default to 50% when there is
no measurement, and 10% when the substance has MW > 500 Da together with a
logP below -1 or above 4. These are dossier defaults, not measurements for a
particular substance. This module also derives a screening DAp scenario from the
Potts-Guy permeability coefficient already implemented in `filters.topical`,
under a finite-dose mass balance, and reports both numbers side by side.

Assumptions, stated rather than buried
--------------------------------------
    relative product exposure   SCCS Table 3A (retention already included)
    product density             1 g/cm3, so C% -> 10*C mg/cm3 in the vehicle
    vehicle                     aqueous; Potts-Guy was fitted to aqueous
                                vehicles and real emulsions, glycols and
                                ethanol shift Kp by up to an order of magnitude
    contact time                24 h for a leave-on product, 15 min rinse-off
    absorption ceiling          the applied amount itself. Flux is steady-state
                                and a thin film is not, so the model is capped
                                at mass balance within this model scenario
    lag time                    ignored, which makes the estimate conservative
                                for slow permeants

CYP and hERG live here rather than in a pharmacology module because for a
topical they are conditional findings: a CYP3A4 inactivator that never reaches
plasma is an academic observation. They are reported with their relevance
not established by an unmeasured SED scenario alone.
"""
from __future__ import annotations

import logging
import math
from typing import Any, Optional

from rdkit import Chem

from . import admet_alerts as al
from . import admet_cmr
from . import admet_photo as photo
# `_MW_SOFT` (500.0) is the Bos-Meinardi cut-off and `topical` owns it.
# This module carried its own bare `500.0` in two places: the same number
# today, two numbers the first time either is tuned. The Potts-Guy
# coefficients were already single-sourced this way; the MW cut-off was
# the copy that got left behind.
from .topical import _MW_SOFT, potts_guy_log_kp

log = logging.getLogger(__name__)

BODY_WEIGHT_KG = 60.0
MOS_TARGET = 100.0
PRODUCT_DENSITY_G_CM3 = 1.0

# name -> (daily amount g, retention factor, treated area cm2, contact hours)
# Daily amounts from SCCS/1647/22 Table 3A; area/contact times are scenario
# assumptions. A is tabulated separately: the study uses individual body
# weights, so the published relative exposure is not daily grams / 60 kg.
PRODUCT_USE: dict[str, tuple[float, float, float, float]] = {
    "body_lotion": (7.82, 1.00, 15670.0, 24.0),
    "face_cream":  (1.54, 1.00, 565.0, 24.0),
    "hand_cream":  (2.16, 1.00, 860.0, 24.0),
    "lip_product": (0.057, 1.00, 4.0, 24.0),
    "body_wash":   (18.67, 0.01, 17500.0, 0.25),
    "shampoo":     (10.46, 0.01, 580.0, 0.25),
}
SCCS_RELATIVE_EXPOSURE_MG_KG_DAY = {
    "body_lotion": 123.20, "face_cream": 24.14, "hand_cream": 32.70,
    "lip_product": 0.90, "body_wash": 2.79, "shampoo": 1.51,
}
EXPOSURE_POLICY = "sccs-table3a-clogp-scenario-v2"
SCCS_SOURCE_URL = "https://health.ec.europa.eu/document/download/32a999f7-d820-496a-b659-d8c296cc99c1_en?filename=sccs_o_273_final.pdf"

# SED bands. These are interpretive, not regulatory: the lowest one is anchored
# on the Cramer class III Threshold of Toxicological Concern (1.5 ug/kg bw/day).
_SED_BANDS: tuple[tuple[float, str, str], ...] = (
    (0.0015, "negligible",
     "low modeled dose; this does not establish TTC eligibility or safety"),
    (0.05, "low",
     "low modeled dose; measured absorption and substance-specific evidence are needed"),
    (0.5, "material",
     "material modeled dose; investigate absorption before choosing a toxicology assessment"),
)
_SED_TOP = ("high", "high modeled dose; actual systemic exposure and toxicology requirements are not established")

#: The use level C is a percentage of the formulation, so the values that
#: describe a formulation at all are the finite ones from 0 to 100. Both ends
#: are inclusive and both are real answers: 0% is "this active is not in this
#: product", which someone may legitimately price, and 100% is the neat active.
APPLIED_PERCENT_MIN = 0.0
APPLIED_PERCENT_MAX = 100.0


def validate_applied_percent(applied_percent: Any) -> float:
    """Return the use level as a float, or refuse it with a diagnostic.

    SED is linear in C, and nothing downstream re-examines its sign or its
    finiteness, so a use level that is not a concentration propagates into a
    complete, confident-looking assessment. Ethanol at `applied_percent=-1`
    returned SED = -1.303333 mg/kg bw/day, a required NOAEL of -130.3333 and the
    `negligible` band -- because "below the Cramer class III TTC" is written as
    `sed < 0.0015`, and every negative number satisfies it. A negative systemic
    dose was therefore reported as the safest result the module can produce.

    The finiteness test comes first, and it is not interchangeable with the
    range test. A range written as `pct < 0 or pct > 100` cannot fire on NaN,
    because every comparison against NaN is False -- and that is exactly what
    `click.FloatRange(min=0, max=100)` compiles to, so a bare range option on
    the CLI would let `--applied-percent nan` through untouched. NaN used to
    reach `_band`, fall off the end of the loop, and be reported as the top band
    with the verdict "SED at most nan mg/kg bw/day".

    This is the one validator for the whole contract. `dermal_absorption`,
    `systemic_exposure`, `admet.systemic_exposure`, `admet.summarise` and the
    CLI's `--applied-percent` all call it, so the Python API and the command
    line cannot come to different answers about the same number.
    """
    if isinstance(applied_percent, bool):
        raise ValueError("applied_percent must be a concentration, not a bool")
    try:
        pct = float(applied_percent)
    except (TypeError, ValueError):
        raise ValueError(
            f"applied_percent must be a number between {APPLIED_PERCENT_MIN:g} "
            f"and {APPLIED_PERCENT_MAX:g} (percent of the formulation the active "
            f"makes up); got {applied_percent!r}") from None
    if not math.isfinite(pct):
        raise ValueError(
            f"applied_percent must be a finite number between "
            f"{APPLIED_PERCENT_MIN:g} and {APPLIED_PERCENT_MAX:g}; got "
            f"{applied_percent!r}. A non-finite use level yields a non-finite "
            f"SED, which compares False against every band threshold and is "
            f"reported as an assessment rather than as the failure it is.")
    if not APPLIED_PERCENT_MIN <= pct <= APPLIED_PERCENT_MAX:
        raise ValueError(
            f"applied_percent must be between {APPLIED_PERCENT_MIN:g} and "
            f"{APPLIED_PERCENT_MAX:g}; got {applied_percent!r}. It is the "
            f"percentage of the formulation the active makes up: a negative use "
            f"level gives a negative SED and a negative required NOAEL, and a "
            f"formulation cannot be more than 100% one ingredient.")
    return pct


def dermal_absorption(logp: float, mw: float, applied_percent: float,
                      product: str) -> dict[str, Any]:
    """Finite-dose dermal absorption, as a percent of the applied active.

    A cosmetic film is finite: roughly 0.5 mg/cm2 of product for a body lotion,
    about 5 micrometres thick. Steady-state flux applied to an infinite dose
    would "absorb" many times that, so what is used here is the depleting-donor
    mass balance, in which the donor concentration falls as the film is drained:

        fraction absorbed = 1 - exp(-Kp * t / h_film)

    with h_film the applied film thickness. Two consequences worth stating.
    The result does not depend on the use level, because the system is linear -
    dermal absorption is a property of the molecule and the film, not of the
    concentration, which is exactly how SCCS treats DAp. And the result is an
    conditional model scenario: omitted lag time and binding can reduce
    absorption, but vehicle effects and estimated logP can move it either way.
    It is not a demonstrated upper bound on measured human absorption.

    Potts-Guy was fitted to polar and moderately lipophilic solutes from
    aqueous vehicles and is known to overestimate Kp for large polyphenols, so
    compounds such as EGCG will read higher here than in-vitro data support.
    """
    applied_percent = validate_applied_percent(applied_percent)
    amount_g, _retention, area_cm2, hours = PRODUCT_USE[product]
    log_kp = potts_guy_log_kp(logp, mw)
    kp = 10.0 ** log_kp                                        # cm/h
    c_vehicle = applied_percent * 10.0 * PRODUCT_DENSITY_G_CM3  # mg/cm3
    flux = kp * c_vehicle                                       # mg/cm2/h
    applied_per_area = (amount_g * 1000.0 / area_cm2
                        * applied_percent / 100.0)              # mg active/cm2
    # This is the applied donor film, before rinsing. Retention belongs once
    # in A (systemic_exposure), not also in DAp's denominator: doing both
    # makes rinse-off absorption larger simply because less product remains.
    film_cm = amount_g / area_cm2 / PRODUCT_DENSITY_G_CM3
    pct = 100.0 * (1.0 - math.exp(-kp * hours / film_cm)) if film_cm else 0.0

    notes: list[str] = ["MODEL_ONLY: RDKit cLogP substitutes for measured log Kow; "
                        "vehicle and absorption are unmeasured. Caps are scenario heuristics, "
                        "not validated human absorption bounds."]
    if _retention != 1.0:
        notes.append("RINSE_OFF_SCENARIO: DAp is conditional on the full pre-rinse "
                     "donor film; SED multiplies that fraction by retention-adjusted "
                     "Table 3A exposure A. This factorisation does not model "
                     "absorption before rinsing separately and may underestimate "
                     "total uptake. It is not a validated rinse-off mass-balance model.")
    saturated = pct > 90.0
    if saturated:
        notes.append("the film is essentially exhausted over the contact time, so "
                     "this saturates near 100%; the real figure is lower because "
                     "lag time, vehicle evaporation and stratum-corneum binding "
                     "are all ignored")
    # SCCS size and partition cut-offs, applied as ceilings on the flux estimate.
    # The two notes are different kinds of statement and only one of them is
    # conditional. "Essentially excluded by an intact barrier" is a fact about
    # a 1000 Da molecule, true whether or not the cap bound, and dropping it
    # would lose the reason the number is small. "The rule *caps* plausible
    # absorption near 10%" asserts that a cap replaced the estimate, and that
    # was appended even when the Potts-Guy figure was already under it -- which
    # tells the reader a number was overridden when nothing touched it.
    # Which cap *bound* the answer, not merely which ones exist. Both can
    # apply: at MW 1100 and logP 12 the film saturates first and the 1000 Da
    # ceiling then fixes the result at 1%, and the verdict called that "the
    # mass-balance cap" because `saturated` was true. It names the last cap
    # that lowered the number instead.
    # What the model said before any ceiling replaced it. The ceilings are
    # step functions on MW, so the reported percent jumps the whole gap at the
    # cut-off: at logP 2 in a body lotion, 499 Da reports 89.78% and 501 Da
    # reports 10.00%. Neither number is wrong for what it is, but only one of
    # them is a model output, and without this field a reader cannot tell
    # which -- nor how far the estimate was from the ceiling that replaced it.
    # The 500 Dalton rule is a screening heuristic about what is *likely* to
    # cross an intact barrier, not a measured bound on every molecule above it.
    model_pct = round(min(pct, 100.0), 2)
    ceiling_applied = False
    binding_cap = "mass_balance" if saturated else None
    if mw > 1000.0:
        ceiling_applied = pct > 1.0
        if ceiling_applied:
            binding_cap = "mw_1000_ceiling"
        pct = min(pct, 1.0)
        notes.append("MW > 1000 Da: the screening scenario uses a 1% heuristic cap; "
                     "this is not an SCCS default or a validated human absorption ceiling")
    elif mw > _MW_SOFT:
        ceiling_applied = pct > 10.0
        if ceiling_applied:
            binding_cap = "mw_500_rule"
            notes.append(
                f"MW > 500 Da: the flux model estimated {model_pct:g}%; the "
                "Bos-Meinardi 500 Dalton rule replaces it with a conservative "
                "10% default. That is a screening heuristic about what is "
                "likely to cross an intact barrier, not a measured ceiling on "
                "this compound -- a molecule 2 Da lighter is reported from the "
                "model instead. Treat it as a bound to be tested, not a value")
        pct = min(pct, 10.0)

    sccs_default = 10.0 if (mw > _MW_SOFT and (logp < -1.0 or logp > 4.0)) else 50.0
    # `upper_bound_only` is the contract term: True whenever `percent` is a cap
    # rather than a model output. Two conditions qualify, and they are the same
    # kind of statement. Mass-balance saturation means Potts-Guy asked for more
    # flux than the film contains; an MW ceiling means the Bos-Meinardi rule
    # overrode the flux estimate outright. In both cases the number reported is
    # the ceiling, so quoting it as an estimate would be wrong in the same way.
    # `saturated` is kept as the narrower mass-balance flag for diagnostics.
    #
    # Three other modules act on this: `systemic_exposure` below copies it onto
    # the exposure dict, from which `admet.summarise` reads it to suppress CYP
    # gating and `admet_report.overall` reads it to quote the SED as a bound.
    # The key name here and the key they read must stay identical; a previous
    # revision emitted `saturated` and read `mass_balance_capped`, which left
    # all four guards silently inert.
    return {"percent": round(min(pct, 100.0), 2),
            "assessment_status": "MODEL_ONLY",
            "measured_absorption": False,
            #: The flux model's own estimate, before any ceiling. Equal to
            #: `percent` when nothing bound. Keeps "what the model computed"
            #: separable from "the conservative default we substituted".
            "model_percent": model_pct,
            "sccs_default_percent": sccs_default,
            "log_kp": round(log_kp, 2),
            "flux_mg_cm2_h": round(flux, 6),
            "applied_mg_cm2": round(applied_per_area, 4),
            "film_thickness_um": round(film_cm * 1e4, 2),
            "film_mass_basis": "full_pre_rinse_daily_product",
            "absorption_fraction_basis": "full_pre_rinse_donor_film",
            "sed_mass_basis": "retention_adjusted_sccs_table3a_relative_exposure",
            "saturated": saturated,
            "upper_bound_only": saturated or ceiling_applied,
            #: Which cap the reported percent actually is, or None when the
            #: model's own estimate came through. Consumers describe the bound
            #: from this rather than inferring it from `saturated`, which is
            #: true whenever mass balance bound *at any point* and not
            #: necessarily at the end.
            "binding_cap": binding_cap,
            "caps_applied": [c for c in ("mass_balance" if saturated else None,
                                         "mw_ceiling" if ceiling_applied else None)
                             if c],
            "notes": notes}


def _band(sed: float) -> tuple[str, str]:
    for threshold, label, note in _SED_BANDS:
        if sed < threshold:
            return label, note
    return _SED_TOP


def systemic_exposure(mol: Chem.Mol, d: dict, *, applied_percent: float,
                      product: str) -> dict[str, Any]:
    """SCCS Systemic Exposure Dosage and the NOAEL it would demand."""
    if product not in PRODUCT_USE:
        raise ValueError(f"unknown product format {product!r}; "
                         f"expected one of {sorted(PRODUCT_USE)}")
    # Both halves of the exposure contract are checked here, before any
    # arithmetic. The product format was already a closed vocabulary; the use
    # level was not, and an open one produced signed and non-finite SEDs that
    # every consumer downstream treated as a completed assessment.
    applied_percent = validate_applied_percent(applied_percent)
    amount_g, retention, area_cm2, hours = PRODUCT_USE[product]
    # Table 3A includes retention. Re-applying RF or dividing daily grams by
    # the default 60 kg would change the published relative-exposure values.
    a_mg_per_kg = SCCS_RELATIVE_EXPOSURE_MG_KG_DAY[product]

    absorption = dermal_absorption(d["logp"], d["mw"], applied_percent, product)
    sed = a_mg_per_kg * applied_percent / 100.0 * absorption["percent"] / 100.0
    sed_sccs = (a_mg_per_kg * applied_percent / 100.0
                * absorption["sccs_default_percent"] / 100.0)

    band, band_note = _band(sed)
    concern = {"negligible": "low", "low": "low",
               "material": "moderate", "high": "high"}[band]

    # Potts-Guy is a steady-state, infinite-dose model. Applied to a finite film
    # of a small hydrophilic active it predicts a flux larger than the amount
    # applied, which is physically impossible, so absorption is capped at 100%.
    # This is only a bound inside the chosen scenario. Uncapped modeled doses
    # are likewise unmeasured and cannot alone establish systemic toxicity.
    upper_bound_only = bool(absorption["upper_bound_only"])
    if concern == "high":
        concern = "moderate"

    # The cap that fixed the number, from the data rather than re-derived.
    # Inferring it from `saturated` described a 1% MW ceiling as a mass-balance
    # cap whenever the film had also saturated on the way there.
    _CAP_WORDS = {
        "mass_balance": "the mass-balance cap",
        "mw_1000_ceiling": "the screening 1000 Dalton heuristic cap",
        "mw_500_rule": "the Bos-Meinardi 500 Dalton rule",
    }
    cap = _CAP_WORDS.get(absorption.get("binding_cap") or "",
                         "a screening molecular-weight heuristic cap")
    both = len(absorption.get("caps_applied") or []) > 1
    verdict = (f"Modeled SED {sed:.4g} mg/kg bw/day at {applied_percent:g}% in a "
               f"{product.replace('_', ' ')} ({band}); {band_note}. "
               "A model scenario, not measured systemic exposure or a safety finding.")
    if upper_bound_only:
        verdict = (f"SED at most {sed:.4g} mg/kg bw/day at {applied_percent:g}% in a "
                   f"{product.replace('_', ' ')}. Absorption was fixed by {cap}"
                   + (" (the film also saturated on the way there)" if both else "")
                   + ", so this is a capped model scenario, not an estimate "
                   "or a validated human upper bound. Settle it with measured dermal "
                   "absorption (OECD TG 428) before treating it as a finding.")

    return {"concern": concern,
            "assessment_status": "MODEL_ONLY",
            "measured_absorption": False,
            "exposure_policy": EXPOSURE_POLICY,
            "upper_bound_only": upper_bound_only,
            "sed_mg_kg_day": round(sed, 6),
            "sed_sccs_default_mg_kg_day": round(sed_sccs, 6),
            "band": band,
            "applied_percent": applied_percent,
            "product": product,
            "relative_daily_exposure_mg_kg_day": round(a_mg_per_kg, 2),
            "dermal_absorption": absorption,
            "noael_required_for_mos_100": round(sed * MOS_TARGET, 4),
            "assumptions": {
                "body_weight_kg": BODY_WEIGHT_KG,
                "relative_exposure_source": "SCCS/1647/22 Table 3A",
                "relative_exposure_source_url": SCCS_SOURCE_URL,
                "body_weight_role": "reference only; published A uses study participants' body weights",
                "daily_product_g": amount_g,
                "retention_factor": retention,
                "area_cm2": area_cm2,
                "contact_hours": hours,
                "vehicle": "aqueous (Potts-Guy fitting domain)"},
            "verdict": verdict,
            "basis": "SCCS Notes of Guidance SCCS/1647/22 section 3-4.5, with "
                     "Table 3A relative exposure and dermal absorption derived "
                     "from Potts-Guy (1992), substituting RDKit cLogP for measured log Kow, rather "
                     "than taken from the 50% default"}


def cyp_liability(mol: Chem.Mol, d: dict,
                  sed_mg_kg_day: Optional[float] = None) -> dict[str, Any]:
    """CYP inhibition and hERG alerts, with relevance gated on systemic exposure.

    For a leave-on cosmetic these are almost always footnotes. They stop being
    footnotes when SED climbs into the range where a plasma concentration could
    approach an inhibition constant, which is why the gate exists rather than a
    blanket dismissal: a small lipophilic amine applied over the whole body at
    a high use level genuinely can reach a systemic dose.
    """
    alerts = al.match_alerts(mol, al.CYP_COMPILED)
    if photo.is_furocoumarin(mol):
        # Detected in code, not SMARTS, for the reason given in admet_photo.
        # Worth carrying here because it is the same chemotype as the
        # phototoxicity flag: bergamottin is both the berloque-dermatitis agent
        # and the grapefruit-juice drug interaction.
        alerts.insert(0, {
            "name": "furocoumarin_cyp3a4", "severity": "high", "n_matches": 1,
            "mechanism": "furocoumarin; the furan epoxidises and alkylates the "
                         "CYP3A4 apoprotein. This is the grapefruit-juice "
                         "interaction, and it is irreversible - activity returns "
                         "only as new enzyme is synthesised",
            "action": "CYP3A4 time-dependent-inhibition assay if SED is material"})
    herg = _herg(mol, d)
    if herg:
        alerts.append(herg)

    # `_SED_BANDS[1][0]` is this same edge; it was written out again here,
    # with nothing holding the two together.
    relevant = sed_mg_kg_day is not None and sed_mg_kg_day >= _SED_BANDS[1][0]
    order = ("critical", "high", "moderate", "low", "context")
    worst = next((lv for lv in order if any(a["severity"] == lv for a in alerts)),
                 "none")
    if not alerts:
        verdict = "no CYP or hERG structural alert"
        concern = "none"
    elif worst == "low" and not herg:
        # Only note-level alerts fired. Escalating these on exposure alone would
        # make niacinamide, whose pyridine is the entire "alert", look like a
        # drug-interaction risk.
        verdict = ("note-level CYP alerts only (" +
                   ", ".join(a["name"] for a in alerts) +
                   "); nothing here would change a formulation decision")
        concern = "low"
    elif relevant:
        verdict = (f"CYP/hERG alerts present AND SED is {sed_mg_kg_day:.3g} "
                   "mg/kg bw/day, which is high enough to matter; run a CYP "
                   "inhibition panel and, for the hERG alert, a patch-clamp assay")
        concern = "moderate"
    elif sed_mg_kg_day is None:
        verdict = ("CYP/hERG alerts present; systemic exposure relevance is "
                   "NOT_EVALUATED because no exposure estimate was supplied. "
                   "An upper-bound-only absorption model is not an exposure "
                   "estimate and does not establish negligible exposure. "
                   "Retain these alerts and evaluate dermal absorption before "
                   "interpreting their systemic relevance")
        concern = "context"
    else:
        verdict = (f"CYP/hERG alerts present; supplied SED {sed_mg_kg_day:.3g} "
                   f"mg/kg bw/day is below the {_SED_BANDS[1][0]:g} mg/kg "
                   "screening threshold, so the alerts are recorded rather "
                   "than escalated by this rule. This threshold does not "
                   "establish safety; revisit if use level or product format "
                   "changes")
        concern = "context"

    return {"concern": concern, "alerts": alerts,
            "sed_mg_kg_day": sed_mg_kg_day, "gated_as_relevant": relevant,
            "verdict": verdict,
            "basis": "structural alerts only; in-vitro CYP inhibition and hERG "
                     "patch clamp are the reference methods"}


def _herg(mol: Chem.Mol, d: dict) -> Optional[dict]:
    """The classic hERG pharmacophore: a basic centre, lipophilicity, aromatics.

    A protonatable nitrogen makes the cation-pi contact with Tyr652/Phe656 in
    the channel pore, and lipophilic aromatic bulk supplies the rest. The rule
    of thumb is cLogP above ~3.7 with a basic amine; it is coarse, and it is
    also the reason most hERG binders look alike.
    """
    basic = Chem.MolFromSmarts(
        "[NX3;H0,H1,H2;!$(N[CX3]=[OX1,SX1,NX2]);!$(N[SX4]);!$(Na);!$(N=*);!$(NN);"
        "!$(N[OX2])]")
    if basic is None or not mol.GetSubstructMatches(basic):
        return None
    if d["logp"] < 3.7 or d["n_aromatic_rings"] < 2 or not 250 <= d["mw"] <= 700:
        return None
    return {"name": "herg_pharmacophore", "severity": "low", "n_matches": 1,
            "mechanism": f"basic nitrogen with logP {d['logp']:.1f} and "
            f"{int(d['n_aromatic_rings'])} aromatic rings; the cation-pi "
            "pharmacophore that binds Tyr652/Phe656 in the hERG pore",
            "action": "hERG patch clamp, but only if SED is material"}


def cmr_alerts(mol: Chem.Mol) -> dict[str, Any]:
    """CMR screening features and bounded reference matches, not classification."""
    return admet_cmr.assess(mol)
