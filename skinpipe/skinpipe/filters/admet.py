"""Topical ADMET and toxicity profiling for cosmetic ingredient discovery.

Why this is not an oral-drug ADMET module
-----------------------------------------
The exposure being modelled is a leave-on product spread on skin, often on the
face and forearms, every day, for decades, by a population that includes
pregnant women and children. Very little of what an oral drug filter cares
about survives that change of context:

*   Oral bioavailability is irrelevant. The question is how much of the applied
    dose crosses the stratum corneum into circulation - a flux problem, not a
    first-pass problem - which `admet_exposure` answers the way the SCCS Notes
    of Guidance does, as a Systemic Exposure Dosage in mg/kg bw/day.
*   CYP inhibition and hERG, which dominate a drug ADMET panel, are close to
    meaningless when SED is micrograms per kilogram, and real the moment it is
    not. They are reported gated on the computed SED, not dropped.
*   Phototoxicity is a relevant endpoint for a sun-exposed topical. Structural
    screening can return `critical` in several domains, including CMR and
    genotoxicity; that priority is not a measured toxicity or legal result.

Domains, and the reference method each one stands in for
--------------------------------------------------------
    phototoxicity     OECD TG 432 (3T3 NRU), gated by the ICH S10 / SCCS
                      light-absorption screen; engine in `admet_photo`
    skin_irritation   OECD TG 439 (reconstructed human epidermis);
                      engine in `admet_irritation`
    genotoxicity      OECD TG 471 (Ames) and TG 487 (micronucleus), alerts
                      after Kazius, McGuire & Bursi (2005) J Med Chem 48:312
                      and the ICH M7(R2) structural-alert concept
    systemic_exposure SCCS Notes of Guidance SCCS/1647/22, SED and MoS
    cyp_liability     in-vitro CYP panel and hERG patch clamp
    cmr_alerts        structural features and limited reference-connectivity
                      warnings; experimental toxicity and legal identity unestablished

Skin sensitisation is deliberately absent: `filters.safety` owns that endpoint
and this module cross-references it from the irritation domain rather than
duplicating it. Permeation is imported from `filters.topical`, not rebuilt.

Nothing here is a regulatory result. Every function returns a graded,
mechanism-annotated prioritisation signal in the style of `filters.safety`:
never a bare pass/fail, and never a silent drop of a whole chemotype. Every
calculation runs offline from the connection table alone; no network call is
made, attempted or required.
"""
from __future__ import annotations

import logging
from typing import Any, Optional

from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors

from ..types import ADMETConcern, ADMETProfile
from . import admet_alerts as al
from . import admet_exposure as expo
from . import admet_irritation as irr
from . import admet_photo as photo
from . import admet_report as report

RDLogger.DisableLog("rdApp.*")
log = logging.getLogger(__name__)

__all__ = ["phototoxicity", "skin_irritation", "genotoxicity",
           "systemic_exposure", "cyp_liability", "cmr_alerts", "summarise"]

def _mol(smiles: str) -> Optional[Chem.Mol]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        log.warning("admet: SMILES could not be parsed: %r", smiles)
        return None
    # RDKit parses "" into a perfectly valid molecule with no atoms, so an
    # empty or hydrogen-only SMILES walked the whole assessment and came back
    # `assessed=True, hazard_severity="none"` -- "we looked and found no
    # structural liability" about a structure that is not there. Every alert
    # below is a substructure search, and nothing matches nothing.
    # `catalog.screen` already guards its *query* SMILES this way and says why;
    # the compounds being screened had no such guard.
    if mol.GetNumHeavyAtoms() == 0:
        log.warning("admet: SMILES has no heavy atoms, nothing to assess: %r",
                    smiles)
        return None
    return mol


def _descriptors(mol: Chem.Mol) -> dict[str, float]:
    return {"mw": Descriptors.MolWt(mol),
            "logp": Crippen.MolLogP(mol),
            "tpsa": rdMolDescriptors.CalcTPSA(mol),
            "hbd": rdMolDescriptors.CalcNumHBD(mol),
            "n_aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(mol)}


def _worst(severities: list[str]) -> str:
    for level in ADMETProfile.SEVERITY_ORDER:
        if level in severities:
            return level
    return "none"


def _unparseable(domain: str) -> dict[str, Any]:
    return {"concern": "unknown", "alerts": [],
            "verdict": f"SMILES could not be parsed; {domain} not assessed"}


# --------------------------------------------------------- 1. phototoxicity

def phototoxicity(smiles: str) -> dict[str, Any]:
    """Screen against the ICH S10 / SCCS light-absorption criterion.

    The rule is an exclusion, and only an exclusion. Stated as the guideline
    and `admet_photo` state it: a substance that does *not* absorb between 290
    and 700 nm with a molar extinction coefficient of at least
    1000 L mol-1 cm-1 is not expected to be directly phototoxic and needs no
    further photosafety testing. Absorbing above that threshold establishes
    nothing; it removes the exclusion, which reopens the question and leaves
    OECD TG 432 to answer it. Absorption is necessary for direct
    phototoxicity, never sufficient, and reading it as sufficient is affirming
    the consequent
    (Henry et al. 2009 Regul Toxicol Pharmacol 55:e1; ICH S10 (2013);
    SCCS Notes of Guidance SCCS/1647/22 section 3-3.9).

    A molar extinction coefficient is an experimental quantity and cannot be
    computed from a connection table, so this returns an explicitly labelled
    structural proxy - conjugation extent plus membership of a documented
    phototoxic chemotype - as a band, never a number. The recommendation is
    always "measure the spectrum and run TG 432", never "this is phototoxic".
    """
    mol = _mol(smiles)
    if mol is None:
        return _unparseable("phototoxicity")

    score, feats, reasons = photo.extinction_proxy(mol)
    concern, epsilon = photo.band(score)
    chemotypes = al.match_alerts(mol, al.PHOTO_COMPILED)

    if feats["furocoumarin"]:
        chemotypes.insert(0, {
            "name": "furocoumarin", "severity": "critical", "n_matches": 1,
            "mechanism": "psoralen/bergapten class; intercalates DNA and forms "
                         "[2+2] cyclobutane adducts with thymine under UVA, then "
                         "interstrand cross-links. Phototoxic and photogenotoxic, "
                         "which is why bergapten in expressed citrus oils is "
                         "capped at 1 ppm in EU leave-on cosmetics",
            "action": "quantify furocoumarins in the extract; EU Annex II/III"})
    if feats["tetrapyrrole"]:
        chemotypes.insert(0, {
            "name": "porphyrin", "severity": "critical", "n_matches": 1,
            "mechanism": "tetrapyrrole; Soret band near 400 nm with epsilon of "
                         "order 1e5 and a near-unit singlet-oxygen yield. This is "
                         "the photodynamic-therapy chromophore",
            "action": "not usable in a leave-on product"})

    # A documented chemotype outranks the proxy in both directions: the proxy is
    # an inference about absorption, the chemotype is an observation about people.
    worst_chemotype = _worst([c["severity"] for c in chemotypes])
    if worst_chemotype == "critical":
        concern = "critical"
    elif worst_chemotype == "high":
        concern = "high"
    elif worst_chemotype == "moderate" and concern in ("none", "low"):
        concern = "moderate"

    verdict = {
        "critical": "documented photosensitising chemotype; do not take into a "
                    "sun-exposed leave-on product without regulatory clearance",
        "high": "absorbs in the 290-700 nm sunlight window; treat as a "
                "phototoxicity concern and run OECD TG 432 (3T3 NRU) before any "
                "sun-exposed use",
        "moderate": "borderline chromophore; measure the UV/Vis spectrum in the "
                    "intended vehicle and apply the 1000 L/mol/cm rule before "
                    "committing to TG 432",
        "low": "weak chromophore; likely exempt from photosafety testing under "
               "ICH S10, but confirm with a measured spectrum",
    }.get(concern, "no chromophore reaching 290 nm; no photosafety testing indicated")

    return {"concern": concern, "proxy_score": score, "estimated_mec": epsilon,
            "chemotypes": chemotypes, "chromophore": feats, "reasons": reasons,
            "verdict": verdict,
            "basis": "structure-based proxy for the ICH S10 / SCCS 1000 L/mol/cm "
                     "screen; NOT a measured spectrum"}


# ------------------------------------------------------- 2. skin irritation

def skin_irritation(smiles: str) -> dict[str, Any]:
    """Irritation and corrosion triage; see `admet_irritation` for the method.

    Acid/base class, amphiphile architecture and lipophilicity, plus a
    cross-reference to the electrophiles `filters.safety` already tracks.
    OECD TG 439 is the reference method and this does not replace it.
    """
    mol = _mol(smiles)
    if mol is None:
        return _unparseable("skin irritation")
    return irr.assess(mol, smiles, _descriptors(mol), _worst)


# --------------------------------------------------------- 3. genotoxicity

def genotoxicity(smiles: str) -> dict[str, Any]:
    """Ames/mutagenicity structural alerts, after Kazius-Bursi and ICH M7.

    Kazius, McGuire & Bursi (2005) J Med Chem 48:312 derived general
    toxicophores from 4337 Ames-tested compounds; ICH M7(R2) formalised the
    same idea as the alert that puts an impurity in the cohort-of-concern.
    Every alert names its activation route to a DNA-reactive species, because
    a bare alert name is not actionable and a chemist can argue with a
    mechanism. Fused polycyclic aromatics are detected by ring analysis, not
    SMARTS: the class is topological and no pattern set enumerates it.
    """
    mol = _mol(smiles)
    if mol is None:
        return _unparseable("genotoxicity")

    alerts = al.match_alerts(mol, al.GENOTOX_COMPILED)

    n_fused, n_fused_arom = photo.max_fused_rings(mol)
    if n_fused_arom >= 3 and _all_carbocyclic(mol):
        alerts.append({
            "name": "polycyclic_aromatic", "severity": "high",
            "n_matches": 1,
            "mechanism": f"{n_fused_arom} fused aromatic carbocycles; CYP1A1/1B1 "
            "epoxidation gives a bay- or fjord-region diol epoxide that alkylates "
            "guanine N2. This is the benzo[a]pyrene mechanism and the reason "
            "PAH content is regulated in mineral-oil and carbon-black grades",
            "action": "Ames (OECD TG 471) plus a PAH specification on the raw material"})

    concern = _worst([a["severity"] for a in alerts])
    verdict = {
        "critical": "cohort-of-concern chemotype; not acceptable in a cosmetic "
                    "at any level without removal or structural change",
        "high": "DNA-reactive structural alert; an Ames test (OECD TG 471) is "
                "required before this goes near a product",
        "moderate": "indirect or weak genotoxicity alert; the in-vitro "
                    "micronucleus test (OECD TG 487) reads this better than Ames",
        "low": "weak alert; note it and move on",
    }.get(concern, "no Ames structural alert")

    return {"concern": concern, "alerts": alerts,
            "n_fused_aromatic_rings": n_fused_arom, "verdict": verdict,
            "basis": "Kazius/Bursi toxicophores and ICH M7(R2) structural alerts; "
                     "structural alerts have roughly 80% sensitivity against the "
                     "Ames test and are not a substitute for running it"}


def _all_carbocyclic(mol: Chem.Mol) -> bool:
    """True when some fused aromatic system is made only of carbon.

    Guards the PAH alert against heteroaromatics: a pteridine or a flavonoid is
    tricyclic and aromatic but does not form a bay-region diol epoxide.
    """
    rings = mol.GetRingInfo().AtomRings()
    for system in photo.ring_systems(mol):
        members = [r for r in rings if set(r) <= system]
        aromatic = [r for r in members
                    if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in r)]
        if len(aromatic) >= 3 and all(mol.GetAtomWithIdx(i).GetAtomicNum() == 6
                                      for r in aromatic for i in r):
            return True
    return False


# ----------------------------------------------------- 4. systemic exposure

def systemic_exposure(smiles: str, *, applied_percent: float = 1.0,
                      product: str = "body_lotion") -> dict[str, Any]:
    """SCCS Systemic Exposure Dosage in mg/kg bw/day; see `admet_exposure`.

    `applied_percent` is exposed so a formulator can test a use level rather
    than receive one. `product` selects the SCCS use scenario; body lotion is
    the default because it is the worst case for a leave-on active (largest
    area, largest daily amount, retention factor 1).

    The use level is validated before the structure is parsed, so an impossible
    concentration is refused the same way whatever molecule it was asked about:
    an unparseable SMILES must not be what decides whether -1% is acceptable.
    """
    expo.validate_applied_percent(applied_percent)
    mol = _mol(smiles)
    if mol is None:
        return _unparseable("systemic exposure")
    return expo.systemic_exposure(mol, _descriptors(mol),
                                  applied_percent=applied_percent, product=product)


# ----------------------------------------------------------- 5. CYP / hERG

def cyp_liability(smiles: str, *,
                  sed_mg_kg_day: Optional[float] = None) -> dict[str, Any]:
    """CYP inhibition and hERG alerts, gated on systemic exposure.

    Low priority for a topical with negligible SED, and explicitly labelled as
    such - but pass `sed_mg_kg_day` from `systemic_exposure` and the same
    alerts are reported as actionable once exposure is material.
    """
    mol = _mol(smiles)
    if mol is None:
        return _unparseable("CYP liability")
    return expo.cyp_liability(mol, _descriptors(mol), sed_mg_kg_day)


def cmr_alerts(smiles: str) -> dict[str, Any]:
    """Carcinogenic / mutagenic / reprotoxic chemotypes; see `admet_exposure`."""
    mol = _mol(smiles)
    if mol is None:
        return _unparseable("CMR")
    return expo.cmr_alerts(mol)


# --------------------------------------------------------------- summarise














MECHANISM_ESCALATION: dict[str, dict[str, str]] = {
    "TYR": {"leukoderma": "high"},
    "TYRP1": {"leukoderma": "high"},
}


def _escalate_for_target(concerns: list[ADMETConcern],
                         target_gene: Optional[str]) -> list[str]:
    """Raise alerts the target itself would activate. Returns what was raised."""
    rules = MECHANISM_ESCALATION.get((target_gene or "").upper())
    if not rules:
        return []
    raised: list[str] = []
    for i, c in enumerate(concerns):
        for motif, severity in rules.items():
            if motif in c.name and c.severity != severity:
                concerns[i] = ADMETConcern(
                    domain=c.domain, name=c.name, severity=severity,
                    mechanism=(f"{c.mechanism} At {target_gene} this is not a "
                               "background alert: the target enzyme is what makes "
                               "the toxic quinone, so the programme's own "
                               "pharmacology drives the hazard."),
                    action=c.action)
                raised.append(c.name)
    return raised


def summarise(smiles: str, *, applied_percent: Optional[float] = None,
              product: str = "body_lotion",
              ligand_id: Optional[str] = None,
              target_gene: Optional[str] = None) -> ADMETProfile:
    """Full topical ADMET profile with one plain-language call at the end.

    `applied_percent` defaults to 1%, a typical leave-on active level, and is
    recorded in the profile so the SED is never read without its use level.

    `target_gene` lets an alert be judged against the programme it belongs to.
    A 4-substituted phenol is a footnote in most screens and the central risk in
    a tyrosinase screen, because tyrosinase is the enzyme that converts it to a
    melanocyte-killing quinone. See `MECHANISM_ESCALATION`.
    """
    # `None` is this function's documented "use the default" and nothing else
    # is: `float(applied_percent)` used to accept -1, 101 and NaN, and the
    # profile that came back carried a negative or non-finite SED with
    # `assessed=True`. Refused before the structure is parsed, and refused with
    # the producer's own validator so the rule cannot differ from the CLI's.
    pct = (1.0 if applied_percent is None
           else expo.validate_applied_percent(applied_percent))
    if _mol(smiles) is None:
        return ADMETProfile(smiles=smiles, ligand_id=ligand_id,
                            overall="SMILES could not be parsed; no assessment made")

    photo_block = phototoxicity(smiles)
    irritation = skin_irritation(smiles)
    geno = genotoxicity(smiles)
    exposure = systemic_exposure(smiles, applied_percent=pct, product=product)
    # Both capped and uncapped cLogP-derived SEDs are unmeasured scenarios.
    # Neither alone establishes circulating exposure or requires a CYP panel.
    _sed_for_gating = (None if exposure.get("upper_bound_only") or not exposure.get("measured_absorption")
                       else exposure.get("sed_mg_kg_day"))
    cyp = cyp_liability(smiles, sed_mg_kg_day=_sed_for_gating)
    cmr = cmr_alerts(smiles)

    concerns: list[ADMETConcern] = []
    report.collect("phototoxicity", photo_block, concerns, "uv_vis_chromophore")
    report.collect("irritation", irritation, concerns, "ph_shift")
    report.collect("genotoxicity", geno, concerns)
    report.collect("cyp", cyp, concerns)
    report.collect("cmr", cmr, concerns)
    # Use the exposure block's own `concern`, not its raw band. When Potts-Guy
    # runs past the mass-balance cap the band is an upper bound from a model
    # outside its domain, and the block downgrades the severity accordingly;
    # reading the band here would undo that and flag niacinamide as high.
    if exposure.get("concern") in ("moderate", "high"):
        concerns.append(ADMETConcern(
            domain="systemic", name="systemic_exposure",
            severity=exposure["concern"],
            mechanism=exposure["verdict"],
            action=f"a repeat-dose NOAEL of at least "
                   f"{exposure['noael_required_for_mos_100']:.3g} mg/kg bw/day is "
                   f"needed to clear a margin of safety of 100"))

    escalated = _escalate_for_target(concerns, target_gene)
    report.rank(concerns)

    verdicts = {"phototoxicity": photo_block["verdict"],
                "skin_irritation": irritation["verdict"],
                "genotoxicity": geno["verdict"],
                "systemic_exposure": exposure["verdict"],
                "cyp_liability": cyp["verdict"],
                "cmr": cmr["verdict"]}

    return ADMETProfile(
        smiles=smiles, ligand_id=ligand_id, phototoxicity=photo_block,
        skin_irritation=irritation, genotoxicity=geno,
        systemic_exposure=exposure, cyp_liability=cyp, cmr=cmr,
        concerns=concerns, verdicts=verdicts,
        overall=report.overall(concerns, photo_block, exposure, pct, product)
        + (f" Escalated for {target_gene}: {', '.join(escalated)} is the mechanism "
           "of harm at this target, not a background alert." if escalated else ""))
