"""Skin irritation and corrosion triage for `filters.admet`.

The reference method is OECD TG 439 (in-vitro skin irritation on a
reconstructed human epidermis model), with OECD TG 431 sitting above it for
corrosion. This module is not that and does not substitute for it. It is a
prioritisation signal built from three things structure does predict:

*   Acid and base strength, estimated by functional-group class. Skin's acid
    mantle sits at pH 4.5-5.5; a substance that drives the surface outside
    roughly pH 3-9 disrupts the lipid-processing enzymes that build the
    barrier. The pKa here is a class default, not a computed value - a
    carboxylic acid on an electron-poor ring is a full unit stronger than the
    4.5 assumed for the class, and this module will not notice.
*   Amphiphile architecture. Surfactant irritation is structural and
    well-ordered: cationic > anionic > non-ionic, scaling with tail length.
    Benzalkonium chloride is the TG 439 positive control precisely because a
    C12 tail on a permanent cation is a reliable irritant.
*   Gross lipophilicity, which decides whether an irritant reaches living
    tissue at all or sits on the surface.

Reactivity is not re-derived here. `filters.safety` already carries the
electrophilic alert set for sensitisation; this module imports it and reads it
for the irritation endpoint, because the chemistry that haptenates keratin also
denatures it at higher concentration.
"""
from __future__ import annotations

import logging
from typing import Any, Optional

from rdkit import Chem

from . import admet_alerts as al

log = logging.getLogger(__name__)

# Acid/base strength by functional-group class. pKa values are class
# representatives, not computed for the molecule in hand.
# (name, SMARTS, direction, representative pKa, note)
ACID_BASE_CLASSES: tuple[tuple[str, str, str, Optional[float], str], ...] = (
    ("sulfonic_acid", "[#6][SX4](=[OX1])(=[OX1])[OX2H1,OX1H0-]", "acid", 0.5,
     "fully dissociated at any formulation pH; corrosive as the free acid"),
    ("sulfate_ester_acid", "[#6][OX2][SX4](=[OX1])(=[OX1])[OX2H1,OX1H0-]",
     "acid", 1.0, "alkyl sulfate; strong acid head group, and an amphiphile"),
    # The SMARTS must require a free P-OH or P-O(-). Written without that last
    # constraint it matched a fully esterified phosphodiester, so a
    # phosphatidylcholine -- the commonest cosmetic emulsifier -- came back
    # "strongly acidic" on the first dissociation of phosphoric acid, a proton
    # it does not have.
    ("phosphonic_phosphoric",
     "[PX4](=[OX1])([OX2H1,OX1H0-])[OX2H1,OX1H0-]", "acid", 2.0,
     "first dissociation is strong; buffers to a low pH"),
    ("alpha_halo_acid", "[CX4]([F,Cl,Br,I])[CX3](=[OX1])[OX2H1,OX1H0-]", "acid",
     2.5, "inductive stabilisation; trichloroacetic-acid class, a peeling agent"),
    ("carboxylic_acid", "[CX3](=[OX1])[OX2H1,OX1H0-]", "acid", 4.5,
     "weak acid; irritation is pH- and concentration-driven, not intrinsic "
     "(azelaic and glycolic acid are both marketed leave-on actives)"),
    ("phenol", "[OX2H][c]", "acid", 10.0,
     "very weak acid; the irritation risk is oxidative, not pH"),
    ("guanidine", "[NX3][CX3](=[NX2])[NX3]", "base", 13.5,
     "protonated at every skin pH; strongly alkaline as the free base"),
    ("amidine", "[NX3][CX3]=[NX2;!$([NX2][NX3])]", "base", 12.0,
     "strongly basic; disrupts the acid mantle"),
    ("aliphatic_amine", "[NX3;H1,H2,H0;!$(N[CX3]=[OX1,SX1,NX2]);!$(N[SX4]);"
                        "!$(N[PX4]);!$(Na);!$(N[OX2]);!$(N=*);!$(NN)]", "base", 10.3,
     "raises formulation pH; the free base also partitions into the "
     "stratum-corneum lipid and fluidises it"),
    # `None`, not a number. This row used to carry 99.0 as a sentinel meaning
    # "never deprotonates", and `ph_class` read it as a measured basicity:
    # 99.0 >= 12.0, so every permanently cationic compound came back "strongly
    # basic" with a sentence about depilatory burns. Choline, betaine and
    # L-carnitine -- humectants in thousands of marketed products -- reached
    # `hazard_severity = high` on that alone, with no alert fired, and
    # `--admet-gate high` deleted them with the reason `ph_shift`. A quaternary
    # ammonium has no conjugate base and cannot titrate; benzalkonium solutions
    # are pH 6-8. The real irritation signal for this class is the cationic
    # surfactant alert below, which fires on a C>=8 tail and is unaffected.
    ("quaternary_ammonium", "[NX4+;!$([NX4+][OX1-])]", "permanent_cation", None,
     "permanently cationic; binds anionic keratin and is the most irritant "
     "surfactant class (benzalkonium is the positive control in TG 439). "
     "It cannot shift pH: there is no proton to lose."),
)

ACID_BASE_COMPILED = tuple(
    (name, Chem.MolFromSmarts(sma), direction, pka, note)
    for name, sma, direction, pka, note in ACID_BASE_CLASSES)


def ph_class(acid_pka: Optional[float],
             base_pka: Optional[float]) -> tuple[str, str]:
    """Coarse pH-shift class from the strongest ionisable group present."""
    if acid_pka is not None and acid_pka <= 2.5:
        return ("strongly acidic",
                f"strongest acidic group has a class pKa near {acid_pka}; as the "
                "free acid this drives skin pH well below the acid mantle (4.5-5.5)")
    if base_pka is not None and base_pka >= 12.0:
        return ("strongly basic",
                f"class pKa near {base_pka}; alkaline pH swells corneocytes and "
                "strips the lipid lamellae, the mechanism behind depilatory burns")
    if acid_pka is not None and acid_pka <= 5.0:
        return ("weakly acidic",
                "weak acid; azelaic and glycolic acid are both marketed leave-on "
                "actives, so this is a pH and concentration question, not a hazard")
    if base_pka is not None and base_pka >= 9.0:
        return ("weakly basic", "buffer to skin pH and the free base is minimised")
    return ("near-neutral", "no ionisable group strong enough to shift skin pH")


def longest_alkyl_chain(mol: Chem.Mol) -> int:
    """Longest run of acyclic sp3 carbons: the hydrophobic tail of a surfactant."""
    keep = {a.GetIdx() for a in mol.GetAtoms()
            if a.GetAtomicNum() == 6 and not a.IsInRing()
            and a.GetHybridization() == Chem.HybridizationType.SP3}
    best = 0

    def walk(node: int, seen: frozenset[int]) -> None:
        nonlocal best
        best = max(best, len(seen))
        for nb in mol.GetAtomWithIdx(node).GetNeighbors():
            if nb.GetIdx() in keep and nb.GetIdx() not in seen:
                walk(nb.GetIdx(), seen | {nb.GetIdx()})

    for start in keep:
        walk(start, frozenset({start}))
    return best


def _amphiphile(chain: int, bases: list[dict], acids: list[dict],
                d: dict) -> Optional[dict]:
    if chain < 8:
        return None
    # Amphoterics first, because they satisfy the cationic test as well and the
    # branch order alone decided which one answered. A betaine carries a
    # permanent quaternary nitrogen *and* its own anion, which is the whole
    # point of the class: the internal charge pair is why cocamidopropyl
    # betaine is the mild substitute for a cationic surfactant rather than an
    # instance of one. Tested below the cationic branch it came back
    # `cationic_surfactant` at high severity, carrying benzalkonium's TG 439
    # mechanism text -- while the same call had already collected
    # `carboxylic_acid` in `acids` and so knew the anion was there.
    if (any(b["name"] == "quaternary_ammonium" for b in bases)
            and any(a["name"] in ("sulfonic_acid", "sulfate_ester_acid",
                                  "phosphonic_phosphoric", "carboxylic_acid")
                    for a in acids)):
        return {"name": "amphoteric_surfactant", "severity": "moderate", "n_matches": 1,
                "mechanism": f"C{chain} tail on a zwitterionic head; the internal "
                "charge pair blunts the keratin binding that makes a permanent "
                "cation harsh. Milder than either parent class, not inert",
                "action": "OECD TG 439; mild but confirm at the use level"}
    if any(b["name"] == "quaternary_ammonium" for b in bases):
        return {"name": "cationic_surfactant", "severity": "high", "n_matches": 1,
                "mechanism": f"C{chain} tail on a permanent cation; binds anionic "
                "keratin and solubilises the lipid lamellae. Benzalkonium "
                "chloride is the TG 439 positive control for exactly this reason",
                "action": "OECD TG 439; expect a very low ceiling concentration"}
    if any(a["name"] in ("sulfonic_acid", "sulfate_ester_acid",
                         "phosphonic_phosphoric", "carboxylic_acid") for a in acids):
        return {"name": "anionic_surfactant", "severity": "moderate", "n_matches": 1,
                "mechanism": f"C{chain} tail on an anionic head; denatures stratum "
                "corneum protein and extracts lipid (the sodium lauryl sulfate "
                "mechanism). Rinse-off formats tolerate this, leave-on do not",
                "action": "OECD TG 439 and a cumulative-irritation patch test"}
    if d["hbd"] >= 2 or d["tpsa"] >= 60:
        return {"name": "nonionic_amphiphile", "severity": "low", "n_matches": 1,
                "mechanism": f"C{chain} tail on a polar head; non-ionic "
                "surfactants are the mildest class but still fluidise lipid",
                "action": "note only"}
    return None


def _reactivity_crossref(smiles: str) -> Optional[dict]:
    """Read `filters.safety` for a different endpoint instead of duplicating it."""
    try:
        from .safety import sensitisation_alerts
    except Exception:                                   # pragma: no cover
        log.warning("admet: filters.safety unavailable; reactivity cross-check skipped")
        return None
    hits = [a for a in sensitisation_alerts(smiles) if a["severity"] == "high"]
    if not hits:
        return None
    return {"name": "protein_reactive", "severity": "moderate", "n_matches": len(hits),
            "mechanism": "covalently protein-reactive ("
                         + ", ".join(sorted({h["name"] for h in hits}))
                         + "); the electrophilicity that haptenates keratin also "
                           "denatures it directly at higher concentration",
            "action": "filters.safety gives the sensitisation reading; TG 439 "
                      "gives the irritation half"}


def assess(mol: Chem.Mol, smiles: str, d: dict, worst: Any) -> dict[str, Any]:
    """Full irritation triage. `worst` is the severity-ranking helper."""
    alerts = al.match_alerts(mol, al.IRRITATION_COMPILED)

    acids: list[dict] = []
    bases: list[dict] = []
    cations: list[dict] = []      # permanently charged: no pKa, no pH shift
    for name, patt, direction, pka, note in ACID_BASE_COMPILED:
        if patt is None or not mol.GetSubstructMatches(patt):
            continue
        row = {"name": name, "class_pka": pka, "note": note}
        if direction == "acid":
            acids.append(row)
        elif direction == "permanent_cation":
            cations.append(row)
        else:
            bases.append(row)

    # A group with no pKa cannot be the strongest anything. Guarding on `is not
    # None` rather than on the direction alone means a future row that forgets
    # to set one degrades to "not considered" instead of raising inside max().
    strongest_acid = min((a["class_pka"] for a in acids
                          if a["class_pka"] is not None), default=None)
    strongest_base = max((b["class_pka"] for b in bases
                          if b["class_pka"] is not None), default=None)
    ph, ph_note = ph_class(strongest_acid, strongest_base)

    chain = longest_alkyl_chain(mol)
    # Cations are passed alongside the bases: they are excluded from the pH
    # arm, not from the surfactant arm, which is where their real hazard is.
    surfactant = _amphiphile(chain, bases + cations, acids, d)
    if surfactant:
        alerts.append(surfactant)

    if d["logp"] > 5.0:
        alerts.append({"name": "very_lipophilic", "severity": "low", "n_matches": 1,
                       "mechanism": f"logP {d['logp']:.1f}; partitions into the "
                       "stratum-corneum lipid and stays there. Acute irritation is "
                       "low, but follicular accumulation and comedogenicity are not",
                       "action": "a comedogenicity screen reads this better than TG 439"})
    if d["logp"] < -1.0 and ph.startswith("strongly"):
        alerts.append({"name": "polar_ionisable", "severity": "moderate", "n_matches": 1,
                       "mechanism": f"logP {d['logp']:.1f} with a {ph} group; it "
                       "barely permeates, so the whole applied dose stays on the "
                       "surface where the pH shift acts on the acid mantle",
                       "action": "buffer the formulation, then confirm with TG 439"})

    reactive = _reactivity_crossref(smiles)
    if reactive:
        alerts.append(reactive)

    concern = worst([a["severity"] for a in alerts]
                    + (["high"] if ph.startswith("strongly") else []))
    verdict = {
        "critical": "corrosive chemotype; OECD TG 431 before anything else",
        "high": "expect irritation at a functional use level; TG 439 (RhE) plus a "
                "human patch test should set the ceiling concentration",
        "moderate": "possible irritant; confirm at the intended use level with "
                    "TG 439 rather than assuming in either direction",
        "low": "no strong irritation alert; the usual formulation controls apply",
    }.get(concern, "no structural irritation alert")

    return {"concern": concern, "alerts": alerts, "ph_class": ph, "ph_note": ph_note,
            "acid_groups": acids, "base_groups": bases,
            "longest_alkyl_chain": chain, "logp": round(d["logp"], 2),
            "verdict": verdict,
            "basis": "structural and physicochemical triage; OECD TG 439 "
                     "(reconstructed human epidermis) is the reference method and "
                     "this does not substitute for it"}
