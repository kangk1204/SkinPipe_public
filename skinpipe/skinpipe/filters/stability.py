"""Formulation stability and colour alerts.

A compound that is potent, permeant and safe can still be unusable if it
degrades in the jar. In cosmetics this is not a footnote: it is the reason
several well-supported actives underperform in finished product relative to
their published cell data.

The three failure modes that matter for a leave-on topical product:

*   Photolability. Extended conjugation absorbs in the near-UV, so the actives
    most often paired with sun exposure are the ones least able to survive it.
    Retinoids isomerise and degrade; resveratrol photoisomerises trans->cis,
    and the cis form is much less active.
*   Oxidation. Polyphenols, catechols, enediols and thiols oxidise in air,
    usually browning as they go. Ascorbic acid is the canonical example.
*   Hydrolysis. Esters, lactones and glycosides cleave at formulation pH,
    which is sometimes intended (a prodrug ester releasing the active) and
    sometimes just loss of material.

Alerts here are reasons to specify packaging, pH, antioxidant, or a stabilised
derivative. They are not reasons to drop a compound.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass

from rdkit import Chem, RDLogger

from . import admet_photo as photo
from .safety import assessable_molecule

RDLogger.DisableLog("rdApp.*")
log = logging.getLogger(__name__)


@dataclass(frozen=True)
class StabilityAlert:
    name: str
    mode: str          # photo | oxidation | hydrolysis | colour
    smarts: str
    mitigation: str


_ALERTS: tuple[StabilityAlert, ...] = (
    StabilityAlert("polyene", "photo", "[CX3]=[CX3][CX3]=[CX3][CX3]=[CX3]",
                   "opaque airless packaging, antioxidant, or an ester prodrug (retinyl/HPR)"),
    StabilityAlert("stilbene", "photo", "[c][CX3]=[CX3][c]",
                   "trans->cis photoisomerisation; UV-opaque packaging or encapsulation"),
    StabilityAlert("aryl_enone", "photo", "[c][CX3]=[CX3][CX3]=[OX1]",
                   "chalcone/curcuminoid motif; photodegrades, use opaque packaging"),
    StabilityAlert("nitroaromatic", "photo", "[c][N+](=[OX1])[O-]", "amber glass"),
    StabilityAlert("azo", "photo", "[#6][NX2]=[NX2][#6]", "amber glass"),

    StabilityAlert("catechol", "oxidation", "[OX2H][c]1[c]([OX2H])[c][c][c][c]1",
                   "chelate trace metals (EDTA/phytate), degas, add a compatible antioxidant"),
    StabilityAlert("hydroquinone", "oxidation", "[OX2H][c]1[c][c][c]([OX2H])[c][c]1",
                   "oxygen-barrier packaging and a reducing agent; expect browning"),
    StabilityAlert("pyrogallol", "oxidation", "[OX2H][c]1[c]([OX2H])[c]([OX2H])[c][c][c]1",
                   "very oxidation-prone; consider a glycoside or ester derivative"),
    StabilityAlert("enediol", "oxidation", "[OX2H][CX3]=[CX3][OX2H]",
                   "ascorbic acid class; anhydrous base or a stabilised derivative"),
    StabilityAlert("thiol", "oxidation", "[SX2H]",
                   "oxidises to disulfide and smells; anhydrous or encapsulated"),
    StabilityAlert("aldehyde", "oxidation", "[CX3H1](=[OX1])[#6]",
                   "oxidises to the acid; also a sensitisation alert"),

    StabilityAlert("ester", "hydrolysis", "[CX3;!R](=[OX1])[OX2][#6;!$([CX3]=[OX1])]",
                   "control pH; may be intentional if this is a prodrug"),
    StabilityAlert("lactone", "hydrolysis", "[#6;R][CX3;R](=[OX1])[OX2;R][#6;!$([CX3]=[OX1])]",
                   "ring-opens at alkaline pH; hold the formulation acidic"),
    StabilityAlert("glycoside", "hydrolysis",
                   "[OX2;R][CX4;R]([OX2][#6])",
                   "acid- and enzyme-labile; can be the intended slow-release form"),
    StabilityAlert("anhydride", "hydrolysis", "[CX3](=[OX1])[OX2][CX3](=[OX1])",
                   "hydrolyses on contact with water; not viable in an emulsion"),
)

_COMPILED = [(a, Chem.MolFromSmarts(a.smarts)) for a in _ALERTS]
_COMPILED = [(a, p) for a, p in _COMPILED if p is not None]


# The -N=N- SMARTS is taken from the photostability alert above rather than
# written out a second time: the two uses must agree about what an azo group
# is, and a second copy is a second thing to forget to update.
_AZO = next(p for a, p in _COMPILED if a.name == "azo")

# Thresholds for `colour_risk`, all read off the largest single conjugated
# system. Each is anchored to a pair of real compounds that straddle it; the
# justification is in the `colour_risk` docstring.
_VISIBLE_FUSED_RINGS = 3        # anthracene 3 / naphthalene 2
_VISIBLE_FUSED_PI_ATOMS = 16    # tetracene 18 orange / anthracene 14 white
_VISIBLE_POLYENE = 7            # beta-carotene 9 orange / retinol 5 pale
_TINT_POLYENE = 5
_VISIBLE_CHAIN_PI_BONDS = 3     # curcumin 4 yellow / avobenzone 2 colourless
_VISIBLE_CHAIN_PI_ATOMS = 20    # curcumin 25 / avobenzone 18
_TINT_FUSED_RINGS = 2           # naphthalene pair, 275 nm on its own
_TINT_FUSED_PI_ATOMS = 20       # quercetin 22 yellow / chlorpromazine 13 white


def stability_alerts(smiles: str) -> list[dict]:
    mol = assessable_molecule(smiles)
    if mol is None:
        return []
    out = []
    for alert, patt in _COMPILED:
        if mol.GetSubstructMatches(patt):
            out.append({"name": alert.name, "mode": alert.mode,
                        "mitigation": alert.mitigation})
    return out


def colour_risk(smiles: str) -> str:
    'Rough call on whether the molecule will tint the product.'
    mol = assessable_molecule(smiles)
    if mol is None:
        return "unknown"
    # Descriptors of the largest single conjugated system, computed once in
    # admet_photo so the colour call and the photosafety call cannot drift.
    f = photo.features(mol)
    n_pi = f["n_pi_atoms"]
    n_fused = f["n_fused_pi_rings"]
    n_ext = f["extending_pi_bonds"]
    n_ene = f["polyene_length"]
    is_azo = bool(mol.GetSubstructMatches(_AZO))

    if (is_azo
            or (n_fused >= _VISIBLE_FUSED_RINGS and n_pi >= _VISIBLE_FUSED_PI_ATOMS)
            or n_ene >= _VISIBLE_POLYENE
            or (n_ext >= _VISIBLE_CHAIN_PI_BONDS and n_pi >= _VISIBLE_CHAIN_PI_ATOMS)):
        return "strong (expect visible colour)"
    if (n_fused >= _VISIBLE_FUSED_RINGS
            or n_ene >= _TINT_POLYENE
            or n_ext >= _VISIBLE_CHAIN_PI_BONDS
            or (n_fused >= _TINT_FUSED_RINGS and n_pi >= _TINT_FUSED_PI_ATOMS)):
        return "moderate (may tint)"
    return "low"


def summarise(smiles: str) -> dict:
    if assessable_molecule(smiles) is None:
        return {"alerts": [], "modes": [], "colour": "unknown",
                "assessment_status": "NOT_EVALUATED",
                "reason": "invalid_or_empty_structure",
                "verdict": "unknown stability risk; structure not assessable"}
    alerts = stability_alerts(smiles)
    modes = sorted({a["mode"] for a in alerts})
    colour = colour_risk(smiles)
    if not alerts and colour == "low":
        verdict = "no stability flag"
    else:
        parts = []
        if "photo" in modes:
            parts.append("light-sensitive")
        if "oxidation" in modes:
            parts.append("oxidises in air")
        if "hydrolysis" in modes:
            parts.append("hydrolysis-prone")
        if colour != "low":
            parts.append(f"colour {colour}")
        verdict = "; ".join(parts) or "minor flags only"
    return {"alerts": alerts, "modes": modes, "colour": colour, "verdict": verdict,
            "assessment_status": "EVALUATED", "reason": None}
