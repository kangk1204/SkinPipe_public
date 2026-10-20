"""A curated bioisostere table, applied as substructure transformations.

Where the idea comes from
-------------------------
Langmuir (1919, J Am Chem Soc 41:1543) noticed that molecules with the same
number of atoms and the same valence-electron count behave alike; Grimm's
hydride displacement law (1925) generalised it, and Friedman (1951) coined
"bioisostere" for the case where the resemblance survives contact with a
protein. Patani & LaVoie's review (1996, Chem Rev 96:3147) is still the
standard catalogue of classical replacements, and Meanwell (2011, J Med Chem
54:2529) covers the non-classical ones that dominate modern practice.

Why bioisosteres and not just "generate analogues"
--------------------------------------------------
A bioisostere is the one kind of structural change that comes with a prior
reason to expect the binding interaction to survive: the replacement group
presents the same hydrogen bonds, the same charge at physiological pH, or the
same aromatic quadrupole from roughly the same geometry. Carboxylic acid ->
tetrazole is the worked example everyone quotes. Tetrazole is about as acidic
(pKa ~4.9 vs ~4.2), so it is ionised in the same compartments, and it is a
planar, delocalised anion of similar footprint -- but it is ~10x more
lipophilic than the carboxylate, and it is a completely different ring system
to a patent examiner. That is precisely the combination this module is hunting
for: same interaction, different chemistry (Ballatore, Huryn & Smith 2013,
ChemMedChem 8:385; Herr 2002, Bioorg Med Chem 10:3379).

What this table can and cannot do
---------------------------------
CAN: apply a documented, literature-sourced replacement wherever its SMARTS
matches, and tell you which rule fired, so every analogue is traceable to a
citation rather than to a black box.

CANNOT: know whether the replacement is right *here*. A bioisostere is a
statement about average behaviour across many series, and every one of these
rules has published failures. The catechol -> chelator rules below are the
sharpest example: on a tyrosinase inhibitor the catechol is often the group
that coordinates the dicopper centre, so swapping it is either the whole point
(a novel chelator with the same metal-binding geometry) or fatal (a group that
no longer reaches the metal). Nothing in this file can tell those apart. The
3D scoring in `hopping.py` is a weak proxy; an assay is the answer.

Also: the RDKit reaction machinery can drop stereochemistry on double bonds
immediately adjacent to a transformed atom. Products are re-standardised on
the way out, but a lost E/Z assignment is a real possibility on rules that
rebuild a ring the alkene hangs off. Check any hit you intend to synthesise.

Skin-specific note
------------------
Several replacements move a molecule the wrong way for topical delivery. Acid
-> tetrazole adds ~40 Da and keeps the molecule anionic, and ionised species
barely cross the stratum corneum; acid -> hydroxamic acid adds hydrogen-bond
donors, which `filters.topical` penalises twice. Rank hops with the topical
filter attached (`rank_hops` does this) rather than on hop quality alone.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Iterable, Optional

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

from .library import standardize

RDLogger.DisableLog("rdApp.*")
log = logging.getLogger(__name__)


@dataclass(frozen=True)
class Bioisostere:
    """One documented replacement, as a reaction SMARTS.

    `klass` separates the classical replacements (same valence-electron count,
    Grimm/Erlenmeyer lineage) from the non-classical ones (no atom-count
    correspondence at all; they work because the *interaction* is conserved,
    not the atoms). Non-classical replacements are where the interesting hops
    live, and also where the failure rate is highest.
    """
    name: str
    group_from: str
    group_to: str
    klass: str            # classical | non-classical | ring | chelator
    smarts: str           # reaction SMARTS, reactant>>product
    citation: str
    note: str = ""


# Ordered roughly by how often each is reached for in practice. Every SMARTS
# here is exercised by tests/test_hopping.py against a molecule it must hit.
_BIOISOSTERES: tuple[Bioisostere, ...] = (
    # --- carboxylic acid surrogates ---------------------------------------
    # The single most-worked problem in the field, because acids are potent
    # and poorly absorbed. Ballatore et al. 2013 tabulate ~20; these four
    # cover the acidity range that matters.
    Bioisostere(
        "acid_to_tetrazole", "-COOH", "1H-tetrazol-5-yl", "non-classical",
        "[#6:1][CX3](=[OX1])[OX2H1]>>[#6:1]c1nn[nH]n1",
        "Herr 2002 Bioorg Med Chem 10:3379; Ballatore et al. 2013 ChemMedChem 8:385",
        "pKa ~4.9 vs ~4.2, so still anionic at skin pH; ~10x more lipophilic "
        "than the carboxylate but ~40 Da heavier and metabolically robust"),
    Bioisostere(
        "acid_to_acylsulfonamide", "-COOH", "-C(=O)NHSO2Me", "non-classical",
        "[#6:1][CX3](=[OX1])[OX2H1]>>[#6:1]C(=O)NS(C)(=O)=O",
        "Ballatore et al. 2013 ChemMedChem 8:385",
        "acidity tunable by the sulfonyl substituent; larger footprint, so it "
        "needs room in the pocket that the carboxylate did not"),
    Bioisostere(
        "acid_to_hydroxamic", "-COOH", "-C(=O)NHOH", "non-classical",
        "[#6:1][CX3](=[OX1])[OX2H1]>>[#6:1]C(=O)NO",
        "Patani & LaVoie 1996 Chem Rev 96:3147",
        "weaker acid (pKa ~9) and a strong bidentate metal chelator -- a "
        "liability at an off-target metalloenzyme, an asset at a real one"),
    Bioisostere(
        "acid_to_hydroxyisoxazole", "-COOH", "3-hydroxyisoxazol-5-yl", "non-classical",
        "[#6:1][CX3](=[OX1])[OX2H1]>>[#6:1]c1cc(O)no1",
        "Ballatore et al. 2013 ChemMedChem 8:385",
        "planar acidic heterocycle; the classic route out of a flat "
        "carboxylate chemotype into unclaimed ring space"),

    # --- amide surrogates --------------------------------------------------
    Bioisostere(
        "amide_to_1_2_4_oxadiazole", "-C(=O)NH-", "1,2,4-oxadiazole", "ring",
        "[#6:1][CX3](=[OX1])[NX3H1][#6:3]>>[#6:1]c1noc([#6:3])n1",
        "Bostrom et al. 2012 J Med Chem 55:1817",
        "keeps the amide dipole and the two acceptors, loses the N-H donor "
        "and the amide's rotational freedom -- rigidification that helps if "
        "the bound conformer is the planar one and hurts otherwise"),
    Bioisostere(
        "amide_to_1_3_4_oxadiazole", "-C(=O)NH-", "1,3,4-oxadiazole", "ring",
        "[#6:1][CX3](=[OX1])[NX3H1][#6:3]>>[#6:1]c1nnc([#6:3])o1",
        "Bostrom et al. 2012 J Med Chem 55:1817",
        "different acceptor vector to the 1,2,4-isomer; worth enumerating both "
        "because they are not interchangeable in practice"),
    Bioisostere(
        "amide_to_triazole", "-C(=O)NH-", "1,2,3-triazol-1,4-diyl", "ring",
        "[#6:1][CX3](=[OX1])[NX3H1][#6:3]>>[#6:1]c1cn([#6:3])nn1",
        "Meanwell 2011 J Med Chem 54:2529",
        "amide mimic with a large dipole and no hydrolysable bond; "
        "synthetically cheap via Huisgen cycloaddition"),
    Bioisostere(
        "amide_to_reverse_amide", "-C(=O)NH-", "-NHC(=O)-", "classical",
        "[#6:1][CX3](=[OX1])[NX3H1][#6:3]>>[#6:1]NC(=O)[#6:3]",
        "Patani & LaVoie 1996 Chem Rev 96:3147",
        "same atoms, donor and acceptor swapped end for end; a minimal change "
        "that nonetheless often walks out of a narrow composition-of-matter claim"),

    # --- ester surrogates --------------------------------------------------
    Bioisostere(
        "ester_to_oxadiazole", "-C(=O)O-", "1,2,4-oxadiazole", "ring",
        "[#6:1][CX3](=[OX1])[OX2][#6;!$([#6]=O):3]>>[#6:1]c1noc([#6:3])n1",
        "Bostrom et al. 2012 J Med Chem 55:1817",
        "the standard fix for an ester that skin esterases cleave before it "
        "reaches the viable epidermis -- relevant for any pro-drug-shaped active"),
    Bioisostere(
        "ester_to_amide", "-C(=O)O-", "-C(=O)NH-", "classical",
        "[#6:1][CX3](=[OX1])[OX2][#6;!$([#6]=O):3]>>[#6:1]C(=O)N[#6:3]",
        "Patani & LaVoie 1996 Chem Rev 96:3147",
        "O -> NH: hydrolytically far more stable, one extra donor, and the "
        "amide prefers a different torsion to the ester"),

    # --- aryl ring replacements -------------------------------------------
    Bioisostere(
        "phenyl_to_pyridine", "phenyl", "pyridyl", "ring",
        "[cH1:1]1[cH1:2][cH1:3][cH1:4][cH1:5][c:6]1[*:7]>>[cH1:1]1[cH1:2][n:3][cH1:4][cH1:5][c:6]1[*:7]",
        "Patani & LaVoie 1996 Chem Rev 96:3147",
        "one CH -> N: ~1 log unit off logP, one acceptor gained, aromatic "
        "footprint unchanged. The cheapest real hop available on most scaffolds"),
    Bioisostere(
        "phenyl_to_thiophene", "phenyl", "thiophen-2-yl", "ring",
        "[cH1:1]1[cH1:2][cH1:3][cH1:4][cH1:5][c:6]1[*:7]>>[cH1:1]1[cH1:2][cH1:3][c:6](s1)[*:7]",
        "Patani & LaVoie 1996 Chem Rev 96:3147",
        "six-ring to five-ring: substituent vectors move by ~30 degrees, which "
        "is either a productive perturbation or a lost contact"),
    Bioisostere(
        "phenyl_to_bicyclopentane", "1,4-phenylene", "bicyclo[1.1.1]pentane", "non-classical",
        "[*;!#1:1][c:2]1[cH][cH][c:3]([*;!#1:4])[cH][cH]1>>[*:1]C12CC(C1)(C2)[*:4]",
        "Stepan et al. 2012 J Med Chem 55:3414",
        "para-phenylene -> BCP: same ~1.9 A exit-vector separation and same "
        "linear geometry, but sp3, so logP and metabolic liability drop and "
        "aqueous solubility rises. Only valid for 1,4-disubstituted rings, and "
        "it deletes the aromatic ring, so any pi-stacking contact is gone"),

    # --- metal-chelating head groups --------------------------------------
    # Directly relevant to this project: tyrosinase is a type-3 dicopper
    # enzyme, and most of its useful inhibitors are copper chelators wearing
    # different rings. Swapping one chelator for another is the classic hop in
    # this target class -- kojic acid's own 3-hydroxy-4-pyranone and tropolone
    # are exactly that relationship to a catechol.
    Bioisostere(
        "catechol_to_hydroxypyridinone", "catechol", "3-hydroxy-4-pyridinone", "chelator",
        "[*;!#1:1][c:2]1[cH:3][cH:4][c:5]([OX2H])[c:6]([OX2H])[cH:7]1"
        ">>[*:1]C1=CC(=O)C(O)=CN1C",
        "Patani & LaVoie 1996 Chem Rev 96:3147 (chelator isosterism); "
        "deferiprone as the reference 3-hydroxy-4-pyridinone",
        "keeps a bidentate O,O donor set with similar bite angle; removes the "
        "catechol's air-oxidation-to-quinone liability that filters.safety "
        "flags as a prohapten"),
    Bioisostere(
        "catechol_to_tropolone", "catechol", "tropolone", "chelator",
        "[*;!#1:1][c:2]1[cH:3][cH:4][c:5]([OX2H])[c:6]([OX2H])[cH:7]1"
        ">>[*:1]C1=CC=CC(=O)C(O)=C1",
        "Patani & LaVoie 1996 Chem Rev 96:3147 (chelator isosterism); "
        "tropolone is a long-known low-micromolar tyrosinase inhibitor",
        "seven-membered aromatic alpha-hydroxyketone; a genuinely different "
        "ring system that presents the same chelating oxygens"),
    Bioisostere(
        "catechol_to_pyranone", "catechol", "5-hydroxy-4H-pyran-4-one", "chelator",
        "[*;!#1:1][c:2]1[cH:3][cH:4][c:5]([OX2H])[c:6]([OX2H])[cH:7]1"
        ">>[*:1]C1=CC(=O)C(O)=CO1",
        "kojic acid as the reference 3-hydroxy-4-pyranone chelator",
        "walks a catechol series onto the kojic acid head group -- useful as a "
        "positive control, useless for novelty, since kojic acid is prior art"),
    # The reverse direction, off the kojic/maltol head group itself. These are
    # the rules that give a 4-pyranone seed anywhere to go: the ring oxygen is
    # the only heteroatom, so swapping it is the one edit that changes the
    # scaffold while leaving the chelating O,O pair in place.
    Bioisostere(
        "pyranone_to_pyridinone", "4H-pyran-4-one", "1-methyl-4-pyridinone", "chelator",
        "[o;r6:1]1[c:2][c:3][c:4](=[O:9])[c:5][c:6]1>>[n:1](C)1[c:2][c:3][c:4](=[O:9])[c:5][c:6]1",
        "Patani & LaVoie 1996 Chem Rev 96:3147 (ring O/NR isosterism); "
        "maltol -> deferiprone is this exact transformation",
        "ring O -> N-Me keeps the O,O bite angle and raises the pKa of the "
        "chelating hydroxyl, which strengthens metal binding at skin pH. "
        "Synthetically real: condense the pyranone with methylamine"),
    Bioisostere(
        "pyranone_to_thiopyranone", "4H-pyran-4-one", "4H-thiopyran-4-one", "classical",
        "[o;r6:1]1[c:2][c:3][c:4](=[O:9])[c:5][c:6]1>>[s:1]1[c:2][c:3][c:4](=[O:9])[c:5][c:6]1",
        "Grimm 1925 Z Elektrochem 31:474 (divalent O/S displacement)",
        "the cheapest scaffold change available on a pyranone, but sulfur "
        "heterocycles of this type are awkward to make -- expect the SA score "
        "to reject most of what this rule produces"),
    Bioisostere(
        "pyranone_to_thionopyridinone", "3-hydroxy-4-pyranone",
        "3-hydroxy-4-pyridinethione", "chelator",
        "[o;r6:1]1[c:2][c:3][c:4](=[O:9])[c:5][c:6]1"
        ">>[n:1](C)1[c:2][c:3][c:4](=[S:9])[c:5][c:6]1",
        "Patani & LaVoie 1996 Chem Rev 96:3147 (C=O/C=S isosterism)",
        "an O,S donor set instead of O,O. Softer, so it prefers Cu(I) over "
        "Cu(II) -- a real change in which oxidation state of the tyrosinase "
        "copper it would engage, not a cosmetic edit"),
    Bioisostere(
        "hydroxypyranone_to_tropolone", "5-hydroxy-4H-pyran-4-one", "tropolone", "chelator",
        "[*;!#1:7][c:6]1[o:1][cH:2][c:3]([OX2H])[c:4](=[O:9])[cH:5]1"
        ">>[*:7]C1=CC=CC(=O)C(O)=C1",
        "tropolone as the reference seven-membered alpha-hydroxyketone "
        "tyrosinase inhibitor",
        "positional, not general: written for the kojic acid substitution "
        "pattern (2-substituted-5-hydroxy-4-pyranone) and it will not match "
        "maltol's 3-hydroxy-2-substituted isomer. Six-ring to seven-ring is "
        "about as large a scaffold change as this table makes"),

    # --- classical single-atom swaps --------------------------------------
    # Grimm's hydride displacement law in its original form. These barely move
    # 2D similarity, so they rarely score as hops; they are here because they
    # are cheap to make and sometimes fix a metabolic or stability problem.
    Bioisostere(
        "phenol_to_fluorine", "Ar-OH", "Ar-F", "classical",
        "[c:1][OX2H]>>[c:1]F",
        "Langmuir 1919 J Am Chem Soc 41:1543; Grimm 1925 Z Elektrochem 31:474",
        "isosteric by electron count and nearly isosteric by size, but it "
        "deletes a hydrogen-bond donor entirely. Fatal whenever that OH was "
        "the binding interaction -- which on a tyrosinase inhibitor it usually is"),
    Bioisostere(
        "ether_to_thioether", "-O-", "-S-", "classical",
        "[#6:1][OX2;!$([OX2][CX3]=O)][CH3:2]>>[#6:1]S[CH3:2]",
        "Grimm 1925 Z Elektrochem 31:474; Patani & LaVoie 1996 Chem Rev 96:3147",
        "divalent-for-divalent; S is larger and far less basic, and it brings "
        "an oxidative-metabolism soft spot the ether did not have"),
    Bioisostere(
        "methyl_to_trifluoromethyl", "-CH3", "-CF3", "classical",
        "[c:1][CH3]>>[c:1]C(F)(F)F",
        "Meanwell 2011 J Med Chem 54:2529",
        "similar steric bulk to an isopropyl, opposite electronics to methyl; "
        "blocks benzylic oxidation and raises logP by ~0.9"),
)


def table() -> tuple[Bioisostere, ...]:
    """The curated replacement set, for inspection or documentation."""
    return _BIOISOSTERES


def _compile() -> list[tuple[Bioisostere, AllChem.ChemicalReaction]]:
    out: list[tuple[Bioisostere, AllChem.ChemicalReaction]] = []
    for b in _BIOISOSTERES:
        try:
            rxn = AllChem.ReactionFromSmarts(b.smarts)
        except Exception as exc:                          # pragma: no cover
            log.error("bad reaction SMARTS for %s: %s (%s)", b.name, b.smarts, exc)
            continue
        if rxn is None:                                   # pragma: no cover
            log.error("bad reaction SMARTS for %s: %s", b.name, b.smarts)
            continue
        rxn.Initialize()
        out.append((b, rxn))
    return out


_COMPILED = _compile()


def applicable(smiles: str) -> list[str]:
    """Names of the rules whose reactant pattern matches this molecule.

    Cheap enough to call before enumerating, e.g. to explain why a seed
    produced nothing.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []
    return [b.name for b, rxn in _COMPILED if rxn.IsMoleculeReactant(mol)]


def bioisostere_replace(smiles: str, rules: Optional[Iterable[str]] = None,
                        max_per_rule: int = 4) -> list[dict]:
    """Apply every matching replacement, one at a time, and report provenance.

    Single-edit only, deliberately. Stacking replacements compounds the risk
    that the binding interaction is gone, and it makes the "which change did
    this?" question unanswerable. Call it again on its own output if you want
    double-edited analogues, and accept that the hit rate drops.

    `max_per_rule` caps products from one rule on one molecule, which matters
    for symmetric aromatics: phenyl -> pyridine on a plain phenyl ring matches
    at several equivalent positions and RDKit returns each match separately.
    Duplicates are removed by canonical SMILES first, so the cap applies to
    genuinely distinct products.

    Returns one dict per analogue with keys `smiles`, `rule`, `group_from`,
    `group_to`, `klass`, `citation`, `note`. The parent is never returned:
    several rules (reverse amide on a symmetric substrate, for instance) can
    regenerate the input, and that is not an analogue.

    >>> hits = bioisostere_replace("O=C(O)c1ccccc1")   # benzoic acid
    >>> any(h["rule"] == "acid_to_tetrazole" for h in hits)
    True
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        log.warning("bioisostere_replace: unparseable SMILES %r", smiles[:80])
        return []
    parent = Chem.MolToSmiles(mol)
    wanted = set(rules) if rules is not None else None

    out: list[dict] = []
    seen: set[str] = {parent}
    for b, rxn in _COMPILED:
        if wanted is not None and b.name not in wanted:
            continue
        try:
            products = rxn.RunReactants((mol,))
        except Exception as exc:                          # pragma: no cover
            log.debug("%s failed on %s: %s", b.name, parent[:60], exc)
            continue
        kept = 0
        for prod_tuple in products:
            if kept >= max_per_rule:
                break
            product = prod_tuple[0]
            try:
                Chem.SanitizeMol(product)
            except Exception:
                continue                       # rule produced an invalid valence
            # Reuse the library's desalt/neutralise step so generated molecules
            # enter the pipeline in exactly the same form as fetched ones.
            std = standardize(Chem.MolToSmiles(product))
            if std is None or std in seen:
                continue
            seen.add(std)
            kept += 1
            out.append({"smiles": std, "rule": b.name, "group_from": b.group_from,
                        "group_to": b.group_to, "klass": b.klass,
                        "citation": b.citation, "note": b.note})
    log.debug("bioisostere_replace(%s): %d analogues from %d rules",
              parent[:40], len(out), len({h["rule"] for h in out}))
    return out
