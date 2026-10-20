"""Structure-based proxy for the OECD TG 432 / SCCS light-absorption screen.

The rule being approximated
---------------------------
ICH S10 (2013) and the SCCS Notes of Guidance (SCCS/1647/22, section 3-3.9)
share one screening principle, first argued quantitatively by Henry et al.
(2009) Regul Toxicol Pharmacol 55:e1:

    A chemical that does not absorb light between 290 and 700 nm with a molar
    extinction coefficient of at least 1000 L mol-1 cm-1 is not expected to be
    directly phototoxic, and needs no further photosafety testing.

290 nm is where terrestrial sunlight starts (the ozone cut-off); 700 nm is
where photochemistry stops being able to populate a reactive triplet. The
threshold is deliberately conservative - measured phototoxic drugs sit two to
three orders of magnitude above it.

Why this module is a proxy, and what that costs
-----------------------------------------------
A molar extinction coefficient is an experimental quantity. It cannot be
derived from a connection table without a quantum-chemical excited-state
calculation, and even TD-DFT gets oscillator strengths wrong by factors of two.
What structure does predict reliably is the *qualitative* presence of a
low-energy, symmetry-allowed pi->pi* transition, which is what a large epsilon
above 290 nm requires. So this module scores:

*   the extent of the largest conjugated pi system (rings and chain),
*   whether a carbonyl extends that system (and is not an amide, whose
    nitrogen lone pair blocks the extension - this is what separates
    anthraquinone from caffeine; an ester is NOT excluded, see
    `extending_pi_bonds`),
*   donor/acceptor and donor-only substitution, which red-shifts an otherwise
    sub-290 nm benzene band into sunlight (engine in `admet_chromophore`),
*   membership of a chemotype whose phototoxicity is already documented.

The output is an estimated *band* for epsilon, never a number, and the verdict
is always "run TG 432", never "this compound is phototoxic".

Known limitation, stated rather than tuned away: this proxy cannot distinguish
a strong absorber that is photostable and harmless from one that photoreacts.
It over-calls polyphenols, which do absorb past 290 nm and do generate singlet
oxygen, and which are nonetheless used at scale in cosmetics. That is the
correct behaviour for a screen whose job is to nominate compounds for TG 432.

The commercial case that makes this worth doing
------------------------------------------------
Bergapten (5-methoxypsoralen) in expressed bergamot and other citrus oils is
capped at 1 ppm in leave-on EU cosmetics (Annex II/III) because furocoumarins
photobind to DNA under UVA and cause berloque dermatitis. A pipeline that ranks
a citrus-derived active highly and never notices its furocoumarin content has
produced a commercially unusable recommendation.
"""
from __future__ import annotations

import logging

from rdkit import Chem, RDLogger

from . import admet_chromophore as chromo

RDLogger.DisableLog("rdApp.*")
log = logging.getLogger(__name__)

# Ring skeletons are written with element-only atoms and "any" bonds. Bond
# orders vary because RDKit aromatises the coumarin pyranone but not the
# anthraquinone central ring; ring-size primitives are avoided because RDKit's
# `r<n>` tests the MINIMUM ring size, which is 5 for every atom psoralen shares
# between its furan and its benzo ring. See admet_alerts for the full note.
_COUMARIN_CORE = Chem.MolFromSmarts(
    "[#8]1~[#6](=[OX1])~[#6]~[#6]~[#6]2~[#6]~[#6]~[#6]~[#6]~[#6]12")
_FURAN_RING = Chem.MolFromSmarts("[#8]1~[#6]~[#6]~[#6]~[#6]1")
_PYRROLE_RING = Chem.MolFromSmarts("[#7]1~[#6]~[#6]~[#6]~[#6]1")

# A 1,3-dicarbonyl enolises, and the enol is a single conjugated system across
# the sp3 bridge. Curcumin and avobenzone are both drawn in the keto form and
# are both fully conjugated in reality - curcumin is yellow (425 nm) and
# avobenzone is a licensed UVA filter (357 nm). Treating the bridge as sp3
# would score both as having no chromophore at all.
_ENOLISABLE_DIKETONE = Chem.MolFromSmarts(
    "[#6X3](=[OX1])[#6X4;H1,H2][#6X3]=[OX1]")

# Donor/acceptor perception and the biaryl charge-transfer test live in
# `admet_chromophore`; see that module for the aromaticity bug that made every
# licensed non-fused UV filter read as having no chromophore.

# A dialkyl ketone carbonyl, reported as a descriptor only - it is not scored.
# `extending_pi_bonds` below is what decides whether a carbonyl extends the
# chromophore, and it applies a deliberately narrower rule than this pattern.
_EXTENDING_CO = Chem.MolFromSmarts("[CX3](=[OX1])([#6])[#6]")
_CONJ_CO_TO_ARYL = Chem.MolFromSmarts("[c][CX3](=[OX1])[#6]")


def ring_systems(mol: Chem.Mol) -> list[set[int]]:
    """Group SSSR rings into fused systems (rings sharing at least one bond)."""
    systems: list[set[int]] = []
    for ring in mol.GetRingInfo().AtomRings():
        atoms = set(ring)
        merged = [s for s in systems if len(s & atoms) >= 2]
        for s in merged:
            systems.remove(s)
            atoms |= s
        systems.append(atoms)
    return systems


def max_fused_rings(mol: Chem.Mol) -> tuple[int, int]:
    """(rings, aromatic rings) in the largest fused ring system."""
    rings = mol.GetRingInfo().AtomRings()
    best = (0, 0)
    for system in ring_systems(mol):
        n = sum(1 for r in rings if set(r) <= system)
        n_arom = sum(1 for r in rings
                     if set(r) <= system
                     and all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in r))
        if (n, n_arom) > best:
            best = (n, n_arom)
    return best


def conjugated_components(mol: Chem.Mol) -> list[set[int]]:
    """Connected components of the subgraph of conjugated/aromatic bonds."""
    parent = {a.GetIdx(): a.GetIdx() for a in mol.GetAtoms()}

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for bond in mol.GetBonds():
        if not (bond.GetIsConjugated() or bond.GetIsAromatic()):
            continue
        a, b = find(bond.GetBeginAtomIdx()), find(bond.GetEndAtomIdx())
        if a != b:
            parent[a] = b

    for match in mol.GetSubstructMatches(_ENOLISABLE_DIKETONE):
        head = find(match[0])
        for idx in match[1:]:
            other = find(idx)
            if other != head:
                parent[other] = head
                head = find(idx)

    groups: dict[int, set[int]] = {}
    for idx in parent:
        groups.setdefault(find(idx), set()).add(idx)
    return [g for g in groups.values() if len(g) > 1]


def longest_polyene(mol: Chem.Mol) -> int:
    """Longest run of conjugated, non-aromatic C=C bonds.

    Retinoic acid has five (four in the chain plus the cyclohexene), which is
    what puts its band at 350 nm; beta-carotene has eleven, which is what makes
    it orange. Aromatic bonds are excluded so a benzene ring is not counted as
    three conjugated alkenes.
    """
    pi = [b for b in mol.GetBonds()
          if b.GetBondType() == Chem.BondType.DOUBLE
          and not b.GetIsAromatic()
          and b.GetBeginAtom().GetAtomicNum() == 6
          and b.GetEndAtom().GetAtomicNum() == 6]
    if not pi:
        return 0
    idx = {b.GetIdx(): i for i, b in enumerate(pi)}
    adj: dict[int, set[int]] = {i: set() for i in range(len(pi))}
    for bond in mol.GetBonds():                 # single bonds that link two pi bonds
        if bond.GetIdx() in idx or not bond.GetIsConjugated():
            continue
        ends = (bond.GetBeginAtom(), bond.GetEndAtom())
        touching = [idx[b.GetIdx()] for atom in ends for b in atom.GetBonds()
                    if b.GetIdx() in idx]
        for i in touching:
            for j in touching:
                if i != j:
                    adj[i].add(j)
                    adj[j].add(i)

    best = 0

    def walk(node: int, seen: frozenset[int]) -> int:
        nonlocal best
        best = max(best, len(seen))
        for nxt in adj[node]:
            if nxt not in seen:
                walk(nxt, seen | {nxt})
        return best

    for start in range(len(pi)):
        walk(start, frozenset({start}))
    return best


def is_furocoumarin(mol: Chem.Mol) -> bool:
    """Coumarin core with a furan ring fused to it (psoralen, bergapten).

    Detected as a ring intersection rather than one SMARTS, because the linear
    (psoralen) and angular (angelicin) isomers share no single tractable
    pattern and a Kekule-form SMARTS would miss whichever form RDKit does not
    perceive.
    """
    coumarins = mol.GetSubstructMatches(_COUMARIN_CORE)
    furans = mol.GetSubstructMatches(_FURAN_RING)
    return any(len(set(c) & set(f)) >= 2 for c in coumarins for f in furans)


def is_tetrapyrrole(mol: Chem.Mol) -> bool:
    """Four pyrrole-type rings inside one conjugated system (porphyrins)."""
    pyrroles = mol.GetSubstructMatches(_PYRROLE_RING)
    if len(pyrroles) < 4:
        return False
    return any(sum(1 for p in pyrroles if set(p) <= comp) >= 4
               for comp in conjugated_components(mol))


def extending_pi_bonds(mol: Chem.Mol, component: set[int]) -> int:
    """Non-aromatic pi bonds inside `component` that genuinely extend it.

    A carbonyl whose carbon also carries a nitrogen is cross-conjugated: the
    amide nitrogen lone pair feeds the same carbon and the band stays in the
    deep UV. That single exclusion is what separates caffeine (two lactam
    carbonyls on a bicyclic aromatic, lambda_max 273 nm, harmless in sunlight)
    from anthraquinone (two aryl ketones, absorbs to 430 nm, phototoxic).

    Nitrogen is the only exclusion, deliberately: an ester alkoxy lone pair
    cross-conjugates far more weakly than an amide nitrogen's - the reason
    amides are planar with a large rotational barrier and esters are not - so
    an ester carbonyl still extends an aryl chromophore. Styrene absorbs at
    248 nm, methyl cinnamate at 280 nm and 2-ethylhexyl 4-methoxycinnamate at
    310 nm, where it is a licensed UV filter. The cost is that methyl benzoate
    is over-called at one extending bond, which for a screen that nominates
    compounds for TG 432 is the right direction to err.
    """
    n = 0
    for bond in mol.GetBonds():
        if bond.GetIsAromatic():
            continue
        if bond.GetBondType() not in (Chem.BondType.DOUBLE, Chem.BondType.TRIPLE):
            continue
        a, b = bond.GetBeginAtom(), bond.GetEndAtom()
        if a.GetIdx() not in component or b.GetIdx() not in component:
            continue
        if {a.GetAtomicNum(), b.GetAtomicNum()} == {6, 8}:
            carbon = a if a.GetAtomicNum() == 6 else b
            if any(nb.GetAtomicNum() == 7 for nb in carbon.GetNeighbors()):
                continue
        n += 1
    return n


def features(mol: Chem.Mol) -> dict:
    """Raw chromophore descriptors feeding the extinction-coefficient proxy."""
    components = conjugated_components(mol)
    largest = max(components, key=len) if components else set()
    rings = mol.GetRingInfo().AtomRings()
    n_rings_conj = sum(1 for r in rings if set(r) <= largest)
    n_fused, n_fused_arom = max_fused_rings(mol)
    # Only rings that are themselves part of the conjugated system count towards
    # the chromophore. Without this restriction a triterpene or a steroid - five
    # fused but entirely saturated rings - scores like an acene, and ursolic
    # acid comes out looking as phototoxic as anthracene.
    n_fused_pi = max((sum(1 for r in rings if set(r) <= system and set(r) <= largest)
                      for system in ring_systems(mol)), default=0)
    donors, acceptors = chromo.donor_acceptor_atoms(mol)
    has_donor = bool(donors & largest)
    push_pull = has_donor and bool(acceptors & largest)
    has_aromatic = any(mol.GetAtomWithIdx(i).GetIsAromatic() for i in largest)
    return {
        "n_pi_atoms": len(largest),
        "n_fused_pi_rings": n_fused_pi,
        "aromatic_pi_system": has_aromatic,
        "extending_pi_bonds": extending_pi_bonds(mol, largest) if has_aromatic else 0,
        "n_rings_in_pi_system": n_rings_conj,
        "n_fused_rings": n_fused,
        "n_fused_aromatic_rings": n_fused_arom,
        "polyene_length": longest_polyene(mol),
        "aryl_ketone": bool(mol.GetSubstructMatches(_CONJ_CO_TO_ARYL)),
        "extending_carbonyl": bool(mol.GetSubstructMatches(_EXTENDING_CO)),
        "push_pull": push_pull,
        # A donor alone on a two-ring-or-larger system is the classic
        # red-shift: naphthalene 275 nm, 2-methoxynaphthalene 326 nm.
        "donor_on_extended_system": has_donor and n_rings_conj >= 2,
        "charge_transfer_biaryl": chromo.charge_transfer_biaryl(
            mol, largest, push_pull),
        "furocoumarin": is_furocoumarin(mol),
        "tetrapyrrole": is_tetrapyrrole(mol),
    }


# Weights for the extinction proxy.
#
# What this actually scores, measured rather than claimed: the proxy rewards
# fused-ring extent and polyene length, so the TG 432 positive controls it
# clears 3.0 on are the condensed polycyclics and the polyenes. It does not
# carry the phenothiazines -- chlorpromazine and promethazine both score 1.30,
# and caffeine, which has no sunlight-window absorption worth the name, scores
# 0.80 rather than 0. Those are not misses: `PHOTO_CHEMOTYPES` recognises the
# phenothiazine, furocoumarin and tetrapyrrole classes by structure and is what
# returns `high` for them, and the proxy is the continuous fallback for
# everything the table does not name.
#
# The comment here used to say the weights were "calibrated so that the OECD TG
# 432 positive controls clear 3.0 and the cosmetic staples with no chromophore
# stay at 0", which is a claim about two numbers this file does not produce.
# Raising the weights until chlorpromazine cleared 3.0 would drag every benign
# bicyclic up with it; the chemotype table already answers that case correctly.
_W_TRICYCLIC_FUSED = 3.0
_W_BICYCLIC_FUSED = 0.8
_W_EXTRA_FUSED_RING = 0.8
_W_POLYENE = {2: 0.5, 3: 1.6, 4: 3.0}
_W_ARYL_KETONE = 2.0
_W_ARYL_CHAIN = {1: 0.8, 2: 2.0}
_W_PUSH_PULL = 1.5
_W_EXTRA_RING_IN_PI = 0.5
# Biaryl extension saturates: p-terphenyl and p-quaterphenyl absorb at
# almost the same wavelength, unlike the acene series where every fused
# ring adds 60-80 nm. Cap the contribution so a polyaryl drug scaffold
# cannot accumulate its way to a phototoxicity call.
_W_EXTRA_RING_CAP = 1.5
# ...but that saturation argument holds only for a symmetric homo-biaryl. A
# donor-acceptor or heteroaryl biaryl carries a charge-transfer band along the
# biaryl axis and does not saturate, which is why drometrizole, ensulizole,
# bemotrizinol and bisoctrizole are licensed UV filters. Under the flat cap
# all four came back "low": no photosafety attention needed, for molecules
# whose regulatory purpose is to absorb in the sunlight window.
_W_CT_BIARYL_RING = 1.5
_W_CT_BIARYL_CAP = 3.0
# A donor-only auxochrome on an extended aromatic system. Smaller than
# push-pull because one substituent shifts the band less than two, but not
# zero: 2-naphthol, 2-methoxynaphthalene and naproxen all reach 326-331 nm.
_W_DONOR_EXTENDED = 0.5

HIGH, MODERATE, LOW = 3.0, 2.0, 1.0


def extinction_proxy(mol: Chem.Mol) -> tuple[float, dict, list[str]]:
    """Score the likelihood that epsilon(290-700 nm) exceeds 1000 L/mol/cm.

    Returns (score, features, human-readable reasons). Bands: >=3.0 very
    likely above the threshold, 2.0-3.0 borderline, 1.0-2.0 unlikely, below
    1.0 no chromophore at all.
    """
    f = features(mol)
    score = 0.0
    why: list[str] = []

    n_fused = f["n_fused_pi_rings"]
    if n_fused >= 3:
        score += _W_TRICYCLIC_FUSED + _W_EXTRA_FUSED_RING * (n_fused - 3)
        why.append(f"{n_fused} fused rings in one system; linear acenes red-shift "
                   "roughly 60-80 nm per added ring")
    elif n_fused == 2:
        score += _W_BICYCLIC_FUSED
        why.append("bicyclic conjugated system; naphthalene- and purine-type "
                   "bands peak near 260-275 nm, below the 290 nm sunlight cut-off")

    extra_rings = max(0, f["n_rings_in_pi_system"] - max(n_fused, 1))
    if extra_rings and f["charge_transfer_biaryl"]:
        score += min(_W_CT_BIARYL_RING * extra_rings, _W_CT_BIARYL_CAP)
        why.append(f"{extra_rings} further ring(s) joined across a "
                   "donor-acceptor or heteroaryl biaryl bond; the "
                   "charge-transfer band along that axis does not saturate "
                   "the way a p-terphenyl chain does")
    elif extra_rings:
        score += min(_W_EXTRA_RING_IN_PI * extra_rings, _W_EXTRA_RING_CAP)
        why.append(f"{extra_rings} further ring(s) conjugated into the same pi "
                   "system (biaryl/flavonoid-type extension)")

    pl = f["polyene_length"]
    if pl >= 2:
        score += _W_POLYENE.get(min(pl, 4), _W_POLYENE[4])
        why.append(f"{pl} conjugated C=C in a chain; each added double bond "
                   "moves the band about 30 nm to the red")

    n_ext = f["extending_pi_bonds"]
    if n_ext >= 1:
        score += _W_ARYL_CHAIN.get(min(n_ext, 2), _W_ARYL_CHAIN[2])
        why.append(f"{n_ext} non-amide pi bond(s) conjugated into the aromatic "
                   "system; cinnamate- and chalcone-type extension moves the "
                   "band from 250 nm to past 300 nm")

    if f["aryl_ketone"]:
        score += _W_ARYL_KETONE
        why.append("carbonyl conjugated to an aromatic ring and not an amide; "
                   "the n->pi* band sits near 330 nm and populates a triplet")

    if f["donor_on_extended_system"]:
        score += _W_DONOR_EXTENDED
        why.append("donor auxochrome on a conjugated system of two or more "
                   "rings; naphthalene peaks at 275 nm and 2-methoxy-, "
                   "2-hydroxy- and 2-alkylnaphthalenes at 326-331 nm")

    if f["push_pull"]:
        score += _W_PUSH_PULL
        why.append("donor and acceptor on the same pi system; the charge-transfer "
                   "band is far to the red of the parent (4-nitroaniline, 380 nm)")

    return round(score, 2), f, why


def band(score: float) -> tuple[str, str]:
    """Map the proxy score to a concern level and an epsilon statement."""
    if score >= HIGH:
        return "high", "epsilon(290-700 nm) very likely above 1000 L/mol/cm"
    if score >= MODERATE:
        return "moderate", "epsilon(290-700 nm) plausibly above 1000 L/mol/cm"
    if score >= LOW:
        return "low", "weak chromophore; epsilon(290-700 nm) likely below 1000 L/mol/cm"
    return "none", "no chromophore reaching 290 nm"
