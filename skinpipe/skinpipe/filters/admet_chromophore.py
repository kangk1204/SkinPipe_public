"""Auxochrome and charge-transfer analysis for the `admet_photo` proxy.

Split out of `admet_photo` for file size; it answers three questions about a
conjugated system that the ring/polyene count alone gets wrong.

1. Donors and acceptors
-----------------------
A ring carrying both a donor and an acceptor has a charge-transfer band tens
of nanometres to the red of the parent. Amide carbonyls are excluded from the
acceptor set: the nitrogen lone pair feeds the same carbon, so a benzamide is
a weak acceptor and its band stays in the deep UV. Without that exclusion
every anilide and every kinase inhibitor with a benzamide reads as a push-pull
chromophore, imatinib among them.

The acceptor set had the aromaticity bug this project documents in
`admet_alerts`: it was written entirely in aliphatic notation, so no acceptor
alternative could ever match a ring nitrogen RDKit had aromatised. Drometrizole
- a licensed UV filter whose phenol donor and benzotriazole acceptor are the
whole point of the molecule - came back with push_pull False for that reason.

2. Donor on an extended system
------------------------------
`push_pull` requires both a donor and an acceptor, so a donor-only auxochrome
earned nothing even though it is the classic red-shift. Naphthalene peaks at
275 nm and is genuinely invisible to sunlight, but 2-methoxynaphthalene sits at
326 nm, 2-naphthol at 330 nm, and naproxen at 331 nm - and naproxen is a
documented photosensitiser. A donor on a two-ring-or-larger conjugated system
gets credit for that shift.

3. Charge-transfer biaryls
--------------------------
Extension across a biaryl single bond saturates: p-terphenyl and p-quaterphenyl
absorb at almost the same wavelength, which is why `admet_photo` caps the
per-ring contribution. That is true of homo-biaryls and wrong for the two cases
that matter commercially. A donor-acceptor biphenyl has a charge-transfer band
along the biaryl axis, and a heteroaryl biaryl has one because the two rings
differ in electron density to begin with: 2-(2-hydroxyphenyl)benzotriazoles
(drometrizole, bisoctrizole), 2-phenylbenzimidazoles (ensulizole) and the
triazine trisaryls (bemotrizinol) are all commercial UV absorbers *because* of
that chromophore, and all four were banded "low" under the flat cap.

Detected by ring analysis rather than SMARTS, following the house rule in
`admet_alerts`: ring topology is not something a pattern set enumerates.
"""
from __future__ import annotations

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

#: Donors push electron density into the ring (auxochromes).
DONOR = Chem.MolFromSmarts(
    "[$([OX2H][c]),$([OX2]([#6])[c]),$([NX3;H2,H1,H0;!$(N[CX3]=[OX1]);"
    "!$(N[SX4]);!$(N=*)][c]),$([SX2][c])]")

#: Acceptors pull it out. The last alternative is the aromatised azole case:
#: an aryl joined across a non-ring bond to a pyrrole-type ring nitrogen that
#: is itself flanked by a pyridine-type nitrogen, i.e. an N-aryl azole that is
#: electron-poor (benzotriazole, triazole, tetrazole). Written `!@` so it can
#: only ever describe a bond BETWEEN two rings: an in-ring match would make
#: every purine its own acceptor and would band adenosine as a chromophore.
ACCEPTOR = Chem.MolFromSmarts(
    "[$([c][$([NX3](=O)=O),$([NX3+](=O)[O-])]),$([c][CX3](=[OX1])[#6,#8]),"
    "$([c][CX2]#[NX1]),$([c][SX4](=[OX1])=[OX1]),$([c][NX2]=[NX2]),"
    "$([c][CX3]=[CX3][CX3](=[OX1])[#6,#8]),"     # vinyl spacer: the cinnamate case
    "$([c]!@[nX3](:[nX2]))]")                    # N-aryl azole: drometrizole


def donor_acceptor_atoms(mol: Chem.Mol) -> tuple[set[int], set[int]]:
    """Indices of the atoms that carry a donor / an acceptor substituent."""
    return ({m[0] for m in mol.GetSubstructMatches(DONOR)},
            {m[0] for m in mol.GetSubstructMatches(ACCEPTOR)})


def aromatic_rings_within(mol: Chem.Mol, component: set[int]) -> list[set[int]]:
    """Fully aromatic SSSR rings lying inside one conjugated component."""
    return [set(r) for r in mol.GetRingInfo().AtomRings()
            if set(r) <= component
            and all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in r)]


def _is_heteroaromatic(mol: Chem.Mol, ring: set[int]) -> bool:
    return any(mol.GetAtomWithIdx(i).GetAtomicNum() != 6 for i in ring)


def charge_transfer_biaryl(mol: Chem.Mol, component: set[int],
                           push_pull: bool) -> bool:
    """True when two aromatic rings are joined by a ring-to-ring single bond
    and the pair is not a symmetric homo-biaryl.

    A fused bond is inside a ring and is skipped, so naphthalene, the purines
    and every acene are untouched; only a genuine inter-ring bond qualifies.
    Biphenyl itself does not: with no heteroaromatic partner and no donor plus
    acceptor, it is the saturating case the cap was written for.
    """
    rings = aromatic_rings_within(mol, component)
    if len(rings) < 2:
        return False
    for bond in mol.GetBonds():
        if bond.IsInRing():
            continue
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        for ring_a in (r for r in rings if a in r):
            for ring_b in (r for r in rings if b in r and r != ring_a):
                if (push_pull or _is_heteroaromatic(mol, ring_a)
                        or _is_heteroaromatic(mol, ring_b)):
                    return True
    return False
