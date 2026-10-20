'The three similarity measurements a scaffold hop is judged on.'
from __future__ import annotations

import logging
from typing import Optional

from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, rdFingerprintGenerator, rdMolAlign, rdShapeHelpers
from rdkit.Chem.Pharm2D import Generate, Gobbi_Pharm2D
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")
log = logging.getLogger(__name__)

CONFORMER_SEED = 0xC05          # fixed so every score in this module is reproducible
N_CONFORMERS_DEFAULT = 10       # justified in the module docstring
_PRUNE_RMS = 0.5                # A; collapses redundant conformers on rigid molecules
_MMFF_ITERS = 200

_MORGAN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
_PHARM_FACTORY = Gobbi_Pharm2D.factory


# --------------------------------------------------------------- scaffolds

def murcko_scaffold(smiles: str, generic: bool = False) -> str:
    """The Bemis-Murcko scaffold: ring systems plus the linkers joining them.

    Bemis & Murcko (1996, J Med Chem 39:2887) defined it by stripping every
    side chain and keeping what is left. It is the standard operational answer
    to "is this the same chemical series?", and it is the first gate in
    `hop_score`: two molecules sharing a Murcko scaffold are analogues of each
    other, whatever their fingerprints say.

    With `generic=True` the framework is returned instead -- every atom made
    carbon and every bond made single, so only the ring-and-linker topology
    survives. That distinguishes a real topological change from an element
    swap: pyridine and benzene are different Murcko scaffolds but the *same*
    generic framework, so a phenyl-to-pyridine edit has not changed the shape
    of the molecule at all, only its electronics.

    Returns "" for an acyclic molecule, which has no Murcko scaffold, and for
    an unparseable one. Callers must not treat "" == "" as a scaffold match;
    `hop_score` handles this explicitly.

    >>> murcko_scaffold("OCC1=CC(=O)C(O)=CO1")            # kojic acid
    'O=c1ccocc1'
    >>> murcko_scaffold("OCC1=CC(=O)C(O)=CO1", generic=True)
    'CC1CCCCC1'
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        log.warning("murcko_scaffold: unparseable SMILES %r", smiles[:80])
        return ""
    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        if generic:
            scaffold = MurckoScaffold.MakeScaffoldGeneric(scaffold)
        return Chem.MolToSmiles(scaffold)
    except Exception as exc:
        log.warning("murcko_scaffold failed on %s: %s", smiles[:60], exc)
        return ""


# --------------------------------------------------------------------- 2D

def tanimoto_2d(smiles_a: str, smiles_b: str) -> float:
    """ECFP4 Tanimoto. LOW is what a scaffold hop needs. See module docstring."""
    a, b = Chem.MolFromSmiles(smiles_a), Chem.MolFromSmiles(smiles_b)
    if a is None or b is None:
        return float("nan")
    return round(DataStructs.TanimotoSimilarity(
        _MORGAN.GetFingerprint(a), _MORGAN.GetFingerprint(b)), 4)


def pharmacophore_similarity(smiles_a: str, smiles_b: str) -> float:
    """Gobbi 2D pharmacophore fingerprint Tanimoto. HIGH is what a hop needs."""
    a, b = Chem.MolFromSmiles(smiles_a), Chem.MolFromSmiles(smiles_b)
    if a is None or b is None:
        return float("nan")
    try:
        fp_a = Generate.Gen2DFingerprint(a, _PHARM_FACTORY)
        fp_b = Generate.Gen2DFingerprint(b, _PHARM_FACTORY)
    except Exception as exc:
        log.warning("pharmacophore fingerprint failed: %s", exc)
        return float("nan")
    if not fp_a.GetNumOnBits() or not fp_b.GetNumOnBits():
        # A molecule with no recognised pharmacophoric feature (a plain
        # hydrocarbon) has nothing to compare. Returning 0.0 here would read
        # as "no pharmacophore in common", which is true but for the wrong
        # reason, so say it out loud.
        log.info("one molecule has no Gobbi pharmacophore features; similarity is 0 by "
                 "absence, not by disagreement")
        return 0.0
    return round(DataStructs.TanimotoSimilarity(fp_a, fp_b), 4)


# --------------------------------------------------------------------- 3D

def embed_conformers(smiles: str, n_confs: int = N_CONFORMERS_DEFAULT,
                     seed: int = CONFORMER_SEED) -> Optional[Chem.Mol]:
    """Hydrogen-added, MMFF-minimised, RMSD-pruned conformer ensemble.

    Hydrogens are added because MMFF minimisation and Open3DAlign's atom typing
    both need them; they are not, however, part of the shape measurement.
    RDKit's ``ShapeTanimotoDist`` defaults to ``ignoreHs=True``, so
    `shape_tanimoto` is a **heavy-atom** volume overlap. That default is kept
    rather than overridden: hydrogen positions on a rotatable group are a
    property of the conformer generator, not of the molecule, and letting them
    into the grid adds sampling noise to a measure that is already a lower
    bound.

    Returns None when embedding fails, which happens on strained cages and on
    some of what BRICS recombination produces -- a legitimate outcome, not an
    error to suppress.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.pruneRmsThresh = _PRUNE_RMS
    try:
        cids = AllChem.EmbedMultipleConfs(mol, numConfs=n_confs, params=params)
    except Exception as exc:
        log.debug("embedding raised on %s: %s", smiles[:60], exc)
        return None
    if not len(cids):
        log.debug("could not embed %s", smiles[:60])
        return None
    try:
        AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=_MMFF_ITERS)
    except Exception:
        # No MMFF parameters for some atom type. The raw ETKDG geometry is
        # still a reasonable conformer, so keep going rather than dropping the
        # molecule -- but O3A also needs MMFF, so shape_tanimoto will fail
        # downstream and be reported as missing.
        log.debug("MMFF minimisation unavailable for %s", smiles[:60])
    return mol


def shape_tanimoto(mol_a: Chem.Mol, mol_b: Chem.Mol) -> Optional[float]:
    """Best shape overlap over all conformer pairs, after Open3DAlign.

    1.0 is identical volume, 0.0 is no overlap. HIGH is what a hop needs: it
    is the claim that the new molecule can occupy the same space in the
    pocket. The cost is O(n_a * n_b) alignments, which is why `rank_hops`
    prefilters on the cheap 2D gates before it ever calls this.

    Returns None if either molecule has no conformers or if O3A cannot be set
    up (no MMFF parameters). None means "not measured", and `hop_score` scores
    it as missing rather than as zero -- a molecule that failed to embed is
    not the same as one with no shape in common.
    """
    if mol_a is None or mol_b is None:
        return None
    if not mol_a.GetNumConformers() or not mol_b.GetNumConformers():
        return None
    props_a = AllChem.MMFFGetMoleculeProperties(mol_a)
    props_b = AllChem.MMFFGetMoleculeProperties(mol_b)
    if props_a is None or props_b is None:
        log.debug("no MMFF parameters; cannot run O3A")
        return None

    best = 0.0
    measured = False
    for cid_a in range(mol_a.GetNumConformers()):
        for cid_b in range(mol_b.GetNumConformers()):
            try:
                o3a = rdMolAlign.GetO3A(mol_b, mol_a, props_b, props_a,
                                        prbCid=cid_b, refCid=cid_a)
                o3a.Align()
                overlap = 1.0 - rdShapeHelpers.ShapeTanimotoDist(
                    mol_b, mol_a, confId1=cid_b, confId2=cid_a)
            except Exception as exc:
                log.debug("O3A failed on conformer pair (%d, %d): %s", cid_a, cid_b, exc)
                continue
            measured = True
            best = max(best, overlap)
    return round(best, 4) if measured else None


def best_shape_from_smiles(smiles_a: str, smiles_b: str,
                           n_confs: int = N_CONFORMERS_DEFAULT) -> Optional[float]:
    """Convenience wrapper: embed both molecules and return the best overlap.

    Prefer `embed_conformers` + `shape_tanimoto` when scoring many candidates
    against one parent, so the parent's ensemble is built once.
    """
    return shape_tanimoto(embed_conformers(smiles_a, n_confs),
                          embed_conformers(smiles_b, n_confs))
