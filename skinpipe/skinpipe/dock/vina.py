"""AutoDock Vina docking: ligand preparation, screening, and redock control.

Scientific choices
------------------
**Conformer generation.** ETKDGv3 with a fixed seed, then MMFF94s (UFF fallback
for atoms MMFF cannot type, e.g. boron or most transition metals). Vina searches
torsions itself, so the input conformer only fixes what Vina cannot change: ring
puckers, bridgehead geometry and stereocentres. Generating several conformers and
keeping the lowest-energy one therefore buys ring quality, not torsional
coverage, which is why ``n_conformers`` defaults to 1.

**One bad SMILES must not kill a screen.** ``prepare_ligand_pdbqt`` returns
``None`` and logs the reason for anything that will not sanitise, embed or
optimise. A virtual screen over a natural-product library will hit these
regularly (unparsable stereo, hypervalent nitrogen in a mis-typed alkaloid,
macrocycles that ETKDG cannot close).

**Ligand efficiency.** ``-affinity / heavy_atom_count``. Vina's score grows
roughly linearly with size, so raw affinity ranks big greasy molecules first -
exactly the compounds that will not cross the stratum corneum. LE is the
size-corrected column to sort on for a topical programme.

**The redock control.** Docking a ligand back into the structure it came from is
the only cheap check that the receptor preparation and the box are sane. The
metric is *in-place* symmetry-corrected RMSD (``rdMolAlign.CalcRMS``), not
``GetBestRMS``: ``GetBestRMS`` superimposes the two molecules before measuring,
which answers "is the conformer right" rather than "is it in the right place",
and reports small values even for a pose docked into the wrong subpocket. Both
numbers are returned so the difference is visible. The conventional pass mark is
RMSD <= 2.0 A on the top-scoring pose.
"""
from __future__ import annotations

from .engine import dock_many, dock_one
from .ligprep import prepare_ligand_pdbqt
from .redock import REDOCK_RMSD_PASS, redock_validation

__all__ = [
    "prepare_ligand_pdbqt",
    "dock_one",
    "dock_many",
    "redock_validation",
    "REDOCK_RMSD_PASS",
]
