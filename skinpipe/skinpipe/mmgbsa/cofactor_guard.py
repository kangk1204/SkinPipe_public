'Refuse an MM/GBSA number for a pocket whose cofactor the model has lost.'
from __future__ import annotations

import logging
import math
from pathlib import Path
from typing import Any, Iterable, Sequence

from ..structure import cofactors as cofactor_sites
from .metal_guard import _finite_coordinates, ligand_coords

log = logging.getLogger(__name__)


class CofactorSiteError(RuntimeError):
    """The pocket these poses were docked into is not the one being scored.

    Mirrors :class:`~skinpipe.mmgbsa.metal_guard.MetalSiteError`: it is not a
    failed pose, it is the wrong receptor, so a caller must not record it per
    row and move on.
    """


#: A cofactor heavy atom this close to any pose heavy atom is a wall of the
#: pocket that pose occupies. Vina's own interaction terms run out by about 8 A;
#: 4.5 is the range over which removing a group changes the shape of the site
#: rather than only its long-range electrostatics.
CONTACT_A = 4.5


def _contacting(cofactors: Sequence["cofactor_sites.CofactorAtom"],
                ligand_xyz: Any, cutoff: float) -> dict[str, float]:
    """Nearest approach per cofactor group, for those within ``cutoff``."""
    import numpy as np

    if not math.isfinite(cutoff) or cutoff < 0:
        raise ValueError("cofactor cutoff must be nonnegative and finite")
    lig = _finite_coordinates(ligand_xyz, "ligand")
    positions = _finite_coordinates([atom.xyz for atom in cofactors], "cofactor")
    nearest: dict[str, float] = {}
    for atom, xyz in zip(cofactors, positions):
        d = float(np.linalg.norm(lig - xyz, axis=1).min())
        if not math.isfinite(d):
            raise ValueError("cofactor/ligand distances must be finite")
        if d <= cutoff and d < nearest.get(atom.label, float("inf")):
            nearest[atom.label] = d
    return nearest


def check_poses(receptor_pdb: Path, pose_sdfs: Iterable[Path],
                cutoff: float = CONTACT_A) -> list[str]:
    'Reasons this batch must not be MM/GBSA rescored. Empty means go ahead.'
    if not math.isfinite(cutoff) or cutoff < 0:
        raise ValueError("cofactor cutoff must be nonnegative and finite")
    try:
        cofactors = cofactor_sites.receptor_cofactors(receptor_pdb)
    except ValueError as exc:
        raise CofactorSiteError(f"cannot determine cofactor site: {exc}") from exc
    if not cofactors:
        return []

    import numpy as np

    try:
        readable = [xyz for xyz in (ligand_coords(p) for p in pose_sdfs)
                    if xyz is not None]
    except ValueError as exc:
        raise CofactorSiteError(f"cannot determine cofactor site: {exc}") from exc
    groups = sorted({a.label for a in cofactors})
    if not readable:
        return [f"{len(groups)} structural cofactor(s) ({', '.join(groups)}) "
                f"are in the docking receptor and not in {Path(receptor_pdb).name}, "
                "and no pose could be read to tell whether they line the pocket. "
                "MM/GBSA would be scored against a different structure from the "
                "one that produced these poses."]

    try:
        near = _contacting(cofactors, np.vstack(readable), cutoff)
    except ValueError as exc:
        raise CofactorSiteError(f"cannot determine cofactor site: {exc}") from exc
    if not near:
        return []
    where = ", ".join(f"{label} at {d:.1f} A" for label, d in sorted(near.items()))
    return [f"{len(near)} structural cofactor(s) form a wall of this pocket "
            f"({where}) and are absent from {Path(receptor_pdb).name}; "
            "`mmgbsa.system.prepare_receptor` removes heterogens, so the "
            "rescore would run against a site the poses were not docked into. "
            "Docking kept them, MM/GBSA cannot represent them, and a number "
            "from the second structure is not a refinement of the first."]


def refuse_poses(receptor_pdb: Path, pose_sdfs: Iterable[Path],
                 cutoff: float = CONTACT_A) -> None:
    'Raise :class:`CofactorSiteError` if this batch must not be rescored.'
    reasons = check_poses(receptor_pdb, pose_sdfs, cutoff)
    if reasons:
        raise CofactorSiteError(" ".join(reasons))
