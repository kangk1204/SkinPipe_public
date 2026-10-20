'Refuse to report an MM/GBSA number for a site whose metal the model has lost.'
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Iterable, Optional, Sequence

from ..structure import metal_sites

LOGGER = logging.getLogger(__name__)

#: Refuse (default): raise :class:`MetalSiteError` instead of returning a dG.
REFUSE = "refuse"
#: Score anyway and stamp every result with the fact that the metal is missing.
#: Legitimate for a congeneric series that does not touch the metal; the name is
#: deliberately long so it cannot be set and then forgotten.
SCORE_WITHOUT_METAL = "score_without_metal"
POLICIES = (REFUSE, SCORE_WITHOUT_METAL)









SITE_CUTOFF_A: float = 6.0


class MetalSiteError(RuntimeError):
    """The pose binds a metal that the MM/GBSA model does not contain.

    Raised rather than returned. A metal-stripped site is not a failed pose, it
    is the wrong question, and a caller that only checks ``dG_bind is not None``
    would otherwise file it next to a ligand that failed to parameterise and
    carry on. :func:`skinpipe.mmgbsa.openmm_gbsa.mmgbsa_single` swallows every
    other exception into ``GBSAResult.notes``; this one it re-raises.
    """


def ligand_coords(sdf_path: Path) -> Optional["Any"]:
    """Heavy-atom coordinates of a pose as an (N, 3) array, or ``None``.

    ``sanitize=False`` on purpose: the guard has to run before anything
    expensive, on a file that may well turn out to be unparameterisable, and all
    it needs is where the atoms are. ``None`` means the file could not be read
    at all, which the caller must treat as "cannot rule the metal out" rather
    than as "no metal nearby".
    """
    from rdkit import Chem, RDLogger

    RDLogger.DisableLog("rdApp.*")
    path = Path(sdf_path)
    if not path.is_file():
        return None
    mol = Chem.MolFromMolFile(str(path), sanitize=False, removeHs=False)
    if mol is None or mol.GetNumConformers() == 0:
        return None
    conf = mol.GetConformer()
    points = [list(conf.GetAtomPosition(atom.GetIdx()))
              for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    if not points:
        return None
    return _finite_coordinates(points, f"pose {path.name}")


def _finite_coordinates(points: "Any", context: str) -> "Any":
    """Geometry must be known before a site can be declared clear."""
    import numpy as np

    try:
        xyz = np.asarray(points, dtype=float).reshape(-1, 3)
    except (ValueError, TypeError) as exc:
        raise ValueError(f"invalid {context} coordinates") from exc
    if not len(xyz) or not np.isfinite(xyz).all():
        raise ValueError(f"{context} coordinates must be nonempty and finite")
    return xyz


def _distances(metals: Sequence[metal_sites.MetalAtom],
               ligand_xyz: "Any") -> list[float]:
    import numpy as np

    lig = _finite_coordinates(ligand_xyz, "ligand")
    positions = _finite_coordinates([m.xyz for m in metals], "metal")
    distances = [float(np.linalg.norm(lig - xyz, axis=1).min()) for xyz in positions]
    if not np.isfinite(distances).all():
        raise ValueError("metal/ligand distances must be finite")
    return distances


def _describe(metals: Sequence[metal_sites.MetalAtom],
              distances: Optional[Sequence[float]]) -> str:
    if distances is None:
        return ", ".join(m.label for m in metals)
    return ", ".join(f"{m.label} at {d:.1f} A" for m, d in zip(metals, distances))


def check_binding_site(receptor_pdb: Path, ligand_xyz: Optional["Any"],
                       policy: str = REFUSE,
                       context: str = "") -> list[str]:
    """Decide whether an MM/GBSA number for this pose would mean anything.

    ``ligand_xyz`` is an (N, 3) array of pose heavy atoms from
    :func:`ligand_coords`, or ``None`` when the pose could not be read -- in
    which case every metal is treated as in-site, because an unreadable pose
    cannot be shown to be clear of one.

    Returns notes to attach to the result. Raises :class:`MetalSiteError` under
    the default ``policy=REFUSE`` when a metal is within :data:`SITE_CUTOFF_A`.
    A receptor with no metals at all returns ``[]`` without reading anything
    else, so this costs nothing on the ordinary path.
    """
    if policy not in POLICIES:
        raise ValueError(f"unknown on_metal_site policy {policy!r}; "
                         f"choose from {list(POLICIES)}")
    try:
        metals = metal_sites.receptor_metals(receptor_pdb)
    except ValueError as exc:
        raise MetalSiteError(f"cannot determine metal site: {exc}") from exc
    if not metals:
        return []

    name = Path(receptor_pdb).name
    where = f" for {context}" if context else ""
    if ligand_xyz is None:
        near, far = list(metals), []
        reason = ("the pose could not be read, so its distance to them is "
                  "unknown and they cannot be ruled out")
    else:
        try:
            pairs = list(zip(metals, _distances(metals, ligand_xyz)))
        except ValueError as exc:
            raise MetalSiteError(f"cannot determine metal site: {exc}") from exc
        near = [m for m, d in pairs if d <= SITE_CUTOFF_A]
        far = [(m, d) for m, d in pairs if d > SITE_CUTOFF_A]
        reason = (f"{len(near)} of them within {SITE_CUTOFF_A:.0f} A of the pose: "
                  + _describe(near, [d for m, d in pairs if d <= SITE_CUTOFF_A]))

    notes: list[str] = []
    if far:
        notes.append(
            f"{len(far)} catalytic metal atom(s) absent from the MM/GBSA model "
            f"but more than {SITE_CUTOFF_A:.0f} A from the pose "
            f"({_describe([m for m, _ in far], [d for _, d in far])}); scored through")

    if not near:
        return notes

    message = (
        f"{name} has {len(metals)} catalytic metal atom(s) that the MM/GBSA model "
        f"does not contain, {reason}. amber14+GBn2 has no parameters for them "
        "(see skinpipe.mmgbsa.metal_guard for the measurements), so a dG computed "
        f"here{where} would be for the apo site -- which for a chelating inhibitor "
        "is a different molecule's question. Refusing to report a number. Pass "
        f"on_metal_site={SCORE_WITHOUT_METAL!r} to get the apo number anyway, and "
        "label it as such wherever it is reported.")

    if policy == REFUSE:
        LOGGER.error("%s", message)
        raise MetalSiteError(message)

    stamp = (f"SCORED WITHOUT METAL: {_describe(near, None)} "
             f"{'is' if len(near) == 1 else 'are'} in the binding site but absent "
             "from the force field; this dG is for the apo site")
    LOGGER.warning("%s%s", stamp, where)
    return notes + [stamp]


def check_poses(receptor_pdb: Path, pose_sdfs: Iterable[Path],
                policy: str = REFUSE) -> list[str]:
    'Guard a whole batch up front: one receptor, every pose, before any work.'
    try:
        metals = metal_sites.receptor_metals(receptor_pdb)
    except ValueError as exc:
        raise MetalSiteError(f"cannot determine metal site: {exc}") from exc
    if not metals:
        return []
    import numpy as np

    try:
        readable = [xyz for xyz in (ligand_coords(p) for p in pose_sdfs)
                    if xyz is not None]
    except ValueError as exc:
        raise MetalSiteError(f"cannot determine metal site: {exc}") from exc
    stacked = np.vstack(readable) if readable else None
    return check_binding_site(
        receptor_pdb, stacked, policy=policy,
        context=f"this batch ({len(readable)} readable pose(s))")
