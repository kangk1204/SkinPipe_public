'Metal-site repair for receptors whose metals were held out of PDBFixer.'
from __future__ import annotations

import logging
from typing import Optional

import numpy as np

log = logging.getLogger(__name__)

#: AutoDock element token per metal, written into PDBQT columns 78-79. Vina
#: 1.2 maps every one of these onto its single internal ``Met_D`` type, so the
#: oxidation-state variants share their element's token.
METAL_AD_TYPE: dict[str, str] = {
    "ZN": "Zn", "CU": "Cu", "FE": "Fe", "MG": "Mg",
    "MN": "Mn", "CA": "Ca", "NI": "Ni", "CO": "Co",
    "CU1": "Cu", "CU3": "Cu", "FE2": "Fe", "MN3": "Mn",
    "3CO": "Co", "3NI": "Ni",
}
#: Formal charge used for the PDBQT partial-charge column. Unused by the Vina
#: scoring function, which is atom-type based, but read by the AD4 one.
METAL_CHARGE: dict[str, float] = {
    "ZN": 2.0, "CU": 2.0, "FE": 3.0, "MG": 2.0,
    "MN": 2.0, "CA": 2.0, "NI": 2.0, "CO": 2.0,
    "CU1": 1.0, "CU3": 3.0, "FE2": 2.0, "MN3": 3.0,
    "3CO": 3.0, "3NI": 3.0,
}


def metal_pdbqt_line(serial: int, src: str) -> str:
    """Render one crystallographic metal HETATM as an AutoDock PDBQT record."""
    resname = src[17:20].strip().upper()
    ad_type = METAL_AD_TYPE[resname]
    charge = METAL_CHARGE[resname]
    name = f"{ad_type.upper():<4}"[:4]
    coords = src[30:54]
    return (
        f"HETATM{serial:5d} {name} {resname:>3s} {src[21]}{src[22:26]}"
        f"{src[26]}   {coords}  1.00  0.00    {charge:+6.3f} {ad_type}"
    )


def render_metal_records(metal_lines: list[str],
                         start_serial: int = 90000) -> tuple[list[str], list[str]]:
    """Render retained metals as PDBQT HETATM records, plus honest notes.

    The note counts what actually reached the file. Reporting
    ``len(metal_lines)`` would claim a metal was kept in exactly the case where
    it was dropped for want of an AutoDock element type -- the one case worth
    knowing about, because a metalloenzyme pocket without its metal scores
    garbage that looks like a score.
    """
    records: list[str] = []
    notes: list[str] = []
    written: list[str] = []
    skipped: list[str] = []
    serial = start_serial
    for src in metal_lines:
        resname = src[17:20].strip().upper()
        if resname not in METAL_AD_TYPE:
            skipped.append(resname)
            continue
        serial += 1
        written.append(resname)
        records.append(metal_pdbqt_line(serial, src))
    if written:
        notes.append(f"kept {len(written)} catalytic metal atom(s): "
                     f"{', '.join(sorted(set(written)))}")
    if skipped:
        notes.append(f"DROPPED {len(skipped)} metal atom(s) with no AutoDock element "
                     f"type: {', '.join(sorted(set(skipped)))} - the pocket is missing "
                     "them and scores from it are not trustworthy")
        log.warning("no AutoDock type for metal(s) %s", sorted(set(skipped)))
    return records, notes


#: Longest N-metal distance still counted as coordination, Angstrom.  His-Cu
#: and His-Zn are 2.0-2.2 A; Mn and Ca reach ~2.4 A.  Beyond 2.6 A the
#: nitrogen is not in the first coordination shell.
COORDINATION_CUTOFF_A: float = 2.6

#: N-H bond length used when a ring proton is rebuilt, Angstrom.
_NH_BOND_A: float = 1.01

#: (ring nitrogen, its proton, the two ring neighbours of that nitrogen)
_RING_N = {
    "ND1": ("HD1", ("CG", "CE1")),
    "NE2": ("HE2", ("CD2", "CE1")),
}


def _xyz(line: str) -> np.ndarray:
    return np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])], dtype=float)


def _with_position(line: str, name: str, point: np.ndarray) -> str:
    """Rewrite an ATOM record's name and coordinates, leaving every other column."""
    padded = f" {name:<3}"[:4]
    return (f"{line[:12]}{padded}{line[16:30]}"
            f"{point[0]:8.3f}{point[1]:8.3f}{point[2]:8.3f}{line[54:]}")


def _rebuilt_proton(n_pos: np.ndarray, neighbours: list[np.ndarray]) -> Optional[np.ndarray]:
    """Place a ring N-H in the imidazole plane, bisecting the external angle."""
    if len(neighbours) != 2:
        return None
    arms = []
    for nb in neighbours:
        v = n_pos - nb
        norm = float(np.linalg.norm(v))
        if norm < 1e-6:
            return None
        arms.append(v / norm)
    direction = arms[0] + arms[1]
    norm = float(np.linalg.norm(direction))
    if norm < 1e-6:
        return None
    return n_pos + _NH_BOND_A * direction / norm


def fix_metal_site_histidines(
    protein_lines: list[str], metal_lines: list[str],
    cutoff: float = COORDINATION_CUTOFF_A,
) -> tuple[list[str], list[str]]:
    """Move ring protons off metal-coordinating histidine nitrogens.

    ``protein_lines`` is the protonated receptor (PDBFixer output, ATOM records
    with hydrogens); ``metal_lines`` are the catalytic metals held out of it.
    Returns ``(lines, notes)`` with one note per residue touched or skipped.
    See the module docstring for why this is necessary and why the proton is
    moved rather than deleted.
    """
    if not metal_lines or not protein_lines:
        return protein_lines, []

    metals = [_xyz(line) for line in metal_lines]
    residues: dict[tuple[int, str, str], list[int]] = {}
    segment = 0
    for index, line in enumerate(protein_lines):
        if line.startswith("TER"):
            segment += 1
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        if line[17:20].strip().upper() != "HIS":
            continue
        residues.setdefault((segment, line[21], line[22:27]), []).append(index)

    drop: set[int] = set()
    notes: list[str] = []
    for (_segment, chain, resid), indices in residues.items():
        by_name = {protein_lines[i][12:16].strip().upper(): i for i in indices}
        coordinating = [
            name for name in _RING_N
            if name in by_name
            and min((float(np.linalg.norm(_xyz(protein_lines[by_name[name]]) - m))
                     for m in metals), default=1e9) <= cutoff
        ]
        if not coordinating:
            continue

        label = f"HIS {chain}:{resid.strip()}"
        present = [n for n in _RING_N if _RING_N[n][0] in by_name]
        if len(coordinating) == 2:
            notes.append(
                f"{label} bridges two metals through both ring nitrogens; the correct "
                "species is the imidazolate, which has no residue template, so it was "
                "left protonated - the metal site is approximate here")
            continue

        bound = coordinating[0]
        bound_h = _RING_N[bound][0]
        if bound_h not in by_name:
            continue                                    # already the right tautomer

        if len(present) == 2:                           # HIP: just drop the bad one
            drop.add(by_name[bound_h])
            notes.append(f"{label} was doubly protonated and coordinates a metal through "
                         f"{bound}; removed {bound_h} (the metal-facing proton)")
            continue

        other = next(n for n in _RING_N if n != bound)
        if other not in by_name:
            continue
        neighbour_names = _RING_N[other][1]
        neighbours = [_xyz(protein_lines[by_name[n]]) for n in neighbour_names if n in by_name]
        point = _rebuilt_proton(_xyz(protein_lines[by_name[other]]), neighbours)
        if point is None:
            notes.append(f"{label} coordinates a metal through {bound} but its imidazole "
                         "geometry could not be read; tautomer left as PDBFixer set it")
            continue

        distance = min(float(np.linalg.norm(_xyz(protein_lines[by_name[bound_h]]) - m))
                       for m in metals)
        protein_lines = list(protein_lines)
        protein_lines[by_name[bound_h]] = _with_position(
            protein_lines[by_name[bound_h]], _RING_N[other][0], point)
        notes.append(
            f"{label} coordinates a metal through {bound} but PDBFixer protonated that "
            f"nitrogen ({bound_h} was {distance:.2f} A from the metal); moved the proton "
            f"to {other} so the coordinating nitrogen is an acceptor")

    if drop:
        protein_lines = [pdb_line for i, pdb_line in enumerate(protein_lines) if i not in drop]
    if notes:
        log.info("metal-site histidine repair: %s", "; ".join(notes))
    return protein_lines, notes


def repair_metal_sites(clean_pdb, metal_lines: list[str]) -> list[str]:
    """Apply :func:`fix_metal_site_histidines` to a PDBFixer output file in place.

    Returns the notes; the file is only rewritten when something changed.
    """
    from pathlib import Path

    path = Path(clean_pdb)
    text = path.read_text()
    lines = text.splitlines()
    # TER records carry the topology's observed fragment boundaries. Keep
    # every non-coordinate record in place, including between repaired atoms.
    fixed, notes = fix_metal_site_histidines(lines, metal_lines)
    if fixed != lines:
        path.write_text("\n".join(fixed) + "\n")
    return notes
