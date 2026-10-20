"""Turn Vina's PDBQT output back into chemically valid SDF.

The awkward part of a Vina run is the exit, not the entry. PDBQT keeps only
AutoDock atom types and merges non-polar hydrogens into their carbons, so a
naive PDBQT -> SDF conversion (OpenBabel's default, for instance) guesses bond
orders from geometry and routinely aromatises or kekulises things wrongly. Any
downstream MM/GBSA stage then parameterises the wrong molecule.

Meeko solves this properly: ``MoleculePreparation`` stores the input RDKit
molecule's bond orders in ``REMARK SMILES`` records inside the PDBQT it writes,
and ``RDKitMolCreate.from_pdbqt_mol`` reads them back. The docked coordinates
are grafted onto the original chemistry instead of being re-perceived.

Every SDF this module writes is read straight back with ``Chem.SDMolSupplier``
and checked for a valid molecule with a matching heavy-atom count, because a
file that does not round-trip is worse than no file at all.
"""
from __future__ import annotations

import logging
import math
from pathlib import Path
from typing import Optional

from rdkit import Chem

log = logging.getLogger(__name__)


def parse_vina_remarks(pdbqt_text: str) -> list[tuple[float, float, float]]:
    """``(affinity, rmsd_lb, rmsd_ub)`` per pose, in pose order."""
    out: list[tuple[float, float, float]] = []
    for line in pdbqt_text.splitlines():
        if not line.startswith("REMARK VINA RESULT:"):
            continue
        parts = line.split(":", 1)[1].split()
        try:
            out.append((float(parts[0]), float(parts[1]), float(parts[2])))
        except (IndexError, ValueError):
            log.warning("unparsable VINA RESULT line: %s", line.strip())
    return out


def scored_pose_remarks(pdbqt_text: str) -> list[tuple[float, float, float]]:
    """Require one score for each complete MODEL, without dropping bad rows.

    ``parse_vina_remarks`` remains a permissive inspection API. Production
    conversion needs a one-to-one mapping: skipping a malformed first score
    would attach the second MODEL's affinity to the first MODEL's coordinates.
    Unwrapped single-pose output is supported, with exactly one score.
    Numeric values are preserved; the operational score window is applied
    separately, after conversion.
    """
    lines = pdbqt_text.splitlines()
    wrapped = any(line[:6].strip() == "MODEL" for line in lines)
    scores: list[tuple[float, float, float]] = []
    model_ids: set[int] = set()
    active: Optional[int] = None
    model_scores: list[tuple[float, float, float]] = []
    for line_number, line in enumerate(lines, start=1):
        record = line[:6].strip()
        if record == "MODEL":
            parts = line.split()
            try:
                model_id = int(parts[1]) if len(parts) == 2 else 0
            except ValueError:
                model_id = 0
            if active is not None or model_id <= 0 or model_id in model_ids:
                raise ValueError(f"invalid or duplicate MODEL at line {line_number}")
            model_ids.add(model_id)
            active, model_scores = model_id, []
        elif record == "ENDMDL":
            if active is None or len(model_scores) != 1:
                raise ValueError(f"MODEL must have exactly one VINA RESULT before line {line_number}")
            scores.append(model_scores[0])
            active = None
        elif line.startswith("REMARK VINA RESULT:"):
            if wrapped and active is None:
                raise ValueError(f"VINA RESULT outside MODEL at line {line_number}")
            parts = line.split(":", 1)[1].split()
            try:
                if len(parts) != 3:
                    raise ValueError("expected three values")
                score = (float(parts[0]), float(parts[1]), float(parts[2]))
            except ValueError as exc:
                raise ValueError(f"invalid VINA RESULT at line {line_number}") from exc
            if (not all(math.isfinite(value) for value in score[1:])
                    or not 0 <= score[1] <= score[2]):
                raise ValueError(f"invalid VINA RMSD bounds at line {line_number}")
            model_scores.append(score)
        elif wrapped and active is None and line.startswith(
                ("ATOM", "HETATM", "ROOT", "ENDROOT", "BRANCH", "ENDBRANCH", "TORSDOF")):
            raise ValueError(f"pose coordinates outside MODEL at line {line_number}")
    if active is not None:
        raise ValueError("MODEL is missing ENDMDL")
    if not wrapped:
        scores = model_scores
        if len(scores) != 1:
            raise ValueError("single-pose output must have exactly one VINA RESULT")
    return scores


def mols_from_docked_pdbqt(pdbqt_text: str) -> Optional[Chem.Mol]:
    """Rebuild the docked ligand, one conformer per pose, or ``None``."""
    from meeko import PDBQTMolecule, RDKitMolCreate

    try:
        pdbqt_mol = PDBQTMolecule(pdbqt_text, is_dlg=False, skip_typing=True)
        mols = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)
    except (ValueError, RuntimeError, KeyError) as exc:
        log.error("meeko could not rebuild molecule from docked PDBQT: %s", exc)
        return None
    mols = [m for m in mols if m is not None]
    if not mols:
        log.error("meeko returned no molecule from docked PDBQT")
        return None
    if len(mols) > 1:
        log.warning("docked PDBQT held %d molecules; using the first", len(mols))
    return mols[0]


def write_pose_sdf(
    mol: Chem.Mol,
    conf_id: int,
    out_sdf: Path,
    props: Optional[dict[str, object]] = None,
) -> Optional[Path]:
    """Write one conformer as SDF and verify it reloads. ``None`` on failure."""
    out_sdf.parent.mkdir(parents=True, exist_ok=True)
    writer = Chem.SDWriter(str(out_sdf))
    try:
        for key, value in (props or {}).items():
            mol.SetProp(str(key), str(value))
        writer.write(mol, confId=conf_id)
    except (ValueError, RuntimeError) as exc:
        log.error("could not write %s: %s", out_sdf, exc)
        return None
    finally:
        writer.close()

    if not sdf_roundtrips(out_sdf, expected_heavy_atoms=mol.GetNumHeavyAtoms(),
                          expected_graph=_flat_smiles(mol)):
        out_sdf.unlink(missing_ok=True)
        return None
    return out_sdf


def _flat_smiles(mol: Chem.Mol) -> Optional[str]:
    """Canonical SMILES without stereochemistry, or ``None`` if it cannot be made.

    Stereo is deliberately dropped from the comparison. RDKit assigns chiral
    tags from 3D coordinates when it reads a molblock back, so a ligand whose
    input SMILES left a centre undefined legitimately reloads *with* that centre
    assigned -- comparing isomeric SMILES would throw away a valid pose. What is
    worth comparing is the bond-order graph, which is exactly what a naive
    PDBQT conversion gets wrong and what an atom count cannot see.
    """
    try:
        return Chem.MolToSmiles(Chem.RemoveHs(Chem.Mol(mol)), isomericSmiles=False)
    except (ValueError, RuntimeError):
        return None


def sdf_roundtrips(sdf_path: Path, expected_heavy_atoms: Optional[int] = None,
                   expected_graph: Optional[str] = None) -> bool:
    """True when the file reloads as the same molecule.

    ``expected_graph`` is a stereo-free canonical SMILES. Checking it is what
    makes this function live up to its own module docstring: an atom count is
    unchanged by a benzene read back as cyclohexatriene, and that is precisely
    the corruption a geometry-based PDBQT conversion produces.
    """
    supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=False, sanitize=True)
    mols = [m for m in supplier if m is not None]
    if not mols:
        log.error("SDF %s does not round-trip: no valid molecule", sdf_path)
        return False
    if expected_heavy_atoms is not None:
        got = mols[0].GetNumHeavyAtoms()
        if got != expected_heavy_atoms:
            log.error(
                "SDF %s round-trips with %d heavy atoms, expected %d",
                sdf_path,
                got,
                expected_heavy_atoms,
            )
            return False
    if mols[0].GetNumConformers() == 0:
        log.error("SDF %s round-trips without coordinates", sdf_path)
        return False
    if expected_graph is not None:
        got = _flat_smiles(mols[0])
        if got != expected_graph:
            log.error("SDF %s round-trips as a different molecule: %s, expected %s",
                      sdf_path, got, expected_graph)
            return False
    return True


def load_single_mol(sdf_path: Path) -> Optional[Chem.Mol]:
    """First sanitised molecule in an SDF, or ``None``."""
    supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=False, sanitize=True)
    for mol in supplier:
        if mol is not None:
            return mol
    log.error("no valid molecule in %s", sdf_path)
    return None
