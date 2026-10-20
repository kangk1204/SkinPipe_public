"""Redock control: does this receptor and box reproduce the crystal pose?

Split out of :mod:`skinpipe.dock.vina` purely for file size; ``redock_validation``
is re-exported there.

The metric is *in-place* symmetry-corrected RMSD (``rdMolAlign.CalcRMS``) rather
than ``GetBestRMS``. ``GetBestRMS`` superimposes probe onto reference before
measuring, so it reports how well the conformer was reproduced and stays small
even when the pose landed in the wrong subpocket - useless as a control. Both
values are returned so the two failure modes (wrong conformer vs wrong place)
can be told apart: a small aligned RMSD with a large in-place RMSD means the
chemistry is fine and the placement is not.
"""
from __future__ import annotations

import logging
import math
from pathlib import Path
from typing import Any, Optional

from rdkit import Chem
from rdkit.Chem import rdMolAlign

from ..types import LigandRecord, Receptor
from . import poses as poseio
from .engine import DockingEngineError, dock_one

log = logging.getLogger(__name__)

#: Conventional redock success threshold, Angstrom.
REDOCK_RMSD_PASS: float = 2.0

#: The `failure_kind` values that mean *no RMSD was obtained*, as opposed to one
#: that was obtained and missed. They are not degrees of the same thing: a
#: measured miss fails this control at the declared search budget, while these
#: leave it unmeasured. Neither establishes the performance of other ligands.
#: Reporting them as "FAIL, in-place RMSD n/a A" told a reader the box had been
#: tested and found wrong when nothing had been compared at all.
NOT_MEASURED_KINDS: frozenset[str] = frozenset({
    "no_reference",   # no crystallographic ligand to redock, or it is unreadable
    "engine_error",   # the docking engine failed -- a transient, not a result
    "no_pose",        # the engine ran and returned nothing
    "unmatchable",    # poses exist, but no finite comparable RMSD was obtained
    "top_pose_missing",  # original engine pose 1 did not survive portable-pose checks
})
















REDOCK_N_POSES: int = 30


def redock_validation(receptor: Receptor, workdir: Path,
                      exhaustiveness: int = 16,
                      n_poses: int = REDOCK_N_POSES,
                      seed: int = 42, engine: str = "vina",
                      unidock_bin: Optional[str] = None,
                      unidock_search_mode: str = "balance") -> dict[str, Any]:
    """Re-dock the receptor's own cocrystal ligand and measure how far it lands.

    Returns ``rmsd`` (in-place, symmetry-corrected), ``rmsd_best_aligned``
    (``GetBestRMS``, i.e. conformer-only error after superposition), ``affinity``
    of the top pose, ``rmsd_all`` for every pose, ``best_pose_rmsd`` with the
    ``best_pose_index`` that achieved it (1-based), and ``pass``.

    The search budget and seed are exposed for repeat controls. A near-reference
    pose below the top rank demonstrates an ordering miss for this control;
    its presence does not validate other ligands or the preparation as a whole.
    Failure to sample one at a finite budget does not identify the cause.
    """
    workdir = Path(workdir)
    execution: dict[str, Any] = {"requested_engine": engine, "actual_engine": None,
                                "state": "not_started", "attempts": []}
    result: dict[str, Any] = {
        "pdb_id": receptor.pdb_id,
        "chain": receptor.chain,
        "ligand": receptor.cocrystal_ligand,
        "pocket_source": receptor.pocket_source,






        "reference_completeness": receptor.cocrystal_completeness,
        "reference_is_fragment": bool(receptor.cocrystal_is_fragment),
        "reference_atoms": receptor.cocrystal_atoms,
        "rmsd": None,
        "rmsd_best_aligned": None,
        "rmsd_all": [],
        # Declared in this function's docstring, so a caller that follows it
        # reads them unconditionally. They were only ever set on the success
        # path, which handed a `KeyError` to anyone inspecting a redock that
        # *failed* -- the one case the control exists to report on.
        "best_pose_rmsd": None,
        "best_pose_index": None,
        "affinity": None,
        "exhaustiveness": exhaustiveness,
        "seed": int(seed),
        "n_poses_requested": int(n_poses),
        "n_poses_returned": 0,
        "pose_indices": [],
        "top_pose_index": None,
        "first_readable_pose_index": None,
        "rmsd_pass_threshold_A": REDOCK_RMSD_PASS,
        "requested_engine": engine,
        "actual_engine": None,
        "execution": execution,
        "interpretation_scope": "this control ligand at the declared search settings",
        "pass": False,
        # The kind of failure, as a value rather than only as prose. Callers
        # have to act differently on each -- a search failure says raise
        # exhaustiveness, a scoring failure says the ordering itself is wrong
        # and no amount of searching fixes it -- and parsing that decision back
        # out of `note` is how a message edit becomes a behaviour change.
        #: None (control never invoked) | "pass" | "search" | "scoring"
        #: | any member of NOT_MEASURED_KINDS
        "failure_kind": None,
        "note": "",
    }

    if not receptor.cocrystal_sdf or not Path(receptor.cocrystal_sdf).is_file():
        # "No ligand" and "no usable ligand this run" are different facts and
        # only the first is a property of the entry. `place_box` writes the
        # reason into the receptor notes; surfacing it here keeps the reader
        # from concluding the structure is unusable when the ligand is in it
        # and a template fetch failed. 1HFC is the worked example: PLH, 25
        # atoms, present -- and reported as unavailable.
        transient = next(
            (n for n in (receptor.notes or [])
             if "could not be fetched from RCSB" in n), None)
        if receptor.cocrystal_ligand and transient:
            result["note"] = (
                f"redock validation could not run this time: {transient} "
                "The box and receptor preparation are unvalidated until it does."
            )
        elif receptor.cocrystal_ligand:
            result["note"] = (
                f"the crystallographic ligand {receptor.cocrystal_ligand} could "
                "not be turned into a comparable molecule, so redock validation "
                "cannot be run and the box and receptor preparation are "
                "unvalidated"
            )
        else:
            result["note"] = (
                "no crystallographic ligand in this entry - redock validation "
                "cannot be run, so the box and receptor preparation are "
                "unvalidated"
            )
        result["failure_kind"] = "no_reference"
        log.warning("%s: %s", receptor.pdb_id, result["note"])
        return result

    reference = poseio.load_single_mol(Path(receptor.cocrystal_sdf))
    if reference is None:
        result["note"] = "crystallographic ligand SDF could not be read"
        result["failure_kind"] = "no_reference"
        return result
    reference = Chem.RemoveHs(reference)

    smiles = Chem.MolToSmiles(reference)
    ligand = LigandRecord(
        ligand_id=f"redock_{receptor.pdb_id}_{receptor.cocrystal_ligand}",
        smiles=smiles,
        name=f"{receptor.cocrystal_ligand} (crystallographic)",
        source="pdb",
    )











    try:
        docked = dock_one(receptor, ligand, workdir,
                          exhaustiveness=exhaustiveness, n_poses=n_poses,
                          seed=seed, engine_name=engine, unidock_bin=unidock_bin,
                          unidock_search_mode=unidock_search_mode, execution=execution)
    except DockingEngineError as exc:
        result["failure_kind"] = "engine_error"
        result["note"] = (
            f"the docking engine failed on the control ligand "
            f"{receptor.cocrystal_ligand}: {exc}. Nothing was measured, so the "
            "box and receptor preparation are neither confirmed nor refuted; "
            "this is a fault on this run, not a property of the target.")
        log.warning("%s redock: %s", receptor.pdb_id, result["note"])
        return result
    finally:
        result["actual_engine"] = execution.get("actual_engine")
    result["n_poses_returned"] = len(docked)
    result["pose_indices"] = [pose.pose_index for pose in docked]
    if not docked:
        result["failure_kind"] = "no_pose"
        result["note"] = (
            f"redock of {receptor.cocrystal_ligand} produced no pose, so no "
            "RMSD was computed and the box is unvalidated rather than refuted")
        log.warning("%s redock: %s", receptor.pdb_id, result["note"])
        return result

    rmsds: list[Optional[float]] = []
    for pose in docked:
        probe = poseio.load_single_mol(Path(pose.sdf_path)) if pose.sdf_path else None
        rmsds.append(_symmetry_rmsd(probe, reference) if probe is not None else None)

    result["rmsd_all"] = rmsds
    result["first_readable_pose_index"] = docked[0].pose_index
    # SDF round-trip or operational-score checks can discard engine pose 1.
    # The first surviving pose is then a lower-ranked pose, not a successful
    # top-pose control. Bind this decision to the original engine ordinal.
    top_matches = [i for i, pose in enumerate(docked)
                   if type(pose.pose_index) is int and pose.pose_index == 1]
    top_at = top_matches[0] if len(top_matches) == 1 else None
    if top_at is not None:
        result["top_pose_index"] = 1
        result["affinity"] = docked[top_at].affinity
        result["rmsd"] = rmsds[top_at]
        result["rmsd_best_aligned"] = _aligned_rmsd(
            poseio.load_single_mol(Path(docked[top_at].sdf_path))
            if docked[top_at].sdf_path else None, reference,
        )
    valid = [r for r in rmsds if r is not None]
    result["best_pose_rmsd"] = min(valid) if valid else None
    # Two different numbers, and they used to be one. `best_at` is the position
    # in `docked`, which is what indexes back into it; the number *reported* is
    # the pose's own Vina index, which is what a reader looks up. `docked` holds
    # only the poses that survived the SDF round-trip, so the moment one is
    # dropped the positions slide and "pose 2 of 9" names a pose that is not the
    # one measured.
    best_at = rmsds.index(result["best_pose_rmsd"]) if valid else None
    result["best_pose_index"] = (docked[best_at].pose_index
                                 if best_at is not None else None)
    result["pass"] = bool(result["rmsd"] is not None
                          and result["rmsd"] <= REDOCK_RMSD_PASS)

    if top_at is None:
        result["failure_kind"] = "top_pose_missing"
        result["note"] = (
            "The original engine pose 1 is absent or not uniquely identified "
            "among readable saved poses. Its top-pose RMSD was not measured; "
            "a surviving lower-ranked pose cannot replace it. Sampled-best "
            "RMSD remains a separate observation, and the control is unvalidated."
        )
    elif result["rmsd"] is None:
        result["failure_kind"] = "unmatchable"
        result["note"] = (
            "RMSD could not be computed as a finite non-negative distance: "
            "the molecular graphs could not be matched or the native geometry "
            "result was unusable. The top-pose control remains unmeasured.")
    elif result["pass"]:
        result["failure_kind"] = "pass"
        result["note"] = (
            f"top pose {result['rmsd']:.2f} A from crystal for this control ligand "
            "at the declared search settings; other ligands and affinity ranking "
            "were not validated")
    elif result["best_pose_rmsd"] is not None and result["best_pose_rmsd"] <= REDOCK_RMSD_PASS:
        # The search found it; the scoring function threw it away.
        result["failure_kind"] = "scoring"
        gap = docked[best_at].affinity - docked[top_at].affinity
        result["note"] = (
            f"SCORING failure on this control: the top pose is {result['rmsd']:.2f} A out, "
            f"but pose {result['best_pose_index']} of {len(docked)} sits "
            f"{result['best_pose_rmsd']:.2f} A from the crystal and loses by only "
            f"{gap:.2f} kcal/mol. The sampled near-reference pose was not ranked "
            "first in this run. This does not validate receptor preparation or "
            "other ligands. Affinity ranking and any MM/GBSA refinement require "
            "independent validation.")
    else:
        result["failure_kind"] = "search"
        result["note"] = (
            f"SEARCH failure: no pose came within {REDOCK_RMSD_PASS} A "
            f"(best {result['best_pose_rmsd']:.2f} A of {len(docked)} poses) at "
            f"exhaustiveness {exhaustiveness}, seed {seed}. This is a miss at the "
            "declared search budget, not a diagnosis of the box or receptor. "
            "Repeat with predeclared seeds and appropriate search settings; "
            "inspect preparation and reference compatibility separately.")
    log.info("%s redock: %s", receptor.pdb_id, result["note"])
    return result


def _strip_for_comparison(mol: Chem.Mol) -> Chem.Mol:
    """Heavy-atom mapping representation with charge and H counts ignored.

    RMSD is a geometric question, but RDKit answers it through a substructure
    match, so any graph difference makes it fail outright. Since ligands are now
    titrated to skin pH before docking, the docked molecule routinely carries a
    charge the crystallographic reference does not -- 4-hydroxytamoxifen docks
    as the ammonium and the PDB ideal SDF is neutral. Without this the redock
    control silently returns "RMSD could not be computed" on exactly the
    ionisable ligands that most needed checking.
    This copy is only for atom correspondence. Sanitising after deleting a
    charged atom's H counts can assign spurious radicals: an ammonium with one
    heavy neighbour becomes a bare nitrogen radical and no longer matches the
    neutral reference. Ignore H counts on both molecules, keep their heavy
    bonds and coordinates, and update caches without inferring new chemistry.
    The actual prepared ligand and crystallographic reference remain intact.
    """
    work = Chem.RWMol(Chem.RemoveHs(mol))
    for atom in work.GetAtoms():
        atom.SetFormalCharge(0)
        atom.SetNumExplicitHs(0)
        atom.SetNoImplicit(True)
    out = work.GetMol()
    out.UpdatePropertyCache(strict=False)
    return out


def _finite_rmsd(value: float) -> Optional[float]:
    """A native numeric failure is unmeasured, not a search/scoring miss."""
    value = float(value)
    return value if math.isfinite(value) and value >= 0 else None


def _symmetry_rmsd(probe: Chem.Mol, reference: Chem.Mol) -> Optional[float]:
    """In-place symmetry-corrected RMSD; no superposition."""
    probe, reference = Chem.RemoveHs(probe), Chem.RemoveHs(reference)
    if (probe.GetNumAtoms(), probe.GetNumBonds()) != (reference.GetNumAtoms(), reference.GetNumBonds()):
        return None
    try:
        return _finite_rmsd(rdMolAlign.CalcRMS(probe, reference))
    except (RuntimeError, ValueError):
        pass
    try:      # graphs differ, most often because the docked form was titrated
        return _finite_rmsd(rdMolAlign.CalcRMS(_strip_for_comparison(probe),
                                              _strip_for_comparison(reference)))
    except (RuntimeError, ValueError) as exc:
        log.warning("in-place RMSD failed even on the neutral skeleton: %s", exc)
        return None


def _aligned_rmsd(probe: Optional[Chem.Mol], reference: Chem.Mol) -> Optional[float]:
    """``GetBestRMS``: superimposes first, so it reports conformer error only."""
    if probe is None:
        return None
    probe, reference = Chem.RemoveHs(probe), Chem.RemoveHs(reference)
    if (probe.GetNumAtoms(), probe.GetNumBonds()) != (reference.GetNumAtoms(), reference.GetNumBonds()):
        return None
    try:
        return _finite_rmsd(rdMolAlign.GetBestRMS(Chem.Mol(probe), Chem.Mol(reference)))
    except (RuntimeError, ValueError):
        pass
    try:
        return _finite_rmsd(rdMolAlign.GetBestRMS(_strip_for_comparison(probe),
                                                 _strip_for_comparison(reference)))
    except (RuntimeError, ValueError) as exc:
        log.warning("aligned RMSD failed even on the neutral skeleton: %s", exc)
        return None
