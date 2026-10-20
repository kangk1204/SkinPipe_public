"""Receptor preparation and OpenMM system assembly for MM/GBSA.

Force-group convention used throughout the module (see
:func:`build_complex_system`)::

    group 0  bonded + NonbondedForce      -> gas-phase MM energy
    group 1  CustomGBForce                -> GB polar (+ ACE nonpolar)
    group 2  positional restraints        -> excluded from every reported energy

Evaluating energies by force group means the restraint bias never leaks into a
reported number and no Context has to be rebuilt between minimisation and
scoring.
"""
from __future__ import annotations

import logging
import math
import random
from dataclasses import dataclass, field
from io import StringIO
from pathlib import Path
from threading import RLock
from typing import Any, Optional

from . import metal_guard

LOGGER = logging.getLogger(__name__)

MM_GROUP = 0
GB_GROUP = 1
RESTRAINT_GROUP = 2

#: ``implicit_solvent`` selector -> OpenMM force-field XML.
IMPLICIT_SOLVENT_FILES: dict[str, str] = {
    "gbn2": "implicit/gbn2.xml",
    "obc2": "implicit/obc2.xml",
    "obc1": "implicit/obc1.xml",
    "hct": "implicit/hct.xml",
}

DEFAULT_PROTEIN_FF = "amber14-all.xml"
DEFAULT_PREPARATION_SEED = 20260929
PREPARATION_PLATFORM = "Reference"
# PDBFixer/OpenMM's hydrogen placement uses Python's process-global RNG.
# Serialize cooperating MM/GBSA preparations while saving/restoring that state.
_PREPARATION_LOCK = RLock()

#: ACE nonpolar surface tension baked into OpenMM's implicit XMLs: the energy
#: term carries 28.3919551 kJ/mol/nm^2 == 4*pi*gamma with
#: gamma = 0.0054 kcal/mol/A^2 (AMBER default ``SURFTEN``; Schaefer, Bartels &
#: Karplus, *J. Mol. Biol.* 1998, 284, 835).
ACE_SURFACE_TENSION_KCAL_PER_A2 = 0.0054

KJ_PER_KCAL = 4.184

#: Amber ion-exclusion factor (0.73) times the A^-1 -> nm^-1 conversion (10).
SALT_KAPPA_CONVERSION = 7.3


class SystemBuildError(RuntimeError):
    """Raised when a receptor/ligand complex cannot be turned into a System."""


def validate_preparation_seed(seed: int) -> int:
    """Require an explicit, nonzero OpenMM random seed (zero means random)."""
    if isinstance(seed, bool) or not isinstance(seed, int) or not 0 < seed < 2**31:
        raise ValueError("preparation_seed must be an integer from 1 to 2147483647")
    return seed


def kappa_from_salt(salt_molar: float, temperature_k: float = 300.0,
                    solvent_dielectric: float = 78.5) -> float:
    """Debye screening parameter in nm^-1 for a 1:1 salt.

    Reproduces exactly what OpenMM does with ``implicitSolventSaltConc`` for
    Amber prmtop files (``openmm/app/amberprmtopfile.py``)::

        kappa[1/A]  = 50.33355 * sqrt(C / eps / T)
        kappa[1/nm] = kappa[1/A] * 0.73 * 10

    ``C`` is mol/L, ``T`` kelvin.  The 0.73 is Amber's ion-exclusion
    correction; the 10 converts A^-1 to OpenMM's nm^-1.  Physiological 0.15 M
    at 300 K gives ~0.927 nm^-1 (Debye length ~10.8 A), and the screened GB
    term is ``exp(-kappa*B)`` inside the CustomGBForce.

    Dropping the 7.3 leaves the screening ten times too weak -- silently, since
    the energies stay finite and plausible.
    """
    for name, value, zero_allowed in (("salt_molar", salt_molar, True),
                                     ("temperature_k", temperature_k, False),
                                     ("solvent_dielectric", solvent_dielectric, False)):
        if (isinstance(value, bool) or not isinstance(value, (int, float)) or
                not math.isfinite(value) or (value < 0 if zero_allowed else value <= 0)):
            raise ValueError(f"{name} must be finite and {'>= 0' if zero_allowed else '> 0'}")
    if salt_molar == 0:
        return 0.0
    kappa_per_angstrom = 50.33355 * math.sqrt(
        salt_molar / solvent_dielectric / temperature_k)
    return kappa_per_angstrom * SALT_KAPPA_CONVERSION


def prepare_receptor(receptor_pdb: Path, out_pdb: Path, ph: float = 7.4,
                     keep_water: bool = False, *,
                     preparation_seed: int = DEFAULT_PREPARATION_SEED) -> Path:
    """Retain supported native polymer chemistry and prepare the receptor.

    Amber14 has no templates for waters-with-no-hydrogens, metals, cofactors or
    the co-crystallised ligand, so anything that is not a standard amino acid
    is removed, except supported native phosphorylated amino acids. Their
    original heavy graph and both peptide links are required, and Amber14
    with the installed phosaa14SB supplement represents the template-defined
    protonation state. Unsupported modified chemistry raises before output.
    Missing side-chain atoms are rebuilt because a truncated
    residue in the binding site silently changes the interaction energy.

    Note that ``removeHeterogens`` also strips ACE/NME capping groups, since
    PDBFixer does not count them as standard residues.  If the receptor's
    biology depends on a cofactor, a structural metal or a capped terminus,
    this function is the wrong preparation step -- prepare the receptor
    yourself and pass ``prepare_receptor_pdb=False``.

    Catalytic metals are the one case this function does not simply forget.
    Whatever it deletes (and whatever was already held out upstream, recorded
    in a sidecar by :mod:`skinpipe.structure.metal_sites`) is carried forward
    into a sidecar beside ``out_pdb``, so the prepared receptor still knows its
    site has lost a metal and :mod:`skinpipe.mmgbsa.metal_guard` can refuse to
    score a pose that binds one.  Without that the loss is invisible: a PDB
    with no copper in it looks exactly like a protein that never had one.

    The empty case is handled too, and is not a no-op: preparing a metal-free
    receptor to a path some holo one used before removes that path's sidecar,
    so the record cannot outlive the preparation that wrote it.

    Both stochastic repair steps use ``preparation_seed`` on the Reference
    platform. Python's RNG state is restored even when repair fails. This
    pins preparation within the same dependency environment; it does not
    promise identical coordinates across OpenMM/PDBFixer version changes.
    """
    from openmm import Platform
    from openmm.app import PDBFile
    from pdbfixer import PDBFixer

    from ..structure import metal_sites

    preparation_seed = validate_preparation_seed(preparation_seed)
    receptor_pdb = Path(receptor_pdb)
    if not receptor_pdb.is_file():
        raise SystemBuildError(f"receptor PDB not found: {receptor_pdb}")
    from ..structure import clean, modified
    input_text = receptor_pdb.read_text()
    native = modified.groups(input_text, known_phospho=True)
    try:
        model = modified.parameterization_descriptor(input_text)
    except (ImportError, ValueError) as exc:
        raise SystemBuildError(f"native modified receptor parameters unavailable: {exc}; output not written") from exc
    out_pdb = Path(out_pdb)
    carried = metal_sites.all_metal_lines(receptor_pdb)

    try:
        platform = Platform.getPlatformByName(PREPARATION_PLATFORM)
    except Exception as exc:
        raise SystemBuildError(
            f"deterministic receptor preparation requires the "
            f"{PREPARATION_PLATFORM} OpenMM platform: {exc}") from exc
    fixer = PDBFixer(filename=str(receptor_pdb), platform=platform)
    n_in = fixer.topology.getNumAtoms()
    if native:
        try:
            modified.register_templates(fixer, native)
        except ValueError as exc:
            raise SystemBuildError(f"native modified polymer connectivity unsupported: {exc}; output not written") from exc
        modified.remove_other_heterogens(fixer, native, keep_water=keep_water)
    else:
        fixer.removeHeterogens(keepWater=keep_water)
    fixer.findMissingResidues()
    # Do not model in entire absent loops: a guessed loop far from the pocket
    # adds noise to the total energy without improving the interaction term.
    fixer.missingResidues = {}
    fixer.findNonstandardResidues()
    if native:
        fixer.nonstandardResidues = [item for item in fixer.nonstandardResidues
                                    if (item[0].chain.id, item[0].id + item[0].insertionCode.strip(),
                                        item[0].name) not in native]
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    with _PREPARATION_LOCK:
        rng_state = random.getstate()
        random.seed(preparation_seed)
        try:
            fixer.addMissingAtoms(seed=preparation_seed)
            if native:
                modified.deduplicate_repair_bonds(fixer)
                modified.validate_native_bonds(fixer.topology, native)
                modified.add_hydrogens(fixer, native, ph)
            else:
                fixer.addMissingHydrogens(ph)
        finally:
            random.setstate(rng_state)

    n_out = fixer.topology.getNumAtoms()
    if n_out == 0:
        raise SystemBuildError(
            f"{receptor_pdb.name} has no standard protein residues left after "
            "removing heterogens")
    # OpenMM's header includes today's date. It changes file identity across
    # days despite unchanged coordinates; keep its version, omit only the date.
    buffer = StringIO()
    PDBFile.writeFile(fixer.topology, fixer.positions, buffer, keepIds=True)
    lines = buffer.getvalue().splitlines(keepends=True)
    if lines and lines[0].startswith("REMARK   1 CREATED WITH OPENMM "):
        lines[0] = f"REMARK   1 CREATED WITH OPENMM {Platform.getOpenMMVersion()}\n"
    text = "".join(lines)
    if native:
        text = "\n".join(modified.declaration_lines(input_text)) + "\n" + text
        observed = {key for key in clean._heavy_atom_keys(input_text) if key[:3] in native}
        clean._validate_observed_heavy_atoms(text, input_text, observed)
        clean._validate_assembled_receptor_records(clean.atom_lines(text))
        modified.validate_native_bonds(fixer.topology, native)
        if modified.parameterization_descriptor(text) != model:
            raise SystemBuildError("native modified receptor model changed during preparation; output not written")
    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    out_pdb.write_text(text)
    LOGGER.info("prepared receptor %s -> %s: %d atoms in, %d out, %d residues",
                receptor_pdb.name, out_pdb.name, n_in, n_out,
                fixer.topology.getNumResidues())
    if carried:
        sidecar = metal_sites.write_sidecar(out_pdb, carried)
        LOGGER.warning("%s has %d catalytic metal atom(s) that amber14/GBn2 cannot "
                       "represent and this preparation does not keep; recorded in %s",
                       receptor_pdb.name, len(carried), sidecar)
    elif metal_sites.clear_sidecar(out_pdb):
        # Not "nothing to do".  ``out_pdb`` is a derived path and callers reuse
        # it: mmgbsa_batch prepares into ``workdir/<stem>_prepared.pdb``, and a
        # workdir is meant to be re-run -- that is how the resume cache works.
        # A sidecar left by whatever was prepared here last describes a
        # structure this file is not, and
        # metal_guard believes it: it would refuse this metal-free receptor, or
        # stamp an apo-site warning on a number that lost nothing.  The record's
        # lifetime belongs to the preparation that wrote it, so an empty metal
        # set has to clear it rather than leave it standing.
        #
        # After the write, not before: if PDBFixer raises, ``out_pdb`` still
        # holds the previous preparation, and its sidecar still describes it.
        LOGGER.warning("removed a stale metal record beside %s: %s has no "
                       "catalytic metals, and the old one would have made "
                       "skinpipe.mmgbsa.metal_guard refuse this receptor",
                       out_pdb.name, receptor_pdb.name)
    return out_pdb


@dataclass
class SpeciesSystems:
    """The two Systems needed to decompose one species' energy.

    ``full`` carries the ACE nonpolar term, ``polar_only`` does not; the
    difference of their GB force groups is the nonpolar contribution.
    """

    topology: Any
    full: Any          # openmm.System, SA='ACE'
    polar_only: Any    # openmm.System, SA=None
    n_atoms: int


@dataclass
class ComplexSystem:
    """Receptor + ligand complex and the two subsystems carved out of it."""

    complex: SpeciesSystems
    receptor: SpeciesSystems
    ligand: SpeciesSystems
    positions: Any                       # openmm.unit.Quantity, complex order
    receptor_slice: slice
    ligand_slice: slice
    restraint_force: Optional[Any] = None
    implicit_solvent: str = "gbn2"
    salt_molar: float = 0.15
    kappa: float = 0.0
    provenance: list[str] = field(default_factory=list)
    parameterization: Optional[dict[str, Any]] = None

    @property
    def n_atoms(self) -> int:
        return self.complex.n_atoms


def _make_forcefield(implicit_solvent: str, generator: Any,
                     protein_ff: str = DEFAULT_PROTEIN_FF, *,
                     modified_files: tuple[Any, ...] = ()) -> Any:
    from openmm.app import ForceField

    key = (implicit_solvent or "gbn2").strip().lower()
    if key not in IMPLICIT_SOLVENT_FILES:
        raise SystemBuildError(
            f"unknown implicit solvent {implicit_solvent!r}; "
            f"choose from {sorted(IMPLICIT_SOLVENT_FILES)}")
    ff = ForceField(protein_ff, *modified_files, IMPLICIT_SOLVENT_FILES[key])
    ff.registerTemplateGenerator(generator.generator)
    return ff


def _tag_force_groups(system: Any) -> None:
    """Put the GB force in :data:`GB_GROUP` and everything else in ``MM_GROUP``."""
    import openmm

    for force in system.getForces():
        if isinstance(force, openmm.CustomGBForce):
            force.setForceGroup(GB_GROUP)
        else:
            force.setForceGroup(MM_GROUP)


def _create_pair(ff: Any, topology: Any, kappa: float, constrain_h: bool,
                 cutoff_nm: Optional[float] = None) -> SpeciesSystems:
    from openmm import app, unit

    common: dict[str, Any] = {
        "constraints": app.HBonds if constrain_h else None,
        "rigidWater": False,
        "removeCMMotion": False,
        "implicitSolventKappa": kappa,
    }
    if cutoff_nm:
        common["nonbondedMethod"] = app.CutoffNonPeriodic
        common["nonbondedCutoff"] = cutoff_nm * unit.nanometer
    else:
        common["nonbondedMethod"] = app.NoCutoff
    try:
        full = ff.createSystem(topology, sasaMethod="ACE", **common)
        polar = ff.createSystem(topology, sasaMethod=None, **common)
    except Exception as exc:  # noqa: BLE001
        raise SystemBuildError(f"createSystem failed: {exc}") from exc
    _tag_force_groups(full)
    _tag_force_groups(polar)
    return SpeciesSystems(topology=topology, full=full, polar_only=polar,
                          n_atoms=topology.getNumAtoms())


def _add_restraints(system: Any, topology: Any, positions: Any,
                    k_kcal_per_mol_a2: float, n_receptor_atoms: int) -> Optional[Any]:
    """Harmonic positional restraints on protein heavy atoms.

    Without these, 500 steps of minimisation on an implicit-solvent protein
    pulls the backbone off the crystal geometry and the receptor internal
    energy swamps the interaction term.  The ligand is left completely free so
    it can relax into the pocket.
    """
    import openmm
    from openmm import unit

    if k_kcal_per_mol_a2 <= 0:
        return None
    force = openmm.CustomExternalForce("k_rest*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    force.addGlobalParameter(
        "k_rest",
        k_kcal_per_mol_a2 * KJ_PER_KCAL * 100.0 * unit.kilojoule_per_mole / unit.nanometer**2,
    )
    for name in ("x0", "y0", "z0"):
        force.addPerParticleParameter(name)

    atoms = list(topology.atoms())
    coords = positions.value_in_unit(unit.nanometer)
    n_restrained = 0
    for atom in atoms:
        if atom.index >= n_receptor_atoms or atom.element is None:
            continue
        if atom.element.symbol == "H":
            continue
        force.addParticle(atom.index, list(coords[atom.index]))
        n_restrained += 1
    if n_restrained == 0:
        return None
    force.setForceGroup(RESTRAINT_GROUP)
    system.addForce(force)
    LOGGER.debug("restrained %d protein heavy atoms at %.1f kcal/mol/A^2",
                 n_restrained, k_kcal_per_mol_a2)
    return force


def build_complex_system(receptor_pdb: Path, ligand_sdf: Path,
                         implicit_solvent: str = "GBn2",
                         ligand_forcefield: str = "gaff-2.11",
                         charge_method: str = "auto",
                         salt_molar: float = 0.15,
                         temperature_k: float = 300.0,
                         restraint_k: float = 5.0,
                         constrain_h_bonds: bool = False,
                         nonbonded_cutoff_nm: Optional[float] = None,
                         prepare: bool = True,
                         on_metal_site: str = metal_guard.REFUSE,
                         workdir: Optional[Path] = None,
                         preparation_seed: int = DEFAULT_PREPARATION_SEED) -> ComplexSystem:
    """Assemble the complex and the receptor/ligand subsystems.

    The receptor and ligand Systems are built on *subsets of the complex
    topology* and are later evaluated on *complex coordinates*.  This is the
    single-trajectory MM/GBSA approximation -- see
    :func:`skinpipe.mmgbsa.openmm_gbsa.mmgbsa_single` for why separating the
    geometries instead is wrong.

    ``nonbonded_cutoff_nm`` trades accuracy for wall time.  ``None`` (default)
    is exact: every pair interacts, which is what MM/GBSA assumes, but it is
    O(N^2) and on a CPU that dominates the run.  A value around 1.6-2.0 nm
    makes a ~300-residue complex several times faster; because the *same*
    cutoff is applied to the complex, the receptor and the ligand, most of the
    truncation error cancels in the difference, but long-range electrostatics
    between distant charged residues does not.  Benchmark it on your own
    system before trusting it for ranking.

    ``constrain_h_bonds`` defaults to False.  Constraints only matter for MD,
    and they make OpenMM's minimiser an order of magnitude slower (see
    :func:`skinpipe.mmgbsa.engine.minimise`).  Potential energies are
    unaffected either way for the reported dG: ``HBonds`` deletes the same
    X-H bond terms from the complex and from both subsystems, so they cancel.

    ``on_metal_site`` decides what happens when the pose sits on a catalytic
    metal this force field cannot represent -- the tyrosinase and MMP case.
    The default :data:`~skinpipe.mmgbsa.metal_guard.REFUSE` raises
    :class:`~skinpipe.mmgbsa.metal_guard.MetalSiteError` before any work is
    done; ``"score_without_metal"`` computes the apo number and stamps every
    result with the fact.  See :mod:`skinpipe.mmgbsa.metal_guard`.
    """
    from openmm.app import Modeller, PDBFile

    from .ligand import parameterise_ligand

    solvent_key = (implicit_solvent or "gbn2").strip().lower()
    if solvent_key not in IMPLICIT_SOLVENT_FILES:
        raise SystemBuildError(f"unknown implicit solvent {implicit_solvent!r}; "
                               f"expected one of {', '.join(IMPLICIT_SOLVENT_FILES)}")
    receptor_pdb = Path(receptor_pdb)
    ligand_sdf = Path(ligand_sdf)
    preparation_seed = validate_preparation_seed(preparation_seed)
    # First, before charge assignment or anything else expensive: a refusal
    # should cost milliseconds, not the minutes it takes to reach the point
    # where the metal's absence would have mattered.
    metal_notes = metal_guard.check_binding_site(
        receptor_pdb, metal_guard.ligand_coords(ligand_sdf),
        policy=on_metal_site, context=ligand_sdf.name)
    from ..structure import modified
    try:
        native_model = modified.parameterization_descriptor(receptor_pdb.read_text())
        if native_model:
            native_input = PDBFile(str(receptor_pdb))
            modified.validate_native_bonds(native_input.topology,
                                          modified.groups(receptor_pdb.read_text(), known_phospho=True))
    except (ImportError, ValueError) as exc:
        raise SystemBuildError(f"native modified receptor cannot be parameterized: {exc}") from exc
    charge_workdir = Path(workdir) / "charge_producers" if workdir else None
    params = parameterise_ligand(ligand_sdf, forcefield=ligand_forcefield,
                                 charge_method=charge_method,
                                 charge_workdir=charge_workdir)

    if prepare:
        workdir = Path(workdir) if workdir else ligand_sdf.parent
        receptor_pdb = prepare_receptor(
            receptor_pdb, workdir / f"{receptor_pdb.stem}_prepared.pdb",
            preparation_seed=preparation_seed)
    pdb = PDBFile(str(receptor_pdb))
    n_receptor = pdb.topology.getNumAtoms()
    if n_receptor == 0:
        raise SystemBuildError(f"{receptor_pdb.name} has no atoms after preparation")

    lig_topology = params.molecule.to_topology().to_openmm()
    lig_positions = params.molecule.conformers[0].to_openmm()

    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.add(lig_topology, lig_positions)
    n_total = modeller.topology.getNumAtoms()
    n_ligand = n_total - n_receptor
    if n_ligand != params.n_atoms:
        raise SystemBuildError(
            f"ligand atom count mismatch: {n_ligand} in complex vs {params.n_atoms}")

    kappa = kappa_from_salt(salt_molar, temperature_k)
    modified_files = modified.amber_supplements(
        {row['name'] for row in native_model['native_residues']}) if native_model else ()
    ff = _make_forcefield(implicit_solvent, params.generator, modified_files=modified_files)

    cutoff = nonbonded_cutoff_nm
    cplx = _create_pair(ff, modeller.topology, kappa, constrain_h_bonds, cutoff)
    receptor = _create_pair(ff, pdb.topology, kappa, constrain_h_bonds, cutoff)
    ligand = _create_pair(ff, lig_topology, kappa, constrain_h_bonds, cutoff)

    restraint = _add_restraints(cplx.full, modeller.topology, modeller.positions,
                                restraint_k, n_receptor)
    parameterization = None
    if params.charge_producer is not None:
        parameterization = _fixed_parameterization_receipt(
            params.charge_producer, cplx, receptor, ligand, n_receptor,
            modeller.positions, Path(workdir) if workdir else None)
    if native_model:
        if modified.parameterization_descriptor(receptor_pdb.read_text()) != native_model:
            raise SystemBuildError("native modified receptor model changed during System construction")
        parameterization = {**(parameterization or {}), "receptor_modified_model": native_model}

    provenance = [params.provenance(),
                  f"solvent={implicit_solvent.lower()}",
                  f"salt={salt_molar:.3f}M(kappa={kappa:.4f}/nm)",
                  f"restraint_k={restraint_k:g} kcal/mol/A^2",
                  f"cutoff={cutoff:g}nm" if cutoff else "cutoff=none",
                  *metal_notes]
    if prepare:
        provenance.append(f"receptor_preparation_seed={preparation_seed}; "
                          f"receptor_preparation_platform={PREPARATION_PLATFORM}")
    if native_model:
        provenance.append("receptor_forcefield=amber14-all.xml+phosaa14SB.xml; "
                          "native modified protonation is template-defined, not a site-specific pKa estimate")
    return ComplexSystem(
        complex=cplx, receptor=receptor, ligand=ligand,
        positions=modeller.positions,
        receptor_slice=slice(0, n_receptor),
        ligand_slice=slice(n_receptor, n_total),
        restraint_force=restraint,
        implicit_solvent=implicit_solvent.lower(), salt_molar=salt_molar,
        kappa=kappa, provenance=provenance,
        parameterization=parameterization,
    )


def _fixed_parameterization_receipt(charge_receipt: dict[str, Any],
                                    complex_systems: SpeciesSystems,
                                    receptor_systems: SpeciesSystems,
                                    ligand_systems: SpeciesSystems,
                                    n_receptor: int, positions: Any,
                                    workdir: Optional[Path]) -> dict[str, Any]:
    """Bind the actual charge vector to all six Systems before minimization."""
    import gzip
    import hashlib
    import tempfile

    import numpy as np
    import openmm
    from openmm import unit

    from .fixed_amber import vector_sha256

    expected = charge_receipt["normalized_charges_e"]
    rows = []
    directory = None
    if workdir is not None:
        base = workdir / "parameterization_systems"
        base.mkdir(parents=True, exist_ok=True)
        directory = Path(tempfile.mkdtemp(prefix="fixed-system-", dir=base))
    for species, pair in (("complex", complex_systems), ("receptor", receptor_systems),
                           ("ligand", ligand_systems)):
        for name, system in (("full", pair.full), ("polar_only", pair.polar_only)):
            forces = [force for force in system.getForces() if isinstance(force, openmm.NonbondedForce)]
            if len(forces) != 1:
                raise SystemBuildError("cannot verify the fixed charge vector in this System")
            force = forces[0]
            charges = [float(force.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge))
                       for i in range(force.getNumParticles())]
            observed = charges[n_receptor:] if species == "complex" else charges
            if species != "receptor" and observed != expected:
                raise SystemBuildError(f"{species}/{name} did not preserve the produced ligand charges")
            xml = openmm.XmlSerializer.serialize(system).encode()
            row = {"species": species, "kind": name, "particles": system.getNumParticles(),
                   "xml_bytes": len(xml), "xml_sha256": hashlib.sha256(xml).hexdigest(),
                   "nonbonded_charge_vector_sha256": vector_sha256(charges)}
            if directory is not None:
                path = directory / (species + "_" + name + ".xml.gz")
                with path.open("wb") as stream:
                    with gzip.GzipFile(filename="", mode="wb", fileobj=stream, mtime=0) as compressed:
                        compressed.write(xml)
                row["saved_path"] = str(path)
                row["saved_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
            rows.append(row)
    coordinates = np.asarray(positions.value_in_unit(unit.nanometer), dtype="<f8")
    return {"charges": charge_receipt, "systems": rows, "ligand_system_charges_verified": True,
            "initial_positions_nm_float64_sha256": hashlib.sha256(coordinates.tobytes(order="C")).hexdigest(),
            "scope": "actual pre-minimization parameterization; no affinity or numerical stability claim"}
