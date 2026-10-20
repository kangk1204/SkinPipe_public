"""Ligand loading and force-field parameterisation for MM/GBSA.

Small molecules are read through the OpenFF ``Molecule`` class rather than
through a bare PDB/mol2 reader, because that is the only route that preserves
bond orders, aromaticity and formal charges -- all three are required to pick
correct atom types and to get the net charge of the ligand right.

Two template-generator backends are supported:

``gaff-*``
    :class:`openmmforcefields.generators.GAFFTemplateGenerator`.  Requires
    AmberTools (``antechamber``/``parmchk2``) on ``PATH``.
``openff-*`` / ``smirnoff``
    :class:`openmmforcefields.generators.SMIRNOFFTemplateGenerator`.  Pure
    Python, no external binaries.

If GAFF is requested but AmberTools is missing, the module falls back to
``SMIRNOFF_FALLBACK`` and records the substitution, rather than failing the
whole run.  The substitution is always logged at WARNING level and surfaced in
the ``notes`` of the resulting ``GBSAResult``.

Requested is not resolved
-------------------------
Nothing here is guaranteed to use what it was asked for.  ``gaff-2.11``
becomes ``openff-2.2.0`` on a machine without AmberTools; ``charge_method=
"auto"`` is whichever of NAGL, AM1BCC and Gasteiger answers first, and an
explicit request falls through the same ladder.  Two runs that asked for the
same thing therefore routinely compute two different numbers.

So the *resolved* path has to be legible to the caller, in two forms:

* after the fact, in :meth:`LigandParameters.provenance` -- and
  :func:`parse_provenance` reads it back out of a finished record's ``notes``,
  which is how :mod:`skinpipe.mmgbsa.cachekey` learns what actually ran;
* before the fact, from :func:`resolve_forcefield`,
  :func:`resolve_charge_method` and :func:`backend_versions`, so the MM/GBSA
  cache key can refuse a record computed down a different path instead of
  serving Gasteiger charges to a run that resolves AM1BCC (N-12).

The prediction and the run walk **one** ladder, :func:`charge_ladder`, rather
than two copies of it; the availability probes are the same functions the run
itself calls.  That is deliberate: the cache key has fallen behind the worker
before, by copying what the worker did instead of asking it.
"""
from __future__ import annotations

import importlib.metadata
import importlib.util
import json
import logging
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from contextlib import ExitStack
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

from .amber_identity import ambertools_identity

LOGGER = logging.getLogger(__name__)

#: Used when a GAFF force field is requested but AmberTools is unavailable.
SMIRNOFF_FALLBACK = "openff-2.2.0"

#: Graph-neural-network reproduction of AM1BCC charges.  ~0.5 s per molecule
#: versus minutes for a real AM1 optimisation, and no AmberTools dependency.
#: Ref: Horton et al., *J. Chem. Inf. Model.* 2022, 62, 5539 (OpenFF NAGL).
NAGL_MODELS = ("openff-gnn-am1bcc-1.0.0.pt", "openff-gnn-am1bcc-0.1.0-rc.3.pt")

#: Preference order used when ``charge_method="auto"``.
CHARGE_LADDER = ("nagl", "am1bcc", "gasteiger")

#: Charge methods that shell out to an external binary and can therefore hang.
_EXTERNAL_CHARGE_METHODS = frozenset({"am1bcc", "am1-mulliken", "am1bccelf10"})
FIXED_AMBER_METHOD = "am1bcc-fixed"


def validate_screen_protocol(charge_method: str, forcefield: str) -> tuple[str, str]:
    """Validate screen-level settings before fetching or creating output."""
    if not isinstance(charge_method, str) or not isinstance(forcefield, str):
        raise TypeError("MM/GBSA charge method and ligand force field must be strings")
    method, ff = charge_method.strip().lower(), forcefield.strip().lower()
    if method not in {"auto", FIXED_AMBER_METHOD, *CHARGE_LADDER, *_EXTERNAL_CHARGE_METHODS}:
        raise ValueError("unsupported MM/GBSA charge method: " + repr(charge_method))
    if re.fullmatch(r"(?:gaff|openff)-\d+(?:\.\d+){1,2}", ff) is None:
        raise ValueError("MM/GBSA ligand force field must name a GAFF or OpenFF version")
    return method, ff

_NET_CHARGE_TOL = 1.0e-4

#: Prefixes :meth:`LigandParameters.provenance` writes and
#: :func:`parse_provenance` reads.  They are constants because a record's
#: ``notes`` is the only place the *resolved* path survives, and the reader
#: must not drift from the writer.
FORCEFIELD_NOTE = "ff="
CHARGE_NOTE = "charges="

#: Which distributions actually produce the numbers, per resolved path.  A
#: version bump in one of these changes the energy; a bump anywhere else does
#: not, so only the ones on the path are recorded -- see
#: :func:`backend_versions`.
_ALWAYS_BACKENDS = ("openff-toolkit", "openmmforcefields")
_CHARGE_BACKENDS: dict[str, tuple[str, ...]] = {
    "nagl": ("openff-nagl", "openff-nagl-models"),
    "am1bcc": ("ambertools",),
    FIXED_AMBER_METHOD: ("ambertools",),
    "am1-mulliken": ("ambertools",),
    "am1bccelf10": ("ambertools",),
    "gasteiger": ("rdkit",),
}

_ABSENT = "absent"


class LigandParameterisationError(RuntimeError):
    """Raised when a ligand cannot be turned into OpenMM parameters."""


@dataclass
class LigandParameters:
    """A parameterised ligand plus provenance for the ``notes`` field."""

    molecule: Any                      # openff.toolkit.Molecule
    generator: Any                     # GAFF/SMIRNOFF template generator
    forcefield: str                    # force field actually used
    charge_method: str                 # charge method actually used
    net_charge: int
    n_atoms: int
    warnings: list[str] = field(default_factory=list)
    charge_producer: Optional[dict[str, Any]] = None
    charge_model_sha256: Optional[str] = None

    def provenance(self) -> str:
        'What actually ran, and every substitution made to get there.'
        bits = [f"{FORCEFIELD_NOTE}{self.forcefield}",
                f"{CHARGE_NOTE}{self.charge_method}",
                f"q={self.net_charge:+d}"]
        if self.charge_model_sha256:
            bits.append(f"charge_model_sha256={self.charge_model_sha256}")
        return "; ".join(bits + self.warnings)


def parse_provenance(notes: str) -> dict[str, str]:
    """The resolved path read back out of a finished record's ``notes``.

    Returns ``{}`` when there is nothing to read -- a pose that failed before
    parameterisation has no force field and no charges, and inventing one for
    it would be worse than admitting it is unknown.

    >>> parse_provenance("ff=openff-2.2.0; charges=nagl:model.pt; q=+0")
    {'ligand_forcefield': 'openff-2.2.0', 'charge_method': 'nagl', \
'charge_model': 'nagl:model.pt'}
    """
    found: dict[str, str] = {}
    for token in str(notes or "").split(";"):
        token = token.strip()
        if token.startswith("charge_model_sha256="):
            found["charge_model_sha256"] = token.split("=", 1)[1]
        if token.startswith(FORCEFIELD_NOTE):
            found["ligand_forcefield"] = token[len(FORCEFIELD_NOTE):].strip()
        elif token.startswith(CHARGE_NOTE):
            used = token[len(CHARGE_NOTE):].strip()
            found["charge_method"] = charge_method_family(used)
            found["charge_model"] = used
    return found


def charge_method_family(used: str) -> str:
    """``"nagl:openff-gnn-am1bcc-1.0.0.pt"`` -> ``"nagl"``.

    The family is what can be predicted before a molecule is charged; the part
    after the colon names the model file and is only knowable afterwards.
    """
    return str(used or "").split(":", 1)[0].strip().lower()


_AMBER_EXES = ("antechamber", "parmchk2")


def ensure_amber_tools_on_path() -> bool:
    """Make AmberTools reachable and say whether it is. Idempotent.

    A conda or micromamba environment invoked by absolute path -- which is how
    this pipeline is normally run, and how its own README says to run it -- has
    its ``bin`` directory populated but **not** exported on ``PATH``. A bare
    ``shutil.which`` therefore misses ``antechamber`` and ``parmchk2`` sitting
    right next to ``python``, and every MM/GBSA run silently swaps the requested
    GAFF force field for OpenFF while reporting in ``notes`` that AmberTools
    "is not installed". :mod:`skinpipe.dock.protonate` already solves this for
    Open Babel by looking beside ``sys.executable``.

    Looking is not sufficient here: ``GAFFTemplateGenerator`` shells out to
    those binaries itself and resolves them through ``PATH``, so finding them
    without exporting them would turn a silent substitution into a hard
    parameterisation failure. The interpreter's own ``bin`` directory is
    therefore prepended to ``PATH`` -- which is exactly what activating the
    environment would have done.
    """
    if all(shutil.which(exe) for exe in _AMBER_EXES):
        return True
    env_bin = os.path.dirname(sys.executable)
    if not all(os.access(os.path.join(env_bin, exe), os.X_OK) for exe in _AMBER_EXES):
        return False
    os.environ["PATH"] = env_bin + os.pathsep + os.environ.get("PATH", "")
    LOGGER.debug("added %s to PATH so AmberTools is reachable", env_bin)
    return all(shutil.which(exe) for exe in _AMBER_EXES)


def amber_tools_available() -> bool:
    """True when ``antechamber`` and ``parmchk2`` are both usable."""
    return ensure_amber_tools_on_path()


def nagl_available() -> bool:
    'True when the NAGL charge backend is installed.'
    try:
        return importlib.util.find_spec("openff.nagl") is not None
    except (ImportError, ValueError):     # namespace package half-installed
        return False


def _version(distribution: str) -> str:
    """An installed distribution's version, or ``"absent"``.

    AmberTools is identified by executable and table content even when Python
    distribution metadata exists. Presence or a package version cannot detect
    replacements made inside the same installation.
    """
    if distribution == "ambertools":
        amber_tools_available()  # expose interpreter-adjacent tools on PATH
        return ambertools_identity()
    try:
        return importlib.metadata.version(distribution)
    except Exception:  # noqa: BLE001 - PackageNotFound and metadata corruption alike
        return _ABSENT


def backend_versions(forcefield: str, charge_method: str) -> str:
    """The producers on this resolved path, with versions/content identities.

    Both arguments are *resolved* values, not requested ones.  Only the
    packages that actually make the numbers for that path are listed, so
    upgrading NAGL does not invalidate a cache full of Gasteiger records --
    and upgrading ``openmmforcefields``, which is what writes the GAFF
    templates, does.

    >>> backend_versions("openff-2.2.0", "gasteiger")   # doctest: +ELLIPSIS
    'openff-toolkit=...; openmmforcefields=...; rdkit=...'
    """
    names = list(_ALWAYS_BACKENDS)
    if str(forcefield or "").strip().lower().startswith("gaff"):
        names.append("ambertools")
    for extra in _CHARGE_BACKENDS.get(charge_method_family(charge_method), ()):
        if extra not in names:
            names.append(extra)
    result = "; ".join(f"{name}={_version(name)}" for name in names)
    if charge_method_family(charge_method) == FIXED_AMBER_METHOD:
        from .fixed_amber import PROTOCOL

        result += f"; charge_protocol={PROTOCOL}"
    return result


def load_molecule(sdf_path: Path, allow_undefined_stereo: bool = True) -> Any:
    """Read one molecule from an SDF, keeping bond orders and formal charges.

    Multi-molecule SDFs are rejected: a pose file must be unambiguous, and
    silently taking record 0 is how the wrong ligand ends up being scored.
    """
    from openff.toolkit import Molecule

    sdf_path = Path(sdf_path)
    if not sdf_path.is_file():
        raise LigandParameterisationError(f"ligand SDF not found: {sdf_path}")
    try:
        loaded = Molecule.from_file(
            str(sdf_path), file_format="SDF",
            allow_undefined_stereo=allow_undefined_stereo,
        )
    except Exception as exc:  # noqa: BLE001 - surface any reader failure uniformly
        raise LigandParameterisationError(
            f"could not parse {sdf_path.name} as SDF: {exc}") from exc

    mols = loaded if isinstance(loaded, list) else [loaded]
    if not mols:
        raise LigandParameterisationError(f"{sdf_path.name} contains no molecules")
    if len(mols) > 1:
        raise LigandParameterisationError(
            f"{sdf_path.name} contains {len(mols)} molecules; pass a single pose")

    mol = mols[0]
    if mol.n_conformers == 0:
        raise LigandParameterisationError(
            f"{sdf_path.name} has no 3D coordinates; MM/GBSA needs a pose")
    if any(atom.atomic_number == 0 for atom in mol.atoms):
        raise LigandParameterisationError(
            f"{sdf_path.name} contains a dummy/query atom")
    return mol


def nagl_model_identity(model: str) -> Optional[str]:
    """Hash an installed model without fetching, loading or accepting ambiguity."""
    from ..dock.software import file_identity
    try:
        from openff.nagl_models import get_nagl_model_dirs_paths
        paths = [Path(directory) / model for directory in get_nagl_model_dirs_paths()]
        digests = {file_identity(path) for path in paths if path.is_file()}
        return next(iter(digests)) if len(digests) == 1 and None not in digests else None
    except (ImportError, OSError, ValueError, RuntimeError):
        return None


def _nagl_charges(mol: Any) -> str:
    from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper

    wrapper = NAGLToolkitWrapper()
    last: Optional[Exception] = None
    for model in NAGL_MODELS:
        try:
            before = nagl_model_identity(model)
            mol.assign_partial_charges(model, toolkit_registry=wrapper)
            after = nagl_model_identity(model)
            if isinstance(getattr(mol, "properties", None), dict):
                mol.properties["skinpipe_charge_model_sha256"] = after if before == after else None
            return f"nagl:{model}"
        except Exception as exc:  # noqa: BLE001 - try the next model
            last = exc
    raise LigandParameterisationError(f"NAGL charge assignment failed: {last}")


def _charge_values(smiles: str, method: str) -> tuple[str, Any]:
    """Run the existing adaptive charge protocol in an owned subprocess.

    The molecule crosses the process boundary as a mapped SMILES, which
    preserves atom order but not the pose, so a fresh conformer is generated.
    Standard AM1BCC uses AM1 Mulliken charges plus bond charge corrections on
    an internally generated conformer. This path does not request ELF10
    conformer selection/averaging and does not charge the docked geometry.
    """
    try:
        from openff.toolkit import Molecule

        mol = Molecule.from_mapped_smiles(smiles, allow_undefined_stereo=True)
        mol.generate_conformers(n_conformers=1)
        mol.assign_partial_charges(method)
        values = [float(q.m) for q in mol.partial_charges]
        if not values or not all(math.isfinite(value) for value in values):
            raise LigandParameterisationError("charge producer returned a non-finite or empty vector")
        return "ok", values
    except Exception as exc:  # noqa: BLE001 - report to the parent, never raise
        return "error", f"{type(exc).__name__}: {exc}"


def _external_charges(mol: Any, method: str, timeout_s: float) -> None:
    """Assign ``method`` charges in a subprocess with a hard timeout.

    The mapped-SMILES/conformer/charge protocol is unchanged. The calculation
    guardian owns the worker and its ordinary external descendants, including
    ``sqm``. A file response avoids joining a producer blocked on a queue pipe
    and the unreliable ``Queue.empty()`` check.
    """
    from openff.units import unit as off_unit
    from ..processes import run_calculation

    with tempfile.TemporaryDirectory(prefix="skinpipe-adaptive-charge-") as directory:
        root = Path(directory)
        request, response = root / "request.json", root / "response.json"
        request.write_text(json.dumps({"smiles": mol.to_smiles(mapped=True), "method": method}))
        try:
            with (root / "worker.log").open("w") as log:
                result = run_calculation(
                    [sys.executable, "-m", "skinpipe.mmgbsa.charge_worker", str(request), str(response)],
                    timeout=timeout_s, stdout=log, stderr=log)
        except subprocess.TimeoutExpired as exc:
            raise LigandParameterisationError(
                f"{method} charge assignment exceeded {timeout_s:.0f} s") from exc
        except OSError as exc:
            raise LigandParameterisationError(f"{method} charge assignment could not start: {exc}") from exc
        if result.returncode != 0:
            raise LigandParameterisationError(
                f"{method} charge assignment exited with status {result.returncode}")
        try:
            record = json.loads(response.read_text())
            if not isinstance(record, dict) or record.get("status") not in ("ok", "error"):
                raise ValueError("invalid worker response")
            status, payload = record["status"], record["payload"]
            if status == "ok" and not isinstance(payload, list):
                raise ValueError("charges are not a list")
        except (OSError, ValueError, KeyError) as exc:
            raise LigandParameterisationError(f"{method} charge assignment died silently: {exc}") from exc
    if status != "ok":
        raise LigandParameterisationError(f"{method} charge assignment failed: {payload}")
    if len(payload) != mol.n_atoms:
        raise LigandParameterisationError(
            f"{method} returned {len(payload)} charges for {mol.n_atoms} atoms")
    if not payload or any(isinstance(value, bool) or not isinstance(value, (int, float))
                          or not math.isfinite(value) for value in payload):
        raise LigandParameterisationError(f"{method} returned non-finite or non-numeric charges")
    mol.partial_charges = payload * off_unit.elementary_charge


def _fixed_amber_charges(mol: Any, timeout_s: float,
                         charge_workdir: Optional[Path] = None) -> None:
    """Explicit prospective protocol; the existing adaptive default is retained.

    Conformers follow the inspected OpenFF charge worker/wrapper sequence.
    External commands are directly owned by the fixed producer, so its total
    command timeout terminates the antechamber/SQM process group. The receipt
    reaches the scored record; a supplied workdir retains all raw files too.
    """
    from openff.toolkit import Molecule
    from openff.toolkit.utils.rdkit_wrapper import RDKitToolkitWrapper
    from openff.units import unit as off_unit

    from .fixed_amber import charge_sdf

    with ExitStack() as stack:
        if charge_workdir is None:
            root = Path(stack.enter_context(tempfile.TemporaryDirectory(prefix="skinpipe-fixed-charge-")))
        else:
            charge_workdir = Path(charge_workdir)
            charge_workdir.mkdir(parents=True, exist_ok=True)
            root = Path(tempfile.mkdtemp(prefix="fixed-charge-", dir=charge_workdir))
        worker = Molecule.from_mapped_smiles(mol.to_smiles(mapped=True), allow_undefined_stereo=True)
        worker.generate_conformers(n_conformers=1)
        generated = Molecule(worker)
        generated.generate_conformers(n_conformers=1, rms_cutoff=0.25 * off_unit.angstrom,
                                     toolkit_registry=RDKitToolkitWrapper())
        sdf = root / "generated.sdf"
        RDKitToolkitWrapper().to_file(generated, str(sdf), file_format="sdf")
        formal = generated.total_charge.m_as(off_unit.elementary_charge)
        if formal != int(formal):
            raise LigandParameterisationError("fixed AM1BCC requires an integer formal charge")
        receipt = charge_sdf(sdf, int(formal), "am1bcc", root / "raw", timeout_s)
        if receipt["n_atoms"] != mol.n_atoms:
            raise LigandParameterisationError("fixed AM1BCC changed the ligand atom count")
        mol.partial_charges = receipt["normalized_charges_e"] * off_unit.elementary_charge
        receipt["raw_directory"] = str(root / "raw") if charge_workdir is not None else None
        mol.properties["skinpipe_charge_producer"] = receipt


def charge_ladder(method: str = "auto") -> list[str]:
    """The order :func:`assign_charges` will try methods in.

    One list, used by the run and by the cache key's prediction, so the two
    cannot disagree about what "auto" means or about whether an explicit
    request falls back. Legacy explicit methods try the remaining ladder;
    ``am1bcc-fixed`` must either execute its recorded protocol or fail.
    """
    requested = (method or "auto").lower()
    if requested == FIXED_AMBER_METHOD:
        return [requested]
    if requested == "auto":
        return list(CHARGE_LADDER)
    return [requested] + [m for m in CHARGE_LADDER if m != requested]


def charge_method_installed(candidate: str) -> bool:
    'Validate the selected structure and retain explicit calculation provenance.'
    if candidate == "nagl":
        return nagl_available()
    if candidate in _EXTERNAL_CHARGE_METHODS or candidate == FIXED_AMBER_METHOD:
        return amber_tools_available()
    return True


def resolve_charge_method(method: str = "auto") -> str:
    'Validate the selected structure and retain explicit calculation provenance.'
    ladder = charge_ladder(method)
    for candidate in ladder:
        if charge_method_installed(candidate):
            return candidate
    return ladder[-1] if ladder else CHARGE_LADDER[-1]


def assign_charges(mol: Any, method: str = "auto",
                   timeout_s: float = 600.0,
                   charge_workdir: Optional[Path] = None) -> tuple[str, list[str]]:
    """Assign partial charges in place.  Returns ``(method_used, warnings)``.

    ``"auto"`` walks :data:`CHARGE_LADDER` and takes the first method that
    works: NAGL (AM1BCC-quality, ~0.5 s), then real AM1BCC if AmberTools is
    present, then Gasteiger as a last resort.  Gasteiger charges are adequate
    for rank-ordering but noticeably degrade absolute electrostatics, so the
    substitution is always reported -- in the returned warnings, which reach
    the caller through ``notes``, and not only in the log, where a batch
    running under a quiet logger would lose it.
    """
    warnings: list[str] = []
    requested = (method or "auto").lower()
    ladder = charge_ladder(requested)
    wanted = ladder[0]

    for candidate in ladder:



        if not charge_method_installed(candidate):
            warnings.append(f"charge method {candidate} unavailable "
                            f"(backend not installed)")
            continue
        try:
            if candidate == "nagl":
                used = _nagl_charges(mol)
            elif candidate == FIXED_AMBER_METHOD:
                _fixed_amber_charges(mol, timeout_s, charge_workdir)
                used = candidate
            elif candidate in _EXTERNAL_CHARGE_METHODS:
                _external_charges(mol, candidate, timeout_s)
                used = candidate
            else:
                mol.assign_partial_charges(candidate)
                used = candidate
        except Exception as exc:  # noqa: BLE001 - fall through the ladder
            LOGGER.debug("charge method %s unavailable: %s", candidate, exc)
            warnings.append(f"charge method {candidate} failed ({exc})")
            continue
        if candidate != wanted:
            msg = f"charge method {wanted} unavailable; substituted {used}"
            LOGGER.warning("%s", msg)
            warnings.append(msg)
        return used, warnings

    raise LigandParameterisationError(
        "no charge method succeeded: " + "; ".join(warnings))


def resolve_forcefield(forcefield: str,
                       allow_fallback: bool = True) -> tuple[str, list[str]]:
    'Validate the selected structure and retain explicit calculation provenance.'
    return _resolve_forcefield(forcefield, allow_fallback)


def _resolve_forcefield(forcefield: str, allow_fallback: bool) -> tuple[str, list[str]]:
    name = (forcefield or "").strip().lower()
    if not name.startswith("gaff"):
        return name, []
    if amber_tools_available():
        return name, []
    if not allow_fallback:
        raise LigandParameterisationError(
            f"{name} requires AmberTools (antechamber/parmchk2) which is not installed")
    msg = (f"{name} needs AmberTools which is not installed; "
           f"substituted {SMIRNOFF_FALLBACK}")
    LOGGER.warning("%s", msg)
    return SMIRNOFF_FALLBACK, [msg]


def _validate_partial_charges(mol: Any) -> None:
    """Refuse an invalid charge vector before any force-field generator sees it.

    This gate never repairs or neutralises charges: each producer keeps its
    own numerical protocol, including the fixed Amber normalization.
    """
    try:
        values = [charge.m for charge in mol.partial_charges]
        formal = mol.total_charge.m
        if (len(values) != mol.n_atoms or not values or
                any(isinstance(value, bool) or not math.isfinite(value) for value in values) or
                isinstance(formal, bool) or not math.isfinite(formal)):
            raise ValueError("charges and formal charge must be finite, with one charge per atom")
        net = float(sum(values))
        if not math.isfinite(net) or abs(net - formal) > _NET_CHARGE_TOL:
            raise ValueError(f"partial charges sum to {net:+.4f}, formal charge {formal:+.0f}")
    except (AttributeError, TypeError, ValueError, OverflowError) as exc:
        raise LigandParameterisationError(f"invalid partial charges: {exc}") from exc


def parameterise_ligand(sdf_path: Path, forcefield: str = "gaff-2.11",
                        charge_method: str = "auto",
                        charge_timeout_s: float = 600.0,
                        allow_forcefield_fallback: bool = True,
                        charge_workdir: Optional[Path] = None) -> LigandParameters:
    """Build an OpenMM residue-template generator for one ligand pose.

    Returns a :class:`LigandParameters` whose ``.generator`` is a
    GAFF/SMIRNOFF template generator; register it on an
    :class:`openmm.app.ForceField` with ``ff.registerTemplateGenerator(
    params.generator.generator)``.

    Parameters
    ----------
    sdf_path:
        Single-molecule SDF with 3D coordinates (a docked pose).
    forcefield:
        ``gaff-2.11`` (default, needs AmberTools) or an OpenFF line such as
        ``openff-2.2.0``.
    charge_method:
        ``auto`` (default) walks :data:`CHARGE_LADDER`.  Explicit values:
        ``nagl``, ``am1bcc``, ``gasteiger``. ``am1bcc-fixed`` selects the
        prospective SQM protocol and records its actual charge producer.
    """
    from openmmforcefields.generators import (GAFFTemplateGenerator,
                                              SMIRNOFFTemplateGenerator)

    mol = load_molecule(Path(sdf_path))
    resolved, warnings = _resolve_forcefield(forcefield, allow_forcefield_fallback)
    used_charges, charge_warnings = assign_charges(mol, charge_method, charge_timeout_s,
                                                charge_workdir=charge_workdir)
    warnings.extend(charge_warnings)
    charge_producer = mol.properties.pop("skinpipe_charge_producer", None)

    _validate_partial_charges(mol)
    formal = mol.total_charge.m

    cls = GAFFTemplateGenerator if resolved.startswith("gaff") else SMIRNOFFTemplateGenerator
    try:
        generator = cls(molecules=[mol], forcefield=resolved)
    except Exception as exc:  # noqa: BLE001
        raise LigandParameterisationError(
            f"{resolved} could not parameterise {Path(sdf_path).name}: {exc}") from exc

    return LigandParameters(
        molecule=mol, generator=generator, forcefield=resolved,
        charge_method=used_charges, net_charge=int(round(formal)),
        n_atoms=mol.n_atoms, warnings=warnings,
        charge_producer=charge_producer,
        charge_model_sha256=getattr(mol, "properties", {}).get("skinpipe_charge_model_sha256"),
    )
