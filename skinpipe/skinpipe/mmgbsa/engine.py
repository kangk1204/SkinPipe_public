"""Energy evaluation, minimisation and MD sampling for MM/GBSA.

Kept apart from the public API so that :mod:`openmm_gbsa` stays readable: this
module knows about Contexts and force groups, that one knows about poses and
DataFrames.
"""
from __future__ import annotations

import logging
import os
from dataclasses import dataclass
from typing import Any, Optional, Sequence

from .system import GB_GROUP, KJ_PER_KCAL, MM_GROUP, ComplexSystem

LOGGER = logging.getLogger(__name__)

#: Imported for the platform warning below, which quotes it.
from ..score import GBSA_RESOLUTION_KCAL  # noqa: E402

#: RMS force below which a restrained minimisation is called converged,
#: kJ/mol/nm.  The conversion is 1 kcal/mol/A = 4.184 / 0.1 = 41.84 kJ/mol/nm,
#: so a tight 0.5 kcal/mol/A bar would be ~20.9 kJ/mol/nm and this 100.0 is
#: ~2.39 kcal/mol/A -- an order of magnitude looser, and deliberately so.  The
#: poses are restrained and the budget is a few hundred L-BFGS steps, so this
#: says "the geometry is not grossly strained", not "the gradient vanished".
CONVERGENCE_RMS_FORCE = 100.0

#: L-BFGS tolerance handed to OpenMM, kJ/mol/nm.
MINIMISE_TOLERANCE = 10.0

#: Langevin timestep.  2 fs needs H-bond constraints; without them 1 fs is the
#: largest safe step for a system containing X-H stretches.
TIMESTEP_CONSTRAINED_PS = 0.002
TIMESTEP_UNCONSTRAINED_PS = 0.001
DEFAULT_FRICTION_PER_PS = 1.0

ENERGY_KEYS = ("dG_bind", "e_gas", "e_solv", "e_polar", "e_nonpolar",
               "e_complex", "e_receptor", "e_ligand")


@dataclass(frozen=True)
class Terms:
    """One species' energy, kcal/mol, split by force group."""

    gas: float
    polar: float
    nonpolar: float

    @property
    def total(self) -> float:
        return self.gas + self.polar + self.nonpolar

    def __sub__(self, other: "Terms") -> "Terms":
        return Terms(self.gas - other.gas, self.polar - other.polar,
                     self.nonpolar - other.nonpolar)


#: Probe result per platform name. A probe costs one tiny Context; the answer
#: cannot change within a process, and `mmgbsa_batch` asks once per pose.
_PLATFORM_PROBE: "dict[str, Optional[str]]" = {}

#: Tried in order when the requested platform will not run. Reference is last
#: and always works; OpenCL before CPU because it is the one that keeps a GPU
#: run a GPU run.
_PLATFORM_FALLBACKS = ("OpenCL", "CPU", "Reference")


def platform_error(name: str) -> Optional[str]:
    'Why this OpenMM platform will not run here, or None if it will.'
    if name in _PLATFORM_PROBE:
        return _PLATFORM_PROBE[name]
    import openmm as mm
    from openmm import unit as u

    error: Optional[str] = None
    try:
        system = mm.System()
        for _ in range(8):
            system.addParticle(1.0 * u.amu)
        force = mm.NonbondedForce()
        for _ in range(8):
            force.addParticle(0.0, 0.3 * u.nanometer, 1.0 * u.kilojoule_per_mole)
        system.addForce(force)
        context = mm.Context(system, mm.VerletIntegrator(1.0 * u.femtosecond),
                             mm.Platform.getPlatformByName(name))
        context.setPositions([(0.1 * i, 0.0, 0.0) for i in range(8)])
        context.getState(getEnergy=True).getPotentialEnergy()
        del context
    except Exception as exc:                      # noqa: BLE001 - probe
        error = f"{type(exc).__name__}: {exc}"
    _PLATFORM_PROBE[name] = error
    return error


def cpu_threads(threads: Optional[int] = None, cpu_profile: str = "auto") -> int:
    """Validate and resolve the CPU execution request independently of OpenMM.

    ``auto`` defaults to one thread and honours an explicit count. ``serial``
    requires one thread; ``parallel`` opts into a supplied count or all cores.
    These names describe execution, not a guarantee of numerical determinism.
    """
    if cpu_profile not in ("auto", "serial", "parallel"):
        raise ValueError("cpu_profile must be auto, serial or parallel")
    if threads is not None:
        if isinstance(threads, bool) or not isinstance(threads, int):
            raise TypeError("threads must be a positive integer or None")
        if threads < 1:
            raise ValueError("threads must be a positive integer")
    if cpu_profile == "serial" and threads not in (None, 1):
        raise ValueError("cpu_profile serial requires threads=1 or None")
    if threads is not None:
        return threads
    return max(1, os.cpu_count() or 1) if cpu_profile == "parallel" else 1


def execution_conditions(platform: Any, properties: dict[str, str],
                         context: Any = None) -> dict[str, Any]:
    """Cache-relevant CPU properties, read from the Context when available.

    GPU device/driver/precision reproducibility is a separate contract. The
    CPU properties here are compared along with the resolved platform and
    OpenMM version; temporary directories never participate in reuse.
    """
    import openmm

    name = str(platform.getName())
    effective: dict[str, str] = {}
    if name == "CPU":
        for key in platform.getPropertyNames():
            effective[str(key)] = str(
                platform.getPropertyValue(context, key) if context is not None
                else properties.get(key, platform.getPropertyDefaultValue(key)))
    count = int(effective["Threads"]) if "Threads" in effective else None
    return {"resolved_platform": name,
            "resolved_platform_properties": effective,
            "resolved_cpu_profile": ("serial" if count == 1 else "parallel")
            if count is not None else None,
            "openmm_version": str(openmm.__version__)}


def make_platform(name: Optional[str] = None,
                  threads: Optional[int] = None,
                  cpu_profile: str = "auto") -> tuple[Any, dict[str, str]]:
    'Pick an OpenMM platform. Defaults to CPU; `name` selects another.'
    count = cpu_threads(threads, cpu_profile)
    from openmm import Platform

    available = {Platform.getPlatform(i).getName()
                 for i in range(Platform.getNumPlatforms())}
    chosen = (name or os.environ.get("SKINPIPE_OPENMM_PLATFORM")
              or ("CPU" if "CPU" in available else "Reference"))
    if chosen not in available:
        # A name OpenMM does not know is a typo or a wrong environment, not a
        # platform that is merely busy -- and falling back put `--platform CUDA`
        # misspelt onto Reference, which is correct arithmetic at a hundredth of
        # the speed and no warning a reader would connect to a stalled run.
        # `registered but unusable` keeps its fallback ladder below; this is the
        # other case and it is refused.
        raise ValueError(
            f"{chosen!r} is not an OpenMM platform in this environment. "
            f"Available: {sorted(available)}. A misspelt name used to fall back "
            "to Reference, which runs correctly and roughly a hundred times "
            "slower.")
    # Registered is not the same as working: OpenMM registers CUDA whenever it
    # finds the library, and the module fails at load if the toolkit's PTX is
    # newer than the driver takes. Ask before the run rather than once per pose.
    failure = platform_error(chosen)
    if failure:
        for alternative in _PLATFORM_FALLBACKS:
            if alternative in available and alternative != chosen \
                    and not platform_error(alternative):
                LOGGER.warning(
                    "OpenMM platform %s is registered but will not run here "
                    "(%s); falling back to %s. This is not a rounding "
                    "difference to assume away. A historical report described "
                    "a CUDA/OpenCL difference of up to 13.3 kcal/mol; that "
                    "comparison has not been independently reproduced. The "
                    "scorer's %.1f kcal/mol resolution does not validate a "
                    "platform change or repeat stability. Do not "
                    "mix this table with one scored on %s; rescore it.",
                    chosen, failure, alternative, GBSA_RESOLUTION_KCAL, chosen)
                chosen = alternative
                break
        else:
            raise RuntimeError(
                f"OpenMM platform {chosen} will not run here ({failure}) and "
                f"no fallback in {_PLATFORM_FALLBACKS} does either")
    properties: dict[str, str] = {}
    platform = Platform.getPlatformByName(chosen)
    if chosen == "CPU":
        properties["Threads"] = str(count)
        # Pin other CPU defaults too. They are otherwise mutable process-wide
        # state, and must be the same in the parent key and worker Contexts.
        for key in platform.getPropertyNames():
            if key != "Threads":
                properties[str(key)] = str(platform.getPropertyDefaultValue(key))
    return platform, properties


class Evaluator:
    """The six Contexts needed to decompose one complex geometry.

    Contexts are expensive to build and cheap to re-position, so they are
    created once and reused for every MD snapshot.  Two Contexts per species:
    one System carries the ACE nonpolar term, one does not, and the difference
    of their GB force groups isolates the nonpolar contribution.
    """

    def __init__(self, cs: ComplexSystem, platform: Any, properties: dict[str, str]):
        import openmm

        self._cs = cs
        self._contexts: dict[str, tuple[Any, Any]] = {}
        for label, species in (("complex", cs.complex), ("receptor", cs.receptor),
                               ("ligand", cs.ligand)):
            self._contexts[label] = (
                openmm.Context(species.full, openmm.VerletIntegrator(0.001),
                               platform, properties),
                openmm.Context(species.polar_only, openmm.VerletIntegrator(0.001),
                               platform, properties),
            )

    @property
    def complex_context(self) -> Any:
        """The Context minimisation runs in; it owns the restraint force."""
        return self._contexts["complex"][0]

    def _terms(self, label: str, positions: Any) -> Terms:
        ctx_full, ctx_polar = self._contexts[label]
        ctx_full.setPositions(positions)
        ctx_polar.setPositions(positions)
        gas = ctx_full.getState(getEnergy=True, groups={MM_GROUP}
                                ).getPotentialEnergy()._value / KJ_PER_KCAL
        gb_ace = ctx_full.getState(getEnergy=True, groups={GB_GROUP}
                                   ).getPotentialEnergy()._value / KJ_PER_KCAL
        gb_polar = ctx_polar.getState(getEnergy=True, groups={GB_GROUP}
                                      ).getPotentialEnergy()._value / KJ_PER_KCAL
        return Terms(gas=gas, polar=gb_polar, nonpolar=gb_ace - gb_polar)

    def decompose(self, positions: Any) -> dict[str, float]:
        """dG and its parts for one complex geometry.

        Receptor and ligand energies are taken from slices of these same
        coordinates -- the single-trajectory approximation.  Restraint energy
        sits in its own force group and is never read here.
        """
        cplx = self._terms("complex", positions)
        rec = self._terms("receptor", positions[self._cs.receptor_slice])
        lig = self._terms("ligand", positions[self._cs.ligand_slice])
        delta = cplx - rec - lig
        # dG is built from the two reported halves rather than from
        # ``delta.total`` so that ``e_gas + e_solv == dG_bind`` holds exactly in
        # floating point, not merely to within an association-order rounding
        # error. The module documents that identity; this is what makes it true.
        e_gas = delta.gas
        e_solv = delta.polar + delta.nonpolar
        return {
            "dG_bind": e_gas + e_solv,
            "e_gas": e_gas,
            "e_solv": e_solv,
            "e_polar": delta.polar,
            "e_nonpolar": delta.nonpolar,
            "e_complex": cplx.total,
            "e_receptor": rec.total,
            "e_ligand": lig.total,
        }


def minimise(context: Any, max_iterations: int) -> dict[str, float]:
    'Restrained L-BFGS minimisation; returns convergence diagnostics.'
    import numpy as np
    import openmm

    e_start = context.getState(getEnergy=True).getPotentialEnergy()._value
    openmm.LocalEnergyMinimizer.minimize(context, MINIMISE_TOLERANCE, max_iterations)
    state = context.getState(getEnergy=True, getForces=True)
    norms = np.linalg.norm(np.asarray(state.getForces(asNumpy=True)._value), axis=1)
    return {
        "e_before_kcal": e_start / KJ_PER_KCAL,
        "e_after_kcal": state.getPotentialEnergy()._value / KJ_PER_KCAL,
        "rms_force": float(np.sqrt((norms ** 2).mean())),
        "max_force": float(norms.max()),
    }


def md_snapshots(cs: ComplexSystem, evaluator: Evaluator, platform: Any,
                 properties: dict[str, str], md_ps: float, n_snapshots: int,
                 temperature_k: float, seed: int) -> list[dict[str, float]]:
    """Langevin dynamics on the restrained complex, sampled evenly.

    Restraints stay on during dynamics: with an implicit solvent and no
    explicit water shell a free protein drifts away from the crystal geometry
    within a few picoseconds, which shows up as receptor-internal noise rather
    than as binding signal.  The first 10% of the run is discarded as
    equilibration.

    ``seed`` follows OpenMM's convention: 0 means "pick a random seed", so
    repeated runs will not agree.  Pass a nonzero seed for a reproducible
    trajectory.
    """
    import openmm
    from openmm import unit

    dt = (TIMESTEP_CONSTRAINED_PS if cs.complex.full.getNumConstraints()
          else TIMESTEP_UNCONSTRAINED_PS)
    n_steps = max(1, int(round(md_ps / dt)))
    LOGGER.info("%.3g ps of Langevin dynamics at %.0f K (%d steps of %.0f fs)",
                md_ps, temperature_k, n_steps, dt * 1000)
    integrator = openmm.LangevinMiddleIntegrator(
        temperature_k * unit.kelvin,
        DEFAULT_FRICTION_PER_PS / unit.picosecond,
        dt * unit.picoseconds)
    integrator.setRandomNumberSeed(seed)
    context = openmm.Context(cs.complex.full, integrator, platform, properties)
    context.setPositions(
        evaluator.complex_context.getState(getPositions=True).getPositions())
    context.setVelocitiesToTemperature(temperature_k * unit.kelvin, seed)

    # The trajectory has to be the length that was asked for. `stride =
    # max(1, (n_steps - equil) // n_snapshots)` stepped `n_snapshots * stride`
    # regardless, so integer division lost the remainder and the `max(1, ...)`
    # floor overshot when the request was short. Measured against a step-counting
    # integrator, dt = 0.001 ps:
    #
    #     md_ps  snapshots  requested  actually ran
    #     0.001      5           1          6
    #     0.01      20          10         21
    #     0.01       5          10          6
    #     0.5        4         500        498
    #
    # The default `md_ps=0` path never reaches here, so this is about the short
    # non-default runs -- where the discrepancy is largest.
    equil = max(1, min(n_steps, n_steps // 10))
    integrator.step(equil)
    sampling = n_steps - equil
    if n_snapshots > max(sampling, 1):
        # More frames than there are steps to put between them. Taking them
        # anyway means `max(1, ...)` invents the steps.
        LOGGER.warning("md_ps %.4g at dt %.4g ps leaves %d step(s) after "
                       "equilibration, which cannot carry %d snapshots; taking "
                       "%d", md_ps, dt, sampling, n_snapshots, max(sampling, 1))
        n_snapshots = max(sampling, 1)
    # Spread the remainder over the first frames so the strides sum to exactly
    # the sampling budget.
    base, extra = divmod(sampling, n_snapshots)
    strides = [base + (1 if i < extra else 0) for i in range(n_snapshots)]

    frames: list[dict[str, float]] = []
    for stride in strides:
        if stride:
            integrator.step(stride)
        frames.append(evaluator.decompose(
            context.getState(getPositions=True).getPositions()))
    return frames


def mean_std(values: Sequence[float]) -> tuple[float, Optional[float]]:
    """Mean and sample standard deviation; ``None`` for a single value."""
    import numpy as np

    arr = np.asarray(values, dtype=float)
    if arr.size < 2:
        return float(arr[0]), None
    return float(arr.mean()), float(arr.std(ddof=1))
