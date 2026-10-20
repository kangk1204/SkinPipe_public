'What must match before a cached MM/GBSA record may be reused.'
from __future__ import annotations

import inspect
import json
import logging
import re
from pathlib import Path
from typing import Any, Callable, Optional

# The same questions -- "is this the same receptor?", "is this the same file?"
# -- with the same answers, so they are imported rather than written twice.
# skinpipe.dock.cachekey pulls in nothing heavier than hashlib and
# skinpipe.dock.__init__ is empty, so this does not drag Vina or RDKit into an
# MM/GBSA-only run.
from ..dock.cachekey import file_sha256
# Imported as a module, not as names, so that a test (or a machine) moving
# ``amber_tools_available``/``nagl_available`` moves what this key predicts
# too: the key must resolve the path the same way the worker will.
from . import ligand as _ligand

LOGGER = logging.getLogger(__name__)
















CONDITIONS_VERSION = 8

#: Key under which the conditions live inside the cached record. The record
#: stays flat so the file is still readable as the row it produced.
CONDITIONS_FIELD = "cache_conditions"

#: Every ``_run_single`` parameter that changes the energy.
PROTOCOL_FIELDS = (
    "minimise_steps", "implicit_solvent", "ligand_forcefield", "n_snapshots",
    "md_ps", "salt_molar", "charge_method", "restraint_k", "temperature_k",
    "nonbonded_cutoff_nm", "constrain_h_bonds", "seed", "platform",
    "preparation_seed", "threads", "cpu_profile",
)

#: Parameters kept out of the key, and why. Anything that is in neither this
#: mapping nor :data:`PROTOCOL_FIELDS` makes :func:`protocol_conditions` raise.
NOT_IN_KEY = {
    "receptor_pdb": "hashed by content instead",
    "ligand_sdf": "hashed by content instead",
    "workdir": "where the answer is written, not what it is",
    "ligand_id": "carried in the key separately",
    "pose_index": "carried in the key separately",
    "prepare_receptor_pdb": "its effect is visible in the receptor content hash",
    # Safe to omit only because mmgbsa_batch runs the metal guard before it
    # reads the cache: under the default "refuse" the batch raises and nothing
    # is served, so a record written under "score_without_metal" can never be
    # handed back to a caller who asked to refuse. It changes no energy either
    # way -- the System is identical; only the note on the result differs.
    "on_metal_site": "guarded before the cache is read; changes notes, not energy",
}

#: What the worker will actually use, as opposed to what was asked for. Part
#: of the key: a machine that resolves Gasteiger must not be handed a record
#: computed on one that resolved AM1BCC.
RESOLVED_FIELDS = ("resolved_ligand_forcefield", "resolved_charge_method",
                   "resolved_charge_model", "charge_model_sha256", "ligand_backends")

#: Reserved for non-identity diagnostics. Charge-model identity is not evidence-only.
EVIDENCE_FIELDS: tuple[str, ...] = ()

#: The two shapes "we do not know which platform ran" arrives in.
_PLATFORM_UNKNOWN = (None, "<absent>")


def _UNKNOWN_PLATFORM(was: Any, now: Any) -> bool:
    return was in _PLATFORM_UNKNOWN or now in _PLATFORM_UNKNOWN

#: ``_run_single`` reads ``constrain_h_bonds`` as ``bool(...)``, so ``None``
#: and ``False`` are the same protocol and must hash the same.
_NORMALISE: dict[str, Callable[[Any], Any]] = {"constrain_h_bonds": bool}

def _plain(value: Any) -> Any:
    """JSON-comparable form: a value read back must equal the one written."""
    if isinstance(value, bool) or value is None:
        return value
    if isinstance(value, int):
        return int(value)
    if isinstance(value, float):
        return round(float(value), 6)
    if isinstance(value, Path):
        return str(value)
    return str(value)


def protocol_conditions(run_single: Callable[..., Any], kwargs: dict[str, Any],
                        receptor_pdb: Path) -> dict[str, Any]:
    """The protocol half of the key, hashing the receptor once for the batch.

    Defaults are read from ``run_single``'s own signature rather than copied
    here, so an omitted keyword records exactly what the worker will use and
    the two cannot drift apart.

    Raises
    ------
    RuntimeError
        If ``run_single`` grew a parameter that is in neither
        :data:`PROTOCOL_FIELDS` nor :data:`NOT_IN_KEY`. A new knob must be
        classified before it can be used, because the alternative is a key
        that quietly ignores it.
    """
    params = inspect.signature(run_single).parameters
    unclassified = sorted(set(params) - set(PROTOCOL_FIELDS) - set(NOT_IN_KEY))
    if unclassified:
        raise RuntimeError(
            f"{run_single.__name__} has parameter(s) {unclassified} that the "
            f"MM/GBSA cache key neither hashes nor excludes. Add them to "
            f"PROTOCOL_FIELDS if they change the energy, or to NOT_IN_KEY with "
            f"the reason they do not.")
    missing = sorted(set(PROTOCOL_FIELDS) - set(params))
    if missing:
        raise RuntimeError(
            f"the MM/GBSA cache key hashes {missing}, which "
            f"{run_single.__name__} no longer accepts")

    conditions: dict[str, Any] = {"version": CONDITIONS_VERSION,
                                  "receptor_atoms_sha256": file_sha256(receptor_pdb)}
    # Native modified receptors add parameters to the default protein model.
    # Bind their actual XML content, even if the package version stays equal.
    # Standard receptors retain their existing keys and physical model.
    if Path(receptor_pdb).is_file():
        from ..structure.modified import parameterization_descriptor
        native_model = parameterization_descriptor(Path(receptor_pdb).read_text())
        if native_model is not None:
            conditions["receptor_modified_model"] = native_model
    for name in PROTOCOL_FIELDS:
        value = kwargs.get(name, params[name].default)
        conditions[name] = _plain(_NORMALISE.get(name, lambda v: v)(value))
    conditions.update(resolved_conditions(conditions["ligand_forcefield"],
                                          conditions["charge_method"]))











    from .engine import cpu_threads
    cpu_threads(kwargs.get("threads"), kwargs.get("cpu_profile", "auto"))
    conditions.update(_resolved_execution(conditions["platform"],
                                         kwargs.get("threads"),
                                         kwargs.get("cpu_profile", "auto")))
    return conditions


def _resolved_execution(requested: Any, threads: Optional[int],
                        cpu_profile: str) -> dict[str, Any]:
    try:
        from .engine import execution_conditions, make_platform
        platform, properties = make_platform(requested, threads, cpu_profile)
        return execution_conditions(platform, properties)
    except Exception:                         # noqa: BLE001 - advisory probe
        return {"resolved_platform": None, "resolved_platform_properties": None,
                "resolved_cpu_profile": None, "openmm_version": None}


def _resolved_platform(requested: Any) -> Optional[str]:
    """The platform OpenMM would really use for this request, or None.

    None when it cannot be determined -- OpenMM missing, or every platform
    unusable. `describe_mismatch` treats None as "unknown", not as a value, so
    a record written where the answer was unknowable is not refused for it.
    """
    try:
        from .engine import make_platform
        platform, _properties = make_platform(requested)
        return str(platform.getName())
    except Exception:                        # noqa: BLE001 - advisory only
        return None


def resolved_conditions(ligand_forcefield: Any, charge_method: Any) -> dict[str, Any]:
    'Validate the selected structure and retain explicit calculation provenance.'
    forcefield, _substitution = _ligand.resolve_forcefield(ligand_forcefield)
    charges = _ligand.resolve_charge_method(charge_method)
    return {"resolved_ligand_forcefield": _plain(forcefield),
            "resolved_charge_method": _plain(charges),
            "ligand_backends": _ligand.backend_versions(forcefield, charges),
            "resolved_charge_model": (f"nagl:{_ligand.NAGL_MODELS[0]}"
                                      if charges == "nagl" else charges),
            "charge_model_sha256": (_ligand.nagl_model_identity(_ligand.NAGL_MODELS[0])
                                    if charges == "nagl" else "not_applicable")}


def observed_conditions(conditions: dict[str, Any],
                        record: dict[str, Any]) -> dict[str, Any]:
    """``conditions`` corrected by what the finished record actually used.

    The predicted path is the right answer for a machine, not for a pose: NAGL
    can be installed and still throw on one molecule, and that pose was then
    computed with Gasteiger charges while the prediction still said NAGL.
    Writing the prediction would make that record a *hit* for a later NAGL run
    -- which is the same defect one level down.

    A record with no provenance to read (a pose that failed before it was
    parameterised) keeps the predicted values: there is no force field and no
    charge method to record, and the record is a cached failure either way.
    """
    observed = dict(conditions)
    if "receptor_modified_model" in conditions:
        parameters = record.get("parameterization")
        actual = parameters.get("receptor_modified_model") if isinstance(parameters, dict) else None
        observed["receptor_modified_model"] = actual or {"status": "unverified: no actual native receptor model receipt"}
    execution = record.get("execution")
    if isinstance(execution, dict):
        for key in ("resolved_platform", "resolved_platform_properties",
                    "resolved_cpu_profile", "openmm_version"):
            if key in execution:
                observed[key] = execution[key]
    used = _ligand.parse_provenance(record.get("notes", ""))
    if not used:
        return observed
    forcefield = used.get("ligand_forcefield",
                          conditions.get("resolved_ligand_forcefield"))
    charges = used.get("charge_method", conditions.get("resolved_charge_method"))
    backends = _ligand.backend_versions(str(forcefield), str(charges))
    if charges == _ligand.FIXED_AMBER_METHOD:
        parameterization = record.get("parameterization")
        receipt = parameterization.get("charges") if isinstance(parameterization, dict) else None
        before = receipt.get("producer_before") if isinstance(receipt, dict) else None
        after = receipt.get("producer_after") if isinstance(receipt, dict) else None
        if (not isinstance(before, dict) or not isinstance(after, dict) or
                not before.get("verified") or before != after or
                not re.fullmatch(r"sha256:[0-9a-f]{64}", str(before.get("identity", "")))):
            backends += "; ligand_backend_identity=unverified:missing-or-changed-charge-producer"
        else:
            # Use the installation inspected around the actual fixed charge
            # commands, rather than relabelling it with the parent's later
            # inspection (including a fallback or installation restoration).
            backends = re.sub(r"\bambertools=[^;]*", "ambertools=" + before["identity"], backends)
            backends = re.sub(r"\bcharge_protocol=[^;]*", "charge_protocol=" + str(receipt.get("protocol")), backends)
    elif (forcefield == conditions.get("resolved_ligand_forcefield") and
            charges == conditions.get("resolved_charge_method") and
            backends != conditions.get("ligand_backends")):
        # A post-run installation fingerprint is not the producer identity if
        # the installation changed since the pre-run inspection. Keep the
        # discrepancy explicit and refuse reuse instead of relabelling it.
        backends += "; ligand_backend_identity=unverified:changed-during-run"
    return {**observed,
            "resolved_ligand_forcefield": _plain(forcefield),
            "resolved_charge_method": _plain(charges),
            "ligand_backends": backends,
            "resolved_charge_model": _plain(used.get("charge_model")),
            "charge_model_sha256": (used.get("charge_model_sha256")
                                     if charges == "nagl" else "not_applicable")}


def pose_conditions(protocol: dict[str, Any], ligand_id: str, pose_index: int,
                    sdf: Path) -> dict[str, Any]:
    """The protocol plus this pose's identity and its SDF's content hash."""
    return {**protocol, "ligand_id": str(ligand_id), "pose_index": int(pose_index),
            "pose_sha256": file_sha256(sdf)}


def describe_mismatch(cached: Any, current: dict[str, Any]) -> list[str]:
    """Every condition that differs, phrased as ``field was -> now``.

    An empty list means the entry may be reused. A record whose receptor or
    pose file is no longer on disk is refused outright: its inputs cannot be
    checked at all, and an unverifiable number is the thing this cache exists
    to stop being served.

    A NAGL record needs an explicit recognised model on both sides. Family
    equality alone cannot establish that its energy belongs to this request.
    """
    if cached is None:
        return ["written before the MM/GBSA protocol was recorded"]
    if not isinstance(cached, dict):






        return [f"the stored conditions are {type(cached).__name__} "
                f"({cached!r}), not the mapping protocol_conditions produces; "
                "the record is corrupt and cannot be verified"]
    for field, what in (("receptor_atoms_sha256", "receptor PDB"),
                        ("pose_sha256", "pose SDF")):
        if cached.get(field) is None or current.get(field) is None:
            return [f"the {what} is not on disk, so the record cannot be verified"]
    changed = []
    for label, conditions in (("cached", cached), ("current", current)):
        if (conditions.get("resolved_charge_method") == "nagl" and
                (conditions.get("resolved_charge_model") not in
                 {f"nagl:{model}" for model in _ligand.NAGL_MODELS}
                 or not re.fullmatch(r"[0-9a-f]{64}", str(conditions.get("charge_model_sha256", ""))))):
            changed.append(f"{label} NAGL model identity is unverified")
        backends = str(conditions.get("ligand_backends", ""))
        entries = dict(token.strip().split("=", 1) for token in backends.split(";")
                       if "=" in token)
        uses_amber = (str(conditions.get("resolved_ligand_forcefield", "")).lower().startswith("gaff") or
                      conditions.get("resolved_charge_method") in
                      {"am1bcc", "am1-mulliken", "am1bccelf10", _ligand.FIXED_AMBER_METHOD})
        if uses_amber and not re.fullmatch(r"sha256:[0-9a-f]{64}",
                                           entries.get("ambertools", "")):
            changed.append(f"{label} AmberTools content identity is unverified")
        if entries.get("ligand_backend_identity", "").startswith("unverified:"):
            changed.append(f"{label} ligand backend identity is unverified")
    for field in sorted((set(cached) | set(current)) - set(EVIDENCE_FIELDS)):
        was = cached.get(field, "<absent>")
        now = current.get(field, "<absent>")










        if field == "resolved_platform" and _UNKNOWN_PLATFORM(was, now):
            continue
        if was != now:
            changed.append(f"{field} {was!r} -> {now!r}")
    return changed


def platform_unverified(cached: Any, current: dict[str, Any]) -> Optional[str]:
    """Why this record's backend cannot be confirmed to match, or None.

    `describe_mismatch` returns `[]` for a record that never wrote
    `resolved_platform`, and returning `[]` is right -- refusing 120 records for
    a field that was unknowable when they were written would discard them for a
    reason that is not true. But `[]` reads as "verified the same", and it is
    not: it is "cannot be checked". The two deserve different words, and a
    caller that wants to say so in a report needs somewhere to read it.

    This is traceability, not arithmetic. Two *known* platforms that differ are
    still a mismatch and still refused by `describe_mismatch`.
    """
    if not isinstance(cached, dict):
        return None
    was = cached.get("resolved_platform", "<absent>")
    now = current.get("resolved_platform", "<absent>")
    if not _UNKNOWN_PLATFORM(was, now):
        return None
    if was == "<absent>":
        return ("this record predates `resolved_platform`, so the backend it "
                "was computed on was never written down and cannot be "
                "compared with this run's")
    if was is None:
        return ("this record's backend could not be determined when it was "
                "written (no OpenMM, or nothing usable), so it cannot be "
                "compared with this run's")
    return ("this run's backend could not be determined, so it cannot be "
            "compared with the record's")


def cache_path(workdir: Path, ligand_id: str, pose_index: int) -> Path:
    """Stable path for one pose's record.

    The conditions are *inside* the file rather than in its name on purpose:
    a stale entry is then found, read, and reported as stale with the field
    that changed, instead of being missed silently as if it had never existed.
    """
    safe = "".join(c if c.isalnum() or c in "-_." else "_" for c in str(ligand_id))
    return Path(workdir) / "cache" / f"{safe}__pose{int(pose_index)}.json"


def load_record(path: Path, conditions: dict[str, Any]) -> Optional[dict[str, Any]]:
    """The cached record if it was computed under ``conditions``, else ``None``.

    Every refusal is logged at WARNING naming the field that changed, so a
    recomputed batch says why it is recomputing.
    """
    path = Path(path)
    if not path.is_file():
        return None
    try:
        payload = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        LOGGER.warning("discarding unreadable cache %s: %s", path, exc)
        return None
    if not isinstance(payload, dict):
        LOGGER.warning("discarding cache %s: not a record", path)
        return None
    changed = describe_mismatch(payload.get(CONDITIONS_FIELD), conditions)
    if changed:
        LOGGER.warning("%s: the cached MM/GBSA record was computed under "
                       "different conditions (%s); recomputing",
                       path.name, "; ".join(changed))
        return None
    return {k: v for k, v in payload.items() if k != CONDITIONS_FIELD}


def save_record(path: Path, record: dict[str, Any],
                conditions: dict[str, Any]) -> None:
    """Write the record with the conditions that actually produced it alongside.

    The conditions are reconciled against the record's own provenance first
    (:func:`observed_conditions`), so the file says which force field and
    which charge method made this number rather than which ones were asked
    for.

    A caller that hands over something that is *not* the conditions mapping
    gets an ERROR naming what arrived. The file is still written, because
    refusing would only trade one unusable record for none, but such a record
    can never be verified again -- ``describe_mismatch`` refuses a non-mapping
    outright -- so every later run recomputes it and the resolved method never
    reaches disk. That is a silent permanent cache miss, and silence is how it
    would survive: nothing downstream fails, the batch just quietly never
    resumes.
    """
    path = Path(path)
    if not isinstance(conditions, dict):
        LOGGER.error("%s is being cached with conditions that arrived as %s "
                     "(%r) instead of the mapping from protocol_conditions; "
                     "this record can never be verified and every later run "
                     "will recompute it", path.name,
                     type(conditions).__name__, conditions)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps({**record, CONDITIONS_FIELD: conditions},
                                   indent=2))
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(
        {**record, CONDITIONS_FIELD: observed_conditions(conditions, record)},
        indent=2))
