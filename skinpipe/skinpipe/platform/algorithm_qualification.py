'Operator-owned chemistry asset admission and CAS-backed runtime qualification.'
from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import re
import time

from .artifacts import ArtifactIntegrityError, atomic_bytes, canonical_json, sha256_file
from .identity import chemical_identity
from .algorithm_manifest import decode_manifest, MAX_MANIFEST_BYTES

SCHEMA = "skinpipe-chemistry-qualification-v1"
MODEL_SCHEMA = "skinpipe-model-qualification-v1"
SOURCE_SCHEMA = "skinpipe-registered-algorithm-source-v1"
DEPENDENCIES = {"dimorphite": {"rdkit", "dimorphite-dl"}, "cdpl": {"rdkit", "CDPKit"},
                "posebusters": {"rdkit", "posebusters"}, "prolif": {"rdkit", "prolif", "MDAnalysis"},
                "fpsim2": {"rdkit", "FPSim2"}}
METHODS = {"dimorphite": "Dimorphite-DL+RDKit-tautomer",
           "cdpl": "CDPKit-CONFORGE+GaussianShapeAlignment+PharmacophoreAlignment",
           "posebusters": "PoseBusters-dock", "prolif": "ProLIF", "fpsim2": "FPSim2-Morgan"}
LICENSES = {"MIT", "Apache-2.0", "BSD-2-Clause", "BSD-3-Clause", "PSF-2.0",
            "LGPL-2.1-only", "LGPL-2.1-or-later", "LGPL-3.0-only", "LGPL-3.0-or-later"}
SHA = re.compile(r"[0-9a-f]{64}\Z")
MAX_BYTES = 64 * 1024 * 1024


class QualificationHeld(ValueError):
    """The optional runtime remains held; source candidates remain unchanged."""


def _hash(value):
    return hashlib.sha256(canonical_json(value)).hexdigest()


def _pairs(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ArtifactIntegrityError("Duplicate JSON object key")
        result[key] = value
    return result


def _json(raw, *, maximum=MAX_BYTES):
    if not isinstance(raw, bytes) or len(raw) > maximum:
        raise ArtifactIntegrityError("Qualification JSON exceeds its byte limit or is not bytes")
    value = json.loads(raw, object_pairs_hook=_pairs,
                       parse_constant=lambda value: (_ for _ in ()).throw(ArtifactIntegrityError("Nonfinite JSON: " + value)))
    if not isinstance(value, dict):
        raise ArtifactIntegrityError("Qualification JSON requires an object")
    return value


def _digest(value):
    if not isinstance(value, str) or not SHA.fullmatch(value):
        raise ValueError("Exact lowercase SHA256 required")
    return value


def _asset(path, expected, *, executable=False):
    if not isinstance(path, str) or not Path(path).is_absolute():
        raise ValueError("Administrator assets require absolute paths")
    supplied = Path(path)
    resolved = supplied.resolve()
    if not resolved.is_file():
        raise QualificationHeld("configured_asset_missing")
    if supplied.is_symlink() and not executable:
        raise QualificationHeld("runtime_asset_symlink_not_admitted")
    if executable and not os.access(resolved, os.X_OK):
        raise QualificationHeld("configured_python_not_executable")
    if sha256_file(resolved) != _digest(expected):
        raise ArtifactIntegrityError("Configured chemistry asset bytes changed")
    return {"path": str(resolved), "sha256": expected, "size_bytes": resolved.stat().st_size}


def _cas(sha256, fetch_artifact):
    _digest(sha256)
    if not callable(fetch_artifact):
        raise QualificationHeld("CAS_verification_callback_required")
    try:
        raw = fetch_artifact(sha256)
    except (KeyError, OSError) as error:
        raise QualificationHeld("qualification_CAS_artifact_missing") from error
    if not isinstance(raw, bytes) or len(raw) > MAX_BYTES or hashlib.sha256(raw).hexdigest() != sha256:
        raise ArtifactIntegrityError("Qualification CAS bytes do not match their identity")
    return raw


def _manifest(source, fetch_artifact):
    raw = decode_manifest(source, lambda sha: _cas(sha, fetch_artifact))
    return _json(raw, maximum=MAX_MANIFEST_BYTES)


def chemistry_asset_identity(provider, profile):
    """Verify exact operator assets and per-component code license evidence."""
    if provider not in DEPENDENCIES or not isinstance(profile, dict) or profile.get("provider") != provider:
        raise ValueError("Chemistry provider/profile mismatch")
    for field in ("version", "preprocessing_version"):
        if not isinstance(profile.get(field), str) or not 0 < len(profile[field]) <= 200:
            raise QualificationHeld("pinned_chemistry_versions_required")
    python = _asset(profile.get("python"), profile.get("python_sha256"), executable=True)
    lock = _asset(profile.get("runtime_lock"), profile.get("runtime_lock_sha256"))
    entries = profile.get("runtime_files")
    if not isinstance(entries, list) or not 2 <= len(entries) <= 100000:
        raise QualificationHeld("implementation_and_native_runtime_manifest_required")
    files, seen = [], set()
    roles = set()
    for entry in entries:
        if not isinstance(entry, dict) or entry.get("role") not in {"implementation", "native", "support"}:
            raise ValueError("Runtime manifest entry requires implementation/native/support role")
        item = _asset(entry.get("path"), entry.get("sha256"))
        if item["path"] in seen:
            raise ValueError("Duplicate runtime asset path")
        seen.add(item["path"])
        roles.add(entry["role"])
        files.append({**item, "role": entry["role"]})
    if not {"implementation", "native"} <= roles:
        raise QualificationHeld("implementation_and_native_runtime_manifest_required")
    expected_module = str(Path(__file__).with_name("algorithm_chemistry.py").resolve())
    if expected_module not in seen:
        raise QualificationHeld("chemistry_implementation_hash_required")
    versions = profile.get("package_versions")
    if not isinstance(versions, dict) or not DEPENDENCIES[provider] <= versions.keys():
        raise QualificationHeld("pinned_provider_and_RDKit_package_versions_required")
    if any(not isinstance(version, str) or not 0 < len(version) <= 100 for version in versions.values()):
        raise ValueError("Invalid pinned package version")
    pinned_lines = {line.lower().replace("_", "-") for line in Path(lock["path"]).read_text().splitlines()}
    for package, version in versions.items():
        if (package + "==" + version).lower().replace("_", "-") not in pinned_lines:
            raise QualificationHeld("runtime_lock_does_not_pin_declared_package")
    licenses, components = [], set()
    entries = profile.get("license_receipts")
    if not isinstance(entries, list) or not 1 <= len(entries) <= 100:
        raise QualificationHeld("commercial_dependency_license_receipts_required")
    for entry in entries:
        if not isinstance(entry, dict):
            raise ValueError("Invalid code license receipt entry")
        asset = _asset(entry.get("path"), entry.get("sha256"))
        if asset["size_bytes"] > 256 * 1024:
            raise ValueError("Code license receipt exceeds256KiB")
        record = _json(Path(asset["path"]).read_bytes())
        if (record.get("scope") != "code" or record.get("license") not in LICENSES
                or not isinstance(record.get("source_url"), str) or not record["source_url"].startswith("https://")
                or not isinstance(record.get("component"), str)):
            raise QualificationHeld("code_license_scope_or_source_unresolved")
        component = record["component"]
        if component in components:
            raise ValueError("Duplicate component license receipt")
        components.add(component)
        text = record.get("license_text")
        if not isinstance(text, dict):
            raise QualificationHeld("exact_code_license_text_required")
        license_text = _asset(text.get("path"), text.get("sha256"))
        if component in versions and record.get("version") != versions[component]:
            raise QualificationHeld("license_receipt_package_version_mismatch")
        licenses.append({"component": component, "license": record["license"],
                         "source_url": record["source_url"], "receipt": asset, "license_text": license_text})
    if not DEPENDENCIES[provider] <= components:
        raise QualificationHeld("separate_provider_and_dependency_code_receipts_required")
    value = {"schema": "skinpipe-chemistry-assets-v1", "provider": provider,
             "version": profile["version"], "preprocessing_version": profile["preprocessing_version"],
             "python": python, "runtime_lock": lock, "runtime_files": sorted(files, key=lambda item: item["path"]),
             "package_versions": dict(sorted(versions.items())), "license_receipts": sorted(licenses, key=lambda item: item["component"])}
    return {**value, "sha256": _hash(value)}


def _source(provider, input_sha256, source_receipt_sha256, fetch_artifact, *, maximum_rows=32):
    payload = _json(_cas(input_sha256, fetch_artifact))
    source = _json(_cas(source_receipt_sha256, fetch_artifact))
    if (source.get("schema") != SOURCE_SCHEMA or source.get("provider") != provider
            or source.get("input_sha256") != input_sha256 or payload.get("provider") != provider
            or source.get("campaign_id") != payload.get("campaign_id")
            or not isinstance(source.get("campaign_id"), str)
            or type(source.get("synthetic")) is not bool
            or source.get("source_type") not in {"frozen_campaign", "frozen_artifact_canary", "prepared_package_fixture", "synthetic"}):
        raise QualificationHeld("registered_source_receipt_not_bound_to_input")
    manifest = _manifest(source, fetch_artifact)
    rows = payload.get("rows")
    if not isinstance(rows, list) or not 1 <= len(rows) <= maximum_rows:
        raise QualificationHeld("bounded_real_input_pilot_required")
    ids = [row.get("candidate_id") for row in rows if isinstance(row, dict)]
    if (len(ids) != len(rows) or any(not isinstance(cid, str) or not 0 < len(cid) <= 200 for cid in ids)
            or len(set(ids)) != len(ids)):
        raise QualificationHeld("distinct_registered_pilot_candidates_required")
    context = payload.get("context", {})
    if not isinstance(context, dict):
        raise QualificationHeld("registered_input_context_required")
    target = context.get("target_id")
    if source.get("source_type") in {"frozen_campaign", "frozen_artifact_canary"}:
        if not isinstance(target, str) or not target or manifest.get("plan", {}).get("target_id") != target:
            raise QualificationHeld("registered_target_manifest_mismatch")
        registered_rows = manifest.get("candidates")
        if (not isinstance(registered_rows, list) or any(not isinstance(row, dict) or not isinstance(row.get("candidate_id"), str)
                                                       for row in registered_rows)):
            raise QualificationHeld("registered_candidate_manifest_required")
        registered = {row["candidate_id"]: row for row in registered_rows}
        if len(registered) != len(registered_rows):
            raise QualificationHeld("registered_manifest_candidate_ids_not_unique")
        for row in rows:
            original = registered.get(row["candidate_id"])
            if (not original or row.get("target_id") != target or original.get("target_id") != target
                    or original.get("smiles") != row.get("smiles")):
                raise QualificationHeld("pilot_candidate_not_in_frozen_target_manifest")
    for row in rows:
        identity = chemical_identity(row.get("smiles"))
        if identity["structure_status"] != "valid":
            raise QualificationHeld("real_input_molecular_identity_required")
        if row.get("exact_structure_id", identity["exact_structure_id"]) != identity["exact_structure_id"]:
            raise QualificationHeld("source_input_exact_identity_mismatch")
    artifacts = source.get("source_artifacts")
    if not isinstance(artifacts, list) or not 1 <= len(artifacts) <= 100:
        raise QualificationHeld("actual_source_artifact_bytes_required")
    verified_artifacts = []
    source_hashes = set()
    for entry in artifacts:
        if not isinstance(entry, dict) or not isinstance(entry.get("role"), str):
            raise ValueError("Source artifact requires role and SHA256")
        raw = _cas(entry.get("sha256"), fetch_artifact)
        if not raw:
            raise QualificationHeld("empty_source_artifact")
        source_hashes.add(entry["sha256"])
        verified_artifacts.append({"sha256": entry["sha256"], "role": entry["role"], "size_bytes": len(raw)})
    required = set(payload.get("context_artifacts", {}).values())
    for row in rows:
        required.update(row[key] for key in ("pose_sdf_sha256", "pose_sha256") if row.get(key))
    if not required <= source_hashes:
        raise QualificationHeld("pilot_uses_artifact_outside_registered_source")
    synthetic = source["synthetic"] or source["source_type"] in {"synthetic", "prepared_package_fixture"}
    return payload, {"sha256": source_receipt_sha256, "manifest_sha256": source["manifest_sha256"],
                     "source_type": source["source_type"], "synthetic": synthetic,
                     "source_artifacts": verified_artifacts}


def _observed(provider, identity, payload, output_sha256, fetch_artifact):
    output = _json(_cas(output_sha256, fetch_artifact))
    expected_ids = [row["candidate_id"] for row in payload["rows"]]
    rows = output.get("candidates")
    if (output.get("provider") != provider or output.get("status") != "complete" or not isinstance(rows, list)
            or any(not isinstance(row, dict) for row in rows)
            or [row.get("candidate_id") for row in rows] != expected_ids
            or any(row.get("status") != "complete" for row in rows)):
        raise QualificationHeld("successful_exact_candidate_output_required")
    for source, row in zip(payload["rows"], rows):
        expected = chemical_identity(source.get("smiles"))["exact_structure_id"]
        if row.get("input_exact_structure_id") != expected:
            raise QualificationHeld("qualification_output_molecular_identity_mismatch")
        if not isinstance(row.get("metrics"), dict) or row["metrics"].get("method") != METHODS[provider]:
            raise QualificationHeld("actual_provider_output_method_required")
    measurement = output.get("measurements", {})
    if not isinstance(measurement, dict):
        raise QualificationHeld("actual_time_and_memory_measurements_required")
    seconds = measurement.get("wall_seconds")
    memory = measurement.get("process_lifetime_peak_rss_bytes")
    if (isinstance(seconds, bool) or not isinstance(seconds, (int, float)) or not math.isfinite(seconds) or seconds <= 0
            or isinstance(memory, bool) or not isinstance(memory, (int, float)) or not math.isfinite(memory) or memory <= 0):
        raise QualificationHeld("actual_time_and_memory_measurements_required")
    if measurement.get("memory_scope") != "process lifetime high-water mark; not per-job incremental peak":
        raise QualificationHeld("chemistry_memory_measurement_scope_required")
    observed = output.get("software", {})
    for component, version in identity["package_versions"].items():
        if observed.get("packages", {}).get(component, {}).get("version") != version:
            raise QualificationHeld("observed_package_version_differs_from_assets")
    runtime = {entry["path"]: entry["sha256"] for entry in identity["runtime_files"]}
    module_hashes = observed.get("module_sha256", {})
    if not isinstance(module_hashes, dict) or not module_hashes:
        raise QualificationHeld("actual_imported_module_hashes_required")
    for entry in module_hashes.values():
        if not isinstance(entry, dict) or runtime.get(str(Path(entry.get("path", "")).resolve())) != entry.get("sha256"):
            raise QualificationHeld("observed_imported_module_not_in_asset_manifest")
    artifacts = [*output.get("artifacts", []), *(artifact for row in rows for artifact in row.get("artifacts", []))]
    if not artifacts:
        raise QualificationHeld("actual_output_artifact_bytes_required")
    for artifact in artifacts:
        if not isinstance(artifact, dict):
            raise ArtifactIntegrityError("Invalid declared chemistry output artifact")
        raw = _cas(artifact.get("sha256"), fetch_artifact)
        if artifact.get("size_bytes", len(raw)) != len(raw):
            raise ArtifactIntegrityError("Declared chemistry output size differs from CAS bytes")
    return output, {"elapsed_seconds": float(seconds), "peak_memory_mib": float(memory) / 1024**2,
                    "memory_scope": measurement.get("memory_scope"), "completed_candidates": len(rows)}


def create_chemistry_qualification(provider, profile, *, input_sha256, output_sha256,
                                    fetch_artifact, destination, source_receipt_sha256):
    """Create only from verified CAS inputs/results and registered source bytes."""
    identity = chemistry_asset_identity(provider, profile)
    payload, source = _source(provider, input_sha256, source_receipt_sha256, fetch_artifact)
    _, measured = _observed(provider, identity, payload, output_sha256, fetch_artifact)
    receipt = {"schema": SCHEMA, "status": "software_only" if source["synthetic"] else "qualified",
               "provider": provider, "asset_identity_sha256": identity["sha256"],
               "input_sha256": input_sha256, "output_sha256": output_sha256,
               "source_receipt_sha256": source_receipt_sha256, "source": source, **measured,
               "observed_at_unix": time.time(), "scientific_validation": "NOT_ESTABLISHED",
               "scientific_promotion_allowed": False}
    path = Path(destination).resolve()
    if path.exists():
        raise ValueError("Qualification receipts are immutable; use a fresh destination")
    atomic_bytes(path, canonical_json(receipt))
    return {"path": str(path), "sha256": sha256_file(path), "size_bytes": path.stat().st_size,
            "status": receipt["status"], "asset_identity_sha256": identity["sha256"]}


def verify_chemistry_qualification(provider, profile, *, fetch_artifact, qualification_mode=False,
                                   input_sha256=None, source_receipt_sha256=None):
    """Validate fresh assets/licenses in both modes; ordinary runs require proof."""
    if type(qualification_mode) is not bool:
        raise ValueError("Qualification mode must be boolean")
    identity = chemistry_asset_identity(provider, profile)
    if qualification_mode:
        if input_sha256 is None or source_receipt_sha256 is None:
            raise QualificationHeld("registered_qualification_pilot_source_required")
        _, source = _source(provider, input_sha256, source_receipt_sha256, fetch_artifact)
        return {**identity, "qualification_status": "pilot_source_verified", "source": source}
    ref = profile.get("qualification_receipt")
    if not isinstance(ref, dict):
        raise QualificationHeld("real_input_runtime_qualification_required")
    asset = _asset(ref.get("path"), ref.get("sha256"))
    if asset["size_bytes"] > 256 * 1024:
        raise ValueError("Qualification receipt exceeds256KiB")
    receipt = _json(Path(asset["path"]).read_bytes())
    if (receipt.get("schema") != SCHEMA or receipt.get("status") != "qualified"
            or receipt.get("provider") != provider or receipt.get("asset_identity_sha256") != identity["sha256"]
            or receipt.get("scientific_validation") != "NOT_ESTABLISHED"
            or receipt.get("scientific_promotion_allowed") is not False):
        raise QualificationHeld("runtime_qualification_not_bound_or_scientific_claim_invalid")
    payload, source = _source(provider, receipt.get("input_sha256"), receipt.get("source_receipt_sha256"), fetch_artifact)
    _, measurement = _observed(provider, identity, payload, receipt.get("output_sha256"), fetch_artifact)
    if source["synthetic"] or receipt.get("source") != source or any(receipt.get(key) != value for key, value in measurement.items()):
        raise QualificationHeld("qualification_source_or_measurements_changed")
    if input_sha256 is not None or source_receipt_sha256 is not None:
        _, current_source = _source(provider, input_sha256, source_receipt_sha256, fetch_artifact, maximum_rows=10000)
        if current_source["synthetic"]:
            raise QualificationHeld("synthetic_current_inputs_require_qualification_mode")
    return {**identity, "qualification_status": "qualified", "qualification_receipt": asset}


def verify_model_qualification(provider, profile, *, fetch_artifact, qualification_mode=False,
                               input_sha256=None, source_receipt_sha256=None):
    """Verify model pilot provenance or a real CAS-bound measured runtime receipt.

    Model command preparation separately enforces code/weight license evidence.
    GPU qualification requires observed process VRAM; RSS cannot substitute for
    it. Missing weights or receipts hold the branch instead of enabling a model.
    """
    from .algorithm_models import ModelHeld, model_asset_identity

    if type(qualification_mode) is not bool:
        raise ValueError("Qualification mode must be boolean")
    if not isinstance(profile, dict) or not profile:
        raise QualificationHeld("operator_model_runtime_profile_missing")
    try:
        identity = model_asset_identity(provider, profile)
    except ModelHeld as error:
        raise QualificationHeld(str(error)) from error
    if qualification_mode:
        if input_sha256 is None or source_receipt_sha256 is None:
            raise QualificationHeld("registered_model_pilot_source_required")
        _, source = _source(provider, input_sha256, source_receipt_sha256, fetch_artifact)
        return {**identity, "qualification_status": "pilot_source_verified", "source": source}
    ref = profile.get("qualification_receipt")
    if not isinstance(ref, dict):
        raise QualificationHeld("real_model_runtime_qualification_required")
    asset = _asset(ref.get("path"), ref.get("sha256"))
    if asset["size_bytes"] > 256 * 1024:
        raise ValueError("Model qualification receipt exceeds256KiB")
    receipt = _json(Path(asset["path"]).read_bytes())
    if (receipt.get("schema") != MODEL_SCHEMA or receipt.get("status") != "qualified" or receipt.get("provider") != provider
            or receipt.get("asset_identity_sha256") != identity["sha256"]
            or receipt.get("scientific_validation") != "NOT_ESTABLISHED"
            or receipt.get("scientific_promotion_allowed") is not False):
        raise QualificationHeld("model_runtime_qualification_not_bound_or_scientific_claim_invalid")
    payload, source = _source(provider, receipt.get("input_sha256"), receipt.get("source_receipt_sha256"), fetch_artifact)
    output = _json(_cas(receipt.get("output_sha256"), fetch_artifact))
    expected_ids = [row["candidate_id"] for row in payload["rows"]]
    rows = output.get("candidates")
    if (source["synthetic"] or output.get("provider") != provider or output.get("status") != "complete"
            or output.get("asset_identity_sha256") != identity["sha256"] or not isinstance(rows, list)
            or any(not isinstance(row, dict) for row in rows)
            or [row.get("candidate_id") for row in rows] != expected_ids
            or any(row.get("status") != "complete" or row.get("asset_identity_sha256") != identity["sha256"] for row in rows)):
        raise QualificationHeld("successful_model_output_not_bound_to_real_assets_and_candidates")
    if receipt.get("source") != source:
        raise QualificationHeld("model_source_receipt_differs_from_CAS_source")
    for source_row, row in zip(payload["rows"], rows):
        if row.get("input_exact_structure_id") != chemical_identity(source_row["smiles"])["exact_structure_id"]:
            raise QualificationHeld("model_output_molecular_identity_mismatch")
    seconds = output.get("compute_seconds")
    if isinstance(seconds, bool) or not isinstance(seconds, (int, float)) or not math.isfinite(seconds) or seconds <= 0:
        raise QualificationHeld("actual_model_compute_measurement_required")
    if receipt.get("elapsed_seconds") != seconds:
        raise QualificationHeld("model_time_measurement_differs_from_CAS_output")
    executions = output.get("executions", [])
    if not isinstance(executions, list) or len(executions) != len(rows) or any(not isinstance(entry, dict) or entry.get("returncode") != 0
                                                                           for entry in executions):
        raise QualificationHeld("successful_model_execution_measurements_required")
    for entry in executions:
        elapsed = entry.get("elapsed_seconds")
        if isinstance(elapsed, bool) or not isinstance(elapsed, (int, float)) or not math.isfinite(elapsed) or elapsed <= 0:
            raise QualificationHeld("actual_model_execution_elapsed_measurements_required")
    peaks = [entry.get("peak_mib") for entry in executions]
    device = receipt.get("runtime_device")
    if device not in {"cpu", "gpu"} or output.get("runtime_device") != device:
        raise QualificationHeld("actual_model_runtime_device_required")
    if provider in {"boltzina", "boltz2", "chai1", "unidock2"} and device != "gpu":
        raise QualificationHeld("GPU_model_runtime_required")
    gpu_provider = device == "gpu"
    if gpu_provider:
        if (not peaks or any(isinstance(peak, bool) or not isinstance(peak, (int, float))
                            or not math.isfinite(peak) or peak <= 0 for peak in peaks)
                or receipt.get("peak_memory_mib") != max(peaks)
                or receipt.get("memory_scope") != "gpu_process_peak"):
            raise QualificationHeld("actual_GPU_process_VRAM_measurements_required")
    else:
        measurement = output.get("measurements", {})
        memory = measurement.get("process_peak_rss_bytes")
        if (isinstance(memory, bool) or not isinstance(memory, (int, float)) or not math.isfinite(memory) or memory <= 0
                or receipt.get("peak_memory_mib") != memory / 1024**2 or receipt.get("memory_scope") != "cpu_process_peak_rss"):
            raise QualificationHeld("actual_CPU_process_RSS_measurements_required")
    artifacts = [artifact for row in rows for artifact in row.get("artifacts", [])]
    if not artifacts:
        raise QualificationHeld("actual_model_output_artifact_bytes_required")
    for artifact in artifacts:
        if not isinstance(artifact, dict):
            raise ArtifactIntegrityError("Invalid declared model output artifact")
        raw = _cas(artifact.get("sha256"), fetch_artifact)
        if artifact.get("size_bytes", len(raw)) != len(raw):
            raise ArtifactIntegrityError("Model output artifact size differs from CAS bytes")
    if input_sha256 is not None or source_receipt_sha256 is not None:
        _, current_source = _source(provider, input_sha256, source_receipt_sha256, fetch_artifact, maximum_rows=10000)
        if current_source["synthetic"]:
            raise QualificationHeld("synthetic_current_inputs_require_qualification_mode")
    return {**identity, "qualification_status": "qualified", "qualification_receipt": asset}
