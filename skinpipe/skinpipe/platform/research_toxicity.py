"""Optional pinned local ADMET-AI inference, separate from published DB values.

Only operator configuration selects runtime/weights. This fixed runner neither
downloads models nor invokes a client-provided program. All outputs remain
auxiliary model predictions with unestablished applicability/biological validity.
"""
from __future__ import annotations

from copy import deepcopy
import csv
import hashlib
import importlib.metadata
import importlib.util
import io
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import time

from .artifacts import atomic_bytes, canonical_json, sha256_file
from .identity import chemical_identity, digest

POLICY = "local-admet-ai-exact-input-cpu-v1"
PACKAGES = {"admet-ai", "torch", "chemprop", "rdkit", "lightning", "numpy"}
RUNNER_FILES = ("__init__.py", "platform/__init__.py", "platform/artifacts.py",
                "platform/identity.py", "platform/research_toxicity.py")


def _finite_number(value):
    """A finite JSON number; oversized integers are invalid, not a crash.

    The runner and the cached-row scan both read JSON that may carry integers
    beyond the float range, where ``math.isfinite`` raises ``OverflowError``.
    Such a value must be treated as unusable output, like ``NaN``.
    """
    if type(value) not in (int, float):
        return False
    try:
        return math.isfinite(value)
    except OverflowError:
        return False


def _runner_sources():
    package = Path(__file__).resolve().parents[1]
    return {name: sha256_file(package / name) for name in RUNNER_FILES}


def _stage_runner(work, identity):
    """Expose only frozen SkinPipe modules, never the server's dependencies."""
    expected = identity.get("runner_sources_sha256")
    if not isinstance(expected, dict) or set(expected) != set(RUNNER_FILES):
        raise ValueError("Pinned isolated runner source inventory required")
    if expected["platform/research_toxicity.py"] != identity.get("runner_sha256"):
        raise ValueError("Isolated runner identity differs")
    source = Path(__file__).resolve().parents[1]
    prefix = work / "runner-code"
    prefix.mkdir(mode=0o700, exist_ok=False)
    for name in RUNNER_FILES:
        raw = (source / name).read_bytes()
        if len(raw) > 4 * 1024**2 or hashlib.sha256(raw).hexdigest() != expected[name]:
            raise ValueError("Pinned isolated runner source changed")
        target = prefix / "skinpipe" / name
        atomic_bytes(target, raw)
        if sha256_file(target) != expected[name] or sha256_file(source / name) != expected[name]:
            raise ValueError("Isolated runner source copy differs")
    atomic_bytes(prefix / "source-manifest.json", canonical_json(expected))
    return prefix


def _file(entry, maximum):
    if not isinstance(entry, dict) or set(entry) != {"path", "sha256"}:
        raise ValueError("Pinned file paths and SHA256 required")
    path = Path(entry["path"])
    if (not path.is_absolute() or not path.is_file() or path.stat().st_size > maximum
            or not re.fullmatch(r"[0-9a-f]{64}", entry["sha256"]) or sha256_file(path) != entry["sha256"]):
        raise ValueError("Pinned model/runtime file unavailable or changed")
    return path.resolve()


def _profile(profile):
    allowed = {"provider", "version", "python", "implementation", "models_dir", "weights", "model_manifest",
               "endpoint_definitions", "packages", "preprocessing", "batch_size", "max_candidates", "budget_seconds"}
    if (not isinstance(profile, dict) or set(profile) - allowed or profile.get("provider") != "admet_ai"
            or profile.get("version") != "2.0.1" or profile.get("preprocessing") != "exact-rdkit-graph-v1"):
        raise ValueError("Only pinned ADMET-AI 2.0.1 exact-graph profiles are supported")
    profile = deepcopy(profile)
    for key, default, maximum in (("batch_size", 16, 32), ("max_candidates", 100, 200), ("budget_seconds", 120, 600)):
        value = profile.setdefault(key, default)
        if isinstance(value, bool) or not isinstance(value, int) or not 1 <= value <= maximum:
            raise ValueError("Invalid bounded model execution policy")
    if not isinstance(profile.get("packages"), dict) or set(profile["packages"]) != PACKAGES:
        raise ValueError("Exact model dependency versions required")
    python = _file(profile["python"], 256 * 1024**2)
    if not os.access(python, os.X_OK):
        raise ValueError("Pinned runtime must be executable")
    implementation = _file(profile["implementation"], 8 * 1024**2)
    manifest_path = _file(profile["model_manifest"], 4 * 1024**2)
    endpoints_path = _file(profile["endpoint_definitions"], 4 * 1024**2)
    manifest = json.loads(manifest_path.read_bytes())
    if (manifest.get("version") != profile["version"] or manifest.get("implementation_sha256") != profile["implementation"]["sha256"]
            or manifest.get("endpoint_metadata_sha256") != profile["endpoint_definitions"]["sha256"]
            or profile["packages"]["admet-ai"] != profile["version"]
            or any(profile["packages"].get(name) != version for name, version in manifest.get("packages", {}).items())):
        raise ValueError("Model manifest/runtime/endpoint identity differs")
    root = Path(profile["models_dir"])
    if not root.is_absolute() or not root.is_dir():
        raise ValueError("Installed model directory required")
    root = root.resolve()
    if not isinstance(profile.get("weights"), list) or len(profile["weights"]) != 10:
        raise ValueError("ADMET-AI requires ten pinned ensemble checkpoints")
    weights = {}
    for entry in profile["weights"]:
        path = _file(entry, 1024**3)
        if not path.is_relative_to(root):
            raise ValueError("Checkpoint outside selected model directory")
        weights[path.relative_to(root).as_posix()] = entry["sha256"]
    expected = {"/".join(Path(entry["path"]).parts[-2:]): entry["sha256"] for entry in manifest.get("weights", [])}
    if (weights != expected or len(weights) != 10
            or {p.resolve() for p in root.rglob("*.pt")} != {root / name for name in weights}):
        raise ValueError("Installed checkpoint population differs from model manifest")
    endpoints = list(csv.DictReader(io.StringIO(endpoints_path.read_text())))
    ids = [entry.get("id") for entry in endpoints]
    if not ids or len(ids) > 200 or len(set(ids)) != len(ids) or any(not isinstance(name, str) or not name for name in ids):
        raise ValueError("Bounded unique endpoint definitions required")
    identity = {"policy": POLICY, "provider": "admet_ai", "version": profile["version"],
                "model_id": manifest["model_id"], "model_manifest_sha256": profile["model_manifest"]["sha256"],
                "endpoint_definitions_sha256": profile["endpoint_definitions"]["sha256"],
                "implementation_sha256": profile["implementation"]["sha256"], "python_sha256": profile["python"]["sha256"],
                "runner_sha256": sha256_file(__file__), "weights": weights, "packages": profile["packages"],
                "runner_sources_sha256": _runner_sources(),
                "preprocessing": profile["preprocessing"], "batch_size": profile["batch_size"],
                "max_candidates": profile["max_candidates"], "budget_seconds": profile["budget_seconds"],
                "device": "cpu", "threads": 1, "drugbank_enabled": False, "percentiles_enabled": False,
                "applicability_domain": "NOT_ESTABLISHED", "scientific_validation": "NOT_ESTABLISHED"}
    return profile, identity, manifest, endpoints, implementation


def freeze_toxicity(config, root):
    """Freeze trusted paths privately; export only the opaque ID and file hashes."""
    selected = config.get("toxicity_profile")
    if selected is None:
        return {"policy": POLICY, "status": "not_configured", "reason": "additional_model_not_configured"}
    try:
        path = Path(selected)
        if not path.is_absolute() or not path.is_file() or path.stat().st_size > 1024**2:
            raise ValueError("Private operator profile unavailable")
        profile, identity, _, _, _ = _profile(json.loads(path.read_bytes()))
        record = {"profile": profile, "identity": identity}
        binding = digest(record)
        private = Path(root) / "toxicity-profile-bindings"
        private.mkdir(parents=True, exist_ok=True, mode=0o700)
        private.chmod(0o700)
        target = private / (binding + ".json")
        if not target.exists():
            atomic_bytes(target, canonical_json(record))
            target.chmod(0o600)
        return {"policy": POLICY, "status": "configured", "binding_id": binding, "identity": identity}
    except (ValueError, KeyError, TypeError, OSError) as error:
        return {"policy": POLICY, "status": "held", "reason": "model_profile_unavailable_or_unsupported", "error_type": type(error).__name__}


def validate_toxicity_profile(path):
    path = Path(path)
    if not path.is_absolute() or not path.is_file() or path.stat().st_size > 1024**2:
        raise ValueError("Private operator profile unavailable")
    profile, _, _, _, _ = _profile(json.loads(path.read_bytes()))
    return profile


def toxicity_profile_status(path=None):
    if path is None:
        return {"state": "unconfigured", "model_version": None, "message": "Additional local model is not configured"}
    try:
        profile = validate_toxicity_profile(path)
        return {"state": "available", "model_version": profile["version"], "message": "Pinned optional CPU model profile is available; runtime qualification required"}
    except (ValueError, KeyError, TypeError, OSError):
        return {"state": "invalid", "model_version": None, "message": "Additional local model profile is unavailable or changed"}


def research_toxicity_capability(config):
    """Read operator readiness without creating bindings or revealing paths.

    A valid pinned profile is installation readiness only. Actual inference
    remains recorded independently for each frozen research request.
    """
    result = {"state": "unconfigured", "provider": "admet_ai", "model_version": None,
              "policy": "missing_exact_db_predictions_only", "max_candidates": None,
              "budget_seconds": None, "scientific_validation": "NOT_ESTABLISHED",
              "runtime_qualification": "NOT_ESTABLISHED"}
    selected = config.get("toxicity_profile") if isinstance(config, dict) else None
    if selected is None:
        return result
    try:
        profile = validate_toxicity_profile(selected)
        result.update(state="available", model_version=profile["version"],
                      max_candidates=profile["max_candidates"], budget_seconds=profile["budget_seconds"])
    except (ValueError, KeyError, TypeError, OSError, AttributeError, csv.Error):
        result["state"] = "held"
    return result


def _run(payload):
    """Fixed numerical adapter. Called inside the pinned, offline CPU runtime."""
    profile, identity, _, endpoints, implementation = _profile(payload["profile"])
    if identity != payload["identity"]:
        raise ValueError("Model profile changed since submission")
    for name, expected in profile["packages"].items():
        if importlib.metadata.version(name) != expected:
            raise ValueError("Model dependency version changed: " + name)
    spec = importlib.util.find_spec("admet_ai.admet_model")
    if spec is None or Path(spec.origin).resolve() != implementation:
        raise ValueError("Imported model implementation differs from the pinned profile")
    import torch
    torch.set_num_threads(1)
    torch.set_num_interop_threads(1)
    if torch.cuda.is_available():
        raise ValueError("Additional toxicity runner must remain CPU-only")
    from admet_ai import ADMETModel
    model = ADMETModel(models_dir=profile["models_dir"], include_physchem=True, drugbank_path=None, num_workers=0)
    if model.device != "cpu" or model.drugbank is not None or model.drugbank_atc_filtered is not None:
        raise ValueError("Model selected an unpermitted device or external reference dataset")
    ids = {entry["id"] for entry in endpoints}
    outputs = []
    def predict(items):
        try:
            for item in items:
                if chemical_identity(item["smiles"])["exact_structure_id"] != item["exact_structure_id"]:
                    raise ValueError("Runtime exact graph differs from the frozen input")
            frame = model.predict(smiles=[item["smiles"] for item in items])
            if (list(frame.index) != [item["smiles"] for item in items] or len(frame) != len(items)
                    or len(set(frame.columns)) != len(frame.columns) or set(frame.columns) != ids):
                raise ValueError("Model input index/cardinality/endpoint schema differs")
            batch_outputs = []
            for item, values in zip(items, frame.to_dict("records")):
                raw = {key: str(value) for key, value in values.items()}
                clean = {key: float(value) if math.isfinite(float(value)) else None for key, value in values.items()}
                invalid = []
                for endpoint in endpoints:
                    value = clean[endpoint["id"]]
                    if endpoint.get("task_type") == "classification" and value is not None and not 0 <= value <= 1:
                        invalid.append(endpoint["id"])
                computed = {entry["id"] for entry in endpoints if entry.get("category") == "Physicochemical"}
                batch_outputs.append({**item, "status": "failed_invalid_probability" if invalid else "success" if all(v is not None for v in clean.values()) else "partial_nonfinite",
                                      "values": {k: v for k, v in clean.items() if k not in computed},
                                      "physicochemical_values": {k: v for k, v in clean.items() if k in computed},
                                      "raw_values": raw, "quality_flags": {"invalid_classifier_probability": invalid}})
            outputs.extend(batch_outputs)
        except Exception as error:
            if len(items) > 1:
                mid = len(items) // 2
                predict(items[:mid])
                predict(items[mid:])
            else:
                outputs.append({**items[0], "status": "failed", "values": {}, "error_type": type(error).__name__})
    for offset in range(0, len(payload["items"]), profile["batch_size"]):
        predict(payload["items"][offset:offset + profile["batch_size"]])
    return {"status": "complete", "identity": identity, "input_sha256": digest(payload["items"]), "outputs": outputs}


def _admit_result(result, selected, identity, endpoints):
    if (result.get("status") != "complete" or result.get("identity") != identity
            or result.get("input_sha256") != digest(selected) or len(result.get("outputs", [])) != len(selected)):
        raise ValueError("Additional model output receipt differs from exact inputs")
    for expected, output in zip(selected, result["outputs"]):
        if any(output.get(key) != value for key, value in expected.items()):
            raise ValueError("Additional model output reordered or changed an input")
        status = output.get("status")
        if status not in {"success", "partial_nonfinite", "failed", "failed_invalid_probability"}:
            raise ValueError("Unsupported additional model output status")
        computed = {entry["id"] for entry in endpoints if entry.get("category") == "Physicochemical"}
        learned = {entry["id"] for entry in endpoints} - computed
        for key, expected_keys in (("values", learned), ("physicochemical_values", computed)):
            values = output.get(key, {})
            if (not isinstance(values, dict) or set(values) - expected_keys
                    or status in {"success", "partial_nonfinite"} and set(values) != expected_keys
                    or any(value is not None and not _finite_number(value) for value in values.values())):
                raise ValueError("Additional model endpoint schema/nonfinite output differs")
            if status == "success" and any(value is None for value in values.values()):
                raise ValueError("Successful additional model result contains missing endpoints")
        if status in {"success", "partial_nonfinite"}:
            for endpoint in endpoints:
                if endpoint.get("task_type") == "classification":
                    value = output.get("values", {}).get(endpoint["id"])
                    if value is not None and not 0 <= value <= 1:
                        raise ValueError("Successful additional classifier probability outside 0..1")


def _execute(payload, work, python, *, stop, budget):
    from .engines import terminate_owned
    stop()
    request_path, result_path = work / "request.json", work / "result.json"
    completion_path = work / "result.complete.json"
    if result_path.exists() or completion_path.exists():
        raise ValueError("Model execution requires a fresh owned attempt directory")
    runner = _stage_runner(work, payload["identity"])
    atomic_bytes(request_path, canonical_json(payload))
    environment = {**{key: os.environ[key] for key in ("PATH", "HOME", "TMPDIR", "LANG", "SYSTEMROOT", "WINDIR") if key in os.environ},
                   "PYTHONPATH": str(runner), "PYTHONNOUSERSITE": "1", "PYTHONSAFEPATH": "1", "CUDA_VISIBLE_DEVICES": "",
                   "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1",
                   "HF_HUB_OFFLINE": "1", "TRANSFORMERS_OFFLINE": "1"}
    started = time.monotonic()
    durable_result = False
    result = None
    endpoint_path = _file(payload["profile"]["endpoint_definitions"], 4 * 1024**2)
    endpoints = list(csv.DictReader(io.StringIO(endpoint_path.read_text())))
    with (work / "execution.log").open("wb") as log:
        process = subprocess.Popen([str(python), "-m", "skinpipe.platform.research_toxicity", str(request_path), str(result_path)],
                                   env=environment, stdout=log, stderr=log, start_new_session=True)
        try:
            while process.poll() is None:
                stop()
                if time.monotonic() - started >= budget:
                    raise RuntimeError("additional_model_stage_budget_exhausted")
                # The fixed module publishes one fsynced atomic result only
                # after every input has a result. Some scientific runtimes
                # retain interpreter-exit threads; reap that owned session.
                if completion_path.is_file():
                    if result_path.stat().st_size > 32 * 1024**2:
                        raise ValueError("Additional model response exceeds its bound")
                    if completion_path.stat().st_size > 4096:
                        raise ValueError("Additional model completion receipt exceeds its bound")
                    completion = json.loads(completion_path.read_bytes())
                    if completion != {"result_sha256": sha256_file(result_path), "input_sha256": digest(payload["items"]),
                                      "model_identity_sha256": digest(payload["identity"])}:
                        raise ValueError("Additional model durable completion receipt differs")
                    result = json.loads(result_path.read_bytes())
                    _admit_result(result, payload["items"], payload["identity"], endpoints)
                    durable_result = True
                    break
                time.sleep(0.05)
            stop()
            if time.monotonic() - started >= budget:
                raise RuntimeError("additional_model_stage_budget_exhausted")
            if result is None and completion_path.is_file() and result_path.is_file() and result_path.stat().st_size <= 32 * 1024**2:
                completion = json.loads(completion_path.read_bytes())
                if completion != {"result_sha256": sha256_file(result_path), "input_sha256": digest(payload["items"]),
                                  "model_identity_sha256": digest(payload["identity"])}:
                    raise ValueError("Additional model durable completion receipt differs")
                result = json.loads(result_path.read_bytes())
                _admit_result(result, payload["items"], payload["identity"], endpoints)
        finally:
            terminate_owned(process)
    if (process.returncode != 0 and not durable_result) or not result_path.is_file() or result_path.stat().st_size > 32 * 1024**2:
        raise RuntimeError("additional_model_runtime_held")
    if result is None:
        raise RuntimeError("additional_model_runtime_held")
    result["runtime_completion"] = {"returncode": process.returncode, "durable_result_observed": durable_result,
                                    "owned_process_cleanup": True}
    return result


def apply_toxicity(rows, service, request, *, stop, progress):
    """Add a separate channel only; cached predictions and priorities stay intact."""
    policy = request.get("toxicity_policy") or {"policy": POLICY, "status": "not_configured"}
    receipt = {"policy": deepcopy(policy), "status": policy["status"], "eligible_count": 0, "selected_count": 0,
               "completed_count": 0, "failed_count": 0, "partial_count": 0, "skipped_cached_count": 0,
               "unsupported_count": 0, "scientific_validation": "NOT_ESTABLISHED",
               "scope": "Separate auxiliary model channel; cached values, measured activity and ranking are not overwritten"}
    eligible = []
    for index, row in enumerate(rows):
        if row.get("additional_toxicity") is not None:
            row.setdefault("source_additional_toxicity", []).append(deepcopy(row["additional_toxicity"]))
        cached = row.get("admet") or []
        cached = [cached] if isinstance(cached, dict) else cached
        # A partial cache row can contain only null endpoints. That records a
        # failed measurement of this channel, not a reusable prediction. Keep
        # a genuine zero as observed, and leave the source cache untouched.
        usable = any(p.get("status") in {"success", "complete", "partial_nonfinite"}
            and isinstance(p.get("values"), dict)
            and any(_finite_number(value) for value in p["values"].values())
            for p in cached if isinstance(p, dict))
        row["additional_toxicity"] = {"status": "not_requested_cached_prediction" if usable else "unassessed_model_not_configured",
                                      "scientific_validation": "NOT_ESTABLISHED", "values": {}}
        if usable:
            row["additional_toxicity"]["cached_prediction_statuses"] = sorted({p.get("status", "unknown") for p in cached if isinstance(p, dict)})
            receipt["skipped_cached_count"] += 1
            continue
        identity = chemical_identity(row.get("smiles"))
        if identity["structure_status"] != "valid" or identity["parent_ambiguous"]:
            row["additional_toxicity"]["status"] = "unsupported_structure"
            receipt["unsupported_count"] += 1
            continue
        eligible.append({"row_index": index, "candidate_id": row["candidate_id"], "smiles": row["smiles"],
                         "exact_structure_id": identity["exact_structure_id"]})
    receipt["eligible_count"] = len(eligible)
    if policy["status"] != "configured" or not eligible:
        if policy["status"] == "held":
            for item in eligible:
                rows[item["row_index"]]["additional_toxicity"].update(status="unassessed_model_profile_held",
                    reason=policy.get("reason", "additional_model_profile_held"))
        return rows, receipt
    selected = []
    work = None
    try:
        binding = policy["binding_id"]
        if not re.fullmatch(r"[0-9a-f]{64}", binding):
            raise ValueError("Invalid private model binding")
        record = json.loads((service.root / "toxicity-profile-bindings" / (binding + ".json")).read_bytes())
        if digest(record) != binding:
            raise ValueError("Private model binding changed")
        profile, identity, manifest, endpoints, _ = _profile(record["profile"])
        if identity != policy["identity"]:
            raise ValueError("Model identity changed since submission")
        selected = eligible[:profile["max_candidates"]]
        receipt["selected_count"] = len(selected)
        for item in eligible[len(selected):]:
            rows[item["row_index"]]["additional_toxicity"]["status"] = "unassessed_model_candidate_limit"
        stop()
        progress({"stage": "additional_toxicity", "selected": len(selected)})
        work = service.root / "additional-toxicity" / digest([request["discovery_run_id"], binding]) / ("attempt-" + str(time.time_ns()))
        work.mkdir(parents=True, exist_ok=True)
        receipt["execution_id"] = work.name
        payload = {"profile": profile, "identity": identity, "items": selected}
        result = _execute(payload, work, profile["python"]["path"], stop=stop, budget=profile["budget_seconds"])
        stop()
        # Preserve a quarantined child response before interpreting it. The
        # parent never promotes malformed/raw nonfinite values into results.
        raw = json.dumps(result, ensure_ascii=False, sort_keys=True, separators=(",", ":"), allow_nan=True).encode()
        if len(raw) > 32 * 1024**2:
            raise ValueError("Additional raw model response exceeds its bound")
        raw_sha = service.artifacts.put_bytes(raw, media_type="application/json")["sha256"]
        receipt.update(raw_predictions_sha256=raw_sha, raw_predictions_admission="quarantined")
        _admit_result(result, selected, identity, endpoints)
        input_sha = service.artifacts.put_bytes(canonical_json(selected), media_type="application/json")["sha256"]
        manifest_sha = service.artifacts.put_bytes(canonical_json(manifest), media_type="application/json")["sha256"]
        endpoints_sha = service.artifacts.put_bytes(canonical_json(endpoints), media_type="application/json")["sha256"]
        for item, output in zip(selected, result["outputs"]):
            admitted = deepcopy(output)
            if output["status"].startswith("failed"):
                admitted["rejected_values"] = admitted.get("values", {})
                admitted["values"] = {}
            rows[item["row_index"]]["additional_toxicity"] = {**admitted, "model_identity": identity,
                "endpoint_definitions": [{key: entry.get(key) if isinstance(entry.get(key), str) and len(entry[key]) <= 512
                                          and not any(ord(c) < 32 for c in entry[key]) else None
                                          for key in ("id", "task_type", "category", "name", "units")} for entry in endpoints],
                "raw_predictions_sha256": raw_sha, "exact_input_receipt_sha256": input_sha,
                "model_manifest_sha256": manifest_sha, "endpoint_definitions_sha256": endpoints_sha,
                "source": "additional_local_cpu_prediction", "activity_inherited": False,
                "scientific_validation": "NOT_ESTABLISHED", "applicability_domain": "NOT_ESTABLISHED",
                "stereochemical_sensitivity": "NOT_ESTABLISHED"}
        receipt.update(status="complete", raw_predictions_sha256=raw_sha, raw_predictions_admission="admitted", exact_input_receipt_sha256=input_sha,
                       model_manifest_sha256=manifest_sha, endpoint_definitions_sha256=endpoints_sha,
                       completed_count=sum(r["status"] == "success" for r in result["outputs"]),
                       failed_count=sum(r["status"].startswith("failed") for r in result["outputs"]),
                       partial_count=sum(r["status"] == "partial_nonfinite" for r in result["outputs"]))
    except TimeoutError:
        raise
    except (ValueError, KeyError, TypeError, OSError, RuntimeError) as error:
        # Lease/budget/shutdown fences must reach the coordinator unchanged.
        from .engines import OwnedProcessCleanupError, ResourceDeferred
        from .queue import LeaseLost, Cancelled
        if isinstance(error, (OwnedProcessCleanupError, ResourceDeferred, LeaseLost, Cancelled)):
            raise
        if work is not None:
            raw_path = work / "result.json"
            if raw_path.is_file() and raw_path.stat().st_size <= 32 * 1024**2:
                receipt.update(raw_predictions_sha256=service.artifacts.put_file(raw_path, media_type="application/json")["sha256"],
                               raw_predictions_admission="quarantined")
        for item in selected or eligible:
            rows[item["row_index"]]["additional_toxicity"].update(status="unassessed_model_stage_held", error_type=type(error).__name__)
        receipt.update(status="held", reason="additional_model_stage_held", error_type=type(error).__name__)
    return rows, receipt


def main():
    # CPU inference is deliberately separate from the donor GPU wrapper/model
    # channel. Disable network APIs even if an imported dependency tries them.
    import socket
    def no_network(*args, **kwargs):
        raise RuntimeError("Additional local model runtime is offline")
    socket.create_connection = no_network
    socket.socket.connect = no_network
    request_path, result_path = map(Path, sys.argv[1:])
    result = _run(json.loads(request_path.read_bytes()))
    atomic_bytes(result_path, canonical_json(result))
    # This second publication proves the result file and its containing
    # directory have finished their first fsync/atomic-rename transaction.
    atomic_bytes(result_path.with_name("result.complete.json"), canonical_json({"result_sha256": sha256_file(result_path),
        "input_sha256": result["input_sha256"], "model_identity_sha256": digest(result["identity"])}))


if __name__ == "__main__":
    main()
