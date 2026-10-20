"""Lease-owned, offline evaluation of an exact, bounded CAS input bundle.

Only operator-pinned Python/runtime assets may execute. The child receives local
hash-addressed bytes and a source snapshot; no API client or credentials cross
the process boundary. Evaluation never activates or promotes a model.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import hashlib
import json
import math
import os
from pathlib import Path
import re
import time
import uuid

from .artifacts import ArtifactIntegrityError, atomic_bytes, canonical_json, sha256_file

MAX_TRANSFER_BYTES = 256 * 1024 * 1024
MAX_JSON_BYTES = 32 * 1024 * 1024
MAX_OUTPUT_BYTES = 32 * 1024 * 1024
SHA = re.compile(r"[0-9a-f]{64}\Z")
SOURCE_FILES = (
    "skinpipe/__init__.py", "skinpipe/platform/__init__.py",
    "skinpipe/platform/evaluation_execution.py", "skinpipe/platform/artifacts.py",
    "skinpipe/platform/algorithm_evaluation.py", "skinpipe/platform/algorithm_pockets.py",
    "skinpipe/platform/evaluation.py", "skinpipe/platform/identity.py",
    "skinpipe/platform/ranking.py",
)
# -I blocks PYTHONPATH, user site and the script directory before this trusted
# bootstrap runs. Clear inherited secrets before importing the owned snapshot.
BOOTSTRAP = (
    "import os,runpy,sys; "
    "keys=('PATH','HOME','TMPDIR','TEMP','TMP','LANG','LC_ALL','SYSTEMROOT',"
    "'OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','CUDA_VISIBLE_DEVICES',"
    "'HF_HUB_OFFLINE','TRANSFORMERS_OFFLINE','HF_DATASETS_OFFLINE'); "
    "clean={k:os.environ[k] for k in keys if k in os.environ}; "
    "os.environ.clear(); os.environ.update(clean); "
    "sys.path.insert(0,sys.argv.pop(1)); "
    "runpy.run_module('skinpipe.platform.evaluation_execution',run_name='__main__')"
)


class EvaluationHeld(ValueError):
    """An optional evaluation prerequisite or bounded input is unavailable."""


def _sha(value):
    if not isinstance(value, str) or not SHA.fullmatch(value):
        raise ArtifactIntegrityError("Exact lowercase SHA256 required")
    return value


def _pairs(pairs):
    value = {}
    for key, item in pairs:
        if key in value:
            raise ArtifactIntegrityError("Duplicate evaluation JSON key")
        value[key] = item
    return value


def _json(raw):
    if not isinstance(raw, bytes) or len(raw) > MAX_JSON_BYTES:
        raise EvaluationHeld("evaluation_JSON_exceeds32MiB")
    return json.loads(raw, object_pairs_hook=_pairs,
                      parse_constant=lambda _: (_ for _ in ()).throw(ArtifactIntegrityError("Nonfinite JSON")))


def _file(path, expected, *, executable=False):
    if not isinstance(path, str) or not Path(path).is_absolute():
        raise EvaluationHeld("operator_runtime_profile_required")
    supplied = Path(path)
    resolved = supplied.resolve()
    if not resolved.is_file() or supplied.is_symlink() and not executable:
        raise EvaluationHeld("operator_runtime_asset_missing_or_unsafe")
    if executable and not os.access(resolved, os.X_OK):
        raise EvaluationHeld("operator_python_not_executable")
    if sha256_file(resolved) != _sha(expected):
        raise ArtifactIntegrityError("Pinned evaluation runtime changed")
    # Launch the configured venv path: resolving its Python symlink would lose
    # pyvenv.cfg discovery and silently execute the base interpreter instead.
    return {"path": str(supplied) if executable else str(resolved),
            "sha256": expected, "size_bytes": resolved.stat().st_size}


class EvaluationRunner:
    def __init__(self, root, profile=None):
        self.root = Path(root).resolve()
        self.root.mkdir(parents=True, exist_ok=True)
        self.profile = deepcopy(profile or {})

    def _runtime(self, work):
        if not isinstance(self.profile, dict):
            raise EvaluationHeld("operator_evaluation_profile_missing")
        python = _file(self.profile.get("python"), self.profile.get("python_sha256"), executable=True)
        lock = _file(self.profile.get("runtime_lock"), self.profile.get("runtime_lock_sha256"))
        extra = self.profile.get("runtime_files", [])
        if (not isinstance(extra, list) or len(extra) > 10000
                or any(not isinstance(item, dict) for item in extra)):
            raise ValueError("Bounded pinned runtime file manifest required")
        runtime_files = [_file(item.get("path"), item.get("sha256")) for item in extra]
        package_root = Path(__file__).resolve().parents[2]
        sources = []
        for name in SOURCE_FILES:
            source = package_root / name
            raw = source.read_bytes()
            if len(raw) > 4 * 1024 * 1024:
                raise EvaluationHeld("evaluation_source_exceeds4MiB")
            atomic_bytes(work / "runtime" / name, raw)
            sources.append({"name": name, "sha256": hashlib.sha256(raw).hexdigest()})
        identity = {"schema": "skinpipe-evaluation-runtime-v1", "python": python,
                    "runtime_lock": lock, "runtime_files": runtime_files, "sources": sources}
        identity["sha256"] = hashlib.sha256(canonical_json(identity)).hexdigest()
        return python, identity

    def execute(self, request, *, fetch_artifact, publish_artifact, cancel_event,
                remaining_seconds=3600, progress=None):
        from .engines import ComputeBudgetExceeded, EngineCancelled, run_owned

        started = time.monotonic()
        result = {"provider": "evaluation", "status": "held", "promoted": False,
                  "ranking_weight": 0, "scientific_validation": "NOT_ESTABLISHED",
                  "population_preserved": True, "compute_seconds": 0}
        work = self.root / uuid.uuid4().hex
        work.mkdir()

        def live():
            if cancel_event.is_set():
                raise EngineCancelled("Evaluation lease lost or cancelled")
            if time.monotonic() - started >= remaining_seconds:
                raise ComputeBudgetExceeded("Evaluation compute budget exhausted")

        def update(stage, seconds=None):
            if progress:
                progress({"stage": stage, "compute_seconds": time.monotonic()-started if seconds is None else seconds})

        try:
            if (not isinstance(request, dict) or set(request) != {"evaluation_manifest_sha256"}
                    or isinstance(remaining_seconds, bool) or not isinstance(remaining_seconds, (int, float))
                    or not math.isfinite(remaining_seconds) or remaining_seconds <= 0):
                raise ValueError("One immutable evaluation manifest and positive finite budget required")
            manifest_sha = _sha(request["evaluation_manifest_sha256"])
            result["evaluation_manifest_sha256"] = manifest_sha
            live()
            python, runtime = self._runtime(work)
            result["runtime_profile_sha256"] = runtime["sha256"]
            result["runtime_identity"] = runtime
            collected, total = {}, 0

            def get(sha256, *, as_json=False):
                nonlocal total
                live()
                sha256 = _sha(sha256)
                if sha256 not in collected:
                    raw = fetch_artifact(sha256)
                    live()
                    if not isinstance(raw, bytes) or hashlib.sha256(raw).hexdigest() != sha256:
                        raise ArtifactIntegrityError("Evaluation CAS bytes differ from SHA256")
                    total += len(raw)
                    if total > MAX_TRANSFER_BYTES:
                        raise EvaluationHeld("evaluation_input_transfer_exceeds256MiB")
                    name = "inputs/" + sha256 + ".blob"
                    atomic_bytes(work / name, raw)
                    collected[sha256] = {"name": name, "size_bytes": len(raw)}
                raw = (work / collected[sha256]["name"]).read_bytes()
                return _json(raw) if as_json else raw

            update("evaluation_materializing")
            manifest = get(manifest_sha, as_json=True)
            if not isinstance(manifest, dict) or manifest.get("schema") != "skinpipe-evaluation-v1":
                raise ValueError("Unsupported evaluation manifest")
            for key in ("snapshot_sha256", "policy_sha256"):
                get(manifest[key], as_json=True)
            bindings = {}
            for key in ("model_binding_sha256", "baseline_model_binding_sha256"):
                binding = get(manifest[key], as_json=True)
                if not isinstance(binding, dict):
                    raise ValueError("Model binding must be a record")
                get(binding["training_rows_sha256"], as_json=True)
                get(binding["model_sha256"])
                bindings[key] = binding
            for key, binding_key in (("predictions_sha256", "model_binding_sha256"),
                                     ("baseline_predictions_sha256", "baseline_model_binding_sha256")):
                prediction = get(manifest[key], as_json=True)
                if not isinstance(prediction, dict):
                    raise ValueError("Prediction bundle must be a record")
                if (prediction.get("model_sha256") != bindings[binding_key]["model_sha256"]
                        or prediction.get("snapshot_sha256") != manifest["snapshot_sha256"]):
                    raise ArtifactIntegrityError("Prediction bundle refers to an unbound model or snapshot")
                get(prediction["model_sha256"])
            bundle = {"schema": "skinpipe-evaluation-local-input-v1", "request": request,
                      "artifacts": collected, "runtime_identity": runtime}
            bundle_path, output_path = work / "bundle.json", work / "output.json"
            atomic_bytes(bundle_path, canonical_json(bundle))
            home = work / "home"
            home.mkdir()
            env = {"PATH": str(Path(python["path"]).parent), "HOME": str(home), "TMPDIR": str(work),
                   "TEMP": str(work), "TMP": str(work), "LANG": "C.UTF-8", "LC_ALL": "C.UTF-8",
                   "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
                   "CUDA_VISIBLE_DEVICES": "", "HF_HUB_OFFLINE": "1", "TRANSFORMERS_OFFLINE": "1",
                   "HF_DATASETS_OFFLINE": "1"}
            live()
            process_started = time.monotonic()
            update("evaluation_compute")
            execution = run_owned([python["path"], "-I", "-B", "-c", BOOTSTRAP, str(work / "runtime"),
                "--bundle", str(bundle_path), "--output", str(output_path)], cwd=work,
                timeout=remaining_seconds-(time.monotonic()-started), cancel_event=cancel_event, env=env,
                on_compute=lambda elapsed: update("evaluation_compute", process_started-started+elapsed))
            live()
            # Recheck mutable administrator files after the isolated execution.
            for asset in [runtime["python"], runtime["runtime_lock"], *runtime["runtime_files"]]:
                _file(asset["path"], asset["sha256"], executable=asset is runtime["python"])
            if execution["returncode"] != 0:
                raise EvaluationHeld("isolated_evaluation_process_failed")
            envelope = _json(_owned_bytes(output_path, work, MAX_OUTPUT_BYTES))
            if (not isinstance(envelope, dict) or envelope.get("schema") != "skinpipe-evaluation-local-result-v1"
                    or envelope.get("evaluation_manifest_sha256") != manifest_sha
                    or envelope.get("runtime_profile_sha256") != runtime["sha256"]):
                raise ArtifactIntegrityError("Evaluation output differs from its frozen bundle")
            evaluation = envelope.get("evaluation")
            if (not isinstance(evaluation, dict) or evaluation.get("status") not in {"held", "evaluated"}
                    or evaluation.get("promoted") is not False):
                raise ArtifactIntegrityError("Evaluation output attempted unsupported promotion")
            published, output_total = {}, 0
            outputs = envelope.get("artifacts")
            if not isinstance(outputs, list) or len(outputs) > 8:
                raise ValueError("Bounded local evaluation output manifest required")
            for item in outputs:
                if not isinstance(item, dict):
                    raise ArtifactIntegrityError("Evaluation output receipt must be a record")
                sha256 = _sha(item["sha256"])
                if item.get("name") != "outputs/" + sha256 + ".blob" or sha256 in published:
                    raise ArtifactIntegrityError("Invalid local evaluation output mapping")
                raw = _owned_bytes(work / item["name"], work, MAX_OUTPUT_BYTES)
                output_total += len(raw)
                if (output_total > MAX_OUTPUT_BYTES or hashlib.sha256(raw).hexdigest() != sha256
                        or type(item.get("size_bytes")) is not int or item["size_bytes"] != len(raw)):
                    raise ArtifactIntegrityError("Evaluation output bytes/size differ from receipt")
                live()
                receipt = publish_artifact(raw)
                if receipt.get("sha256") != sha256 or receipt.get("size_bytes") != len(raw):
                    raise ArtifactIntegrityError("Evaluation CAS publication differs from bytes")
                published[sha256] = receipt
            if evaluation["status"] == "evaluated":
                artifact = evaluation.get("artifact", {})
                expected = hashlib.sha256(canonical_json({k: v for k, v in evaluation.items() if k != "artifact"})).hexdigest()
                if artifact.get("sha256") != expected or expected not in published:
                    raise ArtifactIntegrityError("Evaluation receipt is not its actual CAS output")
                evaluation["artifact"] = published[expected]
            result.update(status="complete" if evaluation["status"] == "evaluated" else "held",
                          evaluation=evaluation, input_transfer_bytes=total,
                          input_artifacts=[{"sha256": sha, "size_bytes": item["size_bytes"]} for sha,item in sorted(collected.items())],
                          execution={key: execution[key] for key in ("returncode", "elapsed_seconds", "peak_mib")})
            if result["status"] == "held":
                result["reason"] = evaluation.get("reason", "evaluation_prerequisite_held")
        except (ValueError, OSError, KeyError, TypeError, EngineCancelled, ComputeBudgetExceeded) as error:
            result.update(status="held", reason=type(error).__name__ + ":" + str(error)[:1000])
        result["compute_seconds"] = time.monotonic() - started
        if not cancel_event.is_set():
            update("evaluation_" + result["status"])
        if not cancel_event.is_set():
            raw = canonical_json(result)
            artifact = publish_artifact(raw)
            if artifact.get("sha256") != hashlib.sha256(raw).hexdigest() or artifact.get("size_bytes") != len(raw):
                raise ArtifactIntegrityError("Evaluation result publication differs from bytes")
            result["result_artifact"] = artifact
        return result


def _owned_bytes(path, root, maximum):
    path, root = Path(path), Path(root).resolve()
    if (path.is_symlink() or not path.is_file() or not path.resolve().is_relative_to(root)
            or path.stat().st_size > maximum):
        raise ArtifactIntegrityError("Evaluation local file is missing, unsafe or excessive")
    return path.read_bytes()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bundle", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    root = Path(args.bundle).resolve().parent
    bundle = _json(_owned_bytes(Path(args.bundle), root, MAX_JSON_BYTES))
    request, references = bundle["request"], bundle["artifacts"]
    _sha(request["evaluation_manifest_sha256"])
    if not isinstance(references, dict) or len(references) > 16:
        raise ValueError("Invalid bounded local CAS mapping")
    artifacts = {}

    def fetch(sha256):
        _sha(sha256)
        item = references[sha256]
        if item.get("name") != "inputs/" + sha256 + ".blob":
            raise ArtifactIntegrityError("Unknown local input mapping")
        raw = _owned_bytes(root / item["name"], root, MAX_TRANSFER_BYTES)
        if hashlib.sha256(raw).hexdigest() != sha256 or item.get("size_bytes") != len(raw):
            raise ArtifactIntegrityError("Local input changed after collection")
        return raw

    def publish(raw):
        if not isinstance(raw, bytes) or len(raw) > MAX_OUTPUT_BYTES:
            raise ValueError("Bounded evaluation output required")
        sha256 = hashlib.sha256(raw).hexdigest()
        name = "outputs/" + sha256 + ".blob"
        atomic_bytes(root / name, raw)
        receipt = {"sha256": sha256, "size_bytes": len(raw)}
        artifacts[sha256] = {**receipt, "name": name}
        return receipt

    from .algorithm_evaluation import evaluate_comparison
    evaluation = evaluate_comparison(request, fetch_artifact=fetch, publish_artifact=publish)
    envelope = {"schema": "skinpipe-evaluation-local-result-v1",
                "evaluation_manifest_sha256": request["evaluation_manifest_sha256"],
                "runtime_profile_sha256": bundle["runtime_identity"]["sha256"],
                "evaluation": evaluation, "artifacts": list(artifacts.values())}
    output = Path(args.output)
    if output.resolve().parent != root or output.is_symlink():
        raise ValueError("Evaluation output must belong to its owned directory")
    atomic_bytes(output, canonical_json(envelope))


if __name__ == "__main__":
    main()
