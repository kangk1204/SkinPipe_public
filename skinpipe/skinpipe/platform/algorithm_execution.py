"""Artifact-bound optional execution under the existing worker lease and budget.

Administrator profiles supply pinned runtimes. Research requests carry only
registered candidates and CAS identifiers, never executable or filesystem paths.
Model outputs retain their units and have zero automatic ranking weight.
"""
from __future__ import annotations

from copy import deepcopy
import hashlib
import json
from pathlib import Path
import time
import uuid

from .algorithm_profiles import PROVIDERS, provider_limit, validate_profile
from .artifacts import ArtifactIntegrityError, atomic_bytes, canonical_json, sha256_file
from .engines import ComputeBudgetExceeded, EngineCancelled, run_owned
from .identity import digest

CHEMISTRY = {"fpsim2", "cdpl", "dimorphite", "posebusters", "prolif"}
POCKETS = {"fpocket", "p2rank"}


class AlgorithmRunner:
    def __init__(self, root, profiles, *, gpu_uuid=None, cpu_threads=1):
        self.root = Path(root).resolve()
        self.root.mkdir(parents=True, exist_ok=True)
        self.profiles = deepcopy(profiles)
        self.gpu_uuid, self.cpu_threads = gpu_uuid, cpu_threads

    @staticmethod
    def _materialize(sha256, path, fetch_artifact):
        if not isinstance(sha256, str) or len(sha256) != 64 or any(c not in "0123456789abcdef" for c in sha256):
            raise ValueError("Exact artifact SHA256 required")
        raw = fetch_artifact(sha256)
        if hashlib.sha256(raw).hexdigest() != sha256:
            raise ArtifactIntegrityError("Algorithm input artifact changed")
        if len(raw) > 64 * 1024 * 1024:
            raise ValueError("Algorithm input exceeds64MiB")
        atomic_bytes(path, raw)
        return str(path)

    @staticmethod
    def _publish_files(value, directory, publish_artifact):
        """Replace verified owned output file receipts with immutable CAS receipts."""
        if isinstance(value, list):
            return [AlgorithmRunner._publish_files(v, directory, publish_artifact) for v in value]
        if not isinstance(value, dict):
            return value
        if "path" in value and "sha256" in value and "module_sha256" not in value:
            path = Path(value["path"])
            if path.is_symlink() or not path.is_file() or not path.resolve().is_relative_to(directory):
                raise ArtifactIntegrityError("Algorithm output is outside its owned directory")
            if path.stat().st_size > 64 * 1024 * 1024 or sha256_file(path) != value["sha256"]:
                raise ArtifactIntegrityError("Algorithm output size/hash mismatch")
            artifact = publish_artifact(path.read_bytes())
            if artifact.get("sha256") != value["sha256"]:
                raise ArtifactIntegrityError("Algorithm publication hash mismatch")
            return {**{k: v for k, v in value.items() if k != "path"}, **artifact}
        return {k: (deepcopy(v) if k == "software" else AlgorithmRunner._publish_files(v, directory, publish_artifact))
                for k, v in value.items() if k not in {"artifact_directory", "log_path"}}

    def execute(self, request, *, fetch_artifact, publish_artifact, cancel_event,
                remaining_seconds=3600, progress=None, reserve_launch=None, prior_result=None):
        started = time.monotonic()
        if remaining_seconds <= 0:
            raise ComputeBudgetExceeded('No optional compute budget remains')
        profile = validate_profile(request["compute_profile"])
        provider = request["provider"]
        if provider not in PROVIDERS:
            raise ValueError("Unsupported research algorithm")
        qualification = request.get("qualification_mode", False)
        if type(qualification) is not bool:
            raise ValueError("Qualification mode must be boolean")
        work = self.root / uuid.uuid4().hex
        work.mkdir()
        source = work / "registered-inputs.json"
        self._materialize(request["input_sha256"], source, fetch_artifact)
        payload = json.loads(source.read_text())
        if payload.get("campaign_id") != request["campaign_id"] or payload.get("provider") != provider:
            raise ArtifactIntegrityError("Algorithm request differs from frozen inputs")
        rows = payload["rows"]
        ids = [r["candidate_id"] for r in rows]
        cap = provider_limit(provider, qualification=qualification)
        if len(rows) > cap or len(set(ids)) != len(ids):
            raise ValueError("Algorithm request exceeds its unique candidate cap")
        if any(k.endswith("_path") or k in {"argv", "command", "shell", "executable"}
               for obj in [payload.get("context", {}), *rows] for k in obj):
            raise ValueError("Research inputs may not supply filesystem paths or commands")
        context = deepcopy(payload.get("context", {}))
        context.update(cpu_threads=self.cpu_threads, gpu_uuid=self.gpu_uuid,
                       qualification_mode=qualification, budget_seconds=min(remaining_seconds, 28800))
        context.update(conformers=profile["multiplicities"]["conformers"],
                       total_states=profile["multiplicities"]["prepared_states_total"],
                       reference_limit=profile["multiplicities"]["alignment_references"], candidate_limit=len(rows))
        result = {"provider": provider, "input_sha256": request["input_sha256"],
                  "compute_profile_sha256": profile["sha256"], "candidates": [],
                  "ranking_weight": 0, "population_preserved": True, "scientific_validation": "NOT_ESTABLISHED"}
        configured = self.profiles.get(provider)
        if not configured:
            result.update(status="held", reason="operator_runtime_profile_missing",
                          candidates=[{"candidate_id": cid, "status": "held", "reason": "operator_runtime_profile_missing"} for cid in ids], compute_seconds=time.monotonic()-started)
            self._bind(result, rows, context)
            result["result_artifact"] = publish_artifact(canonical_json(result))
            return result
        try:
            if provider in CHEMISTRY:
                from .algorithm_qualification import verify_chemistry_qualification
                verify_chemistry_qualification(provider, configured, fetch_artifact=fetch_artifact,
                    qualification_mode=qualification, input_sha256=request["input_sha256"],
                    source_receipt_sha256=request.get("source_receipt_sha256"))
            elif provider in POCKETS:
                from .algorithm_pocket_qualification import verify_pocket_qualification
                admitted=verify_pocket_qualification(provider,configured,fetch_artifact=fetch_artifact,
                    qualification_mode=qualification,input_sha256=request['input_sha256'],
                    source_receipt_sha256=request.get('source_receipt_sha256'))
                result['asset_identity_sha256']=admitted['sha256']
            elif provider not in POCKETS:
                from .algorithm_qualification import verify_model_qualification
                verify_model_qualification(provider, configured, fetch_artifact=fetch_artifact,
                    qualification_mode=qualification, input_sha256=request["input_sha256"],
                    source_receipt_sha256=request.get("source_receipt_sha256"))
        except ValueError as error:
            result.update(status="held", reason=str(error)[:500], compute_seconds=time.monotonic()-started,
                          candidates=[{"candidate_id":cid,"status":"held","reason":str(error)[:500]} for cid in ids])
            self._bind(result,rows,context)
            result["result_artifact"] = publish_artifact(canonical_json(result))
            return result
        # Materialize one receptor/sequence, then copy exact inputs into each model
        # directory. Models deliberately never share an unverified generated pose.
        for name, sha256 in payload.get("context_artifacts", {}).items():
            if name not in {"receptor", "receptor_identity", "sequence_reference", "msa", "reference"}:
                raise ValueError("Unsupported registered context artifact")
            suffix = ".pdb" if name == "receptor" else ".sdf" if name == "reference" else ".json" if name == "sequence_reference" else ".a3m"
            context[name + "_path"] = self._materialize(sha256, work / (name + suffix), fetch_artifact)
            context[name + "_sha256"] = sha256
        def owned(argv, *, cwd, timeout_seconds, env=None):
            if cancel_event.is_set():
                raise EngineCancelled("Algorithm lease lost")
            elapsed = time.monotonic() - started
            def account(current):
                if progress:
                    progress({"stage": "algorithm_" + provider, "compute_seconds": elapsed + current})
            execution = run_owned(argv, cwd=cwd, timeout=min(timeout_seconds, remaining_seconds-elapsed),
                                  cancel_event=cancel_event, gpu_uuid=self.gpu_uuid, env=env, on_compute=account)
            return {**execution, 'env':dict(env or {}), "status": "complete" if execution["returncode"] == 0 else "failed"}

        if provider in POCKETS:
            from .algorithm_pockets import run_pocket_search
            pocket_request = {"method": provider, "receptor_sha256": context.get("receptor_sha256"),
                              "structural_source": context.get("structural_source", "unknown"),
                              "timeout_seconds": min(remaining_seconds, 3600)}
            branch = run_pocket_search(pocket_request, configured, fetch_artifact=fetch_artifact,
                                      publish_artifact=publish_artifact, run_owned=owned)
            if branch.get('result_artifact'):
                branch['native_result_artifact'] = branch.pop('result_artifact')
            result.update(branch)
            result["candidates"] = [{"candidate_id": cid, "status": branch["status"], "reason": branch.get("reason")} for cid in ids]
        elif provider in CHEMISTRY:
            python = Path(configured.get("python", ""))
            lock = Path(configured.get("runtime_lock", ""))
            if (not python.is_absolute() or not python.is_file() or sha256_file(python) != configured.get("python_sha256")
                    or not lock.is_file() or sha256_file(lock) != configured.get("runtime_lock_sha256")):
                raise ArtifactIntegrityError("Chemistry runtime is missing or changed")
            for row in rows:
                if context.get("reference_path"):
                    row["reference_path"] = context["reference_path"]
                    row["reference_sha256"] = context["reference_sha256"]
                artifact = row.pop("pose_sdf_sha256", None)
                if artifact:
                    row["pose_path"] = self._materialize(artifact, work / (digest(row["candidate_id"]) + ".sdf"), fetch_artifact)
                    row["pose_sha256"] = artifact
            command_input = work / "chemistry-input.json"
            output = work / "chemistry-output.json"
            atomic_bytes(command_input, canonical_json({"provider": provider, "rows": rows, "context": context}))
            env = {"PYTHONPATH": str(Path(__file__).resolve().parents[2]),
                   "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"}
            execution = owned([str(python), "-m", "skinpipe.platform.algorithm_entry", "--input", str(command_input), "--output", str(output)], cwd=work, timeout_seconds=remaining_seconds, env=env)
            if execution["returncode"] != 0 or not output.is_file() or output.stat().st_size > 64*1024*1024:
                result.update(status="held", reason="isolated_algorithm_process_failed", candidates=[{"candidate_id": cid, "status": "held"} for cid in ids])
            else:
                branch = json.loads(output.read_text())
                if [r.get("candidate_id") for r in branch.get("candidates", [])] != ids:
                    raise ArtifactIntegrityError("Chemistry output candidate identity/order mismatch")
                result.update(self._publish_files(branch, work, publish_artifact))
                if provider=='dimorphite':
                    for receipt,row in zip(result['candidates'],rows):
                        metrics=receipt.get('metrics',{})
                        states=metrics.get('states',[])
                        original=row.get('prepared_exact_structure_id')
                        if original and row.get('original_prepared_smiles'):
                            metrics['workload_states']=[{'exact_structure_id':original,'smiles':row['original_prepared_smiles'],'role':'frozen_broad_prepared_state'}]
                            metrics['workload_states']+= [{**s,'role':'additional_sensitivity_state'} for s in states if s['exact_structure_id']!=original][:1]
                        else:
                            metrics['workload_states']=states[:2]
                            metrics['original_prepared_state_status']='unavailable; freeze before physical state docking'
                        metrics['automatic_state_docking_count']=0
        else:
            from .algorithm_models import prepare_model_command, parse_model_outputs
            result['runtime_device'] = 'gpu' if self.gpu_uuid else 'cpu'
            previous = {}
            if prior_result:
                if (prior_result.get("input_sha256") != request["input_sha256"]
                        or prior_result.get("runtime_profile_sha256") != digest(configured)):
                    raise ArtifactIntegrityError("Algorithm checkpoint has a different frozen identity")
                for saved in prior_result.get("candidates", []):
                    if saved.get("status") == "complete" and saved.get("candidate_id") in ids:
                        for artifact in saved.get("artifacts", []):
                            raw = fetch_artifact(artifact["sha256"])
                            if hashlib.sha256(raw).hexdigest() != artifact["sha256"]:
                                raise ArtifactIntegrityError("Algorithm checkpoint artifact changed")
                        previous[saved["candidate_id"]] = saved
            def checkpoint():
                result["compute_seconds"] = time.monotonic() - started
                result["runtime_profile_sha256"] = digest(configured)
                self._bind(result, rows, context)
                if progress:
                    progress({"stage": "algorithm_" + provider, "compute_seconds": result["compute_seconds"],
                              "partial_result": deepcopy(result)})
            stopped = None
            for row in rows:
                if cancel_event.is_set():
                    raise EngineCancelled("Algorithm lease lost")
                if row["candidate_id"] in previous:
                    result["candidates"].append(deepcopy(previous[row["candidate_id"]]))
                    continue
                if stopped:
                    result["candidates"].append({"candidate_id": row["candidate_id"], "status": stopped, "reason": "actual_compute_budget_exhausted"})
                    continue
                directory = work / digest(row["candidate_id"])
                directory.mkdir()
                for name in ("receptor", "sequence_reference", "msa"):
                    if context.get(name + "_sha256"):
                        row[name + "_path"] = self._materialize(context[name + "_sha256"], directory / Path(context[name + "_path"]).name, fetch_artifact)
                        row[name + "_sha256"] = context[name + "_sha256"]
                if row.get("pose_sdf_sha256"):
                    row["pose_sdf_path"] = self._materialize(row["pose_sdf_sha256"], directory / "pose.sdf", fetch_artifact)
                try:
                    prepared = prepare_model_command(provider, row, context, configured, directory)
                except (ValueError, OSError, KeyError) as error:
                    result["candidates"].append({"candidate_id": row["candidate_id"], "status": "held", "reason": type(error).__name__ + ":" + str(error)[:500]})
                    checkpoint()
                    continue
                if prepared["status"] != "ready":
                    result["candidates"].append({k: v for k, v in prepared.items() if k != "argv"})
                    continue
                if provider in {"boltz2", "chai1"} and (reserve_launch is None or not reserve_launch(digest(row["candidate_id"]))):
                    result["candidates"].append({"candidate_id": row["candidate_id"], "status": "budget_deferred", "reason": "combined_physical_launch_limit"})
                    continue
                try:
                    execution = owned(prepared["argv"], cwd=directory, timeout_seconds=min(3600, remaining_seconds), env=prepared.get("env"))
                    if execution["returncode"] == 0:
                        receipt = parse_model_outputs(provider, directory, row)
                        artifacts = []
                        for artifact in receipt.get("artifacts", []):
                            path = directory / artifact["name"]
                            artifacts.append(self._publish_files({"path": str(path), **artifact}, work, publish_artifact))
                        receipt["artifacts"] = artifacts
                        result["candidates"].append(receipt)
                    else:
                        result["candidates"].append({"candidate_id": row["candidate_id"], "status": "failed", "reason": "owned_model_process_failed"})
                    result.setdefault("executions", []).append({k: execution.get(k) for k in ("returncode", "elapsed_seconds", "peak_mib", "process_peak_rss_bytes")})
                    peaks = [e.get('process_peak_rss_bytes') or 0 for e in result['executions']]
                    result['measurements'] = {'process_peak_rss_bytes': max(peaks), 'memory_scope':'sampled_process_group_RSS'}
                except ComputeBudgetExceeded:
                    stopped = "budget_deferred"
                    result["candidates"].append({"candidate_id": row["candidate_id"], "status": stopped, "reason": "actual_compute_budget_exhausted"})
                except (ValueError, OSError, KeyError) as error:
                    result["candidates"].append({"candidate_id": row["candidate_id"], "status": "held", "reason": type(error).__name__ + ":" + str(error)[:500]})
                result["asset_identity_sha256"] = prepared["asset_identity"]["sha256"]
                checkpoint()
            statuses = [r["status"] for r in result["candidates"]]
            result["status"] = "complete" if statuses and all(s == "complete" for s in statuses) else "partial" if "complete" in statuses else "held"
        result["compute_seconds"] = time.monotonic() - started
        result["runtime_profile_sha256"] = digest(configured)
        self._bind(result, rows, context)
        raw = canonical_json(result)
        result["result_artifact"] = publish_artifact(raw)
        if result["result_artifact"].get("sha256") != hashlib.sha256(raw).hexdigest():
            raise ArtifactIntegrityError("Algorithm result publication mismatch")
        return result

    @staticmethod
    def _bind(result, rows, context):
        by_id = {r["candidate_id"]: r for r in rows}
        for receipt in result["candidates"]:
            row = by_id[receipt["candidate_id"]]
            receipt["target_id"] = context.get("target_id")
            receipt["input_exact_structure_id"] = row.get("exact_structure_id")
            receipt["prepared_exact_structure_id"] = row.get("prepared_exact_structure_id")
        return result
