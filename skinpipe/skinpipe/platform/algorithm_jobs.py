'Validate the selected structure and retain explicit calculation provenance.'

from __future__ import annotations

from copy import deepcopy
import json

from .algorithm_profiles import PROVIDERS, compute_profile, provider_limit, select_bands
from .artifacts import ArtifactIntegrityError, canonical_json
from .campaign import _sealed_read, _sealed_write
from .identity import digest
from .algorithm_manifest import freeze_manifest

POSE_PROVIDERS = {"posebusters", "prolif", "gnina", "boltzina", "unidock2"}


class AlgorithmJobs:
    def __init__(self, campaigns, artifacts, queue):
        self.campaigns, self.artifacts, self.queue = campaigns, artifacts, queue

    def submit(self, campaign_id, body, *, owner="researcher", automatic=False):
        if not isinstance(body, dict) or set(body) - {"provider", "qualification_mode"}:
            raise ValueError(
                "Algorithm requests accept provider and qualification_mode only"
            )
        provider, qualification = (
            body.get("provider"),
            body.get("qualification_mode", False),
        )
        if (
            not isinstance(provider, str)
            or provider not in PROVIDERS
            or type(qualification) is not bool
        ):
            raise ValueError("Invalid research algorithm or qualification mode")
        if provider == "openfe":
            return {
                "status": "held",
                "reason": "OpenFE requires a separately frozen post-hit transformation",
                "jobs": [],
                "population_preserved": True,
            }
        coordinator = self.campaigns._coordinator(owner)
        manifest = coordinator._manifest(campaign_id)
        profile = manifest["plan"].get("compute_profile") or compute_profile()
        if automatic and (
            provider != "posebusters" or profile["version"] != "generalized-research-v1"
        ):
            raise ValueError(
                "Only representative pose quality is automatic in generalized runs"
            )
        available = [
            node
            for node in self.queue.resources()["nodes"]
            if node.get("fresh")
            and "algorithm_" + provider
            in node.get("telemetry", {}).get("capabilities", [])
        ]
        if not available:
            return {
                "status": "held",
                "reason": "registered_algorithm_worker_unavailable",
                "jobs": [],
                "provider": provider,
                "population_preserved": True,
                "ranking_weight": 0,
            }
        records = []
        with coordinator._lock(campaign_id) as directory:
            state = _sealed_read(directory / "state.json")
            if "broad_docking" in state["stages"]:
                records, missing = coordinator._stage_rows(directory, state, "broad_docking")
                if provider in POSE_PROVIDERS:
                    pending = [cid for cid in missing
                               if manifest["prepared"].get(cid, {}).get("pdbqt_sha256")]
                    failed_gpu = {r["candidate_id"] for r in records
                                  if r.get("engine") == "unidock"
                                  and r.get("status") not in {"complete", "review_required"}}
                    if pending:
                        return {"status": "held", "reason": "broad_pose_receipts_pending",
                                "missing_count": len(pending), "jobs": [],
                                "population_preserved": True}
                    if failed_gpu:
                        if "broad_cpu_fallback" not in state["stages"]:
                            return {"status": "held", "reason": "broad_cpu_fallback_receipts_pending",
                                    "missing_count": len(failed_gpu), "jobs": [],
                                    "population_preserved": True}
                        fallback, missing = coordinator._stage_rows(directory, state, "broad_cpu_fallback")
                        if missing:
                            return {"status": "held", "reason": "broad_cpu_fallback_receipts_pending",
                                    "missing_count": len(missing), "jobs": [],
                                    "population_preserved": True}



                        replacement = {r["candidate_id"]: r for r in fallback
                                       if r["candidate_id"] in failed_gpu}
                        records = [replacement.get(r["candidate_id"], r)
                                   if r["candidate_id"] in failed_gpu else r for r in records]
        ranked = coordinator._rank_results(manifest, records, unique_parents=False) if provider in POSE_PROVIDERS else coordinator._rank_results(manifest, records)
        target = manifest["plan"]["target_id"]
        # Include controls while selecting so their parent chemistry cannot
        # consume experimental candidate quotas even under a different supplier ID.
        controls = {r['candidate_id']: deepcopy(r)
            for r in manifest["candidates"] + ranked
            if r.get("role")
            in {"control", "positive_control", "negative_control", "counter_control"}
            or r.get("opposing_direction")
        }
        controls = list(controls.values())
        selector_rows = {r["candidate_id"]: r for r in ranked}
        for row in controls:
            selector_rows.setdefault(row["candidate_id"], row)
        by_id = {
            r["candidate_id"]: r
            for r in records
            if r.get("status") in {"complete", "review_required"}
            and r.get("sdf_artifact")
        }
        unavailable = []
        if provider in POSE_PROVIDERS:
            # Preserve full-population scores and reserve all control parents,
            # but only a verified pose can represent a parent for this branch.
            for cid in list(by_id):
                try:
                    self.artifacts.verified_path(by_id[cid]["sdf_artifact"]["sha256"])
                except (ArtifactIntegrityError, ValueError, KeyError, TypeError, OSError):
                    del by_id[cid]
            control_ids = {r["candidate_id"] for r in controls}
            unavailable = [{"candidate_id": cid, "reason": "verified_representative_sdf_pose_unavailable"}
                           for cid in sorted(set(selector_rows) - set(by_id) - control_ids)]
            selector_rows = {cid: row for cid, row in selector_rows.items() if cid in by_id or cid in control_ids}
        selection = select_bands(
            list(selector_rows.values()), target, provider, trust_stored_identity=True
        )
        rows, held = [], deepcopy(selection["omitted"]) + unavailable
        for candidate in selection["candidates"]:
            cid = candidate["candidate_id"]
            if provider in POSE_PROVIDERS and cid not in by_id:
                held.append(
                    {
                        "candidate_id": cid,
                        "reason": "verified_representative_sdf_pose_unavailable",
                    }
                )
                continue
            identity = candidate["identity"]
            row = {
                "candidate_id": cid,
                "smiles": candidate["smiles"],
                "target_id": target,
                "exact_structure_id": identity["exact_structure_id"],
                "parent_id": identity["parent_id"],
            }
            prepared = manifest["prepared"].get(cid, {})
            receipt_sha = prepared.get("preparation_receipt_sha256")
            if receipt_sha:
                from .identity import chemical_identity

                preparation = json.loads(self.artifacts.get_bytes(receipt_sha))
                row["original_prepared_smiles"] = preparation.get("effective_smiles")
                row["prepared_exact_structure_id"] = chemical_identity(
                    row["original_prepared_smiles"]
                )["exact_structure_id"]
            if cid in by_id:
                row["pose_sdf_sha256"] = by_id[cid]["sdf_artifact"]["sha256"]
                self.artifacts.verified_path(row["pose_sdf_sha256"])
            rows.append(row)
        rows = rows[: provider_limit(provider, qualification=qualification)]
        receptor = manifest.get("receptor") or {}
        context = {
            "target_id": target,
            "actual_protein_id": receptor.get("actual_protein_id"),
            "seed": 42,
            "center": receptor.get("center"),
            "size": receptor.get("size"),
            "reference_smiles": [
                r["smiles"]
                for r in manifest["plan"]["reference_snapshot"]["records"]
                if r.get("structure_status") == "valid"
            ][:5],
            "structural_source": "unknown",
        }
        context_artifacts = {}
        if receptor.get("uniprot_record_sha256"):
            context_artifacts["sequence_reference"] = receptor["uniprot_record_sha256"]
        if receptor.get("identity_card_sha256"):
            context_artifacts["receptor_identity"] = receptor["identity_card_sha256"]
            card = json.loads(
                self.artifacts.get_bytes(receptor["identity_card_sha256"])
            )
            for key in ("clean_pdb", "clean", "pdb"):
                entry = card.get("artifacts", {}).get(key, {})
                if entry.get("present") is True:
                    context_artifacts["receptor"] = entry["sha256"]
                    break
        for sha256 in context_artifacts.values():
            self.artifacts.verified_path(sha256)
        if provider in {"fpocket", "p2rank"}:
            rows = []
        if provider in {"boltz2", "unidock2"}:
            return {
                "status": "held",
                "reason": "registered_offline_MSA_required"
                if provider == "boltz2"
                else "registered_DMS_receptor_preparation_required",
                "jobs": [],
                "held": held,
                "population_preserved": True,
                "ranking_weight": 0,
            }
        missing = (
            provider in POSE_PROVIDERS | {"fpocket", "p2rank"}
            and "receptor" not in context_artifacts
            or provider in {"boltz2", "chai1"}
            and "sequence_reference" not in context_artifacts
            or provider not in {"fpocket", "p2rank"}
            and not rows
        )
        if missing:
            return {
                "status": "held",
                "reason": "registered_artifact_prerequisite_unavailable",
                "jobs": [],
                "selection": {k: v for k, v in selection.items() if k != "candidates"},
                "held": held,
                "population_preserved": True,
            }



        cohort_key = digest({"provider": provider, "profile": profile["sha256"]})
        with coordinator._lock(campaign_id) as directory:
            state = _sealed_read(directory / "state.json")
            cohorts = state.setdefault("algorithm_execution_cohorts", {})
            cohort = cohorts.setdefault(cohort_key, {"provider": provider,
                "compute_profile_sha256": profile["sha256"], "parents": {},
                "cap": provider_limit(provider)})
            admitted = []
            for row in rows:
                parent = row["parent_id"]
                if parent not in cohort["parents"] and len(cohort["parents"]) >= cohort["cap"]:
                    held.append({"candidate_id": row["candidate_id"], "reason": "campaign_provider_cohort_cap"})
                    continue
                cohort["parents"].setdefault(parent, row["candidate_id"])
                admitted.append(row)
            rows = admitted
            _sealed_write(directory / "state.json", state)
        if not rows and provider not in {"fpocket", "p2rank"}:
            return {"status": "held", "reason": "campaign_provider_cohort_cap", "jobs": [],
                    "held": held, "population_preserved": True, "ranking_weight": 0}
        jobs = []
        manifest_artifact, manifest_transport = freeze_manifest(manifest, self.artifacts)
        # Queue chunks are bounded independently of whole-target cohort caps.
        chunks = [rows[start : start + 16] for start in range(0, len(rows), 16)] or [[]]
        scope = (
            manifest["plan"]["budget_scope"]
            if automatic
            else "optional:"
            + digest(
                {"owner": owner, "campaign": campaign_id, "profile": profile["sha256"]}
            )
        )
        budget = (
            manifest["plan"]["budget_seconds"]
            if automatic
            else profile["budget"]["optional_inference_seconds"]
        )
        for chunk in chunks:
            frozen = {
                "campaign_id": campaign_id,
                "provider": provider,
                "rows": chunk,
                "context": context,
                "context_artifacts": context_artifacts,
                "campaign_manifest_sha256": manifest_artifact["sha256"],
            }
            artifact = self.artifacts.put_bytes(canonical_json(frozen))
            source = self.artifacts.put_bytes(
                canonical_json(
                    {
                        "schema": "skinpipe-registered-algorithm-source-v1",
                        "provider": provider,
                        "input_sha256": artifact["sha256"],
                        "source_type": "frozen_campaign",
                        "synthetic": False,
                        "campaign_id": campaign_id,
                        "manifest_sha256": manifest_artifact["sha256"],
                        **({"manifest_transport": manifest_transport} if manifest_transport else {}),
                        "source_artifacts": [
                            {
                                "sha256": (manifest_transport or manifest_artifact)["sha256"],
                                "role": "frozen_manifest",
                            }
                        ]
                        + [
                            {"sha256": sha, "role": name}
                            for name, sha in context_artifacts.items()
                        ]
                        + [
                            {"sha256": r["pose_sdf_sha256"], "role": "pose"}
                            for r in chunk
                            if r.get("pose_sdf_sha256")
                        ],
                    }
                )
            )
            request = {
                "campaign_id": campaign_id,
                "provider": provider,
                "input_sha256": artifact["sha256"],
                "source_receipt_sha256": source["sha256"],
                "compute_profile": profile,
                "qualification_mode": qualification,
                "budget_scope": scope,
            }
            job = self.queue.enqueue(
                "algorithm_" + provider,
                request,
                owner=owner,
                project="campaign-" + digest(owner)[:32],
                idempotency_key=digest(request),
                budget_seconds=budget,
            )
            jobs.append(
                {
                    "id": job["id"],
                    "status": job["status"],
                    "input_sha256": artifact["sha256"],
                }
            )
        executed_ids = [r["candidate_id"] for r in rows]
        executed_id_set = set(executed_ids)
        return {
            "status": "submitted",
            "provider": provider,
            "jobs": jobs,
            "held": held,
            "selection": {k: v for k, v in selection.items() if k != "candidates"},
            "execution_cohort": {
                "candidate_ids": executed_ids,
                "candidate_count": len(rows),
                "cap": provider_limit(provider, qualification=qualification),
                "campaign_unique_parent_count": len(cohort["parents"]),
                "qualification_mode": qualification,
                "selected_but_not_executed_ids": [
                    cid
                    for cid in selection["candidate_ids"]
                    if cid not in executed_id_set
                ],
            },
            "compute_profile_sha256": profile["sha256"],
            "population_preserved": True,
            "ranking_weight": 0,
        }

    def validate_result(self, request, result):
        frozen = json.loads(self.artifacts.get_bytes(request["input_sha256"]))
        expected = {r["candidate_id"]: r for r in frozen["rows"]}
        if not isinstance(result, dict):
            raise ArtifactIntegrityError("Algorithm result must be an object")
        rows = result.get("candidates", [])
        if (
            not isinstance(rows, list)
            or result.get("provider") != request["provider"]
            or result.get("input_sha256") != request["input_sha256"]
            or result.get("compute_profile_sha256")
            != request["compute_profile"]["sha256"]
            or len(rows) != len(expected)
            or {r.get("candidate_id") for r in rows} != set(expected)
            or result.get("ranking_weight") != 0
        ):
            raise ArtifactIntegrityError(
                "Algorithm output does not bind every frozen candidate"
            )
        for row in rows:
            original = expected[row["candidate_id"]]
            if (
                row.get("target_id") != frozen.get("context", {}).get("target_id")
                or row.get("input_exact_structure_id")
                != original.get("exact_structure_id")
                or row.get("prepared_exact_structure_id")
                != original.get("prepared_exact_structure_id")
                or row.get("status")
                not in {"complete", "held", "failed", "cancelled", "budget_deferred"}
            ):
                raise ArtifactIntegrityError(
                    "Algorithm candidate input identity/status mismatch"
                )
        receipt = result.get("result_artifact")
        if not isinstance(receipt, dict):
            raise ArtifactIntegrityError(
                "Algorithm result needs its immutable publication"
            )
        raw = self.artifacts.get_bytes(receipt["sha256"])
        payload = {
            k: v
            for k, v in result.items()
            if k not in {"result_artifact", "worker_telemetry"}
        }
        if raw != canonical_json(payload):
            raise ArtifactIntegrityError(
                "Algorithm result differs from its published bytes"
            )

        def verify(value):
            if isinstance(value, dict):
                if "path" in value:
                    raise ArtifactIntegrityError(
                        "Published algorithm artifacts cannot contain local paths"
                    )
                if "sha256" in value and "size_bytes" in value:
                    path = self.artifacts.verified_path(value["sha256"])
                    if (
                        type(value["size_bytes"]) is not int
                        or path.stat().st_size != value["size_bytes"]
                    ):
                        raise ArtifactIntegrityError("Algorithm artifact size mismatch")
                for k, v in value.items():
                    if k != "software":
                        verify(v)
            elif isinstance(value, list):
                for v in value:
                    verify(v)

        # Worker telemetry describes the operator's installed engines and may
        # include their executable paths. It is outside the immutable scientific
        # payload, just as it is excluded from the canonical CAS comparison.
        verify(payload)
        return True
