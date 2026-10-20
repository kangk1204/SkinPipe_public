'Validate the selected structure and retain explicit calculation provenance.'
from __future__ import annotations

from copy import deepcopy

from .artifacts import ArtifactIntegrityError
from .campaign import SUCCESS, _sealed_read
from .identity import candidate_identity, digest
from .ranking import COMPARABILITY_VERSION, PRIORITY_POLICY_VERSION, comparable_protocol, rank_candidates


def publish_campaign_results(coordinator, catalog, campaign_id, *, owner=None):
    """Return a content-idempotent, whole-population derived snapshot manifest.

    Retained failed/quarantined observations never become scores. Every accepted
    pose, prepared input and receptor is checksum-verified again at publication.
    The source snapshot and receipt ledger are unchanged, including controls.
    """
    with coordinator._lock(campaign_id) as directory:
        campaign = coordinator._manifest(campaign_id)
        source_id, target = campaign["plan"]["snapshot_id"], campaign["plan"]["target_id"]
        source = catalog.manifest(source_id)
        source_owner = catalog.snapshot_owner(source_id)
        if owner is not None and source_owner != owner:
            raise PermissionError("Campaign publication source belongs to another project")
        if target not in source["target_counts"]:
            raise ValueError("Campaign target is absent from its source snapshot")
        state = _sealed_read(directory / "state.json")
        expected = {r["candidate_id"]: candidate_identity(r, trust_stored=True)["exact_structure_id"]
                    for r in campaign["candidates"]}
        diagnostics = {r["candidate_id"]: deepcopy(r["preparation_diagnostic"])
                       for r in campaign["candidates"] if r.get("preparation_diagnostic")}
        observations, receipts, quarantine = {}, [], []
        expected_observations = 0
        for stage_name in sorted(state["stages"]):
            spec = coordinator._stage_read(directory, state, stage_name)
            expected_observations += len(spec["seeds"]) * len(spec["candidate_ids"] + spec["control_ids"])
            for seed in spec["seeds"]:
                task = coordinator._task_key(spec, seed)
                for cid, item in sorted(state["outcomes"].get(task, {}).items()):
                    if cid not in expected or cid not in spec["candidate_ids"] + spec["control_ids"]:
                        raise ArtifactIntegrityError("Campaign outcome has an unexpected candidate")
                    receipts.append({"candidate_id": cid, "stage": stage_name, "seed": seed,
                                     "receipt_sha256": item["receipt_sha256"]})
                    try:
                        receipt = coordinator._read_outcome(directory, item)
                        if receipt["campaign_id"] != campaign_id or receipt["candidate_id"] != cid or receipt["task_key"] != task:
                            raise ArtifactIntegrityError("Campaign receipt binding mismatch")
                        raw = deepcopy(receipt["row"])
                        protocol = raw["protocol"]
                        if (raw["candidate_id"] != cid or raw["input_sha256"] != campaign["prepared"][cid]["pdbqt_sha256"]
                                or protocol["receptor_sha256"] != spec["receptor"]["receptor_sha256"]):
                            raise ArtifactIntegrityError("Campaign receipt structure or receptor mismatch")
                        successful = raw["status"] in SUCCESS
                        scope = state.get("confirmation_scope")
                        historical_confirmation = bool(stage_name == "confirmation" and scope and (
                            cid not in scope["candidate_ids"] + scope["control_ids"] or seed not in scope["seeds"]))
                        valid_pose = raw["status"] in {"complete", "review_required"} and bool(raw.get("pose_sha256"))
                        observations.setdefault(cid, []).append({
                            **raw, "target_id": target, "receptor_id": protocol["receptor_sha256"],
                            "receptor_sha256": protocol["receptor_sha256"],
                            "receptor_identity_card_sha256": spec["receptor"].get("identity_card_sha256"),
                            "receptor_interpretation": deepcopy(spec["receptor"].get("interpretation")),
                            "pocket_id": digest([protocol["center"], protocol["size"]]),
                            "protocol_id": digest(comparable_protocol(protocol)), "unit": "kcal/mol",
                            "comparison_policy": COMPARABILITY_VERSION,
                            "execution_protocol_id": digest(protocol),
                            "score": raw.get("affinity") if successful else None,
                            "higher_is_better": False, "artifact_valid": valid_pose,
                            "ranking_eligible": successful and not historical_confirmation,
                            "historical_confirmation_policy": historical_confirmation,
                            "source_stage": stage_name, "source_seed": seed, "source_campaign_id": campaign_id,
                            "source_receipt_sha256": item["receipt_sha256"],
                            "source_job_id": receipt["job_id"], "exact_structure_id": expected[cid],
                            "receipt_runtime": receipt.get("runtime"),
                            "receipt_worker_telemetry": receipt.get("worker_telemetry"),
                            "scientific_validation": "NOT_ESTABLISHED"})
                    except (ArtifactIntegrityError, OSError, ValueError, KeyError, TypeError) as error:
                        quarantine.append({"candidate_id": cid, "source_stage": stage_name,
                                           "source_receipt_sha256": item["receipt_sha256"],
                                           "status": "quarantined", "reason": str(error)})
        if not receipts:
            raise ValueError("Campaign has no collected observations to publish")

        pool = {}
        for row in catalog.iter_candidates(source_id, target_id=target):
            cid = row["candidate_id"]
            if cid not in expected:
                continue
            if candidate_identity(row, trust_stored=True)["exact_structure_id"] != expected[cid]:
                raise ArtifactIntegrityError("Frozen candidate structure differs from source snapshot")
            previous = row.get("docking") or []
            previous = [previous] if isinstance(previous, dict) else list(previous)
            row["docking"] = [r for r in previous if r.get("source_campaign_id") != campaign_id] + observations.get(cid, [])
            row["campaign_result_quarantine"] = [q for q in quarantine if q["candidate_id"] == cid]
            if observations.get(cid):
                row["computation_status"] = ("completed" if any(r["status"] == "complete" for r in observations[cid])
                                             else "review_required" if any(r["artifact_valid"] for r in observations[cid])
                                             else "failed")
            elif row["campaign_result_quarantine"]:
                row["computation_status"] = "quarantined"
            elif cid in diagnostics:
                row["preparation_diagnostic"] = diagnostics[cid]
                row["computation_status"] = "preparation_held"
            row["campaign_computation_status"] = (row["computation_status"]
                if observations.get(cid) or row["campaign_result_quarantine"] or cid in diagnostics
                else "not_calculated")
            row["campaign_result_source"] = campaign_id
            broad = [(index, observation) for index, observation in enumerate(row["docking"])
                     if observation.get("source_campaign_id") == campaign_id and observation.get("source_stage") == "broad_docking"]
            if broad:
                index, primary = min(broad, key=lambda item: (item[1]["source_seed"], item[1]["protocol_id"]))
                row["primary_docking"] = {key: primary.get(key) for key in
                    ("score", "status", "engine", "source_stage", "protocol_id", "receptor_sha256", "pose_sha256")}
                row["primary_docking"].update(raw_affinity=primary.get("affinity"), seed=primary["source_seed"],
                                              observation_index=index)
            pool[cid] = row




        control_ids = set(source.get("provenance", {}).get("reference_control_ids", []))
        if (set(pool) != set(expected) or not control_ids <= set(pool)
                or any(pool[cid].get("role") not in {
                    "control", "positive_control", "negative_control", "counter_control", "reference_control"
                } for cid in control_ids)
                or len(set(pool) - control_ids) > 10000):
            raise ArtifactIntegrityError("Campaign source population is missing or exceeds reranking budget")
        if source.get("provenance", {}).get("discovery_run_id"):
            # Select the frozen broad stage and engine, never confirmation or
            # sensitivity protocols. Multiple execution cohorts are explicit;
            # the planned primary engine cannot silently become another engine.
            broad_spec = coordinator._stage_read(directory, state, "broad_docking")
            primary_ids = {o["protocol_id"] for rs in observations.values() for o in rs
                           if o.get("source_stage") == "broad_docking" and o.get("engine") == broad_spec["engine"]}
            if len(primary_ids) > 1:
                raise ArtifactIntegrityError("Primary scientific protocol differs within frozen broad stage")
            primary_id = next(iter(primary_ids), None)
            for row in pool.values():
                row["priority_primary_protocol_id"] = primary_id
        ranked = {r["candidate_id"]: r for r in rank_candidates(list(pool.values()))}
        primary_protocols = {}
        for row in ranked.values():
            primary = row.get("primary_docking")
            if primary is None:
                continue
            entry = primary_protocols.setdefault(primary["protocol_id"], {
                "target_id": target, **{key: primary[key] for key in ("protocol_id", "engine", "seed", "receptor_sha256")},
                "candidate_count": 0, "scored_candidate_count": 0})
            entry["candidate_count"] += 1
            entry["scored_candidate_count"] += int(primary["score"] is not None)
        provenance = {**deepcopy(source.get("provenance", {})), "adapter": "campaign-results-v2",
                      "parent_snapshot_id": source_id, "parent_manifest_sha256": digest(source),
                      "campaign_id": campaign_id, "target_id": target,
                      "confirmation_amendment": deepcopy(state.get("confirmation_amendment")),
                      "confirmation_execution_scope": deepcopy(state.get("confirmation_scope")),
                      "receipts": receipts, "quarantined_receipts": quarantine,
                      "primary_docking_protocols": [primary_protocols[key] for key in sorted(primary_protocols)],
                      "reranking_candidate_count": len(ranked), "policy_version": PRIORITY_POLICY_VERSION,
                      "comparison_policy": COMPARABILITY_VERSION,
                      "scientific_validation": "NOT_ESTABLISHED"}
        received = len(receipts)
        collection_complete = ("confirmation" in state["stages"]
                               and received == expected_observations and not quarantine
                               and not campaign["held"]
                               and all(o["ranking_eligible"] for rs in observations.values() for o in rs
                                       if not o.get("historical_confirmation_policy")))
        provenance["publication_scope"] = {
            "scope": "registered campaign observation collection; not efficacy qualification",
            "status": "COMPLETE_REGISTERED_COLLECTION" if collection_complete else "PARTIAL_OR_HELD",
            "registered_collection_complete": collection_complete,
            "whole_workflow_complete": False,
            "whole_workflow_status": "NOT_ESTABLISHED",
            "registered_stages": sorted(state["stages"]),
            "expected_seeded_observations": expected_observations,
            "collected_observations": received,
            "uncollected_observations": max(0, expected_observations - received),
            "quarantined_observations": len(quarantine),
            "preparation_held_candidates": len(campaign["held"]),
        }
        if source_owner is not None:
            provenance["owner_id"] = source_owner
        identifier = "campaign-" + digest(provenance)[:24]
        if (catalog.root / identifier).exists():
            catalog.verify(identifier)
            return catalog.manifest(identifier)

        def rows():
            for row in catalog.iter_candidates(source_id):
                if row["target_id"] == target and row["candidate_id"] in ranked:
                    row = ranked[row["candidate_id"]]
                    row["priority_score_status"] = "recomputed_in_campaign_pool"
                yield row

        result = catalog.create_snapshot(rows(), snapshot_id=identifier, provenance=provenance)
        if result["candidate_count"] != source["candidate_count"] or result["target_counts"] != source["target_counts"]:
            raise ArtifactIntegrityError("Campaign publication changed candidate population")
        catalog.verify(identifier)
        return result
