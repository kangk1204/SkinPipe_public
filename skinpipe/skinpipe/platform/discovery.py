'Frozen discovery runs backed by the PostgreSQL lease/publication ledger.'

from __future__ import annotations

from copy import deepcopy
import json
import math
import sqlite3
import threading
import time
from pathlib import Path

from psycopg.types.json import Jsonb

from .artifacts import canonical_json
from .identity import chemical_identity, digest
from .procurement import has_purchase_record, procurement_state
from .ranking import bind_target_goals, finite, rank_candidates
from .targets import SourceClient, resolve_target, target_request, CHEMBL, now
from .references import measurement, freeze_references
from .catalog_hopping import catalog_hopping, select_docking
from .algorithm_profiles import compute_profile
from .queue import identifier, LeaseLost, Cancelled
from .engines import OwnedProcessCleanupError, ResourceDeferred
from .admet_search import annotate_risk, search_policy
from .native_prepare import run_native, PreparationUnavailable
from .scout_predictions import pin_prediction_metadata
from .scout import require_complete_offers, offer_evidence_status

DEFAULT_FILTERS = {
    "purchasable_only": True,
    "hide_critical": True,
    "hide_opposing": True,
    "hide_severe_admet": True,
    "min_topical": 0,
    "min_pchembl": 0,
    "include_missing_topical": True,
    "include_missing_pchembl": True,
}
ZERO_REASONS = {
    "zero_after_filters",
    "procurement_unconfirmed",
    "source_unavailable",
    "partial_acquisition",
    "structure_preparation_failed",
    "unsupported_chemistry",
    "reference_controls_only",
    "budget_deferred",
}


def freeze_catalog_limit(value=None, *, default=1000):
    value = default if value is None else value
    if type(value) is not int or not 1 <= value <= 10000:
        raise ValueError("catalog_limit must be an integer within 1..10000")
    return value


def freeze_filters(body):
    if not isinstance(body, dict) or set(body) - set(DEFAULT_FILTERS):
        raise ValueError("Unknown discovery filter")
    result = {**DEFAULT_FILTERS, **body}
    for key, value in result.items():
        if key.startswith("min_"):
            if finite(value) is None or not 0 <= value <= 100:
                raise ValueError("Filter thresholds must be finite and nonnegative")
        elif not isinstance(value, bool):
            raise ValueError("Discovery switches must be boolean")
    return result


def filter_reasons(row, filters):
    reasons = []
    if filters["purchasable_only"] and not has_purchase_record(row):
        reasons.append("purchase_record_missing")
    if filters["hide_critical"] and (
        row.get("critical_warning") is True or row.get("hazard_severity") == "critical"
    ):
        reasons.append("critical_warning")
    if filters["hide_opposing"] and row.get("opposing_direction") is True:
        reasons.append("opposing_direction")
    if (
        filters.get("hide_severe_admet")
        and row.get("admet_search_risk", {}).get("excluded") is True
    ):
        reasons.append("severe_predicted_admet")
    for key, source in (("topical", "topical_score"), ("pchembl", "pchembl")):
        value = finite(row.get(source))
        if value is None and not filters["include_missing_" + key]:
            reasons.append(key + "_missing")
        elif value is not None and value < filters["min_" + key]:
            reasons.append(key + "_below_threshold")
    return reasons


def discovery_route(
    reference_count, receptor_valid, completeness, *, identity_status="verified"
):
    if identity_status != "verified":
        return "identity_or_mechanism_selection_required"
    if reference_count:
        return (
            "reference_hopping_then_catalog_docking"
            if receptor_valid
            else "ligand_based_docking_held"
        )
    if receptor_valid:
        return (
            "catalog_docking"
            if completeness == "complete"
            else "data_incomplete_catalog_fallback"
        )
    return "evidence_and_receptor_insufficient"


class DiscoveryRuns:
    def __init__(self, service, config=None, *, source_factory=None, initialize=True):
        self.service, self.queue, self.config = service, service.queue, config or {}
        self.source_factory = source_factory or (
            lambda: SourceClient(service.artifacts)
        )
        self._local_chembl = None
        # Operator connections to an already migrated live database must not
        # acquire migration locks. Normal installation retains its DDL path.
        if not initialize:
            return
        with self.queue.connection() as conn:
            conn.execute("""CREATE TABLE IF NOT EXISTS platform_discovery_runs (
                id TEXT PRIMARY KEY, owner_id TEXT NOT NULL, project_id TEXT NOT NULL,
                request JSONB NOT NULL, request_sha256 TEXT NOT NULL,
                job_id TEXT NOT NULL REFERENCES platform_jobs(id), campaign_id TEXT,
                created_at TIMESTAMPTZ NOT NULL DEFAULT clock_timestamp(),
                UNIQUE(owner_id,request_sha256))""")
            conn.execute(
                "ALTER TABLE platform_discovery_runs ADD COLUMN IF NOT EXISTS desired_state TEXT NOT NULL DEFAULT 'active'"
            )
            conn.execute("""CREATE TABLE IF NOT EXISTS platform_discovery_controls (
                run_id TEXT PRIMARY KEY, owner_id TEXT NOT NULL, desired_state TEXT NOT NULL,
                updated_at TIMESTAMPTZ NOT NULL DEFAULT clock_timestamp())""")

    def resolve(self, body):
        descriptor = self._resolve_target(body)
        receipt = self.service.artifacts.put_bytes(
            canonical_json(descriptor), media_type="application/json"
        )
        return {**descriptor, "resolution_id": receipt["sha256"]}

    def _catalog_source_identity(self):
        sources = getattr(self.service, "catalog_sources", None)
        if sources is not None:
            return sources.identity()
        from .research_catalog import source_identity
        return source_identity(getattr(self.service, "local_catalog", None), getattr(self.service, "scout", {}))

    def _catalog_source_matches(self, frozen):
        return frozen is None or frozen == self._catalog_source_identity()

    def _resolve_input(self, body, *, owner="researcher"):
        registry = getattr(self.service, "research_inputs", None)
        if registry is not None:
            return registry.resolve_for_run({"input": body}, owner=owner)
        from .research_inputs import resolve_input
        return resolve_input(
            body, artifact_store=self.service.artifacts,
            protein_resolver=self._resolve_target,
            source_get=getattr(self, "research_source_get", None),
        )

    def resolve_input(self, body, *, owner="researcher"):
        registry = getattr(self.service, "research_inputs", None)
        if registry is not None:
            return registry.resolve(body, owner=owner)
        descriptor = self._resolve_input(body, owner=owner)
        receipt = self.service.artifacts.put_bytes(
            canonical_json(descriptor), media_type="application/json"
        )
        return {**descriptor, "resolution_id": receipt["sha256"]}

    def _get_local_chembl(self):
        config = self.config.get('chembl_local')
        if not config:
            return None
        if self._local_chembl is None:
            from .chembl_local import LocalChEMBL
            self._local_chembl = LocalChEMBL(config['database'], source_version=config['version'], source_sha256=config['sha256'])
        return self._local_chembl

    def _resolve_target(self, body):
        request = target_request(body)
        client, local, failed_local = self.source_factory(), None, None
        if request['kind'] != 'endpoint' and self.config.get('chembl_local'):
            try:
                local = self._get_local_chembl()
            except (OSError, ValueError, KeyError, TypeError, sqlite3.Error) as error:
                config = self.config['chembl_local']
                failed_local = {'source': 'ChEMBL SQLite', 'version': config.get('version'),
                    'sha256': config.get('sha256'), 'operation': 'local_source_admission',
                    'retrieved_at': now(), 'status': 'source_unavailable', 'raw_sha256': None,
                    'error_type': type(error).__name__}
                client.receipts.append(failed_local)
        descriptor = resolve_target(body, client, local=local)
        if failed_local:
            descriptor['local_target_mapping'] = {'complete': False, 'status': 'source_unavailable',
                'source_provenance': {k: failed_local[k] for k in ('source', 'version', 'sha256')}}
            if descriptor.get('status') == 'verified' and not descriptor.get('chembl_target_id'):
                descriptor['target_mapping_status'] = 'source_unavailable'
        return descriptor

    def submit(self, body, *, owner="researcher"):
        allowed = {
            "target",
            "input",
            "resolution_id",
            "skin_context",
            "search_mode",
            "generation",
            "ai_permissions",
            "snapshot_id",
            "source_target_id",
            "catalog_query",
            "filters",
            "input_id",
            "receptor_request",
            "project_id",
            "idempotency_key",
            "compute_budget_seconds",
            "catalog_limit",
            "screening",
        }
        if not isinstance(body, dict) or set(body) - allowed:
            raise ValueError(
                "Discovery runs accept target and source references, not supplied candidate scores"
            )
        request = deepcopy(body)
        # Client replay tokens may vary after a lost response. Scientific input,
        # run identity and the shared budget do not vary with that token.
        request.pop("idempotency_key", None)
        new_input = "input" in body or "resolution_id" in body
        if new_input or "catalog_limit" in body:
            request["catalog_limit"] = freeze_catalog_limit(body.get("catalog_limit"),
                                                           default=1000 if new_input else 10000)
        if new_input == ("target" in body):
            raise ValueError("Choose one research input or legacy target")
        if "input" in body:
            from .research_inputs import normalize_input
            request["input"] = normalize_input(body["input"])
        elif not new_input:
            request["target"] = target_request(body["target"])
        if new_input or "skin_context" in body:
            from .catalog_assessment import freeze_skin_context
            request["skin_context"] = freeze_skin_context(body.get("skin_context"))
        if new_input:
            mode = body.get("search_mode", "similarity")
            if not isinstance(mode, str) or mode not in {"similarity", "scaffold_hopping", "assessment"}:
                raise ValueError("search_mode must be similarity, scaffold_hopping or assessment")
            from .research_generation import freeze_generation
            request["search_mode"] = mode
            request["generation"] = freeze_generation(body.get("generation"))
            if mode == "assessment" and request["generation"]["enabled"]:
                raise ValueError("Generation requires similarity or scaffold_hopping mode")
            permissions = body.get("ai_permissions", {})
            if (not isinstance(permissions, dict) or set(permissions) - {"public_literature", "research_data"}
                    or any(not isinstance(v, bool) for v in permissions.values())):
                raise ValueError("AI permissions must be explicit boolean switches")
            request["ai_permissions"] = {"public_literature": True, "research_data": False, **permissions}
            from .research_docking import freeze_docking_policy
            request["research_docking_policy"] = freeze_docking_policy(self.config)
            from .research_toxicity import freeze_toxicity
            request["toxicity_policy"] = freeze_toxicity(self.config, self.service.root)
            sources = getattr(self.service, "catalog_sources", None)
            if sources is not None:
                request["catalog_source_binding"], request["catalog_source_identity"] = sources.bind()
            else:
                request["catalog_source_identity"] = self._catalog_source_identity()
        # Preserve the explicit default-species flag without passing it back to
        # the public target input validator on retries.
        request["filters"] = freeze_filters(body.get("filters", {}))
        if new_input:
            # Freeze acquisition semantics independently of scoring profiles.
            # Older stored requests retain their original strict supplier gate.
            request["supplier_evidence_policy"] = "optional_unless_purchase_filter-v1"
        request["admet_search_policy"] = search_policy()
        request["compute_profile"] = compute_profile("generalized-research")
        from .literature import freeze_policy
        literature_policy = freeze_policy(self.config)
        if literature_policy["enabled"] and request.get("ai_permissions", {}).get("public_literature", True):
            request["literature_policy"] = literature_policy
        budget = body.get("compute_budget_seconds", 7200 if new_input else 28800)
        if finite(budget) is None or not 0 < budget <= 28800:
            raise ValueError("Discovery budget must be within 8h actual processing")
        request["compute_budget_seconds"] = budget
        if new_input:
            from .routine_screening import freeze_screening
            policy = freeze_screening(body.get('screening'), compute_budget_seconds=budget)
            if request['search_mode'] != 'assessment':
                request['screening_policy'] = policy
            request.pop('screening', None)
        elif 'screening' in body:
            raise ValueError('Routine screening requires a resolved research input')
        if body.get("snapshot_id"):
            if self.service.catalog.snapshot_owner(body["snapshot_id"]) != owner:
                raise PermissionError("Source snapshot belongs to another project")
            manifest = self.service.catalog.manifest(body["snapshot_id"])
            request["source_manifest_sha256"] = digest(manifest)
            target = body.get("source_target_id")
            if target is None and len(manifest["target_counts"]) == 1:
                target = next(iter(manifest["target_counts"]))
            if target not in manifest["target_counts"]:
                raise ValueError("Choose a source_target_id within the source snapshot")
            request["source_target_id"] = target
        query = request.get("catalog_query", "")
        if not isinstance(query, str) or len(query) > 1024:
            raise ValueError("Catalog query must be a bounded literal string")
        receptor = request.get("receptor_request")
        if receptor is not None:
            import re

            if (
                not isinstance(receptor, dict)
                or set(receptor) - {"pdb_id", "chain", "ligand_code", "assembly_id", "ph", "explicit_pocket", "structural_component_decisions"}
                or not re.fullmatch(r"[0-9][A-Za-z0-9]{3}", receptor.get("pdb_id", ""))
            ):
                raise ValueError(
                    "Receptor request requires a PDB identifier and optional exact chain/pocket"
                )
            if any(
                v is not None
                and (
                    not isinstance(v, str)
                    or not re.fullmatch(r"[A-Za-z0-9_.-]{1,32}", v)
                )
                for k, v in receptor.items() if k not in {"ph", "explicit_pocket", "structural_component_decisions"}
            ):
                raise ValueError("Invalid receptor identifier")
            if "ph" in receptor:
                ph = receptor["ph"]
                if (isinstance(ph, bool) or not isinstance(ph, (int, float))
                        or not math.isfinite(ph) or not 0 <= ph <= 14):
                    raise ValueError("Receptor pH must be a finite number within 0..14")
            if "explicit_pocket" in receptor:
                from ..structure.box import validate_explicit_pocket
                pocket = validate_explicit_pocket(receptor["explicit_pocket"])
                if receptor.get("ligand_code") is not None:
                    raise ValueError("Choose either a cocrystal ligand or an explicit residue pocket")
                if receptor.get("chain") is not None and receptor["chain"] != pocket["chain"]:
                    raise ValueError("Explicit pocket chain differs from the requested chain")
                request["receptor_request"] = {**receptor, "explicit_pocket": pocket}
            if "structural_component_decisions" in receptor:
                from ..structure.cofactors import normalize_component_decisions
                decisions = normalize_component_decisions(receptor["structural_component_decisions"])
                request["receptor_request"] = {**request["receptor_request"], "structural_component_decisions": decisions}
        identifier(owner, "owner")
        project = body.get("project_id") or "discovery-" + digest(owner)[:24]
        identifier(project, "project")
        if new_input:
            # Retry identity is the normalized user request. Resolve a name
            # once, freeze its receipt, and reuse that exact identity on retry.
            input_spec = request.get("input")
            kind = input_spec.get("kind") if isinstance(input_spec, dict) else None
            if kind is None and request.get("resolution_id") is not None:
                registry = getattr(self.service, "research_inputs", None)
                if registry is not None:
                    kind = registry.get(request["resolution_id"], owner=owner).get("kind")
            if kind == "protein":
                # Versioned only on NEW submissions: stored 0.2.5 and legacy
                # requests retain their previous acquisition on later resume.
                from .catalog_diversity import selection_policy
                request["catalog_selection_policy"] = selection_policy()
            input_signature = digest({"owner": owner, "request": request})
            with self.queue.connection() as conn:
                previous = conn.execute(
                    "SELECT id FROM platform_discovery_runs WHERE owner_id=%s AND request->>'input_request_sha256'=%s ORDER BY created_at LIMIT 1",
                    (owner, input_signature),
                ).fetchone()
            if previous is not None:
                return self.get(previous["id"], owner=owner)
            registry = getattr(self.service, "research_inputs", None)
            descriptor = registry.resolve_for_run(body, owner=owner) if registry is not None else self._resolve_input(request["input"], owner=owner)
            if descriptor.get("status") != "verified":
                raise ValueError("Research input requires identity selection before starting a run")
            request["input_request_sha256"] = input_signature
            request["resolved_input"] = descriptor
            request["input"] = deepcopy(descriptor["input"])
            if descriptor.get("kind") in {"compound", "bound_complex"} and (
                body.get("receptor_request") is not None or body.get("input_id") is not None
            ):
                raise ValueError("Structure-seed research cannot override its frozen receptor context or infer a protein")
            if descriptor.get("kind") not in {"compound", "bound_complex"}:
                request["target"] = descriptor["request"]
                from .research_docking import resolve_protein_context
                request["protein_receptor_context"] = resolve_protein_context(
                    self.service, descriptor, request.get("receptor_request"), owner=owner,
                )
        signature = digest({"owner": owner, "request": request})
        run_id = "run-" + signature[:32]
        request["discovery_run_id"] = run_id
        job = self.queue.enqueue(
            "discovery_prepare",
            request,
            owner=owner,
            project=project,
            idempotency_key="discovery-" + signature,
            budget_scope="discovery-" + signature,
            budget_seconds=budget,
        )
        with self.queue.connection() as conn:
            conn.execute(
                """INSERT INTO platform_discovery_runs (id,owner_id,project_id,request,request_sha256,job_id)
                            VALUES (%s,%s,%s,%s,%s,%s) ON CONFLICT DO NOTHING""",
                (run_id, owner, project, Jsonb(request), signature, job["id"]),
            )
        return self.get(run_id, owner=owner)

    def _record(self, run_id, owner):
        identifier(run_id, "run")
        with self.queue.connection() as conn:
            row = conn.execute(
                "SELECT * FROM platform_discovery_runs WHERE id=%s AND owner_id=%s",
                (run_id, owner),
            ).fetchone()
            control = conn.execute("SELECT desired_state FROM platform_discovery_controls WHERE run_id=%s AND owner_id=%s", (run_id, owner)).fetchone()
            if row is not None and control is not None:
                row["desired_state"] = control["desired_state"]
        if row is None:
            raise KeyError("Unknown discovery run")
        if (
            digest(
                {
                    "owner": owner,
                    "request": {
                        k: v
                        for k, v in row["request"].items()
                        if k != "discovery_run_id"
                    },
                }
            )
            != row["request_sha256"]
        ):
            raise ValueError("Discovery request integrity failed")
        return row

    def get(self, run_id, *, owner="researcher"):
        row = self._record(run_id, owner)
        job = self.queue.get(row["job_id"], owner=owner)
        outcome = job.get("result") if job["status"] == "complete" else None
        result = {
            "run_id": run_id,
            "job_id": job["id"],
            "status": job["status"],
            "progress": job.get("progress"),
            "compute_seconds": job["compute_seconds"],
            "compute_budget_seconds": row["request"]["compute_budget_seconds"],
            "filters": row["request"]["filters"],
            "policy_sha256": row["request"]["compute_profile"]["sha256"],
            "campaign_id": row["campaign_id"],
            "scientific_validation": "NOT_ESTABLISHED",
        }
        if "input" in row["request"]:
            result["research_input"] = deepcopy(row["request"]["input"])
            result["skin_context"] = deepcopy(row["request"]["skin_context"])
        result["lifecycle_state"] = row["desired_state"]
        result["preparation_status"] = job["status"]
        if outcome:
            result.update(outcome)
        if row["campaign_id"]:
            campaign = self.service.campaigns._status(row["campaign_id"], owner)
            result["campaign"] = campaign
            result["pipeline_status"] = campaign["lifecycle"]["status"]
            stages = campaign.get("stages", {})
            result["counts"] = {
                **result.get("counts", {}),
                "docking_completed": stages.get("broad_docking", {}).get(
                    "verified_complete", 0
                ),
                "docking_controls_completed": stages.get("broad_docking", {}).get(
                    "verified_controls_complete", 0
                ),
                "confirmation_completed": stages.get("confirmation", {}).get(
                    "verified_complete", 0
                ),
            }
            with self.queue.connection() as conn:
                times = conn.execute(
                    """SELECT coalesce(request->>'stage',request->>'provider',kind) AS stage,
                    sum(compute_seconds) AS seconds FROM platform_jobs
                    WHERE owner_id=%s AND request->>'campaign_id'=%s GROUP BY 1""",
                    (owner, row["campaign_id"]),
                ).fetchall()
            result["processing_seconds_by_stage"] = {
                "preparation": job["compute_seconds"],
                **{r["stage"]: r["seconds"] for r in times},
            }
            result["compute_seconds"] = sum(
                result["processing_seconds_by_stage"].values()
            )
            if campaign.get("published_snapshot_id"):
                result["result_snapshot_id"] = campaign["published_snapshot_id"]
        if row["desired_state"] == "cancelled":
            with self.queue.connection() as conn:
                jobs = conn.execute("""SELECT id,status,cancel_requested,attempts,progress FROM platform_jobs
                    WHERE owner_id=%s AND (id=%s OR (%s::text IS NOT NULL AND request->>'campaign_id'=%s))""",
                    (owner, row["job_id"], row["campaign_id"], row["campaign_id"])).fetchall()
            active = [item["id"] for item in jobs if item["status"] in {"running", "queued"}]
            unverified = [item["id"] for item in jobs if item["status"] != "complete"
                and item["attempts"] > 0 and (item.get('progress', {}).get('resource_quarantine', {}).get('active') is True or not (
                    item.get("progress", {}).get("resource_cleanup", {}).get("confirmed") is True
                    and item.get("progress", {}).get("resource_cleanup", {}).get("attempt") == item["attempts"]))]
            released = not active and not unverified
            result["cancellation"] = {"requested": True,
                "state": "stopping" if active else "stopped" if released else "cleanup_unverified",
                "resources_released": True if released else False if active else None,
                "active_job_ids": active, "cleanup_unverified_job_ids": unverified,
                "resume_allowed": released and result["compute_seconds"] < result["compute_budget_seconds"]}
            result["status"] = "cancelled" if released else "cancelling"
        else:
            result["cancellation"] = {"requested": False, "state": None, "resources_released": None,
                "active_job_ids": [], "resume_allowed": False}
        return result

    def list(self, *, owner="researcher"):
        with self.queue.connection() as conn:
            ids = [
                r["id"]
                for r in conn.execute(
                    "SELECT id FROM platform_discovery_runs WHERE owner_id=%s ORDER BY created_at DESC LIMIT 100",
                    (owner,),
                )
            ]
        return {"runs": [self.get(id, owner=owner) for id in ids]}

    def events(self, run_id, *, owner="researcher", after=0):
        row = self._record(run_id, owner)
        return {
            "run_id": run_id,
            "events": self.queue.events(row["job_id"], owner=owner, after=after),
        }

    def acquire_references(self, descriptor, client, stop):
        rows, receipts, completeness = [], [], "complete"
        tid = descriptor.get("chembl_target_id")
        if not tid:
            status = descriptor.get('target_mapping_status')
            return (
                freeze_references([], descriptor),
                status if status in {'complete', 'partial_acquisition'} else 'source_unavailable',
                [r for r in descriptor.get('sources', []) if r.get('source') == 'ChEMBL SQLite'],
            )
        local = self.config.get("chembl_local")
        try:
            if local:
                db = self._get_local_chembl()
                mapping = db.target_identity(tid)
                expected = {c["accession"] for c in descriptor["components"]}
                if (
                    mapping["accessions"] != sorted(expected)
                    or mapping["taxon_id"] != descriptor["taxon_id"]
                    or mapping['target_type'] != ('PROTEIN COMPLEX' if descriptor['kind'] == 'complex' else 'SINGLE PROTEIN')
                ):
                    raise ValueError(
                        "Pinned local ChEMBL target mapping differs from the resolved target"
                    )
                offset = 0
                while True:
                    stop()
                    page = db.activities_page(tid, limit=10000, offset=offset)
                    rows.extend(
                        measurement(
                            r,
                            descriptor,
                            {
                                k.removeprefix("assay_"): v
                                for k, v in r.items()
                                if k.startswith("assay_")
                            },
                        )
                        for r in page["rows"]
                    )
                    if page["next_offset"] is None:
                        break
                    offset = page["next_offset"]
                receipts.append(
                    {
                        **page["source_provenance"],
                        "target_chembl_id": tid,
                        "retrieved_at": now(),
                        "complete": True,
                    }
                )
                return freeze_references(rows, descriptor), completeness, receipts
            assays, molecules = {}, {}
            offset, total = 0, None
            while total is None or offset < total:
                stop()
                page = client.request(
                    CHEMBL + "/activity.json",
                    {"target_chembl_id": tid, "limit": 1000, "offset": offset},
                )
                values = page.get("activities")
                if not isinstance(values, list) or not isinstance(
                    page.get("page_meta", {}).get("total_count"), int
                ):
                    raise ValueError("Invalid ChEMBL activity pagination")
                total = page["page_meta"]["total_count"]
                if not values and offset < total:
                    raise ValueError("Incomplete ChEMBL page")
                for r in values:
                    aid, mid = r.get("assay_chembl_id"), r.get("molecule_chembl_id")
                    if aid not in assays:
                        stop()
                        assays[aid] = client.request(
                            CHEMBL + "/assay/" + str(aid) + ".json"
                        )
                    if mid not in molecules:
                        stop()
                        molecules[mid] = client.request(
                            CHEMBL + "/molecule/" + str(mid) + ".json"
                        )
                    record = {
                        **r,
                        "canonical_smiles": (
                            molecules[mid].get("molecule_structures") or {}
                        ).get("canonical_smiles"),
                        "target_organism": descriptor.get("organism"),
                    }
                    rows.append(measurement(record, descriptor, assays[aid]))
                offset += len(values)
        except (OSError, ValueError, KeyError, TypeError, sqlite3.Error) as error:
            completeness = "partial_acquisition" if rows else "source_unavailable"
            receipts.append(
                {
                    "status": completeness,
                    "error_type": type(error).__name__,
                    "retrieved_at": now(),
                }
            )
        return (
            freeze_references(rows, descriptor),
            completeness,
            receipts + deepcopy(client.receipts),
        )

    def _reference_index(self, request, references, stop, progress):
        """Union bounded whole-index queries without dropping observed controls.

        The index admits at most five references per request. A protein study
        can have five measured seeds plus its observed redocking ligand. Each
        batch's top N contains every possible global top-N member under the
        maximum-Tanimoto policy. Index IDs are assigned in compound-ID order.
        """
        limit = freeze_catalog_limit(request.get("catalog_limit"), default=10000)
        deadline = time.monotonic() + 300
        batches = []
        for offset in range(0, len(references), 5):
            stop()
            result = self.service._scout_request("/similarity", body={
                "reference_smiles": [r["smiles"] for r in references[offset:offset + 5]],
                "limit": limit, "purchasable_only": request["filters"]["purchasable_only"],
                "hide_severe_admet": request["filters"]["hide_severe_admet"],
                "admet_policy_sha256": request["admet_search_policy"]["sha256"],
            })
            while result.get("status") == "running":
                stop()
                if time.monotonic() > deadline:
                    raise ValueError("Catalog similarity search exceeded acquisition budget")
                time.sleep(0.5)
                result = self.service._scout_request("/queries/" + result["query_id"],
                                                    response_limit=16 * 1024 * 1024)
            if result.get("status") != "complete" or not isinstance(result.get("rows"), dict):
                raise ValueError("Catalog similarity response incomplete")
            receipt = result["rows"]
            ids = [r["compound_id"] for r in receipt["rows"]]
            if len(ids) != len(set(ids)):
                raise ValueError("Duplicate similarity identity")
            batches.append(receipt)
        if len(batches) == 1:
            return batches[0]  # Preserve the historical single-query contract.
        binding = ("index_manifest_sha256", "indexed_structures", "source_report_sha256",
                   "scope", "search_filters", "search_gate_manifest_sha256", "threshold")
        first = batches[0]
        if any(any(batch.get(k) != first.get(k) for k in binding) for batch in batches[1:]):
            raise ValueError("Catalog similarity batches used different index conditions")
        best = {}
        for batch_number, batch in enumerate(batches):
            offset = batch_number * 5
            for row in batch["rows"]:
                score, ref = finite(row.get("tanimoto_2d")), row.get("reference_index")
                if (score is None or not 0 <= score <= 1 or type(ref) is not int
                        or not 0 <= ref < min(5, len(references) - offset)):
                    raise ValueError("Invalid batched similarity observation")
                item = {**row, "reference_index": ref + offset}
                prior = best.get(item["compound_id"])
                if prior is None or (-score, item["reference_index"]) < (-prior["tanimoto_2d"], prior["reference_index"]):
                    best[item["compound_id"]] = item
        rows = sorted(best.values(), key=lambda r: (-r["tanimoto_2d"], r["compound_id"]))[:limit]
        return {**first, "rows": rows, "reference_count": len(references),
                "reference_batches": batches, "aggregation_policy": "max-tanimoto-global-union-v1"}

    def _catalog_rows(self, request, descriptor, stop, progress, reference_set=None):
        retrieval_limit = freeze_catalog_limit(request.get("catalog_limit"), default=10000)
        if request.get("search_mode") == "assessment":
            return [], {"completeness": "not_requested", "scope": "query assessment only",
                        "retrieval_policy": "not_requested", "similarity_index": None}
        if not self._catalog_source_matches(request.get("catalog_source_identity")):
            return [], {"completeness": "source_unavailable", "reason": "frozen_catalog_source_version_changed",
                        "scope": "Frozen catalog version unavailable; no replacement package used",
                        "retrieval_policy": "frozen_source_required", "similarity_index": None,
                        "frozen_catalog_source_identity": deepcopy(request["catalog_source_identity"])}
        source_id = request.get("snapshot_id")
        if source_id:
            manifest = self.service.catalog.manifest(source_id)
            if digest(manifest) != request["source_manifest_sha256"]:
                raise ValueError("Frozen source snapshot changed")
            # Full source population is streamed into the derived snapshot;
            # only a bounded selection pool enters costly chemistry.
            return list(
                self.service.catalog.iter_candidates(
                    source_id, target_id=request["source_target_id"]
                )
            ), {
                "snapshot_id": source_id,
                "manifest_sha256": digest(manifest),
                "completeness": "complete",
                "scope": "frozen source population",
            }
        rows, source, offset, index_receipt = [], None, 0, None
        catalog_selection = None
        indexed_ids = None
        if reference_set and reference_set["computational_seeds"]:
            try:
                index_receipt = self._reference_index(
                    request, reference_set["computational_seeds"], stop, progress)
                indexed_ids = [r["compound_id"] for r in index_receipt["rows"]]
                if len(indexed_ids) != len(set(indexed_ids)):
                    raise ValueError("Duplicate similarity identity")
            except ValueError as error:
                index_receipt = {
                    "status": "held",
                    "reason": str(error),
                    "catalog_fallback": "input" not in request,
                }
                if "input" in request:
                    return [], {"completeness": "source_unavailable", "scope": "similarity retrieval held",
                                "retrieval_policy": "similarity_index_required", "similarity_index": index_receipt,
                                "sampling_limitations": "The configured similarity source was unavailable; no unrelated catalog window was substituted"}
        elif (request.get("input", {}).get("kind") == "protein" and descriptor.get("kind") == "protein"
              and not request.get("catalog_query") and request.get("catalog_selection_policy") is not None):
            from .catalog_diversity import selection_policy, selection_summary
            policy = request["catalog_selection_policy"]
            if policy != selection_policy():
                raise ValueError("Frozen protein diversity policy is unsupported")
            try:
                result = self.service._scout_request("/diversity", body={
                    "limit": retrieval_limit, "seed": policy["seed"], "algorithm": policy["algorithm"],
                    "purchasable_only": request["filters"]["purchasable_only"],
                    "hide_severe_admet": request["filters"]["hide_severe_admet"],
                    "admet_policy_sha256": request["admet_search_policy"]["sha256"],
                })
                deadline = time.monotonic() + 300
                while result.get("status") == "running":
                    stop()
                    if time.monotonic() >= deadline:
                        raise ValueError("Catalog diversity search exceeded acquisition budget")
                    time.sleep(0.2)
                    result = self.service._scout_request("/queries/" + result["query_id"],
                                                       response_limit=16 * 1024 * 1024)
                if result.get("status") != "complete" or not isinstance(result.get("rows"), dict):
                    raise ValueError("Catalog diversity response incomplete")
                index_receipt = result["rows"]
                indexed_ids = [row["compound_id"] for row in index_receipt["rows"]]
                if len(indexed_ids) > retrieval_limit or len(indexed_ids) != len(set(indexed_ids)):
                    raise ValueError("Diversity identities exceed the requested budget or duplicate")
                if index_receipt.get("strategy") != policy["algorithm"]:
                    raise ValueError("Catalog diversity algorithm differs from the frozen policy")
                catalog_selection = selection_summary(index_receipt)
            except (ValueError, OSError) as error:
                index_receipt = {"status": "held", "reason": str(error), "requested_limit": retrieval_limit}
                catalog_selection = selection_summary(index_receipt, status="held", reason=str(error))
                return [], {"completeness": "source_unavailable", "scope": "protein diversity retrieval held",
                            "retrieval_policy": "global_diversity_required", "diversity_index": index_receipt,
                            "catalog_selection": catalog_selection,
                            "sampling_limitations": "The admitted global diversity source was unavailable; no catalog window was substituted"}
        elif "input" in request and not request.get("catalog_query"):
            return [], {"completeness": "source_unavailable", "scope": "similarity retrieval held",
                        "retrieval_policy": "query_structure_required", "similarity_index": None,
                        "sampling_limitations": "No validated query/reference structure was available for catalog similarity"}
        acquisition_error = None
        component_warnings = []
        optional_offers = ("input" in request and request["filters"]["purchasable_only"] is False
            and request.get("supplier_evidence_policy") == "optional_unless_purchase_filter-v1")
        try:
            while offset < retrieval_limit:
                stop()
                response = (
                    {
                        "status": "complete",
                        "rows": [
                            {"compound_id": cid}
                            for cid in indexed_ids[offset : offset + min(1000, retrieval_limit - offset)]
                        ],
                    }
                    if indexed_ids is not None
                    else self.service._scout_query(
                        {
                            "q": request.get("catalog_query", ""),
                            "limit": min(1000, retrieval_limit - offset),
                            "offset": offset,
                            "purchasable_only": request["filters"]["purchasable_only"],
                            "hide_severe_admet": request["filters"][
                                "hide_severe_admet"
                            ],
                            "admet_policy_sha256": request["admet_search_policy"][
                                "sha256"
                            ],
                        }, owner=request.get("_owner", "researcher")
                    )
                )
                deadline = time.monotonic() + 120
                while response.get("status") == "running":
                    stop()
                    if time.monotonic() > deadline:
                        raise ValueError("Catalog query deadline exceeded")
                    time.sleep(0.2)
                    response = self.service._scout_request(
                        "/queries/" + response["query_id"]
                    )
                if response.get("status") != "complete":
                    raise ValueError("Catalog query not complete")
                ids = [r["compound_id"] for r in response["rows"]]
                if len(ids) > min(1000, retrieval_limit - offset):
                    raise ValueError("Catalog response exceeds the requested candidate budget")
                if not ids:
                    break
                materialization_options = {"require_offers": False} if optional_offers else {}
                detail = self.service._scout_materialization_chunk(
                    ids, progress=progress, **materialization_options
                )
                if source is not None and detail["source"] != source:
                    raise ValueError("Catalog source changed during discovery")
                source = detail["source"]
                receipt = self.service.artifacts.put_bytes(canonical_json(detail))[
                    "sha256"
                ]
                for item in detail["rows"]:
                    compound = item["compound"]
                    smiles = compound.get("identity_smiles") or compound.get(
                        "canonical_smiles"
                    )
                    if not optional_offers:
                        require_complete_offers(item)
                    supplier_status = offer_evidence_status(item)
                    if supplier_status != "complete":
                        component_warnings.append({"candidate_id": compound["compound_id"],
                            "component": "offers", "status": supplier_status,
                            "reason_code": "catalog_offers_source_unavailable"})
                    offers = item.get("offers", [])
                    for prediction in item.get("admet", []):
                        for key in ("model_manifest", "endpoint_definitions"):
                            if key in prediction:
                                pin_prediction_metadata(prediction, self.service.artifacts, key, pop_inline=True)
                    rows.append(
                        {
                            "candidate_id": compound["compound_id"],
                            "smiles": smiles,
                            "target_id": descriptor["target_id"],
                            "name": compound.get("name"),
                            "offers": offers,
                            "offers_status": supplier_status,
                            "purchasable": None if supplier_status != "complete" else any(
                                o.get("available") is True
                                or o.get("availability")
                                in {"catalog_reported_in_stock", "catalog_available"}
                                for o in offers
                            ),
                            "catalog_evidence": item,
                            "source_chunk_sha256": receipt,
                            "admet": item.get("admet"),
                            "admet_status": item.get("admet_status"),
                            "admet_recalculation_required": not bool(item.get("admet")),
                        }
                    )
                offset += len(ids)
                progress({"stage": "catalog_acquisition", "acquired": offset})
                if len(ids) < 1000:
                    break
        except (ValueError, OSError) as error:
            acquisition_error = {
                "type": type(error).__name__,
                "acquired": len(rows),
                "reason_code": getattr(error, "reason_code", "catalog_acquisition_failed"),
            }
        if catalog_selection is not None:
            catalog_selection["selected_count"] = len(rows)
            catalog_selection["requested_limit"] = retrieval_limit
            if acquisition_error and not rows:
                catalog_selection.update(status="held", reason="diversity_materialization_unavailable")
        return rows, {
            "source": source,
            "completeness": ("partial_acquisition" if rows else "source_unavailable")
            if acquisition_error
            else "partial_acquisition" if component_warnings
            else "bounded_retrieval"
            if offset >= retrieval_limit
            else "complete",
            "acquisition_error": acquisition_error,
            "component_warnings": component_warnings,
            "scope": "catalog retrieval pool; not whole shared catalog",
            "retrieval_limit": retrieval_limit,
            "similarity_index": index_receipt if catalog_selection is None else None,
            **({"catalog_selection": catalog_selection, "diversity_index": index_receipt}
               if catalog_selection is not None else {}),
            "retrieval_policy": "global_diversity" if catalog_selection is not None else "reference_similarity"
            if indexed_ids is not None
            else "deterministic_catalog_window",
            "sampling_limitations": None
            if indexed_ids is not None
            else "No whole-catalog similarity index; this catalog window is not an exhaustive target-relevant retrieval",
        }

    def execute(self, request, *, owner, stop, progress):
        from contextlib import nullcontext
        from .catalog_cancellation import owned_catalog_scope
        context = getattr(self.service, "catalog_source_context", None)
        if context is not None and request.get("catalog_source_binding") is not None:
            with context(request["catalog_source_binding"]):
                with owned_catalog_scope(self.service, stop=stop) if request.get("_owned_compute") else nullcontext():
                    return self._execute(request, owner=owner, stop=stop, progress=progress)
        with owned_catalog_scope(self.service, stop=stop) if request.get("_owned_compute") else nullcontext():
            return self._execute(request, owner=owner, stop=stop, progress=progress)

    def _execute(self, request, *, owner, stop, progress):
        if "resolved_input" in request:
            descriptor = deepcopy(request["resolved_input"])
        else:
            descriptor_body = {
                k: v for k, v in request["target"].items() if k != "species_defaulted"
            }
            descriptor = self._resolve_target(descriptor_body)
            descriptor["request"]["species_defaulted"] = request["target"][
                "species_defaulted"
            ]
        bound_complex = descriptor.get("kind") == "bound_complex"
        compound_only = descriptor.get("kind") in {"compound", "bound_complex"}
        protein_cpu = "input" in request and descriptor.get("kind") == "protein"
        protein_context = deepcopy(request.get("protein_receptor_context") or {})
        activity_status = "NOT_ESTABLISHED_STRUCTURE_ONLY" if bound_complex else "NOT_ESTABLISHED_NO_PROTEIN_SELECTED"
        descriptor_sha = self.service.artifacts.put_bytes(canonical_json(descriptor))[
            "sha256"
        ]
        if descriptor["status"] != "verified":
            return {
                "route": "identity_or_mechanism_selection_required",
                "target_descriptor": descriptor,
                "target_descriptor_sha256": descriptor_sha,
                "zero_reasons": ["source_unavailable"]
                if descriptor["status"] == "source_unavailable"
                else [],
                "outcome_status": "held",
                "verified_at": now(),
            }
        query_seed = None
        if compound_only:
            source_seed = descriptor.get("query_seed", {})
            seed_smiles = source_seed.get("smiles") or descriptor.get("exact_smiles") or descriptor.get("identity", {}).get("identity_smiles") or descriptor.get("canonical_smiles")
            identity = chemical_identity(seed_smiles)
            query_seed = {
                **deepcopy(source_seed),
                "candidate_id": "query-" + identity["exact_structure_id"].split(":", 1)[1][:32],
                "target_id": descriptor["target_id"],
                "smiles": seed_smiles,
                "name": source_seed.get("name") or descriptor.get("name"),
                "identity": identity,
                "role": "query_seed",
                "measurements": [], "docking": [], "target_association": False,
                "query_seed_status": "observed_complex_ligand_not_activity_reference" if bound_complex else "user_supplied_structure_not_measured_reference",
                "positive_control_status": "not_established",
            }
            policy = {"policy": "user-query-seed-v1", "exact_structure_id": identity["exact_structure_id"]}
            reference_set = {
                "population": [deepcopy(query_seed)],
                "computational_seeds": [deepcopy(query_seed)],
                "observations": [], "sha256": digest(policy), "policy": policy,
            }
            completeness, receipts = "complete", deepcopy(descriptor.get("source_receipts", []))
            literature_evidence = {
                "status": "not_requested", "scope": "compound_only_structure_search",
                "reason": "The query ligand is not a measured activity reference; no protein activity or literature evidence is inferred",
            }
            progress({"stage": "query_seed_frozen", "target_id": descriptor["target_id"]})
        else:
            progress({"stage": "reference_acquisition", "target_id": descriptor["target_id"]})
            from .literature import collect_for_run
            progress({"stage": "literature_evidence", "target_id": descriptor["target_id"]})
            literature_evidence = collect_for_run(
                self.service, descriptor, request, getattr(self, "config", {}), owner=owner, stop=stop
            )
            stop()
            progress({"stage": "reference_acquisition", "target_id": descriptor["target_id"]})
            reference_set, completeness, receipts = self.acquire_references(
                descriptor, self.source_factory(), stop
            )
        comparison_references = deepcopy(reference_set)
        observed_control = None
        if protein_cpu and protein_context.get("status") == "verified":
            bound = protein_context["bound_descriptor"]
            observed_control = {
                **deepcopy(bound["query_seed"]),
                "candidate_id": "observed-redock-" + bound["query_seed"]["identity"]["exact_structure_id"].split(":", 1)[1][:24],
                "role": "control", "control_kind": "observed_redock_reference", "redock_reference": True,
                "target_id": descriptor["target_id"], "measurements": [], "docking": [],
                "target_association": False, "positive_control_status": "not_established",
                "discovery_origin": {"route": "observed_complex_redock_reference", "activity_inherited": False},
            }
            comparison_references["population"].append(deepcopy(observed_control))
            comparison_references["computational_seeds"].append(deepcopy(observed_control))
            comparison_references["sha256"] = digest({"measured_reference_policy": reference_set["sha256"],
                                                        "observed_redock_exact_structure_id": observed_control["identity"]["exact_structure_id"],
                                                        "observed_sdf_sha256": observed_control["observed_sdf_sha256"]})
        ref_sha = self.service.artifacts.put_bytes(canonical_json(reference_set))[
            "sha256"
        ]
        comparison_ref_sha = self.service.artifacts.put_bytes(canonical_json(comparison_references))["sha256"]
        progress({"stage": "catalog_acquisition", "target_id": descriptor["target_id"]})
        try:
            rows, catalog_receipt = self._catalog_rows(
                request, descriptor, stop, progress, comparison_references
            )
        except (ValueError, OSError) as error:
            rows, catalog_receipt = (
                [],
                {
                    "completeness": "source_unavailable",
                    "error_type": type(error).__name__,
                },
            )
        target_id = descriptor["target_id"]
        catalog_count = len(rows)
        from .research_generation import freeze_generation, generate_research
        generation_policy = request.get("generation") or freeze_generation()
        generation_seed = query_seed or next(iter(comparison_references["computational_seeds"]), None)
        if generation_seed is not None and request.get("_owned_compute"):
            from .research_compute import run_compute
            generated_rows, generation_receipt = run_compute("generation", {
                "seed_smiles": generation_seed["smiles"], "policy": generation_policy, "target_id": target_id},
                self.service.root, stop=stop, progress=progress)
        else:
            generated_rows, generation_receipt = generate_research(
                generation_seed["smiles"], generation_policy, target_id, stop=stop,
            ) if generation_seed is not None else ([], {"enabled": generation_policy["enabled"], "status": "held" if generation_policy["enabled"] else "not_requested", "raw_proposal_count": 0, "unique_structure_count": 0, "reason": "explicit_ligand_query_required"})
        rows.extend(generated_rows)
        if query_seed is not None:
            rows.append(query_seed)
            self._enrich_query_rows(generated_rows + [query_seed], stop=stop, source_identity=request.get("catalog_source_identity"))
        elif generated_rows:
            self._enrich_query_rows(generated_rows, stop=stop, source_identity=request.get("catalog_source_identity"))
        if observed_control is not None:
            rows.append(observed_control)
            self._enrich_query_rows([observed_control], stop=stop, source_identity=request.get("catalog_source_identity"))
        for row in rows:
            if row.get("target_id") != target_id or compound_only:
                # A catalog structure can be reused across targets; activity,
                # docking scores and old target hypotheses cannot.
                row["source_target_observations"] = {
                    "target_id": row.get("target_id"),
                    "measurements": row.get("measurements"),
                    "docking": row.get("docking"),
                    "role": row.get("role"),
                    "control_kind": row.get("control_kind"),
                    "positive_control_status": row.get("positive_control_status"),
                }
                row["measurements"], row["docking"] = [], []
                for key in (
                    "primary_docking",
                    "pchembl",
                    "pose",
                    "dock_affinity",
                    "docking_execution",
                    "opposing_direction",
                    "goal_support",
                    "requested_action",
                    "direction",
                    "direction_verdict",
                    "activity_status",
                    "curated_claim_status",
                    "computation_status",
                    "target_goal_binding",
                ):
                    row.pop(key, None)
                row["target_association"] = False
                row["target_evidence"] = []
                if compound_only and row.get("role") != "query_seed":
                    row["role"] = "candidate"
                    row.pop("control_kind", None)
                    row.pop("positive_control_status", None)
            row["target_id"] = target_id
            row["identity"] = chemical_identity(row.get("smiles"))
            row["discovery_run_id"] = request["discovery_run_id"]
            row["priority_primary_protocol_id"] = None
            row["research_search_mode"] = request.get("search_mode")
            row["procurement"] = procurement_state(row)
            if compound_only:
                row["research_input"] = deepcopy(request["input"])
                row["activity_evidence_status"] = activity_status
                row["query_identity_match"] = row["identity"]["exact_structure_id"] == query_seed["identity"]["exact_structure_id"]
        measured_by_exact = {}
        for record in reference_set["observations"]:
            exact = record.get("identity", {}).get("exact_structure_id")
            if exact:
                measured_by_exact.setdefault(exact, []).append(record)
        for row in rows:
            exact = chemical_identity(row.get("smiles"))["exact_structure_id"]
            if exact in measured_by_exact:
                row["measurements"] = measured_by_exact[exact]
                row["target_association"] = True
        from ..discovery import reweight_records
        from .catalog_assessment import assess_catalog

        def assess(rows):
            if not request.get("_owned_compute"):
                return assess_catalog(rows, skin_context=request.get("skin_context"), stop=stop, progress=progress)
            from .research_compute import run_compute
            assessed = run_compute("assessment", {"rows": rows, "skin_context": request.get("skin_context")},
                self.service.root, stop=stop, progress=progress)
            # Preserve shared control-row references held by the coordinator.
            for original, updated in zip(rows, assessed, strict=True):
                original.clear()
                original.update(updated)
            return rows

        rows = assess(rows)
        rows = rank_candidates(
            bind_target_goals(
                reweight_records(rows), desired_action=descriptor.get("desired_action", "any")
            )
        )
        for row in rows:
            annotate_risk(row, policy=request["admet_search_policy"])
            reasons = filter_reasons(row, request["filters"])
            if "input" in request:
                if row["identity"]["structure_status"] != "valid":
                    reasons.append("invalid_structure")
                elif row["identity"]["parent_ambiguous"]:
                    reasons.append("ambiguous_multicomponent_structure")
            reasons.extend(row.get("generation_filter_reasons", []))
            row["discovery_filter_reasons"], row["discovery_filter_pass"] = (
                reasons,
                not reasons,
            )
            row["discovery_origin"] = {
                "route": "generated_proposal" if row.get("derivation") else "query_seed" if row.get("role") == "query_seed" else "general_catalog",
                "reference_policy_sha256": reference_set["sha256"],
                "activity_inherited": False,
            }
        # The full retrieved population is preserved; only the cheapest top
        # 10,000 eligible rows enter expensive reference comparisons.
        pool = [r for r in rows if r["discovery_filter_pass"] and r.get("role", "candidate") == "candidate"][:10000]
        hop_settings = dict(shape_limit=(generation_policy["max_3d"] if request.get("search_mode") == "scaffold_hopping" else 0)
            if "input" in request else 2000,
            conformers=generation_policy["conformers"],
            conformer_timeout_seconds=10 if "input" in request else 0)
        if request.get("_owned_compute"):
            from .research_compute import run_compute
            pool, hop_receipt = run_compute("hopping", {"rows": pool, "reference_set": comparison_references,
                "artifact_root": str(self.service.artifacts.root), "settings": hop_settings},
                self.service.root, stop=stop, progress=progress)
        else:
            pool, hop_receipt = catalog_hopping(pool, comparison_references, self.service.artifacts,
                stop=stop, progress=progress, **hop_settings)
        by_id = {r["candidate_id"]: r for r in pool}
        rows = [by_id.get(r["candidate_id"], r) for r in rows]
        if catalog_receipt.get("catalog_selection") is not None:
            # Hopping intentionally resets imported origin claims. Restore
            # acquisition provenance only from THIS run's server-owned policy
            # and receipt, after that sanitizing boundary and before saving.
            from .catalog_diversity import selection_policy
            if request.get("catalog_selection_policy") != selection_policy():
                raise ValueError("Diversity origin requires the frozen server selection policy")
            selection = catalog_receipt["catalog_selection"]
            for row in rows:
                if row.get("derivation") or row.get("role", "candidate") != "candidate":
                    continue
                row["discovery_origin"] = {
                    "route": "general_catalog_diversity", "reference_policy_sha256": reference_set["sha256"],
                    "activity_inherited": False, "catalog_selection_policy": request["catalog_selection_policy"]["version"],
                    "target_activity": "NOT_ESTABLISHED_NO_REFERENCE",
                    "source_population": selection.get("source_population"),
                    "source_pool_size": selection.get("source_pool_size"),
                    "index_manifest_sha256": selection.get("index_manifest_sha256"),
                }
                row["activity_evidence_status"] = "NOT_ESTABLISHED_NO_REFERENCE"
        control_ids, reference_controls = [], []
        if observed_control is not None:
            observed_control = next(r for r in rows if r.get("redock_reference") is True)
            reference_controls.append(observed_control)
            control_ids.append(observed_control["candidate_id"])
        for reference in ([] if compound_only else reference_set["computational_seeds"]):
            cid = (
                "ref-"
                + reference["identity"]["exact_structure_id"].split(":", 1)[1][:32]
            )
            control = {
                **deepcopy(reference),
                "candidate_id": cid,
                "role": "control",
                "control_kind": "measured_reference",
                "positive_control_status": "not_established",
                "discovery_run_id": request["discovery_run_id"],
                "discovery_origin": {
                    "route": "frozen_measured_reference",
                    "activity_inherited": False,
                },
                "discovery_filter_pass": True,
                "discovery_filter_reasons": [],
                "priority_primary_protocol_id": None,
                "procurement": procurement_state(reference),
            }
            rows.append(control)
            reference_controls.append(control)
            control_ids.append(cid)
        measured_controls = [r for r in reference_controls if r.get("control_kind") == "measured_reference"]
        if "input" in request and measured_controls:
            self._enrich_query_rows(measured_controls, stop=stop, source_identity=request.get("catalog_source_identity"))
            assess(measured_controls)
            for control in measured_controls:
                annotate_risk(control, policy=request["admet_search_policy"])
        toxicity_receipt = None
        if "input" in request:
            from .research_toxicity import apply_toxicity
            rows, toxicity_receipt = apply_toxicity(rows, self.service, request, stop=stop, progress=progress)
        # Reserve measured-reference parents before spending the candidate
        # quota. Catalog aliases of a reference remain in the saved population.
        selection_rows = pool + reference_controls + ([query_seed] if query_seed else [])
        if request.get('screening_policy') is not None:
            from .routine_screening import eligibility, validate_policy
            routine_policy = validate_policy(request['screening_policy'])
            # Additional prediction attaches to new row objects. Use the
            # current source rows, while preserving the full saved population.
            current = {r['candidate_id']: r for r in rows}
            selection_rows = [current.get(r['candidate_id'], r) for r in selection_rows]
            selection_rows = [r for r in selection_rows if eligibility(r)[0] != 'held'
                              and 'exact_model_predictions_required' not in eligibility(r)[1]]
            selected = select_docking(selection_rows, maximum=routine_policy['max_new_docking'])
        else:
            selected = select_docking(selection_rows)
        docking_receipt = None
        if bound_complex and request.get("search_mode") != "assessment":
            from .research_docking import evaluate_bound_complex
            rows, docking_receipt = evaluate_bound_complex(
                rows, descriptor, selected, self.service, self.config, request,
                stop=stop, progress=progress,
            )
        elif protein_cpu:
            if request.get("search_mode") == "assessment":
                docking_receipt = {"status": "not_requested", "receptor_valid": None,
                                   "prepared_count": 0, "selected_count": 0, "completed_count": 0}
            elif protein_context.get("status") != "verified":
                docking_receipt = {"status": "held", "receptor_valid": False,
                                   "hold_reasons": protein_context.get("hold_reasons") or ["explicit_pdb_and_observed_ligand_selection_required"],
                                   "prepared_count": 0, "selected_count": 0, "completed_count": 0}
            else:
                from .research_docking import evaluate_bound_complex
                rows, docking_receipt = evaluate_bound_complex(
                    rows, protein_context["bound_descriptor"], selected,
                    self.service, self.config, request, stop=stop, progress=progress,
                )
        routine_selection = None
        if request.get('screening_policy') is not None:
            from .routine_screening import select_representatives
            if descriptor.get('kind') == 'compound':
                for row in rows:
                    row['routine_route'] = 'ligand_only'
            routine_selection = select_representatives(rows, request['screening_policy'])
            for row in rows:
                row['routine_disposition'] = routine_selection['decisions'][row['candidate_id']]
        snapshot_id = (
            "discovery-"
            + digest(
                {
                    "run": request["discovery_run_id"],
                    "descriptor": descriptor_sha,
                    "references": ref_sha,
                    "rows": rows,
                    "source": catalog_receipt,
                }
            )[:24]
        )
        provenance = {
            "adapter": "generalized-discovery-v1",
            "owner_id": owner,
            "discovery_run_id": request["discovery_run_id"],
            "target_descriptor": descriptor,
            "target_descriptor_sha256": descriptor_sha,
            "actual_protein_id": descriptor.get("actual_protein_id") if compound_only else target_id,
            "taxon_id": descriptor.get("taxon_id"),
            "target_id": target_id,
            "references": reference_set["population"],
            "reference_set_sha256": ref_sha,
            "comparison_reference_set_sha256": comparison_ref_sha,
            "reference_control_ids": control_ids,
            "selected_candidate_ids": selected,
            "reference_policy_sha256": reference_set["sha256"],
            "source_receipts": receipts,
            "source_completeness": completeness,
            "catalog_receipt": catalog_receipt,
            "catalog_selection": deepcopy(catalog_receipt.get("catalog_selection")),
            "catalog_selection_policy": deepcopy(request.get("catalog_selection_policy")),
            "catalog_source_identity": deepcopy(request.get("catalog_source_identity")),
            "catalog_source_binding": deepcopy(request.get("catalog_source_binding")),
            "filters": request["filters"],
            "admet_search_policy": request["admet_search_policy"],
            "compute_profile": request["compute_profile"],
            **({'routine_selection': routine_selection} if routine_selection is not None else {}),
            "hop_receipt": hop_receipt,
            "literature_evidence": deepcopy(literature_evidence),
            "research_input": deepcopy(request.get("input")),
            "skin_context": deepcopy(request.get("skin_context")),
            "analysis_mode": "bound_complex" if bound_complex else "compound_only" if compound_only else "protein_target",
            "query_seed_id": query_seed["candidate_id"] if query_seed else None,
            "activity_inherited": False,
            "search_mode": request.get("search_mode"),
            "generation_receipt": generation_receipt,
            "additional_toxicity_receipt": toxicity_receipt,
            "bound_complex_docking": docking_receipt,
            "protein_receptor_context": protein_context if protein_cpu else None,
            "observed_redock_control_id": observed_control["candidate_id"] if observed_control else None,
            "measured_reference_control_ids": [r["candidate_id"] for r in reference_controls if r.get("control_kind") == "measured_reference"],
        }
        stop()
        try:
            manifest = self.service.catalog.create_snapshot(
                rows, snapshot_id, provenance
            )
        except FileExistsError:
            manifest = self.service.catalog.manifest(snapshot_id)
        progress(
            {
                "stage": "candidate_population_frozen",
                "snapshot_id": manifest["snapshot_id"],
                "catalog_count": catalog_count,
            }
        )
        if compound_only or protein_cpu:
            input_id, validation = None, {
                "receptor_valid": (docking_receipt or {}).get("receptor_valid"),
                "prepared_count": (docking_receipt or {}).get("prepared_count", 0),
                "status": "explicit_cpu_research_path" if protein_cpu else "not_requested_compound_only",
                "reason": "Explicit CPU research uses the frozen observed complex and fixed pocket; no legacy GPU campaign is activated"
                if protein_cpu else "Compound structure search and skin assessment require no surrogate protein or docking receptor",
            }
        else:
            input_id, validation = self.prepare_inputs(
                request, descriptor, manifest,
                selected + control_ids if selected else [], stop, progress,
            )
        literature_evidence["run_verification_links"] = {
            "scope": "Actual checks already performed by this new run; separate from literature interpretation",
            "target_descriptor_sha256": descriptor_sha,
            "target_identity_status": descriptor["status"],
            "reference_set_sha256": ref_sha,
            "reference_acquisition_completeness": completeness,
            "input_id": input_id,
            "input_validation": deepcopy(validation),
            "pose_geometry": "NOT_EVALUATED_BY_LITERATURE_REVIEW",
            "ranking_performance": "NOT_EVALUATED_BY_LITERATURE_REVIEW",
            "experimental_validation": "NOT_EVALUATED",
        }
        catalog_rows = [r for r in rows if r.get("role", "candidate") == "candidate"]
        receptor_valid = validation.get("receptor_valid") is True
        reasons = []
        if not catalog_count and not generated_rows and request.get("search_mode") != "assessment":
            reasons.append(
                "source_unavailable"
                if catalog_receipt["completeness"] != "complete"
                else "zero_after_filters"
            )
        elif catalog_rows and not any(r["discovery_filter_pass"] for r in catalog_rows):
            reasons.append("zero_after_filters")
        elif catalog_rows and not selected:
            control_parents = {r["identity"]["parent_id"] for r in reference_controls}
            eligible = [r for r in pool if r.get("role", "candidate") == "candidate"]
            reasons.append(
                "reference_controls_only"
                if eligible and all(r["identity"]["parent_id"] in control_parents for r in eligible)
                else "unsupported_chemistry"
            )
        if catalog_rows and not any(
            r["procurement"]["documented_exact_form"] for r in catalog_rows
        ):
            reasons.append("procurement_unconfirmed")
        if completeness != "complete":
            reasons.append(completeness)
        if catalog_receipt["completeness"] in {"source_unavailable", "partial_acquisition", "bounded_retrieval"}:
            reasons.append(catalog_receipt["completeness"])
        if not receptor_valid and selected and not (compound_only or protein_cpu):
            reasons.append("structure_preparation_failed")
        result = {
            "additional_toxicity_receipt": toxicity_receipt,
            "target_descriptor": descriptor,
            "target_descriptor_sha256": descriptor_sha,
            "snapshot_id": manifest["snapshot_id"],
            "literature_evidence": literature_evidence,
            "target_id": target_id,
            "input_id": input_id,
            "input_validation": validation,
            "route": (("bound_complex_catalog_hopping" if request.get("search_mode") == "scaffold_hopping" else "bound_complex_catalog_similarity")
                      if bound_complex and (catalog_count or generated_rows) else "bound_complex_seed_assessment" if bound_complex else
                      ("compound_only_catalog_hopping" if request.get("search_mode") == "scaffold_hopping" else "compound_only_catalog_similarity")
                      if catalog_count or generated_rows else "compound_only_seed_assessment") if compound_only else (
                      "protein_catalog_research_cpu_docking" if protein_context.get("status") == "verified"
                      else "protein_catalog_research_docking_held") if protein_cpu else discovery_route(
                len(reference_set["population"]), receptor_valid, completeness
            ),
            "source_completeness": completeness,
            "catalog_completeness": catalog_receipt["completeness"],
            "catalog_selection": deepcopy(catalog_receipt.get("catalog_selection")),
            "reference_set_sha256": ref_sha,
            "reference_policy_sha256": reference_set["sha256"],
            "counts": {
                **({'representative_goal': routine_selection['policy']['target_candidates'],
                    'representative_selected': routine_selection['selected_count'],
                    'representative_shortfall': routine_selection['shortfall'],
                    'routine_held': routine_selection['counts'].get('held', 0),
                    'routine_pending': routine_selection['counts'].get('pending', 0),
                    'routine_alternatives': routine_selection['counts'].get('alternative', 0)}
                   if routine_selection is not None else {}),
                "collected": catalog_count,
                "generated": len(generated_rows),
                "raw_generated_proposals": generation_receipt["raw_proposal_count"],
                "saved_rows": len(rows),
                "new_unique_candidates": len({r.get("identity", {}).get("exact_structure_id")
                    for r in catalog_rows if r.get("identity", {}).get("exact_structure_id")
                    and (query_seed is None or not r.get("query_identity_match"))}),
                "reference_controls": sum(r.get("control_kind") == "measured_reference" for r in reference_controls),
                "observed_redock_controls": int(observed_control is not None),
                "purchase_record": sum(has_purchase_record(r) for r in catalog_rows),
                "filter_pass": sum(r["discovery_filter_pass"] for r in catalog_rows),
                "docking_selected": (docking_receipt or {}).get("selected_count", 0) if compound_only or protein_cpu else len(selected),
                "comparison_selected": len(selected) if "input" in request else 0,
                "query_seeds": int(compound_only),
                "severe_predicted_admet": sum(
                    r.get("admet_search_risk", {}).get("excluded") is True
                    for r in catalog_rows
                ),
                "prepared": validation.get("prepared_count", 0),
                "docking_completed": (docking_receipt or {}).get("completed_count", 0),
                "query_redocking_completed": int((docking_receipt or {}).get("redocking", {}).get("status") == "complete"),
                "documented_procurement": sum(
                    r["procurement"]["documented_exact_form"] for r in catalog_rows
                ),
                "references": 0 if compound_only else len(reference_set["population"]),
                "reference_seeds": 0 if compound_only else len(reference_set["computational_seeds"]),
            },
            "selected_candidate_ids": selected,
            **({'routine_selection': routine_selection} if routine_selection is not None else {}),
            "hop_receipt": hop_receipt,
            "generation_receipt": generation_receipt,
            "bound_complex_docking": docking_receipt,
            "protein_receptor_context": protein_context if protein_cpu else None,
            "zero_reasons": sorted(set(reasons)),
            "outcome_status": "completed_zero"
            if catalog_receipt["completeness"] == "complete"
            and not any(r["discovery_filter_pass"] for r in catalog_rows)
            else "review_available",
            "verified_at": now(),
            "scientific_validation": "NOT_ESTABLISHED",
        }
        if compound_only:
            result["research_input"] = deepcopy(request["input"])
            result["skin_context"] = deepcopy(request["skin_context"])
            result["seed_assessment"] = next(r for r in rows if r.get("role") == "query_seed")
            result["outcome_status"] = "bound_complex_review_available" if bound_complex else "compound_only_review_available" if catalog_count or generated_rows else "query_assessment_available"
            result["activity_evidence_status"] = activity_status
        return result

    def _enrich_query_rows(self, rows, *, stop, source_identity=None):
        """Attach only exact DB matches. Missing predictions remain missing."""
        lookup = getattr(self.service, "lookup_compounds", None)
        if not self._catalog_source_matches(source_identity):
            for row in rows:
                row["admet_status"] = "unassessed_frozen_database_version_changed"
                row["exact_database_lookup"] = {"status": "unavailable", "matches": [], "reason": "frozen_catalog_source_version_changed"}
            return
        if lookup is None:
            for row in rows:
                row.setdefault("admet_status", "unassessed_database_not_configured")
            return
        stop()
        try:
            results = lookup([r["smiles"] for r in rows])
        except (ValueError, OSError):
            results = []
        indexed = {r["query_smiles"]: r for r in results}
        for row in rows:
            stop()
            result = indexed.get(row["smiles"], {"status": "unavailable", "matches": []})
            row["exact_database_lookup"] = deepcopy(result)
            exact = chemical_identity(row["smiles"])["exact_structure_id"]
            matches = []
            for match in result.get("matches", []):
                compound = match.get("compound", {})
                smiles = compound.get("identity_smiles") or compound.get("canonical_smiles")
                if smiles and chemical_identity(smiles)["exact_structure_id"] == exact:
                    matches.append(match)
            if not matches:
                row["admet_status"] = "unassessed_database_miss" if result.get("status") == "not_found" else "unassessed_database_unavailable"
                row.setdefault("admet", [])
                continue
            row["offers"] = [offer for match in matches for offer in match.get("offers", [])]
            row["admet"] = [prediction for match in matches for prediction in match.get("admet", [])]
            row["admet_status"] = "available_exact_database_match" if row["admet"] else "unassessed_model_not_configured"
            row["admet_recalculation_required"] = not bool(row["admet"])
            row["purchasable"] = any(o.get("available") is True or o.get("availability") in {"catalog_reported_in_stock", "catalog_available"} for o in row["offers"])

    def prepare_inputs(self, request, descriptor, manifest, selected, stop, progress):
        # Registration reruns native structure/provenance checks. Receptor
        # context can be transferred only by exact target/taxon, never label.
        api, owner = (
            self.service.campaigns,
            request["_owner"] if "_owner" in request else "researcher",
        )
        if not selected:
            return None, {
                "receptor_valid": False,
                "reason": "no_selected_candidates",
                "prepared_count": 0,
            }
        if descriptor["kind"] == "complex":
            return None, {
                "receptor_valid": False,
                "reason": "full_complex_assembly_and_engine_support_required",
                "prepared_count": 0,
            }
        receptor_record = None
        if (request.get("input_id") and request.get("snapshot_id")
                and request.get("receptor_request") is None):
            prior = api._inputs(
                request["input_id"],
                request["snapshot_id"],
                request["source_target_id"],
                owner,
            )
            receptor = prior.get("receptor")
            if (
                receptor
                and receptor["actual_protein_id"] == descriptor["target_id"]
                and receptor["taxon_id"] == descriptor["taxon_id"]
            ):
                receptor_record = {
                    k: receptor[k]
                    for k in ("identity_card_sha256", "uniprot_record_sha256")
                }
        receptor_settings = request.get("receptor_request")
        if receptor_settings is None and descriptor.get("structure_choices"):
            receptor_settings = {"pdb_id": descriptor["structure_choices"][0]["pdb_id"]}
        if receptor_record is None and receptor_settings:
            accession = descriptor["components"][0]["accession"]
            work = self.service.root / "discovery-work" / request["discovery_run_id"]
            work.mkdir(parents=True, exist_ok=True)
            settings = dict(receptor_settings)
            try:
                stop()
                rec = run_native(
                    {"kind": "receptor", "settings": settings, "accession": accession},
                    work,
                    python=self.config.get("preparation_python"),
                    timeout=360,
                    stop=stop,
                )
                card_path = Path(rec["card_path"])
                card = json.loads(card_path.read_bytes())
                fields = card["receptor"]
                if fields.get("pocket_source") == "centroid" or any(
                    "BOX LIGAND AMBIGUOUS" in n for n in fields.get("notes", [])
                ):
                    return None, {
                        "receptor_valid": False,
                        "reason": "explicit_pocket_selection_required",
                        "prepared_count": 0,
                        "notes": fields.get("notes", []),
                    }
                for artifact in card["artifacts"].values():
                    if artifact.get("present") is True:
                        raw = (card_path.parent / artifact["path"]).read_bytes()
                        if (
                            self.service.artifacts.put_bytes(raw)["sha256"]
                            != artifact["sha256"]
                        ):
                            raise ValueError("Native receptor artifact changed")
                receptor_record = {
                    "identity_card_sha256": self.service.artifacts.put_bytes(
                        card_path.read_bytes()
                    )["sha256"],
                    "uniprot_record_sha256": descriptor["components"][0][
                        "source_record_sha256"
                    ],
                }
            except (ValueError, OSError, PreparationUnavailable) as error:
                return None, {
                    "receptor_valid": False,
                    "reason": str(error)
                    if isinstance(error, PreparationUnavailable)
                    else type(error).__name__,
                    "prepared_count": 0,
                }
        if receptor_record is None:
            return None, {
                "receptor_valid": False,
                "reason": "explicit_verified_receptor_or_pocket_selection_required",
                "prepared_count": 0,
            }
        try:
            checked_receptor = api._receptor(receptor_record, manifest, descriptor["target_id"])
        except ValueError as error:
            return None, {
                "receptor_valid": False,
                "reason": str(error),
                "prepared_count": 0,
            }
        prepared, held = {}, []
        preparation_ph = checked_receptor.get("preparation_ph", 7.4)
        work = self.service.root / "discovery-work" / request["discovery_run_id"]
        for index in range(0, len(selected), 32):
            stop()
            ids = selected[index : index + 32]
            candidates = self.service.catalog.get_candidates(
                manifest["snapshot_id"], ids, descriptor["target_id"]
            )
            items = [
                {
                    "candidate_id": cid,
                    "smiles": candidate["smiles"],
                }
                for cid, candidate in zip(ids, candidates)
            ]
            try:
                ligand_request = {"kind": "ligands", "ligands": items}
                if preparation_ph != 7.4:
                    ligand_request["ph"] = preparation_ph
                result = run_native(
                    ligand_request,
                    work,
                    python=self.config.get("preparation_python"),
                    timeout=360,
                    stop=stop,
                )
                diagnostic_log_sha = None
                if any(row.get("reason") for row in result["ligands"]) and result.get("diagnostic_log"):
                    log = result["diagnostic_log"]
                    diagnostic_log_sha = self.service.artifacts.put_file(log["path"])["sha256"]
                    if diagnostic_log_sha != log["sha256"]:
                        raise ValueError("Native preparation diagnostic log changed")
                for row in result["ligands"]:
                    cid = row["candidate_id"]
                    if row.get("reason"):
                        held.append({"candidate_id": cid, "reason": row["reason"],
                            "error_type": row.get("error_type", "PreparationUnavailable"),
                            "original_structure_id": row.get("original_structure_id") or "raw:" + digest(next(item["smiles"] for item in items if item["candidate_id"] == cid)),
                            "diagnostic_log_sha256": diagnostic_log_sha})
                        continue
                    self.service.artifacts.put_file(row["path"])
                    prepared[cid] = self.service.artifacts.put_file(
                        row["receipt_path"]
                    )["sha256"]
            except (ValueError, OSError, PreparationUnavailable) as error:
                held.extend(
                    {"candidate_id": item["candidate_id"], "reason": "native_preparation_runtime_failed",
                     "error_type": type(error).__name__, "original_structure_id": "raw:" + digest(item["smiles"]),
                     "diagnostic_log_sha256": None} for item in items
                )
            progress(
                {
                    "stage": "ligand_preparation",
                    "processed": min(index + 32, len(selected)),
                    "requested": len(selected),
                    "complete": len(prepared),
                }
            )
        native = {
            "schema": "skinpipe-campaign-inputs-v1",
            "snapshot_id": manifest["snapshot_id"],
            "target_id": descriptor["target_id"],
            "prepared": prepared,
            "preparation_holds": held,
            "receptor": receptor_record,
        }
        sha = self.service.artifacts.put_bytes(canonical_json(native))["sha256"]
        checked = api.register_inputs({"manifest_sha256": sha}, owner=owner)
        return checked["input_id"], {
            "receptor_valid": True,
            "prepared_count": len(prepared),
            "held": held,
            "receptor_interpretation": checked_receptor.get("interpretation"),
        }

    def continue_run(self, run_id, *, owner="researcher"):
        self._record(run_id, owner)


        with self.queue.connection() as conn:
            row = conn.execute(
                "SELECT * FROM platform_discovery_runs WHERE id=%s AND owner_id=%s FOR UPDATE",
                (run_id, owner),
            ).fetchone()
            control = conn.execute("SELECT desired_state FROM platform_discovery_controls WHERE run_id=%s AND owner_id=%s", (run_id, owner)).fetchone()
            if control:
                row["desired_state"] = control["desired_state"]
            if row["desired_state"] == "active":
                job = self.queue.get(row["job_id"], owner=owner)
                if job["status"] == "complete" and job.get("result", {}).get(
                    "input_id"
                ):
                    if not row["campaign_id"]:
                        result = job["result"]
                        remaining = (
                            row["request"]["compute_budget_seconds"]
                            - job["compute_seconds"]
                        )
                        if remaining > 0:
                            campaign = self.service.campaigns.post(
                                "/v2/campaigns",
                                {
                                    "snapshot_id": result["snapshot_id"],
                                    "target_id": result["target_id"],
                                    "input_id": result["input_id"],
                                    "compute_profile": "generalized-research",
                                    "desired_action": row["request"].get("target", {}).get("desired_action", "any"),
                                    "compute_budget_seconds": remaining,
                                },
                                owner=owner,
                            )
                            row["campaign_id"] = campaign["campaign_id"]
                            conn.execute(
                                "UPDATE platform_discovery_runs SET campaign_id=%s WHERE id=%s",
                                (row["campaign_id"], run_id),
                            )
                    if row["campaign_id"]:
                        control = conn.execute("SELECT desired_state FROM platform_discovery_controls WHERE run_id=%s AND owner_id=%s", (run_id, owner)).fetchone()
                        if control is None or control["desired_state"] == "active":
                            self.service.campaigns.post(
                                "/v2/campaigns/" + row["campaign_id"] + "/run", {}, owner=owner)
                        # An intent arriving during activation is durably
                        # fenced before this coordinator returns.
                        control = conn.execute("SELECT desired_state FROM platform_discovery_controls WHERE run_id=%s AND owner_id=%s", (run_id, owner)).fetchone()
                        if control and control["desired_state"] == "cancelled":
                            self.service.campaigns.cancel_run_jobs(row["campaign_id"], owner=owner)
        return self.get(run_id, owner=owner)

    def set_lifecycle(self, run_id, state, *, owner="researcher"):
        self._record(run_id, owner)
        if state not in {"active", "cancelled"}:
            raise ValueError("Unknown run lifecycle")
        with self.queue.connection() as conn:
            # Control intent uses an independent row: slow preparation/activation
            # cannot turn a valid cancellation into a five-second lock timeout.
            conn.execute("""INSERT INTO platform_discovery_controls(run_id,owner_id,desired_state) VALUES(%s,%s,%s)
                ON CONFLICT(run_id) DO UPDATE SET desired_state=EXCLUDED.desired_state,updated_at=clock_timestamp()
                WHERE platform_discovery_controls.owner_id=EXCLUDED.owner_id""", (run_id, owner, state))
        return self._record(run_id, owner)

    def run_once(self, *, shutdown_event=None):
        from .resources import probe

        telemetry = probe(directory=self.service.root)
        capabilities = ["discovery_prepare", "result_export"]
        if getattr(getattr(self.service, "personal_literature", None), "configured", False):
            capabilities.append("literature_review")
        telemetry["capabilities"] = capabilities
        self.queue.register_node("discovery-coordinator", telemetry)
        job = self.queue.claim(
            "discovery-coordinator", capabilities
        )
        if not job:
            return None
        if not self.queue.reserve_cpu(job["id"], job["lease_token"], 1):
            return self.queue.defer(job["id"], job["lease_token"], "cpu_capacity_wait")
        started, initial = time.monotonic(), job["compute_seconds"]
        ending, lost = threading.Event(), threading.Event()
        heartbeat_state = {"error_type": None, "terminal": None, "cancel_requested": False}
        progress = {"stage": "identity_resolution"}

        def cleanup_diagnostic(stage, error_type):
            # Diagnostics must not replace the original hold or certify exit.
            try:
                record = getattr(self.queue, 'record_cleanup_diagnostic', None)
                if record is not None:
                    record(job['id'], job['lease_token'], worker='discovery-coordinator',
                           stage=stage, error_type=error_type)
            except Exception:
                pass

        def elapsed():
            return initial + time.monotonic() - started

        def stop():
            if shutdown_event is not None and shutdown_event.is_set():
                raise ResourceDeferred("personal_service_shutdown")
            if lost.is_set():
                if heartbeat_state["cancel_requested"]:
                    raise Cancelled("Cancellation requested")
                if heartbeat_state["error_type"]:
                    raise LeaseLost(
                        "Discovery ownership unverified after heartbeat "
                        + heartbeat_state["error_type"]
                    )
                raise LeaseLost("Discovery ownership lost")
            if elapsed() >= job["budget_seconds"]:
                raise TimeoutError("budget_deferred")

        def heartbeat():
            while not ending.wait(1):
                try:
                    state = self.queue.heartbeat(
                        job["id"],
                        job["lease_token"],
                        progress=dict(progress),
                        compute_seconds=elapsed(),
                    )
                    if state.get("status") != "running":
                        # Budget accounting may atomically close the lease.
                        # Its terminal receipt is authoritative; do not write
                        # another heartbeat/failure using that closed token.
                        if state.get("status") == "budget_deferred":
                            heartbeat_state["terminal"] = dict(state)
                        lost.set()
                        return
                except Cancelled:
                    heartbeat_state["cancel_requested"] = True
                    lost.set()
                    return
                except LeaseLost:
                    lost.set()
                    return
                except Exception as error:
                    # A dead keeper cannot certify ownership. Fence the local
                    # actor at its next checkpoint; retain only the error class
                    # because connection exceptions can contain private DSNs.
                    heartbeat_state["error_type"] = type(error).__name__
                    lost.set()
                    return

        thread = threading.Thread(target=heartbeat, daemon=True)
        thread.start()
        cleanup_known = True
        try:
            request = {**job["request"], "_owner": job["owner_id"], "_owned_compute": True}
            if job["kind"] == "literature_review":
                result = self.service.personal_literature.execute(
                    request, owner=job["owner_id"], stop=stop, progress=progress.update)
            elif job["kind"] == "result_export":
                result = self.service.exports.execute(request, stop=stop, progress=progress.update)
            else:
                result = self.execute(request, owner=job["owner_id"], stop=stop, progress=progress.update)
            stop()
            ending.set()
            thread.join()
            # An in-flight heartbeat can fail during keeper teardown, after
            # the first checkpoint. Recheck its fence before publication.
            stop()
            from .engines import cleanup_receipt
            publication = self.queue.finish(
                job["id"], job["lease_token"], result, compute_seconds=elapsed(), cleanup=cleanup_receipt()
            )
            return publication
        except Cancelled as error:
            ending.set()
            thread.join()
            from .engines import cleanup_receipt, release_stopped_exception_frames
            release_stopped_exception_frames(error)
            progress.clear()
            try:
                return self.queue.acknowledge_cancel(job["id"], job["lease_token"],
                    compute_seconds=elapsed(), cleanup=cleanup_receipt())
            except (LeaseLost, ValueError):
                return {"status": "cleanup_unverified"}
        except LeaseLost:
            return heartbeat_state["terminal"] or {"status": "lease_lost"}
        except OwnedProcessCleanupError:
            cleanup_known = False
            # An unknown native/query lifetime must never trigger an automatic
            # duplicate calculation or a cleanup acknowledgement.
            ending.set()
            thread.join()
            cleanup_diagnostic('owned_cleanup', 'OwnedProcessCleanupError')
            try:
                self.queue.fail(job["id"], job["lease_token"],
                    "owned_resource_cleanup_unverified", status="needs_data")
            except (LeaseLost, Cancelled):
                pass
            return {"status": "cleanup_unverified"}
        except ResourceDeferred as error:
            ending.set()
            thread.join()
            if lost.is_set():
                return heartbeat_state["terminal"] or {"status": "lease_lost"}
            reason = str(error)
            try:
                state = self.queue.heartbeat(
                    job["id"], job["lease_token"],
                    progress={"stage": "resource_wait", "reason": reason},
                    compute_seconds=elapsed(),
                )
                if state.get("status") != "running":
                    return state
                self.queue.defer(job["id"], job["lease_token"], reason)
                return {"status": "resource_deferred", "reason": reason}
            except (LeaseLost, Cancelled):
                return {"status": "lease_lost"}
        except TimeoutError:
            ending.set()
            thread.join()
            if lost.is_set():
                return heartbeat_state["terminal"] or {"status": "lease_lost"}
            try:
                state = self.queue.heartbeat(
                    job["id"],
                    job["lease_token"],
                    progress={"stage": "budget_deferred"},
                    compute_seconds=elapsed(),
                )
                if state.get("status") != "running":
                    return state
                return self.queue.fail(
                    job["id"],
                    job["lease_token"],
                    "budget_deferred",
                    status="budget_deferred",
                )
            except (LeaseLost, Cancelled):
                return {"status": "lease_lost"}
        except Exception as error:
            ending.set()
            thread.join()
            if lost.is_set():
                return heartbeat_state["terminal"] or {"status": "lease_lost"}
            try:
                state = self.queue.heartbeat(
                    job["id"],
                    job["lease_token"],
                    progress={"stage": "failed", "error_type": type(error).__name__},
                    compute_seconds=elapsed(),
                )
                if state.get("status") != "running":
                    return state
                return self.queue.fail(
                    job["id"], job["lease_token"], type(error).__name__,
                    status="failed" if job["kind"] == "literature_review" else "retryable"
                )
            except (LeaseLost, Cancelled):
                return {"status": "lease_lost"}
        finally:
            ending.set()
            thread.join()
            acknowledge = getattr(self.queue, 'acknowledge_cleanup', None)
            if cleanup_known and acknowledge is not None:
                from .engines import cleanup_receipt, release_stopped_exception_frames
                import sys
                stopped_error = sys.exc_info()[1]
                if stopped_error is not None:
                    release_stopped_exception_frames(stopped_error)
                try:
                    acknowledge(job['id'], job['lease_token'], worker='discovery-coordinator', cleanup=cleanup_receipt())
                except (LeaseLost, ValueError) as error:
                    cleanup_diagnostic('cleanup_ack', 'LeaseLost' if isinstance(error, LeaseLost) else 'ValueError')

    def reconcile(self):
        'Restart-safe continuation is independent of preparation publication.'
        outcomes = []



        with self.queue.connection() as conn:
            stopping = conn.execute("""SELECT r.id,r.owner_id,r.job_id,r.campaign_id
                FROM platform_discovery_runs r
                LEFT JOIN platform_discovery_controls c ON c.run_id=r.id AND c.owner_id=r.owner_id
                WHERE COALESCE(c.desired_state,r.desired_state)='cancelled' AND EXISTS(
                    SELECT 1 FROM platform_jobs j WHERE j.owner_id=r.owner_id
                    AND (j.id=r.job_id OR (r.campaign_id IS NOT NULL AND j.request->>'campaign_id'=r.campaign_id))
                    AND j.status IN ('queued','running') AND NOT j.cancel_requested)
                ORDER BY r.created_at LIMIT 10""").fetchall()
        for row in stopping:
            with self.queue.connection() as conn:
                unfenced = conn.execute("""SELECT id FROM platform_jobs WHERE owner_id=%s
                    AND (id=%s OR (%s::text IS NOT NULL AND request->>'campaign_id'=%s))
                    AND status IN ('queued','running') AND NOT cancel_requested ORDER BY created_at,id LIMIT 100""",
                    (row['owner_id'], row['job_id'], row['campaign_id'], row['campaign_id'])).fetchall()
            failures = []
            for item in unfenced:
                try:
                    self.queue.cancel(item['id'], owner=row['owner_id'])
                except Exception as error:
                    failures.append({'job_id': item['id'], 'error_type': type(error).__name__})
            if row['campaign_id']:
                try:
                    pause = self.service.campaigns.cancel_run_jobs(row['campaign_id'], owner=row['owner_id'])
                    failures.extend(pause.get('errors', []))
                except Exception as error:
                    failures.append({'step': 'campaign_pause', 'error_type': type(error).__name__})
            outcomes.append({'run_id': row['id'], 'cancellation_reconciled': True,
                             'status': 'retry_pending' if failures else 'cancel_requested', 'errors': failures})
        with self.queue.connection() as conn:
            rows = conn.execute("""SELECT r.id,r.owner_id FROM platform_discovery_runs r
                JOIN platform_jobs j ON j.id=r.job_id
                LEFT JOIN platform_discovery_controls c ON c.run_id=r.id AND c.owner_id=r.owner_id
                WHERE j.status='complete'
                AND j.result->>'input_id' IS NOT NULL AND r.campaign_id IS NULL
                AND COALESCE(c.desired_state,r.desired_state)='active'
                ORDER BY r.created_at LIMIT 10""").fetchall()
        for row in rows:
            try:
                outcomes.append(self.continue_run(row["id"], owner=row["owner_id"]))
            except Exception as error:
                outcomes.append(
                    {
                        "run_id": row["id"],
                        "continuation": "retry_pending",
                        "error_type": type(error).__name__,
                    }
                )
        with self.queue.connection() as conn:
            active = conn.execute("""SELECT r.id,r.owner_id,r.campaign_id FROM platform_discovery_runs r
                LEFT JOIN platform_discovery_controls c ON c.run_id=r.id AND c.owner_id=r.owner_id
                WHERE r.campaign_id IS NOT NULL AND COALESCE(c.desired_state,r.desired_state)='active'
                ORDER BY r.created_at LIMIT 100""").fetchall()
        for row in active:
            try:
                if not self._reconcile_active(row['id'], row['owner_id']):
                    continue
                coordinator = self.service.campaigns._coordinator(row["owner_id"])
                broad = (
                    coordinator.status(row["campaign_id"])
                    .get("stages", {})
                    .get("broad_docking")
                )
                if (
                    broad
                    and broad.get(
                        "computable_missing_or_corrupt", broad["missing_or_corrupt"]
                    )
                    == 0
                ):
                    # Cancellation is independent of slow status/asset checks.
                    # Recheck immediately before enqueue and fence any jobs
                    # created by an intent arriving during that enqueue.
                    if not self._reconcile_active(row['id'], row['owner_id']):
                        continue
                    result = self.service.algorithm_jobs.submit(
                        row["campaign_id"],
                        {"provider": "posebusters"},
                        owner=row["owner_id"],
                        automatic=True,
                    )
                    if not self._reconcile_active(row['id'], row['owner_id']):
                        for job in result.get('jobs', []):
                            self.queue.cancel(job['id'], owner=row['owner_id'])
                        outcomes.append({'run_id': row['id'], 'pose_quality': 'cancelled'})
                        continue
                    outcomes.append(
                        {"run_id": row["id"], "pose_quality": result["status"]}
                    )
            except (ValueError, KeyError, OSError) as error:
                outcomes.append(
                    {
                        "run_id": row["id"],
                        "pose_quality": "held",
                        "error_type": type(error).__name__,
                    }
                )
        return outcomes

    def _reconcile_active(self, run_id, owner):
        with self.queue.connection() as conn:
            control = conn.execute(
                'SELECT desired_state FROM platform_discovery_controls WHERE run_id=%s AND owner_id=%s',
                (run_id, owner),
            ).fetchone()
        return control is None or control['desired_state'] == 'active'


def main():
    import argparse
    import signal
    from .cli import read_config, service_from_config

    parser = argparse.ArgumentParser(
        description="Resume SkinPipe generalized discovery jobs"
    )
    parser.add_argument("--config", required=True)
    parser.add_argument("--once", action="store_true")
    args = parser.parse_args()
    service = service_from_config(read_config(args.config))
    end = threading.Event()
    for sig in (signal.SIGINT, signal.SIGTERM):
        signal.signal(sig, lambda *_: end.set())
    while not end.is_set():
        try:
            service.discovery.run_once()
            service.discovery.reconcile()
        except Exception as error:
            print(
                json.dumps(
                    {"status": "retry_pending", "error_type": type(error).__name__}
                ),
                flush=True,
            )
        if args.once:
            break
        end.wait(2)


if __name__ == "__main__":
    main()
