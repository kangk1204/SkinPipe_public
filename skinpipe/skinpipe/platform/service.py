"""Version two application service shared by the legacy HTTP bridge and CLI.

The bridge authenticates the research workspace before dispatch. Worker tokens
are scoped to worker routes and never grant access to experiment uploads.
"""

from __future__ import annotations

import base64
import hashlib
import ipaddress
import json
import re
from pathlib import Path
import threading
import time
from urllib.parse import parse_qs, quote, unquote, urlparse
from urllib import error as urlerror, request as urlrequest
from copy import deepcopy

from .artifacts import ArtifactIntegrityError, ArtifactStore, atomic_bytes, canonical_json, sha256_file
from .catalog import Catalog
from .campaign_api import CampaignAPI
from .experiments import ExperimentStore, MAX_UPLOAD_BYTES
from .feedback import measured_feedback
from .identity import chemical_identity
from .queue import Queue
from .policy import assess_candidate
from .ranking import (
    PRIORITY_POLICY_VERSION,
    bind_target_goals,
    rank_candidates,
    validate_requested_action,
)
from .selection import select_experiment_batch
from .scorers import OptionalScorers
from .correspondence import compare_structures
from .algorithm_jobs import AlgorithmJobs
from .algorithm_profiles import PROVIDERS, compute_profile
from .scout_predictions import pin_prediction_metadata
from .scout import require_complete_offers


class _CatalogRedirectRejected(urlerror.URLError):
    """A credential boundary rejection is not a transient transport failure."""


class _NoRedirect(urlrequest.HTTPRedirectHandler):
    def redirect_request(self, request, fp, code, msg, headers, newurl):
        raise _CatalogRedirectRejected("Scout gateway redirects are disabled")


class _CatalogTransportError(ValueError):
    """A bounded read may be retried without creating another gateway job."""

    reason_code = "gateway_transport_unavailable"


def _snapshot_list_summary(snapshot):
    """Project verified list metadata without mutating the sealed manifest."""
    result = {
        key: deepcopy(snapshot[key])
        for key in (
            "schema_version", "snapshot_id", "created_at", "candidate_count",
            "target_counts", "status", "reason", "error", "ownership_status",
        )
        if key in snapshot
    }
    provenance = snapshot.get("provenance")
    if isinstance(provenance, dict):
        result["provenance"] = {
            key: deepcopy(provenance[key])
            for key in (
                "marker", "target_id", "screen_id", "search_mode",
                "primary_docking_protocols",
            )
            if key in provenance
        }
    return result


class PlatformService:
    def __init__(
        self,
        root,
        dsn,
        *,
        bridge=None,
        worker_tokens=None,
        scout=None,
        discovery_config=None,
        access_config=None,
    ):
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)
        self.catalog = Catalog(self.root / "catalog")
        self.artifacts = ArtifactStore(self.root / "artifacts")
        self.queue = Queue(dsn)
        self.experiments = ExperimentStore(dsn, artifact_root=self.root / "experiments")
        self.campaigns = CampaignAPI(
            self.root / "campaign-api", self.catalog, self.artifacts, self.queue
        )
        self.algorithm_jobs = AlgorithmJobs(self.campaigns, self.artifacts, self.queue)
        self.optional_scorers = OptionalScorers(self.root / "optional", {})
        self.bridge = bridge
        self.worker_tokens = worker_tokens or {}
        self.scout = scout or {}
        self.local_catalog = None
        if self.scout.get("kind") == "local":
            from .local_catalog import LocalCatalog
            self.local_catalog = LocalCatalog(self.scout)
        elif self.scout:
            parsed = urlparse(self.scout.get("url", ""))
            try:
                address = ipaddress.ip_address(parsed.hostname or "")
                private_node = address in ipaddress.ip_network("100.64.0.0/10")
                parsed.port
            except ValueError:
                private_node = False
            if (
                parsed.scheme not in {"http", "https"}
                or not private_node
                or parsed.username
                or parsed.password
                or parsed.path not in {"", "/"}
                or parsed.query
                or parsed.fragment
                or not isinstance(self.scout.get("token"), str)
                or not self.scout["token"]
            ):
                raise ValueError(
                    "Scout gateway requires a configured Tailscale IP URL and private bearer token"
                )
            self.scout = {
                "url": self.scout["url"].rstrip("/"),
                "token": self.scout["token"],
            }
        self.lock = threading.Lock()
        from .research_catalog import CatalogSources
        self.catalog_sources = CatalogSources(self)
        self._materializations = {}
        self._legacy_imports = {}
        from .discovery import DiscoveryRuns

        self.discovery = DiscoveryRuns(self, discovery_config)
        from .personal_literature import PersonalLiterature

        self.personal_literature = PersonalLiterature(self, discovery_config)
        from .access import ProjectAccess

        self.access = ProjectAccess(self.queue, access_config)
        from .research_structures import ResearchStructures
        from .research_registry import ResearchInputs
        self.research_structures = ResearchStructures(self.root / "research-structures", artifacts=self.artifacts)
        self.research_inputs = ResearchInputs(self)
        from .research_assistance import ResearchAssistance
        self.research_assistance = ResearchAssistance(self)
        # The adapter's registered paths are not a project ownership grant.
        # Operators bind allowed source IDs to owner namespaces explicitly.
        self.legacy_source_grants = (access_config or {}).get("legacy_sources", {})
        if (not isinstance(self.legacy_source_grants, dict) or any(
            not isinstance(owner, str) or not owner or not isinstance(ids, list)
            or any(not isinstance(sid, str) or not 1 <= len(sid) <= 256 for sid in ids)
            for owner, ids in self.legacy_source_grants.items())):
            raise ValueError("Legacy sources require explicit owner-to-screen-ID grants")
        from .exports import RunExports

        self.exports = RunExports(self)

    @staticmethod
    def _open_scout(request):
        # Neither HTTP redirects nor inherited proxies may forward the private
        # node credential to a different destination.
        opener = urlrequest.build_opener(urlrequest.ProxyHandler({}), _NoRedirect())
        return opener.open(request, timeout=25)

    def _scout_request(self, path, *, body=None, response_limit=4 * 1024 * 1024):
        from .catalog_cancellation import note_catalog_query, note_catalog_transport_failure
        local, scout = self._catalog_runtime()
        try:
            result = self._scout_request_impl(path, body=body, response_limit=response_limit)
        except BaseException:
            # A failed submission response may still have admitted a query.
            # Without its ID, the research actor cannot certify cleanup or
            # retry the computation as if nothing had been started.
            if body is not None and (local is not None or scout):
                note_catalog_transport_failure(path)
            raise
        note_catalog_query(path, result)
        return result

    def _scout_request_impl(self, path, *, body=None, response_limit=4 * 1024 * 1024):
        local, scout = self._catalog_runtime()
        if local is not None:
            result = local.request(path, body)
            if len(canonical_json(result)) > response_limit:
                raise ValueError("Local catalog response exceeds its bounded transfer limit")
            return result
        if not scout:
            raise ValueError("Global catalog gateway is not configured")
        request = urlrequest.Request(
            scout["url"] + path,
            data=canonical_json(body) if body is not None else None,
            headers={
                "Authorization": "Bearer " + scout["token"],
                "Content-Type": "application/json",
                "Accept": "application/json",
            },
        )
        try:
            with self._open_scout(request) as response:
                raw = response.read(response_limit + 1)
            if len(raw) > response_limit:
                raise ValueError(
                    "Global catalog response exceeds its bounded transfer limit"
                )

            def reject_constant(value):
                raise ValueError("Global catalog returned a nonfinite JSON value")

            return json.loads(raw, parse_constant=reject_constant)
        except urlerror.HTTPError as exc:
            raise ValueError(
                f"Global catalog gateway returned HTTP {exc.code}"
            ) from None
        except _CatalogRedirectRejected:
            raise ValueError("Scout gateway redirects are disabled") from None
        except (urlerror.URLError, TimeoutError, OSError):
            raise _CatalogTransportError(
                "Global catalog gateway is temporarily unavailable"
            ) from None

    def lookup_compounds(self, smiles_list):
        local, scout = self._catalog_runtime()
        if local is not None:
            return self._pin_local_lookup_predictions(local.lookup(smiles_list))
        if not isinstance(smiles_list, list) or len(smiles_list) > 1000:
            raise ValueError("At most 1000 exact compound structures per lookup")
        if any(not isinstance(value, str) or not 1 <= len(value) <= 4096
               or any(ord(char) < 32 for char in value) for value in smiles_list):
            raise ValueError("Exact lookup requires bounded structure strings")
        if not scout:
            return [{"query_smiles": value, "matches": [], "status": "unavailable"} for value in smiles_list]
        rows, binding = [], None
        # Model cards and endpoint definitions can be large. Keep each private
        # transfer bounded while preserving the original query population.
        for offset in range(0, len(smiles_list), 32):
            result, source = self._remote_exact_lookup_chunk(smiles_list[offset:offset + 32])
            if binding is not None and source != binding:
                raise ValueError("Exact lookup source changed between chunks")
            binding = source
            rows.extend(result)
        return rows

    def _pin_local_lookup_predictions(self, lookup):
        """Pin validated local lookup metadata without changing the donor DTO.

        Exact query seeds have no catalog source chunk. Their metadata must be
        explicitly pinned here, just like remote exact lookups, before saving
        the predictions. Shared cards may come only from this lookup/model.
        """
        from .scout_predictions import normalize_prediction_values

        rows = deepcopy(lookup)
        predictions, metadata = [], {}
        for result in rows:
            for detail in result["matches"]:
                for prediction in detail.get("admet", []):
                    model_id = prediction.get("model_id")
                    if not isinstance(model_id, str) or not 1 <= len(model_id) <= 256:
                        raise ValueError("Exact lookup prediction model is invalid")
                    predictions.append(prediction)
                    for key in ("model_manifest", "endpoint_definitions"):
                        supplied = prediction.get(key + "_sha256")
                        if not isinstance(supplied, str) or not re.fullmatch(r"[0-9a-f]{64}", supplied):
                            raise ValueError("Exact lookup requires prediction metadata digests")
                        if key in prediction:
                            original = deepcopy(prediction[key])
                            pin_prediction_metadata(prediction, self.artifacts, key)
                            metadata[(model_id, key, supplied)] = original
                            metadata[(model_id, key, prediction[key + "_sha256"])] = original
        for prediction in predictions:
            for key in ("model_manifest", "endpoint_definitions"):
                if key not in prediction:
                    value = metadata.get((prediction["model_id"], key, prediction.get(key + "_sha256")))
                    if value is None:
                        raise ValueError("Exact lookup prediction metadata is unavailable")
                    prediction[key] = deepcopy(value)
                    pin_prediction_metadata(prediction, self.artifacts, key)
            model, endpoints = prediction["model_manifest"], prediction["endpoint_definitions"]
            if (not isinstance(model, dict) or model.get("model_id") != prediction["model_id"]
                    or not isinstance(endpoints, list) or not 1 <= len(endpoints) <= 1000
                    or any(not isinstance(entry, dict) or not isinstance(entry.get("id"), str) for entry in endpoints)
                    or len({entry["id"] for entry in endpoints}) != len(endpoints)):
                raise ValueError("Exact lookup model/endpoint metadata is invalid")
            normalize_prediction_values({**prediction, "values_json": json.dumps(prediction.get("values"), allow_nan=False)}, endpoints)
        return rows

    def _remote_exact_lookup_chunk(self, smiles_list):
        from .identity import digest
        from .scout_predictions import normalize_prediction_values
        request_sha = digest({"smiles": smiles_list})
        response = self._scout_request("/exact-structures", body={"smiles": smiles_list}, response_limit=64 * 1024 * 1024)
        deadline, query_id = time.monotonic() + 120, None
        while True:
            if (not isinstance(response, dict) or response.get("operation") != "exact_structures"
                    or response.get("request_sha256") != request_sha):
                raise ValueError("Exact lookup response is not bound to the requested structures")
            returned_id = response.get("query_id")
            if (not isinstance(returned_id, str) or not re.fullmatch(r"[0-9a-f]{32}", returned_id)
                    or query_id is not None and returned_id != query_id):
                raise ValueError("Exact lookup returned an invalid asynchronous query ID")
            query_id = returned_id
            if response.get("status") != "running":
                break
            if time.monotonic() >= deadline:
                raise ValueError("Exact lookup exceeded its 120 second budget")
            time.sleep(0.5)
            response = self._scout_request("/queries/" + query_id, response_limit=64 * 1024 * 1024)
        source, binding = response.get("source"), response.get("source_binding")
        if (response.get("status") != "complete" or response.get("query_scope") != "full_source_catalog"
                or not isinstance(source, dict) or not isinstance(binding, dict)
                or set(binding) != {"schema", "report_sha256", "compound_count", "database_identity"}
                or binding.get("schema") != "scout-exact-source-v1"
                or source.get("structure_admission") != "reported_verified"
                or source.get("report_sha256") != binding.get("report_sha256")
                or not isinstance(binding.get("report_sha256"), str)
                or not re.fullmatch(r"[0-9a-f]{64}", binding["report_sha256"])
                or type(binding.get("compound_count")) is not int or binding["compound_count"] < 0
                or source.get("compound_count") != binding["compound_count"]):
            raise ValueError("Exact lookup requires a verified catalog population binding")
        database = binding["database_identity"]
        if (not isinstance(database, dict) or set(database) != {"sha256", "size_bytes", "structure_dependencies"}
                or not isinstance(database.get("sha256"), str) or not re.fullmatch(r"[0-9a-f]{64}", database["sha256"])
                or type(database.get("size_bytes")) is not int or database["size_bytes"] < 1
                or not isinstance(database.get("structure_dependencies"), list)):
            raise ValueError("Exact lookup database identity is unavailable")
        for dependency in database["structure_dependencies"]:
            if (not isinstance(dependency, dict) or set(dependency) != {"path", "sha256", "size_bytes"}
                    or not isinstance(dependency.get("path"), str)
                    or not isinstance(dependency.get("sha256"), str)
                    or not re.fullmatch(r"[0-9a-f]{64}", dependency["sha256"])
                    or type(dependency.get("size_bytes")) is not int or dependency["size_bytes"] < 0):
                raise ValueError("Exact lookup dependency identity is unavailable")
        rows = deepcopy(response.get("rows"))
        if not isinstance(rows, list) or len(rows) != len(smiles_list):
            raise ValueError("Exact lookup query population differs from its request")
        predictions, metadata = [], {}
        for value, result in zip(smiles_list, rows):
            identity = chemical_identity(value)
            if (not isinstance(result, dict) or set(result) != {"query_smiles", "identity", "matches", "status"}
                    or result.get("query_smiles") != value or not isinstance(result.get("identity"), dict)
                    or any(result["identity"].get(key) != identity.get(key) for key in
                           ("structure_status", "exact_structure_id", "identity_smiles", "identity_version", "normalizer_runtime"))
                    or not isinstance(result.get("matches"), list) or len(result["matches"]) > 1000):
                raise ValueError("Exact lookup structure identity/order differs from its request")
            expected = "invalid_structure" if identity["structure_status"] != "valid" else "matched" if result["matches"] else "not_found"
            if result.get("status") != expected:
                raise ValueError("Exact lookup status differs from its recorded matches")
            seen = set()
            for detail in result["matches"]:
                compound = detail.get("compound") if isinstance(detail, dict) else None
                cid = compound.get("compound_id") if isinstance(compound, dict) else None
                actual = chemical_identity((compound.get("identity_smiles") or compound.get("canonical_smiles")) if compound else None)
                if (not isinstance(cid, str) or not 1 <= len(cid) <= 256 or cid in seen
                        or detail.get("compound_id") != cid or detail.get("missing") is not False
                        or identity["structure_status"] != "valid" or actual.get("exact_structure_id") != identity["exact_structure_id"]
                        or not isinstance(detail.get("admet"), list) or len(detail["admet"]) > 1000
                        or not isinstance(detail.get("offers"), list) or len(detail["offers"]) > 1000):
                    raise ValueError("Exact lookup contains an unrelated or duplicated compound")
                seen.add(cid)
                if any(not isinstance(offer, dict) or offer.get("compound_id") != cid for offer in detail["offers"]):
                    raise ValueError("Exact lookup supplier rows differ from their compound")
                for prediction in detail["admet"]:
                    if (not isinstance(prediction, dict) or prediction.get("compound_id") != cid
                            or not isinstance(prediction.get("model_id"), str) or not 1 <= len(prediction["model_id"]) <= 256
                            or prediction.get("input_smiles") != compound.get("canonical_smiles")
                            or chemical_identity(prediction.get("input_smiles")).get("exact_structure_id") != identity["exact_structure_id"]):
                        raise ValueError("Exact lookup prediction input differs from its source structure")
                    predictions.append(prediction)
                    for key in ("model_manifest", "endpoint_definitions"):
                        supplied = prediction.get(key + "_sha256")
                        if not isinstance(supplied, str) or not re.fullmatch(r"[0-9a-f]{64}", supplied):
                            raise ValueError("Exact lookup requires prediction metadata digests")
                        if key in prediction:
                            original = deepcopy(prediction[key])
                            pin_prediction_metadata(prediction, self.artifacts, key)
                            metadata[(key, supplied)] = original
                            metadata[(key, prediction[key + "_sha256"])] = original
            result["source_binding"] = deepcopy(binding)
        for prediction in predictions:
            for key in ("model_manifest", "endpoint_definitions"):
                if key not in prediction:
                    value = metadata.get((key, prediction.get(key + "_sha256")))
                    if value is None:
                        raise ValueError("Exact lookup prediction metadata is unavailable")
                    prediction[key] = deepcopy(value)
                    pin_prediction_metadata(prediction, self.artifacts, key)
            model, endpoints = prediction["model_manifest"], prediction["endpoint_definitions"]
            if (not isinstance(model, dict) or model.get("model_id", prediction["model_id"]) != prediction["model_id"]
                    or not isinstance(endpoints, list) or not 1 <= len(endpoints) <= 1000
                    or any(not isinstance(entry, dict) or not isinstance(entry.get("id"), str) for entry in endpoints)
                    or len({entry["id"] for entry in endpoints}) != len(endpoints)):
                raise ValueError("Exact lookup model/endpoint metadata is invalid")
            normalize_prediction_values({**prediction, "values_json": json.dumps(prediction.get("values"), allow_nan=False)}, prediction["endpoint_definitions"])
        return rows, binding

    def _catalog_runtime(self):
        sources = getattr(self, "catalog_sources", None)
        if sources is not None:
            return sources.runtime()
        return getattr(self, "local_catalog", None), getattr(self, "scout", {})

    def catalog_source_context(self, binding):
        return self.catalog_sources.context(binding)

    def close_catalog_sources(self):
        self.catalog_sources.close()

    def _scout_query(self, body, *, owner="researcher"):
        if not isinstance(body, dict):
            raise ValueError("Global catalog query must be an object")
        if set(body) - {
            "q",
            "limit",
            "offset",
            "purchasable_only",
            "hide_severe_admet",
            "admet_policy_sha256",
        }:
            raise ValueError("Unsupported global catalog query fields")
        q, limit, offset = (
            body.get("q", ""),
            body.get("limit", 50),
            body.get("offset", 0),
        )
        purchase = body.get("purchasable_only", False)
        risk = body.get("hide_severe_admet")
        if risk is not None and type(risk) is not bool:
            raise ValueError("Invalid ADMET search switch")
        if not isinstance(q, str) or len(q) > 1024 or not isinstance(purchase, bool):
            raise ValueError("Global catalog search and purchase filter are invalid")
        if (
            isinstance(limit, bool)
            or not isinstance(limit, int)
            or not 1 <= limit <= 1000
        ):
            raise ValueError("Global catalog limit must be 1..1000")
        if (
            isinstance(offset, bool)
            or not isinstance(offset, int)
            or not 0 <= offset <= 100000
        ):
            raise ValueError("Global catalog offset must be 0..100000")
        result = self._scout_request(
            "/query",
            body={
                "q": q,
                "limit": limit,
                "offset": offset,
                "purchasable_only": purchase,
                **({"hide_severe_admet": risk} if risk is not None else {}),
                **(
                    {"admet_policy_sha256": body["admet_policy_sha256"]}
                    if "admet_policy_sha256" in body
                    else {}
                ),
            },
        )
        if result.get("query_id"):
            query_id = result["query_id"]
            self._materialization_path(query_id)  # identical bounded ID contract
            directory = self.root / "catalog-queries"
            directory.mkdir(exist_ok=True)
            with self.lock:
                registry = directory / (query_id + ".json")
                if registry.exists() and json.loads(registry.read_text()).get("owner_id") != owner:
                    raise ValueError("Catalog query identifier already belongs to another owner")
                atomic_bytes(registry, canonical_json({"owner_id": owner, "query_id": query_id}))
        return result

    def _materialization_path(self, materialization_id):
        if (
            not isinstance(materialization_id, str)
            or len(materialization_id) != 32
            or any(char not in "0123456789abcdef" for char in materialization_id)
        ):
            raise ValueError("Invalid catalog materialization ID")
        return self.root / "catalog-materializations" / (materialization_id + ".json")

    def _materialization_state(self, materialization_id, owner):
        path = self._materialization_path(materialization_id)
        if not path.is_file():
            raise KeyError("Unknown catalog materialization")
        state = json.loads(path.read_text())
        if state["owner"] != owner:
            raise KeyError("Unknown catalog materialization")
        return state

    def _materialization_public(self, state):
        public = {
            key: state[key]
            for key in (
                "status",
                "materialization_id",
                "stage",
                "processed",
                "requested",
            )
            if key in state
        }
        if state.get("snapshot_id"):
            public["snapshot"] = self.catalog.manifest(state["snapshot_id"])
        if state.get("error"):
            public["error"] = state["error"]
        return public

    def catalog_snapshot(self, body, *, owner="researcher"):
        """Asynchronously union authoritative Scout records with a frozen pool.

        A registered target pool establishes a destination, not protein or
        receptor validity. Browser-supplied structures, scores and assertions
        never enter this adapter. Restarted pending reads can safely resume.
        """
        if set(body) != {"compound_ids", "target_snapshot_id", "target_id"}:
            raise ValueError(
                "Catalog snapshots require compound_ids, target_snapshot_id and target_id only"
            )
        ids = body["compound_ids"]
        if (
            not isinstance(ids, list)
            or not 1 <= len(ids) <= 10000
            or any(not isinstance(cid, str) or not cid or len(cid) > 256 for cid in ids)
            or len(set(ids)) != len(ids)
        ):
            raise ValueError(
                "Catalog snapshots require 1..10000 distinct bounded compound IDs"
            )
        if not self.scout:
            raise ValueError("Global catalog gateway is not configured")
        manifest = self.catalog.manifest(body["target_snapshot_id"])
        if (
            not isinstance(body["target_id"], str)
            or body["target_id"] not in manifest["target_counts"]
        ):
            raise ValueError("Target does not belong to the registered frozen snapshot")
        request = deepcopy(body)
        request["target_manifest_sha256"] = hashlib.sha256(
            canonical_json(manifest)
        ).hexdigest()
        identifier = hashlib.sha256(
            canonical_json({"owner": owner, "request": request})
        ).hexdigest()[:32]
        path = self._materialization_path(identifier)
        with self.lock:
            previous = json.loads(path.read_text()) if path.exists() else None
            if previous is None or previous.get("status") == "failed":
                if (
                    sum(thread.is_alive() for thread in self._materializations.values())
                    >= 2
                ):
                    raise ValueError(
                        "Two catalog materializations are active; retry after one completes"
                    )
                state = {
                    "materialization_id": identifier,
                    "owner": owner,
                    "request": request,
                    "status": "running",
                    "stage": "authoritative_catalog_fetch",
                    "processed": 0,
                    "requested": len(ids),
                    "chunks": [],
                }
                if previous is not None:
                    state["previous_attempt_sha256"] = self.artifacts.put_bytes(
                        canonical_json(previous)
                    )["sha256"]
                atomic_bytes(path, canonical_json(state))
        return self.catalog_snapshot_status(identifier, owner=owner)

    def catalog_snapshot_status(self, materialization_id, *, owner="researcher"):
        with self.lock:
            state = self._materialization_state(materialization_id, owner)
            if state["status"] == "running" and not self._materializations.get(
                materialization_id, None
            ):
                if (
                    sum(thread.is_alive() for thread in self._materializations.values())
                    >= 2
                ):
                    return {
                        **self._materialization_public(state),
                        "stage": "waiting_for_catalog_slot",
                    }
                thread = threading.Thread(
                    target=self._run_catalog_materialization,
                    args=(materialization_id, owner),
                    daemon=True,
                )
                self._materializations[materialization_id] = thread
                thread.start()
        return self._materialization_public(state)

    def _scout_materialization_chunk(self, ids, *, progress=None, require_offers=True):
        # A full top-K population still has to fit each bounded gateway transfer.
        # Materialize 50 identities at a time without changing retrieval or QC.
        if type(require_offers) is not bool:
            raise ValueError("Supplier evidence requirement must be boolean")
        if (not isinstance(ids, list) or not 1 <= len(ids) <= 1000
                or any(not isinstance(cid, str) or not cid or len(cid) > 256 for cid in ids)
                or len(set(ids)) != len(ids)):
            raise ValueError("Materialization requires 1..1000 unique bounded compound IDs")
        if len(ids) <= 50:
            return self._scout_materialization_batch(ids, progress=progress, require_offers=require_offers)
        from .catalog_cancellation import catalog_checkpoint
        started = time.monotonic()
        combined, rows, chunks = None, [], []
        for offset in range(0, len(ids), 50):
            catalog_checkpoint()
            part = self._scout_materialization_batch(
                ids[offset:offset + 50], progress=progress, require_offers=require_offers)
            if combined is None:
                combined = dict(part)
            elif part["source"] != combined["source"]:
                raise ValueError("Catalog source changed between materialization batches")
            rows.extend(part["rows"])
            chunks.append({"query_id": part.get("query_id"), "candidate_count": len(part["rows"]),
                           "source_sha256": hashlib.sha256(canonical_json(part["source"])).hexdigest()})
        return {**combined, "rows": rows, "materialization_chunks": chunks,
                "elapsed_seconds": round(time.monotonic() - started, 3)}

    def _scout_materialization_batch(self, ids, *, progress=None, require_offers=True):
        if type(require_offers) is not bool:
            raise ValueError("Supplier evidence requirement must be boolean")
        response = self._scout_request(
            "/materialize", body={"compound_ids": ids}, response_limit=64 * 1024 * 1024
        )
        deadline = time.monotonic() + 120
        reported_query_id = None
        poll_retries = 0
        while response.get("status") == "running":
            query_id = response.get("query_id")
            if (
                not isinstance(query_id, str)
                or len(query_id) != 32
                or any(c not in "0123456789abcdef" for c in query_id)
            ):
                raise ValueError(
                    "Catalog gateway returned an invalid asynchronous query ID"
                )
            if progress is not None and query_id != reported_query_id:
                progress({"gateway_query_id": query_id, "gateway_status": "running"})
                reported_query_id = query_id
            if time.monotonic() >= deadline:
                raise ValueError(
                    "Catalog chunk exceeded its 120 second budget; no snapshot was published"
                )
            time.sleep(0.5)
            try:
                response = self._scout_request(
                    "/queries/" + query_id, response_limit=64 * 1024 * 1024
                )
            except _CatalogTransportError:
                # Retry only the existing, read-only query. A failed POST has
                # uncertain admission and must never silently create a second
                # job. Integrity/identity/HTTP failures remain fail-closed.
                if poll_retries >= 2 or time.monotonic() >= deadline:
                    raise
                poll_retries += 1
                if progress is not None:
                    progress({
                        "gateway_query_id": query_id,
                        "gateway_status": "transport_retry",
                        "gateway_poll_retries": poll_retries,
                    })
        if response.get("status") != "complete" or not isinstance(
            response.get("rows"), list
        ):
            raise ValueError(
                "Catalog materialization did not return completed authoritative rows"
            )
        rows = response["rows"]
        returned = [
            row.get("compound", {}).get("compound_id")
            if isinstance(row.get("compound"), dict)
            else row.get("compound_id")
            for row in rows
            if isinstance(row, dict)
        ]
        if returned != ids or any(
            not isinstance(row.get("compound"), dict) or row.get("missing")
            for row in rows
        ):
            raise ValueError(
                "Unknown, missing or reordered catalog IDs; no partial snapshot was published"
            )
        if not isinstance(response.get("source"), dict) or not response["source"].get(
            "report_sha256"
        ):
            raise ValueError(
                "Catalog source audit fingerprint is required for a frozen materialization"
            )
        for detail in rows:
            if require_offers:
                require_complete_offers(detail)
            for prediction in detail.get("admet", []):
                for key in ("model_manifest", "endpoint_definitions"):
                    if key in prediction:
                        pin_prediction_metadata(prediction, self.artifacts, key)
        return response

    def _run_catalog_materialization(self, materialization_id, owner):
        state = self._materialization_state(materialization_id, owner)
        path = self._materialization_path(materialization_id)
        try:
            request = state["request"]
            source = self.catalog.manifest(request["target_snapshot_id"])
            if (
                hashlib.sha256(canonical_json(source)).hexdigest()
                != request["target_manifest_sha256"]
            ):
                raise ValueError(
                    "Target snapshot changed after materialization was requested"
                )
            self.catalog.verify(request["target_snapshot_id"])
            ids = request["compound_ids"]
            chunks = state.get("chunks", [])

            def gateway_progress(event):
                state.update(event)
                atomic_bytes(path, canonical_json(state))

            for start in range(len(chunks) * 1000, len(ids), 1000):
                response = self._scout_materialization_chunk(
                    ids[start : start + 1000], progress=gateway_progress
                )
                artifact = self.artifacts.put_bytes(canonical_json(response))
                if chunks and response["source"] != chunks[0]["source"]:
                    raise ValueError(
                        "Catalog source changed between batches; no mixed snapshot was published"
                    )
                chunks.append(
                    {
                        "sha256": artifact["sha256"],
                        "source": response["source"],
                        "count": len(response["rows"]),
                    }
                )
                state.update(
                    chunks=chunks,
                    processed=sum(chunk["count"] for chunk in chunks),
                    gateway_status="complete",
                )
                atomic_bytes(path, canonical_json(state))
            state["stage"] = "merging_verified_records"
            atomic_bytes(path, canonical_json(state))
            manifest = self._publish_catalog_union(request, chunks, owner=owner)
            state.update(
                status="complete", stage="verified", snapshot_id=manifest["snapshot_id"]
            )
        except Exception as exc:
            # External response bodies and connection details never become UI errors.
            state.update(
                status="failed",
                stage="held",
                error=str(exc)
                if isinstance(exc, ValueError)
                else "Catalog materialization failed: " + type(exc).__name__,
            )
        finally:
            atomic_bytes(path, canonical_json(state))
            with self.lock:
                self._materializations.pop(materialization_id, None)

    def _publish_catalog_union(self, request, chunks, *, owner="researcher"):
        snapshot, target = request["target_snapshot_id"], request["target_id"]
        source = self.catalog.manifest(snapshot)
        if self.catalog.snapshot_owner(snapshot) != owner:
            raise PermissionError("Union parent snapshot belongs to another project")
        details = []
        metadata_artifacts = {}
        for chunk in chunks:
            response = json.loads(self.artifacts.get_bytes(chunk["sha256"]))
            for detail in response["rows"]:
                require_complete_offers(detail)
                detail["source_chunk_sha256"] = chunk["sha256"]
                # The raw response remains intact in CAS. Repeated model and
                # endpoint metadata become exact content references so 10,000
                # compounds do not duplicate the same large model card.
                for prediction in detail.get("admet", []):
                    for key in ("model_manifest", "endpoint_definitions"):
                        if key in prediction:
                            artifact = pin_prediction_metadata(prediction, self.artifacts, key, pop_inline=True)
                            metadata_artifacts[artifact['sha256']] = artifact
            details.extend(response["rows"])
        if len(details) != len(request["compound_ids"]):
            raise ValueError(
                "Materialization receipt count differs from the requested IDs"
            )
        new_groups = {}
        for detail in details:
            compound = detail["compound"]
            smiles = (
                compound.get("identity_smiles")
                or compound.get("canonical_smiles")
                or compound.get("smiles")
            )
            identity = chemical_identity(smiles)
            key = (
                identity["exact_structure_id"] or "unparsed:" + compound["compound_id"]
            )
            new_groups.setdefault(
                key, {"smiles": smiles, "identity": identity, "details": []}
            )["details"].append(detail)
        matched, existing_ids, updates = {}, set(), {}
        for row in self.catalog.iter_candidates(snapshot, target_id=target):
            existing_ids.add(row["candidate_id"])
            identity = row.get("identity") or chemical_identity(row.get("smiles"))
            exact = identity.get("exact_structure_id")
            if exact and exact in new_groups:
                matched.setdefault(exact, []).append(row["candidate_id"])
                updates[row["candidate_id"]] = row

        def annotate(row, group):
            row = deepcopy(row)
            prior = row.get("catalog_sources", [])
            by_source = {
                (item["source"], item["compound"]["compound_id"]): item
                for item in prior
            }
            for detail in group["details"]:
                by_source[("SkinScout", detail["compound"]["compound_id"])] = {
                    "source": "SkinScout",
                    **detail,
                }
            row["catalog_sources"] = list(by_source.values())
            offers = [
                offer
                for item in row["catalog_sources"]
                for offer in item.get("offers", [])
            ]
            row["supplier_offers"] = offers
            row["offers"] = offers
            row["admet"] = [
                prediction
                for item in row["catalog_sources"]
                for prediction in item.get("admet", [])
            ]
            from .admet_search import annotate_risk

            annotate_risk(row)
            if any(
                offer.get("availability") == "catalog_reported_in_stock"
                for offer in offers
            ):
                row["purchasable"] = True
            if isinstance(row.get("states"), dict) and isinstance(
                row["states"].get("experiment"), dict
            ):
                row["states"]["experiment"]["purchasable"] = row.get("purchasable")
            # Skin assay annotations and verified prediction partitions remain
            # source-scoped; they are not automatically target assay evidence.
            row["catalog_evidence_scope"] = (
                "source assay context; target relevance requires review"
            )
            return row

        additions, bindings = [], []
        for key, group in new_groups.items():
            ids = matched.get(key)
            if ids:
                for cid in ids:
                    updates[cid] = annotate(updates[cid], group)
            else:
                cid = (
                    "SCOUT-"
                    + hashlib.sha256(
                        canonical_json(
                            sorted(
                                detail["compound"]["compound_id"]
                                for detail in group["details"]
                            )
                        )
                    ).hexdigest()[:32]
                )
                if cid in existing_ids:
                    raise ValueError(
                        "Scout candidate namespace collides with a different source structure"
                    )
                ids = [cid]
                row = {
                    "candidate_id": cid,
                    "target_id": target,
                    "smiles": group["smiles"],
                    "identity": group["identity"],
                    "parent_id": group["identity"].get("parent_id"),
                    "name": group["details"][0]["compound"].get("name"),
                    "source": "SkinScout",
                    "discovery_score": None,
                    "measurements": [],
                    "docking": [],
                    "target_pool_membership": True,
                    "target_association": False,
                    "target_identity_valid": False,
                    "receptor_available": False,
                    "target_context_status": "registered_pool_without_prepared_receptor",
                    "purchasable": None,
                    "admet_status": "source_annotations_only",
                }
                additions.append(
                    assess_candidate(annotate(row, group), target_identity_valid=False)
                )
            for detail in group["details"]:
                bindings.append(
                    {
                        "compound_id": detail["compound"]["compound_id"],
                        "candidate_ids": ids,
                        "exact_structure_id": group["identity"]["exact_structure_id"],
                    }
                )
        pool_limit = 10000 - len(additions)
        prior_pool = (
            self.catalog.selection_pool(snapshot, target_id=target, limit=pool_limit)
            if pool_limit
            else []
        )
        reranked = {
            row["candidate_id"]: row
            for row in rank_candidates(
                [updates.get(row["candidate_id"], row) for row in prior_pool]
                + additions
            )
        }
        from .planning import target_relations

        marker = source.get("provenance", {}).get("marker")
        relations = [
            relation
            for relation in target_relations()
            if relation["gene"] == target or relation["marker_id"] in {marker, target}
        ]
        provenance = {
            "owner_id": owner,
            "adapter": "scout-union-v1",
            "parent_snapshot_id": snapshot,
            "parent_manifest_sha256": request["target_manifest_sha256"],
            "target_id": target,
            "target_binding_scope": "registered target pool; protein and receptor preparation remain independent",
            "target_relations": relations,
            "metadata_artifact_sha256": sorted(metadata_artifacts),
            "catalog_chunks": chunks,
            "compound_bindings": bindings,
            "added_candidates": len(additions),
            "original_candidate_count": source["candidate_count"],
            "reranking_limit": 10000,
            "references": source.get("provenance", {}).get("references", []),
            "scientific_validation": "NOT_ESTABLISHED",
        }
        identifier = (
            "catalog-" + hashlib.sha256(canonical_json(provenance)).hexdigest()[:24]
        )
        if (self.catalog.root / identifier).exists():
            return self.catalog.manifest(identifier)

        def rows():
            for row in self.catalog.iter_candidates(snapshot):
                if row["target_id"] == target:
                    cid = row["candidate_id"]
                    row = reranked.get(cid, updates.get(cid, row))
                    row["priority_score_status"] = (
                        "recomputed_in_catalog_union_pool"
                        if cid in reranked
                        else "preserved_outside_reranking_budget"
                    )
                yield row
            for row in additions:
                item = reranked[row["candidate_id"]]
                item["priority_score_status"] = "recomputed_in_catalog_union_pool"
                yield item

        manifest = self.catalog.create_snapshot(
            rows(), snapshot_id=identifier, provenance=provenance
        )
        self.catalog.verify(identifier)
        if manifest["candidate_count"] != source["candidate_count"] + len(additions):
            raise ValueError(
                "Catalog union did not preserve the original candidate population"
            )
        return manifest

    def snapshots(self):
        return {
            "snapshots": self.catalog.list_snapshots(),
            "policy_version": PRIORITY_POLICY_VERSION,
        }

    def _legacy_import_path(self, import_id):
        self._materialization_path(import_id)  # Same bounded hexadecimal ID contract.
        return self.root / "legacy-imports" / (import_id + ".json")

    def legacy_import(self, body, *, owner="researcher"):
        if (
            set(body) != {"screen_id"}
            or not isinstance(body["screen_id"], str)
            or not 1 <= len(body["screen_id"]) <= 256
        ):
            raise ValueError("A single registered screen_id is required")
        if self.bridge is None:
            raise ValueError("Legacy source adapter is not configured")
        if body["screen_id"] not in self.legacy_source_grants.get(owner, []):
            # Authorize before source lookup: another tenant cannot use errors
            # to determine whether a guessed marker/result exists.
            raise PermissionError("Legacy source is not assigned to this project")
        digest = sha256_file(self.bridge.screen_file(body["screen_id"]))
        identifier = hashlib.sha256(
            canonical_json(
                {
                    "owner": owner,
                    "screen_id": body["screen_id"],
                    "source_sha256": digest,
                }
            )
        ).hexdigest()[:32]
        path = self._legacy_import_path(identifier)
        with self.lock:
            if not path.exists():
                state = {
                    "import_id": identifier,
                    "owner": owner,
                    "screen_id": body["screen_id"],
                    "source_sha256": digest,
                    "status": "running",
                    "stage": "verifying_source",
                }
                atomic_bytes(path, canonical_json(state))
        return self.legacy_import_status(identifier, owner=owner)

    def legacy_import_status(self, import_id, *, owner="researcher"):
        path = self._legacy_import_path(import_id)
        with self.lock:
            if not path.is_file():
                raise KeyError("Unknown legacy import")
            state = json.loads(path.read_text())
            if state["owner"] != owner:
                raise KeyError("Unknown legacy import")
            if state["status"] == "running" and import_id not in self._legacy_imports:
                if (
                    sum(thread.is_alive() for thread in self._legacy_imports.values())
                    < 2
                ):
                    thread = threading.Thread(
                        target=self._run_legacy_import,
                        args=(import_id, owner),
                        daemon=True,
                    )
                    self._legacy_imports[import_id] = thread
                    thread.start()
                else:
                    state["stage"] = "waiting_for_import_slot"
        result = {
            key: state[key]
            for key in ("import_id", "screen_id", "status", "stage", "error")
            if key in state
        }
        if state.get("snapshot_id"):
            result["snapshot"] = self.catalog.manifest(state["snapshot_id"])
        return result

    def _run_legacy_import(self, import_id, owner):
        path = self._legacy_import_path(import_id)
        state = json.loads(path.read_text())
        try:
            if state["owner"] != owner:
                raise ValueError("Legacy import owner changed")
            manifest = self.import_legacy(
                state["screen_id"],
                activate=False,
                expected_source_sha256=state["source_sha256"],
                owner=owner,
            )
            self.catalog.verify(manifest["snapshot_id"])
            if (
                sha256_file(self.bridge.screen_file(state["screen_id"]))
                != state["source_sha256"]
            ):
                raise ValueError(
                    "Legacy source changed before activation; existing pointer is preserved"
                )
            if owner == "researcher":
                self._activate_snapshot(state["screen_id"], manifest["snapshot_id"])
            state.update(
                status="complete", stage="verified", snapshot_id=manifest["snapshot_id"]
            )
        except Exception as exc:
            state.update(
                status="failed",
                stage="held",
                error=str(exc)
                if isinstance(exc, ValueError)
                else "Legacy import failed: " + type(exc).__name__,
            )
        finally:
            atomic_bytes(path, canonical_json(state))
            with self.lock:
                self._legacy_imports.pop(import_id, None)

    def _population(self, snapshot_id, target_id):
        # Selection is bounded at the approved cheap reranking stage. Full
        # catalog exploration and filters remain SQL operations over every row.
        query = self.catalog.query(
            snapshot_id,
            target_id=target_id,
            sort="priority_score_v2",
            limit=1000,
            page_size=50,
        )
        total = query["counts"]["population"]
        rows = self.catalog.selection_pool(
            snapshot_id, target_id=target_id, limit=10000
        )
        existing = {r["candidate_id"] for r in rows}
        for cid in (
            self.catalog.manifest(snapshot_id)
            .get("provenance", {})
            .get("reference_control_ids", [])
        ):
            if cid not in existing:
                rows.append(self.catalog.get_candidate(snapshot_id, cid, target_id))
        return rows, total

    def _batch(self, body, *, owner="researcher"):
        snapshot, target = body["snapshot_id"], body["target_id"]
        manifest = self.catalog.manifest(snapshot)
        rows, total = self._population(snapshot, target)
        references = manifest.get("provenance", {}).get("references", [])
        requested = (
            validate_requested_action(body["desired_action"])
            if "desired_action" in body
            else None
        )
        rows = bind_target_goals(rows, desired_action=requested)
        generalized = bool(manifest.get("provenance", {}).get("discovery_run_id"))
        selection = select_experiment_batch(
            rows,
            target,
            references=references,
            trust_stored_identity=True,
            require_documented_procurement=generalized,
            preserve_frozen_ranking=generalized,
        )
        for row in selection["candidates"] + selection["controls"]:
            row["selection_lane"] = row.get(
                "selection_group", row.get("role", "control")
            )
        policy = {
            "version": PRIORITY_POLICY_VERSION,
            "weights": {"measured": 0.4, "structural": 0.3, "discovery": 0.3},
            "owner_id": owner,
            "procurement_requirement": "dated-exact-form-v1"
            if generalized
            else "catalog-record-v1",
            "purchase_preference": True,
            "target_population": total,
            "selection_pool": len(rows),
            "desired_action": requested,
            "goal_interpretation": "user hypothesis; not verified skin biology",
            "reference_set": references,
            "selection_report": {
                k: v
                for k, v in selection.items()
                if k not in {"candidates", "controls"}
            },
            "interpretation": "experimental priority heuristic; not efficacy or binding probability",
        }
        if generalized and not selection["candidates"]:
            return {
                "status": "completed_zero",
                "reason": "procurement_unconfirmed",
                "selection": selection,
                "snapshot_id": snapshot,
                "target_id": target,
                "candidate_count": 0,
                "control_count": len(selection["controls"]),
            }
        return self.experiments.create_batch(
            target, selection["candidates"], selection["controls"], snapshot, policy
        )

    def _feedback_records(self, target_id, *, batch_ids=None, owner="researcher"):
        records, receipts, selected = {}, [], set()
        if batch_ids is not None and (
            not isinstance(batch_ids, list)
            or any(not isinstance(x, str) for x in batch_ids)
            or len(set(batch_ids)) != len(batch_ids)
        ):
            raise ValueError("batch_ids must be distinct strings")
        requested = set(batch_ids) if batch_ids is not None else None
        for batch in self.experiments.list_batches(target_id):
            if batch.get("policy", {}).get("owner_id") != owner:
                continue
            if requested is not None and batch.get("batch_id") not in requested:
                continue
            if batch.get("status") != "frozen":
                if requested is not None:
                    raise ValueError(
                        "Feedback requires frozen intact experimental batches"
                    )
                continue
            selected.add(batch["batch_id"])
            bindings = {
                r["candidate_id"]: chemical_identity(r["smiles"])["exact_structure_id"]
                for r in batch["candidates"] + batch["controls"]
            }
            observations = self.experiments.observations(batch["batch_id"])
            receipts.append(
                {
                    "batch_id": batch["batch_id"],
                    "batch_sha256": batch["payload_sha256"],
                    "observations_sha256": hashlib.sha256(
                        canonical_json(observations)
                    ).hexdigest(),
                    "observation_count": len(observations),
                }
            )
            for observation in observations:
                cid = observation["normalized"]["candidate_id"]
                if cid in bindings:
                    item = deepcopy(observation)
                    item["experiment_batch_id"] = batch["batch_id"]
                    item["experiment_batch_sha256"] = batch["payload_sha256"]
                    records.setdefault((cid, bindings[cid]), []).append(item)
        if requested is not None and selected != requested:
            raise ValueError("Unknown or different-target feedback batch")
        return records, sorted(receipts, key=lambda item: item["batch_id"])

    def _candidate_feedback(self, candidate, *, owner="researcher"):
        records, _ = self._feedback_records(candidate["target_id"], owner=owner)
        identity = chemical_identity(candidate.get("smiles"))
        observations = records.get(
            (candidate["candidate_id"], identity["exact_structure_id"]), []
        )
        return measured_feedback(observations)

    def feedback_snapshot(self, body, *, owner="researcher"):
        """Publish a new immutable snapshot; never rewrite its parent or batch.

        Numeric ranking requires explicit reviewed assay semantics. Observations
        lacking those semantics remain visible measured evidence. Only the cheap
        reranking pool is rescored; every other source candidate is preserved.
        """
        snapshot, target = body["snapshot_id"], body["target_id"]
        manifest = self.catalog.manifest(snapshot)
        if target not in manifest["target_counts"]:
            raise ValueError("Target does not exist in this snapshot")
        requested_action = (
            validate_requested_action(body["desired_action"])
            if "desired_action" in body
            else None
        )
        contexts = body.get("assay_contexts", {})
        if not isinstance(contexts, dict):
            raise ValueError("assay_contexts must be an object keyed by assay_id")
        for assay, context in contexts.items():
            if not isinstance(assay, str) or not assay or not isinstance(context, dict):
                raise ValueError(
                    "Each reviewed assay context requires an assay_id and object"
                )
            if any(
                not isinstance(context.get(key), str) or not context[key].strip()
                for key in ("species", "assay_type", "action")
            ) or not isinstance(context.get("higher_is_better"), bool):
                raise ValueError(
                    "Reviewed contexts require species, assay_type, action and higher_is_better"
                )
            if context.get("reviewed") is not True:
                raise ValueError(
                    "Assay interpretation must be explicitly reviewed before ranking"
                )
        records, receipts = self._feedback_records(
            target, batch_ids=body.get("batch_ids"), owner=owner
        )
        if not receipts or not any(records.values()):
            raise ValueError("No registered experimental observations are available")
        observed_ids = {cid for cid, _ in records}
        provenance = {
            "parent_snapshot_id": snapshot,
            "owner_id": owner,
            "parent_manifest_sha256": hashlib.sha256(
                canonical_json(manifest)
            ).hexdigest(),
            "target_id": target,
            "feedback_receipts": receipts,
            "assay_contexts": contexts,
            "target_goal_request": {
                "desired_action": requested_action,
                "origin": "operator_request"
                if "desired_action" in body
                else "preserved_source",
            },
            "policy_version": PRIORITY_POLICY_VERSION,
            "reranking_limit": 10000,
            "references": manifest.get("provenance", {}).get("references", []),
            "interpretation": "Measured assay feedback; missing context does not establish comparable ranks",
        }
        snapshot_id = (
            "feedback-" + hashlib.sha256(canonical_json(provenance)).hexdigest()[:24]
        )
        if (self.catalog.root / snapshot_id).exists():
            return self.catalog.manifest(snapshot_id)

        def enrich(row):
            row = deepcopy(row)
            if "desired_action" in body:
                row = bind_target_goals([row], desired_action=requested_action)[0]
            if row["candidate_id"] not in observed_ids:
                return row
            identity = chemical_identity(row.get("smiles"))
            observations = records.get(
                (row["candidate_id"], identity["exact_structure_id"]), []
            )
            if not observations:
                return row
            prior = row.get("experimental_observations", [])
            by_key = {
                (r.get("experiment_batch_id"), r.get("observation_key")): r
                for r in prior
            }
            by_key.update(
                {
                    (r["experiment_batch_id"], r["observation_key"]): r
                    for r in observations
                }
            )
            row["experimental_observations"] = list(by_key.values())
            row["experimental_feedback"] = measured_feedback(
                row["experimental_observations"]
            )
            original = row.get("measurements") or []
            measurements = [original] if isinstance(original, dict) else list(original)
            # Replace only this adapter's previous derived records; preserve all
            # independently sourced legacy measurements and their null values.
            measurements = [
                m
                for m in measurements
                if not isinstance(m, dict)
                or m.get("adapter") != "experimental-feedback-v1"
            ]
            for observation in row["experimental_observations"]:
                value = observation["normalized"]
                if value.get("status") != "measured":
                    continue
                measured = {
                    **value,
                    "is_measured": True,
                    "adapter": "experimental-feedback-v1",
                    "observation_key": observation["observation_key"],
                    "experiment_batch_id": observation["experiment_batch_id"],
                }
                context = contexts.get(value.get("assay_id"))
                if context:
                    measured.update(
                        {
                            key: context[key]
                            for key in (
                                "species",
                                "assay_type",
                                "action",
                                "higher_is_better",
                            )
                        }
                    )
                    measured["assay_context"] = canonical_json(
                        {
                            key: value.get(key)
                            for key in (
                                "assay_id",
                                "conditions",
                                "biological_context",
                                "concentration",
                                "concentration_unit",
                                "source_id",
                            )
                        }
                    ).decode()
                    measured["context_provenance"] = (
                        "operator-reviewed frozen assay_contexts"
                    )
                measurements.append(measured)
            row["measurements"] = measurements
            return row

        pool = self.catalog.selection_pool(snapshot, target_id=target, limit=10000)
        by_id = {row["candidate_id"]: row for row in pool}
        # Ensure uploaded candidates participate before filling the remaining
        # bounded pool from the original ranking.
        observed = {}
        for cid, exact_id in records:
            row = self.catalog.get_candidate(snapshot, cid, target_id=target)
            if (
                row
                and chemical_identity(row.get("smiles"))["exact_structure_id"]
                == exact_id
            ):
                observed[cid] = row
        if not observed:
            raise ValueError(
                "No observations match the exact candidate structures in this snapshot"
            )
        if len(observed) > 10000:
            raise ValueError(
                "Observed cohort exceeds the reranking budget; choose explicit batch_ids"
            )
        combined = list(observed.values()) + [
            row for cid, row in by_id.items() if cid not in observed
        ]
        reranked = rank_candidates([enrich(row) for row in combined[:10000]])
        scores = {row["candidate_id"]: row for row in reranked}
        provenance["reranking_candidate_count"] = len(reranked)
        provenance["reranking_candidate_ids_sha256"] = hashlib.sha256(
            canonical_json(sorted(scores))
        ).hexdigest()

        def rows():
            for row in self.catalog.iter_candidates(snapshot):
                if row["target_id"] == target:
                    if row["candidate_id"] in scores:
                        row = scores[row["candidate_id"]]
                        row["priority_score_status"] = (
                            "recomputed_in_frozen_feedback_cohort"
                        )
                    else:
                        row = enrich(row)
                        row["priority_score_status"] = (
                            "preserved_outside_reranking_budget"
                        )
                yield row

        try:
            return self.catalog.create_snapshot(
                rows(), snapshot_id=snapshot_id, provenance=provenance
            )
        except FileExistsError:
            return self.catalog.manifest(snapshot_id)

    @staticmethod
    def _decode_upload(body, *, limit=12 * 1024 * 1024):
        value = body.get("base64")
        if not isinstance(value, str) or len(value) > (limit * 4 // 3 + 4):
            raise ValueError("Upload exceeds size limit")
        raw = base64.b64decode(value, validate=True)
        if len(raw) > limit:
            raise ValueError("Upload exceeds size limit")
        return raw

    def get(self, url, *, owner="researcher"):
        self._authorize_data(url, owner=owner)
        result = self._get(url, owner=owner)
        if urlparse(url).path == "/v2/snapshots":
            result = {
                "snapshots": [
                    s
                    for s in result["snapshots"]
                    if self.catalog.snapshot_owner(s["snapshot_id"]) == owner
                ]
            }
        if urlparse(url).path == "/v2/experiment-batches":
            result["batches"] = [
                b
                for b in result["batches"]
                if b.get("policy", {}).get("owner_id") == owner
            ]
        return result

    def _authorize_data(self, url, *, owner, body=None):
        parts = [unquote(p) for p in urlparse(url).path.strip("/").split("/")]
        query = {k: v[-1] for k, v in parse_qs(urlparse(url).query).items()}
        body = body or {}
        for key in ("snapshot_id", "target_snapshot_id"):
            snapshot = body.get(key) or query.get(key)
            if parts[:2] == ["v2", "snapshots"] and len(parts) == 3:
                snapshot = parts[2]
            if (
                snapshot
                and self.catalog.snapshot_owner(snapshot) != owner
            ):
                raise PermissionError("Snapshot belongs to another project")
        batch = body.get("batch_id")
        if parts[:2] == ["v2", "experiment-batches"] and len(parts) >= 3:
            batch = parts[2]
        if parts[:2] == ["v2", "result-imports"] and len(parts) >= 3:
            batch = self.experiments.get_import(parts[2])["batch_id"]
        if (
            batch
            and self.experiments.get_batch(batch)
            .get("policy", {})
            .get("owner_id")
            != owner
        ):
            raise PermissionError("Experimental batch belongs to another project")

    def _get(self, url, *, owner="researcher"):
        parsed = urlparse(url)
        parts = [unquote(p) for p in parsed.path.strip("/").split("/")]
        query = {k: v[-1] for k, v in parse_qs(parsed.query).items()}
        if parts == ["v2", "research-capabilities"]:
            from .research_toxicity import research_toxicity_capability
            return {"additional_toxicity": research_toxicity_capability(self.discovery.config)}
        if len(parts) == 4 and parts[:2] == ["v2", "research-inputs"] and parts[3] == "structure":
            return self.research_inputs.structure(parts[2], owner=owner)
        if len(parts) == 3 and parts[:2] == ["v2", "research-inputs"]:
            return {**self.research_inputs.get(parts[2], owner=owner), "resolution_id": parts[2]}
        if len(parts) == 4 and parts[:2] == ["v2", "research-inputs"] and parts[3] == "image":
            return self.research_inputs.image(parts[2], owner=owner)
        if len(parts) == 3 and parts[:2] in (["v2", "exports"], ["v2", "result-exports"]):
            return self.exports.status(parts[2], owner=owner)
        personal = getattr(self, "personal_literature", None)
        if (parts[:2] == ["v2", "literature"] and personal is not None
                and (personal.configured or parts != ["v2", "literature", "usage"])):
            return personal.get(url, owner=owner)
        if parts == ["v2", "literature", "usage"]:
            # Site identity grants no plan credential. Never borrow a host owner's
            # loopback account, Codex account, or other visitor's private ledger.
            return {
                "schema": "skinpipe-chatgpt-usage-v1",
                "scope": "shared_host_unconfigured",
                "account": {"connected": False, "plan_permission": False},
                "provider_plan": {"status": "UNKNOWN", "remaining_percent": None,
                                  "remaining_tokens": None, "reset_at": None,
                                  "manage_url": "https://chatgpt.com/settings/usage"},
                "app_usage": {"accounting_status": "unavailable"},
                "last_limit_error": None,
            }
        if parts == ["v2", "discovery-runs"]:
            if query.get("view") == "status":
                return self.discovery.list_status(owner=owner)
            return self.discovery.list(owner=owner)
        if len(parts) == 4 and parts[:2] == ["v2", "discovery-runs"] and parts[3] == "summary":
            return self.research_assistance.summary(parts[2], owner=owner)
        if len(parts) in {3, 4} and parts[:2] == ["v2", "discovery-runs"]:
            if len(parts) == 3:
                if query.get("view") == "status":
                    return self.discovery.status(parts[2], owner=owner)
                return self.discovery.get(parts[2], owner=owner)
            if parts[3] == "events":
                return self.discovery.events(
                    parts[2], owner=owner, after=int(query.get("after", "0"))
                )
        if parts == ["v2", "algorithms"]:
            return {
                "compute_profile": compute_profile(),
                "providers": list(PROVIDERS),
                "resources": self.queue.resources(),
                "ranking_weight": 0,
                "scientific_validation": "NOT_ESTABLISHED",
            }
        if parts[:2] in (["v2", "campaigns"], ["v2", "campaign-inputs"]):
            return self.campaigns.get(url, owner=owner)
        if parts == ["v2", "health"]:
            return {
                "service": "skinpipe",
                "api_version": 2,
                "operational_store": "postgresql",
                "catalog_store": "duckdb-parquet",
                "scientific_validation": "not_established",
            }
        if parts == ["v2", "optional-scorers"]:
            return {
                "providers": self.optional_scorers.status(),
                "execution_via_research_api": False,
                "automatic_candidates": 0,
                "maximum_opt_in_candidates": 12,
                "ranking_weight": 0,
            }
        if parts == ["v2", "snapshots"]:
            # Only an explicit summary view changes the legacy full response.
            views = parse_qs(parsed.query, keep_blank_values=True).get("view", ["full"])
            if len(views) != 1 or views[0] not in {"full", "summary"}:
                raise ValueError("Snapshot list view must be full or summary exactly once")
            result = self.snapshots()
            if views[0] == "summary":
                result["snapshots"] = [_snapshot_list_summary(s) for s in result["snapshots"]]
            return result
        if len(parts) == 3 and parts[:2] == ["v2", "legacy-imports"]:
            return self.legacy_import_status(parts[2], owner=owner)
        if parts == ["v2", "catalog-status"]:
            if not self.scout:
                return {
                    "configured": False,
                    "status": "unavailable",
                    "reason": "not_configured",
                    "scientific_validation": "NOT_ESTABLISHED",
                }
            return self._scout_request("/status")
        if len(parts) == 3 and parts[:2] == ["v2", "catalog-snapshots"]:
            return self.catalog_snapshot_status(parts[2], owner=owner)
        if len(parts) == 3 and parts[:2] == ["v2", "catalog-queries"]:
            query_id = parts[2]
            if len(query_id) != 32 or any(
                char not in "0123456789abcdef" for char in query_id
            ):
                raise ValueError("Invalid global catalog query ID")
            path = self.root / "catalog-queries" / (query_id + ".json")
            if not path.is_file() or json.loads(path.read_bytes()).get("owner_id") != owner:
                raise KeyError("Unknown global catalog query")
            return self._scout_request("/queries/" + query_id)
        if len(parts) == 3 and parts[:2] == ["v2", "catalog-candidates"]:
            compound_id = parts[2]
            if (
                not compound_id
                or len(compound_id) > 256
                or any(not (char.isalnum() or char in "_.:-") for char in compound_id)
            ):
                raise ValueError("Invalid global catalog compound ID")
            return self._scout_request("/compound/" + quote(compound_id, safe=""))
        if len(parts) == 3 and parts[:2] == ["v2", "snapshots"]:
            return self.catalog.manifest(parts[2])
        if len(parts) == 3 and parts[:2] == ["v2", "candidates"]:
            candidate = self.catalog.get_candidate(
                query["snapshot_id"], parts[2], target_id=query.get("target_id")
            )
            if candidate is None:
                raise KeyError("Unknown candidate")
            candidate["experimental_feedback"] = self._candidate_feedback(
                candidate, owner=owner
            )
            from .cached_prediction_views import CachedPredictionView
            candidate["admet"] = CachedPredictionView(self.artifacts).for_row(candidate)
            return candidate
        if (
            len(parts) == 4
            and parts[:2] == ["v2", "candidates"]
            and parts[3] == "alignment"
        ):
            from .catalog_hopping import stored_alignment

            candidate = self.catalog.get_candidate(
                query["snapshot_id"], parts[2], target_id=query["target_id"]
            )
            if candidate is None:
                raise KeyError("Unknown candidate")
            return stored_alignment(candidate, self.artifacts)
        if (
            len(parts) == 4
            and parts[:2] == ["v2", "candidates"]
            and parts[3] == "docking"
        ):
            from .pose_view import candidate_docking_view

            try:
                index = int(query.get("observation_index", "0"))
            except (TypeError, ValueError):
                raise ValueError(
                    "Docking observation_index must be a nonnegative integer"
                ) from None
            if not 0 <= index <= 100:
                raise ValueError("Docking observation_index must be between 0 and 100")
            return candidate_docking_view(
                self.catalog,
                self.artifacts,
                parts[2],
                snapshot_id=query["snapshot_id"],
                target_id=query["target_id"],
                observation_index=index,
            )
        if parts == ["v2", "resources"]:
            return self.queue.resources()
        if parts == ["v2", "jobs"]:
            return {"jobs": self._jobs_budget(self.queue.list(owner=owner))}
        if len(parts) in {3, 4} and parts[:2] == ["v2", "jobs"]:
            if len(parts) == 3:
                return self._job_budget(self.queue.get(parts[2], owner=owner))
            if parts[3] == "events":
                return {
                    "events": self.queue.events(
                        parts[2], owner=owner, after=query.get("after", 0)
                    )
                }
        if parts == ["v2", "experiment-batches"]:
            return {"batches": self.experiments.list_batches(query.get("target_id"))}
        if len(parts) >= 3 and parts[:2] == ["v2", "experiment-batches"]:
            if len(parts) == 3:
                return self.experiments.get_batch(parts[2])
            if parts[3:] == ["observations"]:
                observations = self.experiments.observations(parts[2])
                return {
                    "observations": observations,
                    "feedback": measured_feedback(observations),
                }
            if parts[3:] == ["export"]:
                files = self.experiments.export_batch(parts[2])
                return {
                    "files": [
                        {
                            "name": name,
                            "sha256": hashlib.sha256(raw).hexdigest(),
                            "base64": base64.b64encode(raw).decode(),
                        }
                        for name, raw in files.items()
                    ]
                }
        if len(parts) == 3 and parts[:2] == ["v2", "result-imports"]:
            return self.experiments.get_import(parts[2])
        raise KeyError("Unknown v2 route")

    def post(self, url, body, *, owner="researcher"):
        if not isinstance(body, dict):
            raise ValueError("JSON object required")
        if urlparse(url).path == "/v2/result-imports":
            value = body.get("base64")
            if not isinstance(value, str) or len(value) > MAX_UPLOAD_BYTES * 4 // 3 + 4:
                raise ValueError("Upload exceeds size limit")
        self._authorize_data(url, owner=owner, body=body)
        return self._post(url, body, owner=owner)

    def _post(self, url, body, *, owner="researcher"):
        if not isinstance(body, dict):
            raise ValueError("JSON object required")
        parts = [unquote(p) for p in urlparse(url).path.strip("/").split("/")]
        if parts == ["v2", "literature", "reviews"]:
            return self.personal_literature.submit(body, owner=owner)
        if (
            len(parts) == 4
            and parts[:2] == ["v2", "campaigns"]
            and parts[3] == "algorithm-jobs"
        ):
            return self._research_campaign_action(parts[2], owner=owner,
                action=lambda: self.algorithm_jobs.submit(parts[2], body, owner=owner))
        if parts == ["v2", "algorithm-evaluations"]:
            if set(body) != {"evaluation_manifest_sha256"}:
                raise ValueError("Evaluation accepts an immutable manifest SHA256 only")
            manifest = self.campaigns._json_artifact(body["evaluation_manifest_sha256"])
            campaign_id = manifest.get("campaign_id")
            profile = compute_profile()
            if campaign_id is not None:
                campaign = self.campaigns._coordinator(owner)._manifest(campaign_id)
                profile = campaign["plan"].get("compute_profile") or profile
                if manifest.get("compute_profile_sha256") != profile["sha256"]:
                    raise ValueError("Evaluation must bind the registered campaign policy")
            # Legacy independent evaluation remains supported, but all such
            # requests share one owner/profile allowance instead of refreshing
            # one hour for each prediction/model/seed manifest.
            budget_scope = "optional:" + hashlib.sha256(canonical_json({
                "owner": owner, "campaign": campaign_id or "independent-evaluation",
                "profile": profile["sha256"],
            })).hexdigest()
            available = [
                node
                for node in self.queue.resources()["nodes"]
                if node.get("fresh")
                and "algorithm_evaluation"
                in node.get("telemetry", {}).get("capabilities", [])
            ]
            if not available:
                return {
                    "status": "held",
                    "reason": "registered_evaluation_worker_unavailable",
                    "promoted": False,
                }
            identity = hashlib.sha256(
                canonical_json({"owner": owner, "evaluation": body, "budget_scope": budget_scope,
                                "admission_version": "shared-evaluation-v3"})
            ).hexdigest()
            def enqueue_evaluation():
                return self.queue.enqueue(
                    "algorithm_evaluation",
                    {**body, **({"campaign_id": campaign_id} if campaign_id else {})},
                    owner=owner,
                    project="evaluation-" + hashlib.sha256(owner.encode()).hexdigest()[:32],
                    idempotency_key="evaluation-" + identity,
                    budget_seconds=3600,
                    budget_scope=budget_scope,
                )
            job = self._research_campaign_action(campaign_id, owner=owner,
                action=enqueue_evaluation) if campaign_id else enqueue_evaluation()
            return {"job_id": job["id"], "status": job["status"], "promoted": False}
        if parts[:2] in (["v2", "campaigns"], ["v2", "campaign-inputs"]):
            if len(parts) == 4 and parts[:2] == ["v2", "campaigns"] and parts[3] in {
                    "run", "dispatch", "advance", "confirmation-policy"}:
                return self._research_campaign_action(parts[2], owner=owner,
                    action=lambda: self.campaigns.post(url, body, owner=owner))
            return self.campaigns.post(url, body, owner=owner)
        if parts == ["v2", "targets", "resolve"]:
            return self.discovery.resolve(body)
        if parts == ["v2", "research-inputs", "resolve"]:
            if not isinstance(body, dict) or set(body) != {"input"}:
                raise ValueError("Research input resolution requires an input object")
            return self.research_inputs.resolve(body["input"], owner=owner)
        if parts == ["v2", "research-inputs", "uploads"]:
            return self.research_structures.upload(body, owner=owner)
        if parts == ["v2", "research-assistance"]:
            return self.research_assistance.assist(body, owner=owner)
        if parts in (["v2", "exports"], ["v2", "result-exports"]):
            return self.exports.submit(body, owner=owner)
        if parts == ["v2", "discovery-runs"]:
            return self.discovery.submit(body, owner=owner)
        if len(parts) == 4 and parts[:2] == ["v2", "discovery-runs"]:
            if parts[3] in {"cancel", "resume", "continue"}:
                if not isinstance(body, dict) or body:
                    raise ValueError("Research controls do not change frozen inputs or budgets")
                return self._control_discovery_run(parts[2], parts[3], owner=owner)
        if parts == ["v2", "candidate-queries"]:
            allowed = {
                "snapshot_id",
                "target_id",
                "q",
                "filters",
                "sort",
                "descending",
                "limit",
                "page_size",
                "offset",
                "purchase_first",
            }
            if set(body) - allowed:
                raise ValueError("Unsupported candidate query fields")
            return self.catalog.query(**body)
        if parts == ["v2", "catalog-queries"]:
            return self._scout_query(body, owner=owner)
        if parts == ["v2", "catalog-snapshots"]:
            return self.catalog_snapshot(body, owner=owner)
        if parts == ["v2", "legacy-imports"]:
            return self.legacy_import(body, owner=owner)
        if parts == ["v2", "correspondence"]:
            reference, candidate = (
                body.get("reference_smiles"),
                body.get("candidate_smiles"),
            )
            if not all(
                isinstance(x, str) and 0 < len(x) <= 4096
                for x in (reference, candidate)
            ):
                raise ValueError("Two bounded structure strings required")
            return compare_structures(reference, candidate)
        if parts == ["v2", "experiment-batches"]:
            return self._batch(body, owner=owner)
        if parts == ["v2", "feedback-snapshots"]:
            return self.feedback_snapshot(body, owner=owner)
        if len(parts) == 4 and parts[:2] == ["v2", "experiment-batches"]:
            if parts[3] == "freeze":
                return self.experiments.freeze_batch(parts[2])
            if parts[3] == "update":
                # Only reordering/removing catalog-bound rows is accepted from
                # the browser; supplied structures/scores never replace source.
                batch = self.experiments.get_batch(parts[2])
                candidate_rows = {r["candidate_id"]: r for r in batch["candidates"]}
                control_rows = {r["candidate_id"]: r for r in batch["controls"]}
                for key, allowed in (
                    ("candidate_ids", candidate_rows),
                    ("control_ids", control_rows),
                ):
                    ids = body.get(key)
                    if (
                        not isinstance(ids, list)
                        or any(not isinstance(cid, str) for cid in ids)
                        or len(set(ids)) != len(ids)
                        or not set(ids).issubset(allowed)
                    ):
                        raise ValueError(
                            "Draft edits may only reorder/remove rows within their original candidate or control role"
                        )
                candidates = [candidate_rows[cid] for cid in body["candidate_ids"]]
                controls = [control_rows[cid] for cid in body["control_ids"]]
                return self.experiments.update_draft(parts[2], candidates, controls)
        if parts == ["v2", "result-imports"]:
            return self.experiments.preview_import(
                body["batch_id"],
                self._decode_upload(body, limit=MAX_UPLOAD_BYTES),
                body["filename"],
            )
        if (
            len(parts) == 4
            and parts[:2] == ["v2", "result-imports"]
            and parts[3] == "confirm"
        ):
            return self.experiments.confirm_import(parts[2], body["idempotency_key"])
        if parts == ["v2", "jobs"]:
            if set(body) - {"kind", "request", "idempotency_key", "budget_seconds"}:
                raise ValueError("Unsupported job submission fields")
            kind = body.get("kind")
            if kind not in {"probe", "dock_batch"}:
                raise ValueError("Unsupported worker task kind")
            request = body.get("request", {})
            if not isinstance(request, dict) or "budget_scope" in request:
                raise ValueError("Computation budgets are assigned by the server")
            budget = body.get("budget_seconds", 28800)
            from .ranking import finite
            if finite(budget) is None or not 0 < budget <= 28800:
                raise ValueError("Worker computation budget must be within the server 8h cap")
            if kind == "dock_batch":
                if set(request) - {"campaign_id", "stage", "chunk_size", "retry_failed"} or not request.get("campaign_id"):
                    raise ValueError("Docking jobs require a registered campaign; raw preparations cannot bypass admission")
                return self._research_campaign_action(request["campaign_id"], owner=owner,
                    action=lambda: self.campaigns.post(
                        "/v2/campaigns/" + request["campaign_id"] + "/dispatch",
                        {k: v for k, v in request.items() if k != "campaign_id"}, owner=owner,
                    ))
            return self.queue.enqueue(
                kind,
                request,
                owner=owner,
                project="probes-" + hashlib.sha256(owner.encode()).hexdigest()[:24],
                idempotency_key=body["idempotency_key"],
                budget_seconds=budget,
                budget_scope="probe:" + hashlib.sha256(owner.encode()).hexdigest(),
            )
        if len(parts) == 4 and parts[:2] == ["v2", "jobs"]:
            if parts[3] == "cancel":
                return self.queue.cancel(parts[2], owner=owner)
            if parts[3] == "resume":
                job = self._job_budget(self.queue.get(parts[2], owner=owner))
                self._require_whole_study_resume(job, owner=owner)
                if job.get("kind") == "literature_review":
                    raise ValueError("Literature review attempts cannot be resumed or replayed")
                if job['budget']['status'] != 'verified':
                    raise ValueError('Resume requires a verified frozen server budget')
                requested = body.get("budget_seconds")
                if requested is not None:
                    from .ranking import finite
                    if finite(requested) is None or requested > job["budget"]["maximum_seconds"]:
                        raise ValueError("Resume cannot exceed the frozen server budget; create a new run")
                return self.queue.resume(
                    parts[2], owner=owner, budget_seconds=body.get("budget_seconds")
                )
        if parts == ["v2", "plans"]:
            return self.campaigns.plan(body, owner=owner)
        raise KeyError("Unknown v2 route")

    def _research_campaign_action(self, campaign_id, *, owner, action):
        with self.queue.connection() as conn:
            linked = conn.execute("SELECT id FROM platform_discovery_runs WHERE owner_id=%s AND campaign_id=%s LIMIT 1",
                                  (owner, campaign_id)).fetchone()
        if not linked:
            return action()
        run_id = linked["id"]
        with self.queue.connection() as conn:
            conn.execute("SELECT pg_advisory_xact_lock(hashtextextended(%s,0))",
                         ("research-control:" + owner + ":" + run_id,))
            if self.discovery._record(run_id, owner)["desired_state"] != "active":
                raise ValueError("Cancelled research cannot dispatch downstream work; use whole-study resume")
            return action()

    def _require_whole_study_resume(self, job, *, owner):
        # Advanced job controls must not bypass a study's cancellation fence.
        # Always use the whole-study action for linked jobs, even if active;
        # this also avoids a check-then-cancel race between separate routes.
        request = job.get("request") or {}
        if job.get("kind") != "discovery_prepare" and not request.get("campaign_id"):
            return
        with self.queue.connection() as conn:
            linked = conn.execute("""SELECT id FROM platform_discovery_runs
                WHERE owner_id=%s AND (job_id=%s OR campaign_id=%s) LIMIT 1""",
                (owner, job["id"], request.get("campaign_id"))).fetchone()
        if linked:
            raise ValueError("Resume this research through its whole-study action after resource cleanup")

    def _control_discovery_run(self, run_id, action, *, owner):
        # This lock serializes user control requests across HTTP processes. It
        # does not share the preparation/activation row lock, so requesting
        # cancellation does not wait for a long scientific stage to finish.
        self.discovery._record(run_id, owner)
        with self.queue.connection() as conn:
            conn.execute("SELECT pg_advisory_xact_lock(hashtextextended(%s,0))",
                         ("research-control:" + owner + ":" + run_id,))
            row = self.discovery._record(run_id, owner)
            if action == "cancel":
                self.discovery.set_lifecycle(run_id, "cancelled", owner=owner)



                self.queue.cancel(row["job_id"], owner=owner)
                if row["campaign_id"]:
                    self.campaigns.cancel_run_jobs(row["campaign_id"], owner=owner)
                return self.discovery.get(run_id, owner=owner)
            if action == "continue":
                if row["desired_state"] != "active":
                    raise ValueError("Cancelled research requires explicit resume after resource cleanup")
                return self.discovery.continue_run(run_id, owner=owner)
            if action != "resume":
                raise ValueError("Unknown research control")
            current = self.discovery.get(run_id, owner=owner)
            cancellation = current.get("cancellation") or {}
            if row["desired_state"] == "cancelled" and (
                    cancellation.get("resources_released") is not True
                    or cancellation.get("resume_allowed") is not True):
                raise ValueError("Resume requires confirmed resource cleanup and remaining frozen budget")
            job = self.queue.get(row["job_id"], owner=owner)
            if row["campaign_id"]:
                attempt = self.campaigns.resume_run_jobs(row["campaign_id"], owner=owner)
                if attempt.get("status") == "resume_incomplete":
                    # A partial retry must not leave newly queued work running
                    # behind a cancelled study or erase the original holds.
                    self.campaigns.cancel_run_jobs(row["campaign_id"], owner=owner)
                    return {**self.discovery.get(run_id, owner=owner),
                            "resume_status": "resume_incomplete", "resume_attempt": attempt}
            elif job["status"] != "complete":
                self.queue.resume(row["job_id"], owner=owner)
            # A rejected queue/budget resume must never clear the durable
            # cancellation intent and reactivate downstream stages.
            self.discovery.set_lifecycle(run_id, "active", owner=owner)
            return self.discovery.continue_run(run_id, owner=owner)

    def _job_budget(self, job):
        """Expose server-authoritative allowances without allowing scope changes."""
        return self._jobs_budget([job])[0]

    def _jobs_budget(self, jobs):
        'Validate the selected structure and retain explicit calculation provenance.'
        if not jobs:
            return []
        with self.queue.connection() as conn:
            rows = conn.execute(
                'SELECT id,owner_id,limit_seconds,used_seconds FROM platform_compute_budgets '
                'WHERE id=ANY(%s) AND owner_id=ANY(%s)',
                (list({j.get('budget_scope') for j in jobs if j.get('budget_scope')}),
                 list({j.get('owner_id') for j in jobs if j.get('owner_id')})),
            ).fetchall()
        budgets = {(r['id'], r['owner_id']): r for r in rows}
        manifests, result = {}, []
        for job in jobs:
            request = job.get('request') or {}
            maximum = 3600 if job['kind'].startswith('algorithm_') or job['kind'] == 'result_export' else 28800
            if job['kind'] == 'literature_review':
                maximum = 900
            status, error_type = 'verified', None
            if request.get('campaign_id') and job.get('owner_id') and (job['kind'] == 'dock_batch' or job['kind'].startswith('algorithm_')):
                key = (job.get('owner_id'), request['campaign_id'])
                if key not in manifests:
                    try:
                        manifest = self.campaigns._coordinator(key[0])._manifest(key[1])
                        # Shape validation is isolated with the seal read too.
                        manifest['plan']['budget_scope']
                        manifest['plan']['budget_seconds']
                        manifests[key] = manifest
                    except (ArtifactIntegrityError, ValueError, KeyError, OSError) as error:
                        manifests[key] = error
                manifest = manifests[key]
                if isinstance(manifest, Exception):
                    status, error_type, maximum = 'manifest_quarantined', type(manifest).__name__, None
                else:
                    bound_scope = str(job['owner_id']) + ':' + manifest['plan']['budget_scope']
                    if job.get('budget_scope') == bound_scope:
                        maximum = manifest['plan']['budget_seconds']
                    elif job['kind'] == 'dock_batch':


                        status, maximum = 'budget_binding_conflict', None
            if job['kind'] == 'discovery_prepare':
                maximum = request.get('compute_budget_seconds')
            from .ranking import finite
            if status == 'verified' and (finite(maximum) is None or maximum <= 0):
                status, maximum = 'budget_unavailable', None
            budget = budgets.get((job.get('budget_scope'), job.get('owner_id')))
            if budget is None and status == 'verified':
                status, maximum = 'budget_unavailable', None
            result.append({**job, 'budget': {
                'limit_seconds': budget['limit_seconds'] if budget else None,
                'used_seconds': budget['used_seconds'] if budget else None,
                'maximum_seconds': maximum, 'status': status,
                'resume_allowed': job['kind'] != 'literature_review' and status == 'verified' and maximum is not None
                    and (job.get('compute_seconds', 0) < maximum)
                    and (budget['used_seconds'] < maximum),
                **({'error_type': error_type} if error_type else {}),
            }})
        return result

    def worker_get(self, url, node_id):
        parts = urlparse(url).path.strip("/").split("/")
        if len(parts) == 4 and parts[:3] == ["v2", "worker", "artifacts"]:
            raw = self.artifacts.get_bytes(parts[3])
            return {"sha256": parts[3], "base64": base64.b64encode(raw).decode()}
        raise KeyError("Unknown worker route")

    def worker_post(self, url, body, node_id):
        action = urlparse(url).path.removeprefix("/v2/worker/")
        if action == "register":
            self.queue.register_node(node_id, body["telemetry"])
            return {"registered": node_id, "claim_contracts": ["engine_hashes_v1", "cancellation_ack_v1", "owned_cleanup_ack_v1"]}
        if action == "claim":
            kinds = body.get("kinds", [])
            if (
                not isinstance(kinds, list)
                or any(not isinstance(k, str) for k in kinds)
                or not set(kinds)
                <= {
                    "probe",
                    "dock_batch",
                    "algorithm_evaluation",
                    *("algorithm_" + p for p in PROVIDERS),
                }
            ):
                raise ValueError("Unsupported worker kinds")
            # Operator-owned marker drains at a job boundary. It does not
            # cancel a lease, stop an engine, or alter completed checkpoints.
            marker = (
                self.root
                / "worker-drain"
                / hashlib.sha256(node_id.encode()).hexdigest()
            )
            if marker.exists():
                return {"job": None, "status": "maintenance_drain"}
            if "engine_hashes" in body:
                from .engines import validate_engine_hashes

                validate_engine_hashes(body["engine_hashes"])
                return {"job": self.queue.claim(node_id, kinds, engine_hashes=body["engine_hashes"])}
            return {"job": self.queue.claim(node_id, kinds)}
        if action == "artifacts":
            raw = self._decode_upload(body)
            if hashlib.sha256(raw).hexdigest() != body["sha256"]:
                raise ValueError("Artifact content hash mismatch")
            return self.artifacts.put_bytes(raw)
        job_id, token = body["job_id"], body["lease_token"]
        if action == 'cleanup-ack':
            return self.queue.acknowledge_cleanup(job_id, token, worker=node_id, cleanup=body.get('cleanup'))
        job = self.queue.get(job_id)
        if job.get("worker") != node_id:
            raise ValueError("Worker does not own job")
        if action == "heartbeat":
            return self.queue.heartbeat(
                job_id,
                token,
                progress=body.get("progress"),
                compute_seconds=body.get("compute_seconds"),
            )
        if action == "cancel-ack":
            cleanup = body.get("cleanup")
            if (not isinstance(cleanup, dict)
                    or cleanup.get("scope") != "owned_job_processes"
                    or cleanup.get("owned_processes_stopped") is not True):
                raise ValueError("Cancellation acknowledgement requires owned-process cleanup evidence")
            return self.queue.acknowledge_cancel(
                job_id, token, compute_seconds=body.get("compute_seconds"), cleanup=cleanup)
        if action == "finish":
            result = body["result"]
            if job["kind"] == "algorithm_evaluation":
                if (
                    not isinstance(result, dict)
                    or result.get("provider") != "evaluation"
                    or result.get("promoted") is not False
                    or result.get("ranking_weight") != 0
                    or result.get("evaluation_manifest_sha256")
                    != job["request"]["evaluation_manifest_sha256"]
                ):
                    raise ValueError(
                        "Evaluation result identity or promotion state mismatch"
                    )
                raw = self.artifacts.get_bytes(result["result_artifact"]["sha256"])
                if raw != canonical_json(
                    {
                        k: v
                        for k, v in result.items()
                        if k not in {"result_artifact", "worker_telemetry"}
                    }
                ):
                    raise ValueError(
                        "Evaluation result differs from its immutable publication"
                    )
            elif job["kind"].startswith("algorithm_"):
                self.algorithm_jobs.validate_result(job["request"], result)
            if job["kind"] == "dock_batch":
                expected = {r["candidate_id"] for r in job["request"]["ligands"]}
                rows = result.get("candidates", [])
                if (
                    len(rows) != len(expected)
                    or {r.get("candidate_id") for r in rows} != expected
                ):
                    raise ValueError(
                        "Docking receipt must account for every input candidate"
                    )
                for row in rows:
                    if row.get("pose_sha256"):
                        self.artifacts.verified_path(row["pose_sha256"])
            return self.queue.finish(
                job_id, token, result, compute_seconds=body.get("compute_seconds"),
                **({"cleanup": body["cleanup"]} if "cleanup" in body else {}),
            )
        if action == "fail":
            self.queue.fail(
                job_id, token, body["error"], status=body.get("status", "failed")
            )
            return {"recorded": True}
        if action == "defer":
            self.queue.defer(job_id, token, body["reason"])
            return {"deferred": True}
        if action == "reserve":
            return {"reserved": self.queue.reserve_gpu(job_id, token, body["gpu_uuid"])}
        if action == "reserve-cpu":
            return {"reserved": self.queue.reserve_cpu(job_id, token, body["cores"])}
        if action == "reserve-model-launch":
            return {
                "reserved": self.queue.reserve_model_launch(
                    job_id, token, body["work_key"]
                )
            }
        if action == "start-model-launch":
            return {
                "authorized": self.queue.start_model_launch(
                    job_id, token, body["work_key"]
                )
            }
        raise KeyError("Unknown worker route")

    def import_legacy(self, screen_id, *, activate=True, expected_source_sha256=None, owner="researcher"):
        started = time.monotonic()

        def progress(phase, **details):
            callback = getattr(self, "import_progress", None)
            if callback:
                callback(
                    {
                        "screen_id": screen_id,
                        "phase": phase,
                        "elapsed_seconds": round(time.monotonic() - started, 3),
                        **details,
                    }
                )

        if self.bridge is None:
            raise ValueError("Legacy bridge unavailable")
        progress("source_read_started")
        source = self.bridge.screen_file(screen_id)
        source_hash = sha256_file(source)
        if expected_source_sha256 is not None and source_hash != expected_source_sha256:
            raise ValueError("Legacy source differs from the frozen inventory")
        _, data = self.bridge.results.read(source)
        package = self.bridge.results.package(source)
        progress("source_read_complete")
        scientific_hash = (
            package["manifest"]["source_sha256"] if package else source_hash
        )
        snapshot = (
            "legacy-"
            + hashlib.sha256(
                ("legacy-adapter-v2:" + owner + ":" + screen_id + source_hash).encode()
            ).hexdigest()[:24]
        )
        if (self.catalog.root / snapshot).exists():
            manifest = self.catalog.manifest(snapshot)
            if self.catalog.snapshot_owner(snapshot) != owner:
                raise PermissionError("Existing legacy snapshot belongs to another project")
            if manifest.get("provenance", {}).get("source_sha256") != scientific_hash:
                raise ValueError(
                    "Existing legacy snapshot has a conflicting scientific source"
                )
            if activate and owner == "researcher":
                self._activate_snapshot(screen_id, snapshot)
            progress(
                "verified_snapshot_reused", candidate_count=manifest["candidate_count"]
            )
            return manifest
        triage = data.get("triage") or {}
        node = triage.get("node") if isinstance(triage, dict) else None
        if isinstance(node, dict):
            target = node.get("gene") or node.get("uniprot") or node.get("name")
        else:
            target = node
        target = target or data.get("marker") or screen_id
        if not isinstance(target, str) or not target:
            raise ValueError("Legacy target identity must be an explicit string")
        if isinstance(data.get("discovery_candidates"), list):
            originals = data["discovery_candidates"]
        else:
            # Older results may record excluded rows only under dropped. Keep
            # every recoverable ID; never replace its missing values with zeros.
            by_id = {}
            for old in data.get("candidates", []):
                cid = old.get("ligand_id")
                if cid in by_id:
                    raise ValueError("Duplicate legacy candidate identity")
                by_id[cid] = old
            for old in data.get("dropped", []):
                cid = old.get("ligand_id")
                if cid and cid not in by_id:
                    by_id[cid] = old
            originals = list(by_id.values())
        rows = []
        for old in originals:
            row = dict(old)
            cid = old.get("ligand_id")
            if not isinstance(cid, str) or not cid:
                raise ValueError("Legacy candidate has no recoverable ligand ID")
            row.update(
                candidate_id=old.get("ligand_id"),
                target_id=target,
                screen_id=screen_id,
                marker=data.get("marker"),
                target_pool_membership=True,
                target_association=old.get("target_association") is True
                or bool(old.get("target_evidence")),
                legacy_detail_url=f"/v1/screens/{screen_id}/candidates/{old.get('ligand_id')}",
                legacy_docking_url=f"/v1/screens/{screen_id}/docking/{old.get('ligand_id')}",
                source_snapshot_sha256=scientific_hash,
                legacy_source_annotation=dict(old),
            )
            row["opposing_direction"] = row.get("direction_verdict") == "opposing"
            if row["opposing_direction"]:
                row["role"] = "counter_control"
            # Original values remain available. Missing assay/engine context
            # never gets invented to turn a legacy score into ranked evidence.
            # An old scalar docking field is retained as an annotation, not
            # silently reinterpreted as an admitted v2 structural observation.
            if not isinstance(row.get("measurements"), (list, dict)):
                row["measurements"] = []
            if not isinstance(row.get("docking"), (list, dict)):
                row["docking"] = []
            rows.append(row)
        progress("rows_constructed", candidate_count=len(rows))
        ranked = rank_candidates(rows)
        progress("ranking_complete")
        for index, row in enumerate(ranked):
            ranked[index] = assess_candidate(
                row, target_identity_valid=row.get("target_identity_valid") is True
            )
            if index % 100 == 99:
                progress(
                    "structure_assessment",
                    assessed=index + 1,
                    candidate_id=row["candidate_id"],
                )
        progress("structure_assessment_complete")
        if sha256_file(source) != source_hash:
            raise ValueError("Legacy source changed during import")
        if package:
            self.bridge.results.package(source)
        progress("source_reverified")
        manifest = self.catalog.create_snapshot(
            ranked,
            snapshot_id=snapshot,
            provenance={
                "owner_id": owner,
                "screen_id": screen_id,
                "marker": data.get("marker"),
                "target_id": target,
                "source_sha256": scientific_hash,
                "source_file_sha256": source_hash,
                "source_candidate_count": len(rows),
                "source_candidate_ids_sha256": hashlib.sha256(
                    canonical_json(sorted(r["candidate_id"] for r in rows))
                ).hexdigest(),
                "package_manifest_sha256": hashlib.sha256(
                    canonical_json(package["manifest"])
                ).hexdigest()
                if package
                else None,
                "legacy": True,
                "adapter_version": "legacy-adapter-v2",
                "references": [],
                "ranking_note": "Legacy scores lacking complete comparability context remain source annotations.",
            },
        )
        progress("snapshot_published", candidate_count=manifest["candidate_count"])
        if activate and owner == "researcher":
            self._activate_snapshot(screen_id, snapshot)
        return manifest

    def _activate_snapshot(self, screen_id, snapshot):
        self.catalog.manifest(snapshot)
        self._update_active_snapshots({screen_id: snapshot})

    def _update_active_snapshots(self, updates):
        import fcntl

        with self.lock:
            with (self.root / ".active-snapshots.lock").open("a") as lock:
                fcntl.flock(lock, fcntl.LOCK_EX)
                path = self.root / "active-snapshots.json"
                mapping = json.loads(path.read_text()) if path.exists() else {}
                atomic_bytes(path, canonical_json({**mapping, **updates}))
                return mapping
