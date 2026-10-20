'Owner-only literature jobs with a private, explicitly pinned registration.'
from __future__ import annotations

from copy import deepcopy
import hashlib
from pathlib import PurePosixPath
import re
from urllib.parse import parse_qs, unquote, urlparse

from .artifacts import ArtifactIntegrityError, canonical_json
from .literature import collect_for_run, freeze_policy
from .queue import identifier
from .targets import target_request

KIND = "literature_review"
DEFAULT_MODEL = "gpt-6-luna"
DEFAULT_EFFORT = "max"
MAX_SOURCE_BYTES = 4 * 1024 * 1024
_MODEL = re.compile(r"[A-Za-z0-9][A-Za-z0-9._:-]{0,127}\Z")
_REVIEW = re.compile(r"[0-9a-f]{32}\Z")
_PARAMETERS = ("query", "kind", "taxon_id", "action", "max_papers", "model", "reasoning_effort", "public_topic")


def _safe_error(error):
    # Reuse the collector's allowlisted quota/timing metadata, never str(error).
    from ..literature.collect import _review_failure
    return _review_failure(error)["error"]


class PersonalLiterature:
    def __init__(self, service, config=None):
        self.service, self.queue = service, service.queue
        self.config = deepcopy(config or {})
        value = self.config.get("literature", {})
        self.configured = (isinstance(value, dict) and value.get("enabled") is True
            and value.get("review_provider") == "chatgpt"
            and value.get("deployment_mode") in {"personal_self_hosted", "local_single_user"})
        if self.configured:
            value.setdefault("model", DEFAULT_MODEL)
            value.setdefault("reasoning_effort", DEFAULT_EFFORT)
            value.setdefault("max_papers", 2)
            freeze_policy(self.config)
            identifier(value["owner_id"], "literature owner")

    def _guard(self, owner):
        if not self.configured:
            raise PermissionError("Personal literature review is not configured")
        if owner != self.config["literature"]["owner_id"]:
            raise PermissionError("Personal literature belongs to another owner")

    def _client(self, *, model=None, reasoning_effort=None):
        from ..literature.chatgpt import ChatGPTClient
        value = self.config["literature"]
        options = {"model": model or value["model"],
                   "reasoning_effort": reasoning_effort or value["reasoning_effort"]}
        if value["deployment_mode"] == "personal_self_hosted":
            options["account_id"] = value["account_id"]
        return ChatGPTClient(value["state_dir"], **options)

    def status(self, *, owner):
        value = self.config.get("literature", {}) if self.configured else {}
        result = {"configured": self.configured, "connected": False, "plan_enabled": False,
                  "mode": value.get("deployment_mode", "unconfigured"),
                  "default_model": value.get("model", DEFAULT_MODEL),
                  "default_reasoning_effort": value.get("reasoning_effort", DEFAULT_EFFORT),
                  "default_max_papers": 2, "max_papers": 8,
                  "review_timeout_seconds": value.get("review_timeout", 180),
                  "status_scope": "saved selected registration; actual request decides eligibility",
                  "scientific_status": "NOT_ESTABLISHED"}
        if not self.configured:
            return result
        self._guard(owner)
        from ..literature.chatgpt import ChatGPTError
        try:
            with self._client() as client:
                saved = client.status()
            result.update(connected=saved.get("connected") is True,
                          plan_enabled=saved.get("plan_enabled") is True)
        except (OSError, ValueError, ChatGPTError) as error:
            result["error"] = _safe_error(error)
        return result

    def usage(self, *, owner):
        self._guard(owner)
        with self._client() as client:
            usage = client.usage()
        # The client is registration-pinned. Its accounting is local records,
        # never a provider plan balance or account-wide consumption estimate.
        account = usage.get("account", {})
        return {key: deepcopy(usage.get(key)) for key in (
            "schema", "scope", "provider_plan", "app_usage", "last_limit_error", "caveat")} | {
                "account": {key: deepcopy(account.get(key)) for key in (
                    "connected", "plan_permission", "status_scope")}}

    def models(self, *, owner):
        self._guard(owner)
        from ..literature.chatgpt import ChatGPTError
        try:
            with self._client() as client:
                catalog = client.list_models()
            models = [{key: row[key] for key in ("id", "slug", "display_name", "supported_reasoning_efforts", "default_reasoning_effort", "is_default") if key in row}
                      for row in catalog["models"] if isinstance(row, dict)
                      and isinstance(row.get("id"), str) and _MODEL.fullmatch(row["id"])]
            result = {"models": models, "entitlement_established": False,
                    "note": "Catalog visibility is not an inference entitlement"}
            if isinstance(catalog.get("default_model"), str) and catalog["default_model"] in {row["id"] for row in models}:
                result["default_model"] = catalog["default_model"]
            return result
        except (OSError, ValueError, ChatGPTError) as error:
            return {"models": [], "status": "HOLD", "error": _safe_error(error),
                    "entitlement_established": False}

    def _effective_config(self, parameters):
        config = deepcopy(self.config)
        config["literature"].update({key: parameters[key] for key in (
            "max_papers", "model", "reasoning_effort")})
        return config

    def submit(self, body, *, owner):
        self._guard(owner)
        if not isinstance(body, dict) or set(body) - set(_PARAMETERS) - {"idempotency_key"}:
            raise ValueError("Literature review accepts target search options only")
        key = body.get("idempotency_key")
        if not isinstance(key, str) or not re.fullmatch(r"[A-Za-z0-9_.:-]{1,128}", key):
            raise ValueError("An explicit bounded idempotency_key is required")
        compound = body.get("kind") == "compound"
        if compound:
            from ..literature.collect import public_compound_name, public_skin_topic
            query = public_compound_name(body.get("query"))
            topic = public_skin_topic(body.get("public_topic", ""))
            taxon = body.get("taxon_id", 9606)
            if type(taxon) is not int or not 1 <= taxon <= 100_000_000:
                raise ValueError("taxon_id must be a positive integer")
            target = {"query": query, "kind": "compound", "taxon_id": taxon}
        else:
            if "public_topic" in body:
                raise ValueError("Public skin topic is available for compound literature only")
            target = target_request({"query": body.get("query"), "kind": body.get("kind", "protein"),
                                     "taxon_id": body.get("taxon_id", 9606)})
        action = body.get("action", "any")
        if not isinstance(action, str) or action not in {
            "any", "inhibit", "activate", "increase", "decrease", "suppress", "unknown"}:
            raise ValueError("Unsupported literature action")
        count = body.get("max_papers", 2)
        if type(count) is not int or not 1 <= count <= 8:
            raise ValueError("max_papers must be within 1..8")
        model, effort = body.get("model", DEFAULT_MODEL), body.get("reasoning_effort", DEFAULT_EFFORT)
        if not isinstance(model, str) or not _MODEL.fullmatch(model):
            raise ValueError("A bounded model identifier is required")
        if not isinstance(effort, str) or effort not in {"none", "low", "medium", "high", "xhigh", "max"}:
            raise ValueError("Unsupported reasoning_effort")
        parameters = {"query": target["query"], "kind": target["kind"], "taxon_id": target["taxon_id"],
                      "action": action, "max_papers": count, "model": model, "reasoning_effort": effort}
        signature = hashlib.sha256(canonical_json({"owner": owner, "key": key})).hexdigest()
        public_request_sha = None
        if compound:
            from ..literature.collect import PUBLIC_COMPOUND_SCHEMA
            public_request_sha = hashlib.sha256(canonical_json({**parameters, "public_topic": topic})).hexdigest()
            # A retried public request reuses its first immutable identity
            # evidence, even if the name provider is now unavailable or its
            # response formatting changed. No second model job is submitted.
            with self.queue.connection() as connection:
                previous = connection.execute("SELECT id,request FROM platform_jobs WHERE owner_id=%s AND idempotency_key=%s AND kind=%s",
                    (owner, "literature-" + signature, KIND)).fetchone()
            if previous is not None:
                if previous["request"].get("public_request_sha256") != public_request_sha:
                    raise ValueError("Idempotency key already bound to different public literature input")
                return self.get_review(previous["id"], owner=owner)
            resolved = self.service.research_inputs.resolve({"kind": "compound", "format": "name", "value": query}, owner=owner)
            cid = resolved.get("pubchem_cid")
            if (resolved.get("status") != "verified" or resolved.get("kind") != "compound"
                    or type(cid) is not int or cid < 1
                    or (resolved.get("provenance") or {}).get("identity_basis") != "provider_name_to_structure"):
                raise ValueError("화합물 이름을 하나의 PubChem 구조로 확인하지 못했습니다. 구체적인 개별 물질 이름을 확인해 주세요.")
            name = public_compound_name(resolved.get("name"))
            receipts = [r["raw_sha256"] for r in resolved.get("source_receipts", [])
                        if isinstance(r, dict) and r.get("status") == "complete"
                        and isinstance(r.get("raw_sha256"), str) and re.fullmatch(r"[0-9a-f]{64}", r["raw_sha256"])]
            if not receipts:
                raise ValueError("화합물 이름 확인의 원본 조회 근거가 없습니다. 다시 확인해 주세요.")
            target = {"schema": PUBLIC_COMPOUND_SCHEMA, "kind": "compound", "query": name, "pubchem_cid": cid,
                      "public_topic": topic, "taxon_id": taxon, "identity_status": "NAME_RESOLVED_BY_PUBCHEM",
                      "identity_source_sha256": sorted(set(receipts))}
            parameters.update(query=name, public_topic=topic, pubchem_cid=cid)
        target["action"] = action
        config = self._effective_config(parameters)
        request = {"schema": "skinpipe-personal-literature-request-v1", "parameters": parameters,
                   "target": target, "literature_policy": freeze_policy(config)}
        if compound:
            request["public_request_sha256"] = public_request_sha
        job = self.queue.enqueue(KIND, request, owner=owner,
            project="literature-" + hashlib.sha256(owner.encode()).hexdigest()[:24],
            idempotency_key="literature-" + signature, budget_seconds=900, max_attempts=1)
        return self.get_review(job["id"], owner=owner)

    def _job(self, review_id, owner):
        self._guard(owner)
        if not isinstance(review_id, str) or not _REVIEW.fullmatch(review_id):
            raise ValueError("Invalid literature review ID")
        job = self.queue.get(review_id, owner=owner)
        if job.get("kind") != KIND:
            raise KeyError("Unknown literature review")
        return job

    def _result(self, job):
        result = job.get("result")
        if not isinstance(result, dict) or not isinstance(result.get("dossier"), dict):
            return None
        raw = self.service.artifacts.get_bytes(result["dossier_sha256"])
        if raw != canonical_json(result["dossier"]):
            raise ArtifactIntegrityError("Literature dossier binding mismatch")
        return result

    def get_review(self, review_id, *, owner):
        job = self._job(review_id, owner)
        request = job.get("request", {}).get("parameters", {})
        result = {"review_id": job["id"], "job_id": job["id"],
                  "request": {key: request[key] for key in (*_PARAMETERS, "pubchem_cid") if key in request},
                  "result": None, "scope": "literature_only", "scientific_status": "NOT_ESTABLISHED"}
        for key in ("status", "created_at", "updated_at", "compute_seconds", "budget_seconds", "progress"):
            result[key] = deepcopy(job.get(key))
        if job.get("error"):
            result["error"] = {"code": "literature_job_failed"}
        if job.get("status") == "complete":
            try:
                result["result"] = deepcopy(self._result(job))
            except (KeyError, OSError, ValueError):
                result["error"] = {"code": "literature_dossier_unavailable"}
        return result

    def list(self, *, owner):
        self._guard(owner)
        with self.queue.connection() as connection:
            rows = connection.execute("SELECT id FROM platform_jobs WHERE owner_id=%s AND kind=%s "
                "ORDER BY created_at DESC,id DESC LIMIT 100", (owner, KIND)).fetchall()
        return {"reviews": [self.get_review(row["id"], owner=owner) for row in rows], "scope": "literature_only"}

    def source(self, review_id, path, *, owner):
        job = self._job(review_id, owner)
        if job["status"] != "complete":
            raise ValueError("Literature review has no completed dossier")
        result = self._result(job)
        if result is None:
            raise ValueError("Literature dossier unavailable")
        if (not isinstance(path, str) or not path or len(path) > 512 or "\\" in path
                or PurePosixPath(path).is_absolute() or any(p in {".", ".."} for p in path.split("/"))):
            raise ValueError("An exact bound relative source path is required")
        matches = [item for item in result["dossier"].get("literature_artifacts", []) if item.get("path") == path]
        if len(matches) != 1:
            raise ValueError("Source path is not uniquely bound to this dossier")
        binding = matches[0]
        size = binding.get("size_bytes")
        if type(size) is not int or not 0 <= size <= MAX_SOURCE_BYTES:
            raise ValueError("Source exceeds the 4 MiB display limit")
        digest = binding["literature_artifact_sha256"]
        blob = self.service.artifacts.verified_path(digest)
        if blob.stat().st_size != size:
            raise ValueError("Source size differs from its dossier binding")
        with blob.open("rb") as stream:
            raw = stream.read(MAX_SOURCE_BYTES + 1)
        if len(raw) != size or hashlib.sha256(raw).hexdigest() != digest:
            raise ValueError("Source size differs from its dossier binding")
        try:
            text = raw.decode("utf-8")
        except UnicodeDecodeError:
            raise ValueError("Source is not UTF-8 text") from None
        return {"text": text, "sha256": binding["literature_artifact_sha256"], "size_bytes": len(raw), "path": path}

    def get(self, url, *, owner):
        parsed = urlparse(url)
        parts = [unquote(p) for p in parsed.path.strip("/").split("/")]
        if parts == ["v2", "literature", "status"]:
            return self.status(owner=owner)
        if parts == ["v2", "literature", "usage"]:
            return self.usage(owner=owner)
        if parts == ["v2", "literature", "models"]:
            return self.models(owner=owner)
        if parts == ["v2", "literature", "reviews"]:
            return self.list(owner=owner)
        if len(parts) == 4 and parts[:3] == ["v2", "literature", "reviews"]:
            return self.get_review(parts[3], owner=owner)
        if len(parts) == 5 and parts[:3] == ["v2", "literature", "reviews"]:
            if parts[4] == "source":
                query = parse_qs(parsed.query, keep_blank_values=True)
                if set(query) != {"path"} or len(query["path"]) != 1:
                    raise ValueError("One exact bound source path is required")
                return self.source(parts[3], query["path"][0], owner=owner)
            if parts[4] == "dossier":
                job = self._job(parts[3], owner)
                result = self._result(job) if job["status"] == "complete" else None
                if result is None:
                    raise ValueError("Literature dossier unavailable")
                return deepcopy(result["dossier"])
        raise KeyError("Unknown personal literature route")

    def execute(self, request, *, owner, stop, progress):
        self._guard(owner)
        stop()
        progress({"stage": "literature_evidence"})
        config = self._effective_config(request["parameters"])
        dossier = collect_for_run(self.service, request["target"], request, config,
                                  owner=owner, stop=stop, progress=progress)
        stop()
        sha = self.service.artifacts.put_bytes(canonical_json(dossier), media_type="application/json")["sha256"]
        return {"dossier": dossier, "dossier_sha256": sha, "scope": "literature_only",
                "scientific_status": "NOT_ESTABLISHED"}
