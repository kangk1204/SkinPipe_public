"""Optional, source-bound literature sidecar for newly submitted runs.

The sidecar never changes candidate scores, target selection or scientific
validation. Historical requests with no policy remain historical requests.
"""
from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import uuid

from .artifacts import canonical_json

POLICY = "skinpipe-literature-policy-v1"


def freeze_policy(config):
    value = (config or {}).get("literature", {})
    if not isinstance(value, dict) or set(value) - {
        "enabled", "max_papers", "timeout", "review_provider", "model",
        "state_dir", "owner_id", "deployment_mode",
        "reasoning_effort", "account_id", "binding_id", "review_timeout",
    }:
        raise ValueError("Unsupported literature configuration")
    if not isinstance(value.get("enabled", False), bool):
        raise ValueError("Literature enabled must be a boolean")
    if value.get("enabled", False) is not True:
        return {"schema": POLICY, "enabled": False}
    if value.get("review_provider") == "chatgpt" and value.get("deployment_mode") == "personal_self_hosted":
        value = {"model": "gpt-6-luna", "reasoning_effort": "max", "max_papers": 2, **value}
    count, timeout = value.get("max_papers", 8), value.get("timeout", 20)
    if isinstance(count, bool) or not isinstance(count, int) or not 1 <= count <= 8:
        raise ValueError("Literature max_papers must be within 1..8")
    if isinstance(timeout, bool) or not isinstance(timeout, int) or not 1 <= timeout <= 60:
        raise ValueError("Literature timeout must be within 1..60 seconds")
    review_timeout = value.get("review_timeout", 180)
    if type(review_timeout) is not int or not 1 <= review_timeout <= 300:
        raise ValueError("Literature review_timeout must be within 1..300 seconds")
    provider = value.get("review_provider", "none")
    if provider not in {"none", "chatgpt"}:
        raise ValueError("Literature review_provider must be none or chatgpt")
    effort = value.get("reasoning_effort")
    if effort is not None and (not isinstance(effort, str) or effort not in {"none", "low", "medium", "high", "xhigh", "max"}):
        raise ValueError("Unsupported literature reasoning_effort")
    if provider == "chatgpt":
        if value.get("deployment_mode") not in {"local_single_user", "personal_self_hosted"} or not value.get("owner_id"):
            raise ValueError("ChatGPT plan use requires a bound local_single_user or personal_self_hosted owner")
        if not isinstance(value.get("model"), str) or not value["model"].strip():
            raise ValueError("Choose a ChatGPT literature review model explicitly")
        if not value.get("state_dir"):
            raise ValueError("Configure a private ChatGPT state_dir outside Git")
        if value.get("deployment_mode") == "personal_self_hosted":
            import re
            if (not isinstance(value.get("account_id"), str) or not value["account_id"]
                    or len(value["account_id"]) > 256):
                raise ValueError("Personal hosting requires an explicit private account_id")
            if not isinstance(value.get("binding_id"), str) or not re.fullmatch(r"[0-9a-f]{32,64}", value["binding_id"]):
                raise ValueError("Personal hosting requires a random nonidentity binding_id")
    # Credential paths are runtime details, never part of an exported request.
    policy = {"schema": POLICY, "enabled": True, "max_papers": count,
            "timeout": timeout, "review_provider": provider,
            "model": value.get("model") if provider == "chatgpt" else None,
            "owner_id": value.get("owner_id") if provider == "chatgpt" else None,
            "deployment_mode": value.get("deployment_mode") if provider == "chatgpt" else None}
    if "reasoning_effort" in value:
        policy["reasoning_effort"] = effort
    # Missing fields keep historical policy bytes and the old 180s behavior.
    # New installations explicitly choose 300s independently of retrieval.
    if "review_timeout" in value:
        policy["review_timeout"] = review_timeout
    if provider == "chatgpt" and value.get("deployment_mode") == "personal_self_hosted":
        policy["binding_id"] = value["binding_id"]
    return policy


def collect_for_run(service, descriptor, request, config, *, owner, stop, collector=None, progress=None):
    policy = request.get("literature_policy")
    if not policy or policy.get("enabled") is not True:
        return {"schema": "skinpipe-literature-evidence-v1", "status": "NOT_RUN",
                "reason": "No literature policy was frozen with this request",
                "human_review_required": True,
                "human_review_checklist": ["Literature has not been collected for this run"],
                "experimental_validation": "NOT_EVALUATED"}
    if policy != freeze_policy(config):
        return {"schema": "skinpipe-literature-evidence-v1", "status": "HOLD",
                "reason": "Literature policy changed after submission; resubmit a new run",
                "human_review_required": True,
                "human_review_checklist": ["Operator must restore the frozen policy or create a new run"],
                "experimental_validation": "NOT_EVALUATED"}
    if collector is None:
        from ..literature.collect import collect_evidence
        collector = collect_evidence
    reviewer = None
    if policy["review_provider"] == "chatgpt":
        def reviewer(packet):
            if policy["owner_id"] != owner:
                raise ValueError("ChatGPT account is not bound to this run owner")
            stop()
            if progress is not None:
                progress({"stage": "literature_review",
                          "review_timeout_seconds": policy.get("review_timeout", 180)})
            from ..literature.chatgpt import ChatGPTClient
            options = {"model": policy["model"]}
            if policy.get("reasoning_effort") is not None:
                options["reasoning_effort"] = policy["reasoning_effort"]
            if policy.get("deployment_mode") == "personal_self_hosted":
                options["account_id"] = config["literature"]["account_id"]
            with ChatGPTClient(config["literature"]["state_dir"], **options) as client:
                options = {"timeout": policy["review_timeout"]} if "review_timeout" in policy else {}
                result = client.review(packet, **options)
            stop()
            return result
    stop()
    directory = Path(service.root) / "literature" / ("run-" + uuid.uuid4().hex)
    try:
        dossier = collector(descriptor, directory, reviewer=reviewer,
                            max_papers=policy["max_papers"], timeout=policy["timeout"])
    except (ValueError, OSError) as error:
        return {"schema": "skinpipe-literature-evidence-v1", "status": "HOLD",
                "reason": "Literature collection unavailable: " + type(error).__name__,
                "human_review_required": True,
                "human_review_checklist": ["Inspect the preserved collection attempt"],
                "experimental_validation": "NOT_EVALUATED"}
    stop()
    dossier = deepcopy(dossier)
    dossier.pop("artifact_dir", None)
    # Import the entire source/decision receipt closure into the existing CAS.
    # These references are recognized by the ordinary immutable ZIP exporter.
    artifacts = []
    try:
        # Check the complete packet before importing any bytes. The root must
        # also be checked: rglob does not yield its starting directory.
        if directory.is_symlink():
            raise ValueError("Literature artifact symlinks are forbidden")
        paths = sorted(directory.rglob("*"))
        if any(path.is_symlink() for path in paths):
            raise ValueError("Literature artifact symlinks are forbidden")
        for path in paths:
            if path.is_file():
                stored = service.artifacts.put_file(path, media_type="application/octet-stream")
                artifacts.append({"path": path.relative_to(directory).as_posix(),
                                  "literature_artifact_sha256": stored["sha256"],
                                  "size_bytes": stored["size_bytes"]})
    except (ValueError, OSError) as error:
        # A failed optional sidecar must not discard successful candidate work.
        # Preserve the local packet and any already-published CAS bytes; expose
        # no claims as source-bound when the receipt closure is incomplete.
        dossier["collection_status"] = dossier.get("status", "UNKNOWN")
        dossier.update(status="HOLD", claims=[], human_review_required=True,
                       experimental_validation="NOT_EVALUATED",
                       scientific_status="NOT_ESTABLISHED",
                       reason="Literature artifact publication unavailable: " + type(error).__name__)
        if isinstance(dossier.get("review"), dict):
            dossier["review"]["claims"] = []
        dossier.setdefault("checks", {})["artifact_publication"] = {
            "status": "HOLD", "error_type": type(error).__name__,
            "published_files": len(artifacts),
        }
    else:
        dossier.setdefault("checks", {})["artifact_publication"] = {
            "status": "COMPLETE", "published_files": len(artifacts),
        }
    dossier["literature_artifacts"] = artifacts
    dossier["policy"] = policy
    # Validate serializability before the catalog snapshot becomes visible.
    json.loads(canonical_json(dossier))
    return dossier
