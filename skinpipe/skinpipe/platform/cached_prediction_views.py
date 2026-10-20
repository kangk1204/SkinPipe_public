"""Bounded read-only display metadata for an already owner-authorized snapshot.

Stored rows deliberately deduplicate model cards into CAS. Restore them only
for selected summary cards or a single candidate, never by changing source
predictions, snapshots, exports, rankings or the additional-model policy.
"""
from copy import deepcopy
import hashlib
import json
import re

from .artifacts import canonical_json
from .identity import chemical_identity

_SHA = re.compile(r"[0-9a-f]{64}")
_FIELDS = ("model_manifest", "endpoint_definitions")


class CachedPredictionView:
    """One response owns its bounded caches; no cross-owner/model fallback."""

    def __init__(self, artifacts):
        self.artifacts = artifacts
        self.cache = {}
        self.bytes_read = 0

    def _json(self, sha, limit):
        if not isinstance(sha, str) or not _SHA.fullmatch(sha):
            raise ValueError("Invalid cached metadata reference")
        if sha not in self.cache:
            # Check size before ArtifactStore performs its full byte hash.
            blob, _ = self.artifacts._paths(sha)
            size = blob.stat().st_size
            if size > limit or self.bytes_read + size > 64 * 1024**2 or len(self.cache) >= 64:
                raise ValueError("Cached metadata display budget exceeded")
            self.bytes_read += size
            self.cache[sha] = json.loads(self.artifacts.get_bytes(sha))
        return self.cache[sha]

    @staticmethod
    def _reference(prediction, field):
        sha = prediction.get(field + "_sha256")
        refs = prediction.get("metadata_refs")
        ref = refs.get(field + "_sha256") if isinstance(refs, dict) else None
        if (not isinstance(sha, str) or not _SHA.fullmatch(sha) or not isinstance(ref, dict)
                or ref.get("sha256") != sha or ref.get("source") != "SkinScout model registry"
                or ref.get("storage") not in {"platform_cas", "external_registry"}
                or ref.get("digest_format") not in {"canonical_json", "canonical_json_utf8_v1", "canonical_json_ascii_v1"}):
            raise ValueError("Unverified cached metadata provenance")
        return sha, ref

    @staticmethod
    def _content_matches(value, sha, reference):
        raw = (json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
            allow_nan=False).encode() if reference["digest_format"] == "canonical_json_ascii_v1"
            else canonical_json(value))
        return hashlib.sha256(raw).hexdigest() == sha

    def _source(self, row, prediction):
        chunk_sha = row.get("source_chunk_sha256")
        if not chunk_sha:
            # Exact-lookup query seeds retain their full, server-pinned CAS
            # definitions. An external hash alone cannot borrow a global card.
            if not all(self._reference(prediction, field)[1]["storage"] == "platform_cas"
                    for field in _FIELDS):
                raise ValueError("Shared external definitions require the saved source chunk")
            return None
        source = self._json(chunk_sha, 16 * 1024**2)
        if (not isinstance(source, dict) or source.get("status") != "complete"
                or not isinstance(source.get("source"), dict)
                or not isinstance(source["source"].get("report_sha256"), str)
                or not _SHA.fullmatch(source["source"]["report_sha256"])
                or source["source"].get("structure_admission") != "reported_verified"
                or source["source"].get("report_sha256") != prediction.get("source_audit_sha256")
                or not isinstance(source.get("rows"), list) or len(source["rows"]) > 1000):
            raise ValueError("Cached prediction source binding unavailable")
        for detail in source["rows"]:
            if not isinstance(detail, dict) or detail.get("compound_id") != prediction.get("compound_id"):
                continue
            compound = detail.get("compound")
            if not isinstance(compound, dict) or compound.get("compound_id") != prediction.get("compound_id"):
                continue
            for original in detail.get("admet", []):
                if not isinstance(original, dict):
                    continue
                keys = ("model_id", "compound_id", "input_smiles", "input_hash", "status", "values",
                        "model_manifest_sha256", "endpoint_definitions_sha256", "source_audit_sha256")
                if all(original.get(key) == prediction.get(key) for key in keys):
                    return source
        raise ValueError("Cached prediction differs from its saved source")

    def _prediction(self, row, prediction):
        identity = chemical_identity(row.get("smiles"))
        original = chemical_identity(prediction.get("input_smiles"))
        if (identity["structure_status"] != "valid" or original["structure_status"] != "valid"
                or identity["exact_structure_id"] != original["exact_structure_id"]):
            raise ValueError("Cached prediction input differs from the displayed molecule")
        source = self._source(row, prediction)
        values = {}
        for field in _FIELDS:
            sha, reference = self._reference(prediction, field)
            value = self._json(sha, 4 * 1024**2)
            if not self._content_matches(value, sha, reference):
                raise ValueError("Cached metadata content/encoding mismatch")
            # CAS-only external rows are allowed only when the same immutable
            # source response explicitly binds these definitions to this model.
            if reference["storage"] == "external_registry":
                matched = False
                for detail in source["rows"]:
                    for shared in detail.get("admet", []) if isinstance(detail, dict) else []:
                        if (isinstance(shared, dict) and shared.get("model_id") == prediction.get("model_id")
                                and shared.get(field + "_sha256") == sha and field in shared
                                and self._content_matches(shared[field], sha, reference)
                                and shared[field] == value):
                            matched = True
                if not matched:
                    raise ValueError("Shared definition belongs to another response/model")
            values[field] = value
        model, endpoints = values["model_manifest"], values["endpoint_definitions"]
        if (not isinstance(model, dict) or model.get("model_id") != prediction.get("model_id")
                or not isinstance(endpoints, list) or not 1 <= len(endpoints) <= 200
                or any(not isinstance(entry, dict) or not isinstance(entry.get("id"), str)
                    or not 1 <= len(entry["id"]) <= 200 or entry.get("task_type") not in {"classification", "regression"}
                    for entry in endpoints)
                or len({entry["id"] for entry in endpoints}) != len(endpoints)):
            raise ValueError("Cached model and endpoint definitions differ")
        output_values = prediction.get("values")
        if (prediction.get("status") in {"success", "complete", "partial_nonfinite"}
                and (not isinstance(output_values, dict) or set(output_values) != {entry["id"] for entry in endpoints})):
            raise ValueError("Cached outputs do not match the defined endpoints")
        return values

    def for_row(self, row, *, key="admet"):
        cached = row.get(key)
        records = [cached] if isinstance(cached, dict) else cached
        if not isinstance(records, list):
            return deepcopy(cached)
        result = []
        for index, item in enumerate(records):
            if not isinstance(item, dict):
                result.append(deepcopy(item))
                continue
            prediction = deepcopy(item)
            # A display response never inherits an unchecked inline model card.
            for field in _FIELDS:
                prediction.pop(field, None)
            if index >= 100:
                prediction["display_metadata_status"] = "unavailable_or_unverified"
                result.append(prediction)
                continue
            try:
                prediction.update(deepcopy(self._prediction(row, item)))
                prediction["display_metadata_status"] = "verified_source_cas"
            except (ValueError, TypeError, KeyError, OSError, AttributeError, RecursionError):
                prediction["display_metadata_status"] = "unavailable_or_unverified"
            result.append(prediction)
        return result[0] if isinstance(cached, dict) and result else result
