'Validate the selected structure and retain explicit calculation provenance.'

from __future__ import annotations

from contextlib import contextmanager
from collections import OrderedDict, defaultdict
from copy import deepcopy
import fcntl
import hashlib
import json
import math
from pathlib import Path, PurePosixPath
import re
from statistics import median
from threading import Lock
import time
from typing import Any

from .artifacts import ArtifactIntegrityError, atomic_bytes, canonical_json
from .engines import inspect_output, validate_request, validate_engine_hashes
from .identity import chemical_identity, candidate_identity, digest
from .confirmation import confirmation_policy, select_confirmation
from .ranking import (
    COMPARABILITY_VERSION,
    PRIORITY_POLICY_VERSION,
    _percentiles,
    comparable_protocol,
    finite,
    rank_candidates,
)
from .procurement import has_purchase_record
from .algorithm_profiles import validate_profile

SCHEMA = "skinpipe-campaign-v2.1"
LEGACY_SETTINGS = {
    "engine": "hybrid",
    "seed": 42,
    "n_poses": 9,
    "exhaustiveness": 16,
    "search_mode": "balance",
    "batch_size": 300,
}
SUCCESS = {"complete"}
TERMINAL = {"complete", "failed", "cancelled", "needs_data", "budget_deferred"}
_STATUS_PROJECTIONS = OrderedDict()
_STATUS_PROJECTION_LOCK = Lock()
_STATUS_PROJECTION_BYTES = 0
_STATUS_PROJECTION_MAX_ENTRIES = 32768
_STATUS_PROJECTION_MAX_BYTES = 16 * 1024 * 1024
_STATUS_PROJECTION_MAX_ENTRY_BYTES = 128 * 1024
_MANIFEST_VALIDATIONS = OrderedDict()
_MANIFEST_VALIDATION_LOCK = Lock()
_MANIFEST_VALIDATION_MAX_ENTRIES = 256


def _sha(raw):
    return hashlib.sha256(raw).hexdigest()


def _identifier(value):
    if not isinstance(value, str) or not re.fullmatch(r"[a-f0-9]{64}", value):
        raise ValueError("campaign identifier must be a SHA256 digest")
    return value


def _sealed_write(path, payload):
    atomic_bytes(
        path,
        canonical_json({"payload": payload, "sha256": _sha(canonical_json(payload))}),
    )


def _sealed_payload(raw):
    try:
        envelope = json.loads(raw)
        if _sha(canonical_json(envelope["payload"])) != envelope["sha256"]:
            raise ArtifactIntegrityError("campaign checkpoint checksum mismatch")
        return envelope["payload"]
    except (KeyError, TypeError, json.JSONDecodeError) as error:
        raise ArtifactIntegrityError("campaign checkpoint unreadable") from error


def _sealed_read(path):
    try:
        return _sealed_payload(Path(path).read_bytes())
    except OSError as error:
        raise ArtifactIntegrityError("campaign checkpoint unreadable") from error


def _verified_manifest(path, expected_id):
    'Read fresh bytes, retaining only bounded proofs of seal/identity checks.'
    path = Path(path)
    _identifier(expected_id)
    try:
        raw = path.read_bytes()
    except OSError as error:
        raise ArtifactIntegrityError("campaign checkpoint unreadable") from error
    return _manifest_payload(raw, path, expected_id)


def _manifest_payload(raw, path, expected_id):
    """Verify bytes just read by the caller; cache no population or file handle."""
    key = (_sha(str(path.absolute()).encode()), _sha(raw), expected_id)
    with _MANIFEST_VALIDATION_LOCK:
        known = key in _MANIFEST_VALIDATIONS
        if known:
            _MANIFEST_VALIDATIONS.move_to_end(key)
    if known:
        return json.loads(raw)["payload"]
    manifest = _sealed_payload(raw)
    if digest(manifest) != expected_id:
        raise ArtifactIntegrityError("campaign manifest identity mismatch")
    with _MANIFEST_VALIDATION_LOCK:
        _MANIFEST_VALIDATIONS[key] = None
        _MANIFEST_VALIDATIONS.move_to_end(key)
        while len(_MANIFEST_VALIDATIONS) > _MANIFEST_VALIDATION_MAX_ENTRIES:
            _MANIFEST_VALIDATIONS.popitem(last=False)
    return manifest


def _status_projection(path, expected_id, kind):
    """Cache small facts only after current full bytes and semantic ID validation."""
    global _STATUS_PROJECTION_BYTES
    path = Path(path)
    try:
        raw = path.read_bytes()
    except OSError as error:
        raise ArtifactIntegrityError("campaign checkpoint unreadable") from error
    # Path namespaces owners; neither raw bytes nor complete inputs are retained.
    key = (str(path.absolute()), _sha(raw), expected_id, kind)
    with _STATUS_PROJECTION_LOCK:
        encoded = _STATUS_PROJECTIONS.get(key)
        if encoded is not None:
            _STATUS_PROJECTIONS.move_to_end(key)
    if encoded is not None:
        return json.loads(encoded)
    if kind == "manifest":
        payload = _manifest_payload(raw, path, expected_id)
    else:
        payload = _sealed_payload(raw)
        if digest(payload) != expected_id:
            raise ArtifactIntegrityError("candidate receipt identity mismatch")
    if kind == "manifest":
        plan = payload["plan"]
        projection = {
            "mode": payload["mode"],
            "population_count": len(payload["candidates"]),
            "held": payload["held"],
            "plan": {
                key: plan[key]
                for key in (
                    "snapshot_id",
                    "target_id",
                    "input_id",
                    "input_validation",
                    "target_population_count",
                    "population_count",
                    "plan_id",
                    "budget_seconds",
                )
                if key in plan
            },
        }
    else:
        row = payload["row"]
        projection = {
            "candidate_id": row["candidate_id"],
            "status": row["status"],
            "input_sha256": row["input_sha256"],
            "protocol": {
                "seed": row["protocol"].get("seed"),
                "receptor_sha256": row["protocol"]["receptor_sha256"],
            },
            "pose_sha256": row.get("pose_sha256"),
            "sdf_artifact": (
                {"sha256": row["sdf_artifact"]["sha256"]}
                if row.get("sdf_artifact")
                else None
            ),
        }
    encoded = canonical_json(projection)
    if (
        len(encoded) > _STATUS_PROJECTION_MAX_ENTRY_BYTES
        or len(encoded) > _STATUS_PROJECTION_MAX_BYTES
    ):
        return projection
    with _STATUS_PROJECTION_LOCK:
        previous = _STATUS_PROJECTIONS.pop(key, None)
        if previous is not None:
            _STATUS_PROJECTION_BYTES -= len(previous)
        _STATUS_PROJECTIONS[key] = encoded
        _STATUS_PROJECTION_BYTES += len(encoded)
        while (
            len(_STATUS_PROJECTIONS) > _STATUS_PROJECTION_MAX_ENTRIES
            or _STATUS_PROJECTION_BYTES > _STATUS_PROJECTION_MAX_BYTES
        ):
            _, evicted = _STATUS_PROJECTIONS.popitem(last=False)
            _STATUS_PROJECTION_BYTES -= len(evicted)
    return json.loads(encoded)


def _safe_relative(value):
    path = PurePosixPath(value)
    if path.is_absolute() or ".." in path.parts:
        raise ValueError("legacy manifest path escapes campaign")
    return path.as_posix()


def import_legacy_pool(
    manifest_bytes: bytes,
    expected_sha256: str,
    pool_key: str,
    *,
    read_file,
    preparation_receipt_paths=(),
) -> dict[str, Any]:
    'Read exact phase98/101 input bytes through a caller-scoped read function.'
    if _sha(manifest_bytes) != expected_sha256:
        raise ArtifactIntegrityError("legacy manifest does not match pinned checksum")
    manifest = json.loads(manifest_bytes)
    compute = manifest.get("compute", {})
    if any(compute.get(key) != value for key, value in LEGACY_SETTINGS.items()):
        raise ValueError("legacy scientific settings differ from phase98/101 contract")
    pools = [pool for pool in manifest.get("pools", []) if pool.get("key") == pool_key]
    if len(pools) != 1:
        raise ValueError("legacy pool is absent or ambiguous")
    pool = pools[0]
    selected_path = _safe_relative("pools/" + pool_key + "/selected.json")
    selected_bytes = read_file(selected_path)
    if _sha(selected_bytes) != pool["selected_sha256"]:
        raise ArtifactIntegrityError("legacy selected candidates changed")
    selected = json.loads(selected_bytes)
    if len(selected) != pool["n_selected"] or len(
        {r["ligand_id"] for r in selected}
    ) != len(selected):
        raise ValueError("legacy selected count or candidate uniqueness mismatch")
    target_id = pool["node"]
    candidates = [
        {
            **row,
            "candidate_id": row["ligand_id"],
            "target_id": target_id,
            "legacy_selection_index": i,
        }
        for i, row in enumerate(selected)
    ]
    by_id = {row["candidate_id"]: row for row in candidates}
    prepared, artifacts, rejected, ambiguous = {}, {}, [], set()
    for receipt_path in sorted(preparation_receipt_paths):
        try:
            receipt_bytes = read_file(str(receipt_path))
            receipt = json.loads(receipt_bytes)
            cid = receipt.get("ligand_id")
            if cid not in by_id:
                continue
            if cid in ambiguous:
                raise ArtifactIntegrityError("ambiguous prepared states remain held")
            request = receipt["request"]
            if (
                receipt.get("state") != "completed"
                or receipt.get("reusable") is not True
                or receipt.get("software_snapshot_stable") is not True
                or request.get("smiles") != by_id[cid]["smiles"]
                or receipt.get("input_smiles") != by_id[cid]["smiles"]
                or receipt.get("request_sha256") != digest(request)
            ):
                raise ArtifactIntegrityError(
                    "legacy preparation does not bind frozen source identity"
                )
            name = receipt["pdbqt_file"]
            if PurePosixPath(name).name != name:
                raise ValueError("unsafe preparation filename")
            pdbqt_path = str(PurePosixPath(receipt_path).parent / name)
            pdbqt_bytes = read_file(pdbqt_path)
            if _sha(pdbqt_bytes) != receipt["pdbqt_sha256"]:
                raise ArtifactIntegrityError("legacy prepared ligand checksum mismatch")
            item = {
                "pdbqt_sha256": _sha(pdbqt_bytes),
                "preparation_receipt_sha256": _sha(receipt_bytes),
                "exact_structure_id": chemical_identity(by_id[cid]["smiles"])[
                    "exact_structure_id"
                ],
                "prepared_structure_id": "legacy-prepared:" + receipt["request_sha256"],
                "source_path": pdbqt_path,
                "source_receipt_path": str(receipt_path),
            }
            if (
                cid in prepared
                and prepared[cid]["pdbqt_sha256"] != item["pdbqt_sha256"]
            ):
                ambiguous.add(cid)
                prepared.pop(cid)
                raise ArtifactIntegrityError(
                    "multiple prepared states for one frozen candidate"
                )
            prepared[cid] = item
            artifacts[_sha(pdbqt_bytes)] = pdbqt_bytes
            artifacts[_sha(receipt_bytes)] = receipt_bytes
        except (ValueError, KeyError, OSError, TypeError) as error:
            rejected.append({"receipt_path": str(receipt_path), "reason": str(error)})
    receptor = None
    if pool.get("receptor"):
        frozen = pool["receptor"]
        card_path = _safe_relative(frozen["card"])
        card_bytes = read_file(card_path)
        if _sha(card_bytes) != frozen["card_sha256"]:
            raise ArtifactIntegrityError("legacy receptor identity card changed")
        card = json.loads(card_bytes)
        pdbqt_path = _safe_relative(
            str(PurePosixPath(card_path).parent / card["pdbqt"])
        )
        pdbqt_bytes = read_file(pdbqt_path)
        if _sha(pdbqt_bytes) != frozen["pdbqt_sha256"]:
            raise ArtifactIntegrityError("legacy receptor bytes changed")
        accession = frozen.get("target_model", {}).get("accession")
        if not accession:
            raise ValueError("legacy receptor lacks exact target accession")
        receptor = {
            "receptor_sha256": _sha(pdbqt_bytes),
            "center": card["receptor"]["center"],
            "size": card["receptor"]["size"],
            "actual_protein_id": "uniprot:" + accession,
            "target_identity_valid": True,
            "card_sha256": _sha(card_bytes),
            "functional_admissibility": frozen.get(
                "functional_admissibility", "NOT_ESTABLISHED"
            ),
        }
        artifacts[_sha(pdbqt_bytes)] = pdbqt_bytes
        artifacts[_sha(card_bytes)] = card_bytes
    for row in candidates:
        if row["candidate_id"] not in prepared:
            prepared[row["candidate_id"]] = {
                "status": "held",
                "reason": "verified_prepared_artifact_not_supplied",
            }
    artifacts[expected_sha256] = manifest_bytes
    artifacts[_sha(selected_bytes)] = selected_bytes
    return {
        "candidates": candidates,
        "prepared": prepared,
        "receptor": receptor,
        "artifacts": artifacts,
        "source": {
            "manifest_sha256": expected_sha256,
            "selected_sha256": _sha(selected_bytes),
            "source_commit": manifest.get("source_commit"),
            "pool_key": pool_key,
            "settings": deepcopy(compute),
            "selection_policy": manifest.get("selection_policy"),
            "original_state": pool.get("state"),
            "original_reason": pool.get("reason"),
            "parent_manifest_sha256": manifest.get("parent_manifest_sha256"),
            "unidock_binary_sha256": manifest.get("unidock_binary_sha256"),
            "batch_semantics": "legacy native ordered batches; results depend on co-batched identities",
            "ownership_transferred": False,
            "execution_compatible_with_adaptive_v2": False,
        },
        "rejected_preparations": rejected,
        "status": "imported_inputs_only_not_submitted",
    }


class CampaignCoordinator:
    """One API-host coordinator, restartable across submissions and publications."""

    def __init__(self, root, artifacts, queue=None):
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)
        self.artifacts = artifacts
        self.queue = queue

    def _directory(self, campaign_id):
        return self.root / _identifier(campaign_id)

    @contextmanager
    def _lock(self, campaign_id):
        directory = self._directory(campaign_id)
        if not directory.is_dir():
            raise KeyError("campaign not found")
        with (directory / ".lock").open("a+b") as handle:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
            try:
                yield directory
            finally:
                fcntl.flock(handle.fileno(), fcntl.LOCK_UN)

    def _manifest(self, campaign_id):
        return _verified_manifest(
            self._directory(campaign_id) / "manifest.json", campaign_id
        )

    def create(
        self,
        plan,
        candidates,
        prepared,
        receptor,
        *,
        mode="upgraded",
        legacy_source=None,
        trust_stored_identity=False,
        preparation_holds=None,
    ):
        if "engine_hashes" in plan:
            validate_engine_hashes(plan["engine_hashes"])
        if mode not in {"upgraded", "legacy_campaign"}:
            raise ValueError("unknown campaign mode")
        if mode == "legacy_campaign" and not legacy_source:
            raise ValueError("legacy mode requires pinned legacy source metadata")
        rows = [
            deepcopy(r) for r in candidates if r.get("target_id") == plan["target_id"]
        ]
        rows.sort(key=lambda r: r["candidate_id"])
        if len({r["candidate_id"] for r in rows}) != len(rows):
            raise ValueError("candidate IDs must be unique in campaign")
        mapping, held = {}, []
        for row in rows:
            cid = row["candidate_id"]
            identity = candidate_identity(row, trust_stored=trust_stored_identity)
            supplied = deepcopy(prepared.get(cid, {}))
            row["identity"] = identity
            mapping[cid] = supplied
            if supplied.get("pdbqt_sha256"):
                if (
                    supplied.get("exact_structure_id") != identity["exact_structure_id"]
                    or identity["exact_structure_id"] is None
                ):
                    raise ValueError(
                        "prepared artifact must bind the candidate exact_structure_id"
                    )
                self.artifacts.verified_path(supplied["pdbqt_sha256"])
                if supplied.get("preparation_receipt_sha256"):
                    self.artifacts.verified_path(supplied["preparation_receipt_sha256"])
            else:
                diagnostic = deepcopy((preparation_holds or {}).get(cid))
                if diagnostic:
                    if diagnostic.get("original_structure_id") != identity["original_structure_id"]:
                        raise ArtifactIntegrityError("Preparation diagnostic source structure changed")
                    if diagnostic.get("diagnostic_log_sha256"):
                        self.artifacts.verified_path(diagnostic["diagnostic_log_sha256"])
                    row["preparation_diagnostic"] = diagnostic
                    row["computation_status"] = "preparation_held"
                held.append({"candidate_id": cid, "reason": diagnostic["reason"] if diagnostic else supplied.get("reason", "preparation_required"),
                             **({"preparation_diagnostic": diagnostic} if diagnostic else {})})
        if receptor:
            context = plan.get("target_context") or next(
                (
                    job["request"].get("target_context", {})
                    for job in plan.get("jobs", [])
                    if job["request"].get("stage") == "broad_docking"
                ),
                {},
            )
            self._validate_receptor(receptor, context.get("actual_protein_id"))
        if set(plan.get("population_ids", [])) - {r["candidate_id"] for r in rows}:
            raise ValueError("plan candidates missing from frozen campaign population")
        broad_ids = (
            plan.get("stages", {}).get("broad_docking", {}).get("candidate_ids", [])
        )
        if mode == "legacy_campaign":
            broad_ids = [
                r["candidate_id"]
                for r in sorted(rows, key=lambda r: r.get("legacy_selection_index", 0))
            ]
        if (
            len(broad_ids) > 1000
            or len(set(broad_ids)) != len(broad_ids)
            or not set(broad_ids) <= set(mapping)
        ):
            raise ValueError(
                "broad stage requires at most1000 distinct frozen candidate IDs"
            )
        manifest = {
            "schema": SCHEMA,
            "plan": deepcopy(plan),
            "candidates": rows,
            "prepared": mapping,
            "receptor": deepcopy(receptor),
            "mode": mode,
            "legacy_source": deepcopy(legacy_source),
            "held": held,
            "protocol_version": "adaptive-batching-v2.1"
            if mode == "upgraded"
            else "legacy-native-frozen",
            "claim": "frozen inputs; worker execution and biological validation not established",
        }
        cid = digest(manifest)
        directory = self._directory(cid)
        directory.mkdir(exist_ok=True)
        with self._lock(cid):
            if (directory / "manifest.json").exists():
                self._manifest(cid)
                if (directory / "state.json").exists():
                    return self.status(cid, _locked=True)
            _sealed_write(directory / "manifest.json", manifest)
            state = {"submissions": {}, "outcomes": {}, "stages": {}, "quarantined": []}
            if receptor:
                control_ids = (
                    [
                        r["candidate_id"]
                        for r in rows
                        if r.get("role")
                        in {
                            "control",
                            "positive_control",
                            "negative_control",
                            "counter_control",
                        }
                        and mapping[r["candidate_id"]].get("pdbqt_sha256")
                    ]
                    if (plan.get("compute_profile") or {}).get("version")
                    == "generalized-research-v1"
                    else []
                )
                stage = self._stage(
                    manifest,
                    "broad_docking",
                    broad_ids,
                    control_ids,
                    receptor,
                    [42],
                    "unidock",
                )
                _sealed_write(directory / "stage-broad_docking.json", stage)
                state["stages"]["broad_docking"] = digest(stage)
            _sealed_write(directory / "state.json", state)
        return self.status(cid)

    def _validate_receptor(self, receptor, actual_protein_id=None):
        if (
            receptor.get("target_identity_valid") is not True
            or not receptor.get("actual_protein_id")
            or actual_protein_id
            and receptor["actual_protein_id"] != actual_protein_id
        ):
            raise ValueError(
                "receptor actual protein identity is not verified or mismatches target"
            )
        self.artifacts.verified_path(receptor["receptor_sha256"])
        for key in ("center", "size"):
            values = receptor.get(key)
            if (
                not isinstance(values, list)
                or len(values) != 3
                or any(
                    isinstance(v, bool)
                    or not isinstance(v, (int, float))
                    or not math.isfinite(v)
                    for v in values
                )
            ):
                raise ValueError("receptor box requires finite triples")
        if any(not 0 < v <= 100 for v in receptor["size"]):
            raise ValueError("receptor box size out of range")

    @staticmethod
    def _stage(manifest, stage, candidate_ids, control_ids, receptor, seeds, engine):
        return {
            **({"engine_hashes": deepcopy(manifest["plan"]["engine_hashes"])}
               if "engine_hashes" in manifest["plan"] else {}),
            "stage": stage,
            "candidate_ids": list(candidate_ids),
            "control_ids": list(control_ids),
            "receptor": deepcopy(receptor),
            "seeds": list(seeds),
            "engine": engine,
            "settings": {"n_poses": 9, "exhaustiveness": 16, "search_mode": "balance"},
            "snapshot_id": manifest["plan"]["snapshot_id"],
            "policy_id": manifest["plan"]["policy_id"],
            "reference_sha256": manifest["plan"]["reference_snapshot"]["sha256"],
            "protocol_version": manifest["protocol_version"],
        }

    def _stage_read(self, directory, state, name):
        if name not in state["stages"]:
            raise ValueError("stage not frozen; call advance after upstream receipts")
        stage = _sealed_read(directory / ("stage-" + name + ".json"))
        if digest(stage) != state["stages"][name]:
            raise ArtifactIntegrityError("stage checkpoint hash mismatch")
        return stage

    @staticmethod
    def _stage_members(state, name, spec):
        scope = state.get("confirmation_scope") if name == "confirmation" else None
        if not scope:
            return spec["candidate_ids"], spec["control_ids"], spec["seeds"]
        if (
            scope["source_stage_sha256"] != digest(spec)
            or not set(scope["candidate_ids"]) <= set(spec["candidate_ids"])
            or scope["control_ids"] != spec["control_ids"]
            or not set(scope["seeds"]) <= set(spec["seeds"])
        ):
            raise ArtifactIntegrityError(
                "Confirmation execution scope differs from its frozen source stage"
            )
        return scope["candidate_ids"], scope["control_ids"], scope["seeds"]

    def amend_confirmation(self, campaign_id, policy, *, reason, owner="researcher"):
        """Record one explicit execution amendment; original manifests/stages stay sealed.

        Queued chunks are retired only once a verified score selection is frozen
        by advance. Claimed chunks finish normally and their receipts are kept.
        """
        config = confirmation_policy(policy)
        if not isinstance(reason, str) or not reason.strip() or len(reason) > 2000:
            raise ValueError(
                "Confirmation amendment requires a bounded nonempty reason"
            )
        manifest = self._manifest(campaign_id)
        if manifest["mode"] != "upgraded":
            raise ValueError("Legacy scientific settings cannot be amended")
        with self._lock(campaign_id) as directory:
            state = _sealed_read(directory / "state.json")
            if "confirmation" in state["stages"]:
                spec = self._stage_read(directory, state, "confirmation")
                if not set(config["seeds"]) <= set(spec["seeds"]):
                    raise ValueError(
                        "Amendment seeds must be a subset of the frozen confirmation stage"
                    )



                self._amended_confirmation_ids(manifest, directory, state, spec, config)
            previous = state.get("confirmation_policy_override")
            if previous:
                if previous != config:
                    raise ValueError(
                        "A different confirmation amendment requires a new campaign"
                    )
                return deepcopy(state["confirmation_amendment"])
            amendment = {
                "schema": "skinpipe-confirmation-amendment-v1",
                "campaign_id": campaign_id,
                "policy": config,
                "source_plan_id": manifest["plan"]["plan_id"],
                "reason": reason.strip(),
                "owner": owner,
                "recorded_at": time.time(),
                "preserve_running_jobs": True,
                "original_inputs_unchanged": True,
            }
            amendment["amendment_id"] = digest(amendment)
            amendment_path = directory / ('confirmation-amendment-' + amendment['amendment_id'] + '.json')
            state_path = directory / 'state.json'
            original_state = state_path.read_bytes()
            if amendment_path.exists():
                raise ArtifactIntegrityError('Unactivated confirmation amendment already exists; review before retry')
            state["confirmation_policy_override"] = config
            state["confirmation_amendment"] = amendment
            try:
                _sealed_write(amendment_path, amendment)
                _sealed_write(state_path, state)
            except Exception:
                # An ordinary write failure is not permission to leave a new
                # irreversible override partially activated. Never delete a
                # prior amendment; the path above belongs only to this attempt.
                if state_path.read_bytes() != original_state:
                    atomic_bytes(state_path, original_state)
                amendment_path.unlink(missing_ok=True)
                raise
            return deepcopy(amendment)

    def _amended_confirmation_ids(self, manifest, directory, state, spec, config):
        broad_rows, missing = self._stage_rows(directory, state, 'broad_docking')
        if any(manifest['prepared'].get(cid, {}).get('pdbqt_sha256') for cid in missing):
            raise ValueError('Confirmation amendment requires intact frozen broad receipts')
        if 'broad_cpu_fallback' in state['stages']:
            fallback_rows, missing = self._stage_rows(directory, state, 'broad_cpu_fallback')
            if missing:
                raise ValueError('Confirmation amendment requires intact frozen fallback receipts')
            by_candidate = {row['candidate_id']: row for row in broad_rows}
            by_candidate.update({row['candidate_id']: row for row in fallback_rows})
            broad_rows = list(by_candidate.values())
        ids = select_confirmation(self._rank_results(manifest, broad_rows, unique_parents=False), broad_rows, config)
        if not set(ids) <= set(spec['candidate_ids']):
            raise ValueError('Score selection extends the frozen confirmation pool; create a new campaign')
        return ids

    def _retire_confirmation_chunks(self, campaign_id, state, spec):
        candidates, controls, seeds = self._stage_members(state, "confirmation", spec)
        allowed = set(candidates + controls)
        for record in state["submissions"].values():
            request = record["request"]
            if record.get("collected") or record.get("retired_by_confirmation_policy"):
                continue
            if request["stage"] != "confirmation" or not record.get("job_id"):
                continue
            if request["seed"] in seeds and set(record["candidate_ids"]) <= allowed:
                continue
            result = self.queue.cancel_if_queued(
                record["job_id"],
                owner=record["owner"],
                campaign_id=campaign_id,
                stage="confirmation",
                reason="superseded_confirmation_policy",
            )
            if (
                result["status"] == "cancelled"
                and result.get("error") == "superseded_confirmation_policy"
            ):
                record["retired_by_confirmation_policy"] = state[
                    "confirmation_amendment"
                ]["amendment_id"]

    @staticmethod
    def _task_key(stage, seed):
        return digest({"stage": stage, "seed": seed})

    def _read_outcome(self, directory, item, *, verified_assets=None):
        receipt = _sealed_read(
            directory / "receipts" / (item["receipt_sha256"] + ".json")
        )
        if digest(receipt) != item["receipt_sha256"]:
            raise ArtifactIntegrityError("candidate receipt identity mismatch")
        row = receipt["row"]
        self._verify_outcome_assets(row, verified_assets)
        return receipt

    def _verify_outcome_assets(self, row, verified_assets):
        def verify(sha256):
            key = (self.artifacts, sha256)
            if verified_assets is None or key not in verified_assets:
                self.artifacts.verified_path(sha256)
                if verified_assets is not None:
                    verified_assets.add(key)

        verify(row["input_sha256"])
        verify(row["protocol"]["receptor_sha256"])
        for key in ("pose_sha256",):
            if row.get(key):
                verify(row[key])
        if row.get("sdf_artifact"):
            verify(row["sdf_artifact"]["sha256"])

    def _read_status_outcome(self, directory, item, verified_assets):
        row = _status_projection(
            directory / "receipts" / (item["receipt_sha256"] + ".json"),
            item["receipt_sha256"],
            "receipt",
        )
        self._verify_outcome_assets(row, verified_assets)
        return row

    def dispatch(
        self,
        campaign_id,
        *,
        stage="broad_docking",
        chunk_size=32,
        owner="researcher",
        retry_failed=False,
    ):
        if self.queue is None:
            raise ValueError("queue connection required for dispatch")
        if (
            isinstance(chunk_size, bool)
            or not isinstance(chunk_size, int)
            or not 1 <= chunk_size <= 1000
        ):
            raise ValueError("chunk size must be 1..1000")
        manifest = self._manifest(campaign_id)
        profile = manifest["plan"].get("policy", {}).get("compute_profile")
        if profile:
            validate_profile(profile)
        if manifest["mode"] == "legacy_campaign":
            raise ValueError(
                "legacy handoff held: adaptive batches do not preserve native frozen batch semantics; ownership and compatible native engine required"
            )
        with self._lock(campaign_id) as directory:
            state = _sealed_read(directory / "state.json")
            spec = self._stage_read(directory, state, stage)
            self._validate_receptor(
                spec["receptor"], manifest["receptor"]["actual_protein_id"]
            )
            submissions = state["submissions"]
            # Recover an ambiguous enqueue using its original idempotency key.
            for record in submissions.values():
                if record.get("job_id") is None:
                    job = self.queue.enqueue(
                        "dock_batch",
                        record["request"],
                        owner=record["owner"],
                        idempotency_key=record["idempotency_key"],
                        budget_seconds=manifest["plan"]["budget_seconds"],
                    )
                    record["job_id"] = job["id"]
                    _sealed_write(directory / "state.json", state)
            created = []
            candidate_ids, control_ids, seeds = self._stage_members(state, stage, spec)
            if profile:
                key = (
                    "receptor_sensitivity"
                    if stage.startswith("receptor_sensitivity_")
                    else stage
                )
                if (
                    key in profile["counts"]
                    and len(candidate_ids) > profile["counts"][key]
                ):
                    raise ValueError(
                        "Frozen stage exceeds its compute profile candidate cap"
                    )
                if (
                    stage == "confirmation"
                    and seeds != profile["multiplicities"]["confirmation_seeds"]
                ):
                    raise ValueError(
                        "Frozen confirmation exceeds its seed workload contract"
                    )
            all_ids = candidate_ids + control_ids
            for seed in seeds:
                task_key = self._task_key(spec, seed)
                outcomes = state["outcomes"].setdefault(task_key, {})
                active_ids = set()
                generations = {}
                latest = {}
                for record in submissions.values():
                    if record["task_key"] != task_key:
                        continue
                    job = self.queue.get(record["job_id"], owner=record["owner"])
                    for cid in record["candidate_ids"]:
                        generations[cid] = max(
                            generations.get(cid, 0), record["generation"]
                        )
                        if cid not in latest or record["generation"] >= latest[cid][0]:
                            latest[cid] = (
                                record["generation"],
                                job["status"],
                                bool(record.get("retired_by_confirmation_policy")),
                            )
                    # Completed jobs must be collected before any resubmission.
                    if job["status"] not in TERMINAL or (
                        job["status"] == "complete" and not record.get("collected")
                    ):
                        active_ids.update(record["candidate_ids"])
                pending = []
                for cid in all_ids:
                    if cid in active_ids or not manifest["prepared"].get(cid, {}).get(
                        "pdbqt_sha256"
                    ):
                        continue
                    try:
                        self.artifacts.verified_path(
                            manifest["prepared"][cid]["pdbqt_sha256"]
                        )
                    except ArtifactIntegrityError as error:
                        state["quarantined"].append(
                            {
                                "candidate_id": cid,
                                "reason": str(error),
                                "at": time.time(),
                            }
                        )
                        continue
                    outcome = outcomes.get(cid)
                    if outcome:
                        try:
                            row = self._read_outcome(directory, outcome)["row"]
                            if row["status"] in SUCCESS or not retry_failed:
                                continue
                        except ArtifactIntegrityError as error:
                            state["quarantined"].append(
                                {
                                    "candidate_id": cid,
                                    "receipt": outcome,
                                    "reason": str(error),
                                    "at": time.time(),
                                }
                            )
                            outcomes.pop(cid)
                    elif (
                        generations.get(cid)
                        and not retry_failed
                        and not (
                            latest.get(cid, (0, None, False))[1:] == ("cancelled", True)
                        )
                    ):
                        continue
                    pending.append(cid)
                for offset in range(0, len(pending), chunk_size):
                    ids = pending[offset : offset + chunk_size]
                    generation = 1 + max(
                        (generations.get(cid, 0) for cid in ids), default=0
                    )
                    receptor = spec["receptor"]
                    request = {
                        **({"engine_hashes": deepcopy(spec["engine_hashes"])}
                           if "engine_hashes" in spec else {}),
                        **spec["settings"],
                        "engine": spec["engine"],
                        "seed": seed,
                        "receptor_sha256": receptor["receptor_sha256"],
                        "center": receptor["center"],
                        "size": receptor["size"],
                        "ligands": [
                            {
                                "candidate_id": cid,
                                "pdbqt_sha256": manifest["prepared"][cid][
                                    "pdbqt_sha256"
                                ],
                            }
                            for cid in ids
                        ],
                        "campaign_id": campaign_id,
                        "stage": stage,
                        "task_key": task_key,
                        "protocol_version": spec["protocol_version"],
                        "target_id": manifest["plan"]["target_id"],
                        "actual_protein_id": receptor["actual_protein_id"],
                        "snapshot_id": spec["snapshot_id"],
                        "policy_id": spec["policy_id"],
                        "reference_sha256": spec["reference_sha256"],
                        "budget_scope": manifest["plan"]["budget_scope"],
                        "target_budget_seconds": manifest["plan"]["budget_seconds"],
                        "batch_membership_sha256": digest(ids),
                        "canary": False,
                    }
                    validate_request(request)
                    key = digest({"request": request, "generation": generation})
                    record = {
                        "request": request,
                        "idempotency_key": key,
                        "task_key": task_key,
                        "candidate_ids": ids,
                        "generation": generation,
                        "owner": owner,
                        "job_id": None,
                        "collected": False,
                    }
                    submissions[key] = record
                    _sealed_write(directory / "state.json", state)
                    job = self.queue.enqueue(
                        "dock_batch",
                        request,
                        owner=owner,
                        idempotency_key=key,
                        budget_seconds=manifest["plan"]["budget_seconds"],
                    )
                    record["job_id"] = job["id"]
                    created.append(job)
                    _sealed_write(directory / "state.json", state)
            _sealed_write(directory / "state.json", state)
            return {
                "campaign_id": campaign_id,
                "stage": stage,
                "jobs": created,
                "submitted_count": len(created),
                "execution_completed": False,
            }

    def collect(self, campaign_id):
        manifest = self._manifest(campaign_id)
        if self.queue is None:
            raise ValueError("queue connection required")
        accepted, quarantined = 0, 0
        with self._lock(campaign_id) as directory:
            state = _sealed_read(directory / "state.json")
            for submission in state["submissions"].values():
                if submission.get("collected") or not submission.get("job_id"):
                    continue
                job = self.queue.get(submission["job_id"], owner=submission["owner"])
                request = submission["request"]
                complete_job = job["status"] == "complete"
                result = (
                    (job.get("result") or {})
                    if complete_job
                    else (job.get("progress") or {}).get("partial_result", {})
                )
                try:
                    if job.get("request") is not None and digest(
                        job["request"]
                    ) != digest(request):
                        raise ArtifactIntegrityError(
                            "queue job request changed from frozen submission"
                        )
                    if not result and not complete_job:
                        continue
                    if not isinstance(result, dict):
                        raise ArtifactIntegrityError("worker result must be an object")
                    rows = result.get("candidates", [])
                    if not isinstance(rows, list) or any(
                        not isinstance(r, dict)
                        or not isinstance(r.get("candidate_id"), str)
                        for r in rows
                    ):
                        raise ArtifactIntegrityError(
                            "worker candidates must be identified object rows"
                        )
                    reported_ids = {r["candidate_id"] for r in rows}
                    expected_ids = set(submission["candidate_ids"])
                    if (
                        len(reported_ids) != len(rows)
                        or not reported_ids <= expected_ids
                        or complete_job
                        and reported_ids != expected_ids
                    ):
                        raise ArtifactIntegrityError(
                            "worker receipt does not account for exact submitted candidates"
                        )
                except (ValueError, TypeError) as error:
                    quarantine_id = digest(
                        {"job_id": job["id"], "result": result, "reason": str(error)}
                    )
                    recorded = submission.setdefault("quarantine_ids", [])
                    if quarantine_id not in recorded:
                        state["quarantined"].append(
                            {
                                "quarantine_id": quarantine_id,
                                "job_id": job["id"],
                                "candidate_ids": submission["candidate_ids"],
                                "result": result,
                                "reason": str(error),
                                "at": time.time(),
                            }
                        )
                        recorded.append(quarantine_id)
                        quarantined += len(submission["candidate_ids"])
                    submission["collected"] = complete_job
                    submission["receipt_quarantined"] = True
                    _sealed_write(directory / "state.json", state)
                    continue
                outcomes = state["outcomes"].setdefault(submission["task_key"], {})
                for row in rows:
                    cid = row["candidate_id"]
                    try:
                        self._verify_worker_row(
                            row, request, manifest["prepared"][cid], result
                        )
                        row = deepcopy(row)
                        row.setdefault("protocol", result.get("protocol"))
                        receipt = {
                            "campaign_id": campaign_id,
                            "task_key": submission["task_key"],
                            "job_id": job["id"],
                            "candidate_id": cid,
                            "request_sha256": digest(request),
                            "row": row,
                            "runtime": result.get("runtime"),
                            "worker_telemetry": result.get("worker_telemetry"),
                            "scientific_validation": "NOT_ESTABLISHED",
                        }
                        receipt_id = digest(receipt)
                        # A completed verified candidate cannot be replaced by a later failure.
                        # Repeated partial failures are observations too: counting the
                        # same immutable receipt again would inflate collection totals.
                        previous = outcomes.get(cid)
                        if previous:
                            try:
                                previous_receipt = self._read_outcome(
                                    directory, previous
                                )
                                if (
                                    previous["receipt_sha256"] == receipt_id
                                    or previous_receipt["row"]["status"] in SUCCESS
                                ):
                                    continue
                            except ArtifactIntegrityError:
                                pass
                        _sealed_write(
                            directory / "receipts" / (receipt_id + ".json"), receipt
                        )
                        outcomes[cid] = {
                            "receipt_sha256": receipt_id,
                            "status": row["status"],
                        }
                        accepted += 1
                    except (ValueError, KeyError, TypeError, OSError) as error:
                        state["quarantined"].append(
                            {
                                "candidate_id": cid,
                                "job_id": job["id"],
                                "reason": str(error),
                                "row": row,
                                "at": time.time(),
                            }
                        )
                        quarantined += 1
                submission["collected"] = complete_job
                _sealed_write(directory / "state.json", state)
        return {
            "accepted_candidate_receipts": accepted,
            "quarantined_candidate_receipts": quarantined,
        }

    def reuse_broad_results(self, campaign_id, source_id, *, engine_hashes,
                           include_terminal_failures=False):
        """Reuse verified successes under exact preparation/receptor/engine settings.

        Selection weights may differ. Executing runtime and batch membership
        remain those of the original calculation, with an explicit source link.
        Caller supplies freshly verified engine binary hashes; no engine name
        alone authorizes reuse. Source receipts and failures remain unchanged.
        An explicit execution-recovery caller can also preserve terminal
        chemical/pose failures under the exact same executing binary/protocol;
        this avoids rerunning a known failure merely to rebuild a cohort.
        """
        if type(include_terminal_failures) is not bool:
            raise ValueError("Terminal failure reuse requires an explicit boolean")
        if campaign_id == source_id:
            return {"reused": 0, "held": 0}
        source_manifest = self._manifest(source_id)
        source_dir = self._directory(source_id)
        source_state = _sealed_read(source_dir / "state.json")
        source_spec = self._stage_read(source_dir, source_state, "broad_docking")
        reused = held = 0
        with self._lock(campaign_id) as directory:
            manifest = self._manifest(campaign_id)
            state = _sealed_read(directory / "state.json")
            spec = self._stage_read(directory, state, "broad_docking")
            if (
                manifest["plan"]["target_id"] != source_manifest["plan"]["target_id"]
                or manifest["receptor"]["actual_protein_id"]
                != source_manifest["receptor"]["actual_protein_id"]
                or spec["protocol_version"] != source_spec["protocol_version"]
            ):
                raise ValueError("cached result target identity differs")
            source_outcomes = source_state["outcomes"].get(
                self._task_key(source_spec, 42), {}
            )
            outcomes = state["outcomes"].setdefault(self._task_key(spec, 42), {})
            reusable = SUCCESS | ({"failed", "invalid_pose", "unsupported_chemistry"}
                                  if include_terminal_failures else set())
            for cid in spec["candidate_ids"] + spec["control_ids"]:
                if cid in outcomes or cid not in source_outcomes:
                    continue
                item = source_outcomes[cid]
                if item.get("status") not in reusable:
                    continue
                try:
                    old = self._read_outcome(source_dir, item)
                    row = old["row"]
                    engine = row["engine"]
                    if row["protocol"]["engine_identity"][
                        "sha256"
                    ] != engine_hashes.get(engine):
                        raise ValueError("cached engine binary differs")
                    prepared = manifest["prepared"][cid]
                    source_prepared = source_manifest["prepared"][cid]
                    if (
                        prepared.get("exact_structure_id")
                        != source_prepared.get("exact_structure_id")
                        or prepared.get("prepared_structure_id")
                        != source_prepared.get("prepared_structure_id")
                        or prepared.get("pdbqt_sha256")
                        != source_prepared.get("pdbqt_sha256")
                    ):
                        raise ValueError("cached exact preparation differs")
                    request = {
                        **({"engine_hashes": deepcopy(spec["engine_hashes"])}
                           if "engine_hashes" in spec else {}),
                        **spec["settings"],
                        "seed": 42,
                        "engine": spec["engine"],
                        "center": spec["receptor"]["center"],
                        "size": spec["receptor"]["size"],
                        "receptor_sha256": spec["receptor"]["receptor_sha256"],
                    }
                    self._verify_worker_row(row, request, prepared, old)
                    receipt = {
                        **deepcopy(old),
                        "campaign_id": campaign_id,
                        "task_key": self._task_key(spec, 42),
                        "reused_from": {
                            "campaign_id": source_id,
                            "receipt_sha256": item["receipt_sha256"],
                        },
                    }
                    receipt_id = digest(receipt)
                    _sealed_write(
                        directory / "receipts" / (receipt_id + ".json"), receipt
                    )
                    outcomes[cid] = {
                        "receipt_sha256": receipt_id,
                        "status": row["status"],
                    }
                    reused += 1
                except (ValueError, KeyError, TypeError, OSError):
                    held += 1
            _sealed_write(directory / "state.json", state)
        return {"reused": reused, "held": held}

    def _verify_worker_row(self, row, request, prepared, result):
        engine = row.get("engine")
        if engine not in {request["engine"], "vina"}:
            raise ValueError("receipt engine differs from requested engine/fallback")
        protocol = row.get("protocol") or result.get("protocol") or {}
        for key in (
            "seed",
            "n_poses",
            "exhaustiveness",
            "search_mode",
            "center",
            "size",
            "receptor_sha256",
        ):
            if protocol.get(key) != request[key]:
                raise ValueError("receipt protocol mismatch: " + key)
        if protocol.get("engine") != engine or not protocol.get(
            "engine_identity", {}
        ).get("sha256"):
            raise ValueError("receipt lacks matching executing engine identity")
        if "engine_hashes" in request:
            validate_engine_hashes(request["engine_hashes"])
            if protocol["engine_identity"]["sha256"] != request["engine_hashes"][engine]:
                raise ArtifactIntegrityError("Receipt executing binary differs from frozen cohort")
        if row.get("input_sha256") != prepared["pdbqt_sha256"]:
            raise ValueError("receipt prepared ligand identity mismatch")
        self.artifacts.verified_path(prepared["pdbqt_sha256"])
        self.artifacts.verified_path(request["receptor_sha256"])
        if row.get("pose_sha256"):
            raw = self.artifacts.get_bytes(row["pose_sha256"])
            checked = inspect_output(
                raw.decode(),
                self.artifacts.get_bytes(prepared["pdbqt_sha256"]).decode(),
            )
            if row.get("status") in {"complete", "review_required"} and any(
                row.get(key) != checked[key]
                for key in ("status", "affinity", "pose_scores", "n_poses")
            ):
                raise ValueError("reported score/status/pose metadata differs from actual pose bytes")
        elif row.get("status") in {"complete", "review_required"}:
            raise ValueError("successful receipt lacks pose artifact")
        if row.get("sdf_artifact"):
            self.artifacts.verified_path(row["sdf_artifact"]["sha256"])
        if row.get("status") not in {
            "complete",
            "review_required",
            "failed",
            "invalid_pose",
            "unsupported_chemistry",
        }:
            raise ValueError("unknown candidate computation state")

    def _stage_rows(
        self,
        directory,
        state,
        name,
        *,
        spec=None,
        verified_assets=None,
        status_only=False,
    ):
        if spec is None:
            spec = self._stage_read(directory, state, name)
        candidate_ids, control_ids, seeds = self._stage_members(state, name, spec)
        records, missing = [], []
        for seed in seeds:
            outcomes = state["outcomes"].get(self._task_key(spec, seed), {})
            for cid in candidate_ids + control_ids:
                if cid not in outcomes:
                    missing.append(cid)
                    continue
                try:
                    if status_only:
                        records.append(
                            self._read_status_outcome(
                                directory, outcomes[cid], verified_assets
                            )
                        )
                    else:
                        records.append(
                            self._read_outcome(
                                directory,
                                outcomes[cid],
                                verified_assets=verified_assets,
                            )["row"]
                        )
                except ArtifactIntegrityError:
                    missing.append(cid)
        return records, missing

    def advance(self, campaign_id):
        """Freeze next-stage choices only after all computable upstream rows account.

        Failed rows remain terminal evidence, not imputed scores. CPU/GPU and
        receptor/protocol cohorts stay separate in rank_candidates. Advanced
        model/MMGBSA execution is never automatically scheduled.
        """
        manifest = self._manifest(campaign_id)
        if manifest["mode"] == "legacy_campaign":
            return {"status": "legacy_broad_only", "created_stages": []}
        created = []
        profile = manifest["plan"].get("policy", {}).get("compute_profile")
        if profile:
            validate_profile(profile)
        with self._lock(campaign_id) as directory:
            state = _sealed_read(directory / "state.json")
            if "broad_docking" not in state["stages"]:
                return {"status": "receptor_held", "created_stages": []}
            broad_rows, missing = self._stage_rows(directory, state, "broad_docking")
            missing = [
                cid
                for cid in missing
                if manifest["prepared"].get(cid, {}).get("pdbqt_sha256")
            ]
            if missing:
                return {
                    "status": "awaiting_verified_broad_receipts",
                    "missing_count": len(missing),
                    "created_stages": [],
                }
            fallback_ids = [
                r["candidate_id"]
                for r in broad_rows
                if r["status"] not in SUCCESS and r.get("engine") == "unidock"
            ]
            if fallback_ids:
                if "broad_cpu_fallback" not in state["stages"]:
                    spec = self._stage(
                        manifest,
                        "broad_cpu_fallback",
                        sorted(set(fallback_ids)),
                        [],
                        manifest["receptor"],
                        [42],
                        "vina",
                    )
                    spec["upstream_receipts_sha256"] = digest(state["outcomes"])
                    _sealed_write(directory / "stage-broad_cpu_fallback.json", spec)
                    state["stages"]["broad_cpu_fallback"] = digest(spec)
                    _sealed_write(directory / "state.json", state)
                    return {
                        "status": "cpu_fallback_frozen_not_submitted",
                        "created_stages": ["broad_cpu_fallback"],
                    }
                fallback_rows, missing = self._stage_rows(
                    directory, state, "broad_cpu_fallback"
                )
                if missing:
                    return {
                        "status": "awaiting_verified_cpu_fallback_receipts",
                        "missing_count": len(missing),
                        "created_stages": [],
                    }
                # Preserve original failed rows in their receipts; only the
                # separate CPU result can become a final broad-stage success.
                by_candidate = {r["candidate_id"]: r for r in broad_rows}
                for row in fallback_rows:
                    by_candidate[row["candidate_id"]] = row
                broad_rows = list(by_candidate.values())
            if "confirmation" not in state["stages"]:
                confirmation_config = state.get(
                    "confirmation_policy_override"
                ) or manifest["plan"].get("policy", {}).get("confirmation")
                # A missing-score sibling must not reserve the parent before
                # confirmation checks the actual verified GPU eligibility.
                ranked = self._rank_results(manifest, broad_rows, unique_parents=not bool(confirmation_config))
                successful_ids = {
                    r["candidate_id"] for r in broad_rows if r["status"] in SUCCESS
                }
                ids = [
                    r["candidate_id"]
                    for r in ranked
                    if r["candidate_id"] in successful_ids
                ][:288]  # Explicit historical policy; new plans carry a frozen confirmation profile.
                if confirmation_config:
                    ids = select_confirmation(ranked, broad_rows, confirmation_config)
                controls = [
                    r["candidate_id"]
                    for r in manifest["candidates"]
                    if r.get("role")
                    in {
                        "control",
                        "positive_control",
                        "negative_control",
                        "counter_control",
                    }
                    and manifest["prepared"]
                    .get(r["candidate_id"], {})
                    .get("pdbqt_sha256")
                ]
                spec = self._stage(
                    manifest,
                    "confirmation",
                    ids,
                    controls,
                    manifest["receptor"],
                    confirmation_config["seeds"] if confirmation_config else [137, 271],
                    "vina",
                )
                spec["selection_policy"] = {
                    "priority": PRIORITY_POLICY_VERSION,
                    "comparability": COMPARABILITY_VERSION,
                    "ranking_population": len(ranked),
                }
                if confirmation_config:
                    spec["selection_policy"]["confirmation"] = deepcopy(
                        confirmation_config
                    )
                spec["upstream_receipts_sha256"] = digest(state["outcomes"])
                _sealed_write(directory / "stage-confirmation.json", spec)
                state["stages"]["confirmation"] = digest(spec)
                created.append("confirmation")
            if state.get("confirmation_policy_override"):
                spec = self._stage_read(directory, state, "confirmation")
                if "confirmation_scope" not in state:
                    config = state["confirmation_policy_override"]
                    ids = self._amended_confirmation_ids(manifest, directory, state, spec, config)
                    state["confirmation_scope"] = {
                        "candidate_ids": ids,
                        "control_ids": spec["control_ids"],
                        "seeds": config["seeds"],
                        "policy": deepcopy(config),
                        "source_stage_sha256": digest(spec),
                        "upstream_receipts_sha256": digest(state["outcomes"]),
                        "amendment_id": state["confirmation_amendment"]["amendment_id"],
                    }
                    _sealed_write(directory / "state.json", state)
                self._retire_confirmation_chunks(campaign_id, state, spec)
                _sealed_write(directory / "state.json", state)
            confirmation_rows, missing = self._stage_rows(
                directory, state, "confirmation"
            )
            if not missing:
                ranked = self._rank_sensitivity(manifest, broad_rows, confirmation_rows)
                successful = {
                    r["candidate_id"]
                    for r in confirmation_rows
                    if r["status"] in SUCCESS
                }

                sensitivity_limit = (
                    profile["counts"]["receptor_sensitivity"] if profile else 48
                )
                ids = [r["candidate_id"] for r in self._unique_ranked_parents(
                    ranked, eligible_ids=successful)][:sensitivity_limit]
                extras = (
                    manifest["plan"]
                    .get("target_context", {})
                    .get("additional_receptors", [])
                )
                # Planning stores context in job requests in the public contract.
                if not extras:
                    for planned in manifest["plan"].get("jobs", []):
                        if planned["request"].get("stage") == "receptor_sensitivity":
                            extras = planned["request"].get("additional_receptors", [])
                additional_limit = (
                    profile["counts"]["additional_receptors"] if profile else 2
                )
                for index, receptor in enumerate(extras[:additional_limit]):
                    name = "receptor_sensitivity_" + str(index + 1)
                    if name in state["stages"]:
                        continue
                    self._validate_receptor(
                        receptor, manifest["receptor"]["actual_protein_id"]
                    )
                    spec = self._stage(manifest, name, ids, [], receptor, [42], "vina")
                    spec["selection_policy"] = {
                        "version": "paired-protocol-rank-change-v1",
                        "priority_score_changed": False,
                        "candidates": {
                            row["candidate_id"]: row["receptor_sensitivity_selection"]
                            for row in ranked if row["candidate_id"] in ids
                        },
                    }
                    spec["upstream_receipts_sha256"] = digest(state["outcomes"])
                    _sealed_write(directory / ("stage-" + name + ".json"), spec)
                    state["stages"][name] = digest(spec)
                    created.append(name)
            _sealed_write(directory / "state.json", state)
        return {
            "status": "stages_frozen_not_submitted",
            "created_stages": created,
            "automatic_deep_count": 0,
        }

    @staticmethod
    def _rank_sensitivity(manifest, broad_results, confirmation_results):
        """Compare ranks over the same candidates, keeping engine scores separate.

        Each protocol pair uses only its shared, successful candidate population;
        comparing the full broad population with a selected confirmation subset
        would manufacture rank changes. Execution batch membership does not form
        a scientific cohort. Missing/singleton comparisons remain unmeasured.
        This changes sensitivity selection only, never the primary priority score.
        """
        ranked = CampaignCoordinator._rank_results(
            manifest, broad_results, unique_parents=False
        )
        candidate_ids = {
            row["candidate_id"] for row in ranked
            if row.get("role", "candidate") == "candidate"
            and not row.get("opposing_direction")
        }

        def groups(results):
            grouped = defaultdict(lambda: defaultdict(list))
            for result in results:
                cid, value = result["candidate_id"], finite(result.get("affinity"))
                if result["status"] not in SUCCESS or value is None or cid not in candidate_ids:
                    continue
                protocol = comparable_protocol(result["protocol"])
                grouped[digest(protocol)][cid].append(-value)
            return grouped

        comparisons = defaultdict(list)
        broad, confirmation = groups(broad_results), groups(confirmation_results)
        for broad_id, broad_values in sorted(broad.items()):
            for confirmation_id, confirmation_values in sorted(confirmation.items()):
                shared = sorted(set(broad_values) & set(confirmation_values))
                if len(shared) < 2:
                    continue
                broad_ranks = _percentiles({cid: median(broad_values[cid]) for cid in shared})
                confirmation_ranks = _percentiles({cid: median(confirmation_values[cid]) for cid in shared})
                population_hash = digest(shared)
                for cid in shared:
                    comparisons[cid].append({
                        "broad_protocol_id": broad_id,
                        "confirmation_protocol_id": confirmation_id,
                        "population_sha256": population_hash,
                        "cohort_size": len(shared),
                        "broad_percentile": broad_ranks[cid],
                        "confirmation_percentile": confirmation_ranks[cid],
                        "absolute_rank_change": abs(broad_ranks[cid] - confirmation_ranks[cid]),
                    })
        for row in ranked:
            values = comparisons[row["candidate_id"]]
            row["receptor_sensitivity_selection"] = {
                "status": "compared" if values else "not_evaluated",
                "maximum_rank_change": max((v["absolute_rank_change"] for v in values), default=None),
                "comparisons": values,
                "interpretation": "paired within-protocol rank change; not affinity or a statistical interval",
            }

        def order(row):
            change = row["receptor_sensitivity_selection"]["maximum_rank_change"]
            # Observed disagreement first, then unresolved comparisons, then
            # rank-stable candidates. Frozen primary priority breaks ties.
            group = 0 if change is not None and change > 0 else 1 if change is None else 2
            return (group, -(change or 0), not has_purchase_record(row),
                    -row["priority_score_v2"], row["candidate_id"])

        ranked.sort(key=order)
        return ranked

    @staticmethod
    def _rank_results(manifest, results, *, unique_parents=True):
        rows = deepcopy(manifest["candidates"])
        by_id = {row["candidate_id"]: row for row in rows}
        generalized = (
            manifest["plan"].get("policy", {}).get("compute_profile") or {}
        ).get("version") == "generalized-research-v1"
        primary = set()
        if generalized:
            for result in results:
                if result["status"] in SUCCESS and result.get("engine") == "unidock":
                    primary.add(digest(comparable_protocol(result["protocol"])))
            if len(primary) > 1:
                raise ValueError(
                    "Primary GPU results span incompatible frozen protocols"
                )
        for row in rows:
            row["docking"] = []
            if generalized:
                row["priority_primary_protocol_id"] = next(iter(primary), None)
        for result in results:
            if result["status"] not in SUCCESS:
                continue
            protocol = result["protocol"]
            by_id[result["candidate_id"]]["docking"].append(
                {
                    "target_id": manifest["plan"]["target_id"],
                    "receptor_id": protocol["receptor_sha256"],
                    "pocket_id": digest([protocol["center"], protocol["size"]]),
                    "engine": result["engine"],
                    "protocol_id": digest(comparable_protocol(protocol)),
                    "unit": "kcal/mol",
                    "score": result["affinity"],
                    "higher_is_better": False,
                    "artifact_valid": True,
                }
            )
        ranked = rank_candidates(rows)
        ranked.sort(
            key=lambda row: (
                not has_purchase_record(row),
                -row["priority_score_v2"],
                row["candidate_id"],
            )
        )
        if not unique_parents:
            # Scores keep their original full frozen comparison population.
            # The consumer applies execution eligibility before parent choice.
            return ranked
        return CampaignCoordinator._unique_ranked_parents(ranked)

    @staticmethod
    def _unique_ranked_parents(ranked, *, eligible_ids=None):
        """Keep score/order fixed; apply computation eligibility before parent choice."""
        # Controls disappear from the candidate ranking, but their chemistry
        # remains reserved before the ranked pool reaches confirmation.
        seen = {
            candidate_identity(row, trust_stored=True)["parent_id"] for row in ranked
            if row.get("role") in {"control", "positive_control", "negative_control", "counter_control"}
            or row.get("opposing_direction")
        }
        out = []
        for row in ranked:
            if eligible_ids is not None and row["candidate_id"] not in eligible_ids:
                continue
            parent = candidate_identity(row, trust_stored=True)["parent_id"]
            if (
                parent
                and parent not in seen
                and row.get("role")
                not in {
                    "control",
                    "positive_control",
                    "negative_control",
                    "counter_control",
                }
            ):
                out.append(row)
                seen.add(parent)
        return out

    def status(self, campaign_id, _locked=False):
        return self._status_and_plan(campaign_id, _locked=_locked)[0]

    def _status_and_plan(self, campaign_id, _locked=False, *, verified_assets=None):
        """Read one verified snapshot; the API also needs its small plan fields."""
        if not _locked:
            with self._lock(campaign_id):
                return self._status_and_plan(
                    campaign_id, _locked=True, verified_assets=verified_assets
                )
        directory = self._directory(campaign_id)
        manifest = _status_projection(
            directory / "manifest.json", campaign_id, "manifest"
        )
        state = _sealed_read(directory / "state.json")
        stages = {}
        held_ids = {item["candidate_id"] for item in manifest["held"]}
        job_statuses = {}
        pending_stages = {}
        pending_owners = {}
        if verified_assets is None:
            verified_assets = (
                set()
            )  # Per request only; the next status verifies bytes again.
        submissions_by_task = {}
        if not isinstance(state.get("submissions"), dict):
            raise ArtifactIntegrityError("Campaign checkpoint submissions are quarantined: invalid registry")
        for submission in state["submissions"].values():
            if not isinstance(submission, dict):
                raise ArtifactIntegrityError("Campaign checkpoint submissions are quarantined: invalid record")
            if submission.get("job_id"):
                if (not isinstance(submission.get("task_key"), str) or not submission["task_key"]
                        or not isinstance(submission.get("candidate_ids"), list)
                        or not isinstance(submission.get("owner"), str)
                        or type(submission.get("generation")) is not int):
                    raise ArtifactIntegrityError("Campaign checkpoint submissions are quarantined: missing execution identity")
                submissions_by_task.setdefault(submission["task_key"], []).append(
                    submission
                )
        for name in state["stages"]:
            spec = self._stage_read(directory, state, name)
            candidate_ids, control_ids, seeds = self._stage_members(state, name, spec)
            members = set(candidate_ids) | set(control_ids)
            rows, missing = self._stage_rows(
                directory,
                state,
                name,
                spec=spec,
                verified_assets=verified_assets,
                status_only=True,
            )
            recorded_keys = {
                (row.get("protocol", {}).get("seed"), row["candidate_id"])
                for row in rows
            }
            pending_stages[name] = []
            for seed in seeds:
                task_key = self._task_key(spec, seed)
                outcomes = state["outcomes"].get(task_key, {})
                latest = {}
                for submission in submissions_by_task.get(task_key, ()):
                    for cid in submission["candidate_ids"]:
                        if cid not in members:
                            continue
                        if (
                            cid not in latest
                            or submission["generation"] > latest[cid]["generation"]
                        ):
                            latest[cid] = submission
                for cid, submission in latest.items():
                    # Failed/corrupt receipts need an explicit retry; never a fabricated success.
                    if cid in outcomes and (seed, cid) in recorded_keys:
                        continue
                    job_id = submission["job_id"]
                    pending_stages[name].append((seed, cid, job_id))
                    pending_owners.setdefault(submission["owner"], set()).add(job_id)
            successes = {}
            for row in rows:
                if row["status"] in SUCCESS:
                    successes[row["candidate_id"]] = (
                        successes.get(row["candidate_id"], 0) + 1
                    )
            stages[name] = {
                "verified_complete": sum(
                    successes.get(cid, 0) == len(seeds) for cid in candidate_ids
                ),
                "verified_controls_complete": sum(
                    successes.get(cid, 0) == len(seeds) for cid in control_ids
                ),
                "verified_execution_receipts": sum(successes.values()),
                "planned_candidate_count": len(candidate_ids),
                "planned_control_count": len(control_ids),
                "required_seed_count": len(seeds),
                "terminal_missing_count": 0,
                "terminal_missing_statuses": [],
                "recorded_failures_or_review": sum(
                    r["status"] not in SUCCESS for r in rows
                ),
                "missing_or_corrupt": len(missing),
                "computable_missing_or_corrupt": sum(
                    cid not in held_ids for cid in missing
                ),
                "unprepared_missing_count": sum(cid in held_ids for cid in missing),
                "candidate_receipt_count": len(rows),
            }
            if name == "confirmation" and state.get("confirmation_scope"):
                stages[name]["execution_policy"] = deepcopy(
                    state["confirmation_scope"]["policy"]
                )
                stages[name]["original_frozen_candidate_count"] = len(
                    spec["candidate_ids"]
                )
                stages[name]["original_frozen_seed_count"] = len(spec["seeds"])
        if self.queue:
            for owner, job_ids in pending_owners.items():
                if hasattr(self.queue, "get_statuses"):
                    job_statuses.update(
                        self.queue.get_statuses(sorted(job_ids), owner=owner)
                    )
                else:
                    job_statuses.update(
                        {
                            job_id: self.queue.get(job_id, owner=owner)["status"]
                            for job_id in job_ids
                        }
                    )
        for name, pending in pending_stages.items():
            terminal = {
                (seed, cid)
                for seed, cid, job_id in pending
                if job_statuses.get(job_id) in TERMINAL
            }
            stages[name]["terminal_missing_count"] = len(terminal)
            stages[name]["terminal_missing_statuses"] = sorted(
                {
                    job_statuses[job_id]
                    for _, _, job_id in pending
                    if job_statuses.get(job_id) in TERMINAL
                }
            )
        result = {
            "campaign_id": campaign_id,
            "mode": manifest["mode"],
            "status": "legacy_handoff_held"
            if manifest["mode"] == "legacy_campaign"
            else "frozen_campaign",
            "population_count": manifest["population_count"],
            "held_preparations": manifest["held"],
            "confirmation_amendment": deepcopy(state.get("confirmation_amendment")),
            "stages": stages,
            "submitted_chunks": len(state["submissions"]),
            "quarantined_receipts": len(state["quarantined"]),
            "automatic_deep_count": 0,
            "retry_required": any(s["terminal_missing_count"] for s in stages.values()),
            "retry_action": "dispatch the affected stage with retry_failed=true; original failures remain retained",
            "scientific_validation": "NOT_ESTABLISHED",
        }
        return result, manifest["plan"]
