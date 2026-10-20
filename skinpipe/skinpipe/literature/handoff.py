"""Selected-account handoff for a personally controlled self-hosted VM.

Official procedure:
https://developers.openai.com/siwc/token-sharing-open-source/self-hosted-vms

The caller must transport PreparedHandoff.payload through a protected channel
such as an SSH pipe, keep the VM provider disabled until local finalization,
and never print, log, persist or export the payload. This module performs no
network calls, revocation, OAuth, subprocess execution or secret-file export.
"""
from __future__ import annotations

from copy import deepcopy
from contextlib import contextmanager
from dataclasses import dataclass, field
import hmac
import json
import math
import os
import re
import stat
import uuid

from .chatgpt import ChatGPTClient, ChatGPTError, ISSUER, PLAN_SCOPE, _canonical, _sha
from . import usage as private_usage

SCHEMA = "skinpipe-chatgpt-selected-handoff-v1"
ACK_SCHEMA = "skinpipe-chatgpt-handoff-ack-v1"
MAX_HANDOFF_BYTES = 131072
_TOKENS = {"access_token", "refresh_token", "id_token", "expires_at", "earliest_refresh_at"}
_FIELDS = {"issuer", "subject", "client_id", "email", "name", "scopes", "access_token", "refresh_token",
           "id_token", "expires_at", "earliest_refresh_at", "token_type", "updated_at"}
_ACK_FIELDS = {"schema", "transfer_id", "imported", "remote_host_preserved", "selected_account_count"}


@dataclass(frozen=True, slots=True, repr=False)
class PreparedHandoff:
    """In-memory SECRET container; its repr never exposes credentials or IDs."""

    payload: bytes = field(repr=False)
    account_id: str = field(repr=False)
    transfer_id: str
    _record_sha256: str = field(repr=False)
    _source_directory: tuple[int, int] = field(repr=False)

    def __repr__(self):
        return "<PreparedHandoff SECRET; use protected transport only>"


def _error(code):
    raise ChatGPTError(code)


def _key(value):
    if not isinstance(value, str) or not re.fullmatch(r"[0-9a-f]{64}", value):
        _error("invalid_account_registration")
    return value


def _number(value):
    return type(value) in (int, float) and 0 <= value <= 253_402_300_799 and math.isfinite(value)


def _record(key, value):
    if (not isinstance(value, dict) or set(value) - _FIELDS or value.get("issuer") != ISSUER
            or not isinstance(value.get("subject"), str) or not 1 <= len(value["subject"]) <= 512
            or not isinstance(value.get("client_id"), str)
            or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]{0,255}", value["client_id"])
            or value["client_id"] == "dynamic_agent_client"
            or _key(key) != _sha(_canonical([ISSUER, value["client_id"], value["subject"]]))):
        _error("invalid_handoff_registration")
    scopes = value.get("scopes")
    if (not isinstance(scopes, list) or not 1 <= len(scopes) <= 32 or PLAN_SCOPE not in scopes
            or any(not isinstance(s, str) or not 1 <= len(s) <= 128 or any(c.isspace() for c in s) for s in scopes)):
        _error("chatgpt_plan_permission_required")
    for name in ("access_token", "refresh_token", "id_token"):
        if name in value or name != "id_token":
            token = value.get(name)
            if not isinstance(token, str) or not 1 <= len(token) <= 32768 or any(ord(c) < 32 for c in token):
                _error("invalid_handoff_credentials")
    if (value.get("token_type") != "Bearer" or not _number(value.get("expires_at"))
            or not _number(value.get("updated_at"))):
        _error("invalid_handoff_credentials")
    for name, maximum in (("email", 320), ("name", 256)):
        item = value.get(name)
        if item is not None and (not isinstance(item, str) or len(item) > maximum or any(ord(c) < 32 for c in item)):
            _error("invalid_handoff_registration")
    earliest = value.get("earliest_refresh_at")
    if earliest is not None and not (_number(earliest) or isinstance(earliest, str)
                                      and len(earliest) <= 128 and not any(ord(c) < 32 for c in earliest)):
        _error("invalid_handoff_credentials")
    return deepcopy(value)


def _existing(directory):
    try:
        metadata = os.stat("credentials.json", dir_fd=directory, follow_symlinks=False)
    except FileNotFoundError:
        _error("handoff_state_not_initialized")
    if not stat.S_ISREG(metadata.st_mode):
        _error("credential_file_not_private")


@contextmanager
def _transfer_lock(client):
    # A pending local OAuth must not install a new session just after the CAS.
    # Token refresh/profile selection already serialize on credentials.lock.
    with client._locked("sign-in.lock", blocking=False):
        with client._locked() as directory:
            yield directory


def _empty(client, directory):
    try:
        with client._locked("usage.lock") as usage_directory:
            if private_usage._read(usage_directory)["accounts"]:
                _error("handoff_destination_not_empty")
    except private_usage.UsageUnavailable:
        _error("handoff_destination_unavailable")
    state = client._read(directory)
    if (set(state) - {"schema", "ext_agent_host_id", "accounts", "active_account", "pending_client_id", "handoff_receipts"}
            or state["accounts"] or state.get("active_account") is not None or state.get("pending_client_id") is not None
            or state.get("handoff_receipts")):
        _error("handoff_destination_not_empty")
    return state


def initialize_empty_state(state_dir):
    """Create the VM's own fresh host ID; return no host ID or credentials."""
    with ChatGPTClient(state_dir) as client:
        with _transfer_lock(client) as directory:
            _empty(client, directory)
    return {"schema": "skinpipe-chatgpt-handoff-initialize-v1", "initialized": True,
            "empty_registration_state": True}


def prepare_selected_handoff(state_dir, *, account_id=None):
    """Read one connected registration, without refreshing or changing it."""
    with ChatGPTClient(state_dir, account_id=account_id) as client:
        with _transfer_lock(client) as directory:
            _existing(directory)
            state = client._read(directory)
            key = client._selected_key(state)
            if key is None:
                _error("sign_in_required")
            record = _record(key, state["accounts"].get(key))
            snapshot = _sha(_canonical(record))
            transfer_id = uuid.uuid4().hex
            payload = _canonical({"schema": SCHEMA, "transfer_id": transfer_id, "account_id": key,
                                  "record": record, "record_sha256": snapshot})
            if len(payload) > MAX_HANDOFF_BYTES:
                _error("handoff_payload_too_large")
            metadata = os.fstat(directory)
            return PreparedHandoff(payload, key, transfer_id, snapshot, (metadata.st_dev, metadata.st_ino))


def _decode(payload):
    if not isinstance(payload, bytes) or not 0 < len(payload) <= MAX_HANDOFF_BYTES:
        _error("invalid_handoff_payload")
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError("duplicate key")
            result[key] = value
        return result
    try:
        value = json.loads(payload, object_pairs_hook=pairs,
                           parse_constant=lambda _: (_ for _ in ()).throw(ValueError("nonfinite")))
        if (not isinstance(value, dict) or set(value) != {"schema", "transfer_id", "account_id", "record", "record_sha256"}
                or value["schema"] != SCHEMA or not isinstance(value["transfer_id"], str)
                or not re.fullmatch(r"[0-9a-f]{32}", value["transfer_id"])):
            raise ValueError("invalid envelope")
        key = _key(value["account_id"])
        record = _record(key, value["record"])
        if not isinstance(value["record_sha256"], str) or not hmac.compare_digest(value["record_sha256"], _sha(_canonical(record))):
            raise ValueError("changed record")
        return value, record
    except (ValueError, TypeError, UnicodeError, OverflowError, RecursionError):
        _error("invalid_handoff_payload")


def import_selected_handoff(state_dir, payload, *, expected_account_id=None):
    """Atomically import ONLY the selected registration into an initialized VM.

    The safe acknowledgement may be returned through SSH. It is not proof of
    current provider eligibility. The runtime must remain disabled until source
    finalization succeeds, so the VM alone will own future rotating refreshes.
    """
    value, record = _decode(payload)
    if expected_account_id is not None and _key(expected_account_id) != value["account_id"]:
        _error("account_registration_mismatch")
    acknowledgement = {"schema": ACK_SCHEMA, "transfer_id": value["transfer_id"], "imported": True,
                       "remote_host_preserved": True, "selected_account_count": 1}
    with ChatGPTClient(state_dir) as client:
        with _transfer_lock(client) as directory:
            _existing(directory)
            prior = _empty(client, directory)
            candidate = deepcopy(prior)
            candidate["accounts"] = {value["account_id"]: record}
            candidate["active_account"] = value["account_id"]
            receipts = candidate.setdefault("handoff_receipts", {})
            if not isinstance(receipts, dict) or value["transfer_id"] in receipts:
                _error("handoff_receipt_conflict")
            receipts[value["transfer_id"]] = {**acknowledgement, "direction": "imported"}
            try:
                client._write(directory, candidate)
                verified = client._read(directory)
                if verified != candidate or verified["ext_agent_host_id"] != prior["ext_agent_host_id"]:
                    _error("handoff_import_verification_failed")
            except Exception:
                try:
                    client._write(directory, prior)
                except Exception:
                    _error("handoff_import_rollback_failed")
                _error("handoff_import_failed")
    return acknowledgement


def finalize_local_handoff(state_dir, prepared, acknowledgement):
    """CAS the original selected snapshot, then remove local tokens WITHOUT revoke.

    Never finalize on a failure/unknown acknowledgement or after local refresh,
    reauthorization, account modification or state-directory replacement.
    """
    if (not isinstance(prepared, PreparedHandoff) or not isinstance(acknowledgement, dict)
            or set(acknowledgement) != _ACK_FIELDS or acknowledgement.get("schema") != ACK_SCHEMA
            or acknowledgement.get("transfer_id") != prepared.transfer_id
            or acknowledgement.get("imported") is not True or acknowledgement.get("remote_host_preserved") is not True
            or type(acknowledgement.get("selected_account_count")) is not int
            or acknowledgement["selected_account_count"] != 1):
        _error("handoff_acknowledgement_required")
    # Re-validate the in-memory secret envelope before using any private key.
    value, _ = _decode(prepared.payload)
    if (value["account_id"] != prepared.account_id or value["transfer_id"] != prepared.transfer_id
            or not hmac.compare_digest(value["record_sha256"], prepared._record_sha256)):
        _error("invalid_handoff_payload")
    with ChatGPTClient(state_dir, account_id=prepared.account_id) as client:
        with _transfer_lock(client) as directory:
            _existing(directory)
            metadata = os.fstat(directory)
            if (metadata.st_dev, metadata.st_ino) != prepared._source_directory:
                _error("handoff_source_changed")
            state = client._read(directory)
            current = state["accounts"].get(prepared.account_id)
            receipt = {"schema": "skinpipe-chatgpt-handoff-finalize-v1", "transfer_id": prepared.transfer_id,
                       "local_tokens_removed": True, "remote_revocation_performed": False,
                       "registration_mappings_preserved": True}
            previous_receipts = state.get("handoff_receipts", {})
            if not isinstance(previous_receipts, dict):
                _error("handoff_receipt_conflict")
            if prepared.transfer_id in previous_receipts:
                expected_record = {key: item for key, item in value["record"].items() if key not in _TOKENS}
                if current == expected_record and previous_receipts[prepared.transfer_id] == {**receipt, "direction": "exported"}:
                    return receipt  # Confirm an uncertain prior commit without repeating any mutation.
                _error("handoff_receipt_conflict")
            if not isinstance(current, dict) or not hmac.compare_digest(_sha(_canonical(current)), prepared._record_sha256):
                _error("handoff_source_changed")
            receipts = state.setdefault("handoff_receipts", {})
            if not isinstance(receipts, dict) or prepared.transfer_id in receipts:
                _error("handoff_receipt_conflict")
            for field_name in _TOKENS:
                current.pop(field_name, None)
            receipts[prepared.transfer_id] = {**receipt, "direction": "exported"}
            try:
                client._write(directory, state)
            except (OSError, ChatGPTError):
                _error("handoff_finalize_unconfirmed")
    return receipt
