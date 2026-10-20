"""Private, installation-local inference accounting; never a plan quota API.

Callers hold a cross-process usage lock and supply an already checked private
directory FD. Only registration keys, local UUIDs, counters and safe metadata
belong here. No credentials, prompts, article text or generated prose are stored.
"""
from __future__ import annotations

from datetime import datetime, timezone
import json
import math
import os
import re
import secrets
import stat
from types import SimpleNamespace

SCHEMA = "skinpipe-chatgpt-accounting-v1"
MAX_BYTES = 4 * 1024 * 1024
MAX_INTEGER = 9_007_199_254_740_991
MANAGE_URL = "https://chatgpt.com/settings/usage"
_CORE = ("input_tokens", "output_tokens", "total_tokens")
_QUOTA = {"subscription_sharing_usage_limit_exceeded", "insufficient_quota", "credit_balance_exhausted",
          "organization_spend_limit_exceeded", "project_spend_limit_exceeded", "organization_usage_limit_exceeded"}
_RATE = {"rate_limit_exceeded", "slow_down"}


class UsageUnavailable(RuntimeError):
    """A safe accounting failure which must never cause inference replay."""


def _integer(value):
    return value if type(value) is int and 0 <= value <= MAX_INTEGER else None


def timestamp(value):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not 0 <= value <= 253_402_300_799 or not math.isfinite(value):
        raise UsageUnavailable("invalid_accounting_time")
    return datetime.fromtimestamp(value, timezone.utc).isoformat().replace("+00:00", "Z")


def sanitize_usage(value):
    value = value if isinstance(value, dict) else {}
    result = {name: _integer(value.get(name)) for name in _CORE}
    input_details = value.get("input_tokens_details")
    output_details = value.get("output_tokens_details")
    input_details = input_details if isinstance(input_details, dict) else {}
    output_details = output_details if isinstance(output_details, dict) else {}
    result["input_tokens_details"] = {name: _integer(input_details.get(name)) for name in ("cached_tokens", "cache_write_tokens")}
    result["output_tokens_details"] = {"reasoning_tokens": _integer(output_details.get("reasoning_tokens"))}
    for name in ("cached_tokens", "cache_write_tokens"):
        if result["input_tokens"] is not None and result["input_tokens_details"][name] is not None and result["input_tokens_details"][name] > result["input_tokens"]:
            result["input_tokens_details"][name] = None
    if result["output_tokens"] is not None and result["output_tokens_details"]["reasoning_tokens"] is not None and result["output_tokens_details"]["reasoning_tokens"] > result["output_tokens"]:
        result["output_tokens_details"]["reasoning_tokens"] = None
    if all(result[name] is not None for name in _CORE) and result["total_tokens"] != result["input_tokens"] + result["output_tokens"]:
        result["total_tokens"] = None  # Never silently repair an inconsistent provider total.
    present = sum(result[name] is not None for name in _CORE)
    result["coverage"] = "REPORTED" if present == 3 else "PARTIAL" if present else "UNKNOWN"
    return result


def safe_limit(error, now):
    code = getattr(error, "code", None)
    if code not in {"usage_limit_exceeded", "rate_limited"}:
        return None
    provider = getattr(error, "provider_code", None)
    status = getattr(error, "http_status", None)
    delay = getattr(error, "retry_after_seconds", None)
    reset = getattr(error, "reset_at", None)
    if isinstance(delay, bool) or not isinstance(delay, (int, float)) or not 0 <= delay <= 31_536_000 or not math.isfinite(delay):
        delay = None
    if not isinstance(reset, str) or not re.fullmatch(r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(?:\.\d{1,6})?Z", reset):
        reset = None
    else:
        try:
            datetime.fromisoformat(reset.replace("Z", "+00:00"))
        except ValueError:
            reset = None
    return {"code": code, "usage_exhausted": code == "usage_limit_exceeded",
            "provider_code": provider if isinstance(provider, str) and provider in _QUOTA | _RATE else None,
            "http_status": status if type(status) is int and 100 <= status <= 599 else None,
            "retry_after_seconds": delay, "reset_at": reset, "observed_at": timestamp(now)}


def _is_timestamp(value):
    if not isinstance(value, str) or not re.fullmatch(r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(?:\.\d{1,6})?Z", value):
        return False
    try:
        datetime.fromisoformat(value.replace("Z", "+00:00"))
        return True
    except ValueError:
        return False


def _checked_file(descriptor):
    metadata = os.fstat(descriptor)
    if (not stat.S_ISREG(metadata.st_mode) or metadata.st_uid != os.getuid()
            or metadata.st_nlink != 1 or stat.S_IMODE(metadata.st_mode) & 0o077):
        raise UsageUnavailable("accounting_file_not_private")


def _checked_directory(directory):
    metadata = os.fstat(directory)
    if not stat.S_ISDIR(metadata.st_mode) or metadata.st_uid != os.getuid() or stat.S_IMODE(metadata.st_mode) & 0o077:
        raise UsageUnavailable("accounting_directory_not_private")


def _read(directory):
    _checked_directory(directory)
    try:
        descriptor = os.open("usage.json", os.O_RDONLY | os.O_NOFOLLOW, dir_fd=directory)
    except FileNotFoundError:
        return {"schema": SCHEMA, "accounts": {}}
    with os.fdopen(descriptor, "rb") as stream:
        _checked_file(stream.fileno())
        raw = stream.read(MAX_BYTES + 1)
    if len(raw) > MAX_BYTES:
        raise UsageUnavailable("accounting_size_limit")
    try:
        ledger = json.loads(raw)
        if not isinstance(ledger, dict) or ledger.get("schema") != SCHEMA or not isinstance(ledger.get("accounts"), dict):
            raise ValueError()
        for key, account in ledger["accounts"].items():
            if (not re.fullmatch(r"[0-9a-f]{64}", key) or not isinstance(account, dict)
                    or not isinstance(account.get("attempts"), dict) or not _is_timestamp(account.get("tracking_since"))):
                raise ValueError()
            if "last_limit_error" in account:
                previous = account["last_limit_error"]
                if not isinstance(previous, dict) or not _is_timestamp(previous.get("observed_at")):
                    raise ValueError()
                normalized = safe_limit(SimpleNamespace(**previous), 0)
                if normalized is None:
                    raise ValueError()
                normalized["observed_at"] = previous["observed_at"]
                account["last_limit_error"] = normalized
            for attempt_id, entry in account["attempts"].items():
                if (not re.fullmatch(r"[0-9a-f]{32}", attempt_id) or not isinstance(entry, dict)
                        or entry.get("status") not in {"started", "completed", "failed"}
                        or not _is_timestamp(entry.get("started_at"))
                        or entry["status"] != "started" and not _is_timestamp(entry.get("finished_at"))):
                    raise ValueError()
        return ledger
    except (ValueError, TypeError):
        raise UsageUnavailable("invalid_accounting_record") from None


def _write(directory, ledger):
    _checked_directory(directory)
    raw = json.dumps(ledger, sort_keys=True, ensure_ascii=False, allow_nan=False, separators=(",", ":")).encode()
    if len(raw) > MAX_BYTES:
        raise UsageUnavailable("accounting_size_limit")
    temporary = ".usage-" + secrets.token_hex(16)
    descriptor = os.open(temporary, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600, dir_fd=directory)
    try:
        with os.fdopen(descriptor, "wb") as stream:
            stream.write(raw)
            stream.flush()
            os.fsync(stream.fileno())
        try:
            existing = os.open("usage.json", os.O_RDONLY | os.O_NOFOLLOW, dir_fd=directory)
        except FileNotFoundError:
            pass
        else:
            try:
                _checked_file(existing)
            finally:
                os.close(existing)
        os.replace(temporary, "usage.json", src_dir_fd=directory, dst_dir_fd=directory)
        os.fsync(directory)
    finally:
        try:
            os.unlink(temporary, dir_fd=directory)
        except FileNotFoundError:
            pass


def begin(directory, account_key, attempt_id, now):
    if not isinstance(account_key, str) or not re.fullmatch(r"[0-9a-f]{64}", account_key) or not isinstance(attempt_id, str) or not re.fullmatch(r"[0-9a-f]{32}", attempt_id):
        raise UsageUnavailable("invalid_accounting_key")
    ledger = _read(directory)
    account = ledger["accounts"].setdefault(account_key, {"tracking_since": timestamp(now), "attempts": {}})
    if attempt_id not in account["attempts"]:
        account["attempts"][attempt_id] = {"started_at": timestamp(now), "status": "started"}
        _write(directory, ledger)


def finish(directory, account_key, attempt_id, now, *, completed, usage=None, error=None):
    ledger = _read(directory)
    account = ledger["accounts"].get(account_key)
    entry = account.get("attempts", {}).get(attempt_id) if account else None
    if entry is None:
        raise UsageUnavailable("accounting_attempt_missing")
    if entry["status"] != "started":
        return  # Idempotent finalization never doubles counted tokens/requests.
    entry.update(status="completed" if completed else "failed", finished_at=timestamp(now), usage=sanitize_usage(usage))
    limit = safe_limit(error, now)
    if limit:
        account["last_limit_error"] = limit
    _write(directory, ledger)


def summary(directory, account_key):
    ledger = _read(directory)
    account = ledger["accounts"].get(account_key, {})
    entries = list(account.get("attempts", {}).values())
    usages = [sanitize_usage(entry.get("usage")) for entry in entries]
    missing = sum(item["coverage"] != "REPORTED" for item in usages)
    def total(name, parent=None):
        values = [item[parent].get(name) if parent else item.get(name) for item in usages]
        known = [value for value in values if value is not None]
        return _integer(sum(known)) if known else None
    totals = {name: total(name) for name in _CORE}
    coverage = "COMPLETE_FOR_RECORDED_REQUESTS" if entries and not missing and all(value is not None for value in totals.values()) else "PARTIAL" if any(value is not None for value in totals.values()) else "UNKNOWN"
    return {"tracking_since": account.get("tracking_since"), "requests_started": len(entries),
            "responses_completed": sum(entry["status"] == "completed" for entry in entries),
            "failed_requests": sum(entry["status"] == "failed" for entry in entries),
            "pending_requests": sum(entry["status"] == "started" for entry in entries),
            "usage_missing_requests": missing, **totals,
            "input_tokens_details": {name: total(name, "input_tokens_details") for name in ("cached_tokens", "cache_write_tokens")},
            "output_tokens_details": {"reasoning_tokens": total("reasoning_tokens", "output_tokens_details")},
            "coverage": coverage, "accounting_status": "available"}, account.get("last_limit_error")


def unavailable_summary():
    return {"tracking_since": None, "requests_started": None, "responses_completed": None,
            "failed_requests": None, "pending_requests": None, "usage_missing_requests": None,
            **{name: None for name in _CORE}, "input_tokens_details": {"cached_tokens": None, "cache_write_tokens": None},
            "output_tokens_details": {"reasoning_tokens": None}, "coverage": "UNKNOWN", "accounting_status": "unavailable"}
