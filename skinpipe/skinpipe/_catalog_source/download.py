"""Resumable official-source downloads. Admission requires integrity evidence."""

import concurrent.futures
import gzip
import hashlib
import json
from pathlib import Path
import re
import subprocess
import uuid
import zipfile
from .common import atomic_json, lock, now, sha256


def validate_container(path):
    if str(path).endswith(".gz"):
        with gzip.open(path, "rb") as f:
            for _ in iter(lambda: f.read(8 << 20), b""):
                pass
    elif str(path).endswith((".zip", ".xlsx")):
        with zipfile.ZipFile(path) as z:
            bad = z.testzip()
            if bad:
                raise ValueError(f"CRC failed: {bad}")


def _curl(cmd, cancel_event):
    """Curl inherits the stage process group; stop only this owned child."""
    with subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    ) as process:
        try:
            while True:
                if cancel_event is not None and cancel_event.is_set():
                    raise concurrent.futures.CancelledError(
                        "Another source part failed"
                    )
                try:
                    stdout, stderr = process.communicate(timeout=0.25)
                    return subprocess.CompletedProcess(
                        cmd, process.returncode, stdout, stderr
                    )
                except subprocess.TimeoutExpired:
                    pass
        except BaseException:
            process.terminate()
            try:
                process.communicate(timeout=5)
            except subprocess.TimeoutExpired:
                process.kill()
                process.communicate()
            raise


def _range_rejected(headers):
    if not headers.exists():
        return False
    with headers.open("rb") as stream:
        stream.seek(max(0, headers.stat().st_size - 65536))
        return bool(re.search(rb"(?m)^HTTP/\S+\s+416(?:\s|$)", stream.read()))


def _preserve_download_attempt(root, part, headers, previous_receipt, result, reason):
    """Retain rejected bytes before one empty-file retry; never relax identity."""
    folder = part.parent / ("rejected-attempt-" + uuid.uuid4().hex)
    folder.mkdir(mode=0o700)
    observed = sha256(part) if part.exists() else None
    for file in (part, headers):
        if file.exists():
            file.replace(folder / file.name)
    if previous_receipt is not None:
        (folder / "previous-receipt.json").write_bytes(previous_receipt)
    attempt = {
        "directory": str(folder.relative_to(root)),
        "reason": reason,
        "observed_sha256": observed,
        "expected_sha256": result.get("sha256"),
        "preserved_at": now(),
    }
    atomic_json(folder / "attempt.json", {**result, **attempt, "status": "rejected"})
    result.setdefault("recovery_attempts", []).append(attempt)


def fetch(root, spec, cancel_event=None):
    root = Path(root)
    sid = spec["id"]
    if any(
        not isinstance(value, str)
        or value in {"", ".", ".."}
        or Path(value).name != value
        or "/" in value
        or "\\" in value
        for value in (sid, spec["filename"])
    ):
        raise ValueError("Source ID and filename must be single path components")
    output = root / "raw" / sid / spec["filename"]
    record = root / "manifests" / "downloads" / f"{sid}.json"
    if output.resolve() != output.absolute() or record.resolve() != record.absolute():
        raise ValueError("Source paths must not follow symlinks")
    output.parent.mkdir(parents=True, exist_ok=True)
    with lock(root / "logs" / f"download-{sid}.lock"):
        previous_receipt = record.read_bytes() if record.exists() else None
        if "sha256" in spec and not re.fullmatch(r"[0-9a-f]{64}", str(spec["sha256"])):
            result = {
                **spec,
                "status": "failed",
                "error": "Invalid frozen source SHA-256",
                "failed_at": now(),
            }
            atomic_json(record, result)
            return result
        if cancel_event is not None and cancel_event.is_set():
            return {
                **spec,
                "status": "cancelled",
                "error": "Another source part failed",
            }
        if record.exists() and output.exists():
            old = json.loads(record.read_text())
            if (
                old.get("status") == "verified"
                and old.get("url") == spec["url"]
                and old.get("id") == sid
                and old.get("filename") == spec["filename"]
                and old.get("path") == str(output.relative_to(root))
                and sha256(output) == old["sha256"]
            ):
                if spec.get("sha256") and old["sha256"] != spec["sha256"]:
                    raise ValueError(
                        "Cached artifact differs from frozen source SHA-256: " + sid
                    )
                return old
        result = {
            **spec,
            "started_at": now(),
            "path": str(output.relative_to(root)),
            "status": "downloading",
        }
        atomic_json(record, result)
        part = output.with_name(output.name + ".part")
        headers = output.with_name(output.name + ".headers")
        try:
            if any(p.resolve() != p.absolute() for p in (part, headers)):
                raise ValueError("Partial download paths must not follow symlinks")
            cmd = [
                "curl",
                "--proto",
                "=https",
                "--proto-redir",
                "=https",
                "--fail",
                "--location",
                "--retry",
                "6",
                "--retry-delay",
                "5",
                "--connect-timeout",
                "30",
                "--speed-time",
                "120",
                "--speed-limit",
                "1024",
                "--silent",
                "--show-error",
                "--dump-header",
                str(headers),
                "--continue-at",
                "-",
                "--output",
                str(part),
                spec["url"],
            ]
            resume_option = cmd.index("--continue-at")
            fresh_cmd = cmd[:resume_option] + cmd[resume_option + 2 :]
            for attempt in range(2):
                if cancel_event is not None and cancel_event.is_set():
                    raise concurrent.futures.CancelledError("Download was cancelled")
                p = _curl(cmd if attempt == 0 else fresh_cmd, cancel_event)
                actual = sha256(part) if part.exists() else None
                range_rejected = p.returncode == 33 or _range_rejected(headers)
                verified_complete_partial = (
                    range_rejected
                    and bool(spec.get("sha256"))
                    and actual == spec["sha256"]
                )
                if verified_complete_partial:
                    result["complete_partial_admitted_after_range_rejection"] = True
                if range_rejected and not verified_complete_partial:
                    if attempt == 0:
                        _preserve_download_attempt(
                            root, part, headers, previous_receipt, result,
                            "range_rejected",
                        )
                        continue
                    raise RuntimeError(
                        f"curl {p.returncode}: range rejected after one fresh retry"
                    )
                if p.returncode and not verified_complete_partial:
                    raise RuntimeError(f"curl {p.returncode}: {p.stderr[-2000:]}")
                if spec.get("sha256") and actual != spec["sha256"]:
                    if attempt == 0:
                        _preserve_download_attempt(
                            root, part, headers, previous_receipt, result,
                            "frozen_sha256_mismatch",
                        )
                        continue
                    raise ValueError(
                        "Frozen source SHA-256 mismatch after one fresh retry. "
                        "Downloaded bytes may be corrupted or the provider snapshot "
                        "may have changed; the frozen identity was not relaxed."
                    )
                break
            expected = spec.get("checksum")
            if expected:
                algo, wanted = expected.split(":", 1)
                h = hashlib.new(algo)
                with part.open("rb") as f:
                    for b in iter(lambda: f.read(8 << 20), b""):
                        h.update(b)
                if h.hexdigest().lower() != wanted.lower():
                    raise ValueError(f"{algo} mismatch: {h.hexdigest()} != {wanted}")
            # Validate the format before exposing the final filename.
            if spec["filename"].endswith((".zip", ".xlsx")):
                with zipfile.ZipFile(part) as z:
                    if z.testzip():
                        raise ValueError("ZIP CRC failure")
            elif spec["filename"].endswith(".gz"):
                with gzip.open(part, "rb") as f:
                    for _ in iter(lambda: f.read(8 << 20), b""):
                        pass
            part.replace(output)
            result.update(
                status="verified",
                sha256=actual,
                bytes=output.stat().st_size,
                verified_at=now(),
                provider_checksum_available=bool(expected),
            )
        except concurrent.futures.CancelledError as e:
            result.update(status="cancelled", error=str(e), cancelled_at=now())
        except Exception as e:
            result.update(status="failed", error=str(e), failed_at=now())
        atomic_json(record, result)
        print(
            json.dumps({k: result.get(k) for k in ["id", "status", "bytes", "error"]}),
            flush=True,
        )
        return result


def download_all(root, only=None):
    specs = json.loads((root / "config" / "sources.json").read_text())
    specs = [
        s
        for s in specs
        if s.get("admission") == "allowed"
        and s.get("url")
        and (not only or s["id"] in only)
    ]
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        return list(pool.map(lambda s: fetch(root, s), specs))
