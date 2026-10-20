"""Measured admission; resource scarcity never changes candidate validity."""

from __future__ import annotations

import math
import hashlib
import os
import platform
import shutil
import subprocess
import time

from .ranking import finite


def cpu_limit(available=None):
    available = (
        len(os.sched_getaffinity(0))
        if available is None and hasattr(os, "sched_getaffinity")
        else (available or os.cpu_count() or 1)
    )
    return max(1, min(8, int(available) // 4))


def gpu_admission(free_mib, peak_mib, *, reserve_mib=1536, safety_factor=1.25):
    if (
        isinstance(safety_factor, bool)
        or not isinstance(safety_factor, (int, float))
        or finite(safety_factor) is None
        or safety_factor < 1
    ):
        raise ValueError("Safety factor must be finite and at least one")
    if peak_mib is None:
        return {"admitted": False, "reason": "canary_required", "required_mib": None}
    if any(
        isinstance(x, bool)
        or not isinstance(x, (int, float))
        or finite(x) is None
        or x < 0
        for x in (free_mib, peak_mib, reserve_mib)
    ):
        raise ValueError("Memory measurements must be nonnegative and finite")
    required = math.ceil(peak_mib * safety_factor + reserve_mib)
    return {
        "admitted": free_mib >= required,
        "required_mib": required,
        "reason": "ready" if free_mib >= required else "insufficient_free_vram",
    }


def next_batch_size(current, *, failure=None, maximum=32):
    if any(type(value) is not int or value < 1 for value in (current, maximum)):
        raise ValueError("Positive batch sizes required")
    if failure == "oom":
        return max(1, current // 2)
    if failure:
        return current  # chemistry errors are not a memory signal
    return min(maximum, current * 2)


def classify_failure(message):
    text = str(message).lower()
    if any(
        x in text
        for x in ("out of memory", "cudaerrormemoryallocation", "cuda error: 2")
    ):
        return "oom"
    if any(
        x in text
        for x in ("atom type", "parse error", "invalid smiles", "unsupported chemistry")
    ):
        return "chemistry"
    if "timeout" in text or "timed out" in text:
        return "timeout"
    return "engine_failure"


def probe(*, directory=".", engine_paths=None):
    result = {
        "observed_at": time.time(),
        "hostname": platform.node(),
        "python": platform.python_version(),
        "cpu_count": os.cpu_count(),
        "cpu_available": len(os.sched_getaffinity(0)) if hasattr(os,'sched_getaffinity') else (os.cpu_count() or 1),
        "cpu_limit": cpu_limit(),
        "disk_free_bytes": shutil.disk_usage(directory).free,
        "gpus": [],
        "engines": {},
    }
    try:
        from pathlib import Path
        machine = Path('/etc/machine-id').read_text().strip()
    except OSError:
        machine = platform.node()
    result['physical_host_id'] = hashlib.sha256(machine.encode()).hexdigest()
    binary = shutil.which("nvidia-smi")
    if binary:
        try:
            p = subprocess.run(
            [
                binary,
                "--query-gpu=uuid,name,memory.total,memory.free,utilization.gpu,driver_version",
                "--format=csv,noheader,nounits",
            ],
            capture_output=True,
            text=True,
            timeout=10,
            )
        except (OSError, subprocess.TimeoutExpired) as error:
            p = None
            result["gpu_error"] = type(error).__name__
        if p is not None and p.returncode == 0:
            for line in p.stdout.splitlines():
                fields = [v.strip() for v in line.split(",")]
                if len(fields) == 6:
                    try:
                        measurements = [int(field) for field in fields[2:5]]
                    except ValueError:
                        result["gpu_error"] = "invalid GPU telemetry"
                        continue
                    result["gpus"].append(
                        dict(
                            zip(
                                (
                                    "uuid",
                                    "name",
                                    "total_mib",
                                    "free_mib",
                                    "utilization",
                                    "driver",
                                ),
                                (
                                    fields[0],
                                    fields[1],
                                    *measurements,
                                    fields[5],
                                ),
                            )
                        )
                    )
        elif p is not None:
            result["gpu_error"] = "nvidia-smi failed"
    for engine, configured in (engine_paths or {}).items():
        path = configured or shutil.which(engine)
        if not path:
            result["engines"][engine] = {"available": False}
            continue
        try:
            p = subprocess.run(
                [path, "--version"], capture_output=True, text=True, timeout=15
            )
            result["engines"][engine] = {
                "available": p.returncode == 0,
                "path": path,
                "version": (p.stdout + p.stderr).strip()[:500],
            }
        except (OSError, subprocess.TimeoutExpired) as exc:
            result["engines"][engine] = {
                "available": False,
                "error": type(exc).__name__,
            }
    return result
