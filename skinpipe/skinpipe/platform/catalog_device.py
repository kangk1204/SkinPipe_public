"""Choose one numerical execution channel per catalog, before committing rows.

CUDA availability alone is insufficient: an allocation and kernel must work,
then the caller loads and warms up the actual pinned model. An existing channel
never changes devices silently, including when a previously absent GPU appears.
"""

from __future__ import annotations

import gc


class CatalogDeviceError(RuntimeError):
    """A runtime failure, not a prediction about a compound."""


def release_cuda(torch):
    """Release this process's unused CUDA allocations; never touch other jobs."""
    gc.collect()
    try:
        if torch.cuda.is_initialized():
            torch.cuda.empty_cache()
    except Exception:
        # A failed CUDA context cannot always be flushed; process teardown is
        # still owned by the bootstrap supervisor.
        pass


def _runtime(torch, device):
    value = {
        "device": device,
        "torch": str(torch.__version__),
        "cuda_runtime": torch.version.cuda,
        "cpu_threads": 2,
        "interop_threads": 1,
        "inference_batch_size": 512 if device == "cuda" else 128,
        "tf32": False,
        "float32_matmul_precision": "highest",
    }
    if device == "cuda":
        try:
            if not torch.cuda.is_available():
                raise CatalogDeviceError("CUDA is unavailable")
            properties = torch.cuda.get_device_properties(0)
            # Run a real CUDA kernel, not just a driver-presence check.
            sample = torch.ones((8, 8), device="cuda:0")
            if float((sample @ sample).sum().item()) != 512.0:
                raise CatalogDeviceError("CUDA kernel result differs")
            torch.cuda.synchronize(0)
            del sample
            value["gpu"] = {
                "name": properties.name,
                "capability": [properties.major, properties.minor],
                "total_memory": properties.total_memory,
                "logical_device": 0,
            }
            value["cudnn_version"] = torch.backends.cudnn.version()
        except Exception as error:
            reason = (
                str(error)
                if isinstance(error, CatalogDeviceError)
                else type(error).__name__
            )
            raise CatalogDeviceError("CUDA probe failed: " + reason) from None
    return value


def select_model(torch, factory, policy="auto", existing=None):
    """Return (model, immutable channel), with fallback only before first use.

    ``factory(device, inference_batch_size)`` must load and warm up the real model.
    On resume ``existing`` is the persisted return value. Failure at that point
    preserves existing rows and asks for the same runtime, rather than recording a
    different device under the original model ID.
    """
    if policy not in {"auto", "cpu", "cuda"}:
        raise ValueError("Device must be auto, cpu, or cuda")
    if existing is not None:
        if existing.get("policy") != policy or existing.get("device") not in {
            "cpu",
            "cuda",
        }:
            raise CatalogDeviceError(
                "Saved execution policy differs; use a new directory"
            )
        runtime = _runtime(torch, existing["device"])
        if any(existing.get(key) != value for key, value in runtime.items()):
            raise CatalogDeviceError(
                "Saved numerical runtime differs; restore it or use a new directory"
            )
        try:
            return factory(runtime["device"], runtime["inference_batch_size"]), existing
        except Exception:
            release_cuda(torch)
            raise CatalogDeviceError(
                "Saved execution channel cannot run; committed predictions are preserved. "
                "Restore that device/runtime or build in a new directory."
            ) from None
    reason = None
    if policy != "cpu":
        try:
            runtime = _runtime(torch, "cuda")
            model = factory("cuda", runtime["inference_batch_size"])
            return model, {**runtime, "policy": policy, "fallback_reason": None}
        except Exception as error:
            reason = (
                str(error)
                if isinstance(error, CatalogDeviceError)
                else "CUDA model initialization/warmup failed: " + type(error).__name__
            )
        # The exception traceback is out of scope before GC, so a failed model's
        # own tensors can be collected before CPU fallback.
        release_cuda(torch)
        if policy == "cuda":
            raise CatalogDeviceError(reason + "; explicit CUDA mode does not fall back")
    runtime = _runtime(torch, "cpu")
    try:
        model = factory("cpu", runtime["inference_batch_size"])
    except Exception as error:
        raise CatalogDeviceError(
            "CPU model initialization/warmup failed: " + type(error).__name__
        ) from None
    return model, {**runtime, "policy": policy, "fallback_reason": reason}


def is_cuda_oom(error):
    return (
        "out of memory" in str(error).lower()
        or type(error).__name__ == "OutOfMemoryError"
    )


def is_cuda_runtime_failure(error):
    message = str(error).lower()
    return any(
        token in message
        for token in (
            "cuda",
            "cublas",
            "cudnn",
            "device-side assert",
            "illegal memory access",
            "driver",
            "no kernel image",
            "gpu",
            "nccl",
        )
    )
