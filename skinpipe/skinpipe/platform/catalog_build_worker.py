"""Owned subprocess stages for catalog bootstrap; no deployment/service actions."""

from __future__ import annotations

from collections import Counter
from contextlib import contextmanager
import gzip
import importlib.metadata
import json
import math
import os
from pathlib import Path
import shutil
import sys
import time
import uuid

from .artifacts import atomic_bytes, sha256_file
from .identity import digest as identity_digest
from .._catalog_source.common import atomic_json, digest, now


def download(root, plan):
    from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
    from threading import Event

    from .._catalog_source.download import fetch
    from .catalog_bootstrap import ASSETS

    root = Path(root)
    progress = {
        "status": "running",
        "completed_sources": 0,
        "total_sources": len(plan["sources"]),
        "sources": {},
    }

    def report(source, **details):
        progress["current_source"] = source
        progress["sources"].setdefault(source, {}).update(details)
        progress["updated_at"] = now()
        atomic_json(root / "reports/download-progress.json", progress)

    for spec in plan["sources"]:
        sid = spec["id"]
        report(sid, status="downloading", completed_parts=0, total_parts=1)
        try:
            if sid != "zinc_instock":
                result = fetch(root, spec)
                if result["status"] != "verified":
                    raise ValueError(
                        "Official source unavailable or changed: "
                        + sid
                        + ": "
                        + result.get("error", result["status"])
                    )
                progress["completed_sources"] += 1
                report(sid, status="verified", completed_parts=1)
                continue
            inventory = ASSETS / "zinc-artifacts.json"
            if sha256_file(inventory) != plan["zinc_artifacts_sha256"]:
                raise ValueError("Pinned ZINC inventory changed")
            artifacts = json.loads(inventory.read_bytes())
            if not artifacts or len({p["id"] for p in artifacts}) != len(artifacts):
                raise ValueError("ZINC inventory is empty or has duplicate part IDs")
            report(sid, total_parts=len(artifacts), parallel_downloads=4)
            # Keep both concurrency and queued work bounded. Results are placed by
            # inventory index, so network completion order cannot change assembly.
            results = [None] * len(artifacts)
            stop = Event()
            pool = ThreadPoolExecutor(max_workers=4)
            pending = {}
            next_index = 0
            completed = 0
            try:
                while next_index < len(artifacts) or pending:
                    while next_index < len(artifacts) and len(pending) < 4:
                        part = artifacts[next_index]
                        pending[pool.submit(fetch, root, part, stop)] = next_index
                        next_index += 1
                    done, _ = wait(pending, return_when=FIRST_COMPLETED)
                    for future in sorted(done, key=pending.__getitem__):
                        index = pending.pop(future)
                        result = future.result()
                        if result["status"] != "verified":
                            report(sid, failed_part=artifacts[index]["id"])
                            raise ValueError(
                                "ZINC official artifact unavailable: "
                                + artifacts[index]["id"]
                                + ": "
                                + result.get("error", result["status"])
                            )
                        results[index] = result
                        completed += 1
                        report(sid, completed_parts=completed)
                        if completed % 50 == 0 or completed == len(artifacts):
                            print(
                                json.dumps(
                                    {
                                        "source": sid,
                                        "downloaded": completed,
                                        "total": len(artifacts),
                                    }
                                ),
                                flush=True,
                            )
            finally:
                # An error stops only this stage's downloads. Curl stays in the
                # stage process group for normal parent cancellation and reaping.
                stop.set()
                pool.shutdown(wait=True, cancel_futures=True)
            report(sid, status="assembling", completed_parts=completed)
            path = root / "raw" / sid / spec["filename"]
            path.parent.mkdir(parents=True, exist_ok=True)
            temp = path.with_suffix(path.suffix + ".part")
            count = 0
            with temp.open("wb") as output:
                for part in results:
                    source = root / part["path"]
                    # Receipts cannot redirect assembly outside the owned root.
                    if (
                        source.resolve() != source.absolute()
                        or root.resolve() not in source.resolve().parents
                    ):
                        raise ValueError("ZINC receipt path leaves the owned catalog")
                    last = ""
                    with gzip.open(source, "rt") as stream:
                        for line in stream:
                            last = line
                            count += bool(line.strip())
                    with source.open("rb") as stream:
                        shutil.copyfileobj(stream, output, 8 << 20)
                    if last and not last.endswith("\n"):
                        output.write(gzip.compress(b"\n", mtime=0))
                output.flush()
                os.fsync(output.fileno())
            if count != spec["expected_records"]:
                raise ValueError("ZINC source record count differs")
            with gzip.open(temp, "rt") as stream:
                if sum(bool(line.strip()) for line in stream) != count:
                    raise ValueError("ZINC assembled record boundaries differ")
            temp.replace(path)
            atomic_json(
                root / "manifests/downloads/zinc_instock.json",
                {
                    **spec,
                    "status": "verified",
                    "path": str(path.relative_to(root)),
                    "sha256": sha256_file(path),
                    "bytes": path.stat().st_size,
                    "derived_from": [
                        {"id": r["id"], "sha256": r["sha256"]} for r in results
                    ],
                    "verified_at": now(),
                    "transformation": "ordered gzip members with explicit record boundary v2",
                },
            )
            progress["completed_sources"] += 1
            report(sid, status="verified", records=count)
        except BaseException as error:
            progress["status"] = (
                "interrupted"
                if isinstance(error, (KeyboardInterrupt, SystemExit))
                else "failed"
            )
            report(sid, status=progress["status"], error=str(error))
            raise
    progress["status"] = "complete"
    atomic_json(
        root / "reports/download-progress.json", {**progress, "updated_at": now()}
    )


def predict_batch(model, items, endpoints):
    """Isolate molecule failures; a broken GPU is never a toxicity result."""
    from .catalog_device import (
        CatalogDeviceError,
        is_cuda_oom,
        is_cuda_runtime_failure,
        release_cuda,
    )

    if not items:
        return []
    cuda = getattr(model, "device", "cpu") == "cuda"
    oom = False
    try:
        smiles = [r["canonical_smiles"] for r in items]
        return _checked_predictions(model.predict(smiles=smiles), items, endpoints)
    except Exception as error:
        error_type = type(error).__name__
        oom = cuda and is_cuda_oom(error)
        if isinstance(error, CatalogDeviceError) or isinstance(error, MemoryError):
            raise CatalogDeviceError(
                "Inference runtime failed; committed batches are preserved"
            ) from None
        if (
            cuda
            and not oom
            and (isinstance(error, RuntimeError) or is_cuda_runtime_failure(error))
        ):
            raise CatalogDeviceError(
                "GPU inference failed; restore the saved runtime and resume"
            ) from None
    # Clear exception tracebacks before releasing CUDA buffers and retrying.
    if oom:
        import torch

        release_cuda(torch)
        if len(items) == 1:
            raise CatalogDeviceError(
                "GPU memory is insufficient for one input; existing batches are preserved. "
                "Free this job's resources and resume, or build a separate CPU catalog."
            )
    if len(items) == 1:
        return [("failed", {}, error_type)]
    mid = len(items) // 2
    return predict_batch(model, items[:mid], endpoints) + predict_batch(
        model, items[mid:], endpoints
    )


def _checked_predictions(frame, items, endpoints):
    smiles = [r["canonical_smiles"] for r in items]
    if (
        list(frame.index) != smiles
        or len(frame) != len(items)
        or len(frame.columns) != len(endpoints)
        or set(frame.columns) != set(endpoints)
    ):
        raise ValueError("Prediction output schema/order differs")
    results = []
    for values in frame.to_dict("records"):
        values = {
            k: float(v) if math.isfinite(float(v)) else None for k, v in values.items()
        }
        if any(
            v is not None and not 0 <= v <= 1
            for k, v in values.items()
            if endpoints[k]["task_type"] == "classification"
        ):
            results.append(("failed", {}, "invalid_classifier_probability"))
        else:
            results.append(
                (
                    "success"
                    if all(v is not None for v in values.values())
                    else "partial_nonfinite",
                    values,
                    None,
                )
            )
    return results


def _make_model(profile, endpoints, device, batch_size):
    from admet_ai import ADMETModel
    from chemprop.data import MoleculeDatapoint, MoleculeDataset
    from chemprop.data.dataloader import build_dataloader

    class CatalogADMETModel(ADMETModel):
        def _build_dataloader(self, mols):
            dataset = MoleculeDataset(data=[MoleculeDatapoint(mol=m) for m in mols])
            return build_dataloader(
                dataset=dataset, batch_size=batch_size, num_workers=0, shuffle=False
            )

    model = CatalogADMETModel(
        models_dir=profile["models_dir"],
        include_physchem=True,
        drugbank_path=None,
        num_workers=0,
    )
    # The pinned upstream constructor has no device parameter. It loads weights
    # on CPU, then its Trainer uses this attribute for every prediction.
    model.device = device
    if (
        model.drugbank is not None
        or model.drugbank_atc_filtered is not None
        or sum(map(len, model.task_lists)) != 41
    ):
        raise ValueError("Unexpected model tasks/reference dataset")
    warmup = [{"canonical_smiles": "CCO"}]
    result = _checked_predictions(model.predict(smiles=["CCO"]), warmup, endpoints)
    if any(status != "success" for status, _, _ in result):
        raise ValueError("Pinned model warmup returned invalid endpoint values")
    return model


@contextmanager
def _input_batches(con, batch_size):
    """Sort once and stream; receipt validation must use the other connection."""
    with con.cursor() as reader:
        with reader.execute(
            "SELECT compound_id,canonical_smiles,inference_status FROM compounds ORDER BY compound_id"
        ).fetch_record_batch(rows_per_batch=batch_size) as batches:
            yield (table.to_pylist() for table in batches)


def admet(root, plan):
    import duckdb
    import pyarrow as pa
    import pyarrow.parquet as pq
    from .research_toxicity import _profile
    from .._catalog_source.model_registry import archive_model
    from .._catalog_source.prediction_integrity import validate_batch
    from .._catalog_source.prediction_schema import PREDICTION_SCHEMA

    profile, identity, manifest, definitions, _ = _profile(
        json.loads(Path(plan["profile_path"]).read_bytes())
    )
    if identity != plan["model_identity"]:
        raise ValueError("Frozen model identity changed")
    for name, version in profile["packages"].items():
        if importlib.metadata.version(name) != version:
            raise ValueError("Runtime package differs: " + name)
    # Block all network before importing/running the model. Data downloads are
    # a separate stage and no user structures are transmitted by prediction.
    import socket

    def offline(*_args, **_kwargs):
        raise RuntimeError("ADMET inference network disabled")

    socket.create_connection = offline
    socket.socket.connect = offline
    socket.socket.connect_ex = offline
    import torch
    from .catalog_device import select_model

    torch.set_num_threads(2)
    torch.set_num_interop_threads(1)
    torch.set_float32_matmul_precision("highest")
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False
    endpoints = {r["id"]: r for r in definitions}
    if len(endpoints) != 52:
        raise ValueError("Expected 52 official endpoint definitions")
    execution_path = root / "manifests/catalog-execution.json"
    execution_binding = {
        "schema": "skinpipe-catalog-execution-v1",
        "plan_sha256": identity_digest(plan),
        "profile_identity": identity,
    }
    saved = None
    if execution_path.exists():
        execution = json.loads(execution_path.read_bytes())
        if set(execution) != {*execution_binding, "channel"} or any(
            execution.get(key) != value for key, value in execution_binding.items()
        ):
            raise ValueError(
                "Saved execution binding differs; preserve this catalog and use a new directory"
            )
        saved = execution["channel"]
    elif any((root / "data/predictions").rglob("*.parquet")):
        raise ValueError(
            "Predictions exist without their execution identity; refusing to overwrite"
        )
    model, channel = select_model(
        torch,
        lambda device, batch_size: _make_model(profile, endpoints, device, batch_size),
        policy=plan.get("device_policy", "auto"),
        existing=saved,
    )
    if saved is None:
        atomic_json(execution_path, {**execution_binding, "channel": channel})
    device_label = (
        "GPU (" + channel["gpu"]["name"] + ")"
        if channel["device"] == "cuda"
        else "CPU (기존 실행 조건 유지)"
        if saved is not None
        else "CPU (GPU를 사용할 수 없어 자동 전환)"
        if channel["fallback_reason"]
        else "CPU (직접 선택)"
    )
    print("ADMET 실행 장치: " + device_label, flush=True)
    print(json.dumps({"execution_channel": channel}), flush=True)
    manifest = {
        **manifest,
        "catalog_runner_sha256": sha256_file(__file__),
        "catalog_model_identity": identity,
        "preprocessing": "exact-rdkit-graph-v1",
        "device": channel["device"],
        "runtime_channel": "owner_local_" + channel["device"],
        "catalog_execution": channel,
        "catalog_batch_size": plan["batch_size"],
    }
    # Batch geometry/runner are numerical execution provenance, never inherited
    # from the original GPU database or an unrelated local installation.
    manifest["model_id"] = "admet-ai-2.0.1-local-" + identity_digest(manifest)[:20]
    mid = manifest["model_id"]
    current = root / "models/admet-ai-2.0.1"
    current_manifest = current / "manifest.json"
    if (
        current_manifest.exists()
        and json.loads(current_manifest.read_bytes()) != manifest
    ):
        raise ValueError(
            "Saved catalog model/channel differs; refusing to mix prediction identities"
        )
    archive_model(
        root, manifest, Path(profile["endpoint_definitions"]["path"]).read_bytes()
    )
    current.mkdir(parents=True, exist_ok=True)
    atomic_json(current / "manifest.json", manifest)
    atomic_bytes(
        current / "endpoints.csv",
        Path(profile["endpoint_definitions"]["path"]).read_bytes(),
    )
    folder = root / "data/predictions" / mid
    folder.mkdir(parents=True, exist_ok=True)
    done = reused = computed = 0
    counts = Counter()
    started = time.monotonic()
    with duckdb.connect(str(root / "data/skinscout.duckdb"), read_only=True) as con:
        con.execute("SET threads=2; SET memory_limit='2GB'")
        total = con.execute("SELECT count(*) FROM compounds").fetchone()[0]
        with _input_batches(con, plan["batch_size"]) as batches:
            for inputs in batches:
                ih = digest(json.dumps(inputs, sort_keys=True))
                file = folder / (ih + ".parquet")
                receipt = file.with_suffix(".json")
                if file.exists() and receipt.exists():
                    result = validate_batch(root, file, expected_inputs=inputs)
                    counts.update(result["status_counts"])
                    reused += len(inputs)
                else:
                    # A killed pre-commit file is retained before a new attempt.
                    if file.exists() or receipt.exists():
                        orphan = root / "data/uncommitted" / uuid.uuid4().hex
                        orphan.mkdir(parents=True)
                        for p in [file, receipt]:
                            if p.exists():
                                p.replace(orphan / p.name)
                    eligible = [
                        r for r in inputs if r["inference_status"] == "eligible"
                    ]
                    predictions = []
                    # A GPU shard uses one Trainer invocation per ensemble; the
                    # dataloader remains bounded at 512 and OOM retries bisect the
                    # shard. CPU calls stay small for timely cancellation.
                    prediction_group_size = (
                        max(1, len(eligible))
                        if channel["device"] == "cuda"
                        else channel["inference_batch_size"]
                    )
                    for offset in range(0, len(eligible), prediction_group_size):
                        predictions.extend(
                            predict_batch(
                                model,
                                eligible[offset : offset + prediction_group_size],
                                endpoints,
                            )
                        )
                    predicted = iter(predictions)
                    output = []
                    for row in inputs:
                        status, values, error = (
                            next(predicted)
                            if row["inference_status"] == "eligible"
                            else (row["inference_status"], {}, None)
                        )
                        output.append(
                            {
                                "compound_id": row["compound_id"],
                                "model_id": mid,
                                "input_smiles": row["canonical_smiles"],
                                "input_hash": ih,
                                "computed_at": now(),
                                "status": status,
                                "error": error,
                                "values_json": json.dumps(values, allow_nan=False),
                                "quality_flags_json": json.dumps(
                                    {"exact_structure_only": True}
                                ),
                                "domain_status": "NOT_ESTABLISHED",
                            }
                        )
                    temp = file.with_suffix(".parquet.part")
                    pq.write_table(
                        pa.Table.from_pylist(output, schema=PREDICTION_SCHEMA),
                        temp,
                        compression="zstd",
                    )
                    with temp.open("rb") as stream:
                        os.fsync(stream.fileno())
                    temp.replace(file)
                    atomic_json(
                        receipt,
                        {
                            "rows": len(output),
                            "model_id": mid,
                            "input_hash": ih,
                            "sha256": sha256_file(file),
                        },
                    )
                    result = validate_batch(root, file, expected_inputs=inputs)
                    counts.update(result["status_counts"])
                    computed += len(inputs)
                done += len(inputs)
                seconds = time.monotonic() - started
                progress = {
                    "total": total,
                    "completed": done,
                    "newly_computed": computed,
                    "reused": reused,
                    "percent": round(done * 100 / total, 2),
                    "status_counts": dict(counts),
                    "execution_device": channel["device"],
                    "execution_channel": channel,
                    "elapsed_seconds": round(seconds, 1),
                    "remaining_seconds_estimate": round(
                        (total - done) * seconds / computed
                    )
                    if computed
                    else None,
                }
                atomic_json(root / "reports/admet-local-progress.json", progress)
                print(json.dumps(progress), flush=True)
    atomic_json(
        root / "reports/admet-local.json",
        {
            **progress,
            "model_id": mid,
            "status": "complete",
            "scientific_validation": "NOT_ESTABLISHED",
        },
    )


def indexes(root):
    from .catalog_index import build_index, build_filtered_index
    from .scout_predictions import build_cache, PublishedPredictionCache
    from .admet_search import build_search_gate
    from .scout_gateway import Gateway

    # Incomplete index output remains private and is never registered. Preserve
    # it for diagnosis; build a fresh generation on resume without new inference.
    serving = root / "serving"
    if serving.exists():
        serving.replace(root / ("incomplete-serving-" + uuid.uuid4().hex))
    serving.mkdir()
    build_cache(root, serving / "predictions")
    build_index(root, serving / "similarity")
    build_search_gate(root, serving / "predictions", serving / "admet-gate")
    config = {
        "kind": "local",
        "scout_root": str(root),
        "published_prediction_cache": str(serving / "predictions"),
        "similarity_index": str(serving / "similarity"),
        "admet_search_gate": str(serving / "admet-gate"),
        "filtered_similarity_indices": {},
    }
    for hide, purchase, key in [
        (True, False, "risk"),
        (False, True, "purchase"),
        (True, True, "risk_purchase"),
    ]:
        dest = serving / key
        build_filtered_index(
            serving / "similarity",
            serving / "admet-gate",
            dest,
            hide_severe_admet=hide,
            purchasable_only=purchase,
        )
        config["filtered_similarity_indices"][key] = str(dest)
    gateway = Gateway({**config, "token": "x" * 32})
    try:
        if gateway.admission_errors:
            raise ValueError(
                "Catalog admission failed: " + json.dumps(gateway.admission_errors)
            )
        PublishedPredictionCache(serving / "predictions").verify()
    finally:
        gateway.executor.shutdown(wait=True, cancel_futures=True)
    atomic_json(root / "local-catalog.json", config)


def main():
    root, stage = Path(sys.argv[1]).resolve(), sys.argv[2]
    plan = json.loads((root / "owner-plan.json").read_bytes())
    from .catalog_bootstrap import code_identity

    if code_identity() != plan["code"]:
        raise ValueError("Installed worker code differs from the frozen build plan")
    if stage == "download":
        download(root, plan)
    elif stage == "ingest":
        from .._catalog_source.ingest import ingest

        for spec in plan["sources"]:
            result = ingest(root, spec["id"], workers=2)
            if result["status"] != "complete":
                raise ValueError("Incomplete ingestion: " + spec["id"])
    elif stage == "database":
        from .._catalog_source.database import build

        result = build(root)
        if not result["admitted_structure_sources_complete"] or not result["compounds"]:
            raise ValueError("Source population incomplete or empty")
    elif stage == "admet":
        admet(root, plan)
    elif stage == "audit":
        from .._catalog_source.publish import publish
        from .._catalog_source.audit import audit

        publish(root)
        result = audit(root, verify_hashes=True)
        if not result["auxiliary_admet_db_complete"]:
            raise ValueError("Catalog audit failed: " + json.dumps(result["checks"]))
    elif stage == "indexes":
        indexes(root)
    else:
        raise ValueError("Unsupported stage")
    # Model libraries sometimes retain exit threads. Parent observes durable
    # completion then terminates and reaps only this owned process group.
    atomic_json(
        root / "reports" / (stage + "-worker-complete.json"),
        {"plan_sha256": identity_digest(plan), "stage": stage},
    )


if __name__ == "__main__":
    main()
