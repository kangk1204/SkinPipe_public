"""Immutable model generations; absent historical provenance is never invented."""

import csv
import io
import json
from pathlib import Path

from .common import atomic_json, sha256
from .policy import load_policy


class UnknownModelError(ValueError):
    pass


def archive_model(root, spec, endpoint_bytes):
    mid = spec["model_id"]
    if Path(mid).name != mid or mid in [".", ".."]:
        raise ValueError("Invalid model ID")
    expected = spec.get("endpoint_metadata_sha256")
    if expected:
        import hashlib

        if hashlib.sha256(endpoint_bytes).hexdigest() != expected:
            raise ValueError("Endpoint metadata hash mismatch for " + mid)
    out = root / "models/registry" / mid
    out.mkdir(parents=True, exist_ok=True)
    manifest = out / "manifest.json"
    endpoints = out / "endpoints.csv"
    if manifest.exists() and json.loads(manifest.read_text()) != spec:
        raise ValueError("Immutable model manifest conflict for " + mid)
    if endpoints.exists() and endpoints.read_bytes() != endpoint_bytes:
        raise ValueError("Immutable endpoint metadata conflict for " + mid)
    if not endpoints.exists():
        tmp = out / "endpoints.csv.part"
        tmp.write_bytes(endpoint_bytes)
        tmp.replace(endpoints)
    if not manifest.exists():
        atomic_json(manifest, spec)
    return out


def snapshot_models(root):
    for file in sorted((root / "models").glob("*/manifest.json")):
        spec = json.loads(file.read_text())
        if "model_id" not in spec:
            continue
        endpointfile = file.parent / "endpoints.csv"
        if endpointfile.exists():
            archive_model(root, spec, endpointfile.read_bytes())


def model_metadata(root, mid):
    folder = root / "models/registry" / mid
    manifest = folder / "manifest.json"
    endpointfile = folder / "endpoints.csv"
    if not manifest.exists() or not endpointfile.exists():
        raise UnknownModelError("Missing immutable model/endpoint provenance: " + mid)
    spec = json.loads(manifest.read_text())
    if spec.get("model_id") != mid:
        raise ValueError("Model manifest identity mismatch")
    if (
        spec.get("endpoint_metadata_sha256")
        and sha256(endpointfile) != spec["endpoint_metadata_sha256"]
    ):
        raise ValueError("Changed immutable endpoint metadata: " + mid)
    rows = list(csv.DictReader(io.StringIO(endpointfile.read_text())))
    ids = [r["id"] for r in rows]
    if not rows or len(set(ids)) != len(ids):
        raise ValueError("Missing/duplicate endpoint definitions: " + mid)
    if spec.get("name") == "ADMET-AI" and len(ids) != 52:
        raise ValueError("Expected 52 ADMET and physicochemical endpoint definitions")
    return spec, rows


def collect_models(root):
    snapshot_models(root)
    models = []
    endpoints = []
    known = set()
    for file in sorted((root / "models/registry").glob("*/manifest.json")):
        mid = file.parent.name
        spec, rows = model_metadata(root, mid)
        known.add(mid)
        models.append(
            dict(
                model_id=mid,
                metadata_json=json.dumps(
                    {
                        **spec,
                        "metadata_status": "verified",
                        "use_policy": load_policy(root),
                    }
                ),
            )
        )
        for row in rows:
            endpoints.append(
                dict(
                    model_id=mid,
                    endpoint=row["id"],
                    name=row["name"],
                    unit=row["units"],
                    category=row["category"],
                    task_type=row["task_type"],
                    species=row["species"],
                    source_url=row["url"],
                    published_metrics_not_local_validation=json.dumps(
                        {k: row.get(k) for k in ["AUPRC", "AUROC", "R^2", "MAE"]}
                    ),
                )
            )
    observed = (
        {p.name for p in (root / "data/predictions").iterdir() if p.is_dir()}
        if (root / "data/predictions").exists()
        else set()
    )
    unresolved = [
        {
            "model_id": mid,
            "metadata_status": "unresolved",
            "reason": "Historical manifest and endpoint provenance unavailable; batches quarantined",
        }
        for mid in sorted(observed - known)
    ]
    # Explicit records describe the missing provenance without borrowing another model's metadata.
    models.extend(
        dict(model_id=r["model_id"], metadata_json=json.dumps(r)) for r in unresolved
    )
    return models, endpoints, unresolved
