"""Validate committed prediction bytes and their model/input/output semantics."""

import json
import math
from collections import Counter
import uuid

import pyarrow as pa
import pyarrow.parquet as pq

from .common import atomic_json, digest, lock, now, sha256
from .model_registry import model_metadata, UnknownModelError

UNSUPPORTED = {
    "unsupported_query_or_polymer",
    "unsupported_multi_component",
    "unsupported_element",
    "outside_operational_size_limit",
    "unsupported_enhanced_stereo",
}
TERMINAL = {"success", "failed", "partial_nonfinite"} | UNSUPPORTED


def quarantine_batch(root, file, reason):
    """Preserve invalid bytes/receipt before a resumable job replaces that batch."""
    receipt = file.with_suffix(".json")
    fingerprint = sha256(file)
    target = root / "data/quarantined_predictions" / file.parent.name / fingerprint
    target.mkdir(parents=True, exist_ok=True)
    with lock(root / "logs/publish.lock", blocking=True):
        # A unique event folder preserves receipts even when the same raw bytes
        # have appeared before with a different or damaged receipt.
        event = target / uuid.uuid4().hex
        event.mkdir()
        original = str(file.relative_to(root))
        file.replace(event / file.name)
        if receipt.exists():
            receipt.replace(event / receipt.name)
        link = (
            root / "data/published_predictions" / (file.parent.name + "--" + file.name)
        )
        if link.exists():
            link.unlink()
        state = {
            "original_path": original,
            "archived_path": str((event / file.name).relative_to(root)),
            "sha256": fingerprint,
            "reason": str(reason),
            "quarantined_at": now(),
        }
        atomic_json(event / "quarantine.json", state)
    return state


def validate_values(values, status, expected_endpoints):
    if not isinstance(values, dict):
        raise ValueError("Prediction values must be an object")
    if status in ["success", "partial_nonfinite"]:
        if set(values) != set(expected_endpoints):
            raise ValueError("Prediction endpoint schema mismatch")
        finite = []
        for value in values.values():
            if value is None:
                finite.append(False)
                continue
            if (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(value)
            ):
                raise ValueError("Invalid prediction quantity")
            finite.append(True)
        if status == "success" and not all(finite):
            raise ValueError("Success contains missing endpoint values")
        if status == "partial_nonfinite" and all(finite):
            raise ValueError("Partial status without a missing value")
    elif values:
        raise ValueError("Failed/unsupported row contains prediction values")


def validate_inputs(con, rows):
    """One set-wise join validates identity and eligibility, including enhanced stereo."""
    name = "prediction_inputs_" + uuid.uuid4().hex
    table = pa.Table.from_pylist(
        [
            {
                "compound_id": r["compound_id"],
                "input_smiles": r["input_smiles"],
                "inference_status": "eligible"
                if r["status"] in ["success", "failed", "partial_nonfinite"]
                else r["status"],
            }
            for r in rows
        ]
    )
    con.register(name, table)
    try:
        bad = con.execute(f"""SELECT count(*) FROM {name} p LEFT JOIN compounds c USING(compound_id)
            WHERE c.compound_id IS NULL OR p.input_smiles IS DISTINCT FROM c.canonical_smiles
              OR p.inference_status IS DISTINCT FROM c.inference_status""").fetchone()[
            0
        ]
        if bad:
            raise ValueError(
                f"Prediction input identity/eligibility mismatch: {bad} rows"
            )
    finally:
        con.unregister(name)


def validate_batch(
    root, file, receipt=None, *, require_metadata=True, con=None, expected_inputs=None
):
    receipt = (
        receipt
        if receipt is not None
        else json.loads(file.with_suffix(".json").read_text())
    )
    actual_sha = sha256(file)
    if actual_sha != receipt.get("sha256"):
        raise ValueError("Prediction checksum mismatch: " + str(file))
    mid = file.parent.name
    if receipt.get("model_id") != mid:
        raise ValueError("Prediction receipt model mismatch")
    try:
        spec, definitions = model_metadata(root, mid)
        endpoints = {r["id"] for r in definitions}
        metadata_status = "verified"
    except UnknownModelError:
        if require_metadata:
            raise
        endpoints = None
        metadata_status = "unresolved"
    rows = pq.read_table(file).to_pylist()
    if not rows or len(rows) != receipt.get("rows"):
        raise ValueError("Prediction receipt row count mismatch")
    ids = set()
    inputs = []
    required = {
        "compound_id",
        "model_id",
        "input_smiles",
        "input_hash",
        "computed_at",
        "status",
        "error",
        "values_json",
        "quality_flags_json",
        "domain_status",
    }
    for row in rows:
        if not required <= row.keys():
            raise ValueError("Missing prediction provenance columns")
        if not row["compound_id"] or row["compound_id"] in ids:
            raise ValueError("Missing/duplicate prediction compound ID")
        ids.add(row["compound_id"])
        if row["model_id"] != mid:
            raise ValueError("Prediction row model mismatch")
        if row["status"] not in TERMINAL:
            raise ValueError("Nonterminal/unknown prediction status")
        if not isinstance(row["input_smiles"], str) or not row["input_smiles"]:
            raise ValueError("Missing prediction input")
        if not row["computed_at"] or not row["domain_status"]:
            raise ValueError("Missing prediction provenance value")
        if row["status"] == "failed" and not row["error"]:
            raise ValueError("Failed row lacks explicit failure reason")
        values = json.loads(row["values_json"])
        json.loads(row["quality_flags_json"])
        if endpoints is not None:
            validate_values(values, row["status"], endpoints)
        elif not isinstance(values, dict):
            raise ValueError("Prediction values must be an object")
        inputs.append(
            {
                "compound_id": row["compound_id"],
                "canonical_smiles": row["input_smiles"],
                "inference_status": "eligible"
                if row["status"] in ["success", "failed", "partial_nonfinite"]
                else row["status"],
            }
        )
    input_hash = digest(json.dumps(inputs, sort_keys=True))
    if (
        receipt.get("input_hash") != input_hash
        or any(r["input_hash"] != input_hash for r in rows)
        or file.stem != input_hash
    ):
        raise ValueError("Prediction batch input hash mismatch")
    if expected_inputs is not None and inputs != expected_inputs:
        # The builder already streams these exact rows from its read-only
        # catalog cursor. Comparing them directly avoids a 17M-row join for each
        # 10k shard while still checking order, identity, SMILES and eligibility.
        raise ValueError(
            "Prediction input identity/eligibility differs from catalog batch"
        )
    if con is not None:
        validate_inputs(con, rows)
    return {
        "rows": len(rows),
        "model_id": mid,
        "input_hash": input_hash,
        "sha256": actual_sha,
        "status_counts": dict(Counter(r["status"] for r in rows)),
        "metadata_status": metadata_status,
    }
