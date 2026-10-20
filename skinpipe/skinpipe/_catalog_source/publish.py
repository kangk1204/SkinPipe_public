"""Expose only committed, semantically validated predictions to live SQL readers."""

import json
import os
from pathlib import Path
import pyarrow as pa
import pyarrow.parquet as pq
from .common import atomic_json, digest, lock, now, sha256
from .model_registry import collect_models
from .prediction_integrity import validate_batch
from .prediction_schema import ensure_empty_relations, MODEL_SCHEMA, ENDPOINT_SCHEMA


def file_identity(path):
    st = path.stat()
    return [st.st_dev, st.st_ino, st.st_size, st.st_mtime_ns, st.st_ctime_ns]


def validation_key(root, file, receipt, catalog_identity):
    """Stat identity is only a reuse shortcut after full hash/semantic validation."""
    model = root / "models/registry" / file.parent.name
    metadata = [
        sha256(p) if p.exists() else None
        for p in [model / "manifest.json", model / "endpoints.csv"]
    ]
    validators = [
        Path(__file__),
        Path(__file__).with_name("prediction_integrity.py"),
        Path(__file__).with_name("model_registry.py"),
    ]
    return digest(
        json.dumps(
            [
                file_identity(file),
                sha256(receipt),
                metadata,
                catalog_identity,
                [sha256(p) for p in validators],
            ],
            sort_keys=True,
        )
    )


def atomic_link(file, dest, expected_identity):
    tmp = dest.with_name(dest.name + f".{os.getpid()}.part")
    try:
        if tmp.exists():
            tmp.unlink()
        os.link(file, tmp)
        # Creating a hard link changes ctime, but not inode/size/mtime.
        if file_identity(tmp)[:4] != expected_identity[:4]:
            raise ValueError("Prediction changed during publication")
        tmp.replace(dest)
    finally:
        if tmp.exists():
            tmp.unlink()


def catalog_input_mismatches(con, files):
    """One disk-spillable catalog join for all semantically valid shards.

    Return only per-file mismatch counts, never all 17M input rows in Python.
    The publisher calls this before linking any new batch into the served view.
    """
    if not files:
        return {}
    return dict(
        con.execute(
            """SELECT p.filename, count(*) FROM
        read_parquet(?,union_by_name=true,filename=true) p
        LEFT JOIN compounds c USING(compound_id)
        WHERE c.compound_id IS NULL OR p.input_smiles IS DISTINCT FROM c.canonical_smiles
            OR (CASE WHEN p.status IN ('success','failed','partial_nonfinite') THEN 'eligible'
                     ELSE p.status END) IS DISTINCT FROM c.inference_status
        GROUP BY p.filename""",
            [[str(file) for file in files]],
        ).fetchall()
    )


def publish(root):
    out = root / "data/published_predictions"
    out.mkdir(parents=True, exist_ok=True)
    added = 0
    replaced = 0
    removed = 0
    invalid = []
    accepted = set()

    def reject(item, reason):
        nonlocal removed
        # Keep raw data and its receipt. Only withdraw an invalid serving link.
        if item["dest"].exists():
            item["dest"].unlink()
            removed += 1
        detail = {
            "source": str(item["file"].relative_to(root)),
            "reason": str(reason),
            "status": "quarantined",
        }
        invalid.append(detail)
        atomic_json(
            item["ledger"],
            {
                **detail,
                "validation_key": item["key"],
                "checked_at": now(),
            },
        )

    with lock(root / "logs/publish.lock", blocking=True):
        ensure_empty_relations(root)
        models, endpoints, unresolved = collect_models(root)
        # No source result is exposed until input identity has a catalog to bind to.
        catalog = root / "data/skinscout.duckdb"
        con = None
        catalog_identity = file_identity(catalog) if catalog.exists() else None
        if catalog.exists():
            from .database import connect

            con = connect(root, True)
        try:
            prepared = []
            for receipt in sorted((root / "data/predictions").glob("*/*.json")):
                file = receipt.with_suffix(".parquet")
                dest = out / (file.parent.name + "--" + file.name)
                ledger = (
                    root / "manifests/published_predictions" / (dest.stem + ".json")
                )
                if not file.exists():
                    continue
                item = {
                    "file": file,
                    "receipt": receipt,
                    "dest": dest,
                    "ledger": ledger,
                    "key": None,
                }
                try:
                    if con is None:
                        raise ValueError(
                            "Pending publication: catalog input validation unavailable"
                        )
                    key = validation_key(root, file, receipt, catalog_identity)
                    item["key"] = key
                    previous = json.loads(ledger.read_text()) if ledger.exists() else {}
                    if (
                        previous.get("validation_key") == key
                        and previous.get("status") == "quarantined"
                    ):
                        raise ValueError(previous["reason"])
                    cached = (
                        previous.get("validation_key") == key
                        and previous.get("status") == "validated"
                        and dest.exists()
                        and os.path.samefile(file, dest)
                    )
                    validated_identity = file_identity(file)
                    result = previous if cached else validate_batch(root, file)
                    if file_identity(file) != validated_identity:
                        raise ValueError("Prediction changed during validation")
                    prepared.append(
                        {
                            **item,
                            "cached": cached,
                            "identity": validated_identity,
                            "result": result,
                        }
                    )
                except (ValueError, OSError, KeyError, TypeError) as exc:
                    reject(item, exc)
            mismatches = catalog_input_mismatches(
                con, [item["file"] for item in prepared]
            )
            if prepared and file_identity(catalog) != catalog_identity:
                raise ValueError(
                    "Catalog changed during input validation; publication held"
                )
            # Verify all identities again before exposing any new batch. A query
            # execution failure above leaves existing valid links untouched.
            for item in prepared:
                if (
                    file_identity(item["file"]) != item["identity"]
                    or validation_key(
                        root, item["file"], item["receipt"], catalog_identity
                    )
                    != item["key"]
                ):
                    raise ValueError(
                        "Prediction changed during catalog binding; publication held"
                    )
            for item in prepared:
                file, dest = item["file"], item["dest"]
                if bad := mismatches.get(str(file)):
                    reject(
                        item,
                        f"Prediction input identity/eligibility mismatch: {bad} rows",
                    )
                    continue
                if not item["cached"]:
                    if dest.exists():
                        if not os.path.samefile(file, dest):
                            atomic_link(file, dest, item["identity"])
                            replaced += 1
                    else:
                        atomic_link(file, dest, item["identity"])
                        added += 1
                    key = validation_key(root, file, item["receipt"], catalog_identity)
                    atomic_json(
                        item["ledger"],
                        {
                            **item["result"],
                            "status": "validated",
                            "validation_key": key,
                            "validated_at": now(),
                            "source": str(file.relative_to(root)),
                        },
                    )
                accepted.add(dest.name)
        finally:
            if con is not None:
                con.close()
        # Receipts removed during repair must not leave indefinitely visible results.
        for dest in out.glob("*.parquet"):
            if dest.name != "_empty.parquet" and dest.name not in accepted:
                dest.unlink()
                removed += 1
        for name, rows, schema in [
            ("models", models, MODEL_SCHEMA),
            ("endpoints", endpoints, ENDPOINT_SCHEMA),
        ]:
            folder = root / "data" / name
            folder.mkdir(parents=True, exist_ok=True)
            tmp = folder / "registry.part"
            pq.write_table(
                pa.Table.from_pylist(rows, schema=schema), tmp, compression="zstd"
            )
            tmp.replace(folder / "registry.parquet")
        result = {
            "published_at": now(),
            "new_batches": added,
            "replaced_batches": replaced,
            "removed_links": removed,
            "committed_batches": len(accepted),
            "quarantined_batches": invalid,
            "unresolved_models": unresolved,
            "validation": "committed bytes, immutable model metadata, output schema, input hash and catalog identity",
        }
        atomic_json(root / "reports/published_predictions.json", result)
        return result
