"""Coverage audit validates committed semantics and the SQL-visible catalog."""

import json
from pathlib import Path
import pyarrow as pa
from .common import atomic_json, now, sha256
from .database import connect, quote
from .policy import load_policy, auxiliary_complete
from .model_registry import snapshot_models, model_metadata
from .prediction_integrity import validate_batch

STRUCTURE_KINDS = {
    "smiles",
    "chembl_tsv",
    "coconut_csv",
    "lotus_metadata",
    "npatlas_tsv",
    "sdf",
    "ice_xlsx",
    "huskin_csv",
    "skinpix_xlsx",
    "normalized_jsonl",
}


def record_execution_error(result, check, error):
    """A query that did not finish is not evidence of a semantic mismatch."""
    result["checks"][check] = None
    result["execution_errors"].append(
        {
            "check": check,
            "exception_type": type(error).__name__,
            "message": str(error),
            "recorded_at": now(),
        }
    )


def checked_query(result, con, check, sql, parameters=None):
    try:
        result["checks"][check] = con.execute(sql, parameters or []).fetchone()[0]
    except Exception as error:
        record_execution_error(result, check, error)


def published_endpoint_check(con, batch_size=4096, progress=None):
    """Examine every served endpoint without expanding predictions into 52x rows.

    DuckDB projects one key list/type per prediction; Arrow and Python retain only a
    bounded batch. Membership is exact, including duplicate-key occurrences and
    unknown models for JSON objects. Other non-SQL-NULL JSON shapes are counted
    separately as invalid. Values are not cached or sampled. Definitions are small.
    """
    if batch_size < 1:
        raise ValueError("Endpoint audit batch size must be positive")
    definitions = {}
    for mid, endpoint in con.execute(
        "SELECT model_id,endpoint FROM endpoints"
    ).fetchall():
        # SQL equality never matches NULL, including NULL against NULL.
        if mid is not None and endpoint is not None:
            definitions.setdefault(mid, set()).add(endpoint)
    reader = con.execute(
        "SELECT model_id,json_keys(values_json) AS endpoint_keys,json_type(values_json) AS value_type FROM predictions"
    ).fetch_record_batch(rows_per_batch=batch_size)
    result = {
        "missing_definitions": 0,
        "invalid_value_objects": 0,
        "prediction_rows_examined": 0,
        "endpoint_values_examined": 0,
        "batches": 0,
        "max_batch_rows": 0,
        "requested_batch_rows": batch_size,
        "method": "exact served-SQL endpoint keys in bounded Arrow batches",
    }
    empty = frozenset()
    try:
        for batch in reader:
            result["batches"] += 1
            result["max_batch_rows"] = max(result["max_batch_rows"], batch.num_rows)
            result["prediction_rows_examined"] += batch.num_rows
            for mid, keys, kind in zip(
                batch.column(0).to_pylist(),
                batch.column(1).to_pylist(),
                batch.column(2).to_pylist(),
            ):
                result["invalid_value_objects"] += int(
                    kind is not None and kind != "OBJECT"
                )
                keys = keys or []
                expected = definitions.get(mid, empty)
                result["endpoint_values_examined"] += len(keys)
                result["missing_definitions"] += sum(
                    key not in expected for key in keys
                )
            if progress is not None and result["batches"] % 25 == 0:
                progress(dict(result))
    finally:
        reader.close()
    if progress is not None:
        progress(dict(result))
    return result


def catalog_checks(root, con, checks, result=None):
    from .ingest import expected_runhash

    config = root / "config/sources.json"
    build = root / "reports/db_build.json"
    checks["missing_or_stale_build_report"] = 0
    checks["incomplete_admitted_structure_sources"] = 0
    if not config.exists() or not build.exists():
        checks["missing_or_stale_build_report"] = 1
        return
    specs = json.loads(config.read_text())
    report = json.loads(build.read_text())
    if not report.get("admitted_structure_sources_complete"):
        checks["incomplete_admitted_structure_sources"] += 1
    actual = con.execute("SELECT count(*) FROM compounds").fetchone()[0]
    if report.get("compounds") != actual:
        checks["missing_or_stale_build_report"] += 1
    try:
        stored = {
            sid: json.loads(metadata)
            for sid, metadata in con.execute(
                "SELECT source_id,metadata FROM sources"
            ).fetchall()
        }
        if stored != {s["id"]: s for s in specs}:
            checks["missing_or_stale_build_report"] += 1
    except Exception as error:
        if result is None:
            raise
        record_execution_error(result, "missing_or_stale_build_report", error)
    for spec in specs:
        if (
            spec.get("admission") != "allowed"
            or spec.get("kind") not in STRUCTURE_KINDS
        ):
            continue
        sid = spec["id"]
        ingest = root / "manifests/ingest" / (sid + ".json")
        download = root / "manifests/downloads" / (sid + ".json")
        if not ingest.exists() or not download.exists():
            checks["incomplete_admitted_structure_sources"] += 1
            continue
        state = json.loads(ingest.read_text())
        source = json.loads(download.read_text())
        valid = (
            sid in report.get("included_ingestions", [])
            and state.get("status") == "complete"
            and source.get("status") == "verified"
            and state.get("input_sha256") == source.get("sha256")
            and state.get("runhash") == expected_runhash(spec, source.get("sha256"))
        )
        if not valid:
            checks["incomplete_admitted_structure_sources"] += 1


def audit(root, verify_hashes=False):
    """Run a fresh audit; the progress file is observational, never a resume cache."""
    started = now()

    def progress(stage, **state):
        atomic_json(
            root / "reports/audit_progress.json",
            {
                "started_at": started,
                "updated_at": now(),
                "stage": stage,
                "status": "running",
                **state,
            },
        )

    progress("starting", verify_hashes=verify_hashes)
    try:
        result = _audit(root, verify_hashes, progress)
    except Exception as error:
        progress(
            "execution_failed",
            status="execution_failed",
            exception_type=type(error).__name__,
            message=str(error),
        )
        raise
    progress(
        "finished",
        status="complete"
        if result["complete"]
        else "execution_failed"
        if result["execution_errors"]
        else "needs_attention",
        complete=result["complete"],
        audit_execution_passed=result["audit_execution_passed"],
        execution_error_count=len(result["execution_errors"]),
    )
    return result


def _audit(root, verify_hashes, progress):
    policy = load_policy(root)
    result = {
        "audited_at": now(),
        "database": str(root / "data/skinscout.duckdb"),
        "checks": {},
        "use_policy": policy,
        "hash_verification_performed": verify_hashes,
        "prediction_validation_errors": [],
        "execution_errors": [],
        "audit_code_sha256": sha256(Path(__file__)),
    }
    checks = result["checks"]
    with connect(root, True) as con:
        progress("catalog_checks")
        catalog_checks(root, con, checks, result)
        if verify_hashes:
            progress("ingestion_hash_validation")
            checks["bad_ingest_hashes"] = 0
            for receipt in (root / "manifests/ingest").glob("*.json"):
                state = json.loads(receipt.read_text())
                if state.get("status") != "complete":
                    continue
                for name in state["shards"]:
                    file = root / name
                    meta = file.with_suffix(".json")
                    checks["bad_ingest_hashes"] += int(
                        not file.exists()
                        or not meta.exists()
                        or sha256(file) != json.loads(meta.read_text())["sha256"]
                    )
        checks["duplicate_compound_ids"] = con.execute(
            "SELECT count(*)-count(DISTINCT compound_id) FROM compounds"
        ).fetchone()[0]
        checks["orphan_structure_links"] = con.execute(
            "SELECT count(*) FROM source_records r ANTI JOIN compounds c USING(compound_id) WHERE r.compound_id IS NOT NULL"
        ).fetchone()[0]
        checks["duplicate_source_records"] = con.execute(
            "SELECT count(*)-count(DISTINCT record_id) FROM source_records"
        ).fetchone()[0]
        result["compound_count"] = con.execute(
            "SELECT count(*) FROM compounds"
        ).fetchone()[0]
        result["source_counts"] = (
            con.execute("SELECT * FROM source_counts").fetchdf().to_dict("records")
        )
        result["unsupported_by_reason"] = dict(
            con.execute(
                "SELECT inference_status,count(*) FROM compounds WHERE inference_status<>'eligible' GROUP BY inference_status"
            ).fetchall()
        )
        snapshot_models(root)
        progress("model_metadata")
        model = root / "models/admet-ai-2.0.1/manifest.json"
        checks["invalid_prediction_batches"] = 0
        checks["missing_required_model_metadata"] = 0
        checks["sql_prediction_mismatch"] = 0
        result["prediction_missing"] = result["compound_count"]
        result["prediction_counts"] = {}
        if model.exists():
            mid = json.loads(model.read_text())["model_id"]
            folder = root / "data/predictions" / mid
            try:
                model_spec, definitions = model_metadata(root, mid)
                expected = [
                    dict(
                        model_id=mid,
                        endpoint=r["id"],
                        unit=r["units"],
                        name=r["name"],
                        category=r["category"],
                        task_type=r["task_type"],
                        species=r["species"],
                        source_url=r["url"],
                    )
                    for r in definitions
                ]
                con.register(
                    "expected_current_endpoints", pa.Table.from_pylist(expected)
                )
                checks["served_endpoint_metadata_mismatch"] = con.execute(
                    """SELECT count(*) FROM expected_current_endpoints x
                    FULL JOIN (SELECT * FROM endpoints WHERE model_id=?) e USING(model_id,endpoint)
                    WHERE x.endpoint IS NULL OR e.endpoint IS NULL OR x.unit IS DISTINCT FROM e.unit
                       OR x.name IS DISTINCT FROM e.name OR x.category IS DISTINCT FROM e.category
                       OR x.task_type IS DISTINCT FROM e.task_type OR x.species IS DISTINCT FROM e.species
                       OR x.source_url IS DISTINCT FROM e.source_url""",
                    [mid],
                ).fetchone()[0]
                checks["duplicate_endpoint_definitions"] = con.execute(
                    "SELECT count(*)-count(DISTINCT (model_id,endpoint)) FROM endpoints"
                ).fetchone()[0]
                served = con.execute(
                    "SELECT metadata_json FROM models WHERE model_id=?", [mid]
                ).fetchall()
                checks["served_model_metadata_mismatch"] = int(
                    len(served) != 1
                    or {
                        k: v
                        for k, v in json.loads(served[0][0]).items()
                        if k not in ["metadata_status", "use_policy"]
                    }
                    != model_spec
                )
                con.unregister("expected_current_endpoints")
            except (ValueError, KeyError) as error:
                checks["missing_required_model_metadata"] += 1
                result["served_metadata_error"] = str(error)
                result["served_metadata_error_type"] = type(error).__name__
            except Exception as error:
                record_execution_error(result, "served_metadata_checks", error)
                result["served_metadata_error"] = str(error)
            valid = []
            files = [
                file
                for file in sorted(folder.glob("*.parquet"))
                if file.with_suffix(".json").exists()
            ]
            progress(
                "prediction_batch_validation",
                model_id=mid,
                checked_batches=0,
                total_batches=len(files),
            )
            for index, file in enumerate(files, 1):
                try:
                    # A completion audit always validates semantics and hashes. The option
                    # additionally enables original ingestion hash validation above.
                    # Catalog binding is performed once below, preserving the
                    # exact validate_inputs predicate without 1,708 catalog scans.
                    validate_batch(root, file)
                    valid.append(file)
                except (ValueError, KeyError, TypeError, OSError) as error:
                    checks["invalid_prediction_batches"] += 1
                    result["prediction_validation_errors"].append(
                        {"file": str(file.relative_to(root)), "error": str(error)}
                    )
                if index % 25 == 0 or index == len(files):
                    progress(
                        "prediction_batch_validation",
                        model_id=mid,
                        checked_batches=index,
                        total_batches=len(files),
                        valid_batches=len(valid),
                        invalid_batches=checks["invalid_prediction_batches"],
                    )
            if valid:
                progress(
                    "prediction_counts_and_coverage",
                    model_id=mid,
                    validated_batches=len(valid),
                )
                plist = "[" + ",".join(quote(p) for p in valid) + "]"
                con.execute(
                    f"CREATE TEMP VIEW current_predictions AS SELECT * FROM read_parquet({plist})"
                )
                checks["duplicate_prediction_rows"] = con.execute(
                    "SELECT count(*)-count(DISTINCT compound_id) FROM current_predictions"
                ).fetchone()[0]
                result["prediction_counts"] = dict(
                    con.execute(
                        "SELECT status,count(*) FROM current_predictions GROUP BY status"
                    ).fetchall()
                )
                result["prediction_missing"] = con.execute(
                    "SELECT count(*) FROM compounds ANTI JOIN current_predictions USING(compound_id)"
                ).fetchone()[0]
                checks["forbidden_drugbank_outputs"] = con.execute(
                    "SELECT count(*) FROM current_predictions WHERE contains(lower(values_json),'drugbank') OR contains(lower(values_json),'percentile')"
                ).fetchone()[0]
                progress("prediction_catalog_binding", model_id=mid)
                checked_query(
                    result,
                    con,
                    "invalid_prediction_inputs",
                    """SELECT count(*) FROM current_predictions p
                    LEFT JOIN compounds c USING(compound_id)
                    WHERE c.compound_id IS NULL OR p.input_smiles IS DISTINCT FROM c.canonical_smiles
                        OR (CASE WHEN p.status IN ('success','failed','partial_nonfinite') THEN 'eligible'
                                 ELSE p.status END) IS DISTINCT FROM c.inference_status""",
                )
                progress("sql_prediction_parity", model_id=mid)
                checked_query(
                    result,
                    con,
                    "sql_prediction_mismatch",
                    """SELECT count(*) FROM current_predictions c
                        FULL JOIN (SELECT * FROM predictions WHERE model_id=?) p USING(compound_id)
                        WHERE c.compound_id IS NULL OR p.compound_id IS NULL
                           OR c.model_id IS DISTINCT FROM p.model_id OR c.input_smiles IS DISTINCT FROM p.input_smiles
                           OR c.status IS DISTINCT FROM p.status OR c.values_json IS DISTINCT FROM p.values_json
                           OR c.input_hash IS DISTINCT FROM p.input_hash OR c.computed_at IS DISTINCT FROM p.computed_at
                           OR c.error IS DISTINCT FROM p.error OR c.quality_flags_json IS DISTINCT FROM p.quality_flags_json
                           OR c.domain_status IS DISTINCT FROM p.domain_status""",
                    [mid],
                )
            else:
                checked_query(
                    result,
                    con,
                    "sql_prediction_mismatch",
                    "SELECT count(*) FROM predictions WHERE model_id=?",
                    [mid],
                )
        else:
            checks["missing_required_model_metadata"] += 1
        result["prediction_coverage"] = (
            1 - result["prediction_missing"] / result["compound_count"]
            if result["compound_count"]
            else 0.0
        )
        progress("published_model_and_duplicate_checks")
        checked_query(
            result,
            con,
            "published_orphan_models",
            "SELECT count(*) FROM predictions p ANTI JOIN models m USING(model_id)",
        )
        checked_query(
            result,
            con,
            "duplicate_published_predictions",
            "SELECT count(*)-count(DISTINCT (compound_id,model_id)) FROM predictions",
        )
        try:
            progress("published_endpoint_definitions", prediction_rows_examined=0)
            scan = published_endpoint_check(
                con,
                progress=lambda state: progress(
                    "published_endpoint_definitions", **state
                ),
            )
            checks["published_missing_endpoint_definitions"] = scan[
                "missing_definitions"
            ]
            checks["published_invalid_value_objects"] = scan["invalid_value_objects"]
            result["published_endpoint_scan"] = scan
        except Exception as error:
            record_execution_error(
                result, "published_missing_endpoint_definitions", error
            )
            checks["published_invalid_value_objects"] = None
        result["base_prediction_complete"] = (
            result["compound_count"] > 0 and result["prediction_missing"] == 0
        )
        result["unresolved_source_records"] = con.execute(
            "SELECT count(*) FROM source_records WHERE compound_id IS NULL"
        ).fetchone()[0]
        result["unresolved_policy"] = (
            "Retained in ingredients_unresolved with structure_status; no artificial single molecule or safe value assigned"
        )
        checks["unresolved_without_explicit_status"] = con.execute(
            "SELECT count(*) FROM source_records WHERE compound_id IS NULL AND (structure_status IS NULL OR structure_status='valid')"
        ).fetchone()[0]
    result["audit_execution_passed"] = not result["execution_errors"]
    result["software_integrity_passed"] = result["audit_execution_passed"] and all(
        v == 0 for v in checks.values()
    )
    result["sota_status"] = "NOT_ESTABLISHED"
    result["auxiliary_admet_db_complete"] = auxiliary_complete(result, policy)
    result["complete"] = result["auxiliary_admet_db_complete"]
    result["completion_scope"] = (
        "Auxiliary ADMET for the current admitted structure catalog; committed semantics, SQL visibility and source snapshot completeness are required. Scientific validation is separate."
    )
    result["failed_prediction_rows"] = sum(
        result["prediction_counts"].get(k, 0) for k in ["failed", "partial_nonfinite"]
    )
    atomic_json(root / "reports/audit.json", result)
    return result
