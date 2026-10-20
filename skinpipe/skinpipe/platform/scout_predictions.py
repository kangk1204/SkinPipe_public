"""Immutable indexed reuse of Scout's audited, published ADMET predictions.

The donor catalog is read-only. A separate serving copy keeps bounded lookups
from scanning every Parquet batch and never executes the prediction model.
"""

from collections import Counter
from copy import deepcopy
from pathlib import Path
import csv
import hashlib
import io
import json
import math
import shutil
import time

from .artifacts import atomic_bytes, canonical_json, sha256_file, ArtifactIntegrityError
from .identity import digest

TERMINAL = {
    "success",
    "partial_nonfinite",
    "failed",
    "unsupported_query_or_polymer",
    "unsupported_multi_component",
    "unsupported_element",
    "outside_operational_size_limit",
    "unsupported_enhanced_stereo",
}


def pin_prediction_metadata(prediction, artifacts, key, *, pop_inline=False):
    """Verify donor JSON semantics, preserving its digest when stored in CAS.

    Partition identity JSON escapes Unicode; cache/platform JSON uses UTF-8.
    These are two distinct byte contracts even for the same semantic metadata.
    """
    field = key + '_sha256'
    value = prediction[key]
    utf8 = canonical_json(value)
    ascii_raw = json.dumps(value, sort_keys=True, separators=(',', ':'),
        ensure_ascii=True, allow_nan=False).encode()
    utf8_sha, ascii_sha = hashlib.sha256(utf8).hexdigest(), hashlib.sha256(ascii_raw).hexdigest()
    supplied = prediction.get(field)
    refs = prediction.setdefault('metadata_refs', {})
    reference = refs.get(field)
    if reference is not None:
        if (not isinstance(reference, dict) or reference.get('sha256') != supplied
                or reference.get('storage') not in {'external_registry', 'platform_cas'}
                or reference.get('source') != 'SkinScout model registry'):
            raise ArtifactIntegrityError('ADMET metadata reference/content mismatch')
        encoding = reference.get('digest_format')
        expected = (ascii_sha if encoding == 'canonical_json_ascii_v1' else
            utf8_sha if encoding in {'canonical_json_utf8_v1', 'canonical_json'} else None)
        if expected is None or supplied != expected:
            raise ArtifactIntegrityError('ADMET metadata encoding/content mismatch')
    elif supplied is not None:
        if supplied not in {utf8_sha, ascii_sha}:
            raise ArtifactIntegrityError('ADMET shared metadata content mismatch')
        reference = {'storage': 'external_registry', 'sha256': supplied,
            'source': 'SkinScout model registry',
            'digest_format': 'canonical_json_utf8_v1' if supplied == utf8_sha else 'canonical_json_ascii_v1'}
    if reference is not None and reference.get('storage') == 'external_registry':
        preserved = prediction.setdefault('donor_metadata_refs', {})
        if field in preserved and preserved[field] != reference:
            raise ArtifactIntegrityError('Immutable donor metadata provenance changed')
        preserved[field] = deepcopy(reference)
    artifact = artifacts.put_bytes(utf8)
    prediction[field] = artifact['sha256']
    refs[field] = {'storage': 'platform_cas', 'sha256': artifact['sha256'],
        'source': 'SkinScout model registry', 'digest_format': 'canonical_json_utf8_v1'}
    if pop_inline:
        prediction.pop(key)
    return artifact


def _finite_json_number(value):
    """A JSON number that is finite; an oversized integer is invalid, not a crash.

    ``int`` values beyond the float range make ``math.isfinite`` raise
    ``OverflowError``. An unreadable value must fail the admission contract the
    same way ``NaN`` does, or one corrupt row turns the whole cache scan into
    an unexplained crash.
    """
    if type(value) not in (int, float):
        return False
    try:
        return math.isfinite(value)
    except OverflowError:
        return False


def normalize_prediction_values(row, endpoints):
    """One numeric contract for partition and indexed-cache admission."""
    try:
        values = json.loads(row["values_json"])
    except (KeyError, TypeError, ValueError) as error:
        raise ArtifactIntegrityError("Invalid ADMET values JSON") from error
    expected = {entry["id"] for entry in endpoints}
    status = row.get("status")
    if not isinstance(values, dict) or status not in TERMINAL:
        raise ArtifactIntegrityError("Invalid ADMET value/status shape")
    if status in {"success", "partial_nonfinite"}:
        finite = [_finite_json_number(v) for v in values.values()]
        if (
            not expected or set(values) != expected
            or any(v is not None and not valid for v, valid in zip(values.values(), finite))
            or (status == "success") != all(finite)
        ):
            raise ArtifactIntegrityError("Invalid cached ADMET endpoint values")
        for endpoint in endpoints:
            value = values[endpoint["id"]]
            if (endpoint.get("task_type") == "classification"
                    or endpoint.get("units", endpoint.get("unit")) == "probability"):
                if value is not None and not 0 <= value <= 1:
                    raise ArtifactIntegrityError("Cached ADMET classifier probability outside 0..1")
    elif values:
        raise ArtifactIntegrityError("Unsupported ADMET cannot contain numerical values")
    return values


class PublishedPredictionCache:
    def __init__(self, root):
        self.root = Path(root)
        self.manifest = json.loads((self.root / "manifest.json").read_bytes())
        if self.manifest.get("schema") != "scout-published-admet-cache-v1":
            raise ValueError("Unknown ADMET cache schema")
        self.database = self.root / "predictions.duckdb"
        self._state = None
        self.verify()

    def verify(self):
        s = self.database.stat()
        state = (s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns, s.st_ctime_ns)
        if state != self._state:
            if sha256_file(self.database) != self.manifest["database_sha256"]:
                raise ArtifactIntegrityError("ADMET serving cache changed")
            after = self.database.stat()
            if state != (
                after.st_dev,
                after.st_ino,
                after.st_size,
                after.st_mtime_ns,
                after.st_ctime_ns,
            ):
                raise ArtifactIntegrityError("ADMET cache changed during verification")
            self._state = state
        return state

    def query(self, ids):
        import duckdb
        from .catalog_query_lifetime import interruptible_query, query_checkpoint

        query_checkpoint()
        if (
            not isinstance(ids, list)
            or not 1 <= len(ids) <= 1000
            or len(set(ids)) != len(ids)
        ):
            raise ValueError("ADMET cache requires 1..1000 distinct IDs")
        before = self.verify()
        rows = []
        # A 1,000-ID ART lookup on the full serving DB exceeds 512 MiB.
        # Close each small read-only connection so index allocations cannot
        # accumulate across chunks. Preserve the same global result ordering.
        for offset in range(0, len(ids), 16):
            query_checkpoint()
            chunk = ids[offset : offset + 16]
            with duckdb.connect(str(self.database), read_only=True) as conn:
                conn.execute("SET threads=1; SET memory_limit='512MB'")
                with interruptible_query(conn):
                    try:
                        c = conn.execute(
                            "SELECT * FROM predictions WHERE compound_id IN ("
                            + ",".join("?" for _ in chunk)
                            + ") ORDER BY compound_id,model_id",
                            chunk,
                        )
                        rows.extend(dict(zip([x[0] for x in c.description], r)) for r in c.fetchall())
                    except Exception:
                        query_checkpoint()
                        raise
                query_checkpoint()
        rows.sort(key=lambda row: (row["compound_id"], row["model_id"]))
        if self.verify() != before:
            raise ArtifactIntegrityError("ADMET cache changed during lookup")
        seen = set()
        for row in rows:
            model = self.manifest["models"].get(row["model_id"])
            if model is None:
                raise ArtifactIntegrityError("Unknown cached ADMET model")
            values = normalize_prediction_values(row, model["endpoints"])
            row.update(
                values=values,
                quality_flags=json.loads(row["quality_flags_json"]),
                model_manifest_sha256=model["manifest_sha256"],
                endpoint_definitions_sha256=model["endpoints_sha256"],
                metadata_refs={
                    'model_manifest_sha256': {'storage': 'external_registry',
                        'sha256': model['manifest_sha256'], 'source': 'SkinScout model registry',
                        'digest_format': 'canonical_json_utf8_v1'},
                    'endpoint_definitions_sha256': {'storage': 'external_registry',
                        'sha256': model['endpoints_sha256'], 'source': 'SkinScout model registry',
                        'digest_format': 'canonical_json_utf8_v1'},
                },
                cache_manifest_sha256=digest(self.manifest),
                cache_reused=True,
                ranking_weight=0,
                applicability_domain=row["domain_status"],
                source_audit_sha256=self.manifest["source_audit_sha256"],
            )
            # Shared metadata appears once per response; every other row has its
            # exact content hashes, preserving units without megabytes of repeats.
            if row["model_id"] not in seen:
                row.update(
                    model_manifest=model["manifest"],
                    endpoint_definitions=model["endpoints"],
                )
                seen.add(row["model_id"])
        return rows


def build_cache(scout_root, output):
    import duckdb

    root, output = Path(scout_root), Path(output)
    output.mkdir(parents=True, exist_ok=False)
    audit_path = root / "reports/audit.json"
    audit = json.loads(audit_path.read_bytes())
    audit_sha = sha256_file(audit_path)
    if not (
        audit.get("audit_execution_passed")
        and audit.get("software_integrity_passed")
        and audit.get("auxiliary_admet_db_complete")
        and audit.get("hash_verification_performed")
    ):
        raise ValueError(
            "Complete published-source integrity audit required; scientific validation is independent"
        )
    if (
        any(v != 0 for v in audit.get("checks", {}).values())
        or audit.get("execution_errors")
        or audit.get("prediction_validation_errors")
    ):
        raise ValueError("Published ADMET integrity unresolved")
    code = root / "src/skinscout/audit.py"
    if sha256_file(code) != audit.get("audit_code_sha256"):
        raise ArtifactIntegrityError(
            "Donor audit implementation differs from its receipt"
        )
    started = time.monotonic()
    files = []
    models = {}
    counts = Counter()
    for index, path in enumerate(
        sorted((root / "data/published_predictions").glob("*.parquet"))
    ):
        if path.name == "_empty.parquet":
            continue
        model_id, filename = path.name.split("--", 1)
        source = root / "data/predictions" / model_id / filename
        receipt = json.loads(source.with_suffix(".json").read_bytes())
        sha = sha256_file(path)
        if (
            receipt.get("sha256") != sha
            or receipt.get("model_id") != model_id
            or source.stat().st_ino != path.stat().st_ino
        ):
            # Non-hardlinked deployments can still prove identical bytes.
            if (
                receipt.get("sha256") != sha
                or receipt.get("model_id") != model_id
                or sha256_file(source) != sha
            ):
                raise ArtifactIntegrityError(
                    "Published batch differs from committed prediction receipt"
                )
        files.append(
            {
                "file": str(path.relative_to(root)),
                "sha256": sha,
                "rows": receipt["rows"],
                "model_id": model_id,
                "input_hash": receipt["input_hash"],
                "receipt_sha256": sha256_file(source.with_suffix(".json")),
                "source_state": [
                    path.stat().st_dev,
                    path.stat().st_ino,
                    path.stat().st_size,
                    path.stat().st_mtime_ns,
                    path.stat().st_ctime_ns,
                ],
            }
        )
        if model_id not in models:
            registry = root / "models/registry" / model_id
            spec = json.loads((registry / "manifest.json").read_bytes())
            raw = (registry / "endpoints.csv").read_bytes()
            endpoints = list(csv.DictReader(io.StringIO(raw.decode())))
            if (
                spec.get("model_id") != model_id
                or hashlib.sha256(raw).hexdigest()
                != spec.get("endpoint_metadata_sha256")
                or spec.get("drugbank_enabled") is not False
            ):
                raise ArtifactIntegrityError(
                    "ADMET model/endpoint source identity mismatch"
                )
            if len(endpoints) != 52 or len({e["id"] for e in endpoints}) != 52:
                raise ValueError("ADMET-AI endpoint definitions incomplete")
            models[model_id] = {
                "manifest": spec,
                "manifest_sha256": hashlib.sha256(canonical_json(spec)).hexdigest(),
                "endpoints": endpoints,
                "endpoints_sha256": hashlib.sha256(
                    canonical_json(endpoints)
                ).hexdigest(),
                "endpoint_file_sha256": hashlib.sha256(raw).hexdigest(),
            }
        if index % 32 == 0:
            atomic_bytes(
                output / "progress.json",
                canonical_json(
                    {
                        "stage": "verify_published_bytes",
                        "batches": index + 1,
                        "elapsed_seconds": time.monotonic() - started,
                    }
                ),
            )
    if not files:
        raise ValueError("No committed published ADMET batches")
    database = output / "predictions.duckdb"
    with duckdb.connect(str(database)) as conn:
        conn.execute("SET threads=2; SET memory_limit='6GB'")
        atomic_bytes(
            output / "progress.json",
            canonical_json(
                {
                    "stage": "build_indexed_serving_copy",
                    "batches": len(files),
                    "elapsed_seconds": time.monotonic() - started,
                }
            ),
        )
        paths = [str(root / f["file"]) for f in files]
        conn.execute(
            "CREATE TABLE predictions AS SELECT * FROM read_parquet(?,union_by_name=true)",
            [paths],
        )
        count = conn.execute("SELECT count(*) FROM predictions").fetchone()[0]
        if count != sum(f["rows"] for f in files) or count != audit["compound_count"]:
            raise ValueError("ADMET serving population differs from audited source")
        conn.execute(
            "CREATE UNIQUE INDEX exact_model_identity ON predictions(compound_id,model_id)"
        )
        conn.execute("CREATE INDEX compound_lookup ON predictions(compound_id)")
        counts.update(
            dict(
                conn.execute(
                    "SELECT status,count(*) FROM predictions GROUP BY status"
                ).fetchall()
            )
        )
        if dict(counts) != audit["prediction_counts"]:
            raise ValueError("Prediction status counts differ from full audit")
        if not set(counts) <= TERMINAL:
            raise ValueError("Unrecognized prediction status")
        conn.execute("CHECKPOINT")
    if sha256_file(audit_path) != audit_sha:
        raise ArtifactIntegrityError("Source audit changed during serving copy")
    for f in files:
        path = root / f["file"]
        # The exact source bytes were hashed above. A changing source is held,
        # rather than declaring the last copied shard current.
        receipt = (
            root
            / "data/predictions"
            / f["model_id"]
            / Path(f["file"]).name.split("--", 1)[1]
        ).with_suffix(".json")
        if sha256_file(receipt) != f["receipt_sha256"]:
            raise ArtifactIntegrityError(
                "Source batch receipt changed during serving copy"
            )
        s = path.stat()
        if [s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns, s.st_ctime_ns] != f[
            "source_state"
        ]:
            raise ArtifactIntegrityError(
                "Source prediction bytes changed during serving copy"
            )
    manifest = {
        "schema": "scout-published-admet-cache-v1",
        "source_audit_sha256": audit_sha,
        "source_audit": audit,
        "database_sha256": sha256_file(database),
        "record_count": count,
        "status_counts": dict(counts),
        "models": models,
        "source_files": files,
        "elapsed_seconds": time.monotonic() - started,
        "scientific_validation": "NOT_ESTABLISHED",
        "use": "auxiliary_metrics",
        "automatic_safety_decision": False,
        "commercial_scope": "Code/weight terms recorded in original model manifest; dataset terms retain original qualification status",
    }
    atomic_bytes(output / "manifest.json", canonical_json(manifest))
    atomic_bytes(
        output / "progress.json",
        canonical_json(
            {
                "status": "complete",
                "rows": count,
                "elapsed_seconds": time.monotonic() - started,
            }
        ),
    )
    return {
        k: v
        for k, v in manifest.items()
        if k not in {"source_audit", "source_files", "models"}
    }


def add_point_lookup_index(cache_root, output):
    """New sealed generation; preserve the original published serving copy."""
    import duckdb

    source = PublishedPredictionCache(cache_root)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    database = output / "predictions.duckdb"
    started = time.monotonic()
    atomic_bytes(
        output / "progress.json",
        canonical_json({"stage": "copy_preserved_serving_data"}),
    )
    shutil.copyfile(source.database, database)
    if database.stat().st_ino == source.database.stat().st_ino:
        raise ArtifactIntegrityError("Serving upgrade must not mutate original inode")
    with duckdb.connect(str(database)) as conn:
        conn.execute("SET threads=2; SET memory_limit='4GB'")
        atomic_bytes(
            output / "progress.json",
            canonical_json({"stage": "build_single_column_point_lookup"}),
        )
        conn.execute(
            "CREATE INDEX IF NOT EXISTS compound_lookup ON predictions(compound_id)"
        )
        counts = dict(
            conn.execute(
                "SELECT status,count(*) FROM predictions GROUP BY status"
            ).fetchall()
        )
        if counts != source.manifest["status_counts"]:
            raise ArtifactIntegrityError("Serving upgrade population changed")
        conn.execute("CHECKPOINT")
    source.verify()
    manifest = {
        **source.manifest,
        "database_sha256": sha256_file(database),
        "previous_generation_sha256": digest(source.manifest),
        "lookup_contract": "single-column compound_id ART; bound IN list; exact compound/model uniqueness",
        "upgrade_elapsed_seconds": time.monotonic() - started,
    }
    atomic_bytes(output / "manifest.json", canonical_json(manifest))
    atomic_bytes(
        output / "progress.json",
        canonical_json(
            {
                "stage": "complete",
                "rows": manifest["record_count"],
                "elapsed_seconds": time.monotonic() - started,
            }
        ),
    )
    return {
        k: v
        for k, v in manifest.items()
        if k not in {"source_audit", "source_files", "models"}
    }
