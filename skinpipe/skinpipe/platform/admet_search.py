"""Reversible research search admission from a pinned auxiliary ADMET model.

High model outputs are research preferences, not calibrated toxicity estimates.
No donor rows are deleted, and unmeasured/unsupported predictions stay visible.
"""

from __future__ import annotations

from copy import deepcopy
import json
import math
from pathlib import Path
import re
import time

from .artifacts import ArtifactIntegrityError, atomic_bytes, canonical_json, sha256_file
from .identity import digest

RISK_ENDPOINTS = (
    "Skin_Reaction",
    "AMES",
    "Carcinogens_Lagunin",
    "ClinTox",
    "DILI",
    "hERG",
)
SEARCH_POLICY = {
    "version": "skin-toxicity-search-v2",
    "primary_endpoints": ["Skin_Reaction", "AMES"],
    "primary_threshold": 0.90,
    "combined_endpoints": ["Carcinogens_Lagunin", "ClinTox", "DILI", "hERG"],
    "combined_threshold": 0.90,
    "combined_minimum": 2,
    "missing": "retain",
    "interpretation": "uncalibrated model-output cutoff; not a safety determination",
    "excluded_scope": "search and new computation only; source records preserved",
    "endpoint_source": "https://tdcommons.ai/single_pred_tasks/tox/",
}


def search_policy():
    return {**deepcopy(SEARCH_POLICY), "sha256": digest(SEARCH_POLICY)}


def evaluate_risk(predictions, *, policy=None):
    """Keep model channels separate; never count missing values as zero/safe."""
    policy = deepcopy(policy) if policy is not None else search_policy()
    spec = {k: v for k, v in policy.items() if k != "sha256"}
    if (
        policy.get("sha256") != digest(spec)
        or policy.get("primary_endpoints") != SEARCH_POLICY["primary_endpoints"]
        or policy.get("combined_endpoints") != SEARCH_POLICY["combined_endpoints"]
        or policy.get("combined_minimum") != 2
        or any(
            type(policy.get(k)) not in (int, float) or not 0 < policy[k] <= 1
            for k in ("primary_threshold", "combined_threshold")
        )
    ):
        raise ValueError("Unverified frozen ADMET search policy")
    assessments = []
    for row in predictions or []:
        manifest = row.get("model_manifest") or {}
        mid = row.get("model_id") or manifest.get("model_id")
        # Other model endpoint names alone cannot establish equivalent semantics.
        if not isinstance(mid, str) or not mid.startswith("admet-ai-"):
            continue
        if row.get("status") not in {"success", "partial_nonfinite"}:
            continue
        values = row.get("values")
        if not isinstance(values, dict):
            continue
        known = {
            k: v
            for k in RISK_ENDPOINTS
            if type(v := values.get(k)) in (int, float)
            and 0 <= v <= 1
            and math.isfinite(v)
        }
        primary = [
            k
            for k in policy["primary_endpoints"]
            if known.get(k, -1) >= policy["primary_threshold"]
        ]
        combined = [
            k
            for k in policy["combined_endpoints"]
            if known.get(k, -1) >= policy["combined_threshold"]
        ]
        reasons = primary + (
            ["multiple_toxicity_predictions"]
            if len(combined) >= policy["combined_minimum"]
            else []
        )
        assessments.append(
            {
                "model_id": mid,
                "model_manifest_sha256": row.get("model_manifest_sha256"),
                "endpoint_definitions_sha256": row.get("endpoint_definitions_sha256"),
                "values": known,
                "missing_endpoints": [k for k in RISK_ENDPOINTS if k not in known],
                "strong_primary_endpoints": primary,
                "strong_combined_endpoints": combined,
                "excluded": bool(reasons),
                "reasons": reasons,
                "applicability_domain": row.get(
                    "applicability_domain", row.get("domain_status")
                ),
            }
        )
    return {
        "policy": policy,
        "assessments": assessments,
        "excluded": any(a["excluded"] for a in assessments),
        "status": "assessed" if assessments else "unassessed",
        "safety_status": "NOT_ESTABLISHED",
    }


def annotate_risk(row, *, policy=None):
    # Always derive from the original model rows rather than trusting a stale flag.
    row["admet_search_risk"] = evaluate_risk(row.get("admet"), policy=policy)
    return row


def _attach(connection, path, alias):
    # Paths are administrator-configured local files; no SQL is accepted by APIs.
    connection.execute(
        "ATTACH '"
        + str(Path(path)).replace("'", "''")
        + "' AS "
        + alias
        + " (READ_ONLY)"
    )


def build_search_gate(scout_root, cache_root, output):
    """Materialize cheap eligibility once, independently of target execution."""
    import duckdb
    from .scout import ScoutCatalog
    from .scout_predictions import PublishedPredictionCache

    source = ScoutCatalog(scout_root)
    cache = PublishedPredictionCache(cache_root)
    model_ids = list(cache.manifest["models"])
    if len(model_ids) != 1 or not model_ids[0].startswith("admet-ai-"):
        raise ValueError(
            "Exactly one pinned ADMET-AI model required for search admission"
        )
    model_id = model_ids[0]
    model = cache.manifest["models"][model_id]
    definitions = {e["id"]: e for e in model["endpoints"]}
    if any(
        k not in definitions or definitions[k].get("task_type") != "classification"
        for k in RISK_ENDPOINTS
    ):
        raise ValueError(
            "Classification endpoint definitions required before risk admission"
        )
    source_seal = source.structure_fingerprint()
    state = source.database.stat()
    source_state = (state.st_size, state.st_mtime_ns, state.st_ctime_ns)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    database = output / "catalog.duckdb"
    started = time.monotonic()
    # The read-only donor is attached; only the separate output is writable.
    with duckdb.connect(str(database)) as conn:
        conn.execute("SET threads=2; SET memory_limit='6GB'")
        _attach(conn, source.database, "donor")
        _attach(conn, cache.database, "admet_cache")
        paths = ["$." + k for k in RISK_ENDPOINTS]
        extraction = ",".join(
            "TRY_CAST(risk_json[" + str(i + 1) + "] AS DOUBLE) AS risk_" + str(i)
            for i in range(len(paths))
        )
        atomic_bytes(
            output / "progress.json",
            canonical_json({"stage": "extract_risk_and_purchase_population"}),
        )
        conn.execute(
            """CREATE TABLE catalog AS WITH prediction AS (
            SELECT compound_id,status,json_extract(values_json,?) AS risk_json
            FROM admet_cache.predictions WHERE model_id=?), purchase AS (
            SELECT DISTINCT compound_id FROM donor.offers
            WHERE availability='catalog_reported_in_stock'
            AND try_cast(evidence_date AS DATE) IS NOT NULL
            AND purchase_route IS NOT NULL AND purchase_route<>'')
            SELECT c.*, p.status AS admet_status, o.compound_id IS NOT NULL AS purchase_record,
            """
            + extraction
            + """ FROM donor.compounds c
            LEFT JOIN prediction p USING(compound_id) LEFT JOIN purchase o USING(compound_id)""",
            [paths, model_id],
        )
        count = conn.execute("SELECT count(*) FROM catalog").fetchone()[0]
        if (
            count != cache.manifest["record_count"]
            or conn.execute(
                "SELECT count(*) FROM catalog WHERE admet_status IS NULL"
            ).fetchone()[0]
        ):
            raise ArtifactIntegrityError("Prediction and source populations differ")
        risk_cols = ["risk_" + str(i) for i in range(6)]
        if conn.execute(
            "SELECT count(*) FROM catalog WHERE "
            + " OR ".join(
                f"({c} IS NOT NULL AND (NOT isfinite({c}) OR {c}<0 OR {c}>1))"
                for c in risk_cols
            )
        ).fetchone()[0]:
            raise ArtifactIntegrityError(
                "Risk probabilities outside documented endpoint range"
            )
        qualified = "admet_status IN ('success','partial_nonfinite')"
        primary = " OR ".join(f"coalesce(risk_{i}>=0.90,false)" for i in (0, 1))
        combined = "+".join(
            f"CAST(coalesce(risk_{i}>=0.90,false) AS INTEGER)" for i in range(2, 6)
        )
        conn.execute("ALTER TABLE catalog ADD COLUMN severe_admet BOOLEAN")
        conn.execute(
            "UPDATE catalog SET severe_admet=("
            + qualified
            + ") AND (("
            + primary
            + ") OR ("
            + combined
            + ")>=2)"
        )
        conn.execute("CREATE UNIQUE INDEX exact_identity ON catalog(compound_id)")
        metrics = {}
        for cutoff in (0.80, 0.90, 0.95, 0.99):
            skin = f"coalesce(risk_0>={cutoff},false)"
            ames = f"coalesce(risk_1>={cutoff},false)"
            many = "+".join(
                f"CAST(coalesce(risk_{i}>={cutoff},false) AS INTEGER)"
                for i in range(2, 6)
            )
            rule = f"({qualified}) AND ({skin} OR {ames} OR ({many})>=2)"
            names = (
                "population",
                "purchase_population",
                "success",
                "excluded",
                "purchase_excluded",
                "skin_flag",
                "ames_flag",
                "combined_flag",
            )
            values = conn.execute(
                "SELECT count(*),count(*) FILTER(WHERE purchase_record),"
                "count(*) FILTER(WHERE admet_status='success'),count(*) FILTER(WHERE "
                + rule
                + "),"
                "count(*) FILTER(WHERE purchase_record AND "
                + rule
                + "),count(*) FILTER(WHERE "
                + qualified
                + " AND "
                + skin
                + "),"
                "count(*) FILTER(WHERE "
                + qualified
                + " AND "
                + ames
                + "),count(*) FILTER(WHERE "
                + qualified
                + " AND ("
                + many
                + ")>=2) FROM catalog"
            ).fetchone()
            metrics[str(cutoff)] = dict(zip(names, values))
        conn.execute("CHECKPOINT")
    state = source.database.stat()
    if source_seal != source.structure_fingerprint() or source_state != (
        state.st_size,
        state.st_mtime_ns,
        state.st_ctime_ns,
    ):
        raise ArtifactIntegrityError("Source changed during search snapshot build")
    cache.verify()
    current = metrics["0.9"]
    manifest = {
        "schema": "scout-admet-search-gate-v1",
        "policy": search_policy(),
        "model_id": model_id,
        "model_manifest_sha256": model["manifest_sha256"],
        "endpoint_definitions_sha256": model["endpoints_sha256"],
        "source_structure_fingerprint": source_seal,
        "source_database_identity": source.database_identity(),
        "purchase_record_policy": "dated-supplier-stock-report-v2",
        "source_audit_sha256": cache.manifest["source_audit_sha256"],
        "prediction_cache_manifest_sha256": digest(cache.manifest),
        "database_sha256": sha256_file(database),
        "record_count": count,
        "retained_count": count - current["excluded"],
        "purchase_retained_count": current["purchase_population"]
        - current["purchase_excluded"],
        "threshold_sensitivity": metrics,
        "elapsed_seconds": time.monotonic() - started,
        "retained_interpretation": "eligible for research search; safety NOT_ESTABLISHED",
    }
    atomic_bytes(output / "manifest.json", canonical_json(manifest))
    atomic_bytes(
        output / "progress.json",
        canonical_json(
            {
                "stage": "complete",
                **current,
                "elapsed_seconds": time.monotonic() - started,
            }
        ),
    )
    return manifest


class ADMETSearchGate:
    def __init__(self, root):
        self.root = Path(root)
        self.database = self.root / "catalog.duckdb"
        self.manifest = json.loads((self.root / "manifest.json").read_bytes())
        self._state = None
        legacy = {
            **deepcopy(SEARCH_POLICY),
            "version": "skin-toxicity-search-v1",
            "primary_threshold": 0.95,
            "combined_threshold": 0.95,
        }
        supported = [search_policy(), {**legacy, "sha256": digest(legacy)}]
        if (
            self.manifest.get("schema") != "scout-admet-search-gate-v1"
            or self.manifest.get("policy") not in supported
        ):
            raise ValueError("Unsupported or changed search risk policy")
        self.verify()

    def verify(self):
        s = self.database.stat()
        state = (s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns, s.st_ctime_ns)
        if state != self._state:
            if sha256_file(self.database) != self.manifest["database_sha256"]:
                raise ArtifactIntegrityError("Search eligibility snapshot changed")
            s = self.database.stat()
            if state != (s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns, s.st_ctime_ns):
                raise ArtifactIntegrityError(
                    "Search snapshot changed during verification"
                )
            self._state = state
        return state

    def query(
        self,
        *,
        q="",
        limit=50,
        offset=0,
        purchasable_only=False,
        hide_severe_admet=True,
    ):
        import duckdb
        from .catalog_query_lifetime import interruptible_query, query_checkpoint

        query_checkpoint()
        if (
            not isinstance(q, str)
            or len(q) > 1024
            or type(limit) is not int
            or not 1 <= limit <= 1000
            or type(offset) is not int
            or not 0 <= offset <= 100000
        ):
            raise ValueError("Invalid bounded search request")
        if type(purchasable_only) is not bool or type(hide_severe_admet) is not bool:
            raise ValueError("Search switches must be boolean")
        clauses, values = [], []
        if q:
            if re.fullmatch(r"[0-9a-fA-F]{64}", q):
                clauses.append("compound_id=?")
                values.append(q.lower())
            elif re.fullmatch(r"[A-Z]{14}-[A-Z]{10}-[A-Z]", q.upper()):
                clauses.append("inchikey=?")
                values.append(q.upper())
            else:
                clauses.append(
                    "(compound_id=? OR inchikey=? OR contains(identity_smiles,?))"
                )
                values.extend([q, q, q])
        if purchasable_only:
            clauses.append("purchase_record")
        if hide_severe_admet:
            clauses.append("NOT severe_admet")
        before = self.verify()
        with duckdb.connect(str(self.database), read_only=True) as conn:
            conn.execute("SET threads=1; SET memory_limit='512MB'")
            with interruptible_query(conn):
                try:
                    c = conn.execute(
                        "SELECT * FROM catalog WHERE "
                        + (" AND ".join(clauses) or "TRUE")
                        + " ORDER BY compound_id LIMIT ? OFFSET ?",
                        [*values, limit, offset],
                    )
                    rows = [dict(zip([x[0] for x in c.description], r)) for r in c.fetchall()]
                except Exception:
                    query_checkpoint()
                    raise
            query_checkpoint()
        if self.verify() != before:
            raise ArtifactIntegrityError("Search snapshot changed during query")
        return rows


def update_search_policy(source_root, output):
    """New immutable generation from preserved endpoint outputs; no inference."""
    import duckdb

    original = ADMETSearchGate(source_root)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    database = output / "catalog.duckdb"
    with duckdb.connect(str(database)) as conn:
        conn.execute("SET threads=2; SET memory_limit='4GB'")
        _attach(conn, original.database, "previous")
        conn.execute("CREATE TABLE catalog AS SELECT * FROM previous.catalog")
        primary = " OR ".join(
            f"coalesce(risk_{i}>={SEARCH_POLICY['primary_threshold']},false)"
            for i in (0, 1)
        )
        combined = "+".join(
            f"CAST(coalesce(risk_{i}>={SEARCH_POLICY['combined_threshold']},false) AS INTEGER)"
            for i in range(2, 6)
        )
        conn.execute(
            "UPDATE catalog SET severe_admet=(admet_status IN ('success','partial_nonfinite')) AND (("
            + primary
            + ") OR ("
            + combined
            + ")>=2)"
        )
        conn.execute("CREATE UNIQUE INDEX exact_identity ON catalog(compound_id)")
        counts = conn.execute(
            "SELECT count(*),count(*) FILTER(WHERE severe_admet),count(*) FILTER(WHERE purchase_record AND severe_admet) FROM catalog"
        ).fetchone()
        conn.execute("CHECKPOINT")
    metrics = original.manifest["threshold_sensitivity"][
        str(SEARCH_POLICY["primary_threshold"])
    ]
    if counts != (
        metrics["population"],
        metrics["excluded"],
        metrics["purchase_excluded"],
    ):
        raise ArtifactIntegrityError(
            "Revised policy differs from measured sensitivity population"
        )
    original.verify()
    manifest = {
        **original.manifest,
        "policy": search_policy(),
        "database_sha256": sha256_file(database),
        "previous_generation_sha256": digest(original.manifest),
        "retained_count": counts[0] - counts[1],
        "purchase_retained_count": metrics["purchase_population"]
        - metrics["purchase_excluded"],
        "elapsed_seconds": time.monotonic() - started,
    }
    atomic_bytes(output / "manifest.json", canonical_json(manifest))
    atomic_bytes(
        output / "progress.json",
        canonical_json(
            {
                "stage": "complete",
                "retained_count": manifest["retained_count"],
                "purchase_retained_count": manifest["purchase_retained_count"],
                "elapsed_seconds": time.monotonic() - started,
            }
        ),
    )
    return manifest
