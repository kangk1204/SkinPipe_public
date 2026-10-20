"""Read-only SkinScout adapter with independent source and prediction admission.

The donor database references local Parquet files: this adapter runs beside
those files on the registered worker. It does not copy the full catalog or
interpret a global auxiliary-prediction audit failure as a structure failure.
"""

from __future__ import annotations

import hashlib
import csv
import io
import json
import math
import re
from pathlib import Path
import threading

from .artifacts import ArtifactIntegrityError, sha256_file
from .identity import digest as content_digest
from .scout_predictions import normalize_prediction_values


class OffersSourceUnavailable(ValueError):
    """An incomplete supplier lookup cannot establish a negative purchase fact."""

    reason_code = "catalog_offers_source_unavailable"


def offer_evidence_status(detail):
    # Older gateways do not emit offers_status, but do retain component errors.
    # Missing supplier evidence is unknown, not a demonstrated negative.
    errors = detail.get("component_errors") or {}
    if "offers" in errors:
        return "source_unavailable"
    if "offer_exact_form" in errors:
        return "partial_acquisition"
    status = detail.get("offers_status", "complete")
    return status if status in {"complete", "source_unavailable", "partial_acquisition"} else "source_unavailable"


def require_complete_offers(detail):
    if offer_evidence_status(detail) != "complete":
        raise OffersSourceUnavailable("Catalog supplier evidence acquisition incomplete")


def audit_diagnosis(report):
    """Keep failed SQL execution distinct from a demonstrated row mismatch."""
    checks = report.get("checks", {})
    value = checks.get("sql_prediction_mismatch")
    endpoint_query_recorded = "published_missing_endpoint_definitions" in checks
    if value and not endpoint_query_recorded:
        state = "query_execution_error_recorded"
        reason = "Legacy audit increments sql_prediction_mismatch on an exception; the endpoint query produced no count"
    elif value:
        state = "mismatch_or_query_error_unresolved"
        reason = "Legacy audit combines row mismatches with swallowed SQL exceptions"
    elif value == 0 and endpoint_query_recorded:
        state = "sql_checks_passed"
        reason = "Reported SQL checks completed; partition hashes and model identity still require admission"
    else:
        state = "not_audited"
        reason = "No completed SQL prediction audit available"
    return {
        "status": state,
        "reported_value": value,
        "actual_mismatching_rows": None,
        "endpoint_query_count_recorded": endpoint_query_recorded,
        "reason": reason,
    }


class ScoutCatalog:
    def __init__(self, root, *, query_timeout=20):
        self.root = Path(root)
        self.database = self.root / "data/skinscout.duckdb"
        self.query_timeout = max(1, min(float(query_timeout), 120))
        self._partitions = {}
        self.published_cache = None
        self.search_gate = None
        self._database_identity = None
        self._structure_dependencies = None

    def database_identity(self):
        """Pin DB and external structure/offer bytes; cache only unchanged files.

        Prediction/assay views are independent and do not enter this seal.
        DuckDB's structure table is local, while supplier records may be views.
        """
        def state(path):
            s = path.stat()
            return (s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns, s.st_ctime_ns)

        db_state = state(self.database)
        if self._structure_dependencies is None or self._structure_dependencies[0] != db_state:
            with self._connect() as connection:
                views = connection.execute(
                    "SELECT sql FROM duckdb_views() WHERE view_name IN ('compounds','offers','source_records')"
                ).fetchall()
            paths = set()
            for (sql,) in views:
                for name in re.findall(r"'([^']+\.parquet)'", sql or ""):
                    path = Path(name).resolve()
                    if not path.is_relative_to(self.root.resolve() / "data"):
                        raise ArtifactIntegrityError("Structure view references external donor path")
                    paths.add(path)
            if state(self.database) != db_state:
                raise ArtifactIntegrityError("Donor changed while resolving structure views")
            self._structure_dependencies = (db_state, tuple(sorted(paths)))
        paths = (self.database, *self._structure_dependencies[1])
        before = tuple(state(path) for path in paths)
        if self._database_identity is None or self._database_identity[0] != before:
            identity = {"sha256": sha256_file(self.database), "size_bytes": before[0][2],
                        "structure_dependencies": [
                            {"path": str(path.relative_to(self.root.resolve())),
                             "sha256": sha256_file(path), "size_bytes": entry[2]}
                            for path, entry in zip(paths[1:], before[1:])]}
            if tuple(state(path) for path in paths) != before:
                raise ArtifactIntegrityError("Donor database changed during verification")
            self._database_identity = (before, identity)
        return self._database_identity[1]

    def _connect(self):
        import duckdb

        return duckdb.connect(str(self.database), read_only=True,
                              config={"threads": 2, "memory_limit": "1GB"})

    def _query(self, sql, values=()):
        from .catalog_query_lifetime import interruptible_query, query_checkpoint
        query_checkpoint()
        # Cancellation bounds an accidental expensive scan of a donor view.
        with self._connect() as connection:
            timer = threading.Timer(self.query_timeout, connection.interrupt)
            timer.start()
            try:
                with interruptible_query(connection):
                    try:
                        cursor = connection.execute(sql, list(values))
                        columns = [item[0] for item in cursor.description]
                        rows = [dict(zip(columns, row)) for row in cursor.fetchall()]
                        query_checkpoint()
                        return rows
                    except Exception:
                        query_checkpoint()
                        raise
            finally:
                timer.cancel()
                timer.join()

    def status(self):
        audit_path, build_path = (
            self.root / "reports/audit.json",
            self.root / "reports/db_build.json",
        )
        audit = json.loads(audit_path.read_text()) if audit_path.exists() else {}
        build = json.loads(build_path.read_text()) if build_path.exists() else {}
        checks = audit.get("checks", {})
        structure_checks = (
            "missing_or_stale_build_report",
            "incomplete_admitted_structure_sources",
            "duplicate_compound_ids",
            "orphan_structure_links",
            "duplicate_source_records",
        )
        structure_verified = bool(
            build.get("admitted_structure_sources_complete")
        ) and all(checks.get(key) == 0 for key in structure_checks)
        return {
            "database": str(self.database),
            "compound_count": audit.get("compound_count", build.get("compounds")),
            "audited_at": audit.get("audited_at"),
            "structure_admission": "reported_verified"
            if structure_verified
            else "unverified",
            "offers_and_assays": "independent_of_prediction_admission",
            "prediction_audit": audit_diagnosis(audit),
            "prediction_admitted_partitions": len(self._partitions),
            "admet_search_gate": self.search_gate.manifest
            if self.search_gate
            else None,
            "published_admet_cache": (
                {
                    "record_count": self.published_cache.manifest["record_count"],
                    "status_counts": self.published_cache.manifest["status_counts"],
                    "source_audit_sha256": self.published_cache.manifest[
                        "source_audit_sha256"
                    ],
                }
                if self.published_cache
                else None
            ),
            "report_sha256": sha256_file(audit_path) if audit_path.exists() else None,
            "scientific_validation": "NOT_ESTABLISHED",
        }

    def structure_fingerprint(self):
        """Structure-only seal; unrelated prediction audits cannot invalidate it."""
        from .identity import digest

        audit_path, build_path = (
            self.root / "reports/audit.json",
            self.root / "reports/db_build.json",
        )
        audit = json.loads(audit_path.read_bytes())
        build = json.loads(build_path.read_bytes())
        keys = (
            "missing_or_stale_build_report",
            "incomplete_admitted_structure_sources",
            "duplicate_compound_ids",
            "orphan_structure_links",
            "duplicate_source_records",
        )
        checks = {k: audit.get("checks", {}).get(k) for k in keys}
        if not build.get("admitted_structure_sources_complete") or any(
            v != 0 for v in checks.values()
        ):
            raise ValueError("Structure admission unresolved")
        return digest(
            {
                "build_sha256": sha256_file(build_path),
                "structure_checks": checks,
                "compound_count": audit.get("compound_count", build.get("compounds")),
                "database_identity": self.database_identity(),
                "seal_version": 2,
            }
        )

    def query(
        self,
        *,
        q="",
        limit=50,
        offset=0,
        purchasable_only=False,
        hide_severe_admet=None,
    ):
        """Bounded live structure lookup; opaque bulk catalog stays on its host.

        Empty search scans the ordered compound ID index. Search covers exact
        IDs, InChIKeys and literal structure substrings. Name/CAS resolution
        remains in source identifier tables, never an inferred identity join.
        """
        if hide_severe_admet is not None and type(hide_severe_admet) is not bool:
            raise ValueError("ADMET search switch must be boolean")
        if self.search_gate:
            if (
                self.search_gate.manifest["source_structure_fingerprint"]
                != self.structure_fingerprint()
            ):
                raise ValueError("Search eligibility source changed; rebuild required")
            return self.search_gate.query(
                q=q,
                limit=limit,
                offset=offset,
                purchasable_only=purchasable_only,
                hide_severe_admet=True
                if hide_severe_admet is None
                else hide_severe_admet,
            )
        if hide_severe_admet is True:
            raise ValueError("ADMET search admission not configured")
        if not isinstance(q, str) or len(q) > 1024:
            raise ValueError("Invalid Scout search")
        if (
            isinstance(limit, bool)
            or not isinstance(limit, int)
            or not 1 <= limit <= 1000
        ):
            raise ValueError("Scout limit must be 1..1000")
        if (
            isinstance(offset, bool)
            or not isinstance(offset, int)
            or not 0 <= offset <= 100000
        ):
            raise ValueError(
                "Scout offset must be 0..100000; use a materialized snapshot for deeper queries"
            )
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
            clauses.append(
                "compound_id IN (SELECT compound_id FROM offers WHERE "
                "availability='catalog_reported_in_stock' "
                "AND try_cast(evidence_date AS DATE) IS NOT NULL "
                "AND purchase_route IS NOT NULL AND purchase_route<>'')"
            )
        where = " AND ".join(clauses) if clauses else "TRUE"
        return self._query(
            "SELECT * FROM compounds WHERE "
            + where
            + " ORDER BY compound_id LIMIT ? OFFSET ?",
            [*values, limit, offset],
        )

    def get_compound(self, compound_id, *, evidence_limit=1000):
        row = self.materialize([compound_id], evidence_limit=evidence_limit)[0]
        return None if row["missing"] else row

    def materialize(self, compound_ids, *, evidence_limit=100):
        """Read up to 1,000 exact identities with bounded supplier queries.

        The requested order and unknown IDs are preserved. Independent evidence
        errors and per-compound truncation never silently delete a structure.
        A bounded evidence budget prevents a bulk request from expanding into
        millions of source records. All source files remain read-only.
        """
        if (
            not isinstance(compound_ids, list)
            or not 1 <= len(compound_ids) <= 1000
            or any(
                not isinstance(cid, str) or not cid or len(cid) > 256
                for cid in compound_ids
            )
            or len(set(compound_ids)) != len(compound_ids)
        ):
            raise ValueError("compound_ids must contain 1..1000 distinct bounded IDs")
        if (
            isinstance(evidence_limit, bool)
            or not isinstance(evidence_limit, int)
            or not 1 <= evidence_limit <= 1000
            or evidence_limit * len(compound_ids) > 100000
        ):
            raise ValueError(
                "evidence_limit must be 1..1000 within a 100000-record component budget"
            )
        placeholders = ",".join("?" for _ in compound_ids)
        compounds = self._query(
            "SELECT * FROM compounds WHERE compound_id IN (" + placeholders + ")",
            compound_ids,
        )
        by_id = {row["compound_id"]: row for row in compounds}
        if len(by_id) != len(compounds):
            raise ArtifactIntegrityError("Duplicate Scout compound identity")
        results = {
            cid: {
                "compound_id": cid,
                "compound": by_id.get(cid),
                "missing": cid not in by_id,
                "evidence_limit": evidence_limit,
                "admet": [],
                "component_errors": {},
            }
            for cid in compound_ids
        }
        for name, table in (
            ("offers", "offers"),
            ("assays", "skin_assays"),
            ("natural_products", "natural_products"),
        ):
            for result in results.values():
                result[name], result[name + "_truncated"] = [], False
            # Small independent read-only queries also bound index/window
            # allocation. A failed chunk keeps its error and does not mark
            # successfully read sibling chunks unavailable.
            chunk_size = 16 if name == "offers" else len(compound_ids)
            for start in range(0, len(compound_ids), chunk_size):
                ids = compound_ids[start : start + chunk_size]
                try:
                    rows = self._query(
                        f"""SELECT *, count(*) OVER (PARTITION BY compound_id) AS _skinpipe_source_count,
                        row_number() OVER (PARTITION BY compound_id ORDER BY to_json(e)) AS _skinpipe_row_index
                        FROM {table} e WHERE compound_id IN ({','.join('?' for _ in ids)})
                        QUALIFY _skinpipe_row_index <= ?""",
                        [*ids, evidence_limit],
                    )
                    for row in rows:
                        count = row.pop("_skinpipe_source_count")
                        row.pop("_skinpipe_row_index")
                        result = results[row["compound_id"]]
                        result[name].append(row)
                        result[name + "_truncated"] = count > evidence_limit
                except Exception as error:
                    for cid in ids:
                        results[cid]["component_errors"][name] = {
                            "type": type(error).__name__,
                            "message": str(error),
                        }
        for start in range(0, len(compound_ids), 16):
            ids = compound_ids[start : start + 16]
            offer_ids = sorted({offer["record_id"] for cid in ids
                                for offer in results[cid]["offers"] if offer.get("record_id")})
            if not offer_ids:
                continue
            try:
                source_rows = self._query(
                    "SELECT record_id,raw_smiles,source_version,rights_status FROM source_records WHERE compound_id IN ("
                    + ",".join("?" for _ in ids)
                    + ") AND record_id IN (SELECT unnest(?))",
                    [*ids, offer_ids],
                )
                source_by_id = {r["record_id"]: r for r in source_rows}
                for cid in ids:
                    result = results[cid]
                    for offer in result["offers"]:
                        raw = source_by_id.get(offer.get("record_id"))
                        if raw is None:
                            continue
                        offer.update(
                            exact_smiles=raw["raw_smiles"],
                            checked_at=str(offer.get("evidence_date") or ""),
                            product_url=offer.get("purchase_route"),
                            source=offer.get("source_id"),
                            available=offer.get("availability")
                            == "catalog_reported_in_stock",
                            exact_form_evidence={
                                "policy": "raw_supplier_record_link-v1",
                                **raw,
                            },
                            stock_status="unconfirmed",
                        )
            except Exception as error:
                for cid in ids:
                    results[cid]["component_errors"]["offer_exact_form"] = {
                        "type": type(error).__name__
                    }
        if self.published_cache is not None:
            try:
                for prediction in self.published_cache.query(compound_ids):
                    result = results[prediction["compound_id"]]
                    compound = result["compound"] or {}
                    if compound.get("canonical_smiles") != prediction["input_smiles"]:
                        result["component_errors"]["published_admet_identity"] = {
                            "type": "ArtifactIntegrityError"
                        }
                        continue
                    result["admet"].append(prediction)
            except Exception as error:
                for result in results.values():
                    result["component_errors"]["published_admet_cache"] = {
                        "type": type(error).__name__
                    }
        for digest, partition in list(self._partitions.items()):
            try:
                if sha256_file(partition["file"]) != digest:
                    raise ArtifactIntegrityError(
                        "Admitted prediction partition changed"
                    )
                values = self._query(
                    "SELECT * FROM read_parquet(?) WHERE compound_id IN ("
                    + placeholders
                    + ")",
                    [str(partition["file"]), *compound_ids],
                )
                for value in values:
                    result = results[value["compound_id"]]
                    if any(p["model_id"] == value["model_id"] for p in result["admet"]):
                        continue
                    if not result["compound"] or result["compound"].get(
                        "canonical_smiles"
                    ) != value.get("input_smiles"):
                        result["component_errors"]["admet:" + digest] = {
                            "type": "ArtifactIntegrityError",
                            "message": "Prediction exact structure no longer matches source",
                        }
                        continue
                    result["admet"].append(
                        {
                            **value,
                            "values": normalize_prediction_values(
                                value, partition["endpoint_definitions"]
                            ),
                            "model_manifest": partition["model_manifest"],
                            "model_manifest_sha256": partition["model_manifest_sha256"],
                            "endpoint_definitions": partition["endpoint_definitions"],
                            "endpoint_definitions_sha256": content_digest(partition["endpoint_definitions"]),
                            'metadata_refs': {
                                'model_manifest_sha256': {'storage': 'external_registry',
                                    'sha256': partition['model_manifest_sha256'],
                                    'source': 'SkinScout model registry', 'digest_format': 'canonical_json_ascii_v1'},
                                'endpoint_definitions_sha256': {'storage': 'external_registry',
                                    'sha256': content_digest(partition['endpoint_definitions']),
                                    'source': 'SkinScout model registry', 'digest_format': 'canonical_json_ascii_v1'},
                            },
                            "ranking_weight": 0,
                            "applicability_domain": value.get("domain_status"),
                            "partition_sha256": digest,
                        }
                    )
            except Exception as error:
                for result in results.values():
                    result["component_errors"]["admet:" + digest] = {
                        "type": type(error).__name__,
                        "message": str(error),
                    }
        for result in results.values():
            offer_error = any(key in result["component_errors"]
                              for key in ("offers", "offer_exact_form"))
            result["offers_status"] = (
                "partial_acquisition" if offer_error and result["offers"]
                else "source_unavailable" if offer_error else "complete"
            )
            result["purchasable"] = None if offer_error else any(
                offer.get("available") is True
                or offer.get("availability") in {"catalog_reported_in_stock", "catalog_available"}
                for offer in result["offers"]
            )
            admet_error = any(
                key.startswith(("published_admet", "admet:"))
                for key in result["component_errors"]
            )
            result["admet_status"] = (
                "partial_acquisition"
                if admet_error and result["admet"]
                else "source_unavailable"
                if admet_error
                else
                "cache_reused"
                if any(p.get("cache_reused") for p in result["admet"])
                else "admitted_partition"
                if result["admet"]
                else "not_admitted"
            )
        return [results[cid] for cid in compound_ids]

    def admit_prediction_partition(
        self, file, *, model_manifest, endpoint_definitions,
        model_manifest_sha256=None,
    ):
        """Validate a bounded donor batch before exposing its numeric values."""
        file = Path(file).resolve()
        if not file.is_relative_to(self.root.resolve() / "data/predictions"):
            raise ValueError("Prediction partition must belong to the configured donor")
        receipt_path = file.with_suffix(".json")
        receipt = json.loads(receipt_path.read_text())
        digest = sha256_file(file)
        if receipt.get("sha256") != digest:
            raise ArtifactIntegrityError("Prediction partition checksum mismatch")
        count = receipt.get("rows")
        if not isinstance(count, int) or not 1 <= count <= 100000:
            raise ValueError("Prediction admission requires a bounded batch")
        model = (
            json.loads(Path(model_manifest).read_text())
            if not isinstance(model_manifest, dict)
            else model_manifest
        )
        if (
            not isinstance(model_manifest_sha256, str)
            or not re.fullmatch(r"[0-9a-f]{64}", model_manifest_sha256)
            or content_digest(model) != model_manifest_sha256
        ):
            raise ArtifactIntegrityError("Operator-pinned model manifest checksum required")
        if not re.fullmatch(r"[0-9a-f]{64}", str(model.get("endpoint_metadata_sha256", ""))):
            raise ArtifactIntegrityError("Pinned endpoint metadata checksum required")
        if isinstance(endpoint_definitions, list):
            endpoints = endpoint_definitions
            if model.get("endpoint_metadata_sha256"):
                raise ValueError(
                    "This model requires the exact endpoint metadata file for checksum validation"
                )
        else:
            endpoint_file = Path(endpoint_definitions)
            raw = endpoint_file.read_bytes()
            if (
                model.get("endpoint_metadata_sha256")
                and hashlib.sha256(raw).hexdigest() != model["endpoint_metadata_sha256"]
            ):
                raise ArtifactIntegrityError(
                    "Prediction endpoint metadata checksum mismatch"
                )
            endpoints = (
                list(csv.DictReader(io.StringIO(raw.decode())))
                if endpoint_file.suffix == ".csv"
                else json.loads(raw)
            )
        expected = {entry["id"] for entry in endpoints}
        mid = receipt.get("model_id")
        if (
            not expected
            or len(expected) != len(endpoints)
            or not mid
            or model.get("model_id") != mid
            or file.parent.name != mid
        ):
            raise ValueError("Prediction model/endpoint identity mismatch")
        if len(expected) != 52:
            raise ValueError("ADMET-AI requires its 52 endpoint definitions")
        if model.get("drugbank_enabled") is not False:
            raise ArtifactIntegrityError("ADMET source must explicitly disable DrugBank enrichment")
        rows = self._query(
            "SELECT * FROM read_parquet(?) LIMIT ?", [str(file), count + 1]
        )
        if len(rows) != count:
            raise ValueError("Prediction receipt row count mismatch")
        ids, inputs = set(), []
        terminal = {
            "success",
            "failed",
            "partial_nonfinite",
            "unsupported_query_or_polymer",
            "unsupported_multi_component",
            "unsupported_element",
            "outside_operational_size_limit",
            "unsupported_enhanced_stereo",
        }
        for row in rows:
            identity, status = row["compound_id"], row["status"]
            if (
                not identity
                or identity in ids
                or row["model_id"] != mid
                or status not in terminal
            ):
                raise ValueError("Prediction row identity/status mismatch")
            ids.add(identity)
            if not row.get("computed_at") or not row.get("domain_status"):
                raise ValueError("Prediction provenance incomplete")
            json.loads(row["quality_flags_json"])
            values = json.loads(row["values_json"])
            if status in {"success", "partial_nonfinite"}:
                if set(values) != expected:
                    raise ValueError("Prediction endpoints mismatch")
                finite = [
                    isinstance(value, (int, float))
                    and not isinstance(value, bool)
                    and math.isfinite(value)
                    for value in values.values()
                ]
                if any(
                    value is not None and not valid
                    for value, valid in zip(values.values(), finite)
                ):
                    raise ValueError("Invalid prediction values")
                if (status == "success") != all(finite):
                    raise ValueError("Prediction finite status mismatch")
            elif values:
                raise ValueError(
                    "Unsupported/failed prediction contains numerical values"
                )
            if status == "failed" and not row.get("error"):
                raise ValueError("Failed prediction lacks error")
            inputs.append(
                {
                    "compound_id": identity,
                    "canonical_smiles": row["input_smiles"],
                    "inference_status": "eligible"
                    if status in {"success", "failed", "partial_nonfinite"}
                    else status,
                }
            )
        input_hash = hashlib.sha256(
            json.dumps(inputs, sort_keys=True).encode()
        ).hexdigest()
        if (
            receipt.get("input_hash") != input_hash
            or file.stem != input_hash
            or any(row["input_hash"] != input_hash for row in rows)
        ):
            raise ValueError("Prediction input hash mismatch")
        # Validate against exact structure rows in bounded chunks; never accept
        # a parent/molecule-name match for salt/stereochemistry identity.
        for start in range(0, len(inputs), 500):
            chunk = inputs[start : start + 500]
            actual = self._query(
                "SELECT compound_id,canonical_smiles,inference_status FROM compounds WHERE compound_id IN ("
                + ",".join("?" for _ in chunk)
                + ")",
                [row["compound_id"] for row in chunk],
            )
            by_id = {row["compound_id"]: row for row in actual}
            if any(by_id.get(row["compound_id"]) != row for row in chunk):
                raise ValueError("Prediction exact input identity/eligibility mismatch")
        result = {
            "file": str(file),
            "sha256": digest,
            "rows": count,
            "model_id": mid,
            "receipt_sha256": sha256_file(receipt_path),
            "model_manifest": model,
            "model_manifest_sha256": model_manifest_sha256,
            "endpoint_definitions": endpoints,
            "scientific_validation": "NOT_ESTABLISHED",
        }
        self._partitions[digest] = result
        return result

    def diagnose_endpoint_query(self, partition):
        """Reproduce the formerly swallowed SQL failure on one bounded shard."""
        try:
            partition = Path(partition).resolve()
            if not partition.is_relative_to((self.root / "data/predictions").resolve()):
                raise ArtifactIntegrityError("Diagnostic partition must belong to donor")
            receipt = json.loads(partition.with_suffix(".json").read_bytes())
            if receipt.get("sha256") != sha256_file(partition):
                raise ArtifactIntegrityError("Diagnostic partition checksum mismatch")
            rows = self._query(
                """SELECT count(*) AS mismatching_endpoint_values FROM read_parquet(?) p,
                LATERAL json_each(p.values_json) j WHERE NOT EXISTS
                (SELECT 1 FROM endpoints e WHERE e.model_id=p.model_id AND e.endpoint=j.key)""",
                [str(partition)],
            )
            return {"status": "completed", "scope": "one_partition", **rows[0]}
        except Exception as error:
            return {
                "status": "query_error",
                "scope": "one_partition",
                "error_type": type(error).__name__,
                "error": str(error),
            }
