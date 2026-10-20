"""Immutable DuckDB/Parquet snapshots and bounded, whole-population SQL queries.

Publication streams rows in chunks. Query filtering and ranking happen in SQL
before the top-1,000 display cap; only one page reaches Python. The large Scout
catalog is queried in place by ``ScoutCatalog`` rather than copied into Git.
"""
from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import threading
import uuid

from .artifacts import ArtifactIntegrityError, atomic_bytes, canonical_json, fsync_directory, sha256_file

SCHEMA_VERSION = "skinpipe-catalog-v2"
_ID = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,127}\Z")
_COLUMNS = """candidate_id VARCHAR NOT NULL, target_id VARCHAR NOT NULL,
    name VARCHAR, cas VARCHAR, smiles VARCHAR, pchembl DOUBLE, topical_score DOUBLE,
    discovery_score DOUBLE, priority_score_v2 DOUBLE, mw DOUBLE, dock_affinity DOUBLE,
    purchasable BOOLEAN, critical_warning BOOLEAN, opposing_direction BOOLEAN,
    provenance_pass BOOLEAN, calculated BOOLEAN, computation_status VARCHAR,
    admet_status VARCHAR, measured BOOLEAN, payload VARCHAR NOT NULL, PRIMARY KEY(target_id,candidate_id)"""
SORTS = {key: key for key in ("candidate_id", "target_id", "name", "cas", "pchembl", "topical_score",
                            "discovery_score", "priority_score_v2", "mw", "dock_affinity", "purchasable")}
SORTS.update({"ligand_id": "candidate_id", "purchase_available": "purchasable", "docking_score": "dock_affinity"})
SORTS["primary_docking_score"] = "try_cast(json_extract_string(payload, '$.primary_docking.score') AS DOUBLE)"
FILTERS = {"hide_critical", "hide_opposing", "hide_severe_admet", "measured_only", "purchasable_only", "provenance_only",
           "include_missing_topical", "include_missing_pchembl", "min_topical", "min_pchembl",
           "calculated_only", "admet_status", "computation_status", "primary_docking_protocol"}


def _number(value):
    if isinstance(value, bool) or value is None:
        return None
    try:
        number = float(value)
        return number if math.isfinite(number) else None
    except (ValueError, TypeError, OverflowError):
        return None


def _dict(value):
    return value if isinstance(value, dict) else {}


def _row(candidate):
    row = dict(candidate)
    ligand = _dict(row.get("ligand"))
    candidate_id = row.get("candidate_id") or row.get("ligand_id") or row.get("compound_id") or ligand.get("ligand_id")
    target_id = row.get("target_id") or row.get("screen_id") or row.get("target")
    if not isinstance(candidate_id, str) or not candidate_id or not isinstance(target_id, str) or not target_id:
        raise ValueError("Each candidate needs a string candidate_id and target_id/screen_id")
    row.setdefault("candidate_id", candidate_id)
    row.setdefault("target_id", target_id)
    skin = _dict(row.get("skin"))
    pose = _dict(row.get("pose"))
    direction = row.get("direction_verdict") or _dict(row.get("direction")).get("verdict")
    hazard = row.get("admet_hazard") or _dict(row.get("admet")).get("hazard_severity")
    from .procurement import has_purchase_record
    purchase = True if has_purchase_record(row) else row.get("purchasable", row.get("purchase_available"))
    if not isinstance(purchase, bool):
        purchase = None
    affinity = _number(row.get("dock_affinity", pose.get("affinity")))
    execution = _dict(row.get("docking_execution"))
    status = row.get("computation_status") or execution.get("status") or ("completed" if affinity is not None else "not_calculated")
    calculated = affinity is not None and status not in {"failed", "invalid", "quarantined"}
    docking = row.get("docking") or []
    if isinstance(docking, dict):
        docking = [docking]
    if isinstance(docking, list) and any(isinstance(item, dict) and item.get("artifact_valid") is True
            and item.get("status") not in {"failed", "invalid", "unsupported", "quarantined"}
            and _number(item.get("score") if item.get("score") is not None else item.get("affinity")) is not None for item in docking):
        calculated = True
        if status in {"not_calculated", "missing", "pending"}:
            status = "completed"
    from ..measurements import activity_status
    observation = dict(row)
    measurements = row.get("measurements") or []
    if isinstance(measurements, dict):
        measurements = [measurements]
    if not isinstance(measurements, list):
        measurements = []
    measurements = [item for item in measurements if isinstance(item, dict)]
    observation["measurements"] = measurements
    if row.get("source") == "cosmetic_active" and not row.get("curated_claim_status") and not row.get("curated_claim_policy"):
        observation["curated_claim_status"] = "UNREVIEWED_TARGET_ENDPOINT"
    if row.get("activity_status") == "CURATED_CLAIM_ONLY":
        observation["curated_claim_status"] = "UNREVIEWED_TARGET_ENDPOINT"
    try:
        measured = activity_status(observation) in {"MEASURED", "MEASURED_QUALITATIVE", "ACTIVITY_VALUE_RECORDED"}
    except OverflowError:
        measured = False
    for item in measurements:
        if not isinstance(item, dict):
            continue
        item = item.get("normalized", item)
        if not isinstance(item, dict):
            continue
        unsuccessful = {"failed", "invalid", "unsupported", "not_evaluated", "unmeasured", "measurement_failed",
                        "measurement_unavailable", "unmeasurable", "not_measured"}
        if (item.get("status") not in unsuccessful and item.get("artifact_valid") is not False
                and (item.get("is_measured") is True or (item.get("status") == "measured" and item.get("kind") == "measured"))
                and _number(item.get("value")) is not None and item.get("endpoint") and item.get("unit")):
            measured = True
    return (candidate_id, target_id, row.get("name", ligand.get("name")), row.get("cas", ligand.get("cas")),
            row.get("smiles", ligand.get("smiles")), _number(row.get("pchembl", ligand.get("pchembl"))),
            _number(row.get("topical_score", skin.get("topical_score"))), _number(row.get("discovery_score")),
            _number(row.get("priority_score_v2")), _number(row.get("mw", ligand.get("mw"))), affinity,
            purchase, row.get("critical_warning") is True or hazard == "critical",
            row.get("opposing_direction") is True or direction == "opposing",
            (row.get("provenance_policy_eligible", row.get("provenance_pass"))
             if isinstance(row.get("provenance_policy_eligible", row.get("provenance_pass")), bool) else None), calculated, status,
            row.get("admet_status", "missing" if row.get("admet") is None else "recorded"), measured, canonical_json(row).decode())


class Catalog:
    def __init__(self, root, *, memory_limit="1GB", threads=2):
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)
        if not re.fullmatch(r"[1-9][0-9]*(MB|GB)", memory_limit):
            raise ValueError("Invalid catalog memory limit")
        self.memory_limit = memory_limit
        self.threads = max(1, min(int(threads), 8))
        self._verified = {}
        self._verification_lock = threading.Lock()

    def _connect(self, path, *, read_only=False):
        import duckdb
        # Apply limits before DuckDB constructs its native worker pool.
        return duckdb.connect(str(path), read_only=read_only,
                              config={"memory_limit": self.memory_limit, "threads": self.threads})

    def _directory(self, snapshot_id):
        if not isinstance(snapshot_id, str) or not _ID.fullmatch(snapshot_id):
            raise ValueError("Invalid snapshot ID")
        return self.root / snapshot_id

    def create_snapshot(self, rows, snapshot_id=None, provenance=None):
        # Generic snapshots may lack project context, but cannot be admitted
        # through the research API until an owner is explicitly bound.
        provenance = dict(provenance or {})
        if provenance.get("owner_id") is not None:
            from .queue import identifier
            identifier(provenance["owner_id"], "owner")
        snapshot_id = snapshot_id or uuid.uuid4().hex
        destination = self._directory(snapshot_id)
        if destination.exists():
            raise FileExistsError("Immutable snapshot already exists: " + snapshot_id)
        staging = self.root / (".building-" + uuid.uuid4().hex)
        staging.mkdir()
        try:
            with self._connect(staging / "catalog.duckdb") as connection:
                connection.execute("CREATE TABLE candidates (" + _COLUMNS + ")")
                batch = []
                for row in rows:
                    batch.append(_row(row))
                    if len(batch) == 2000:
                        self._insert(connection, batch)
                        batch.clear()
                if batch:
                    self._insert(connection, batch)
                connection.execute("CREATE INDEX candidate_name ON candidates(name)")
                connection.execute("CREATE INDEX candidate_target ON candidates(target_id)")
                counts = dict(connection.execute("SELECT target_id,count(*) FROM candidates GROUP BY target_id").fetchall())
                # SQL parameterization is supported for COPY destination too.
                connection.execute("COPY (SELECT * FROM candidates ORDER BY target_id,candidate_id) TO ? (FORMAT PARQUET, COMPRESSION ZSTD)",
                                   [str(staging / "candidates.parquet")])
                connection.execute("CHECKPOINT")
            files = {name: {"sha256": sha256_file(staging / name), "size_bytes": (staging / name).stat().st_size}
                     for name in ("catalog.duckdb", "candidates.parquet")}
            manifest = {"schema_version": SCHEMA_VERSION, "snapshot_id": snapshot_id,
                        "created_at": datetime.now(timezone.utc).isoformat(), "candidate_count": sum(counts.values()),
                        "target_counts": counts, "files": files, "provenance": provenance or {}}
            raw = canonical_json(manifest)
            atomic_bytes(staging / "manifest.json", raw)
            atomic_bytes(staging / "manifest.sha256", hashlib.sha256(raw).hexdigest().encode())
            staging.rename(destination)
            fsync_directory(destination.parent)
            if provenance.get("owner_id") is not None:
                self._record_owner(snapshot_id, provenance["owner_id"], hashlib.sha256(raw).hexdigest(),
                                   {"source": "server_snapshot_publication"})
            return manifest
        except BaseException:
            shutil.rmtree(staging, ignore_errors=True)
            raise

    @staticmethod
    def _insert(connection, rows):
        if not rows:
            return
        # One vectorized INSERT is one atomic statement/transaction. Avoid
        # executemany's per-row execution and possible durable commit cost.
        columns = [list(values) for values in zip(*rows)]
        connection.execute("INSERT INTO candidates SELECT " + ",".join("unnest(?)" for _ in range(20)), columns)

    @staticmethod
    def _fingerprint(path):
        stat = path.stat()
        return stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns

    def _verified_manifest_bytes(self, snapshot_id):
        """Verify the sealed snapshot; row reads do not need its large metadata.

        Cache only immutable bytes. Fingerprint changes still trigger the same
        schema/hash checks, including ctime and inode changes after tampering.
        """
        directory = self._directory(snapshot_id)
        try:
            paths = [directory / name for name in ("manifest.json", "manifest.sha256", "catalog.duckdb", "candidates.parquet")]
            state = tuple(self._fingerprint(path) for path in paths)
            with self._verification_lock:
                cached = self._verified.get(snapshot_id)
                if cached and cached[0] == state:
                    return cached[1]
                raw = paths[0].read_bytes()
                if hashlib.sha256(raw).hexdigest() != paths[1].read_text().strip():
                    raise ArtifactIntegrityError("Snapshot manifest checksum mismatch")
                manifest = json.loads(raw)
                if manifest.get("snapshot_id") != snapshot_id or manifest.get("schema_version") != SCHEMA_VERSION:
                    raise ArtifactIntegrityError("Snapshot identity/schema mismatch")
                for path in paths[2:]:
                    expected = manifest["files"][path.name]
                    if path.stat().st_size != expected["size_bytes"] or sha256_file(path) != expected["sha256"]:
                        raise ArtifactIntegrityError("Snapshot artifact checksum mismatch: " + path.name)
                if tuple(self._fingerprint(path) for path in paths) != state:
                    raise ArtifactIntegrityError("Snapshot changed during verification")
                self._verified[snapshot_id] = (state, raw)
                return raw
        except (OSError, KeyError, json.JSONDecodeError) as error:
            raise ArtifactIntegrityError("Snapshot incomplete or unreadable: " + snapshot_id) from error

    def manifest(self, snapshot_id):
        # Callers receive independent metadata, never a mutable cached object.
        return json.loads(self._verified_manifest_bytes(snapshot_id))

    def list_snapshots(self):
        result = []
        for directory in sorted(self.root.iterdir()):
            if directory.is_dir() and _ID.fullmatch(directory.name):
                try:
                    result.append({**self.manifest(directory.name), "status": "verified"})
                except ArtifactIntegrityError as error:
                    owner = self.snapshot_owner(directory.name)
                    result.append({"snapshot_id": directory.name, "status": "quarantined", "reason": str(error),
                                   "provenance": {"owner_id": owner},
                                   "ownership_status": "verified" if owner else "unresolved"})
        return result

    def snapshot_owner(self, snapshot_id):
        """Read sealed ownership even when payload bytes are quarantined.

        A separately admitted owner receipt can recover a pre-ownership snapshot
        without mutating its original manifest or relaxing artifact verification.
        """
        directory = self._directory(snapshot_id)
        path = self.root / ".owners" / (snapshot_id + ".json")
        admitted = None
        try:
            receipt_raw = path.read_bytes()
            if hashlib.sha256(receipt_raw).hexdigest() == path.with_suffix(".sha256").read_text().strip():
                admitted = json.loads(receipt_raw)
                if admitted.get("snapshot_id") != snapshot_id:
                    admitted = None
        except (OSError, ValueError, TypeError):
            pass
        try:
            raw = (directory / "manifest.json").read_bytes()
            sha = hashlib.sha256(raw).hexdigest()
            if sha != (directory / "manifest.sha256").read_text().strip():
                return admitted.get("owner_id") if admitted else None
            manifest = json.loads(raw)
            if manifest.get("snapshot_id") != snapshot_id or manifest.get("schema_version") != SCHEMA_VERSION:
                return None
            owner = manifest.get("provenance", {}).get("owner_id")
            if isinstance(owner, str) and owner:
                return owner
            if admitted and admitted.get("manifest_sha256") == sha:
                return admitted.get("owner_id")
        except (OSError, KeyError, ValueError, TypeError):
            return admitted.get("owner_id") if admitted else None
        return None

    def register_snapshot_owner(self, snapshot_id, owner, *, expected_manifest_sha256, evidence):
        """Administrator migration only: bind trusted historical ownership once."""
        from .queue import identifier
        identifier(owner, "owner")
        if not isinstance(evidence, dict) or not evidence:
            raise ValueError("Trusted ownership migration evidence is required")
        manifest = self.manifest(snapshot_id)
        if hashlib.sha256(canonical_json(manifest)).hexdigest() != expected_manifest_sha256:
            raise ArtifactIntegrityError("Ownership migration manifest differs from admitted source")
        existing = self.snapshot_owner(snapshot_id)
        if existing not in {None, owner}:
            raise PermissionError("Immutable snapshot is bound to another owner")
        self._record_owner(snapshot_id, owner, expected_manifest_sha256, evidence)
        return {"snapshot_id": snapshot_id, "owner_id": owner, "manifest_preserved": True}

    def _record_owner(self, snapshot_id, owner, manifest_sha256, evidence):
        receipt = canonical_json({"snapshot_id": snapshot_id, "owner_id": owner,
                                  "manifest_sha256": manifest_sha256, "evidence": evidence})
        directory = self.root / ".owners"
        directory.mkdir(exist_ok=True)
        import fcntl
        with (directory / ".lock").open("a+b") as lock:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
            path = directory / (snapshot_id + ".json")
            if path.exists() and path.read_bytes() != receipt:
                raise ValueError("Ownership receipt is already frozen")
            atomic_bytes(path, receipt)
            atomic_bytes(path.with_suffix(".sha256"), hashlib.sha256(receipt).hexdigest().encode())

    @staticmethod
    def _where(target_id, q, filters):
        unknown = set(filters) - FILTERS
        if unknown:
            raise ValueError("Unknown candidate filters: " + ", ".join(sorted(unknown)))
        clauses, values = [], []
        if target_id is not None:
            clauses.append("target_id=?")
            values.append(target_id)
        if q:
            if not isinstance(q, str) or len(q) > 1024:
                raise ValueError("Search must be a string of at most 1024 characters")
            clauses.append("(contains(lower(coalesce(name,'')),?) OR contains(lower(candidate_id),?) OR contains(lower(coalesce(cas,'')),?) OR contains(lower(coalesce(smiles,'')),?))")
            values.extend([q.lower()] * 4)
        for name, clause in (("hide_critical", "NOT critical_warning"), ("hide_opposing", "NOT opposing_direction"),
                             ('hide_severe_admet', "NOT coalesce(try_cast(json_extract_string(payload,'$.admet_search_risk.excluded') AS BOOLEAN),false)"),
                             ("measured_only", "measured"), ("purchasable_only", "purchasable IS TRUE"),
                             ("provenance_only", "provenance_pass IS TRUE"), ("calculated_only", "calculated")):
            if name in filters and not isinstance(filters[name], bool):
                raise ValueError(name + " must be boolean")
            if filters.get(name):
                clauses.append(clause)
        for name, column, include in (("min_topical", "topical_score", "include_missing_topical"),
                                      ("min_pchembl", "pchembl", "include_missing_pchembl")):
            if include in filters and not isinstance(filters[include], bool):
                raise ValueError(include + " must be boolean")
            if name in filters:
                value = _number(filters[name])
                if value is None:
                    raise ValueError(name + " must be finite")
                clauses.append("(" + column + ">=?" + (" OR " + column + " IS NULL" if filters.get(include, True) else "") + ")")
                values.append(value)
        for name in ("admet_status", "computation_status"):
            if name in filters:
                if not isinstance(filters[name], str):
                    raise ValueError(name + " must be a string")
                clauses.append(name + "=?")
                values.append(filters[name])
        if "primary_docking_protocol" in filters:
            protocol = filters["primary_docking_protocol"]
            if not isinstance(protocol, str) or not re.fullmatch(r"[a-f0-9]{64}", protocol):
                raise ValueError("primary_docking_protocol must be an exact protocol SHA256")
            clauses.append("json_extract_string(payload, '$.primary_docking.protocol_id')=?")
            values.append(protocol)
        return " AND ".join(clauses) or "TRUE", values

    def query(self, snapshot_id, target_id=None, q="", filters=None, sort="priority_score_v2", descending=True,
              limit=1000, page_size=50, offset=0, purchase_first=True):
        for name, value, low, high in (("limit", limit, 1, 1000), ("page_size", page_size, 1, 50), ("offset", offset, 0, 1000)):
            if isinstance(value, bool) or not isinstance(value, int) or not low <= value <= high:
                raise ValueError(f"{name} must be between {low} and {high}")
        if sort not in SORTS or not isinstance(descending, bool) or not isinstance(purchase_first, bool):
            raise ValueError("Unsupported sort or direction")
        if sort == "primary_docking_score" and not (filters or {}).get("primary_docking_protocol"):
            raise ValueError("Raw docking score sorting requires one exact primary_docking_protocol")
        self._verified_manifest_bytes(snapshot_id)
        where, params = self._where(target_id, q, filters or {})
        base, base_params = self._where(target_id, "", {})
        with self._connect(self._directory(snapshot_id) / "catalog.duckdb", read_only=True) as connection:
            population = connection.execute("""SELECT count(*),count(*) FILTER(WHERE calculated),
                count(*) FILTER(WHERE critical_warning),count(*) FILTER(WHERE opposing_direction),
                count(*) FILTER(WHERE pchembl IS NULL),count(*) FILTER(WHERE topical_score IS NULL),
                count(*) FILTER(WHERE purchasable IS NULL),count(*) FILTER(WHERE admet_status='missing'),
                count(*) FILTER(WHERE provenance_pass IS NOT TRUE),
                count(*) FILTER(WHERE provenance_pass IS FALSE),
                count(*) FILTER(WHERE provenance_pass IS NULL) FROM candidates WHERE """ + base, base_params).fetchone()
            filtered, calculated = connection.execute("SELECT count(*),count(*) FILTER(WHERE calculated) FROM candidates WHERE " + where, params).fetchone()
            order = ("purchasable DESC NULLS LAST, " if purchase_first and SORTS[sort] != "purchasable" else "")
            order += SORTS[sort] + (" DESC" if descending else " ASC") + " NULLS LAST, candidate_id ASC, target_id ASC"
            length = max(0, min(page_size, limit - offset))
            selected = connection.execute("SELECT payload FROM candidates WHERE " + where + " ORDER BY " + order + " LIMIT ? OFFSET ?",
                                          [*params, length, offset]).fetchall() if length else []
        missing = dict(zip(("pchembl", "topical_score", "purchase", "admet"), population[4:8]))
        reasons = {"critical": population[2], "opposing": population[3], "provenance_unqualified": population[8],
                   "provenance_failed": population[9], "provenance_unknown": population[10]}
        return {"snapshot_id": snapshot_id, "target_id": target_id, "rows": [json.loads(row[0]) for row in selected],
                "counts": {"population": population[0], "filtered": filtered, "displayed": min(filtered, limit),
                           "calculated": calculated, "population_calculated": population[1], "missing": missing,
                           "reasons": reasons, "excluded": population[0] - filtered},
                "limit": limit, "page_size": page_size, "offset": offset,
                "next_offset": offset + len(selected) if offset + len(selected) < min(filtered, limit) else None,
                "sort": sort, "descending": descending, "purchase_first": purchase_first}

    def selection_pool(self, snapshot_id, target_id=None, *, limit=10000):
        """Low-cost SQL pool includes controls and warnings; selection handles roles."""
        if isinstance(limit, bool) or not isinstance(limit, int) or not 1 <= limit <= 10000:
            raise ValueError("Selection pool limit must be1..10000")
        self._verified_manifest_bytes(snapshot_id)
        where, params = self._where(target_id, "", {})
        with self._connect(self._directory(snapshot_id) / "catalog.duckdb", read_only=True) as connection:
            rows = connection.execute("SELECT payload FROM candidates WHERE " + where +
                                      " ORDER BY purchasable DESC NULLS LAST,priority_score_v2 DESC NULLS LAST,candidate_id,target_id LIMIT ?",
                                      [*params, limit]).fetchall()
        return [json.loads(row[0]) for row in rows]

    def iter_candidates(self, snapshot_id, target_id=None, *, batch_size=2000):
        """Stream a complete verified snapshot for derived feedback snapshots."""
        if isinstance(batch_size, bool) or not isinstance(batch_size, int) or not 1 <= batch_size <= 10000:
            raise ValueError("Snapshot iterator batch size must be1..10000")
        self._verified_manifest_bytes(snapshot_id)
        where, params = self._where(target_id, "", {})
        with self._connect(self._directory(snapshot_id) / "catalog.duckdb", read_only=True) as connection:
            cursor = connection.execute("SELECT payload FROM candidates WHERE " + where + " ORDER BY target_id,candidate_id", params)
            while True:
                rows = cursor.fetchmany(batch_size)
                if not rows:
                    return
                for row in rows:
                    yield json.loads(row[0])

    def verify(self, snapshot_id):
        manifest = self.manifest(snapshot_id)
        directory = self._directory(snapshot_id)
        with self._connect(directory / "catalog.duckdb", read_only=True) as connection:
            actual = dict(connection.execute("SELECT target_id,count(*) FROM candidates GROUP BY target_id").fetchall())
            parquet = connection.execute("SELECT count(*) FROM read_parquet(?)", [str(directory / "candidates.parquet")]).fetchone()[0]
        if actual != manifest["target_counts"] or sum(actual.values()) != manifest["candidate_count"] or parquet != manifest["candidate_count"]:
            raise ArtifactIntegrityError("Snapshot catalog/Parquet/manifest count mismatch")
        return {"snapshot_id": snapshot_id, "status": "verified", "candidate_count": parquet, "target_counts": actual}

    def get_candidate(self, snapshot_id, candidate_id, target_id=None):
        self._verified_manifest_bytes(snapshot_id)
        where = "candidate_id=?" + (" AND target_id=?" if target_id is not None else "")
        values = [candidate_id] + ([target_id] if target_id is not None else [])
        with self._connect(self._directory(snapshot_id) / "catalog.duckdb", read_only=True) as connection:
            rows = connection.execute("SELECT payload FROM candidates WHERE " + where + " LIMIT 2", values).fetchall()
        if len(rows) > 1:
            raise ValueError("Candidate exists in multiple target contexts; supply target_id")
        return json.loads(rows[0][0]) if rows else None

    def get_candidates(self, snapshot_id, candidate_ids, target_id=None):
        """Read at most 2,000 frozen rows using one verified connection.

        Match point-read semantics: request order, duplicate IDs, missing rows
        and ambiguous target contexts survive. Return independent decoded rows;
        no molecular identity, selection or scientific fields are rewritten.
        """
        if (not isinstance(candidate_ids, list) or len(candidate_ids) > 2000
                or any(not isinstance(cid, str) or not cid for cid in candidate_ids)):
            raise ValueError("Candidate batch requires at most2000 nonempty string IDs")
        before = self._verified_manifest_bytes(snapshot_id)
        if not candidate_ids:
            return []
        where = "candidate_id IN (SELECT unnest(?))"
        values = [list(dict.fromkeys(candidate_ids))]
        if target_id is not None:
            where += " AND target_id=?"
            values.append(target_id)
        with self._connect(self._directory(snapshot_id) / "catalog.duckdb", read_only=True) as connection:
            rows = connection.execute(
                "SELECT candidate_id,payload FROM candidates WHERE " + where, values
            ).fetchall()
        by_id = {}
        for cid, payload in rows:
            if cid in by_id:
                raise ValueError("Candidate exists in multiple target contexts; supply target_id")
            by_id[cid] = payload
        if self._verified_manifest_bytes(snapshot_id) != before:
            raise ArtifactIntegrityError("Snapshot changed during candidate batch read")
        return [json.loads(by_id[cid]) if cid in by_id else None for cid in candidate_ids]
