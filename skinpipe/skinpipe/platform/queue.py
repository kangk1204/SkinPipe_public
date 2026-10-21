"""PostgreSQL leased work with fenced publication and immutable attempt history.

Adapted mechanisms: SkinKairos queue at 67dff3e. Every state transition and its
event commit together. An expired owner can neither heartbeat nor publish.
Clients may retry submissions/receipts after ambiguous HTTP responses.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
import uuid
from contextlib import contextmanager

import psycopg
from psycopg.rows import dict_row
from psycopg.types.json import Jsonb


class LeaseLost(RuntimeError):
    pass


class Cancelled(RuntimeError):
    pass


class BudgetExhausted(ValueError):
    """Stopped work cannot resume within its unchanged frozen allowance."""


def canonical(value):
    return json.dumps(
        value,
        sort_keys=True,
        ensure_ascii=False,
        separators=(",", ":"),
        allow_nan=False,
    )


def digest(value):
    return hashlib.sha256(canonical(value).encode()).hexdigest()


def identifier(value, name="identifier"):
    if not isinstance(value, str) or not re.fullmatch(r"[A-Za-z0-9_.:-]{1,160}", value):
        raise ValueError(f"Invalid {name}")
    return value


def _cpu_reservation_limits(telemetry):
    """Separate one worker's affinity/limit from its physical host's capacity.

    Older nodes may omit cpu_available, cpu_count, or cpu_limit. A missing or
    invalid host count must never grant more than the old affinity-based cap;
    an explicitly invalid worker availability/limit cannot authorize a job.
    """
    from .resources import cpu_limit

    available = telemetry.get("cpu_available", telemetry.get("cpu_count"))
    if type(available) is not int or available < 1:
        return None
    affinity_limit = cpu_limit(available)
    worker_limit = affinity_limit
    if "cpu_limit" in telemetry:
        explicit = telemetry["cpu_limit"]
        if type(explicit) is not int or not 1 <= explicit <= 8:
            return None
        worker_limit = min(worker_limit, explicit)
    count = telemetry.get("cpu_count")
    host_limit = cpu_limit(count) if type(count) is int and count > 0 else affinity_limit
    return worker_limit, host_limit


SCHEMA = """
CREATE TABLE IF NOT EXISTS platform_users (
 id TEXT PRIMARY KEY, created_at TIMESTAMPTZ NOT NULL DEFAULT clock_timestamp());
CREATE TABLE IF NOT EXISTS platform_projects (
 id TEXT PRIMARY KEY, owner_id TEXT NOT NULL REFERENCES platform_users(id));
CREATE TABLE IF NOT EXISTS platform_compute_budgets (
 id TEXT PRIMARY KEY, owner_id TEXT NOT NULL, limit_seconds DOUBLE PRECISION NOT NULL,
 used_seconds DOUBLE PRECISION NOT NULL DEFAULT 0, active_job TEXT);
CREATE TABLE IF NOT EXISTS platform_jobs (
 id TEXT PRIMARY KEY, owner_id TEXT NOT NULL REFERENCES platform_users(id),
 project_id TEXT NOT NULL REFERENCES platform_projects(id),
 kind TEXT NOT NULL, request JSONB NOT NULL, request_sha256 TEXT NOT NULL,
 idempotency_key TEXT NOT NULL, status TEXT NOT NULL DEFAULT 'queued',
 created_at TIMESTAMPTZ NOT NULL DEFAULT clock_timestamp(),
 updated_at TIMESTAMPTZ NOT NULL DEFAULT clock_timestamp(),
 not_before TIMESTAMPTZ NOT NULL DEFAULT clock_timestamp(),
 worker TEXT, lease_token TEXT, lease_until TIMESTAMPTZ,
 attempts INTEGER NOT NULL DEFAULT 0, max_attempts INTEGER NOT NULL DEFAULT 3,
 cancel_requested BOOLEAN NOT NULL DEFAULT false,
 compute_seconds DOUBLE PRECISION NOT NULL DEFAULT 0,
 budget_seconds DOUBLE PRECISION NOT NULL DEFAULT 28800,
 progress JSONB NOT NULL DEFAULT '{}', result JSONB, error TEXT,
 UNIQUE(owner_id, idempotency_key));
ALTER TABLE platform_jobs ADD COLUMN IF NOT EXISTS budget_scope TEXT;
CREATE INDEX IF NOT EXISTS platform_jobs_claim ON platform_jobs(status,not_before,created_at);
CREATE TABLE IF NOT EXISTS platform_attempts (
 job_id TEXT NOT NULL REFERENCES platform_jobs(id), attempt INTEGER NOT NULL,
 token_sha256 TEXT NOT NULL, worker TEXT NOT NULL,
 started_at TIMESTAMPTZ NOT NULL DEFAULT clock_timestamp(),
 ended_at TIMESTAMPTZ, status TEXT NOT NULL DEFAULT 'running',
 PRIMARY KEY(job_id,attempt));
CREATE TABLE IF NOT EXISTS platform_job_events (
 id BIGSERIAL PRIMARY KEY, job_id TEXT NOT NULL REFERENCES platform_jobs(id),
 created_at TIMESTAMPTZ NOT NULL DEFAULT clock_timestamp(), type TEXT NOT NULL,
 payload JSONB NOT NULL);
CREATE TABLE IF NOT EXISTS platform_publications (
 job_id TEXT PRIMARY KEY REFERENCES platform_jobs(id),
 result_sha256 TEXT NOT NULL, result JSONB NOT NULL,
 created_at TIMESTAMPTZ NOT NULL DEFAULT clock_timestamp());
CREATE TABLE IF NOT EXISTS platform_nodes (
 id TEXT PRIMARY KEY, telemetry JSONB NOT NULL,
 updated_at TIMESTAMPTZ NOT NULL DEFAULT clock_timestamp());
CREATE TABLE IF NOT EXISTS platform_gpu_reservations (
 gpu_uuid TEXT PRIMARY KEY, job_id TEXT NOT NULL REFERENCES platform_jobs(id),
 lease_token TEXT NOT NULL, expires_at TIMESTAMPTZ NOT NULL);
CREATE TABLE IF NOT EXISTS platform_cpu_reservations (
 job_id TEXT PRIMARY KEY REFERENCES platform_jobs(id), node_id TEXT NOT NULL,
 lease_token TEXT NOT NULL, cores INTEGER NOT NULL CHECK(cores>0));
ALTER TABLE platform_cpu_reservations ADD COLUMN IF NOT EXISTS physical_host TEXT;
CREATE TABLE IF NOT EXISTS platform_launch_budgets (
 scope TEXT PRIMARY KEY, launch_limit INTEGER NOT NULL CHECK(launch_limit>0),
 used INTEGER NOT NULL DEFAULT 0);
CREATE TABLE IF NOT EXISTS platform_model_launches (
 scope TEXT NOT NULL REFERENCES platform_launch_budgets(scope),
 job_id TEXT NOT NULL REFERENCES platform_jobs(id), attempt INTEGER NOT NULL,
 work_key TEXT NOT NULL, created_at TIMESTAMPTZ NOT NULL DEFAULT clock_timestamp(),
 PRIMARY KEY(scope,job_id,attempt,work_key));
ALTER TABLE platform_model_launches ADD COLUMN IF NOT EXISTS spawn_authorized BOOLEAN NOT NULL DEFAULT false;
"""


class Queue:
    def __init__(self, dsn, *, lease_seconds=90, initialize=True):
        self.dsn = dsn
        if not 2 <= lease_seconds <= 600:
            raise ValueError("Lease must be between 2 and 600 seconds")
        self.lease_seconds = lease_seconds
        if initialize:
            with self.connection() as conn:
                conn.execute(SCHEMA)

    @contextmanager
    def connection(self):
        with psycopg.connect(
            self.dsn,
            connect_timeout=5,
            row_factory=dict_row,
            options="-c statement_timeout=15000 -c lock_timeout=5000",
        ) as conn:
            yield conn

    @staticmethod
    def _event(conn, job_id, event_type, payload=None):
        conn.execute(
            "INSERT INTO platform_job_events(job_id,type,payload) VALUES(%s,%s,%s)",
            (job_id, event_type, Jsonb(payload or {})),
        )

    def enqueue(
        self,
        kind,
        request,
        *,
        owner="researcher",
        project="default",
        idempotency_key,
        budget_seconds=28800,
        budget_scope=None,
        max_attempts=3,
    ):
        for value in (kind, owner, project, idempotency_key):
            identifier(value)
        if not isinstance(request, dict):
            raise ValueError("Job request must be an object")
        if isinstance(max_attempts, bool) or not isinstance(max_attempts, int) or not 1 <= max_attempts <= 3:
            raise ValueError("Job max_attempts must be within 1..3")
        if (
            isinstance(budget_seconds, bool)
            or not isinstance(budget_seconds, (int, float))
            or not math.isfinite(budget_seconds)
            or not 0 < budget_seconds <= 604800
        ):
            raise ValueError("Invalid computation budget")
        effective_scope = identifier(
            budget_scope if budget_scope is not None else request.get("budget_scope", idempotency_key),
            "budget scope",
        )
        signature = digest(
            {
                "kind": kind,
                "request": request,
                "project": project,
                "budget_seconds": budget_seconds,
                **({"max_attempts": max_attempts} if max_attempts != 3 else {}),
                **({"budget_scope": effective_scope} if budget_scope is not None
                   and effective_scope != request.get("budget_scope", idempotency_key) else {}),
            }
        )
        scope = owner + ":" + effective_scope
        with self.connection() as conn:
            conn.execute(
                "SELECT pg_advisory_xact_lock(hashtextextended(%s,0))",
                (owner + ":" + idempotency_key,),
            )
            conn.execute(
                "INSERT INTO platform_users(id) VALUES(%s) ON CONFLICT DO NOTHING",
                (owner,),
            )
            conn.execute(
                "INSERT INTO platform_projects(id,owner_id) VALUES(%s,%s) ON CONFLICT DO NOTHING",
                (project, owner),
            )
            if (
                conn.execute(
                    "SELECT owner_id FROM platform_projects WHERE id=%s", (project,)
                ).fetchone()["owner_id"]
                != owner
            ):
                raise ValueError("Project belongs to another owner")
            row = conn.execute(
                "SELECT * FROM platform_jobs WHERE owner_id=%s AND idempotency_key=%s",
                (owner, idempotency_key),
            ).fetchone()
            if row:
                if row["request_sha256"] != signature:
                    raise ValueError(
                        "Idempotency key already bound to different request"
                    )
                bound = conn.execute(
                    "SELECT owner_id FROM platform_compute_budgets WHERE id=%s",
                    (row["budget_scope"],),
                ).fetchone()
                if bound is not None and bound["owner_id"] != owner:
                    raise ValueError("Compute budget belongs to another owner")
                return self._public(row)
            job_id = uuid.uuid4().hex
            conn.execute(
                "INSERT INTO platform_compute_budgets(id,owner_id,limit_seconds) VALUES(%s,%s,%s) ON CONFLICT DO NOTHING",
                (scope, owner, budget_seconds),
            )
            # Colons are valid in both identifiers, so concatenation alone
            # does not prove ownership (a:b + c and a + b:c collide). Retain
            # historical scope IDs but never share a ledger across owners.
            bound = conn.execute(
                "SELECT owner_id FROM platform_compute_budgets WHERE id=%s",
                (scope,),
            ).fetchone()
            if bound["owner_id"] != owner:
                raise ValueError("Compute budget belongs to another owner")
            conn.execute(
                """INSERT INTO platform_jobs
                (id,owner_id,project_id,kind,request,request_sha256,idempotency_key,budget_seconds,max_attempts)
                VALUES(%s,%s,%s,%s,%s,%s,%s,%s,%s)""",
                (
                    job_id,
                    owner,
                    project,
                    kind,
                    Jsonb(request),
                    signature,
                    idempotency_key,
                    budget_seconds,
                    max_attempts,
                ),
            )
            self._event(conn, job_id, "submitted", {"request_sha256": signature})
            conn.execute(
                "UPDATE platform_jobs SET budget_scope=%s WHERE id=%s", (scope, job_id)
            )
        return self.get(job_id, owner=owner)

    @staticmethod
    def _public(row):
        result = dict(row)
        result.pop("lease_token", None)
        for key, value in result.items():
            if hasattr(value, "isoformat"):
                result[key] = value.isoformat()
        return result

    def get(self, job_id, *, owner=None):
        with self.connection() as conn:
            row = conn.execute(
                "SELECT * FROM platform_jobs WHERE id=%s", (job_id,)
            ).fetchone()
            if not row or (owner is not None and row["owner_id"] != owner):
                raise KeyError("Job not found")
            return self._public(row)

    def get_statuses(self, job_ids, *, owner):
        """Read bounded, owner-scoped job states in one database round trip."""
        identifier(owner, "owner")
        if not isinstance(job_ids, (list, tuple)) or len(job_ids) > 10000:
            raise ValueError("At most 10000 job IDs required")
        ids = list(dict.fromkeys(identifier(value, "job ID") for value in job_ids))
        if not ids:
            return {}
        with self.connection() as conn:
            rows = conn.execute(
                "SELECT id,status FROM platform_jobs WHERE owner_id=%s AND id=ANY(%s)",
                (owner, ids),
            ).fetchall()
        result = {row["id"]: row["status"] for row in rows}
        if len(result) != len(ids):
            raise KeyError("Job not found")
        return result

    def list(self, *, owner="researcher", limit=100):
        with self.connection() as conn:
            return [
                self._public(r)
                for r in conn.execute(
                    "SELECT * FROM platform_jobs WHERE owner_id=%s ORDER BY created_at DESC LIMIT %s",
                    (owner, min(max(int(limit), 1), 500)),
                ).fetchall()
            ]

    def events(self, job_id, *, owner="researcher", after=0):
        self.get(job_id, owner=owner)
        with self.connection() as conn:
            return [
                self._public(r)
                for r in conn.execute(
                    "SELECT * FROM platform_job_events WHERE job_id=%s AND id>%s ORDER BY id LIMIT 500",
                    (job_id, int(after)),
                ).fetchall()
            ]

    @staticmethod
    def _host(telemetry):
        value = telemetry.get('physical_host_id') or telemetry.get('hostname')
        return value if isinstance(value, str) and 0 < len(value) <= 200 and not any(ord(c) < 32 for c in value) else None

    def _quarantines(self, conn):
        return [{"job_id": row['id'], **row['progress']['resource_quarantine']} for row in conn.execute(
            "SELECT id,progress FROM platform_jobs WHERE progress->'resource_quarantine'->>'active'='true'").fetchall()]

    def _admission_quarantined(self, conn, worker, *, gpu_uuid=None):
        conn.execute("SELECT pg_advisory_xact_lock(hashtext('skinpipe-resource-quarantine-v1'))")
        quarantines = self._quarantines(conn)
        if not quarantines:
            return False
        node = conn.execute('SELECT telemetry FROM platform_nodes WHERE id=%s', (worker,)).fetchone()
        host = self._host(node['telemetry']) if node else None
        return any(not item.get('physical_host') or not host or item['physical_host'] == host
            or (gpu_uuid is not None and gpu_uuid in item.get('gpu_uuids', [])) for item in quarantines)

    def _quarantine(self, conn, row, reason):
        if row['kind'] == 'probe' or not row['attempts']:
            return
        conn.execute("SELECT pg_advisory_xact_lock(hashtext('skinpipe-resource-quarantine-v1'))")
        node = conn.execute('SELECT telemetry FROM platform_nodes WHERE id=%s', (row.get('worker'),)).fetchone()
        host = self._host(node['telemetry']) if node else None
        gpu = conn.execute('SELECT gpu_uuid FROM platform_gpu_reservations WHERE job_id=%s', (row['id'],)).fetchall()
        quarantine = {'active': True, 'reason': reason, 'attempt': row['attempts'],
            'worker': row.get('worker'), 'physical_host': host,
            'host_identity_status': 'reported_host' if host else 'unknown_global_hold',
            'gpu_uuids': [item['gpu_uuid'] for item in gpu]}
        current = conn.execute('SELECT progress FROM platform_jobs WHERE id=%s', (row['id'],)).fetchone()['progress']
        conn.execute('UPDATE platform_jobs SET progress=%s WHERE id=%s',
            (Jsonb({**current, 'resource_quarantine': quarantine}), row['id']))
        self._event(conn, row['id'], 'resource_quarantined', quarantine)

    def claim(self, worker, kinds, *, engine_hashes=None):
        identifier(worker)
        if engine_hashes is not None:
            from .engines import validate_engine_hashes

            validate_engine_hashes(engine_hashes)
        if not kinds:
            return None
        with self.connection() as conn:
            expired = conn.execute("""SELECT * FROM platform_jobs WHERE
                status='running' AND lease_until < clock_timestamp()
                ORDER BY lease_until,id LIMIT 100 FOR UPDATE SKIP LOCKED""").fetchall()
            for old in expired:
                cleanup = old['progress'].get('resource_cleanup', {})
                cleanup_verified = cleanup.get('confirmed') is True and cleanup.get('attempt') == old['attempts']
                native = old['kind'] in {'dock_batch', 'discovery_prepare'} or old['kind'].startswith('algorithm_')
                unknown_native = native and not cleanup_verified
                if (old['cancel_requested'] or unknown_native) and not cleanup_verified:
                    self._quarantine(conn, old, 'cancelled_lease_expired_without_cleanup_ack'
                                     if old['cancel_requested'] else 'lease_expired_without_cleanup_ack')
                # An unreported interval may have performed computation. Charge
                # conservatively through lease expiry, not an unlimited outage.
                extra = conn.execute(
                    "SELECT GREATEST(0,EXTRACT(EPOCH FROM lease_until-updated_at)) AS seconds FROM platform_jobs WHERE id=%s",
                    (old["id"],),
                ).fetchone()["seconds"]
                conn.execute(
                    "UPDATE platform_compute_budgets SET used_seconds=used_seconds+%s,active_job=NULL WHERE id=%s AND active_job=%s",
                    (float(extra), old["budget_scope"], old["id"]),
                )
                conn.execute(
                    "UPDATE platform_jobs SET compute_seconds=compute_seconds+%s WHERE id=%s",
                    (float(extra), old["id"]),
                )
                conn.execute(
                    "UPDATE platform_attempts SET ended_at=clock_timestamp(),status='lease_expired' WHERE job_id=%s AND attempt=%s",
                    (old["id"], old["attempts"]),
                )
                status = (
                    "cancelled"
                    if old["cancel_requested"]
                    else "needs_data"
                    if old["budget_scope"] is None
                    else "failed"
                    if old["attempts"] >= old["max_attempts"]
                    else "retry_required"
                    if unknown_native
                    else "queued"
                )
                conn.execute(
                    "UPDATE platform_jobs SET status=%s,lease_token=NULL,lease_until=NULL,worker=NULL,"
                    "error=CASE WHEN %s THEN 'Native resource cleanup is unverified after lease expiry' ELSE error END,"
                    "updated_at=clock_timestamp() WHERE id=%s",
                    (status, unknown_native, old["id"]),
                )
                conn.execute(
                    "DELETE FROM platform_gpu_reservations WHERE job_id=%s",
                    (old["id"],),
                )
                self._event(
                    conn,
                    old["id"],
                    "lease_expired",
                    {"status": status, "attempt": old["attempts"]},
                )
            if self._admission_quarantined(conn, worker):
                return None
            unbound = conn.execute("""UPDATE platform_jobs SET status='needs_data',
                error='Legacy computation requires an administrator-bound shared budget',
                updated_at=clock_timestamp() WHERE status='queued' AND budget_scope IS NULL
                RETURNING id""").fetchall()
            for item in unbound:
                self._event(conn, item['id'], 'budget_binding_required')
            exhausted = conn.execute("""UPDATE platform_jobs j SET status='budget_deferred',
                progress=j.progress || '{"stage":"budget_deferred","reason":"Shared target budget exhausted"}'::jsonb,
                updated_at=clock_timestamp() FROM platform_compute_budgets b
                WHERE j.budget_scope=b.id AND j.status='queued' AND b.used_seconds>=b.limit_seconds
                RETURNING j.id""").fetchall()
            for item in exhausted:
                self._event(conn, item['id'], 'budget_deferred', {'reason':'shared_target_budget'})
            # Opt-in workers admit only the exact frozen pair. Omitted pins
            # retain the legacy claim contract; requests/results are never
            # rewritten to make an old cohort appear pinned.
            engine_clause = (
                " AND (j.kind<>'dock_batch' OR j.request->'engine_hashes'=%s::jsonb)"
                if engine_hashes is not None else ""
            )
            row = conn.execute(
                """SELECT j.* FROM platform_jobs j JOIN platform_compute_budgets b ON b.id=j.budget_scope
                WHERE j.status='queued' AND b.owner_id=j.owner_id
                AND b.active_job IS NULL AND b.used_seconds<b.limit_seconds
                AND NOT cancel_requested AND not_before<=clock_timestamp()
                AND kind=ANY(%s) AND compute_seconds<budget_seconds"""
                + engine_clause
                + " ORDER BY created_at,j.id FOR UPDATE OF j,b SKIP LOCKED LIMIT 1",
                (list(kinds), Jsonb(engine_hashes)) if engine_hashes is not None else (list(kinds),),
            ).fetchone()
            if not row:
                return None
            budget = conn.execute(
                "SELECT * FROM platform_compute_budgets WHERE id=%s",
                (row["budget_scope"],),
            ).fetchone()
            conn.execute(
                "UPDATE platform_compute_budgets SET active_job=%s WHERE id=%s",
                (row["id"], row["budget_scope"]),
            )
            conn.execute(
                "UPDATE platform_jobs SET budget_seconds=%s WHERE id=%s",
                (
                    row["compute_seconds"]
                    + budget["limit_seconds"]
                    - budget["used_seconds"],
                    row["id"],
                ),
            )
            token = uuid.uuid4().hex
            row = conn.execute(
                """UPDATE platform_jobs SET status='running',worker=%s,
                lease_token=%s,lease_until=clock_timestamp()+(%s * interval '1 second'),
                attempts=attempts+1,progress=progress-'resource_cleanup',updated_at=clock_timestamp() WHERE id=%s RETURNING *""",
                (worker, token, self.lease_seconds, row["id"]),
            ).fetchone()
            conn.execute(
                "INSERT INTO platform_attempts(job_id,attempt,token_sha256,worker) VALUES(%s,%s,%s,%s)",
                (
                    row["id"],
                    row["attempts"],
                    hashlib.sha256(token.encode()).hexdigest(),
                    worker,
                ),
            )
            self._event(
                conn,
                row["id"],
                "claimed",
                {"worker": worker, "attempt": row["attempts"]},
            )
            return {**self._public(row), "lease_token": token}

    def _owned(self, conn, job_id, token):
        # Lock first, then check the current database clock; a pre-lock WHERE
        # evaluation alone can accept a lease which expired while waiting.
        row = conn.execute(
            "SELECT * FROM platform_jobs WHERE id=%s FOR UPDATE", (job_id,)
        ).fetchone()
        if not row or row["status"] != "running" or row["lease_token"] != token:
            raise LeaseLost("Attempt no longer owns this job")
        live = conn.execute(
            "SELECT lease_until>clock_timestamp() AS live FROM platform_jobs WHERE id=%s",
            (job_id,),
        ).fetchone()["live"]
        if not live:
            raise LeaseLost("Attempt lease expired")
        return row

    def heartbeat(self, job_id, token, *, progress=None, compute_seconds=None):
        with self.connection() as conn:
            row = self._owned(conn, job_id, token)
            if row["cancel_requested"]:
                raise Cancelled("Cancellation requested")
            elapsed = (
                row["compute_seconds"]
                if compute_seconds is None
                else float(compute_seconds)
            )
            if not math.isfinite(elapsed) or elapsed < row["compute_seconds"]:
                raise ValueError("Compute time must be finite and monotonic")
            merged_progress = {**row['progress'], **(progress or {})}
            conn.execute('UPDATE platform_jobs SET progress=%s WHERE id=%s',
                         (Jsonb(merged_progress), job_id))
            if elapsed >= row["budget_seconds"]:
                self._close(
                    conn,
                    row,
                    "budget_deferred",
                    error="Actual compute budget exhausted",
                    elapsed=elapsed,
                )
                return {"status": "budget_deferred"}
            conn.execute(
                "UPDATE platform_compute_budgets SET used_seconds=used_seconds+%s WHERE id=%s",
                (elapsed - row["compute_seconds"], row["budget_scope"]),
            )
            conn.execute(
                """UPDATE platform_jobs SET lease_until=clock_timestamp()+(%s * interval '1 second'),
                compute_seconds=%s,progress=%s,updated_at=clock_timestamp() WHERE id=%s""",
                (
                    self.lease_seconds,
                    elapsed,
                    Jsonb(merged_progress),
                    job_id,
                ),
            )
            conn.execute(
                "UPDATE platform_gpu_reservations SET expires_at=clock_timestamp()+(%s * interval '1 second') WHERE job_id=%s AND lease_token=%s",
                (self.lease_seconds, job_id, token),
            )
            return {"status": "running", "cancel_requested": False}

    def _close(self, conn, row, status, *, result=None, error=None, elapsed=None):
        progress = conn.execute('SELECT progress FROM platform_jobs WHERE id=%s', (row['id'],)).fetchone()['progress']
        cleanup = progress.get('resource_cleanup', {})
        verified = cleanup.get('confirmed') is True and cleanup.get('attempt') == row['attempts']
        if error == 'owned_resource_cleanup_unverified' or (status in {'cancelled', 'budget_deferred'} and not verified):
            self._quarantine(conn, row, error or status + '_without_cleanup_ack')
        used = row["compute_seconds"] if elapsed is None else elapsed
        conn.execute(
            "UPDATE platform_compute_budgets SET used_seconds=used_seconds+%s,active_job=NULL WHERE id=%s AND active_job=%s",
            (used - row["compute_seconds"], row["budget_scope"], row["id"]),
        )
        conn.execute(
            """UPDATE platform_jobs SET status=%s,result=%s,error=%s,
            compute_seconds=%s,lease_token=NULL,lease_until=NULL,updated_at=clock_timestamp()
            WHERE id=%s""",
            (
                status,
                Jsonb(result),
                error,
                row["compute_seconds"] if elapsed is None else elapsed,
                row["id"],
            ),
        )
        conn.execute(
            "UPDATE platform_attempts SET ended_at=clock_timestamp(),status=%s WHERE job_id=%s AND attempt=%s",
            (status, row["id"], row["attempts"]),
        )
        conn.execute(
            "DELETE FROM platform_gpu_reservations WHERE job_id=%s", (row["id"],)
        )
        conn.execute("DELETE FROM platform_cpu_reservations WHERE job_id=%s", (row["id"],))
        self._event(
            conn, row["id"], status, {"attempt": row["attempts"], "error": error}
        )

    def finish(self, job_id, token, result, *, compute_seconds=None, cleanup=None):
        if cleanup is not None and (not isinstance(cleanup, dict)
                or cleanup.get("owned_processes_stopped") is not True
                or cleanup.get("scope") != "owned_job_processes"):
            raise ValueError("Verified owned process cleanup receipt is required")
        signature = digest(result)
        with self.connection() as conn:
            # Lock before inspecting publication: concurrent retries must see
            # the committed receipt after the first finisher releases this row.
            conn.execute("SELECT id FROM platform_jobs WHERE id=%s FOR UPDATE", (job_id,))
            published = conn.execute(
                "SELECT result_sha256 FROM platform_publications WHERE job_id=%s",
                (job_id,),
            ).fetchone()
            if published:
                if published["result_sha256"] != signature:
                    raise ValueError("Job already published a different result")
                return {"status": "complete", "replayed": True}
            row = self._owned(conn, job_id, token)
            elapsed = (
                row["compute_seconds"]
                if compute_seconds is None
                else float(compute_seconds)
            )
            if not math.isfinite(elapsed) or elapsed < row["compute_seconds"]:
                raise ValueError("Compute time must be finite and monotonic")
            if cleanup is not None:
                receipt = {**cleanup, "confirmed": True, "attempt": row["attempts"]}
                conn.execute("UPDATE platform_jobs SET progress=%s WHERE id=%s",
                    (Jsonb({**row["progress"], "resource_cleanup": receipt}), job_id))
                self._event(conn, job_id, "owned_cleanup_confirmed", receipt)
            if row["cancel_requested"]:
                self._close(conn, row, "cancelled", elapsed=elapsed)
                return {"status": "cancelled"}
            # A completed artifact within the limit is committed exactly once.
            if elapsed > row["budget_seconds"]:
                self._close(conn, row, "budget_deferred", elapsed=elapsed)
                return {"status": "budget_deferred"}
            conn.execute(
                "INSERT INTO platform_publications(job_id,result_sha256,result) VALUES(%s,%s,%s)",
                (job_id, signature, Jsonb(result)),
            )
            self._close(conn, row, "complete", result=result, elapsed=elapsed)
            return {"status": "complete", "replayed": False}

    def fail(self, job_id, token, error, *, status="failed"):
        if status not in {"failed", "cancelled", "budget_deferred", "needs_data", "retryable"}:
            raise ValueError("Invalid terminal status")
        with self.connection() as conn:
            row = self._owned(conn, job_id, token)
            retry = status == "retryable" and row["attempts"] < row["max_attempts"] and not row["cancel_requested"]
            self._close(
                conn,
                row,
                "cancelled" if row["cancel_requested"] else "queued" if retry else "failed" if status == "retryable" else status,
                error=str(error)[:2000],
            )
            if retry:
                conn.execute("UPDATE platform_jobs SET not_before=clock_timestamp()+(30 * interval '1 second') WHERE id=%s", (job_id,))
                self._event(conn, job_id, "infrastructure_retry", {"attempt": row["attempts"], "max_attempts": row["max_attempts"]})

    def defer(self, job_id, token, reason, *, delay=30):
        if not 1 <= delay <= 3600:
            raise ValueError("Invalid deferral duration")
        with self.connection() as conn:
            row = self._owned(conn, job_id, token)
            if row["cancel_requested"]:
                self._close(conn, row, "cancelled")
                return
            self._close(conn, row, "queued")
            # Keep monotonically numbered attempts, but resource waits do not
            # consume failure retry allowance.
            conn.execute(
                "UPDATE platform_jobs SET max_attempts=max_attempts+1,not_before=clock_timestamp()+(%s * interval '1 second'),progress=%s WHERE id=%s",
                (
                    delay,
                    Jsonb(
                        {
                            **row['progress'],
                            "stage": "resource_wait",
                            "reason": reason,
                            "queue_eta_seconds": None,
                        }
                    ),
                    job_id,
                ),
            )
            self._event(conn, job_id, "resource_deferred", {"reason": reason})

    def cancel_if_queued(self, job_id, *, owner="researcher", campaign_id, stage, reason):
        'Validate the selected structure and retain explicit calculation provenance.'
        with self.connection() as conn:
            row = conn.execute(
                "SELECT * FROM platform_jobs WHERE id=%s AND owner_id=%s FOR UPDATE",
                (job_id, owner),
            ).fetchone()
            if not row:
                raise KeyError("Job not found")
            if (row["request"].get("campaign_id") != campaign_id
                    or row["request"].get("stage") != stage):
                raise ValueError("Queued retirement campaign/stage mismatch")
            if row["status"] != "queued":
                return self._public(row)
            conn.execute(
                "UPDATE platform_jobs SET status='cancelled',cancel_requested=true,error=%s,updated_at=clock_timestamp() WHERE id=%s",
                (str(reason)[:2000], job_id),
            )
            self._event(conn, job_id, "queued_policy_retirement", {"reason": str(reason)[:2000]})
            return self._public(conn.execute("SELECT * FROM platform_jobs WHERE id=%s", (job_id,)).fetchone())

    def cancel(self, job_id, *, owner="researcher"):
        with self.connection() as conn:
            row = conn.execute(
                "SELECT * FROM platform_jobs WHERE id=%s AND owner_id=%s FOR UPDATE",
                (job_id, owner),
            ).fetchone()
            if not row:
                raise KeyError("Job not found")
            if row["status"] in {"complete", "failed", "cancelled"}:
                return self._public(row)
            conn.execute(
                "UPDATE platform_jobs SET cancel_requested=true,status=CASE WHEN status='running' THEN status ELSE 'cancelled' END,updated_at=clock_timestamp() WHERE id=%s",
                (job_id,),
            )
            self._event(conn, job_id, "cancel_requested")
        return self.get(job_id, owner=owner)

    def acknowledge_cancel(self, job_id, token, *, compute_seconds=None, cleanup=None):
        """Close a current cancelled lease only after its actor reaped owned work.

        An expired lease is intentionally not proof of native-process exit.
        The receipt is retained in progress so restart/reconciliation can
        distinguish an acknowledged stop from an unverified timeout.
        """
        if (not isinstance(cleanup, dict) or cleanup.get("owned_processes_stopped") is not True
                or cleanup.get("scope") != "owned_job_processes"):
            raise ValueError("Verified owned process cleanup receipt is required")
        with self.connection() as conn:
            row = self._owned(conn, job_id, token)
            if not row["cancel_requested"]:
                raise ValueError("Job cancellation was not requested")
            elapsed = row["compute_seconds"] if compute_seconds is None else float(compute_seconds)
            if not math.isfinite(elapsed) or elapsed < row["compute_seconds"]:
                raise ValueError("Compute time must be finite and monotonic")
            receipt = {**cleanup, "confirmed": True, "attempt": row["attempts"]}
            conn.execute("UPDATE platform_jobs SET progress=%s WHERE id=%s",
                         (Jsonb({**row["progress"], "resource_cleanup": receipt}), job_id))
            self._event(conn, job_id, "owned_cleanup_confirmed", receipt)
            self._close(conn, row, "cancelled", elapsed=elapsed)
        return {"status": "cancelled", "resource_cleanup": receipt}

    def record_cleanup_diagnostic(self, job_id, token, *, worker, stage, error_type):
        """Retain only a fixed diagnostic for the authenticated current attempt.

        A closed budget lease can still report why cleanup was unverified.
        This event neither certifies cleanup nor changes any job or budget.
        """
        allowed = {'owned_cleanup': {'OwnedProcessCleanupError'},
                   'cleanup_ack': {'LeaseLost', 'ValueError'}}
        if (not isinstance(stage, str) or not isinstance(error_type, str)
                or error_type not in allowed.get(stage, set())):
            raise ValueError('Unsupported cleanup diagnostic category')
        identifier(worker, 'worker')
        if not isinstance(token, str) or not 1 <= len(token) <= 200:
            raise LeaseLost('Cleanup diagnostic does not own this current attempt')
        with self.connection() as conn:
            row = conn.execute('SELECT * FROM platform_jobs WHERE id=%s FOR UPDATE', (job_id,)).fetchone()
            if (not row or row['status'] == 'queued'
                    or (row['status'] == 'running' and row['lease_token'] != token)):
                raise LeaseLost('Cleanup diagnostic does not own this current attempt')
            attempt = conn.execute('SELECT token_sha256,worker FROM platform_attempts WHERE job_id=%s AND attempt=%s',
                (job_id, row['attempts'])).fetchone()
            if (not attempt or attempt['worker'] != worker
                    or attempt['token_sha256'] != hashlib.sha256(token.encode()).hexdigest()):
                raise LeaseLost('Cleanup diagnostic does not own this current attempt')
            self._event(conn, job_id, 'owned_cleanup_diagnostic', {
                'attempt': row['attempts'], 'worker': worker, 'stage': stage, 'error_type': error_type})
            return {'status': 'cleanup_diagnostic_recorded', 'attempt': row['attempts']}

    def acknowledge_cleanup(self, job_id, token, *, worker, cleanup):
        """Record stopped-actor proof or clear its matching terminal quarantine.

        Expiry proves neither exit nor release. The immutable attempt token
        and worker binding permit a late cleanup proof without reopening its
        lease, publishing scientific results, or refunding compute time.
        """
        if (not isinstance(cleanup, dict) or cleanup.get('scope') != 'owned_job_processes'
                or cleanup.get('owned_processes_stopped') is not True):
            raise ValueError('Verified owned process cleanup receipt is required')
        with self.connection() as conn:
            row = conn.execute('SELECT * FROM platform_jobs WHERE id=%s FOR UPDATE', (job_id,)).fetchone()
            if not row or row['status'] == 'queued':
                raise LeaseLost('Only this stopped attempt may acknowledge cleanup')
            pending_reclamation = row['status'] == 'running'
            if pending_reclamation and row['lease_token'] != token:
                raise LeaseLost('Only this stopped attempt may acknowledge cleanup')
            attempt = conn.execute('SELECT token_sha256,worker FROM platform_attempts WHERE job_id=%s AND attempt=%s',
                (job_id, row['attempts'])).fetchone()
            if (not attempt or attempt['worker'] != worker
                    or attempt['token_sha256'] != hashlib.sha256(token.encode()).hexdigest()):
                raise LeaseLost('Cleanup proof does not own this current attempt')
            prior = row['progress'].get('resource_cleanup', {})
            quarantine = row['progress'].get('resource_quarantine', {})
            if prior.get('confirmed') is True and prior.get('attempt') == row['attempts'] and not quarantine.get('active'):
                return {'status': 'cleanup_confirmed', 'replayed': True}
            if pending_reclamation:
                # A stopped actor can report before its lease expires or the
                # next claimant retires it. Save its exact-attempt proof while
                # retaining status, lease, reservations and budget; expiry
                # still fences publication and decides whether to retry.
                receipt = {**cleanup, 'confirmed': True, 'attempt': row['attempts']}
                conn.execute('UPDATE platform_jobs SET progress=%s WHERE id=%s',
                    (Jsonb({**row['progress'], 'resource_cleanup': receipt}), job_id))
                self._event(conn, job_id, 'attempt_cleanup_confirmed', {'attempt': row['attempts'], 'worker': worker})
                return {'status': 'cleanup_confirmed', 'replayed': False, 'resource_cleanup': receipt}
            if quarantine.get('active') is not True or quarantine.get('attempt') != row['attempts']:
                raise ValueError('No matching current resource quarantine')
            receipt = {**cleanup, 'confirmed': True, 'attempt': row['attempts']}
            progress = {**row['progress'], 'resource_cleanup': receipt,
                'resource_quarantine': {**quarantine, 'active': False, 'cleared_by': 'authenticated_attempt_cleanup_ack'}}
            conn.execute('UPDATE platform_jobs SET progress=%s,updated_at=clock_timestamp() WHERE id=%s', (Jsonb(progress), job_id))
            self._event(conn, job_id, 'resource_quarantine_cleared', {'attempt': row['attempts'], 'worker': worker})
            return {'status': 'cleanup_confirmed', 'replayed': False, 'resource_cleanup': receipt}

    def campaign_job_ids(self, campaign_id, *, owner='researcher'):
        'Validate the selected structure and retain explicit calculation provenance.'
        identifier(campaign_id, 'campaign')
        with self.connection() as conn:
            return [row['id'] for row in conn.execute(
                "SELECT id FROM platform_jobs WHERE owner_id=%s AND request->>'campaign_id'=%s "
                "AND status NOT IN ('complete','cancelled') ORDER BY created_at,id",
                (owner, campaign_id),
            ).fetchall()]

    def resume(self, job_id, *, owner="researcher", budget_seconds=None):
        with self.connection() as conn:
            row = conn.execute(
                "SELECT * FROM platform_jobs WHERE id=%s AND owner_id=%s FOR UPDATE",
                (job_id, owner),
            ).fetchone()
            if not row:
                raise KeyError("Job not found")
            if row['progress'].get('resource_quarantine', {}).get('active') is True:
                raise ValueError('Verified owned resource cleanup is required before resume')
            if row["status"] not in {
                "failed",
                "cancelled",
                "needs_data",
                "budget_deferred",
                "retry_required",
            }:
                raise ValueError("Only stopped jobs may resume")
            budget = (
                row["budget_seconds"]
                if budget_seconds is None
                else float(budget_seconds)
            )
            if (
                not math.isfinite(budget)
                or not 0 < budget <= 604800
            ):
                raise ValueError("Resume requires remaining compute budget")
            if row['compute_seconds'] >= budget:
                raise BudgetExhausted("Resume requires remaining compute budget")
            shared = conn.execute(
                "SELECT * FROM platform_compute_budgets WHERE id=%s FOR UPDATE",
                (row["budget_scope"],),
            ).fetchone()
            if shared is None:
                raise ValueError("Legacy computation requires an administrator-bound shared budget")
            if shared["owner_id"] != row["owner_id"]:
                raise ValueError("Compute budget belongs to another owner")
            if (
                budget_seconds is None
                and shared["used_seconds"] >= shared["limit_seconds"]
            ):
                raise BudgetExhausted("Resume requires increased target compute budget")
            if budget_seconds is not None:
                if budget <= shared["used_seconds"]:
                    raise BudgetExhausted("Target compute budget is already consumed")
                if budget < shared["limit_seconds"]:
                    raise ValueError("Resume cannot reduce an existing shared budget")
                conn.execute(
                    "UPDATE platform_compute_budgets SET limit_seconds=%s WHERE id=%s",
                    (budget, row["budget_scope"]),
                )
            conn.execute(
                """UPDATE platform_jobs SET status='queued',cancel_requested=false,
                error=NULL,lease_token=NULL,lease_until=NULL,max_attempts=attempts+3,progress=progress-'resource_cleanup',
                budget_seconds=%s,not_before=clock_timestamp(),updated_at=clock_timestamp() WHERE id=%s""",
                (budget, job_id),
            )
            self._event(conn, job_id, "resumed", {"budget_seconds": budget})
        return self.get(job_id, owner=owner)

    def register_node(self, node_id, telemetry):
        identifier(node_id)
        if not isinstance(telemetry, dict):
            raise ValueError("Telemetry must be an object")
        with self.connection() as conn:
            conn.execute(
                "INSERT INTO platform_nodes(id,telemetry) VALUES(%s,%s) ON CONFLICT(id) DO UPDATE SET telemetry=EXCLUDED.telemetry,updated_at=clock_timestamp()",
                (node_id, Jsonb(telemetry)),
            )

    def resources(self):
        with self.connection() as conn:
            nodes = [
                self._public(r)
                for r in conn.execute(
                    "SELECT *,updated_at>clock_timestamp()-interval '90 seconds' AS fresh FROM platform_nodes ORDER BY id"
                ).fetchall()
            ]
            reservations = [
                self._public(r)
                for r in conn.execute(
                    "SELECT gpu_uuid,job_id,expires_at FROM platform_gpu_reservations WHERE expires_at>clock_timestamp()"
                ).fetchall()
            ]
        with self.connection() as conn:
            cpu = [self._public(r) for r in conn.execute("""SELECT r.node_id,r.physical_host,r.job_id,r.cores FROM platform_cpu_reservations r
                JOIN platform_jobs j ON j.id=r.job_id WHERE j.status='running' AND j.lease_token=r.lease_token AND j.lease_until>clock_timestamp()""").fetchall()]
            quarantines = self._quarantines(conn)
        return {"nodes": nodes, "reservations": reservations,"cpu_reservations":cpu,
                'quarantines': quarantines, 'resource_release_unverified': bool(quarantines)}

    def reserve_gpu(self, job_id, token, gpu_uuid):
        identifier(gpu_uuid, "GPU UUID")
        with self.connection() as conn:
            row = self._owned(conn, job_id, token)
            if self._admission_quarantined(conn, row['worker'], gpu_uuid=gpu_uuid):
                return False
            existing = conn.execute("SELECT * FROM platform_gpu_reservations WHERE gpu_uuid=%s", (gpu_uuid,)).fetchone()
            if existing and existing["job_id"] == job_id and existing["lease_token"] == token:
                return True
            if existing:
                previous = conn.execute('SELECT kind,status,lease_token,progress,attempts FROM platform_jobs WHERE id=%s',
                                        (existing['job_id'],)).fetchone()
                if previous and previous['status'] == 'running' and previous['lease_token'] == existing['lease_token']:
                    native = previous['kind'] in {'dock_batch', 'discovery_prepare'} or previous['kind'].startswith('algorithm_')
                    cleanup = previous['progress'].get('resource_cleanup', {})
                    if native and not (cleanup.get('confirmed') is True and cleanup.get('attempt') == previous['attempts']):
                        # Another already-claimed job may reserve before a
                        # claimant sweeps expired leases into quarantine.
                        return False
            conn.execute(
                "DELETE FROM platform_gpu_reservations WHERE gpu_uuid=%s AND expires_at<clock_timestamp()",
                (gpu_uuid,),
            )
            return (
                conn.execute(
                    """INSERT INTO platform_gpu_reservations(gpu_uuid,job_id,lease_token,expires_at)
                VALUES(%s,%s,%s,%s) ON CONFLICT(gpu_uuid) DO NOTHING RETURNING gpu_uuid""",
                    (gpu_uuid, job_id, token, row["lease_until"]),
                ).fetchone()
                is not None
            )

    def reserve_cpu(self, job_id, token, cores):
        """Fence each job by worker capacity and aggregate use by host capacity."""
        if type(cores) is not int or not 1 <= cores <= 8:
            raise ValueError("CPU allocation must be 1..8 cores")
        with self.connection() as conn:
            job = self._owned(conn, job_id, token)
            if self._admission_quarantined(conn, job['worker']):
                return False
            node = conn.execute("SELECT * FROM platform_nodes WHERE id=%s AND updated_at>clock_timestamp()-interval '90 seconds' FOR UPDATE", (job["worker"],)).fetchone()
            if not node:
                return False
            limits = _cpu_reservation_limits(node["telemetry"])
            if limits is None or cores > limits[0]:
                return False
            _, host_limit = limits
            physical = node['telemetry'].get('physical_host_id') or node['telemetry'].get('hostname') or job['worker']
            identifier(physical,'physical host')
            conn.execute('SELECT pg_advisory_xact_lock(hashtext(%s))',('cpu:'+physical,))
            occupied = conn.execute("""SELECT coalesce(sum(r.cores),0) AS cores FROM platform_cpu_reservations r
                JOIN platform_jobs j ON j.id=r.job_id WHERE coalesce(r.physical_host,r.node_id)=%s AND r.job_id<>%s
                AND j.status='running' AND j.lease_token=r.lease_token AND
                (j.lease_until>clock_timestamp() OR
                 ((j.kind IN ('dock_batch','discovery_prepare') OR left(j.kind,10)='algorithm_')
                  AND NOT coalesce(j.progress->'resource_cleanup'->>'confirmed'='true'
                      AND j.progress->'resource_cleanup'->>'attempt'=j.attempts::text,false)))""", (physical, job_id)).fetchone()["cores"]
            if occupied + cores > host_limit:
                return False
            conn.execute("""INSERT INTO platform_cpu_reservations(job_id,node_id,lease_token,cores,physical_host) VALUES(%s,%s,%s,%s,%s)
                ON CONFLICT(job_id) DO UPDATE SET node_id=excluded.node_id,lease_token=excluded.lease_token,cores=excluded.cores,physical_host=excluded.physical_host""", (job_id, job["worker"], token, cores,physical))
            return True

    def reserve_model_launch(self, job_id, token, work_key):
        """Count physical Boltz/Chai attempts before spawn, including failed retries.

        The scope and limit come from the sealed job, never from worker arguments.
        Retrying an ambiguous reservation in the same attempt is idempotent.
        """
        identifier(work_key, "model work key")
        with self.connection() as conn:
            job = self._owned(conn, job_id, token)
            request = job["request"]
            if job["kind"] not in {"algorithm_boltz2", "algorithm_chai1"}:
                raise ValueError("Launch reservation only applies to complex models")
            from .algorithm_profiles import validate_profile
            profile = validate_profile(request["compute_profile"])
            scope = "complex:" + digest({"owner": job["owner_id"], "campaign": request["campaign_id"], "profile": profile["sha256"]})
            maximum = profile["counts"]["complex_models"]
            conn.execute("INSERT INTO platform_launch_budgets(scope,launch_limit) VALUES(%s,%s) ON CONFLICT DO NOTHING", (scope, maximum))
            budget = conn.execute("SELECT * FROM platform_launch_budgets WHERE scope=%s FOR UPDATE", (scope,)).fetchone()
            if budget["launch_limit"] != maximum:
                raise ValueError("Frozen launch budget differs from request")
            prior = conn.execute("SELECT 1 FROM platform_model_launches WHERE scope=%s AND job_id=%s AND attempt=%s AND work_key=%s", (scope, job_id, job["attempts"], work_key)).fetchone()
            if prior:
                return True
            if budget["used"] >= maximum:
                return False
            conn.execute("INSERT INTO platform_model_launches(scope,job_id,attempt,work_key) VALUES(%s,%s,%s,%s)", (scope, job_id, job["attempts"], work_key))
            conn.execute("UPDATE platform_launch_budgets SET used=used+1 WHERE scope=%s", (scope,))
            self._event(conn, job_id, "physical_model_launch_reserved", {"work_key": work_key, "used": budget["used"]+1, "limit": maximum})
            return True

    def start_model_launch(self, job_id, token, work_key):
        """One physical spawn authorization; ambiguous retries never spawn twice."""
        identifier(work_key,'model work key')
        with self.connection() as conn:
            job = self._owned(conn,job_id,token)
            return conn.execute('''UPDATE platform_model_launches SET spawn_authorized=true
                WHERE job_id=%s AND attempt=%s AND work_key=%s AND spawn_authorized=false
                RETURNING scope''', (job_id,job['attempts'],work_key)).fetchone() is not None
