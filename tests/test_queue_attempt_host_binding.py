"""Current registration must not relocate already-owned native resources."""
import json
import os
from pathlib import Path
import uuid

import psycopg
from psycopg.conninfo import make_conninfo
import pytest

from skinpipe.platform.queue import Queue


@pytest.fixture
def queue():
    dsn = os.environ.get("SKINPIPE_RUNTIME_AUDIT_DSN") or os.environ.get("SKINPIPE_TEST_PG_DSN")
    config = os.environ.get("SKINPIPE_TEST_POSTGRES_CONFIG")
    if not dsn and config:
        dsn = make_conninfo(**json.loads(Path(config).read_text()))
    if not dsn:
        pytest.skip("Explicit isolated PostgreSQL test DSN required")
    with psycopg.connect(dsn) as conn:
        assert conn.execute("SELECT current_database()").fetchone()[0].endswith("_test")
    queue = Queue(dsn, lease_seconds=30)
    def clear_fixture_state():
        with queue.connection() as conn:
            conn.execute("TRUNCATE platform_jobs,platform_attempts,platform_job_events,platform_publications,"
                         "platform_nodes,platform_gpu_reservations,platform_compute_budgets RESTART IDENTITY CASCADE")
    clear_fixture_state()
    try:
        yield queue
    finally:
        # These tests deliberately retain unknown-ownership holds. Do not leak
        # those fixture-only holds into the next module's independent jobs.
        clear_fixture_state()


def register(queue, worker, host):
    queue.register_node(worker, {"physical_host_id": host, "cpu_count": 8,
                                 "cpu_available": 8, "cpu_limit": 2})


def enqueue(queue):
    return queue.enqueue("discovery_prepare", {}, idempotency_key=uuid.uuid4().hex,
                         budget_seconds=60)


@pytest.mark.parametrize("reserved", [False, True])
@pytest.mark.parametrize("termination", ["failed_cleanup", "lease_expired"])
def test_reregistered_worker_cannot_move_previous_attempt_quarantine(queue, reserved, termination):
    register(queue, "coordinator", "host-a")
    register(queue, "host-a-alias", "host-a")
    register(queue, "host-b-worker", "host-b")
    target = enqueue(queue)
    claim = queue.claim("coordinator", ["discovery_prepare"])
    if reserved:
        assert queue.reserve_cpu(target["id"], claim["lease_token"], 1)
    register(queue, "coordinator", "host-b")
    if termination == "failed_cleanup":
        queue.fail(target["id"], claim["lease_token"], "owned_resource_cleanup_unverified", status="needs_data")
    else:
        with queue.connection() as conn:
            conn.execute("UPDATE platform_jobs SET lease_until=clock_timestamp()-interval '1 second' WHERE id=%s", (target["id"],))
        assert queue.claim("host-a-alias", ["discovery_prepare"]) is None
    hold = queue.resources()["quarantines"][0]
    assert hold["physical_host"] == "host-a"
    enqueue(queue)
    assert queue.claim("host-a-alias", ["discovery_prepare"]) is None
    assert queue.claim("host-b-worker", ["discovery_prepare"]) is not None


def test_registration_cannot_supply_host_proof_for_previously_unknown_attempt(queue):
    target = enqueue(queue)
    claim = queue.claim("unknown-coordinator", ["discovery_prepare"])
    register(queue, "unknown-coordinator", "host-b")
    queue.fail(target["id"], claim["lease_token"], "owned_resource_cleanup_unverified", status="needs_data")
    hold = queue.resources()["quarantines"][0]
    assert hold["physical_host"] is None
    assert hold["host_identity_status"] == "unknown_global_hold"


def test_current_attempt_reservations_do_not_migrate_after_reregistration(queue):
    register(queue, "coordinator", "host-a")
    target = enqueue(queue)
    claim = queue.claim("coordinator", ["discovery_prepare"])
    assert queue.reserve_cpu(target["id"], claim["lease_token"], 1)
    assert queue.reserve_gpu(target["id"], claim["lease_token"], "GPU-a")
    register(queue, "coordinator", "host-b")
    assert queue.reserve_cpu(target["id"], claim["lease_token"], 1) is False
    assert queue.reserve_gpu(target["id"], claim["lease_token"], "GPU-b") is False
    with queue.connection() as conn:
        cpu = conn.execute("SELECT physical_host FROM platform_cpu_reservations WHERE job_id=%s", (target["id"],)).fetchone()
        gpu = conn.execute("SELECT gpu_uuid FROM platform_gpu_reservations WHERE job_id=%s", (target["id"],)).fetchone()
    assert cpu["physical_host"] == "host-a" and gpu["gpu_uuid"] == "GPU-a"


def test_legacy_attempt_without_frozen_host_stays_conservatively_unknown(queue):
    register(queue, "coordinator", "host-a")
    target = enqueue(queue)
    claim = queue.claim("coordinator", ["discovery_prepare"])
    with queue.connection() as conn:
        conn.execute("ALTER TABLE platform_attempts DROP COLUMN physical_host")
    Queue(queue.dsn)  # Upgrade a legacy ledger without inventing its old host.
    queue.fail(target["id"], claim["lease_token"], "owned_resource_cleanup_unverified", status="needs_data")
    assert queue.resources()["quarantines"][0]["host_identity_status"] == "unknown_global_hold"


def test_schema_upgrade_preserves_frozen_attempt_host_on_repeated_initialization(queue):
    register(queue, "coordinator", "host-a")
    target = enqueue(queue)
    queue.claim("coordinator", ["discovery_prepare"])
    register(queue, "coordinator", "host-b")
    Queue(queue.dsn)
    Queue(queue.dsn)
    with queue.connection() as conn:
        attempt = conn.execute("SELECT physical_host FROM platform_attempts WHERE job_id=%s", (target["id"],)).fetchone()
    assert attempt["physical_host"] == "host-a"


def test_registration_waits_until_claim_admission_and_host_snapshot_are_atomic(queue, monkeypatch):
    from concurrent.futures import ThreadPoolExecutor
    import threading
    import time

    register(queue, "coordinator", "host-a")
    register(queue, "held-worker", "host-b")
    held = enqueue(queue)
    held_claim = queue.claim("held-worker", ["discovery_prepare"])
    queue.fail(held["id"], held_claim["lease_token"], "owned_resource_cleanup_unverified", status="needs_data")
    pending = enqueue(queue)
    entered, release = threading.Event(), threading.Event()
    original = queue._quarantines

    def during_admission(conn):
        entered.set()
        assert release.wait(4), "Fixture admission barrier timed out"
        return original(conn)

    monkeypatch.setattr(queue, "_quarantines", during_admission)
    application = "host-race-" + uuid.uuid4().hex
    registrar = Queue(make_conninfo(queue.dsn, application_name=application), initialize=False)
    with ThreadPoolExecutor(2) as pool:
        claim = pool.submit(queue.claim, "coordinator", ["discovery_prepare"])
        registration = None
        try:
            assert entered.wait(2)
            registration = pool.submit(register, registrar, "coordinator", "host-b")
            blocked = False
            deadline = time.monotonic() + 2
            while time.monotonic() < deadline and not registration.done():
                with queue.connection() as conn:
                    activity = conn.execute("SELECT wait_event_type FROM pg_stat_activity WHERE application_name=%s", (application,)).fetchone()
                if activity and activity["wait_event_type"] == "Lock":
                    blocked = True
                    break
                time.sleep(.01)
            assert blocked, "Re-registration changed host before admission and snapshot committed"
        finally:
            release.set()
        owned = claim.result(timeout=3)
        registration.result(timeout=3)
    assert owned["id"] == pending["id"]
    with queue.connection() as conn:
        attempt = conn.execute("SELECT physical_host FROM platform_attempts WHERE job_id=%s", (pending["id"],)).fetchone()
    assert attempt["physical_host"] == "host-a"
