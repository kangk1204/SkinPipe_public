"""Shutdown must fence new coordinator work before preparation starts."""
import threading
from types import SimpleNamespace

import pytest

from skinpipe.platform.discovery import DiscoveryRuns


@pytest.mark.parametrize("export_only", [False, True])
def test_stopped_coordinator_does_not_probe_register_or_claim(export_only):
    runs = DiscoveryRuns.__new__(DiscoveryRuns)
    runs.service, runs.queue = object(), object()
    stopped = threading.Event()
    stopped.set()
    assert runs.run_once(shutdown_event=stopped, export_only=export_only) is None


@pytest.mark.parametrize("export_only", [False, True])
def test_shutdown_during_claim_defers_without_entering_preparation(tmp_path, monkeypatch, export_only):
    stopped = threading.Event()
    calls = []

    class Queue:
        def register_node(self, *args):
            pass

        def claim(self, worker, capabilities):
            stopped.set()
            return {"id": "fixture-job", "lease_token": "fixture-token", "compute_seconds": 0,
                    "budget_seconds": 60, "owner_id": "fixture-owner", "request": {},
                    "kind": "result_export" if export_only else "discovery_prepare"}

        def reserve_cpu(self, *args):
            return True

        def heartbeat(self, *args, **kwargs):
            return {"status": "running"}

        def defer(self, *args):
            calls.append("deferred")

        def fail(self, *args, **kwargs):
            calls.append("failed")

        def acknowledge_cleanup(self, *args, **kwargs):
            calls.append("cleanup_ack")

    def execute(*args, **kwargs):
        calls.append("entered_preparation")
        return {"synthetic": True}

    runs = DiscoveryRuns.__new__(DiscoveryRuns)
    runs.queue, runs.execute = Queue(), execute
    runs.service = SimpleNamespace(root=tmp_path, exports=SimpleNamespace(execute=execute))
    monkeypatch.setattr("skinpipe.platform.resources.probe", lambda **_: {})
    assert runs.run_once(shutdown_event=stopped, export_only=export_only) == {
        "status": "resource_deferred", "reason": "personal_service_shutdown"}
    assert calls == ["deferred", "cleanup_ack"]
