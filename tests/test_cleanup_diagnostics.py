"""Cleanup diagnostics preserve terminal holds; subprocesses are harmless fixtures."""
from contextlib import contextmanager
from copy import deepcopy
import hashlib
import subprocess
import sys
import time
from types import SimpleNamespace

import pytest

from skinpipe.platform import engines, research_docking
from skinpipe.platform.artifacts import ArtifactStore
from skinpipe.platform.discovery import DiscoveryRuns
from skinpipe.platform.queue import LeaseLost, Queue


class DiagnosticConnection:
    def __init__(self, row, attempt):
        self.row, self.attempt = row, attempt
        self.statements = []

    def execute(self, statement, values):
        self.statements.append(statement)
        assert statement.startswith("SELECT "), "Diagnostic must not mutate job state"
        value = self.attempt if "platform_attempts" in statement else self.row
        return SimpleNamespace(fetchone=lambda: deepcopy(value))


def diagnostic_queue():
    row = {"id": "fixture-job", "status": "budget_deferred", "attempts": 2,
           "lease_token": None, "progress": {"resource_quarantine": {"active": True}},
           "compute_seconds": 1, "result": None}
    attempt = {"worker": "fixture-worker", "token_sha256": hashlib.sha256(b"fixture-token").hexdigest()}
    connection = DiagnosticConnection(row, attempt)
    queue = object.__new__(Queue)
    @contextmanager
    def connect():
        yield connection
    queue.connection = connect
    queue.recorded_events = []
    queue._event = lambda conn, job, event, payload: queue.recorded_events.append((job, event, payload))
    return queue, connection


def test_terminal_diagnostic_authenticates_without_changing_job_or_quarantine():
    queue, connection = diagnostic_queue()
    before = deepcopy(connection.row)
    result = queue.record_cleanup_diagnostic("fixture-job", "fixture-token", worker="fixture-worker",
        stage="owned_cleanup", error_type="OwnedProcessCleanupError")
    assert result == {"status": "cleanup_diagnostic_recorded", "attempt": 2}
    assert connection.row == before
    assert queue.recorded_events == [("fixture-job", "owned_cleanup_diagnostic", {
        "attempt": 2, "worker": "fixture-worker", "stage": "owned_cleanup", "error_type": "OwnedProcessCleanupError"})]
    assert "fixture-token" not in repr(queue.recorded_events)


@pytest.mark.parametrize("fault", ["wrong_token", "wrong_worker", "stale_attempt", "queued", "running_token"])
def test_cleanup_diagnostic_rejects_non_current_actor(fault):
    queue, connection = diagnostic_queue()
    token, worker = "fixture-token", "fixture-worker"
    if fault == "wrong_token":
        token = "another-token"
    elif fault == "wrong_worker":
        worker = "another-worker"
    elif fault == "stale_attempt":
        connection.row["attempts"] = 3
        connection.attempt["token_sha256"] = hashlib.sha256(b"new-attempt-token").hexdigest()
    elif fault == "queued":
        connection.row["status"] = "queued"
    else:
        connection.row.update(status="running", lease_token="another-token")
    with pytest.raises(LeaseLost):
        queue.record_cleanup_diagnostic("fixture-job", token, worker=worker,
            stage="cleanup_ack", error_type="LeaseLost")
    assert queue.recorded_events == []


@pytest.mark.parametrize("stage,error_type", [
    ("owned_cleanup", "secret exception details"), ("cleanup_ack", "RuntimeError"),
    ("unknown-stage", "ValueError"), ("owned_cleanup", "ValueError"),
    ([], "ValueError"), ("cleanup_ack", {})])
def test_cleanup_diagnostic_accepts_only_fixed_safe_categories(stage, error_type):
    queue, connection = diagnostic_queue()
    with pytest.raises(ValueError):
        queue.record_cleanup_diagnostic("fixture-job", "fixture-token", worker="fixture-worker",
            stage=stage, error_type=error_type)
    assert connection.statements == [] and queue.recorded_events == []


class BudgetQueue:
    def __init__(self, *, ack_error=None, diagnostic_error=False):
        self.job = {"id": "fixture-job", "lease_token": "fixture-token", "kind": "discovery_prepare",
                    "compute_seconds": 0, "budget_seconds": 1, "request": {}, "owner_id": "fixture-owner"}
        self.status, self.quarantined = "running", False
        self.ack_error, self.diagnostic_error = ack_error, diagnostic_error
        self.diagnostics, self.acknowledgements = [], []

    def register_node(self, *args):
        pass

    def claim(self, *args):
        return dict(self.job)

    def reserve_cpu(self, *args):
        return True

    def get(self, *args):
        return {"status": self.status}

    def heartbeat(self, *args, compute_seconds, **kwargs):
        if compute_seconds >= self.job["budget_seconds"]:
            self.status, self.quarantined = "budget_deferred", True
        return {"status": self.status}

    def fail(self, *args, **kwargs):
        if self.status != "running":
            raise LeaseLost("closed fixture lease")
        raise AssertionError("Fixture should close its budget before cleanup failure")

    def finish(self, *args, **kwargs):
        raise AssertionError("Budget fixture must never publish")

    def acknowledge_cleanup(self, job, token, *, worker, cleanup):
        assert (job, token, worker) == ("fixture-job", "fixture-token", "discovery-coordinator")
        self.acknowledgements.append(cleanup)
        if self.ack_error:
            raise self.ack_error("private fixture details must not reach diagnostic")
        self.quarantined = False
        return {"status": "cleanup_confirmed"}

    def record_cleanup_diagnostic(self, job, token, **kwargs):
        assert (job, token) == ("fixture-job", "fixture-token")
        if self.diagnostic_error:
            raise RuntimeError("fixture diagnostic storage unavailable")
        self.diagnostics.append(kwargs)


def run_engine_budget_case(queue, job_id, tmp_path, monkeypatch, *, fail_proof=False):
    """Real EngineRunner/run_owned process lifetime, synthetic engine executable only."""
    binary = tmp_path / "harmless-engine"
    binary.write_text("#!" + sys.executable + "\nimport time\ntime.sleep(30)\n")
    binary.chmod(0o700)
    store = ArtifactStore(tmp_path / "artifacts")
    atom = b"ATOM      1  C   UNL A   1       0.000   0.000   0.000  1.00  0.00     0.000 C\n"
    digest = store.put_bytes(atom)["sha256"]
    request = {"engine": "vina", "receptor_sha256": digest, "center": [0., 0., 0.],
               "size": [10., 10., 10.], "ligands": [{"candidate_id": "fixture-ligand", "pdbqt_sha256": digest}]}
    runner = engines.EngineRunner(tmp_path / "engine", {"vina": str(binary)}, cpu_threads=1)
    children, proofs = [], []
    original_popen, original_terminate = subprocess.Popen, engines.terminate_owned
    original_ack = queue.acknowledge_cleanup

    def launch(command, **kwargs):
        child = original_popen(command, **kwargs)
        if command[0] == str(binary):
            children.append(child)
        return child

    def verify_or_fail(process):
        original_terminate(process)
        if fail_proof and process in children:
            # Real cleanup happens first. Only its proof is deliberately denied.
            deadline = time.monotonic() + 3
            while queue.get(job_id)["status"] != "budget_deferred" and time.monotonic() < deadline:
                time.sleep(.01)
            assert queue.get(job_id)["status"] == "budget_deferred"
            raise engines.OwnedProcessCleanupError("synthetic cleanup proof unavailable")

    def acknowledge(*args, **kwargs):
        assert len(children) == 1 and children[0].poll() is not None
        assert engines._group_extinguished(children[0].pid)
        proofs.append(True)
        return original_ack(*args, **kwargs)

    monkeypatch.setattr(engines.subprocess, "Popen", launch)
    monkeypatch.setattr(engines, "terminate_owned", verify_or_fail)
    monkeypatch.setattr(queue, "acknowledge_cleanup", acknowledge)
    monkeypatch.setattr("skinpipe.platform.resources.probe", lambda **_: {
        "physical_host_id": "fixture-host", "cpu_count": 8, "cpu_available": 8, "cpu_limit": 2})
    runs = DiscoveryRuns.__new__(DiscoveryRuns)
    runs.queue, runs.service = queue, SimpleNamespace(root=tmp_path)
    runs.execute = lambda body, **kwargs: research_docking._execute_engine_fenced(
        runner, request, store, stop=kwargs["stop"], progress=kwargs["progress"], remaining_seconds=20)
    try:
        result = runs.run_once()
        assert len(children) == 1 and children[0].poll() is not None
        assert engines._group_extinguished(children[0].pid)
        return result, proofs
    finally:
        for child in children:
            original_terminate(child)


@pytest.mark.parametrize("fail_proof", [False, True])
def test_engine_budget_cleanup_preserves_diagnostic_and_quarantine(tmp_path, monkeypatch, fail_proof):
    queue = BudgetQueue()
    result, proofs = run_engine_budget_case(queue, "fixture-job", tmp_path, monkeypatch, fail_proof=fail_proof)
    assert queue.status == "budget_deferred"
    assert queue.quarantined is fail_proof
    assert result["status"] == ("cleanup_unverified" if fail_proof else "budget_deferred")
    assert proofs == ([] if fail_proof else [True])
    assert queue.diagnostics == ([{"worker": "discovery-coordinator", "stage": "owned_cleanup",
                                  "error_type": "OwnedProcessCleanupError"}] if fail_proof else [])


@pytest.mark.parametrize("error", [LeaseLost, ValueError])
def test_engine_cleanup_ack_rejection_records_fixed_diagnostic(tmp_path, monkeypatch, error):
    queue = BudgetQueue(ack_error=error)
    result, proofs = run_engine_budget_case(queue, "fixture-job", tmp_path, monkeypatch)
    assert result["status"] == "budget_deferred" and proofs == [True] and queue.quarantined
    assert queue.diagnostics == [{"worker": "discovery-coordinator", "stage": "cleanup_ack", "error_type": error.__name__}]


@pytest.mark.parametrize("fail_proof", [False, True])
def test_diagnostic_failure_never_replaces_original_hold(tmp_path, monkeypatch, fail_proof):
    queue = BudgetQueue(ack_error=None if fail_proof else ValueError, diagnostic_error=True)
    result, _ = run_engine_budget_case(queue, "fixture-job", tmp_path, monkeypatch, fail_proof=fail_proof)
    assert result["status"] == ("cleanup_unverified" if fail_proof else "budget_deferred")
    assert queue.quarantined and not queue.diagnostics
