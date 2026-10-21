"""Harmless local regressions for runtime ownership and input boundaries."""
import os
from pathlib import Path
import tempfile
import unittest
import uuid

from skinpipe.platform.queue import Queue, SCHEMA
from skinpipe.platform.worker import Worker, WorkerHTTPError


class WorkerShutdownTests(unittest.TestCase):
    def test_stopped_worker_does_not_poll_or_claim(self):
        with tempfile.TemporaryDirectory() as directory:
            worker = Worker({"workdir": directory}, client=object())
            worker.shutdown()
            self.assertFalse(worker.run_once())

    def test_shutdown_during_claim_does_not_start_native_work(self):
        calls = []

        class Client:
            def call(self, path, payload=None, **kwargs):
                calls.append(path)
                if path.endswith("/register"):
                    return {}
                if path.endswith("/claim"):
                    worker.shutdown()
                    return {"job": {"id": "synthetic-job", "lease_token": "synthetic-lease",
                                    "kind": "dock_batch", "request": {"engine": "vina"},
                                    "compute_seconds": 0}}
                if path.endswith("/reserve-cpu"):
                    return {"reserved": True}
                if path.endswith("/heartbeat"):
                    return {"status": "running"}
                if path.endswith("/cancel-ack"):
                    # A local shutdown must not impersonate user cancellation.
                    raise WorkerHTTPError(400, "Cancellation was not requested")
                return {}

            def fetch_artifact(self, digest):
                raise AssertionError("Stopped worker must not fetch inputs")

            def publish_artifact(self, data):
                raise AssertionError("Stopped worker must not publish results")

        class Runner:
            def execute(self, *args, **kwargs):
                calls.append("native-execution")
                return {"compute_seconds": 0}

        with tempfile.TemporaryDirectory() as directory:
            worker = Worker({"workdir": directory}, client=Client(), engine_runner=Runner())
            worker.telemetry = lambda: {"capabilities": ["probe"], "gpus": [], "cpu_limit": 1}
            self.assertTrue(worker.run_once())
            self.assertIsNone(worker.active_cancel)
        self.assertNotIn("native-execution", calls)
        self.assertNotIn("/v2/worker/reserve-cpu", calls)
        self.assertNotIn("/v2/worker/finish", calls)
        self.assertIn("/v2/worker/cleanup-ack", calls)


@unittest.skipUnless(os.environ.get("SKINPIPE_RUNTIME_AUDIT_DSN"),
                     "Explicit isolated PostgreSQL audit DSN required")
class QueueBudgetOwnershipTests(unittest.TestCase):
    def setUp(self):
        self.queue = Queue(os.environ["SKINPIPE_RUNTIME_AUDIT_DSN"], initialize=False)
        with self.queue.connection() as conn:
            database = conn.execute("SELECT current_database() AS name").fetchone()["name"]
            if not database.endswith("_test"):
                self.fail("Use a dedicated database ending in _test")
            conn.execute(SCHEMA)
        self.namespace = "audit-" + uuid.uuid4().hex

    def submit(self, owner, scope, key, budget=100):
        return self.queue.enqueue("probe", {}, owner=owner, project=owner,
                                  idempotency_key=key, budget_scope=scope,
                                  budget_seconds=budget)

    def test_colliding_scope_cannot_share_another_owners_budget(self):
        original = self.submit(self.namespace + ":a", "b", "first")
        with self.assertRaisesRegex(ValueError, "another owner"):
            self.submit(self.namespace, "a:b", "second", 200)
        self.assertEqual(self.queue.list(owner=self.namespace), [])
        replay = self.submit(self.namespace + ":a", "b", "first")
        self.assertEqual(replay["id"], original["id"])
        sibling = self.submit(self.namespace + ":a", "b", "sibling")
        self.assertEqual(sibling["budget_scope"], original["budget_scope"])

    def test_historical_foreign_budget_cannot_resume_or_claim(self):
        original = self.submit(self.namespace + ":a", "b", "first")
        foreign = self.submit(self.namespace, "independent", "second", 200)
        with self.queue.connection() as conn:
            conn.execute("UPDATE platform_jobs SET budget_scope=%s WHERE id=%s",
                         (original["budget_scope"], foreign["id"]))
        with self.assertRaisesRegex(ValueError, "another owner"):
            self.submit(self.namespace, "independent", "second", 200)
        # Claim only this synthetic kind so unrelated test records are ignored.
        kind = "probe-" + uuid.uuid4().hex
        with self.queue.connection() as conn:
            conn.execute("UPDATE platform_jobs SET kind=%s WHERE id=%s", (kind, foreign["id"]))
        self.assertIsNone(self.queue.claim(self.namespace, [kind]))
        self.queue.cancel(foreign["id"], owner=self.namespace)
        with self.assertRaisesRegex(ValueError, "another owner"):
            self.queue.resume(foreign["id"], owner=self.namespace, budget_seconds=300)
        with self.queue.connection() as conn:
            bound = conn.execute("SELECT limit_seconds,used_seconds FROM platform_compute_budgets WHERE id=%s",
                                 (original["budget_scope"],)).fetchone()
        self.assertEqual(bound, {"limit_seconds": 100, "used_seconds": 0})


class ObservedStructureStereoTests(unittest.TestCase):
    def test_specified_stereo_requires_coordinate_support(self):
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from skinpipe.platform.research_structures import _map_template

        molecule = Chem.AddHs(Chem.MolFromSmiles("C[C@H](O)F"))
        self.assertEqual(AllChem.EmbedMolecule(molecule, randomSeed=42), 0)
        molecule = Chem.RemoveHs(molecule)
        observed = [
            {"name": "A" + str(index), "serial": str(index + 1),
             "element": atom.GetSymbol().upper(),
             "xyz": list(molecule.GetConformer().GetAtomPosition(index))}
            for index, atom in enumerate(molecule.GetAtoms())
        ]
        mapping = list(range(len(observed)))
        intact = _map_template(molecule, observed, names=mapping)
        self.assertEqual(Chem.MolToSmiles(intact), "C[C@H](O)F")
        planar = [{**atom, "xyz": [*atom["xyz"][:2], 0.0]} for atom in observed]
        with self.assertRaisesRegex(ValueError, "stereochemistry"):
            _map_template(molecule, planar, names=mapping)
        unspecified = Chem.MolFromSmiles("CC(O)F")
        self.assertIsNotNone(_map_template(unspecified, planar, names=mapping))


class CampaignControlTests(unittest.TestCase):
    def test_reference_control_and_its_parent_stay_out_of_candidate_selection(self):
        from skinpipe.platform.campaign import CampaignCoordinator

        ranked = [
            {"candidate_id": "reference", "smiles": "CCO", "role": "reference_control"},
            {"candidate_id": "same-parent", "smiles": "CCO", "role": "candidate"},
            {"candidate_id": "independent", "smiles": "CCN", "role": "candidate"},
        ]
        selected = CampaignCoordinator._unique_ranked_parents(ranked)
        self.assertEqual([row["candidate_id"] for row in selected], ["independent"])


class CampaignInputOwnershipTests(unittest.TestCase):
    def test_manifest_indirection_does_not_grant_foreign_snapshot_access(self):
        from skinpipe.platform.artifacts import ArtifactStore, canonical_json
        from skinpipe.platform.campaign_api import CampaignAPI, INPUT_SCHEMA
        from skinpipe.platform.catalog import Catalog

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            catalog = Catalog(root / "catalog")
            catalog.create_snapshot(
                [{"candidate_id": "fixture-candidate", "target_id": "fixture-target", "smiles": "CCO"}],
                snapshot_id="fixture-snapshot", provenance={"owner_id": "fixture-owner"})
            artifacts = ArtifactStore(root / "artifacts")
            api = CampaignAPI(root / "campaigns", catalog, artifacts, None)
            manifest = {"schema": INPUT_SCHEMA, "snapshot_id": "fixture-snapshot",
                        "target_id": "fixture-target", "prepared": {}}
            receipt = artifacts.put_bytes(canonical_json(manifest))
            body = {"manifest_sha256": receipt["sha256"]}
            with self.assertRaises(PermissionError):
                api.register_inputs(body, owner="different-owner")
            plan = {"snapshot_id": manifest["snapshot_id"], "target_id": manifest["target_id"]}
            with self.assertRaises(PermissionError):
                api.plan(plan, owner="different-owner")
            registered = api.register_inputs(body, owner="fixture-owner")
            loaded = api.get("/v2/campaign-inputs/" + registered["input_id"], owner="fixture-owner")
            self.assertEqual(loaded["snapshot_id"], manifest["snapshot_id"])
            # A historical foreign registration still cannot resolve its inputs.
            from skinpipe.platform.campaign import _sealed_write
            from skinpipe.platform.identity import digest
            record = {"owner": "different-owner", "manifest_sha256": receipt["sha256"]}
            identity = digest(record)
            _sealed_write(api._owner_root("different-owner") / "inputs" / (identity + ".json"), record)
            with self.assertRaises(PermissionError):
                api.get("/v2/campaign-inputs/" + identity, owner="different-owner")


@unittest.skipUnless(os.environ.get("SKINPIPE_RUNTIME_AUDIT_DSN"),
                     "Explicit isolated PostgreSQL audit DSN required")
class ConcurrentExportTests(unittest.TestCase):
    def test_export_completes_while_analysis_remains_owned(self):
        import threading
        from types import SimpleNamespace
        from unittest.mock import patch
        import psycopg
        from psycopg import sql
        from psycopg.conninfo import make_conninfo
        from skinpipe.platform.discovery import DiscoveryRuns
        dsn = os.environ["SKINPIPE_RUNTIME_AUDIT_DSN"]
        namespace = "export_" + uuid.uuid4().hex + "_test"
        with psycopg.connect(dsn, autocommit=True) as conn:
            self.assertTrue(conn.execute("SELECT current_database()").fetchone()[0].endswith("_test"))
            conn.execute(sql.SQL("CREATE DATABASE {}").format(sql.Identifier(namespace)))
        entered, release = threading.Event(), threading.Event()
        thread = None
        try:
            queue = Queue(make_conninfo(dsn, dbname=namespace))
            first = queue.enqueue("discovery_prepare", {}, owner="synthetic", project="synthetic",
                                  idempotency_key="first", budget_scope="first", budget_seconds=60)
            second = queue.enqueue("discovery_prepare", {}, owner="synthetic", project="synthetic",
                                   idempotency_key="second", budget_scope="second", budget_seconds=60)
            export = queue.enqueue("result_export", {}, owner="synthetic", project="synthetic",
                                   idempotency_key="export", budget_scope="export", budget_seconds=60)
            with tempfile.TemporaryDirectory() as directory:
                runs = DiscoveryRuns.__new__(DiscoveryRuns)
                runs.queue = queue
                def analysis(*args, **kwargs):
                    entered.set()
                    if not release.wait(10):
                        raise TimeoutError("Synthetic analysis release timed out")
                    return {"fixture": "analysis"}
                runs.execute = analysis
                runs.service = SimpleNamespace(root=Path(directory), exports=SimpleNamespace(
                    execute=lambda *args, **kwargs: {"fixture": "export"}))
                telemetry = {"physical_host_id": namespace, "cpu_count": 16, "cpu_available": 16, "cpu_limit": 4}
                with patch("skinpipe.platform.resources.probe", return_value=telemetry):
                    thread = threading.Thread(target=runs.run_once)
                    thread.start()
                    self.assertTrue(entered.wait(5))
                    result = runs.run_once(export_only=True)
                    self.assertIsNotNone(result, queue.get(export["id"]))
                    self.assertEqual(queue.get(export["id"])["status"], "complete")
                    self.assertEqual(queue.get(first["id"])["status"], "running")
                    self.assertEqual(queue.get(second["id"])["status"], "queued")
                    self.assertIsNone(runs.run_once(export_only=True))
                    release.set()
                    thread.join(5)
                    self.assertFalse(thread.is_alive())
                    self.assertEqual(queue.get(first["id"])["status"], "complete")
        finally:
            release.set()
            if thread is not None:
                thread.join(10)
                self.assertFalse(thread.is_alive())
            with psycopg.connect(dsn, autocommit=True) as conn:
                conn.execute(sql.SQL("DROP DATABASE {}").format(sql.Identifier(namespace)))


if __name__ == "__main__":
    unittest.main()
