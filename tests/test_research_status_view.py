"""Progress reads preserve state while leaving bulky evidence in full results."""
from copy import deepcopy
import json
from types import SimpleNamespace
import unittest

from skinpipe.platform.discovery import DiscoveryRuns
from skinpipe.platform.service import PlatformService


class ResearchStatusViewTests(unittest.TestCase):
    def test_projection_preserves_lifecycle_and_counters_without_mutation(self):
        run = {"run_id": "fixture", "status": "running", "counts": {"docking_completed": 2},
               "research_input": {"kind": "bound_complex", "format": "pdb_id", "value": "1CBS"},
               "progress": {"stage": "routine_docking", "completed": 4, "total": 32,
                            "qualified": 3, "target": 200, "raw_rows": ["x" * 100000]},
               "seed_assessment": {"raw_model_output": "x" * 1000000}}
        before = deepcopy(run)
        view = DiscoveryRuns._status_view(run)
        self.assertEqual(view["progress"], {k: v for k, v in run["progress"].items() if k != "raw_rows"})
        self.assertEqual(view["counts"], run["counts"])
        self.assertNotIn("seed_assessment", view)
        self.assertLess(len(json.dumps(view)), 1000)
        view["counts"]["docking_completed"] = 99
        self.assertEqual(run, before)

    def test_cleanup_and_cancellation_proof_are_not_synthesized(self):
        for active in (True, False):
            run = {"status": "needs_data", "lifecycle_state": "active",
                   "cancellation": {"requested": False, "resources_released": None},
                   "progress": {"resource_quarantine": {"active": active},
                                "resource_cleanup": {"confirmed": False}}}
            self.assertEqual(DiscoveryRuns._status_view(run), run)

    def test_status_uses_existing_owner_checked_get(self):
        calls = []
        runs = DiscoveryRuns.__new__(DiscoveryRuns)
        runs.get = lambda run_id, *, owner: calls.append((run_id, owner)) or {"status": "queued"}
        self.assertEqual(runs.status("fixture", owner="alice"), {"status": "queued", "progress": None})
        self.assertEqual(calls, [("fixture", "alice")])

    def test_opt_in_routes_leave_original_routes_unchanged(self):
        calls = []
        def call(name, *args, **kwargs):
            calls.append((name, args, kwargs))
            return name
        service = PlatformService.__new__(PlatformService)
        service.discovery = SimpleNamespace(**{name: lambda *a, _name=name, **kw: call(_name, *a, **kw)
                                              for name in ("get", "list", "status", "list_status")})
        for path, expected in [("/v2/discovery-runs", "list"),
                               ("/v2/discovery-runs?view=status", "list_status"),
                               ("/v2/discovery-runs/fixture", "get"),
                               ("/v2/discovery-runs/fixture?view=status", "status")]:
            self.assertEqual(service._get(path, owner="alice"), expected)
            self.assertEqual(calls[-1][2], {"owner": "alice"})
