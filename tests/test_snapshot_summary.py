"""Local regression for bounded snapshot-list metadata and its trust boundaries."""
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from skinpipe.platform.catalog import Catalog
from skinpipe.platform.service import PlatformService


class SnapshotSummaryTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.service = PlatformService.__new__(PlatformService)
        self.service.catalog = Catalog(Path(self.directory.name) / "catalog")
        self.provenance = {
            "owner_id": "alice",
            "marker": "Synthetic study",
            "target_id": "T",
            "screen_id": "screen-fixture",
            "search_mode": "similarity",
            "primary_docking_protocols": [
                {"target_id": "T", "protocol_id": "a" * 64, "engine": "vina",
                 "seed": 42, "candidate_count": 1}
            ],
            "references": {"raw_source": "fixture " * 200_000},
            "catalog_receipt": {"private_fixture_note": "not needed by list UI"},
        }
        self.service.catalog.create_snapshot(
            [{"candidate_id": "candidate", "target_id": "T", "smiles": "CCO"}],
            "owned", provenance=deepcopy(self.provenance),
        )
        self.service.catalog.create_snapshot(
            [], "empty", provenance={"owner_id": "alice", "target_id": "empty-target"},
        )
        self.service.catalog.create_snapshot(
            [], "other", provenance={"owner_id": "bob", "marker": "Other project"},
        )

    def get(self, suffix=""):
        return self.service.get("/v2/snapshots" + suffix, owner="alice")

    def test_default_and_explicit_full_preserve_complete_metadata(self):
        original = {
            "snapshots": [
                item for item in self.service.catalog.list_snapshots()
                if self.service.catalog.snapshot_owner(item["snapshot_id"]) == "alice"
            ]
        }
        self.assertEqual(self.get(), original)
        self.assertEqual(self.get("?view=full"), original)
        # Previously ignored unrelated query parameters keep their old behavior.
        self.assertEqual(self.get("?unused=legacy"), original)
        self.assertEqual(original["snapshots"][1]["provenance"], self.provenance)

    def test_summary_retains_exact_counts_selection_and_protocol_metadata(self):
        full = self.get()
        compact = self.get("?view=summary")
        self.assertEqual({s["snapshot_id"] for s in compact["snapshots"]}, {"owned", "empty"})
        for before, after in zip(full["snapshots"], compact["snapshots"]):
            for key in ("snapshot_id", "schema_version", "created_at", "candidate_count",
                        "target_counts", "status"):
                self.assertEqual(after[key], before[key])
            self.assertNotIn("files", after)
            self.assertNotIn("references", after["provenance"])
            self.assertNotIn("catalog_receipt", after["provenance"])
            self.assertNotIn("owner_id", after["provenance"])
        owned = next(s for s in compact["snapshots"] if s["snapshot_id"] == "owned")
        for key in ("marker", "target_id", "screen_id", "search_mode", "primary_docking_protocols"):
            self.assertEqual(owned["provenance"][key], self.provenance[key])
        self.assertEqual(compact["snapshots"][0]["candidate_count"], 0)
        self.assertEqual(compact["snapshots"][0]["target_counts"], {})
        self.assertLess(len(json.dumps(compact)), 2000)
        self.assertGreater(len(json.dumps(full)), 1_000_000)

    def test_summary_does_not_mutate_or_alias_catalog_metadata(self):
        path = self.service.catalog.root / "owned/manifest.json"
        source_sha = hashlib.sha256(path.read_bytes()).hexdigest()
        original = self.get()
        # Also exercise a catalog adapter that returns the same object again.
        shared = self.service.catalog.list_snapshots()
        pristine = deepcopy(shared)
        with patch.object(self.service.catalog, "list_snapshots", return_value=shared):
            summary = self.get("?view=summary")
            item = next(s for s in summary["snapshots"] if s["snapshot_id"] == "owned")
            item["target_counts"]["T"] = 999
            item["provenance"]["primary_docking_protocols"][0]["seed"] = 0
            self.assertEqual(shared, pristine)
        self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), source_sha)
        self.assertEqual(self.get(), original)

    def test_tampering_remains_quarantined_and_other_owner_is_not_disclosed(self):
        self.get("?view=summary")  # Populate the normal verification cache.
        for snapshot in ("owned", "other"):
            with (self.service.catalog.root / snapshot / "candidates.parquet").open("ab") as stream:
                stream.write(b"synthetic tamper")
        compact = self.get("?view=summary")
        self.assertEqual({s["snapshot_id"] for s in compact["snapshots"]}, {"owned", "empty"})
        failed = next(s for s in compact["snapshots"] if s["snapshot_id"] == "owned")
        self.assertEqual(failed["status"], "quarantined")
        self.assertEqual(failed["ownership_status"], "verified")
        self.assertIn("checksum mismatch", failed["reason"])
        self.assertNotIn("candidate_count", failed)
        full_failed = next(s for s in self.get()["snapshots"] if s["snapshot_id"] == "owned")
        self.assertEqual(failed["reason"], full_failed["reason"])

    def test_summary_preserves_explicit_owner_authorization(self):
        with self.assertRaises(PermissionError):
            self.get("?view=summary&snapshot_id=other")
        self.assertEqual(self.service.get("/v2/snapshots?view=summary", owner="unknown"),
                         {"snapshots": []})

    def test_invalid_or_ambiguous_view_rejected_before_listing(self):
        for query in ("?view=", "?view", "?view=unknown", "?view=Summary",
                      "?view=summary&view=full", "?view=summary&view=summary"):
            with self.subTest(query=query):
                with patch.object(self.service.catalog, "list_snapshots") as listing:
                    with self.assertRaisesRegex(ValueError, "full or summary exactly once"):
                        self.get(query)
                    listing.assert_not_called()
        self.assertEqual(self.get("?view=%73ummary"), self.get("?view=summary"))

    def test_catalog_listing_failure_is_not_reclassified_as_an_empty_success(self):
        with patch.object(self.service.catalog, "list_snapshots",
                          side_effect=OSError("synthetic storage failure")):
            with self.assertRaisesRegex(OSError, "synthetic storage failure"):
                self.get("?view=summary")


if __name__ == "__main__":
    unittest.main()
