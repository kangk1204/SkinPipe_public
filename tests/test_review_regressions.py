"""Synthetic regressions for the public-shell and frozen-input review."""
# ruff: noqa: E402 -- Import the checkout, never a different installed release.
from copy import deepcopy
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "skinpipe" if (ROOT / "skinpipe/skinpipe").is_dir() else ROOT
sys.path.insert(0, str(SOURCE))

from skinpipe.platform.admet_search import annotate_risk, evaluate_risk
from skinpipe.platform.artifacts import ArtifactStore
from skinpipe.platform.exports import purchase_links, review_columns, toxicity_columns
from skinpipe.platform.identity import digest
from skinpipe.platform.ranking import rank_candidates
from skinpipe.platform.research_inputs import normalize_input, resolve_input
from skinpipe.platform.routine_docking import _cache_key, evaluate_routine_complex
from skinpipe.platform.routine_screening import execution_sources, screening_policy, validate_policy
from skinpipe.structure import rcsb


class ReviewRegressionTests(unittest.TestCase):
    def test_purchase_urls_use_stored_products_and_reject_unsafe_destinations(self):
        row = {"offers": [{"product_url": "https://example.org/product", "supplier": "Synthetic"},
                          {"url": "javascript:alert(1)"}, {"url": "https://user:secret@example.org/"}],
               "supplier_offers": [{"url": "https://example.org/product"}]}
        original = deepcopy(row)
        self.assertEqual(purchase_links(row), [{"url": "https://example.org/product", "supplier": "Synthetic"}])
        self.assertEqual(row, original)
    def test_single_prediction_has_identical_filter_and_export_values(self):
        prediction = {"model_id": "admet-ai-synthetic", "status": "success", "values": {"AMES": 0.95}}
        row = {"candidate_id": "synthetic", "admet": prediction}
        original = deepcopy(row)
        self.assertEqual(evaluate_risk(prediction), evaluate_risk([prediction]))
        annotate_risk(row)
        self.assertTrue(row["admet_search_risk"]["excluded"])
        self.assertEqual(row["admet"], original["admet"])
        self.assertEqual(review_columns(row)["review_tox_AMES"], 0.95)
        self.assertEqual(toxicity_columns(row, [])["tox_AMES"][0]["value"], 0.95)

    def test_malformed_prediction_items_remain_unassessed(self):
        for value in (None, "invalid", 7, [None, "invalid", 7], {"model_manifest": "invalid"}):
            with self.subTest(value=value):
                self.assertEqual(evaluate_risk(value)["status"], "unassessed")
                self.assertEqual(review_columns({"admet": value})["review_tox_AMES"], "미평가")
        self.assertFalse(evaluate_risk({"model_id": "admet-ai-synthetic", "status": "success",
                                      "values": {"AMES": None}})["excluded"])

    def test_scores_remain_unchanged_while_missingness_and_role_are_readable(self):
        row = {"candidate_id": "synthetic", "target_id": "synthetic-target", "role": "query_seed",
               "discovery_score": 84.23, "topical_score": 9.18, "discovery_terms": {"topical": 50},
               "discovery_missing": ["topical"], "skin_goal_assessment": {"permeability_ranking": "neutral"}}
        original = deepcopy(row)
        ranked = rank_candidates([row])[0]
        self.assertEqual(row, original)
        self.assertEqual(ranked["priority_score_v2"], 60.269)
        self.assertEqual(ranked["priority_score_observed_channels"], 84.23)
        columns = review_columns(ranked)
        self.assertEqual(columns["review_priority_neutral_percent"], 70)
        self.assertEqual(columns["review_role"], "출발구조")
        self.assertIn("중립 50", columns["review_topical_contribution"])
        self.assertIn("neutral_imputation", columns["priority_score_status"])
        ranked.pop("priority_score_status")
        self.assertIn("neutral_imputation", review_columns(ranked)["priority_score_status"])
        for extra, label in (({"role": "candidate"}, "탐색 후보"),
                             ({"role": "control", "control_kind": "measured_reference"}, "측정 자료 참조"),
                             ({"redock_reference": True}, "관측 리간드 재도킹 대조")):
            self.assertEqual(review_columns({**row, **extra})["review_role"], label)

    def test_compatible_old_policy_preserves_frozen_provenance(self):
        policy = screening_policy(target_candidates=2, compute_budget_seconds=360)
        policy["execution_sources"] = {name: "a" * 64 for name in policy["execution_sources"]}
        policy["sha256"] = digest({k: v for k, v in policy.items() if k != "sha256"})
        original = deepcopy(policy)
        self.assertEqual(validate_policy(policy), original)
        self.assertEqual(policy, original)
        self.assertNotEqual(execution_sources(), policy["execution_sources"])
        rows, receipt = evaluate_routine_complex([], {"docking_hold_reasons": ["synthetic_hold"]}, [],
                                                None, {}, {"screening_policy": policy},
                                                stop=lambda: None, progress=lambda _: None)
        self.assertEqual(rows, [])
        self.assertEqual(receipt["screening_policy"], policy)
        self.assertEqual(receipt["execution_sources"], execution_sources())
        self.assertNotEqual(_cache_key({}, {}, {"smiles": "CCO"}, "b" * 64, execution_sources()),
                            _cache_key({}, {}, {"smiles": "CCO"}, "b" * 64, policy["execution_sources"]))

    def test_policy_semantic_and_integrity_changes_are_rejected(self):
        for key, value in (("max_new_docking", 50000), ("score_fraction", 1.0), ("version", "unknown")):
            policy = screening_policy()
            policy[key] = value
            policy["sha256"] = digest({k: v for k, v in policy.items() if k != "sha256"})
            with self.subTest(key=key), self.assertRaises(ValueError):
                validate_policy(policy)
        policy = screening_policy()
        policy["execution_sources"]["engines.py"] = "b" * 64
        with self.assertRaises(ValueError):
            validate_policy(policy)

    def test_1001_entries_are_batched_and_predicted_entries_stay_separate(self):
        identifiers = [f"A{i:03X}" for i in range(1001)]
        calls = []
        def provider(url, body):
            self.assertEqual(url, rcsb._GRAPHQL_URL)
            chunk = body["variables"]["ids"]
            calls.append(chunk)
            return {"data": {"entries": [{"rcsb_id": identifier, "struct": {"title": "Synthetic"}}
                                          for identifier in chunk]}}
        with patch.object(rcsb, "_post_json", side_effect=provider):
            result = rcsb.entry_metadata(identifiers + [identifiers[0].lower(), "AF_SYNTHETIC"])
        self.assertEqual([len(chunk) for chunk in calls], [500, 500, 1])
        self.assertEqual(set(result), set(identifiers) | {"AF_SYNTHETIC"})
        self.assertEqual(result["AF_SYNTHETIC"]["method"], "predicted")

    def test_incomplete_or_failed_metadata_batch_is_not_a_successful_partial_result(self):
        for payload in ({"data": {"entries": []}}, {"errors": [{"message": "synthetic failure"}]},
                        {"data": {"entries": [None]}}):
            with self.subTest(payload=payload), patch.object(rcsb, "_post_json", return_value=payload):
                with self.assertRaises(rcsb.RCSBError):
                    rcsb.entry_metadata(["1ABC"])

    def test_name_with_more_than_twenty_matches_accepts_verified_selected_cid(self):
        requests = []
        def source(url, params):
            requests.append(url)
            if "/name/" in url:
                return {"IdentifierList": {"CID": list(range(1, 26))}}
            self.assertIn("/cid/25/property/", url)
            return {"PropertyTable": {"Properties": [{"CID": 25, "Title": "Synthetic", "SMILES": "CCO"}]}}
        with tempfile.TemporaryDirectory() as temporary:
            result = resolve_input(normalize_input({"kind": "compound", "format": "name", "value": "Synthetic",
                                                   "pubchem_cid": 25}),
                                   artifact_store=ArtifactStore(Path(temporary)), protein_resolver=None, source_get=source)
        self.assertEqual(result["status"], "verified")
        self.assertEqual(len(requests), 2)
