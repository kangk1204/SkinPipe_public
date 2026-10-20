"""Synthetic data-boundary regressions; no external catalog or model execution."""
# ruff: noqa: E402 -- Resolve the checkout before importing its package.
import json
import os
from pathlib import Path
import sys
import tempfile
import threading
from types import SimpleNamespace
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "skinpipe" if (ROOT / "skinpipe/skinpipe").is_dir() else ROOT
sys.path.insert(0, str(SOURCE))

import duckdb
from openpyxl import load_workbook

from skinpipe.platform.admet_search import ADMETSearchGate, evaluate_risk
from skinpipe.platform.artifacts import ArtifactStore, canonical_json, sha256_file
from skinpipe.platform.catalog import Catalog
from skinpipe.platform.catalog_hopping import select_docking
from skinpipe.platform.catalog_query_lifetime import CatalogQueryCancelled, query_scope
from skinpipe.platform.database_bundle import _sql_paths, install_database, package_database
from skinpipe.platform.exports import RunExports, review_columns, toxicity_columns
from skinpipe.platform.identity import chemical_identity
from skinpipe.platform.scout import ScoutCatalog
from skinpipe.platform.scout_predictions import PublishedPredictionCache


class DataBoundaryTests(unittest.TestCase):
    def test_implicit_file_views_cannot_bypass_dependency_sealing(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            (root / "data").mkdir()
            parquet = root / "data/implicit.parquet"
            with duckdb.connect(str(root / "data/skinscout.duckdb")) as connection:
                connection.execute("COPY (SELECT 'synthetic' AS compound_id) TO ? (FORMAT PARQUET)",
                                   [str(parquet)])
                connection.execute(
                    "CREATE VIEW compounds AS SELECT * FROM '"
                    + str(parquet).replace("'", "''") + "'"
                )
                sql = connection.execute(
                    "SELECT sql FROM duckdb_views() WHERE view_name='compounds'"
                ).fetchone()[0]
            before = sha256_file(parquet)
            with self.assertRaisesRegex(ValueError, "explicit.*Parquet"):
                _sql_paths(sql)
            with self.assertRaisesRegex(ValueError, "explicit.*Parquet"):
                ScoutCatalog(root).database_identity()
            self.assertEqual(sha256_file(parquet), before)

    def test_path_like_identifiers_are_held_but_escaped_literal_data_is_allowed(self):
        for sql in (
            'SELECT * FROM "relative.parquet"',
            'SELECT * FROM "relative/data.parquet"',
            'SELECT * FROM ordinary, "second.parquet"',
            'SELECT * FROM ordinary JOIN "odd""name.parquet" USING(id)',
        ):
            with self.subTest(sql=sql), self.assertRaisesRegex(ValueError, "explicit.*Parquet"):
                _sql_paths(sql)
        self.assertEqual(_sql_paths("""SELECT 'a "named.parquet" value' AS label FROM ordinary"""), [])

    def test_cancel_interrupts_serving_gate_and_prediction_sql(self):
        for component in (ADMETSearchGate, PublishedPredictionCache):
            with self.subTest(component=component.__name__):
                event, emergency = threading.Event(), threading.Event()
                native = duckdb.connect(":memory:")
                timers = []

                class Connection:
                    def __enter__(self):
                        return self

                    def __exit__(self, *_):
                        native.close()

                    def interrupt(self):
                        native.interrupt()

                    def execute(self, sql, *args):
                        if sql.startswith("SELECT"):
                            def stop_safety():
                                emergency.set()
                                native.interrupt()
                            timers.extend([threading.Timer(.05, event.set),
                                           threading.Timer(2, stop_safety)])
                            for timer in timers:
                                timer.start()
                            return native.execute(
                                "SELECT sum(sin(i+j)) FROM range(1000000) a(i)"
                                " CROSS JOIN range(1000000) b(j)"
                            )
                        return native.execute(sql, *args)

                reader = component.__new__(component)
                reader.database = Path(":memory:")
                reader.verify = lambda: ("synthetic-unchanged",)
                try:
                    with patch.object(duckdb, "connect", return_value=Connection()):
                        with self.assertRaises(CatalogQueryCancelled), query_scope(event):
                            if component is ADMETSearchGate:
                                reader.query()
                            else:
                                reader.query(["synthetic"])
                    self.assertFalse(emergency.is_set(), "Owned cancellation must interrupt before the safety timer")
                finally:
                    for timer in timers:
                        timer.cancel()
                        timer.join()
                    native.close()

    def test_structure_seal_parses_escaped_quotes_in_parquet_paths(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve() / "quoted'catalog"
            (root / "data").mkdir(parents=True)
            parquet = root / "data/compounds.parquet"
            with duckdb.connect(str(root / "data/skinscout.duckdb")) as connection:
                connection.execute("COPY (SELECT 'synthetic' AS compound_id) TO ? (FORMAT PARQUET)",
                                   [str(parquet)])
                connection.execute(
                    "CREATE VIEW compounds AS SELECT * FROM read_parquet('"
                    + str(parquet).replace("'", "''") + "')"
                )
            identity = ScoutCatalog(root).database_identity()
            self.assertEqual(identity["structure_dependencies"], [{
                "path": "data/compounds.parquet", "sha256": sha256_file(parquet),
                "size_bytes": parquet.stat().st_size,
            }])

    def test_reference_control_parent_cannot_reenter_candidate_docking(self):
        rows = [{
            "candidate_id": name, "role": role, "identity": chemical_identity(smiles),
            "discovery_filter_pass": True,
        } for name, role, smiles in (
            ("reference", "reference_control", "CCO"),
            ("duplicate", "candidate", "OCC"),
            ("distinct", "candidate", "CCN"),
        )]
        self.assertEqual(select_docking(rows), ["distinct"])

    def test_invalid_prediction_quantity_has_quarantinable_error(self):
        from skinpipe._catalog_source.prediction_integrity import validate_values
        with self.assertRaisesRegex(ValueError, "Invalid prediction quantity"):
            validate_values({"synthetic_endpoint": 10**1000}, "success", {"synthetic_endpoint"})

    def test_relative_database_install_survives_working_directory_change(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            donor = root / "donor"
            (donor / "data").mkdir(parents=True)
            (donor / "reports").mkdir()
            source = donor / "data/skinscout.duckdb"
            with duckdb.connect(str(source)) as connection:
                connection.execute(
                    "CREATE TABLE compounds(compound_id VARCHAR, identity_smiles VARCHAR)"
                )
                connection.execute("INSERT INTO compounds VALUES ('synthetic', 'CCO')")
                connection.execute(
                    "COPY compounds TO ? (FORMAT PARQUET)",
                    [str(donor / "data/compounds.parquet")],
                )
                connection.execute("DROP TABLE compounds")
                connection.execute(
                    "CREATE VIEW compounds AS SELECT * FROM read_parquet('"
                    + str(donor / "data/compounds.parquet").replace("'", "''") + "')"
                )
            checks = dict.fromkeys(
                ("missing_or_stale_build_report", "incomplete_admitted_structure_sources",
                 "duplicate_compound_ids", "orphan_structure_links", "duplicate_source_records"),
                0,
            )
            (donor / "reports/audit.json").write_bytes(
                canonical_json({"checks": checks, "compound_count": 1})
            )
            (donor / "reports/db_build.json").write_bytes(canonical_json({
                "admitted_structure_sources_complete": True, "compounds": 1,
            }))
            source_hash = sha256_file(source)
            archive = root / "catalog.tar.gz"
            package_database({"scout_root": str(donor)}, archive)
            previous = Path.cwd()
            try:
                os.chdir(root)
                receipt = install_database(archive, Path("installed"))
            finally:
                os.chdir(previous)
            installed = Path(receipt["config"]["scout_root"])
            self.assertTrue(installed.is_absolute())
            with duckdb.connect(str(installed / "data/skinscout.duckdb"), read_only=True) as connection:
                self.assertEqual(
                    connection.execute("SELECT * FROM compounds").fetchall(),
                    [("synthetic", "CCO")],
                )
            self.assertEqual(sha256_file(source), source_hash)

    def test_oversized_probability_is_missing_without_losing_valid_risk(self):
        risk = evaluate_risk([{
            "model_id": "admet-ai-synthetic", "status": "success",
            "values": {"AMES": 10**1000, "Skin_Reaction": 0.95},
        }])
        self.assertTrue(risk["excluded"])
        self.assertIn("AMES", risk["assessments"][0]["missing_endpoints"])
        self.assertEqual(risk["assessments"][0]["strong_primary_endpoints"], ["Skin_Reaction"])

    def test_oversized_prediction_is_quarantined_for_display(self):
        row = {"candidate_id": "synthetic", "admet": [{
            "model_id": "synthetic-model", "status": "success",
            "values": {"AMES": 10**1000},
        }]}
        warnings = []
        self.assertIsNone(toxicity_columns(row, warnings)["tox_AMES"])
        self.assertEqual(warnings[0]["reason"], "invalid_toxicity_display_values")
        self.assertEqual(review_columns(row)["review_tox_AMES"], "값 오류 · 원자료 확인")
        self.assertEqual(row["admet"][0]["values"]["AMES"], 10**1000)

    def test_workbook_preserves_invalid_metadata_and_numeric_raw_values(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            row = {
                "candidate_id": "synthetic", "target_id": "fixture", "smiles": "CCO",
                "mw": 10**1000, "admet": [{
                    "model_id": "synthetic\x01model", "status": "success",
                    "values": {"AMES": 10**1000},
                    "endpoint_definitions": [{
                        "id": "AMES", "name": "Synthetic\x01metadata", "units": "probability",
                        "task_type": "classification",
                    }],
                }],
            }
            catalog = Catalog(root / "catalog")
            snapshot = catalog.create_snapshot([row])
            exporter = RunExports(SimpleNamespace(
                root=root, catalog=catalog, artifacts=ArtifactStore(root / "artifacts"),
            ))
            output = root / "export"
            manifest = exporter.build({
                "snapshot_id": snapshot["snapshot_id"], "target_id": "fixture",
                "format": "xlsx", "candidate_details": False,
            }, output)
            saved = json.loads((output / manifest["chunks"][0]["file"]).read_bytes())[0]
            self.assertEqual(saved["admet"], row["admet"])
            self.assertEqual(saved["mw"], row["mw"])
            book = load_workbook(output / manifest["workbooks"][0])
            try:
                self.assertIn("원자료", book["ADMET 지표"]["C2"].value)
                self.assertTrue(any("원자료" in str(cell.value) for cell in book["후보 검토"][2]))
            finally:
                book.close()


if __name__ == "__main__":
    unittest.main()
