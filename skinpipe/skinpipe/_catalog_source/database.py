import json
import os
from pathlib import Path
import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
from .common import atomic_json, lock, now, sha256
from .prediction_schema import ensure_empty_relations


def quote(value):
    return "'" + str(value).replace("'", "''") + "'"


def connect(root, read_only=False):
    con = duckdb.connect(str(root / "data/skinscout.duckdb"), read_only=read_only)
    configure_connection(con)
    return con


def configure_connection(con):
    # Two full-size ART indexes do not fit in 2GB. Keep the independently owned
    # build within half of the documented 16GB minimum RAM, with disk spill for
    # sortable/hash-join intermediates. No original database settings are changed.
    con.execute(
        "SET threads=2; SET memory_limit='8GB'; SET preserve_insertion_order=false"
    )


def build(root):
    """Publish a database snapshot atomically; active readers retain their snapshot."""
    with lock(root / "logs/build.lock"):
        from .publish import publish

        publish(root)
        ensure_empty_relations(root)
        specs = json.loads((root / "config/sources.json").read_text())
        selected = {s["id"]: s for s in specs if s.get("admission") == "allowed"}
        paths = []
        ingests = []
        excluded = []
        from .ingest import expected_runhash, SCHEMA

        for p in sorted((root / "manifests/ingest").glob("*.json")):
            state = json.loads(p.read_text())
            sid = state.get("source_id", p.stem)
            reason = None
            spec = selected.get(sid)
            download_path = root / "manifests/downloads" / f"{sid}.json"
            download = (
                json.loads(download_path.read_text()) if download_path.exists() else {}
            )
            if spec is None:
                reason = "source_not_currently_admitted"
            elif state.get("status") != "complete":
                reason = "ingestion_not_complete"
            elif download.get("status") != "verified":
                reason = "selected_download_not_verified"
            elif spec.get("sha256") and spec["sha256"] != download.get("sha256"):
                reason = "selected_snapshot_mismatch"
            elif state.get("input_sha256") != download.get("sha256"):
                reason = "ingestion_snapshot_mismatch"
            elif state.get("runhash") != expected_runhash(spec, download["sha256"]):
                reason = "ingestion_transform_mismatch"
            elif state.get("source_version", spec.get("version")) != spec.get(
                "version"
            ):
                reason = "ingestion_version_mismatch"
            if reason:
                excluded.append({"source_id": sid, "reason": reason})
                continue
            for name in state["shards"]:
                path = root / name
                receipt = path.with_suffix(".json")
                if (
                    not path.exists()
                    or not receipt.exists()
                    or sha256(path) != json.loads(receipt.read_text()).get("sha256")
                ):
                    raise ValueError("Corrupted admitted ingestion shard: " + name)
                paths.append(str(path))
            ingests.append(state)
        if not paths:
            empty = root / "data/schema/records_empty.parquet"
            empty.parent.mkdir(parents=True, exist_ok=True)
            pq.write_table(pa.Table.from_pylist([], schema=SCHEMA), empty)
            paths = [str(empty)]
        target = root / "data/skinscout.duckdb"
        tmp = root / "data/skinscout.building.duckdb"
        if tmp.exists():
            tmp.unlink()
        con = duckdb.connect(str(tmp))
        configure_connection(con)
        plist = "[" + ",".join(map(quote, paths)) + "]"
        con.execute(
            f"CREATE VIEW source_records AS SELECT * FROM read_parquet({plist},union_by_name=true)"
        )
        con.execute("""CREATE TABLE compounds AS SELECT compound_id,canonical_smiles,identity_smiles,inchikey,
            parent_smiles,parent_id,formula,mw,fragments,formal_charge,heavy_atoms,structure_status,inference_status,normalizer
            FROM source_records WHERE compound_id IS NOT NULL
            QUALIFY row_number() OVER(PARTITION BY compound_id ORDER BY source_id,ordinal)=1""")
        con.execute("CREATE UNIQUE INDEX compound_idx ON compounds(compound_id)")
        con.execute("CREATE INDEX inchikey_idx ON compounds(inchikey)")
        registry = pa.Table.from_pylist(
            [
                dict(
                    source_id=s["id"],
                    version=s.get("version"),
                    license=s.get("license"),
                    admission=s["admission"],
                    url=s.get("url"),
                    kind=s.get("kind"),
                    reason=s.get("reason"),
                    purchasing=s.get("purchasing"),
                    purchasing_date=s.get("purchasing_date"),
                    supplier=s.get("supplier"),
                    metadata=json.dumps(s),
                )
                for s in specs
            ],
            schema=pa.schema(
                [
                    (k, pa.string())
                    for k in [
                        "source_id",
                        "version",
                        "license",
                        "admission",
                        "url",
                        "kind",
                        "reason",
                        "purchasing",
                        "purchasing_date",
                        "supplier",
                        "metadata",
                    ]
                ]
            ),
        )
        con.register("registry_input", registry)
        con.execute("CREATE TABLE sources AS SELECT * FROM registry_input")
        con.execute("""CREATE VIEW offers AS SELECT r.record_id,r.compound_id,r.source_id,r.source_record_id AS catalog_id,
            s.supplier,s.purchasing AS availability,s.purchasing_date AS evidence_date,
            CASE WHEN r.source_id='mcule_instock' THEN 'https://mcule.com/'||r.source_record_id||'/' ELSE s.url END AS purchase_route
            FROM source_records r JOIN sources s USING(source_id) WHERE s.purchasing IS NOT NULL""")
        con.execute("""CREATE VIEW natural_products AS SELECT record_id,compound_id,source_id,source_record_id,raw_metadata,rights_status
            FROM source_records WHERE source_id IN ('coconut','lotus_metadata','npatlas_2024_03')""")
        con.execute(
            "CREATE VIEW ingredients_unresolved AS SELECT * FROM source_records WHERE compound_id IS NULL"
        )
        con.execute(
            "CREATE VIEW purchasable AS SELECT * FROM compounds SEMI JOIN (SELECT compound_id FROM offers WHERE availability='catalog_reported_in_stock') USING(compound_id)"
        )
        con.execute(
            "CREATE VIEW purchasable_natural AS SELECT * FROM purchasable SEMI JOIN natural_products USING(compound_id)"
        )
        con.execute("""CREATE VIEW source_counts AS SELECT source_id,count(*) records,count(DISTINCT compound_id) unique_structures,
            count(*) FILTER(WHERE compound_id IS NULL) unresolved,
            count(compound_id)-count(DISTINCT compound_id) duplicate_structural_records FROM source_records GROUP BY source_id""")
        for name, columns in [
            (
                "skin_assays",
                "assay_id VARCHAR, compound_id VARCHAR, source_id VARCHAR, endpoint VARCHAR, value DOUBLE, raw_value VARCHAR, unit VARCHAR, species VARCHAR, assay VARCHAR, study_id VARCHAR, evidence_type VARCHAR, conditions_json VARCHAR, training_eligible BOOLEAN",
            ),
            (
                "cosmetic_annotations",
                "annotation_id VARCHAR, compound_id VARCHAR, source_id VARCHAR, ingredient_id VARCHAR, name VARCHAR, annotation_json VARCHAR",
            ),
        ]:
            globpath = root / "data" / name / "**/*.parquet"
            if list((root / "data" / name).rglob("*.parquet")):
                con.execute(
                    f"CREATE VIEW {name} AS SELECT r.* FROM read_parquet({quote(globpath)},union_by_name=true) r JOIN sources s USING(source_id) WHERE s.admission='allowed'"
                )
            else:
                con.execute(f"CREATE TABLE {name} ({columns})")
        # Sentinels make these globs valid on a clean first build. Only the publisher
        # may put data in these directories; uncommitted raw batches are never read.
        con.execute(
            f"CREATE VIEW predictions AS SELECT * FROM read_parquet({quote(root / 'data/published_predictions/*.parquet')},union_by_name=true)"
        )
        con.execute(
            f"CREATE VIEW models AS SELECT * FROM read_parquet({quote(root / 'data/models/*.parquet')},union_by_name=true)"
        )
        con.execute(
            f"CREATE VIEW endpoints AS SELECT * FROM read_parquet({quote(root / 'data/endpoints/*.parquet')},union_by_name=true)"
        )
        for file in sorted((root / "data/reference_tables").glob("*.parquet")):
            stem = file.stem
            sid = (
                "cpdat"
                if stem.startswith("cpdat_")
                else "pubchem_inchikey"
                if stem == "pubchem_identifiers"
                else stem
            )
            # Reference inputs have the same admission boundary as structural inputs.
            derived = {
                "structure_identifiers",
                "identifier_conflicts",
                "skin_observation_id_migration",
            }
            if stem not in derived and sid not in selected:
                continue
            name = "ref_" + stem
            con.execute(
                f'CREATE VIEW "{name}" AS SELECT * FROM read_parquet({quote(file)})'
            )
        con.execute("""CREATE VIEW prediction_values AS SELECT p.compound_id,p.model_id,p.status,p.computed_at,
            j.key AS endpoint,try_cast(j.value AS DOUBLE) AS value,e.unit,e.name,e.species,p.quality_flags_json,p.domain_status
            FROM predictions p CROSS JOIN LATERAL json_each(p.values_json) j
            LEFT JOIN endpoints e ON e.model_id=p.model_id AND e.endpoint=j.key""")
        skin_report = root / "reports/skin_nested_validation.json"
        if not skin_report.exists():
            skin_report = root / "reports/skin_model_performance.json"
        if skin_report.exists():
            rs = json.loads(skin_report.read_text())
            con.register(
                "skin_registry",
                pa.Table.from_pylist(
                    [
                        dict(
                            endpoint=k,
                            status=v.get("status", "NOT_ESTABLISHED"),
                            inference_admitted=bool(v.get("inference_admitted", False)),
                            metadata_json=json.dumps(v),
                        )
                        for k, v in rs.items()
                    ]
                ),
            )
            con.execute("CREATE TABLE skin_model_status AS SELECT * FROM skin_registry")
            con.execute("""CREATE VIEW skin_prediction_status AS SELECT compound_id,endpoint,
                CASE WHEN inference_status<>'eligible' THEN inference_status WHEN NOT inference_admitted THEN 'model_not_validated' ELSE 'pending' END status
                FROM compounds CROSS JOIN skin_model_status""")
        chembl_status = root / "reports/chembl_sqlite_import.json"
        if chembl_status.exists() and "chembl_sqlite" in selected:
            report = json.loads(chembl_status.read_text())
            dm = root / "manifests/downloads/chembl_sqlite.json"
            download = json.loads(dm.read_text()) if dm.exists() else {}
            from .references import CHEMBL_TRANSFORM, _chembl_checkpoint

            if (
                report.get("status") == "complete"
                and report.get("lineage_version") == 2
                and report.get("transform_version") == CHEMBL_TRANSFORM
                and report.get("source_sha256") == download.get("sha256")
                and download.get("status") == "verified"
                and selected["chembl_sqlite"].get("sha256", download.get("sha256"))
                == download.get("sha256")
            ):
                _chembl_checkpoint(root, report, download["sha256"])
                for table, files in report.get("output_files", {}).items():
                    if not table.replace("_", "").isalnum():
                        raise ValueError("Invalid ChEMBL table name")
                    parquet_list = "[" + ",".join(quote(root / f) for f in files) + "]"
                    con.execute(
                        f'CREATE VIEW "chembl_{table}" AS SELECT * FROM read_parquet({parquet_list})'
                    )
        con.execute("""CREATE VIEW identifier_conflicts AS SELECT inchikey,count(DISTINCT compound_id) structures
            FROM compounds WHERE inchikey IS NOT NULL AND inchikey<>'' GROUP BY inchikey HAVING count(DISTINCT compound_id)>1""")
        con.execute("""CREATE VIEW structure_purchase_status AS SELECT c.compound_id,
            CASE WHEN o.compound_id IS NOT NULL THEN 'catalog_reported_in_stock' ELSE 'unknown' END purchase_status
            FROM compounds c LEFT JOIN (SELECT DISTINCT compound_id FROM offers WHERE availability='catalog_reported_in_stock') o USING(compound_id)""")
        con.execute(
            "CREATE VIEW skin_experimental AS SELECT * FROM compounds SEMI JOIN (SELECT compound_id FROM skin_assays WHERE evidence_type='experimental') USING(compound_id)"
        )
        con.execute(
            "CREATE VIEW unresolved_skin_ingredients AS SELECT * FROM skin_assays WHERE compound_id IS NULL"
        )
        con.execute(
            "CREATE VIEW cosmetic_related AS SELECT * FROM compounds SEMI JOIN cosmetic_annotations USING(compound_id)"
        )
        counts = dict(
            zip(
                ["compounds", "source_records"],
                [
                    con.execute("SELECT count(*) FROM " + t).fetchone()[0]
                    for t in ["compounds", "source_records"]
                ],
            )
        )
        source_counts = (
            con.execute("SELECT * FROM source_counts ORDER BY source_id")
            .fetchdf()
            .to_dict("records")
        )
        con.execute("CHECKPOINT")
        con.close()
        os.replace(tmp, target)
        # First build may precede any catalog against which batches can be checked.
        # Reconcile publication with this exact snapshot before returning to callers.
        publish(root)
        report = {
            "built_at": now(),
            **counts,
            "source_counts": source_counts,
            "included_ingestions": [i["source_id"] for i in ingests],
            "excluded_ingestions": excluded,
            "admitted_structure_sources_complete": all(
                s["id"] in [i["source_id"] for i in ingests]
                for s in specs
                if s.get("admission") == "allowed"
                and s.get("kind")
                in [
                    "smiles",
                    "chembl_tsv",
                    "coconut_csv",
                    "lotus_metadata",
                    "npatlas_tsv",
                    "sdf",
                    "ice_xlsx",
                    "huskin_csv",
                    "skinpix_xlsx",
                    "normalized_jsonl",
                ]
            ),
            "held_sources": [
                {"source_id": s["id"], "reason": s.get("reason")}
                for s in specs
                if s.get("admission") != "allowed"
            ],
            "scientific_validation": "NOT_ESTABLISHED",
        }
        atomic_json(root / "reports/db_build.json", report)
        return report


def export(root, sql, dest, fmt):
    dest = Path(dest).resolve()
    dest.parent.mkdir(parents=True, exist_ok=True)
    con = connect(root, True)
    if fmt in ["csv", "parquet"]:
        con.execute(
            f"COPY ({sql}) TO {quote(dest)} (FORMAT {fmt}, "
            + ("HEADER true" if fmt == "csv" else "COMPRESSION 'ZSTD'")
            + ")"
        )
    elif fmt == "sdf":
        from rdkit import Chem

        with Chem.SDWriter(str(dest)) as writer:
            cur = con.execute(sql)
            columns = [d[0] for d in cur.description]
            if "identity_smiles" not in columns:
                con.close()
                raise ValueError(
                    "SDF export requires identity_smiles to preserve exact stereochemistry"
                )
            while batch := cur.fetchmany(1000):
                for values in batch:
                    row = dict(zip(columns, values))
                    m = Chem.MolFromSmiles(row["identity_smiles"])
                    if m is None:
                        raise ValueError("Invalid structure in SDF export")
                    for key, value in row.items():
                        if value is not None:
                            m.SetProp(key, str(value))
                    writer.SetForceV3000(bool(m.GetStereoGroups()))
                    writer.write(m)
    con.close()
