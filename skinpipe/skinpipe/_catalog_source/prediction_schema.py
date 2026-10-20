"""Typed empty files keep SQL views live before the first committed batch."""

import os
import pyarrow as pa
import pyarrow.parquet as pq

PREDICTION_SCHEMA = pa.schema(
    [
        (k, pa.string())
        for k in (
            "compound_id",
            "model_id",
            "input_smiles",
            "input_hash",
            "computed_at",
            "status",
            "error",
            "values_json",
            "quality_flags_json",
            "domain_status",
        )
    ]
)
MODEL_SCHEMA = pa.schema([("model_id", pa.string()), ("metadata_json", pa.string())])
ENDPOINT_SCHEMA = pa.schema(
    [
        (k, pa.string())
        for k in (
            "model_id",
            "endpoint",
            "name",
            "unit",
            "category",
            "task_type",
            "species",
            "source_url",
            "published_metrics_not_local_validation",
        )
    ]
)


def ensure_empty_relations(root):
    for name, schema in [
        ("published_predictions", PREDICTION_SCHEMA),
        ("models", MODEL_SCHEMA),
        ("endpoints", ENDPOINT_SCHEMA),
    ]:
        folder = root / "data" / name
        folder.mkdir(parents=True, exist_ok=True)
        path = folder / "_empty.parquet"
        if (
            path.exists()
            and pq.read_schema(path).equals(schema)
            and pq.ParquetFile(path).metadata.num_rows == 0
        ):
            continue
        tmp = folder / f"_empty.{os.getpid()}.part"
        pq.write_table(pa.Table.from_pylist([], schema=schema), tmp)
        tmp.replace(path)
