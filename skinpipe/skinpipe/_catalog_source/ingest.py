import csv
import gzip
import io
import itertools
import json
import multiprocessing as mp
import sys
import zipfile
import pyarrow as pa
import pyarrow.parquet as pq

from rdkit import Chem
from .chemistry import normalize_record, NORMALIZER_VERSION
from .common import atomic_json, digest, lock, now, sha256

csv.field_size_limit(sys.maxsize)

STRING_FIELDS = [
    "source_id",
    "source_version",
    "source_record_id",
    "raw_smiles",
    "name",
    "raw_metadata",
    "record_id",
    "compound_id",
    "canonical_smiles",
    "identity_smiles",
    "inchikey",
    "parent_smiles",
    "parent_id",
    "formula",
    "structure_status",
    "inference_status",
    "normalizer",
    "rights_status",
]
SCHEMA = pa.schema(
    [(x, pa.string()) for x in STRING_FIELDS]
    + [
        (x, pa.int64())
        for x in ["ordinal", "fragments", "formal_charge", "heavy_atoms"]
    ]
    + [("mw", pa.float64())]
)


def text_file(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def records(path, spec):
    kind = spec["kind"]
    if kind == "smiles":
        with text_file(path) as f:
            for line in f:
                vals = line.rstrip().split()
                if not vals:
                    continue
                yield {
                    "raw_smiles": vals[0],
                    "source_record_id": vals[1] if len(vals) > 1 else "",
                    "raw_metadata": json.dumps({"columns": vals[2:]}),
                }
    elif kind == "chembl_tsv":
        with text_file(path) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                yield {
                    "raw_smiles": row["canonical_smiles"],
                    "source_record_id": row["chembl_id"],
                    "raw_metadata": json.dumps(row),
                }
    elif kind == "npatlas_tsv":
        with text_file(path) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                yield {
                    "raw_smiles": row.get("compound_smiles"),
                    "source_record_id": row.get("npaid", ""),
                    "name": row.get("compound_names"),
                    "raw_metadata": json.dumps(row),
                }
    elif kind in ["coconut_csv", "lotus_metadata"]:
        if kind == "coconut_csv":
            with zipfile.ZipFile(path) as z:
                names = [
                    n
                    for n in z.namelist()
                    if n.endswith(".csv") and not n.startswith("__MACOSX")
                ]
                if len(names) != 1:
                    raise ValueError(f"Expected one CSV: {names}")
                with z.open(names[0]) as f:
                    yield from csv_records(
                        io.TextIOWrapper(f, encoding="utf-8-sig"), kind
                    )
        else:
            with text_file(path) as f:
                yield from csv_records(f, kind)
    elif kind == "sdf":
        with (
            gzip.open(path, "rb")
            if str(path).endswith(".gz")
            else open(path, "rb") as f
        ):
            for idx, m in enumerate(Chem.ForwardSDMolSupplier(f, removeHs=False)):
                if m is None:
                    yield {
                        "source_record_id": str(idx + 1),
                        "raw_smiles": "!INVALID_SDF_RECORD",
                        "raw_metadata": "{}",
                    }
                else:
                    props = m.GetPropsAsDict()
                    yield {
                        "source_record_id": str(
                            props.get(
                                "ChEBI ID",
                                props.get("chebi_id", props.get("_Name", idx + 1)),
                            )
                        ),
                        "raw_smiles": Chem.MolToCXSmiles(
                            m,
                            Chem.SmilesWriteParams(),
                            Chem.CXSmilesFields.CX_ENHANCEDSTEREO,
                        ),
                        "name": str(
                            props.get("ChEBI Name", props.get("chebi_name", ""))
                        ),
                        "raw_metadata": json.dumps(props, default=str),
                    }
    elif kind == "normalized_jsonl":
        with open(path) as f:
            for line in f:
                yield json.loads(line)
    elif kind == "ice_xlsx":
        from .skin import ice_rows

        for sheet, row in ice_rows(path):
            # A constituent SMILES is not the identity of a tested formulation.
            mixture = str(row.get("Mixture", "")).lower() == "mixture"
            yield {
                "source_record_id": str(row.get("Record_ID", "")),
                "name": row.get("Formulation_Name")
                if mixture
                else row.get("Chemical_Name"),
                "raw_smiles": None if mixture else row.get("SMILES"),
                "raw_metadata": json.dumps({"sheet": sheet, **row}, default=str),
            }
    elif kind in ["huskin_csv", "skinpix_xlsx"]:
        from .skin import permeability_rows

        for row in permeability_rows(path, kind):
            yield {
                "source_record_id": str(row.get("data ID", row["_row"])),
                "name": row.get("compound name", row.get("Compound name")),
                "raw_smiles": row.get("SMILES", row.get("Smiles")),
                "raw_metadata": json.dumps(row, default=str),
            }
    else:
        raise ValueError(f"No compound adapter for {kind}")


def csv_records(f, kind):
    # LOTUS releases are TSV despite the .csv suffix; detect from the actual header.
    first = f.readline()
    sep = "\t" if first.count("\t") > first.count(",") else ","
    reader = csv.DictReader(itertools.chain([first], f), delimiter=sep)
    for row in reader:
        fields = {k.lower(): k for k in row if k is not None}

        def pick(*keys):
            return next(
                (row[fields[k]] for k in keys if k in fields and row[fields[k]]), None
            )

        smi = pick(
            "canonical_smiles",
            "isomeric_smiles",
            "smiles",
            "smiles_wikidata",
            "structure_smiles",
            "smiles (isomeric)",
        )
        if not any("smiles" in k for k in fields):
            raise ValueError(f"No recognized SMILES column in {list(fields)}")
        sid = pick(
            "identifier",
            "lotus_id",
            "coconut_id",
            "id",
            "wikidata_id",
            "wikidata",
            "structure_wikidata",
        )
        meta = json.dumps(row, ensure_ascii=False)
        rights = "source_terms_admitted"
        if kind == "coconut_csv" and any(
            v in meta.lower() for v in ["npatlas", "natural products atlas"]
        ):
            rights = "upstream_npatlas_version_unresolved"
        yield {
            "raw_smiles": smi,
            "source_record_id": sid or "",
            "name": pick("name", "iupac_name", "compound_name"),
            "raw_metadata": meta,
            "rights_status": rights,
        }


def adapter_version(spec):
    # Bump this semantic version when an adapter changes its emitted records.
    return {"sdf": "sdf-enhanced-stereo-v2", "ice_xlsx": "ice-xlsx-source-row-v2"}.get(
        spec["kind"], "adapter-v1"
    )


def expected_runhash(spec, input_sha256):
    """Stable for unchanged adapters so existing catalogue shards remain reusable."""
    return digest(input_sha256 + "|" + NORMALIZER_VERSION + "|" + adapter_version(spec))


def _batch_hash(batch):
    return digest(
        json.dumps(
            batch,
            sort_keys=True,
            ensure_ascii=False,
            default=str,
            separators=(",", ":"),
        )
    )


def _verify_shard(shard, receipt, batch, sid, spec, input_sha256, batch_size):
    """Bind legacy receipts to the actual source interval before accepting them."""
    stats = json.loads(receipt.read_text())
    if sha256(shard) != stats["sha256"]:
        raise ValueError("Corrupted committed shard: " + str(shard))
    expected = dict(
        input_sha256=input_sha256,
        source_kind=spec["kind"],
        adapter_version=adapter_version(spec),
        normalizer=NORMALIZER_VERSION,
        ordinal_start=batch[0][0],
        ordinal_end=batch[-1][0],
        input_batch_sha256=_batch_hash(batch),
    )
    for key, value in expected.items():
        if key in stats and stats[key] != value:
            raise ValueError("Incompatible ingestion shard " + key + ": " + str(shard))
    if stats["records"] != len(batch):
        raise ValueError("Incompatible ingestion batch geometry: " + str(shard))
    if not all(key in stats for key in expected):
        columns = [
            "ordinal",
            "record_id",
            "source_id",
            "source_version",
            "raw_smiles",
            "source_record_id",
            "name",
            "raw_metadata",
            "rights_status",
            "normalizer",
        ]
        rows = pq.read_table(shard, columns=columns).to_pylist()
        if len(rows) != len(batch):
            raise ValueError("Ingestion shard row count mismatch: " + str(shard))
        for actual, (ordinal, raw) in zip(rows, batch):
            wanted = dict(
                ordinal=ordinal,
                record_id=digest(sid + "|" + input_sha256 + "|" + str(ordinal)),
                source_id=sid,
                source_version=spec["version"],
                normalizer=NORMALIZER_VERSION,
                rights_status=raw.get("rights_status", "source_terms_admitted"),
            )
            wanted.update(
                {
                    key: raw.get(key)
                    for key in [
                        "raw_smiles",
                        "source_record_id",
                        "name",
                        "raw_metadata",
                    ]
                }
            )
            if any(actual[key] != value for key, value in wanted.items()):
                raise ValueError(
                    "Legacy ingestion shard does not match source interval: "
                    + str(shard)
                )
        stats.update(
            expected, batch_size=batch_size, receipt_version=2, migrated_at=now()
        )
        atomic_json(receipt, stats)
    return stats


def ingest(root, sid, workers=6, batch_size=25000):
    if batch_size < 1:
        raise ValueError("batch_size must be positive")
    spec = next(
        s
        for s in json.loads((root / "config/sources.json").read_text())
        if s["id"] == sid
    )
    if spec["admission"] != "allowed":
        raise ValueError("Source is held")
    manifest = json.loads((root / "manifests/downloads" / f"{sid}.json").read_text())
    if manifest["status"] != "verified":
        raise ValueError("Unverified download")
    path = root / manifest["path"]
    if sha256(path) != manifest["sha256"]:
        raise ValueError("Raw file changed since download")
    adapter = adapter_version(spec)
    runhash = expected_runhash(spec, manifest["sha256"])
    out = root / "data/records" / sid / runhash[:16]
    out.mkdir(parents=True, exist_ok=True)
    status_path = root / "manifests/ingest" / f"{sid}.json"
    with lock(root / "logs" / f"ingest-{sid}.lock"):
        if status_path.exists():
            state = json.loads(status_path.read_text())
            if state.get("runhash") == runhash and state.get("status") == "complete":
                iterator = enumerate(records(path, spec), 1)
                recorded_batch_size = state.get("batch_size", batch_size)
                if "batch_size" not in state and len(state["shards"]) > 1:
                    first_receipt = (root / state["shards"][0]).with_suffix(".json")
                    if first_receipt.exists():
                        recorded_batch_size = json.loads(first_receipt.read_text())[
                            "records"
                        ]
                for name in state["shards"]:
                    shard = root / name
                    receipt = shard.with_suffix(".json")
                    if (
                        not receipt.exists()
                        or sha256(shard) != json.loads(receipt.read_text())["sha256"]
                    ):
                        raise ValueError("Corrupted completed ingestion shard: " + name)
                    stats = json.loads(receipt.read_text())
                    batch = list(itertools.islice(iterator, stats["records"]))
                    if not batch:
                        raise ValueError("Completed ingestion has extra rows: " + name)
                    _verify_shard(
                        shard,
                        receipt,
                        batch,
                        sid,
                        spec,
                        manifest["sha256"],
                        recorded_batch_size,
                    )
                if next(iterator, None) is not None:
                    raise ValueError("Completed ingestion has missing source rows")
                state.update(
                    batch_size=recorded_batch_size,
                    normalizer=NORMALIZER_VERSION,
                    adapter_version=adapter,
                    source_kind=spec["kind"],
                    source_version=spec["version"],
                    input_sha256=manifest["sha256"],
                )
                expected = spec.get("expected_records")
                state["expected_records"] = expected
                state["count_reconciliation"] = (
                    "match"
                    if expected == state["records"]
                    else ("not_provided" if expected is None else "mismatch")
                )
                state["status"] = (
                    "complete"
                    if state["count_reconciliation"] != "mismatch"
                    else "count_mismatch"
                )
                atomic_json(status_path, state)
                return state
            if (
                state.get("runhash") == runhash
                and state.get("batch_size", batch_size) != batch_size
            ):
                raise ValueError(
                    "Cannot resume ingestion with changed batch_size; use the recorded batch geometry"
                )
        # Detect incompatible legacy geometry before replacing the running state.
        first = out / "part-000000.parquet"
        if first.exists() and first.with_suffix(".json").exists():
            batch = list(
                itertools.islice(enumerate(records(path, spec), 1), batch_size)
            )
            if not batch:
                raise ValueError("Existing ingestion shard has no source records")
            _verify_shard(
                first,
                first.with_suffix(".json"),
                batch,
                sid,
                spec,
                manifest["sha256"],
                batch_size,
            )
        state = dict(
            source_id=sid,
            runhash=runhash,
            started_at=now(),
            status="running",
            input_sha256=manifest["sha256"],
            records=0,
            valid=0,
            missing_or_invalid=0,
            shards=[],
            relative_directory=str(out.relative_to(root)),
            batch_size=batch_size,
            normalizer=NORMALIZER_VERSION,
            adapter_version=adapter,
            source_kind=spec["kind"],
            source_version=spec["version"],
        )
        atomic_json(status_path, state)
        iterator = enumerate(records(path, spec), 1)
        with mp.get_context("spawn").Pool(workers) as pool:
            for batch_index in itertools.count():
                batch = list(itertools.islice(iterator, batch_size))
                if not batch:
                    break
                shard = out / f"part-{batch_index:06d}.parquet"
                receipt = shard.with_suffix(".json")
                if shard.exists() and receipt.exists():
                    stats = _verify_shard(
                        shard, receipt, batch, sid, spec, manifest["sha256"], batch_size
                    )
                else:
                    rows = []
                    for ordinal, row in batch:
                        rows.append(
                            {
                                **row,
                                "ordinal": ordinal,
                                "source_id": sid,
                                "source_version": spec["version"],
                                "record_id": digest(
                                    sid + "|" + manifest["sha256"] + "|" + str(ordinal)
                                ),
                                "rights_status": row.get(
                                    "rights_status", "source_terms_admitted"
                                ),
                            }
                        )
                    rows = list(pool.imap(normalize_record, rows, chunksize=128))
                    table = pa.Table.from_pylist(rows, schema=SCHEMA)
                    tmp = shard.with_suffix(".parquet.part")
                    pq.write_table(table, tmp, compression="zstd")
                    tmp.replace(shard)
                    stats = dict(
                        records=len(rows),
                        valid=sum(r["structure_status"] == "valid" for r in rows),
                        sha256=sha256(shard),
                        input_sha256=manifest["sha256"],
                        source_kind=spec["kind"],
                        adapter_version=adapter,
                        normalizer=NORMALIZER_VERSION,
                        ordinal_start=batch[0][0],
                        ordinal_end=batch[-1][0],
                        input_batch_sha256=_batch_hash(batch),
                        batch_size=batch_size,
                        receipt_version=2,
                    )
                    atomic_json(receipt, stats)
                state["records"] += stats["records"]
                state["valid"] += stats["valid"]
                state["missing_or_invalid"] = state["records"] - state["valid"]
                state["shards"].append(str(shard.relative_to(root)))
                state["updated_at"] = now()
                atomic_json(status_path, state)
                print(
                    json.dumps(
                        {
                            k: state[k]
                            for k in ["source_id", "records", "valid", "updated_at"]
                        }
                    ),
                    flush=True,
                )
        expected = spec.get("expected_records")
        state["expected_records"] = expected
        state["count_reconciliation"] = (
            "match"
            if expected == state["records"]
            else ("not_provided" if expected is None else "mismatch")
        )
        state["status"] = (
            "complete"
            if state["count_reconciliation"] != "mismatch"
            else "count_mismatch"
        )
        state["finished_at"] = now()
        atomic_json(status_path, state)
        return state
