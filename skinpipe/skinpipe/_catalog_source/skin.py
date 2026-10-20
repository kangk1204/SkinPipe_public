"""Experimental endpoint preparation. Mixtures, censoring and provenance survive."""

import csv
import json
import math
import re
import openpyxl
import pyarrow as pa
import pyarrow.parquet as pq
from .chemistry import normalize
from .common import atomic_json, digest, now, sha256


def ice_rows(path):
    workbook = openpyxl.load_workbook(path, read_only=True, data_only=True)
    try:
        for sheet in workbook:
            if not sheet.title.startswith("Data"):
                continue
            sheet.reset_dimensions()
            it = iter(sheet.values)
            header = next(it)
            for source_row, values in enumerate(it, 2):
                yield (
                    sheet.title,
                    {**dict(zip(header, values)), "_source_row": source_row},
                )
    finally:
        workbook.close()


def permeability_rows(path, kind):
    if kind == "huskin_csv":
        with open(path) as f:
            for i, row in enumerate(csv.DictReader(f), 1):
                yield {**row, "_row": i}
    else:
        w = openpyxl.load_workbook(path, read_only=True, data_only=True)
        try:
            s = w.active
            s.reset_dimensions()  # Published file incorrectly declares A1 dimension.
            it = iter(s.values)
            header = next(it)
            for i, values in enumerate(it, 1):
                yield {**dict(zip(header, values)), "_row": i}
        finally:
            w.close()


def numeric(v):
    try:
        x = float(v)
        return x if math.isfinite(x) else None
    except (TypeError, ValueError):
        return None


def logkp(value, unit):
    value = float(value)
    if unit == "log10(cm/s)":
        return value
    if unit == "log10(cm/h)":
        return value - math.log10(3600)
    if value <= 0:
        raise ValueError("Permeability must be positive")
    if unit == "cm/h":
        return math.log10(value / 3600)
    if unit == "cm/s":
        return math.log10(value)
    raise ValueError("Unrecognized permeability unit")


def study_id(row):
    text = str(row.get("doi", row.get("DOI", row.get("Reference", ""))))
    match = re.search(r"10\.\d{4,9}/[^\s;]+", text)
    return match.group().lower().rstrip(".;,") if match else text.strip().lower()


def prepare_skin(root):
    specs = json.loads((root / "config/sources.json").read_text())
    report = {}
    allrows = []
    cosmetics = []
    for spec in specs:
        kind = spec.get("kind")
        sid = spec["id"]
        if kind not in ["ice_xlsx", "huskin_csv", "skinpix_xlsx"]:
            continue
        receipt = root / "manifests/downloads" / f"{sid}.json"
        if not receipt.exists():
            continue
        manifest = json.loads(receipt.read_text())
        if manifest["status"] != "verified":
            continue
        path = root / manifest["path"]
        if sha256(path) != manifest["sha256"]:
            raise ValueError("Changed source")
        observations = []
        iterator = (
            ice_rows(path)
            if kind == "ice_xlsx"
            else (("", r) for r in permeability_rows(path, kind))
        )
        for source_ordinal, (sheet, row) in enumerate(iterator, 1):
            if sid == "ice_functional_use" and not any(
                x in str(row.get("Response", "")).lower()
                for x in [
                    "skin",
                    "hair",
                    "fragrance",
                    "humectant",
                    "emollient",
                    "sunscreen",
                    "uv absorber",
                    "nail",
                    "cosmetic",
                    "oral care",
                ]
            ):
                continue
            mixture = str(row.get("Mixture", "")).lower() == "mixture"
            raw_smiles = row.get("SMILES", row.get("Smiles"))
            norm = normalize(None if mixture else raw_smiles)
            ordinal = str(row.get("Record_ID", row.get("data ID", row.get("_row"))))
            legacy_id = digest(sid + "|" + manifest["sha256"] + "|" + ordinal)
            physical_row = str(row.get("_source_row", row.get("_row", source_ordinal)))
            rid = digest(
                sid
                + "|"
                + manifest["sha256"]
                + "|skin-observation-v2|"
                + sheet
                + "|"
                + physical_row
            )
            value = numeric(row.get("Response"))
            raw = row.get("Response")
            unit = row.get("Response_Unit")
            endpoint = row.get("Endpoint")
            species = row.get("Species")
            evidence = (
                "experimental"
                if row.get("Data_Type")
                in ["Measured", "In Vivo", "In Vitro", "In Chemico"]
                else "predicted_or_unspecified"
            )
            eligible = (
                evidence == "experimental"
                and not mixture
                and norm["inference_status"] == "eligible"
            )
            relation = row.get("Response_Modifier")
            if relation and str(relation).strip() not in ["=", "=="]:
                eligible = False
            if kind in ["huskin_csv", "skinpix_xlsx"]:
                endpoint = "skin_permeability_logkp"
                species = "Human"
                evidence = "experimental"
                unit = "log10(cm/s)"
                raw = row.get("logkp (cm/s)", row.get("logKp (cm/s) (converted)"))
                value = numeric(raw)
                relation = row.get("Kp relation", "=")
                eligible = (
                    norm["inference_status"] == "eligible"
                    and relation == "="
                    and value is not None
                )
                if kind == "skinpix_xlsx":
                    kp = numeric(row.get("Kp (cm/h)"))
                    if (
                        kp is not None
                        and kp > 0
                        and value is not None
                        and abs(logkp(kp, "cm/h") - value) > 0.015
                    ):
                        eligible = False
                        row["_unit_conversion_conflict"] = True
            if sid == "ice_functional_use":
                cosmetic_functions = [
                    "skin",
                    "hair",
                    "fragrance",
                    "humectant",
                    "emollient",
                    "sunscreen",
                    "uv absorber",
                    "nail",
                    "cosmetic",
                    "oral care",
                ]
                if not any(
                    x in str(row.get("Response", "")).lower()
                    for x in cosmetic_functions
                ):
                    continue
                cosmetics.append(
                    dict(
                        annotation_id=rid,
                        compound_id=norm["compound_id"],
                        source_id=sid,
                        ingredient_id=str(row.get("DTXSID", "")),
                        name=row.get("Chemical_Name"),
                        annotation_json=json.dumps(row, default=str),
                        evidence_type=evidence,
                    )
                )
                continue
            observations.append(
                dict(
                    assay_id=rid,
                    observation_id=rid,
                    legacy_assay_id=legacy_id,
                    observation_id_scheme="source-sha256-sheet-physical-row-v2",
                    source_sheet=sheet,
                    source_row=physical_row,
                    compound_id=norm["compound_id"],
                    source_id=sid,
                    source_record_id=ordinal,
                    canonical_smiles=norm["canonical_smiles"],
                    endpoint=endpoint,
                    value=value,
                    raw_value=None if raw is None else str(raw),
                    unit=unit,
                    species=species,
                    assay=row.get("Assay", "human ex vivo skin permeation"),
                    study_id=study_id(row),
                    evidence_type=evidence,
                    conditions_json=json.dumps(row, default=str, ensure_ascii=False),
                    training_eligible=eligible,
                    relation=relation,
                    mixture=mixture,
                    source_sha256=manifest["sha256"],
                )
            )
        if observations:
            if len({r["observation_id"] for r in observations}) != len(observations):
                raise ValueError("Duplicate skin observation identity: " + sid)
            schema = pa.schema(
                [
                    (
                        k,
                        pa.bool_()
                        if k in ["training_eligible", "mixture"]
                        else pa.float64()
                        if k == "value"
                        else pa.string(),
                    )
                    for k in observations[0]
                ]
            )
            out = root / "data/skin_assays" / f"{sid}.parquet"
            out.parent.mkdir(parents=True, exist_ok=True)
            tmp = out.with_suffix(".part")
            pq.write_table(
                pa.Table.from_pylist(observations, schema=schema),
                tmp,
                compression="zstd",
            )
            tmp.replace(out)
            allrows.extend(
                {
                    k: r[k]
                    for k in [
                        "observation_id",
                        "legacy_assay_id",
                        "source_id",
                        "source_record_id",
                        "source_sha256",
                        "source_sheet",
                        "source_row",
                        "assay",
                        "endpoint",
                        "species",
                    ]
                }
                for r in observations
            )
        report[sid] = dict(
            observations=len(observations),
            eligible=sum(r["training_eligible"] for r in observations),
            mixtures=sum(r["mixture"] for r in observations),
        )
    if any(s["id"] == "ice_functional_use" for s in specs):
        out = root / "data/cosmetic_annotations/ice_functional_use.parquet"
        out.parent.mkdir(parents=True, exist_ok=True)
        schema = pa.schema(
            [
                (k, pa.string())
                for k in [
                    "annotation_id",
                    "compound_id",
                    "source_id",
                    "ingredient_id",
                    "name",
                    "annotation_json",
                    "evidence_type",
                ]
            ]
        )
        tmp = out.with_suffix(".part")
        pq.write_table(
            pa.Table.from_pylist(cosmetics, schema=schema), tmp, compression="zstd"
        )
        tmp.replace(out)
    if allrows:
        migration = root / "data/reference_tables/skin_observation_id_migration.parquet"
        migration.parent.mkdir(parents=True, exist_ok=True)
        tmp = migration.with_suffix(".part")
        pq.write_table(pa.Table.from_pylist(allrows), tmp, compression="zstd")
        tmp.replace(migration)
    atomic_json(
        root / "reports/skin_preparation.json",
        {
            "created_at": now(),
            "sources": report,
            "cosmetic_function_records": len(cosmetics),
        },
    )
    return report
