'Validate the selected structure and retain explicit calculation provenance.'
from __future__ import annotations

import csv
import hashlib
import json
import math
from contextlib import ExitStack
from pathlib import Path

from ijson.common import ObjectBuilder
import ijson

from .saved_parent import _HashReader, _parse_events, fingerprint, MAX_DEPTH, MAX_RECORD_BYTES

EXPORTS = {"complete": "candidates.csv", "complete_alias": "all_candidates.csv",
           "ranked": "ranked_candidates.csv", "summary": "candidate_summary.csv"}
SCOPE = "delivery identity and exported bytes; scientific validity is not certified"


def _sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError("duplicate delivery metadata key")
        result[key] = value
    return result


def _inspect_result(path, on_row):
    """Validate the entire document, retaining only small metadata and one row."""
    before = fingerprint(path)
    metadata, counts, stack, capture = {}, {}, [], None
    keys = {"marker", "triage", "decision", "target_design", "candidate_exports"}
    with path.open("rb") as stream:
        reader = _HashReader(stream)
        try:
            with _parse_events(reader) as events:
                for event, value in events:
                    reader.last_event = reader.length
                    if event == "map_key":
                        frame = stack[-1]
                        if value in frame["keys"]:
                            raise ValueError("duplicate delivery result key")
                        frame["keys"].add(value)
                        frame["key"] = value
                        position = frame["path"]
                    elif event in {"end_map", "end_array"}:
                        position = stack[-1]["path"]
                    else:
                        if not stack:
                            position = ()
                            if event != "start_map":
                                raise ValueError("delivery result must be an object")
                        else:
                            frame = stack[-1]
                            key = frame["key"] if frame["type"] == "map" else frame["index"]
                            position = (*frame["path"], key)
                            if frame["type"] == "array":
                                frame["index"] += 1
                        if len(position) == 1 and position[0] in {"candidates", "discovery_candidates"}:
                            if event != "start_array":
                                raise ValueError("delivery candidate population must be an array")
                            counts[position[0]] = 0
                        row = len(position) == 2 and position[0] in counts
                        meta = len(position) == 1 and position[0] in keys
                        if row or meta:
                            if row and event != "start_map":
                                raise ValueError("delivery candidate must be an object")
                            if event in {"start_map", "start_array"}:
                                capture = {"path": position, "builder": ObjectBuilder(), "bytes": 0}
                            elif meta:
                                metadata[position[0]] = value
                    if event == "number" and isinstance(value, float) and not math.isfinite(value):
                        raise ValueError("nonfinite delivery result value")
                    if capture is not None:
                        capture["bytes"] += len(value.encode("utf-8")) + 8 if isinstance(value, str) else 16
                        if capture["bytes"] > MAX_RECORD_BYTES:
                            raise ValueError("delivery record exceeds the 16 MiB limit")
                        capture["builder"].event(event, value)
                        if event in {"end_map", "end_array"} and position == capture["path"]:
                            record = capture["builder"].value
                            if len(position) == 1:
                                metadata[position[0]] = record
                            else:
                                on_row(position[0], counts[position[0]], record)
                                counts[position[0]] += 1
                            capture = None
                    if event in {"start_map", "start_array"}:
                        if len(stack) >= MAX_DEPTH:
                            raise ValueError("delivery result nesting exceeds the input limit")
                        stack.append({"path": position, "type": "map" if event == "start_map" else "array",
                                      "key": None, "keys": set(), "index": 0})
                    elif event in {"end_map", "end_array"}:
                        stack.pop()
        except (ijson.JSONError, OverflowError) as exc:
            raise ValueError("invalid delivery result JSON") from exc
    if fingerprint(path) != before:
        raise ValueError("delivery result changed during admission")
    return metadata, counts, reader.digest.hexdigest()


def _cell(value, *, spreadsheet_safe=False):
    if value is None:
        return ""
    if isinstance(value, (dict, list)):
        return json.dumps(value, ensure_ascii=False, sort_keys=True, allow_nan=False)
    if spreadsheet_safe and isinstance(value, str) and value.startswith(("=", "+", "-", "@", "\t", "\r", "\n")):
        return "'" + value
    return str(value)


def _compare_row(reader, record, label, *, spreadsheet_safe=False):
    row = next(reader, None)
    if (row is None or None in row or not reader.fieldnames or "ligand_id" not in reader.fieldnames
            or len(reader.fieldnames) != len(set(reader.fieldnames)) or set(reader.fieldnames) != set(record)):
        raise ValueError(f"incomplete or malformed candidate export: {label}")
    if any(key not in record or value != _cell(record[key], spreadsheet_safe=spreadsheet_safe) for key, value in row.items()):
        raise ValueError(f"candidate export differs from scientific JSON: {label}, {record.get('ligand_id')}")


def validate_screen(screen: Path, expected: dict | None = None) -> dict:
    """Refuse contradictions; explicitly retain missing legacy evidence as NE."""
    screen, expected = Path(screen), expected or {}
    manifest_path, result_path = screen / "run_manifest.json", screen / "result.json"
    tracked = {p: fingerprint(p) for p in (manifest_path, result_path)}
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"), object_pairs_hook=_object)
    if not isinstance(manifest, dict):
        raise ValueError("delivery manifest must be an object")
    contract = manifest.get("candidate_exports") or {}
    modern = isinstance(contract, dict) and contract.get("contract") == "candidate-exports-v2"
    if contract and not modern:
        raise ValueError("unsupported declared candidate export contract")
    files, readers = {}, {}
    with ExitStack() as stack:
        if modern:
            if contract.get("status") != "EXPORTED":
                raise ValueError("declared candidate exports are not complete")
            csv.field_size_limit(max(csv.field_size_limit(), MAX_RECORD_BYTES))
            for role, name in EXPORTS.items():
                if contract.get(role) != name:
                    raise ValueError(f"unexpected candidate export path: {role}")
                path = screen / name
                if not path.is_file() or path.is_symlink():
                    raise ValueError(f"missing declared candidate export: {name}")
                tracked[path] = fingerprint(path)
                files[role] = path
                if role != "complete_alias":
                    readers[role] = csv.DictReader(stack.enter_context(path.open(encoding="utf-8-sig", newline="")))
            for key in ("n_acquired", "n_ranked", "n_summary"):
                if type(contract.get(key)) is not int or contract[key] < 0:
                    raise ValueError(f"invalid export count: {key}")
            limit = contract.get("summary_top")
            if limit is not None and (type(limit) is not int or limit < 0):
                raise ValueError("invalid candidate summary limit")
            expected_summary = min(contract["n_ranked"], limit) if limit else contract["n_ranked"]
            if contract["n_summary"] != expected_summary or contract.get("export_limit") != (limit or None):
                raise ValueError("candidate summary scope differs from declared limit")

        def row(pool, index, record):
            if modern:
                if not isinstance(record.get("ligand_id"), str) or not record["ligand_id"]:
                    raise ValueError("candidate ID is missing from delivery")
                if pool == "discovery_candidates":
                    _compare_row(readers["complete"], record, "complete", spreadsheet_safe=contract.get("csv_cell_policy") == "utf8-bom-formula-escaped-v1")
                else:
                    _compare_row(readers["ranked"], record, "ranked", spreadsheet_safe=contract.get("csv_cell_policy") == "utf8-bom-formula-escaped-v1")
                    if index < contract["n_summary"]:
                        _compare_row(readers["summary"], record, "summary", spreadsheet_safe=contract.get("csv_cell_policy") == "utf8-bom-formula-escaped-v1")

        metadata, counts, result_sha = _inspect_result(result_path, row)
        result_contract = metadata.get("candidate_exports") or {}
        if modern and result_contract != contract:
            raise ValueError("result and manifest candidate export contracts differ")
        if not modern and isinstance(result_contract, dict) and result_contract.get("contract"):
            raise ValueError("result declares candidate exports absent from its manifest")
        if modern:
            if counts.get("discovery_candidates") != contract["n_acquired"] or counts.get("candidates") != contract["n_ranked"]:
                raise ValueError("candidate JSON population differs from declared export count")
            if "n_acquired_candidates" in manifest and manifest["n_acquired_candidates"] != contract["n_acquired"]:
                raise ValueError("manifest acquired population differs from exported population")
            for role, reader in readers.items():
                # Consume the header even when the population is empty.
                if not reader.fieldnames or "ligand_id" not in reader.fieldnames or next(reader, None) is not None:
                    raise ValueError(f"candidate export has missing header or extra rows: {role}")
            hashes = {role: _sha(path) for role, path in files.items()}
            if hashes["complete"] != hashes["complete_alias"]:
                raise ValueError("complete candidate export alias differs")
            exports = {"status": "PASS", "n_acquired": counts["discovery_candidates"],
                       "n_ranked": counts["candidates"], "n_summary": contract["n_summary"], "sha256": hashes}
        else:
            exports = {"status": "NOT_EVALUATED", "reason": "legacy candidate export contract absent"}

    marker = metadata.get("marker")
    if modern and (not isinstance(marker, str) or not marker or manifest.get("screen") != marker):
        raise ValueError("modern delivery marker identity is missing or inconsistent")
    for name, value in (("manifest", manifest.get("screen")), ("coverage", expected.get("marker"))):
        if value is not None and marker is not None and value != marker:
            raise ValueError(f"{name} marker differs from delivered result")
    triage = metadata.get("triage") or {}
    design = metadata.get("target_design") or {}
    decision = metadata.get("decision") or {}
    manifest_design = manifest.get("target_design") or {}
    for value in (triage, design, decision, manifest_design):
        if not isinstance(value, dict):
            raise ValueError("delivery target metadata must be an object")
    node = triage.get("node")
    for value in (design.get("node_gene"), decision.get("node"), manifest_design.get("node_gene")):
        if value is not None and node is not None and str(value).upper() != str(node).upper():
            raise ValueError("conflicting target identities within delivered result/manifest")
    parameters = manifest.get("parameters") or {}
    if not isinstance(parameters, dict):
        raise ValueError("delivery parameters must be an object")
    planned = parameters.get("planned_screen") or {}
    if not isinstance(planned, dict):
        raise ValueError("planned screen must be an object")
    expectations = []
    if "node_group" in expected:
        expectations.append(expected["node_group"])
    elif "node" in expected:
        expectations.append(expected["node"])
    if "node" in planned:
        expectations.append(planned["node"])
    if parameters.get("expected_node") is not None:
        expectations.append(parameters["expected_node"])
    known = [v for v in expectations if v != "?unreadable"]
    if modern and known and node is None and any(v is not None for v in known):
        raise ValueError("planned target cannot be bound to delivered target metadata")
    for value in known:
        if node is not None and (value is None or str(value).upper() != str(node).upper()):
            raise ValueError(f"planned target {value!r} differs from delivered target {node!r}")
    for source in (expected, planned):
        if "node_action" in source and "required_action" in design and source["node_action"] != design["required_action"]:
            raise ValueError("planned action differs from delivered target action")
    for path, before in tracked.items():
        if fingerprint(path) != before:
            raise ValueError("delivery inputs changed during admission")
    return {"contract": "screen-publication-admission-v1", "scope": SCOPE,
            "marker": marker, "node": node, "result_sha256": result_sha,
            "manifest_sha256": _sha(manifest_path), "exports": exports,
            "planned_target": {"status": "PASS" if known and node is not None else "NOT_EVALUATED",
                               "expected": expectations,
                               "reason": "recorded plan and delivered target agree" if known and node is not None
                                         else "plan or delivered target identity absent; no original-run provenance inferred"}}
