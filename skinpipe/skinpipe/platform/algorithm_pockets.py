"""Optional pocket algorithms; a detected cavity is not a validated target site.

The caller supplies the platform's owned executor after lease/resource admission.
No shell command, executable path, or config code is accepted from a request.
P2Rank: https://github.com/rdk/p2rank (alphafold excludes B-factor features).
Fpocket: https://github.com/Discngine/fpocket/blob/master/doc/GETTINGSTARTED.md
"""

from __future__ import annotations

import csv
from copy import deepcopy
import hashlib
import io
import math
from pathlib import Path
import re
import tempfile

from .artifacts import ArtifactIntegrityError, canonical_json, sha256_file


def _number(value):
    if isinstance(value, bool):
        raise ValueError("Boolean is not a pocket coordinate/score")
    try:
        value = float(value)
    except OverflowError as error:
        raise ValueError("Pocket number exceeds finite range") from error
    if not math.isfinite(value):
        raise ValueError("Nonfinite pocket output")
    return value


def _read(path, root, maximum=16 * 1024 * 1024):
    path, root = Path(path), Path(root).resolve()
    if (
        path.is_symlink()
        or root not in path.resolve().parents
        or not path.is_file()
        or path.stat().st_size > maximum
    ):
        raise ValueError("Missing, unsafe or excessive algorithm output")
    return path.read_bytes()


def _trusted_profile(profile, method):
    """File checks supplement, and do not replace, the common qualification ledger."""
    if not isinstance(profile, dict):
        raise ValueError("Administrator profile must be a record")
    value = deepcopy(profile)
    if (
        value.get("method") != method
        or not value.get("version")
        or not value.get("license_sources")
    ):
        raise ValueError("Versioned administrator profile and license sources required")
    path = Path(value.get("executable", ""))
    if (
        not path.is_absolute()
        or not path.is_file()
        or sha256_file(path) != value.get("executable_sha256")
    ):
        raise ArtifactIntegrityError("Algorithm executable is missing or changed")
    assets = value.get("model_files", [])
    if not isinstance(assets, list) or len(assets) > 10000:
        raise ValueError("Bounded pinned algorithm assets required")
    for asset in assets:
        item = Path(asset["path"])
        if (
            not item.is_absolute()
            or not item.is_file()
            or sha256_file(item) != asset["sha256"]
        ):
            raise ArtifactIntegrityError(
                "Algorithm model/config asset is missing or changed"
            )
    if method == "p2rank" and not assets:
        raise ValueError("P2Rank config, model and distribution files must be pinned")
    if method == "p2rank":
        configs = value.get("configs", {})
        pinned = {item["path"] for item in assets}
        if any(configs.get(name) not in pinned for name in ("default", "alphafold")):
            raise ValueError(
                "Both native default and B-factor-independent P2Rank configs must be pinned"
            )
        java_home = Path(value.get("java_home", ""))
        java = java_home / "bin" / "java"
        if (
            not java_home.is_absolute()
            or str(java) not in pinned
            or not java.is_file()
            or sha256_file(java) != value.get("java_executable_sha256")
        ):
            raise ArtifactIntegrityError(
                "P2Rank requires a pinned Java home and executable"
            )
    return value


def parse_fpocket(directory):
    """Read native *_info.txt and pocketN_vert.pqr, preserving descriptor names."""
    root = Path(directory)
    paths = list(root.glob("*_info.txt"))
    if len(paths) != 1:
        raise ValueError("Expected one fpocket descriptor file")
    text = _read(paths[0], root).decode("utf-8")
    blocks = re.split(r"(?m)^\s*Pocket\s+(\d+)\s*:\s*$", text)
    if len(blocks) == 1:
        if text.strip():
            raise ValueError("Unrecognized fpocket descriptors")
        return []
    pockets, seen = [], set()
    for index in range(1, len(blocks), 2):
        pid = int(blocks[index])
        if pid < 1 or pid in seen or len(seen) >= 10000:
            raise ValueError("Invalid/duplicate fpocket identifier")
        seen.add(pid)
        descriptors = {}
        for line in blocks[index + 1].splitlines():
            if ":" in line:
                key, value = line.split(":", 1)
                descriptors[key.strip()] = _number(value.strip())
        if "Score" not in descriptors:
            raise ValueError("Fpocket score is missing")
        raw = _read(root / "pockets" / f"pocket{pid}_vert.pqr", root)
        points = []
        for line in raw.decode().splitlines():
            if line.startswith(("ATOM  ", "HETATM")):
                # Native PQR has fixed PDB coordinate columns and sphere radius.
                points.append(
                    [_number(line[start : start + 8]) for start in (30, 38, 46)]
                )
        if not points:
            raise ValueError("Fpocket alpha-sphere coordinates are absent")
        low = [min(point[axis] for point in points) for axis in range(3)]
        high = [max(point[axis] for point in points) for axis in range(3)]
        pockets.append(
            {
                "pocket_id": str(pid),
                "rank": pid,
                "score": descriptors["Score"],
                "center": [(a + b) / 2 for a, b in zip(low, high)],
                "proposed_size": [max(10.0, b - a + 10.0) for a, b in zip(low, high)],
                "box_policy": "alpha_sphere_extent_plus_5_angstrom_each_side",
                "descriptors": descriptors,
                "coordinate_unit": "angstrom",
                "vertices_sha256": hashlib.sha256(raw).hexdigest(),
                "sphere_count": len(points),
                "status": "predicted_pocket",
                "scientific_validation": "NOT_ESTABLISHED",
            }
        )
    return sorted(pockets, key=lambda item: (item["rank"], item["pocket_id"]))


def parse_p2rank(csv_bytes):
    """Read the real predictions CSV by trimmed names, including empty results."""
    raw = bytes(csv_bytes)
    if len(raw) > 16 * 1024 * 1024:
        raise ValueError("P2Rank table exceeds output bound")
    reader = csv.DictReader(io.StringIO(raw.decode("utf-8-sig")))
    required = {
        "name",
        "rank",
        "score",
        "probability",
        "center_x",
        "center_y",
        "center_z",
    }
    if not reader.fieldnames or not required.issubset(
        {key.strip() for key in reader.fieldnames}
    ):
        raise ValueError("Unrecognized P2Rank predictions columns")
    pockets, seen = [], set()
    for raw_row in reader:
        if None in raw_row or any(value is None for value in raw_row.values()):
            raise ValueError("Malformed P2Rank CSV row")
        row = {key.strip(): value.strip() for key, value in raw_row.items()}
        name, rank = row["name"], int(row["rank"])
        probability = _number(row["probability"])
        if (
            not name
            or name in seen
            or rank < 1
            or not 0 <= probability <= 1
            or len(seen) >= 10000
        ):
            raise ValueError("Invalid P2Rank identifier/rank/site probability")
        seen.add(name)
        pockets.append(
            {
                "pocket_id": name,
                "rank": rank,
                "score": _number(row["score"]),
                "site_probability": probability,
                "site_probability_scope": "provider pocket calibration; not ligand binding probability",
                "center": [_number(row["center_" + axis]) for axis in "xyz"],
                "proposed_size": [24.0, 24.0, 24.0],
                "box_policy": "24_angstrom_proposal_requires_pocket_review",
                "residue_ids": row.get("residue_ids", "").split(),
                "surface_atom_ids": row.get("surf_atom_ids", "").split(),
                "raw": row,
                "coordinate_unit": "angstrom",
                "status": "predicted_pocket",
                "scientific_validation": "NOT_ESTABLISHED",
            }
        )
    return sorted(pockets, key=lambda item: (item["rank"], item["pocket_id"]))


def pocket_atom_coverage(method, receptor_bytes, execution_log):
    """Describe native atom accounting without inventing atom-level coverage.

    P2Rank reports a structure set and a protein set, but those counts do not
    identify the excluded atoms. HETATM may include PTMs, cofactors and waters;
    preserve their input inventory separately from the provider's accounting.
    Equal counts, missing logs or a successful exit never validate full coverage.
    """
    records = [
        line
        for line in bytes(receptor_bytes).decode("utf-8").splitlines()
        if line.startswith(("ATOM  ", "HETATM"))
    ]
    record_counts = {
        kind: sum(line.startswith(kind) for line in records)
        for kind in ("ATOM  ", "HETATM")
    }
    hetero = {}
    for line in records:
        if line.startswith("HETATM"):
            key = (
                line[21:22].strip(),
                line[22:26].strip(),
                line[26:27].strip(),
                line[17:20].strip(),
            )
            hetero[key] = hetero.get(key, 0) + 1
    inventory = [
        {
            "chain": key[0],
            "author_residue_number": key[1],
            "insertion_code": key[2],
            "residue_name": key[3],
            "atom_record_count": count,
            "model_inclusion": "not_verified",
        }
        for key, count in sorted(hetero.items())
    ]
    log = execution_log[-1024 * 1024 :] if isinstance(execution_log, str) else ""
    warnings = []
    observed = {}
    for name in ("structure", "protein"):
        matches = re.findall(
            r"(?im)^.*?\b" + name + r"\s+atoms\s*:\s*([0-9]{1,10})\s*$", log
        )
        values = {int(value) for value in matches}
        observed[name] = next(iter(values)) if len(values) == 1 else None
        if len(values) > 1:
            warnings.append(
                {
                    "code": "native_atom_counts_conflict",
                    "field": name,
                    "reported_values": sorted(values),
                }
            )
    structure, protein = observed["structure"], observed["protein"]
    excluded, status = None, "not_reported"
    if structure is not None and protein is not None:
        if 0 <= protein <= structure:
            excluded = structure - protein
            status = (
                "partial_native_protein_set" if excluded else "reported_counts_equal"
            )
            if excluded:
                warnings.append(
                    {
                        "code": "native_protein_atom_set_is_partial",
                        "outside_native_protein_atom_count": excluded,
                        "interpretation": "Count difference only; excluded atom identities are not established",
                    }
                )
        else:
            status = "inconsistent_reported_counts"
            warnings.append({"code": "native_atom_counts_inconsistent"})
    else:
        warnings.append(
            {
                "code": "native_atom_coverage_not_reported",
                "interpretation": "Returned log does not establish complete native atom coverage",
            }
        )
    if structure is not None and structure != len(records):
        warnings.append(
            {
                "code": "native_structure_count_differs_from_input_records",
                "interpretation": "Provider preprocessing, models or alternate locations may change atom accounting",
            }
        )
    if inventory:
        warnings.append(
            {
                "code": "input_hetero_residue_coverage_unverified",
                "residue_count": len(inventory),
                "interpretation": "Input HETATM residues may include PTMs, cofactors or waters; atom-level model inclusion is unverified",
            }
        )
    return {
        "provider": method,
        "coverage_status": status,
        "coverage_validated": False,
        "input_atom_record_count": len(records),
        "input_record_counts": {
            key.strip(): value for key, value in record_counts.items()
        },
        "reported_structure_atom_count": structure,
        "reported_protein_atom_count": protein,
        "outside_native_protein_atom_count": excluded,
        "engine_reports_ignoring_ligands": bool(
            re.search(r"\bignoring\s+ligands\b", log, re.I)
        ),
        "input_hetero_residues": inventory[:64],
        "input_hetero_residue_count": len(inventory),
        "hetero_inventory_truncated": len(inventory) > 64,
        "log_sha256": hashlib.sha256(log.encode()).hexdigest() if log else None,
        "log_scope": "returned execution log, potentially tail-truncated; not an atom-level mapping",
        "warnings": warnings,
        "binding_site_validity": "NOT_ESTABLISHED",
    }


def run_pocket_search(request, profile, *, fetch_artifact, publish_artifact, run_owned):
    """Execute an administrator-pinned fpocket/P2Rank through owned process control.

    run_owned(argv, cwd=Path, timeout_seconds=float, env=dict) returns a JSON-safe receipt
    with status=complete and returncode=0 on success. It must own cancellation,
    process groups, lease heartbeat, time/memory accounting and resource limits.
    Request fields: method,receptor_sha256,structural_source,timeout_seconds.
    The caller retains every candidate; this result changes no ranking.
    """
    result = {
        "status": "held",
        "pockets": [],
        "candidate_population_changed": False,
        "scientific_validation": "NOT_ESTABLISHED",
    }
    try:
        if set(request) - {
            "method",
            "receptor_sha256",
            "structural_source",
            "timeout_seconds",
        }:
            raise ValueError(
                "Pocket requests accept artifact identities and bounded parameters only"
            )
        method = request["method"]
        if method not in {"fpocket", "p2rank"}:
            raise ValueError("Unsupported pocket algorithm")
        selected = _trusted_profile(profile, method)
        timeout = _number(request.get("timeout_seconds", 120))
        if not 0 < timeout <= 3600 or not callable(run_owned):
            raise ValueError("Owned execution and a bounded timeout are required")
        raw = bytes(fetch_artifact(request["receptor_sha256"]))
        if (
            hashlib.sha256(raw).hexdigest() != request["receptor_sha256"]
            or len(raw) > 16 * 1024 * 1024
        ):
            raise ArtifactIntegrityError("Receptor hash or size mismatch")
        atoms = [
            line
            for line in raw.decode().splitlines()
            if line.startswith(("ATOM  ", "HETATM"))
        ]
        if not atoms:
            raise ValueError("An observed PDB coordinate artifact is required")
        for line in atoms:
            for start in (30, 38, 46):
                _number(line[start : start + 8])
        with tempfile.TemporaryDirectory(prefix="skinpipe-pocket-") as temporary:
            work = Path(temporary)
            receptor = work / "receptor.pdb"
            receptor.write_bytes(raw)
            if method == "fpocket":
                argv = [selected["executable"], "-f", str(receptor)]
                config = None
            else:
                source = (
                    str(request.get("structural_source", "unknown"))
                    .lower()
                    .replace("-", "_")
                )
                config = (
                    "default"
                    if source in {"xray", "x_ray", "x_ray_diffraction"}
                    else "alphafold"
                )
                argv = [
                    selected["executable"],
                    "predict",
                    "-f",
                    str(receptor),
                    "-o",
                    str(work / "predictions"),
                    "-c",
                    selected["configs"][config],
                    "-visualizations",
                    "0",
                    "-threads",
                    "1",
                ]
            execution_options = {"cwd": work, "timeout_seconds": timeout}
            if method == "p2rank":
                execution_options["env"] = {
                    "JAVA_HOME": selected["java_home"],
                    "JAVA_OPTS": "",
                }
            execution = run_owned(argv, **execution_options)
            atom_coverage = pocket_atom_coverage(method, raw, execution.get("log"))
            result.update(
                execution=deepcopy(execution),
                argv=argv,
                profile_sha256=hashlib.sha256(canonical_json(selected)).hexdigest(),
                receptor_sha256=request["receptor_sha256"],
                p2rank_config=config,
                atom_coverage=atom_coverage,
            )
            if (
                execution.get("status") != "complete"
                or execution.get("returncode") != 0
            ):
                result.update(
                    status=execution.get("status", "failed"),
                    reason="Owned pocket process did not complete",
                )
                return result
            if _trusted_profile(profile, method) != selected:
                raise ArtifactIntegrityError("Pocket profile changed during execution")
            output = work / ("receptor_out" if method == "fpocket" else "predictions")
            pockets = (
                parse_fpocket(output)
                if method == "fpocket"
                else parse_p2rank(_read(output / "receptor.pdb_predictions.csv", work))
            )
            paths = (
                list(output.glob("*_info.txt"))
                + list((output / "pockets").glob("*_vert.pqr"))
                if method == "fpocket"
                else list(output.glob("*.csv"))
            )
            raw_artifacts, total = [], 0
            for path in sorted(paths):
                raw_output = _read(path, output)
                total += len(raw_output)
                if total > 64 * 1024 * 1024:
                    raise ValueError("Pocket output bundle exceeds64 MiB")
                published = publish_artifact(raw_output)
                if published.get("sha256") != hashlib.sha256(raw_output).hexdigest():
                    raise ArtifactIntegrityError(
                        "Native pocket output publication mismatch"
                    )
                raw_artifacts.append(
                    {"name": path.relative_to(output).as_posix(), **published}
                )
            data = canonical_json(
                {
                    "method": method,
                    "pockets": pockets,
                    "input_sha256": request["receptor_sha256"],
                    "raw_artifacts": raw_artifacts,
                    "atom_coverage": atom_coverage,
                }
            )
            artifact = publish_artifact(data)
            if artifact.get("sha256") != hashlib.sha256(data).hexdigest():
                raise ArtifactIntegrityError("Pocket result publication mismatch")
            result.update(
                status="complete",
                pockets=pockets,
                result_artifact=artifact,
                raw_artifacts=raw_artifacts,
                box_admission="requires_target_site_review",
                ranking_weight=0,
            )
    except (OSError, ValueError, KeyError, TypeError, UnicodeError) as error:
        result.update(status="held", reason=str(error))
    return result
