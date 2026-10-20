"""Scientific validation decisions, independent of processing and candidate rank.

This module is also shipped, byte for byte, with the browser bridge. It uses
only the standard library. Legacy flags, numerical scores and control results
never establish any of these five candidate-level claims.
"""
from __future__ import annotations

from collections import Counter
from copy import deepcopy
import hashlib
import json
import re

CONTRACT = "skinpipe-scientific-validation-v1"
STATUSES = ("PASS", "FAIL", "PARTIAL", "NOT_EVALUATED", "NOT_APPLICABLE")
DEFINITIONS = {
    "geometry_validated": "Geometry checks on the selected saved ligand pose and recorded receptor; not affinity or original scoring coordinates.",
    "ranking_validated": "Independent evaluation of a frozen ranking policy on the stated target and held-out compounds; not redocking or finite energy.",
    "direction_supported": "Reviewed, structure- and target-linked evidence agrees with the requested target action; not skin efficacy.",
    "provenance_status": "Independent verification of source records, compound identity and assay context; not a human-origin label or file hash alone.",
    "experimental_validation": "Project-specific experimental validation for the stated endpoint and protocol; not a published activity value or prediction.",
}
POLICY_SHA256 = hashlib.sha256(json.dumps(DEFINITIONS, sort_keys=True,
    separators=(",", ":")).encode()).hexdigest()
GEOMETRY_CHECKS = (
    "bond_lengths", "bond_angles", "internal_steric_clash",
    "aromatic_ring_flatness", "non-aromatic_ring_non-flatness",
    "double_bond_flatness", "internal_energy", "minimum_distance_to_protein",
    "volume_overlap_with_protein", "protein-ligand_maximum_distance",
    "minimum_distance_to_organic_cofactors", "minimum_distance_to_inorganic_cofactors",
    "minimum_distance_to_waters", "volume_overlap_with_organic_cofactors",
    "volume_overlap_with_inorganic_cofactors", "volume_overlap_with_waters",
)
LOADER_CHECKS = ("mol_pred_loaded", "mol_cond_loaded", "sanitization",
    "inchi_convertible", "all_atoms_connected", "no_radicals")


def state(key, status="NOT_EVALUATED", reason="No qualifying validation receipt is attached.", evidence=None):
    if key not in DEFINITIONS or status not in STATUSES:
        raise ValueError("unknown scientific validation decision")
    return {"status": status, "scope": DEFINITIONS[key], "reason": reason,
            "evidence": list(evidence or [])}


def finish(states, scope, counts=None):
    statuses = [states[k]["status"] for k in DEFINITIONS]
    overall = ("REVIEW_REQUIRED" if "FAIL" in statuses else
               "SCOPED_CHECKS_PASSED" if all(s == "PASS" for s in statuses) else "EXPLORATORY")
    result = {"contract": CONTRACT, "policy_sha256": POLICY_SHA256,
              "scope": scope, "overall_status": overall, **states}
    if counts is not None:
        result["counts"] = counts
    return result


def _count_status(count, n):
    return ("NOT_EVALUATED" if not n else "FAIL" if count["FAIL"] else
            "PASS" if count["PASS"] == n else
            "NOT_APPLICABLE" if count["NOT_APPLICABLE"] == n else
            "PARTIAL" if count["PASS"] or count["PARTIAL"] else "NOT_EVALUATED")


def empty_contract(scope="candidate", reason="No qualifying validation receipt is attached."):
    return finish({k: state(k, reason=reason) for k in DEFINITIONS}, scope)


def validate_contract(value):
    """Validate recorded decisions without inventing evidence for legacy data."""
    if not isinstance(value, dict) or value.get("contract") != CONTRACT:
        raise ValueError("unsupported scientific validation contract")
    if value.get("policy_sha256") != POLICY_SHA256 or value.get("scope") not in ("candidate", "population", "campaign"):
        raise ValueError("scientific validation policy or scope does not match")
    for k in DEFINITIONS:
        v = value.get(k)
        if (not isinstance(v, dict) or v.get("status") not in STATUSES
                or v.get("scope") != DEFINITIONS[k] or not isinstance(v.get("reason"), str)
                or not v["reason"] or not isinstance(v.get("evidence"), list)):
            raise ValueError("malformed scientific validation decision: " + k)
        if v["status"] in ("PASS", "FAIL") and not v["evidence"]:
            raise ValueError("evaluated scientific decision lacks its receipt: " + k)
    expected = finish({k: value[k] for k in DEFINITIONS}, value["scope"])["overall_status"]
    if value.get("overall_status") != expected:
        raise ValueError("scientific validation summary contradicts its decisions")
    if "counts" in value:
        n = value.get("n_candidates")
        if type(n) is not int or n < 0 or not isinstance(value["counts"], dict):
            raise ValueError("scientific validation population count is invalid")
        for k in DEFINITIONS:
            count = value["counts"].get(k)
            if (not isinstance(count, dict) or set(count) != set(STATUSES)
                    or any(type(c) is not int or c < 0 for c in count.values())
                    or sum(count.values()) != n or _count_status(count, n) != value[k]["status"]):
                raise ValueError("scientific validation population counts contradict its status")
    return deepcopy(value)


def recorded_contract(record, scope="candidate"):
    value = record.get("validation_status") if isinstance(record, dict) else None
    if value is None:
        return empty_contract(scope, "Legacy record has no scoped scientific validation receipt.")
    try:
        decision = validate_contract(value)
        if decision["scope"] != scope:
            raise ValueError("validation scope does not match this record")
        return decision
    except ValueError:
        return empty_contract(scope, "Attached validation contract is unsupported or internally inconsistent.")


def candidate_contract(row, reviewed_direction=False):
    """The exporter supplies the result of the existing identity-aware review.

    Never use row.direction_supported alone: in legacy rows it also describes
    evidence for an opposing mechanism, and says nothing about the desired action.
    """
    result = empty_contract()
    evidence = [deepcopy(e) for e in (row.get("direction_evidence") or [])
                if isinstance(e, dict) and e.get("status") == "reviewed"
                and e.get("structure_key") and e.get("node_gene")
                and e.get("endpoint") and e.get("source_url")]
    verdict = row.get("direction_verdict")
    if (reviewed_direction is True and row.get("direction_required")
            and evidence and verdict in ("compatible", "opposing")):
        result["direction_supported"] = state("direction_supported",
            "PASS" if verdict == "compatible" else "FAIL",
            "Reviewed evidence agrees with the requested action." if verdict == "compatible" else
            "Reviewed evidence opposes the requested action; retain as a counter-control.", evidence)
    return finish({k: result[k] for k in DEFINITIONS}, "candidate")


def aggregate(contracts, scope="population"):
    counters = {k: Counter() for k in DEFINITIONS}
    n = 0
    for value in contracts:
        value = validate_contract(value)
        n += 1
        for k in DEFINITIONS:
            counters[k][value[k]["status"]] += 1
    states, counts = {}, {}
    for k in DEFINITIONS:
        count = counters[k]
        counts[k] = {s: count[s] for s in STATUSES}
        status = _count_status(count, n)
        states[k] = state(k, status,
            "Every candidate in the stated result population is counted; missing evidence stays explicit.",
            [{"kind": "candidate_state_counts", "n_candidates": n, "counts": counts[k]}])
    result = finish(states, scope, counts)
    result["n_candidates"] = n
    return result


def _hash(value):
    return isinstance(value, str) and re.fullmatch(r"[a-f0-9]{64}", value) is not None


def geometry_from_verdict(document, screen, row, result_sha256, receipt_sha256, artifact_lookup):
    """Read a PoseBusters v2 receipt after binding its actual input bytes.

    artifact_lookup is an adapter restricted to the caller's approved local
    artifacts. This module never opens paths supplied in an uploaded receipt.
    A positive control or one checked compound never certifies other candidates.
    """
    key = "geometry_validated"
    def unknown(why):
        return state(key, reason=why)
    try:
        meta = document["_validation"]
        item = document[screen][row["ligand_id"]]
        inputs = item["inputs"]
        if (meta.get("schema_version") != 2 or not meta.get("posebusters_version")
                or not isinstance(meta.get("config"), dict) or not _hash(meta.get("coverage_sha256"))
                or not _hash(result_sha256) or not _hash(receipt_sha256)
                or meta["screens"][screen].get("ranked_result_sha256") != result_sha256):
            return unknown("Pose validation receipt is not bound to this result snapshot.")
        if (type(row.get("rank")) is not int or row["rank"] < 1
                or type(row.get("dock_pose_index")) is not int or row["dock_pose_index"] < 1
                or type(inputs.get("rank")) is not int or inputs["rank"] != row["rank"]
                or type(inputs.get("pose_index")) is not int or inputs["pose_index"] != row["dock_pose_index"]
                or inputs.get("pose_artifact_status") != "portable_pose_bytes_verified"):
            return unknown("Selected pose identity or original saved-byte binding is missing.")
        verified, snapshots = {}, {}
        for name in ("sdf", "receptor", "docking_record"):
            descriptor = inputs[name]
            raw = artifact_lookup(descriptor)
            if (not _hash(descriptor.get("sha256")) or not isinstance(raw, bytes)
                    or hashlib.sha256(raw).hexdigest() != descriptor["sha256"]):
                return unknown("Pose validation input bytes are unavailable or have changed: " + name)
            verified[name] = deepcopy(descriptor)
            snapshots[name] = raw
        docking = json.loads(snapshots["docking_record"])
        if docking.get("ligand_id") != row["ligand_id"] or not isinstance(docking.get("poses"), list):
            return unknown("The actual docking record belongs to a different compound.")
        poses = [p for p in docking["poses"] if isinstance(p, dict)
                 and type(p.get("pose_index")) is int and p["pose_index"] == row["dock_pose_index"]]
        if len(poses) != 1:
            return unknown("The actual docking record does not identify one selected pose.")
        artifact = poses[0].get("pose_artifact") or {}
        if (artifact.get("sha256") != inputs["sdf"]["sha256"]
                or type(artifact.get("size_bytes")) is not int
                or artifact["size_bytes"] != len(snapshots["sdf"])):
            return unknown("The actual docking record does not bind these SDF bytes.")
        checks = item["checks"]
        if (any(checks.get(c) is not True for c in LOADER_CHECKS)
                or any(type(checks.get(c)) is not bool for c in GEOMETRY_CHECKS)
                or not isinstance(item.get("negative_control_failed_geometry"), list)
                or not any(c in GEOMETRY_CHECKS for c in item["negative_control_failed_geometry"])):
            return unknown("Complete readable geometry checks and a failing geometry control are required.")
        failed = [c for c in GEOMETRY_CHECKS if checks[c] is False]
        observed = "FAIL" if failed else "PASS"
        if item.get("status") != observed or item.get("pose_valid") is not (not failed):
            return unknown("Recorded pose verdict contradicts its individual checks.")
        return state(key, observed, "Recorded saved-pose geometry checks: " + (", ".join(failed) if failed else "all passed"),
            [{"kind": "posebusters_v2", "receipt_sha256": receipt_sha256,
              "result_sha256": result_sha256, "posebusters_version": meta["posebusters_version"],
              "ligand_id": row["ligand_id"], "pose_index": row["dock_pose_index"],
              "inputs": verified, "checks": deepcopy(checks),
              "original_receptor_bytes_verified": inputs.get("original_receptor_bytes_verified") is True}])
    except (KeyError, TypeError, ValueError, OSError, AttributeError):
        return unknown("Pose validation receipt cannot be evaluated against these artifacts.")


def result_view(result):
    """An API view adds decisions while the stored JSON and downloads stay raw."""
    view = dict(result)
    population = view.get("discovery_candidates")
    if not isinstance(population, list):
        known = {r["ligand_id"]: r for r in view.get("candidates", [])}
        for row in view.get("dropped", []):
            if row.get("ligand_id"):
                known.setdefault(row["ligand_id"], row)
        population = list(known.values())
    view["validation_status"] = aggregate(recorded_contract(r) for r in population)
    view["validation_status"]["population_basis"] = (
        "discovery_candidates" if "discovery_candidates" in view else "legacy_ranked_plus_recoverable_dropped")
    return view
