"""Synthetic CAS fixtures: raw pose QC and reused row evidence survive export."""
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from skinpipe.platform.artifacts import ArtifactStore, canonical_json
from skinpipe.platform.catalog import Catalog
from skinpipe.platform.exports import RunExports


def build_export(tmp_path, row, store):
    catalog = Catalog(tmp_path / "catalog")
    manifest = catalog.create_snapshot([{
        "candidate_id": "synthetic-candidate", "target_id": "compound:synthetic",
        "smiles": "CCO", "role": "candidate", **row,
    }], provenance={"owner_id": "researcher"})
    service = SimpleNamespace(root=tmp_path, catalog=catalog, artifacts=store)
    directory = tmp_path / "export"
    receipt = RunExports(service).build({"snapshot_id": manifest["snapshot_id"],
        "target_id": "compound:synthetic", "format": "xlsx"}, directory)
    exported = json.loads((directory / receipt["chunks"][0]["file"]).read_bytes())[0]
    return directory, receipt, exported


def assert_asset(directory, receipt, digest, raw):
    relative = "assets/" + digest
    assert (directory / relative).read_bytes() == raw
    assert receipt["files"][relative] == {
        "sha256": hashlib.sha256(raw).hexdigest(), "size_bytes": len(raw)}


def test_export_keeps_raw_pose_qc_and_cached_source_row_without_expanding_cache(tmp_path):
    store = ArtifactStore(tmp_path / "cas")
    # Synthetic check names intentionally do not claim an actual PoseBusters run.
    metrics = {"checks": {f"synthetic_check_{i}": True for i in range(22)},
               "pose_check_status": "pass", "pose_checks_pass": True,
               "missing_binary_checks": [], "inputs": {"fixture": True}}
    raw_qc = json.dumps(metrics, ensure_ascii=False, indent=2).encode() + b"\n"
    qc = store.put_bytes(raw_qc, media_type="application/json")
    nested_raw = b"synthetic historical evidence must not be recursively exported"
    nested = store.put_bytes(nested_raw)
    cached_raw = canonical_json({"candidate_id": "synthetic-source", "discovery_run_id": "synthetic-run",
        "smiles": "CCO", "source_record_sha256": nested["sha256"], "fixture": True})
    cached = store.put_bytes(cached_raw, media_type="application/json")
    row = {"routine_pose_qc": "pass", "routine_pose_quality": {
        "status": "complete", "candidate_id": "synthetic-candidate", "metrics": metrics,
        "artifacts": [{**qc, "path": "/untrusted/never-read/qc.json"}]},
        "routine_reuse": {"source_candidate_id": "synthetic-source", "source_discovery_run_id": "synthetic-run",
                          "cached_row_artifact_sha256": cached["sha256"], "exact_structure_and_protocol": True}}
    directory, receipt, exported = build_export(tmp_path, row, store)
    assert_asset(directory, receipt, qc["sha256"], raw_qc)
    assert_asset(directory, receipt, cached["sha256"], cached_raw)
    assert not (directory / "assets" / nested["sha256"]).exists()
    assert exported["routine_pose_quality"] == row["routine_pose_quality"]
    assert exported["routine_reuse"] == row["routine_reuse"]
    assert not [w for w in receipt["warnings"] if w.get("sha256") in {qc["sha256"], cached["sha256"]}]


@pytest.mark.parametrize("failure", ["missing", "corrupt"])
def test_export_warns_and_preserves_missing_or_corrupt_raw_evidence(tmp_path, failure):
    store = ArtifactStore(tmp_path / "cas")
    qc = store.put_bytes(b"synthetic raw pose QC")
    cached = store.put_bytes(b"synthetic cached row")
    for record in [qc, cached]:
        path = store.verified_path(record["sha256"])
        if failure == "missing":
            path.unlink()
        else:
            path.write_bytes(b"synthetic corruption")
    row = {"routine_pose_quality": {"status": "complete", "metrics": {"fixture": True},
                                   "artifacts": [{**qc, "path": "/never/read"}]},
           "routine_reuse": {"cached_row_artifact_sha256": cached["sha256"]}}
    directory, receipt, exported = build_export(tmp_path, row, store)
    for record, component in [(qc, "routine_pose_quality_artifact"), (cached, "cached_row_artifact_sha256")]:
        assert not (directory / "assets" / record["sha256"]).exists()
        assert "assets/" + record["sha256"] not in receipt["files"]
        assert any(w.get("sha256") == record["sha256"] and w.get("reason") == "artifact_quarantined"
                   and w.get("component") == component and w.get("candidate_id") == "synthetic-candidate"
                   for w in receipt["warnings"])
    assert exported["routine_pose_quality"] == row["routine_pose_quality"]
    assert exported["routine_reuse"] == row["routine_reuse"]


@pytest.mark.parametrize("artifact_container", ["list", "malformed_mapping"])
def test_export_raw_qc_never_reads_declared_paths_or_arbitrary_hashes(tmp_path, monkeypatch, artifact_container):
    store = ArtifactStore(tmp_path / "cas")
    unrelated = store.put_bytes(b"synthetic unrelated hash must not be exported")
    external = tmp_path / "outside-cas.txt"
    external.write_bytes(b"synthetic external path must never be read")
    entry = {"path": str(external), "sha256": str(external)}
    artifacts = [entry, None, "not-an-artifact", {"sha256": unrelated["sha256"]}]
    # The valid unrelated CAS digest is deliberately outside routine_pose_quality.
    row = {"unrelated": {"artifacts": [{"sha256": unrelated["sha256"], "path": str(external)}]},
           "routine_pose_quality": {"metrics": {"fixture": True},
                "artifacts": artifacts[:3] if artifact_container == "list" else entry}}
    original_open = Path.open
    def guarded_open(path, *args, **kwargs):
        if path == external:
            raise AssertionError("Exporter followed an untrusted artifact path")
        return original_open(path, *args, **kwargs)
    with monkeypatch.context() as context:
        context.setattr(Path, "open", guarded_open)
        directory, receipt, exported = build_export(tmp_path, row, store)
    assert not (directory / "assets" / unrelated["sha256"]).exists()
    assert exported["routine_pose_quality"] == row["routine_pose_quality"]
    assert exported["unrelated"] == row["unrelated"]
    assert external.read_bytes() == b"synthetic external path must never be read"
