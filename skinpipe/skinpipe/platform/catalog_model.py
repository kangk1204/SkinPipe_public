"""Pin official model assets for catalog execution and CPU research fallback."""

from __future__ import annotations

import importlib.metadata
import json
from pathlib import Path
import sys

from .artifacts import atomic_bytes, canonical_json, sha256_file
from .identity import digest

ASSETS = Path(__file__).with_name("_catalog_assets")


def create_profile(output):
    expected = json.loads((ASSETS / "admet-expected.json").read_bytes())
    dist = importlib.metadata.distribution("admet-ai")
    if dist.version != expected["version"]:
        raise ValueError("Expected official ADMET-AI 2.0.1")
    root = Path(dist.locate_file("admet_ai")).resolve()
    implementation = root / "admet_model.py"
    if sha256_file(implementation) != expected["implementation_sha256"]:
        raise ValueError("ADMET implementation checksum differs")
    endpoints = [
        p
        for p in (root / "resources").rglob("*.csv")
        if sha256_file(p) == expected["endpoint_metadata_sha256"]
    ]
    if len(endpoints) != 1:
        raise ValueError("Expected official 52 endpoint definitions")
    for item in expected["weights"]:
        path = root / item["path"]
        if path.stat().st_size != item["bytes"] or sha256_file(path) != item["sha256"]:
            raise ValueError("Official checkpoint bytes differ")
    packages = {
        name: importlib.metadata.version(name)
        for name in ("admet-ai", "torch", "chemprop", "rdkit", "lightning", "numpy")
    }
    if (
        packages["torch"] not in {"2.8.0+cpu", "2.8.0+cu128"}
        or packages["rdkit"] != "2025.9.5"
    ):
        raise ValueError("Use the pinned Ubuntu24 CPU or CUDA 12.8 model environment")
    basis = {
        **expected,
        "packages": packages,
        "device": "cpu",
        "preprocessing": "exact-rdkit-graph-v1",
        "builder_profile_version": 2,
        "catalog_devices": ["cpu", "cuda"]
        if packages["torch"].endswith("+cu128")
        else ["cpu"],
    }
    manifest = {
        **basis,
        "model_id": "admet-ai-2.0.1-local-" + digest(basis)[:16],
        "drugbank_enabled": False,
        "percentiles_enabled": False,
        "scientific_validation": "NOT_ESTABLISHED",
        "applicability_domain": "NOT_ESTABLISHED",
        "runtime_channel": "owner_local_cpu",
        "training_provenance": "Official endpoint source URLs retained; no local biological validation",
    }
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True, mode=0o700)
    manifest_path = output / "model-manifest.json"
    raw = canonical_json(manifest)
    if manifest_path.exists() and manifest_path.read_bytes() != raw:
        raise ValueError("Existing model profile differs; use a new folder")
    atomic_bytes(manifest_path, raw)

    def entry(path):
        return {"path": str(path), "sha256": sha256_file(path)}

    profile = {
        "provider": "admet_ai",
        "version": "2.0.1",
        "preprocessing": "exact-rdkit-graph-v1",
        "python": entry(Path(sys.executable).absolute()),
        "implementation": entry(implementation),
        "models_dir": str(root / "resources/models"),
        "weights": [entry(root / w["path"]) for w in expected["weights"]],
        "model_manifest": entry(manifest_path.absolute()),
        "endpoint_definitions": entry(endpoints[0]),
        "packages": packages,
        "batch_size": 16,
        "max_candidates": 100,
        "budget_seconds": 600,
    }
    from .research_toxicity import _profile

    _profile(profile)
    atomic_bytes(output / "profile.json", canonical_json(profile))
    (output / "profile.json").chmod(0o600)
    return profile


if __name__ == "__main__":
    create_profile(Path(sys.argv[1]).absolute())
