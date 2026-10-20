'Administrator-owned pocket runtime qualification with exact source/CAS proof.'

from __future__ import annotations

import math
from pathlib import Path, PurePosixPath
import tempfile
import time

from .algorithm_pockets import _trusted_profile, parse_fpocket, parse_p2rank
from .algorithm_qualification import (
    LICENSES,
    SOURCE_SCHEMA,
    QualificationHeld,
    _asset,
    _cas,
    _digest,
    _hash,
    _json,
    _manifest,
)
from .artifacts import ArtifactIntegrityError, atomic_bytes, canonical_json, sha256_file
from .identity import digest

SCHEMA = "skinpipe-pocket-qualification-v1"
ASSET_SCHEMA = "skinpipe-pocket-assets-v1"
PROVIDERS = {"fpocket", "p2rank"}
MEMORY_SCOPES = {
    "sampled_leader_process_RSS_and_GPU_VRAM",
    "sampled_process_tree_RSS",
    "sampled_process_group_RSS_and_GPU_VRAM",
}


def _positive(value):
    try:
        return type(value) in (int, float) and math.isfinite(value) and value > 0
    except OverflowError:
        return False


def pocket_asset_identity(provider, profile):
    """Verify executable, adapter, configs/models, Java and exact license texts.

    Extend a run_pocket_search profile with runtime_files[{path,sha256}],
    license_receipts[{path,sha256}], and preprocessing_version. P2Rank also needs
    java_version and model_weight_sha256[]. License records contain component,
    version, scope (code/model/runtime), SPDX license, source_url, license_text
    {path,sha256}, and covers_assets[] of exact asset SHA256s. These are operator
    assets, never request-supplied acceptance flags.
    """
    if provider not in PROVIDERS:
        raise ValueError("Unsupported pocket provider")
    selected = _trusted_profile(profile, provider)
    if (
        not isinstance(selected.get("preprocessing_version"), str)
        or not selected["preprocessing_version"]
    ):
        raise QualificationHeld("pinned_pocket_preprocessing_version_required")
    executable = _asset(
        selected["executable"], selected["executable_sha256"], executable=True
    )
    runtime = selected.get("runtime_files")
    if (
        not isinstance(runtime, list)
        or not 1 <= len(runtime) <= 10000
        or any(not isinstance(entry, dict) for entry in runtime)
    ):
        raise QualificationHeld("pocket_implementation_runtime_manifest_required")
    files = [_asset(entry.get("path"), entry.get("sha256")) for entry in runtime]
    paths = [entry["path"] for entry in files]
    if len(set(paths)) != len(paths):
        raise ValueError("Duplicate pocket runtime asset")
    adapter = str(Path(__file__).with_name("algorithm_pockets.py").resolve())
    if adapter not in paths:
        raise QualificationHeld("pocket_adapter_implementation_hash_required")
    assets = [
        _asset(
            entry["path"],
            entry["sha256"],
            executable=entry["path"]
            == str(Path(selected.get("java_home", "")) / "bin/java"),
        )
        for entry in selected.get("model_files", [])
    ]
    all_hashes = {executable["sha256"], *(entry["sha256"] for entry in assets)}
    required = (
        {("fpocket", "code")}
        if provider == "fpocket"
        else {("p2rank", "code"), ("p2rank-models", "model"), ("java", "runtime")}
    )
    receipts, covered, components = [], {}, set()
    entries = selected.get("license_receipts")
    if not isinstance(entries, list) or not 1 <= len(entries) <= 50:
        raise QualificationHeld("pocket_component_license_receipts_required")
    for entry in entries:
        if not isinstance(entry, dict):
            raise ValueError("Pocket license receipt must be a record")
        receipt_file = _asset(entry.get("path"), entry.get("sha256"))
        if receipt_file["size_bytes"] > 256 * 1024:
            raise ValueError("Pocket license receipt exceeds256KiB")
        record = _json(Path(receipt_file["path"]).read_bytes())
        key = (record.get("component"), record.get("scope"))
        accepted = LICENSES | {"GPL-2.0-only WITH Classpath-exception-2.0"}
        if (
            key not in required
            or key in components
            or record.get("license") not in accepted
            or not isinstance(record.get("source_url"), str)
            or not record["source_url"].startswith("https://")
        ):
            raise QualificationHeld(
                "pocket_license_component_scope_or_source_unresolved"
            )
        expected_version = (
            selected.get("java_version") if key[0] == "java" else selected["version"]
        )
        if not expected_version or record.get("version") != expected_version:
            raise QualificationHeld("pocket_license_version_mismatch")
        text_ref = record.get("license_text")
        if not isinstance(text_ref, dict):
            raise QualificationHeld("pocket_exact_license_text_required")
        text = _asset(text_ref.get("path"), text_ref.get("sha256"))
        if not 1 <= text["size_bytes"] <= 1024 * 1024:
            raise QualificationHeld("pocket_license_text_empty_or_excessive")
        hashes = record.get("covers_assets")
        if (
            not isinstance(hashes, list)
            or not hashes
            or any(_digest(value) not in all_hashes for value in hashes)
        ):
            raise QualificationHeld("license_must_cover_exact_pocket_assets")
        components.add(key)
        covered[key] = set(hashes)
        receipts.append(
            {
                "component": key[0],
                "scope": key[1],
                "license": record["license"],
                "version": expected_version,
                "source_url": record["source_url"],
                "receipt": receipt_file,
                "license_text": text,
                "covers_assets": sorted(set(hashes)),
            }
        )
    if (
        components != required
        or executable["sha256"] not in covered[(provider, "code")]
    ):
        raise QualificationHeld(
            "separate_pocket_code_model_Java_license_receipts_required"
        )
    if provider == "p2rank":
        weights = selected.get("model_weight_sha256")
        if (
            not isinstance(weights, list)
            or not weights
            or not set(weights) <= covered[("p2rank-models", "model")]
        ):
            raise QualificationHeld("exact_p2rank_model_weight_license_required")
        java_root = Path(selected["java_home"]).resolve()
        java_assets = {
            entry["sha256"]
            for entry in assets
            if java_root in Path(entry["path"]).parents
        }
        if not java_assets or not java_assets <= covered[("java", "runtime")]:
            raise QualificationHeld("exact_Java_runtime_and_notices_coverage_required")
        nonjava = {entry["sha256"] for entry in assets} - java_assets
        if not nonjava <= (
            covered[("p2rank", "code")] | covered[("p2rank-models", "model")]
        ):
            raise QualificationHeld(
                "exact_p2rank_distribution_license_coverage_required"
            )
    value = {
        "schema": ASSET_SCHEMA,
        "provider": provider,
        "version": selected["version"],
        "preprocessing_version": selected["preprocessing_version"],
        "executable": executable,
        "runtime_files": sorted(files, key=lambda item: item["path"]),
        "model_files": sorted(assets, key=lambda item: item["path"]),
        "configs": selected.get("configs", {}),
        "java_home": selected.get("java_home"),
        "java_version": selected.get("java_version"),
        "model_weight_sha256": sorted(selected.get("model_weight_sha256", [])),
        "license_receipts": sorted(
            receipts, key=lambda item: (item["component"], item["scope"])
        ),
    }
    return {**value, "sha256": _hash(value)}


def _source(provider, input_sha256, source_receipt_sha256, fetch_artifact):
    payload = _json(_cas(input_sha256, fetch_artifact))
    source = _json(_cas(source_receipt_sha256, fetch_artifact))
    if (
        source.get("schema") != SOURCE_SCHEMA
        or source.get("provider") != provider
        or payload.get("provider") != provider
        or source.get("input_sha256") != input_sha256
        or source.get("campaign_id") != payload.get("campaign_id")
        or not isinstance(source.get("campaign_id"), str)
        or not source["campaign_id"]
        or source.get("synthetic") is not False
        or source.get("source_type")
        not in {"frozen_campaign", "frozen_artifact_canary"}
    ):
        raise QualificationHeld("registered_actual_pocket_source_required")
    manifest = _manifest(source, fetch_artifact)


    if payload.get("campaign_manifest_sha256") not in {
        source["manifest_sha256"],
        digest(manifest),
    }:
        raise QualificationHeld("pocket_input_manifest_identity_mismatch")
    rows = payload.get("rows")
    if not isinstance(rows, list) or rows:
        raise QualificationHeld(
            "protein_only_pocket_input_requires_empty_candidate_rows"
        )
    context = payload.get("context")
    receptor = manifest.get("receptor")
    if not isinstance(context, dict) or not isinstance(receptor, dict):
        raise QualificationHeld("registered_protein_context_required")
    target, protein = context.get("target_id"), context.get("actual_protein_id")
    if (
        not isinstance(target, str)
        or not target
        or manifest.get("plan", {}).get("target_id") != target
        or not isinstance(protein, str)
        or not protein.startswith("uniprot:")
        or receptor.get("actual_protein_id") != protein
    ):
        raise QualificationHeld("pocket_target_or_actual_protein_manifest_mismatch")
    target_context = manifest.get("plan", {}).get("target_context") or {}
    if target_context.get("actual_protein_id", protein) != protein:
        raise QualificationHeld("pocket_target_context_protein_mismatch")
    artifacts = payload.get("context_artifacts")
    if not isinstance(artifacts, dict) or "receptor" not in artifacts:
        raise QualificationHeld("registered_receptor_coordinate_artifact_required")
    entries = source.get("source_artifacts")
    if not isinstance(entries, list) or not 1 <= len(entries) <= 100:
        raise QualificationHeld("actual_pocket_source_artifact_bytes_required")
    hashes, verified = set(), []
    for entry in entries:
        if not isinstance(entry, dict) or not isinstance(entry.get("role"), str):
            raise ValueError("Source artifact requires role and SHA256")
        raw = _cas(entry.get("sha256"), fetch_artifact)
        if not raw:
            raise QualificationHeld("empty_pocket_source_artifact")
        hashes.add(entry["sha256"])
        verified.append(
            {"role": entry["role"], "sha256": entry["sha256"], "size_bytes": len(raw)}
        )
    if not set(artifacts.values()) <= hashes:
        raise QualificationHeld("pocket_context_artifact_not_in_registered_source")
    card_sha = receptor.get("identity_card_sha256")
    card = _json(_cas(card_sha, fetch_artifact))
    record_sha = receptor.get("uniprot_record_sha256")
    record = _json(_cas(record_sha, fetch_artifact))
    model = card.get("identity", {}).get("target_model") or {}
    accession, chain = model.get("accession"), model.get("chain")
    if (
        protein != "uniprot:" + str(accession)
        or record.get("primaryAccession") != accession
        or not isinstance(chain, str)
        or not chain
        or type(record.get("organism", {}).get("taxonId")) is not int
    ):
        raise QualificationHeld("pocket_source_accession_chain_or_taxon_unverified")
    card_artifacts = card.get("artifacts")
    if not isinstance(card_artifacts, dict):
        raise QualificationHeld("pocket_native_source_identity_card_required")
    clean = [card_artifacts.get(key, {}) for key in ("clean_pdb", "clean", "pdb")]
    if not any(
        entry.get("present") is True and entry.get("sha256") == artifacts["receptor"]
        for entry in clean
    ):
        raise QualificationHeld("pocket_coordinates_do_not_match_registered_receptor")
    source_cif = card_artifacts.get("assembly_source", {})
    if source_cif.get("present") is not True:
        raise QualificationHeld("deposited_source_mapping_required_for_pocket_context")
    raw_cif = _cas(source_cif.get("sha256"), fetch_artifact)
    try:
        from ..structure.identity import source_entities, assess
        from Bio.PDB.MMCIF2Dict import MMCIF2Dict

        with tempfile.TemporaryDirectory(prefix="skinpipe-pocket-source-") as directory:
            path = Path(directory) / "source.cif"
            path.write_bytes(raw_cif)
            entities = source_entities(path, [])
            experiment_methods = MMCIF2Dict(str(path)).get("_exptl.method", [])
            claimed_source = (
                str(context.get("structural_source", "unknown"))
                .lower()
                .replace("-", "_")
            )
            if (
                claimed_source in {"xray", "x_ray", "x_ray_diffraction"}
                and "X-RAY DIFFRACTION" not in experiment_methods
            ):
                raise QualificationHeld(
                    "xray_pocket_profile_requires_observed_diffraction_source"
                )
        matching = [
            entry
            for entry in entities
            if any(
                segment.get("chain") == chain and segment.get("accession") == accession
                for segment in entry.get("reference_segments", [])
            )
        ]
        if len(matching) != 1:
            raise QualificationHeld("pocket_chain_mapping_is_not_uniquely_observed")
        actual = assess(
            matching[0], chain, accession, model.get("target_role", "protein")
        )
    except ImportError as error:
        raise QualificationHeld(
            "receptor_identity_parser_runtime_unavailable"
        ) from error
    if model.get("policy") not in {
        "observed_reference_target_role_v1",
        "observed_reference_target_role_v2",
    }:
        raise QualificationHeld("pocket_source_identity_policy_unrecognized")
    keys = [
        "accession",
        "chain",
        "reference_segments",
        "construct_sequence_length",
        "target_role",
    ]
    if "target_construct_residues" in model:
        keys.append("target_construct_residues")
    if any(model.get(key) != actual.get(key) for key in keys):
        raise QualificationHeld(
            "pocket_target_identity_disagrees_with_deposited_source"
        )
    return payload, {
        "sha256": source_receipt_sha256,
        "manifest_sha256": source["manifest_sha256"],
        "source_type": source["source_type"],
        "synthetic": False,
        "target_id": target,
        "actual_protein_id": protein,
        "chain": chain,
        "identity_card_sha256": card_sha,
        "uniprot_record_sha256": record_sha,
        "source_cif_sha256": source_cif["sha256"],
        "receptor_sha256": artifacts["receptor"],
        "source_artifacts": verified,
        "binding_site_validity": "NOT_ESTABLISHED",
    }


def _observed(provider, identity, payload, output_sha256, fetch_artifact):
    output = _json(_cas(output_sha256, fetch_artifact))
    receptor_sha = payload["context_artifacts"]["receptor"]
    if (
        output.get("provider") != provider
        or output.get("status") != "complete"
        or output.get("input_sha256") != _hash(payload)
        or output.get("receptor_sha256") != receptor_sha
        or output.get("asset_identity_sha256") != identity["sha256"]
        or output.get("candidates") != []
        or output.get("ranking_weight") != 0
        or output.get("scientific_validation") != "NOT_ESTABLISHED"
        or output.get("box_admission") != "requires_target_site_review"
    ):
        raise QualificationHeld(
            "successful_pocket_output_not_bound_to_source_and_assets"
        )
    execution = output.get("execution")
    if (
        not isinstance(execution, dict)
        or execution.get("status") != "complete"
        or type(execution.get("returncode")) is not int
        or execution["returncode"] != 0
        or not _positive(execution.get("elapsed_seconds"))
        or not _positive(execution.get("process_peak_rss_bytes"))
        or execution.get("memory_scope") not in MEMORY_SCOPES
    ):
        raise QualificationHeld("actual_pocket_time_and_RSS_measurements_required")
    if provider == "p2rank" and execution["memory_scope"] not in {
        "sampled_process_tree_RSS",
        "sampled_process_group_RSS_and_GPU_VRAM",
    }:
        raise QualificationHeld(
            "p2rank_Java_child_requires_process_tree_RSS_measurement"
        )
    if (
        execution.get(
            "process_tree_peak_rss_bytes", execution["process_peak_rss_bytes"]
        )
        != execution["process_peak_rss_bytes"]
    ):
        raise QualificationHeld("pocket_memory_measurement_alias_mismatch")
    artifacts = output.get("raw_artifacts")
    if not isinstance(artifacts, list) or not 1 <= len(artifacts) <= 20000:
        raise QualificationHeld("native_pocket_output_artifacts_required")
    with tempfile.TemporaryDirectory(prefix="skinpipe-pocket-output-") as directory:
        root = Path(directory)
        names, total = set(), 0
        for entry in artifacts:
            if not isinstance(entry, dict):
                raise ArtifactIntegrityError("Native pocket artifact must be a record")
            name = entry.get("name")
            relative = PurePosixPath(name) if isinstance(name, str) else None
            if (
                relative is None
                or relative.is_absolute()
                or ".." in relative.parts
                or "\\" in name
                or not name
                or name in names
            ):
                raise ArtifactIntegrityError(
                    "Unsafe/duplicate native pocket output name"
                )
            names.add(name)
            raw = _cas(entry.get("sha256"), fetch_artifact)
            total += len(raw)
            if total > 64 * 1024 * 1024 or entry.get("size_bytes") != len(raw):
                raise ArtifactIntegrityError("Native pocket output size mismatch")
            path = root / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(raw)
        observed = (
            parse_fpocket(root)
            if provider == "fpocket"
            else parse_p2rank((root / "receptor.pdb_predictions.csv").read_bytes())
        )
    if output.get("pockets") != observed:
        raise ArtifactIntegrityError(
            "Published pockets differ from native output reparsing"
        )
    if provider == "p2rank":
        source = (
            str(payload["context"].get("structural_source", "unknown"))
            .lower()
            .replace("-", "_")
        )
        config = (
            "default"
            if source in {"xray", "x_ray", "x_ray_diffraction"}
            else "alphafold"
        )
        if output.get("p2rank_config") != config or execution.get("env") != {
            "JAVA_HOME": identity["java_home"],
            "JAVA_OPTS": "",
        }:
            raise QualificationHeld(
                "observed_p2rank_config_or_Java_environment_mismatch"
            )
    return {
        "elapsed_seconds": execution["elapsed_seconds"],
        "peak_memory_mib": execution["process_peak_rss_bytes"] / 1024**2,
        "memory_scope": execution["memory_scope"],
        "pocket_count": len(observed),
    }


def create_pocket_qualification(
    provider,
    profile,
    *,
    input_sha256,
    output_sha256,
    fetch_artifact,
    destination,
    source_receipt_sha256,
):
    identity = pocket_asset_identity(provider, profile)
    payload, source = _source(
        provider, input_sha256, source_receipt_sha256, fetch_artifact
    )
    measured = _observed(provider, identity, payload, output_sha256, fetch_artifact)
    receipt = {
        "schema": SCHEMA,
        "status": "qualified",
        "provider": provider,
        "asset_identity_sha256": identity["sha256"],
        "input_sha256": input_sha256,
        "output_sha256": output_sha256,
        "source_receipt_sha256": source_receipt_sha256,
        "source": source,
        **measured,
        "observed_at_unix": time.time(),
        "scientific_validation": "NOT_ESTABLISHED",
        "binding_site_validity": "NOT_ESTABLISHED",
        "scientific_promotion_allowed": False,
    }
    path = Path(destination).resolve()
    if path.exists():
        raise ValueError(
            "Qualification receipts are immutable; use a fresh destination"
        )
    atomic_bytes(path, canonical_json(receipt))
    return {
        "path": str(path),
        "sha256": sha256_file(path),
        "size_bytes": path.stat().st_size,
        "status": "qualified",
        "asset_identity_sha256": identity["sha256"],
    }


def verify_pocket_qualification(
    provider,
    profile,
    *,
    fetch_artifact,
    qualification_mode=False,
    input_sha256=None,
    source_receipt_sha256=None,
):
    if type(qualification_mode) is not bool:
        raise ValueError("Qualification mode must be boolean")
    identity = pocket_asset_identity(provider, profile)
    if qualification_mode:
        if input_sha256 is None or source_receipt_sha256 is None:
            raise QualificationHeld("registered_actual_pocket_pilot_source_required")
        _, source = _source(
            provider, input_sha256, source_receipt_sha256, fetch_artifact
        )
        return {
            **identity,
            "qualification_status": "pilot_source_verified",
            "source": source,
        }
    ref = profile.get("qualification_receipt")
    if not isinstance(ref, dict):
        raise QualificationHeld("real_input_pocket_runtime_qualification_required")
    asset = _asset(ref.get("path"), ref.get("sha256"))
    if asset["size_bytes"] > 256 * 1024:
        raise ValueError("Pocket qualification receipt exceeds256KiB")
    receipt = _json(Path(asset["path"]).read_bytes())
    if (
        receipt.get("schema") != SCHEMA
        or receipt.get("status") != "qualified"
        or receipt.get("provider") != provider
        or receipt.get("asset_identity_sha256") != identity["sha256"]
        or receipt.get("scientific_validation") != "NOT_ESTABLISHED"
        or receipt.get("binding_site_validity") != "NOT_ESTABLISHED"
        or receipt.get("scientific_promotion_allowed") is not False
    ):
        raise QualificationHeld("pocket_runtime_identity_or_scientific_claim_invalid")
    payload, source = _source(
        provider,
        receipt.get("input_sha256"),
        receipt.get("source_receipt_sha256"),
        fetch_artifact,
    )
    measured = _observed(
        provider, identity, payload, receipt.get("output_sha256"), fetch_artifact
    )
    if receipt.get("source") != source or any(
        receipt.get(key) != value for key, value in measured.items()
    ):
        raise QualificationHeld("pocket_qualification_source_or_measurements_changed")
    if input_sha256 is not None or source_receipt_sha256 is not None:
        _source(provider, input_sha256, source_receipt_sha256, fetch_artifact)
    return {
        **identity,
        "qualification_status": "qualified",
        "qualification_receipt": asset,
        "binding_site_validity": "NOT_ESTABLISHED",
    }
