"""Provider-specific commands for owned workers, with asset-first admission.

The caller supplies an administrator profile and materializes hash-bound inputs
inside a private work directory. Neither function reads executable paths, model
files, flags or qualification assertions from a candidate. Process launch,
leases, resource reservation and launch counters belong to the worker; pass the
returned argv/env to its owned runner, never to a shell.

CLI sources: github.com/gnina/gnina; github.com/ohuelab/boltzina;
github.com/jwohlwend/boltz; github.com/chaidiscovery/chai-lab;
docs.openadmet.org/en/latest/_api/cli/predict.html;
chemprop.readthedocs.io/en/latest/tutorial/cli/predict.html;
github.com/dptech-corp/Uni-Dock2; github.com/OpenFreeEnergy/openfe.
These interfaces must be qualified against the pinned installed versions.
"""
from __future__ import annotations

from collections import Counter
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import zipfile

from .artifacts import ArtifactIntegrityError, canonical_json, sha256_file

PROVIDERS = frozenset({"gnina", "boltzina", "boltz2", "chai1", "openadmet", "chemeleon", "unidock2", "openfe"})
MODEL_PROVIDERS = PROVIDERS - {"unidock2", "openfe"}
CONTRACT_VERSION = "artifact-model-command-v1"
_SHA = re.compile(r"[0-9a-f]{64}")
_MANIFEST = "model-input-binding.json"
_MAX_INPUT = 64 * 1024 * 1024
_LICENSES = {"MIT", "Apache-2.0", "BSD-2-Clause", "BSD-3-Clause", "GPL-2.0-only", "GPL-2.0-or-later",
             "GPL-3.0-only", "GPL-3.0-or-later", "LGPL-2.1-only", "LGPL-2.1-or-later"}


class ModelHeld(ValueError):
    """An unavailable or unadmitted branch; candidate handoff remains possible."""


def _hash(value):
    return hashlib.sha256(canonical_json(value)).hexdigest()


def _pairs(items):
    result = {}
    for key, value in items:
        if key in result:
            raise ValueError("Duplicate JSON field: " + key)
        result[key] = value
    return result


def _json(path, maximum=_MAX_INPUT):
    if path.stat().st_size > maximum:
        raise ValueError("JSON artifact exceeds size bound")
    return json.loads(path.read_text(), object_pairs_hook=_pairs,
                      parse_constant=lambda value: (_ for _ in ()).throw(ValueError("Nonfinite JSON: " + value)))


def _asset(path, expected, *, executable=False):
    if not isinstance(path, str) or not Path(path).is_absolute() or not isinstance(expected, str) or not _SHA.fullmatch(expected):
        raise ValueError("Administrator asset requires an absolute path and SHA256")
    supplied = Path(path)
    resolved = supplied.resolve()
    if not resolved.is_file():
        raise ModelHeld("configured_asset_missing")
    if not executable and supplied.is_symlink():
        raise ValueError("Model and receipt assets may not be symlinks")
    if executable and not os.access(resolved, os.X_OK):
        raise ModelHeld("configured_executable_not_executable")
    if sha256_file(resolved) != expected:
        raise ArtifactIntegrityError("Configured asset checksum changed")
    return resolved


def _asset_list(entries, label):
    if not isinstance(entries, list) or len(entries) > 100000:
        raise ValueError("Bounded " + label + " manifest required")
    verified = {}
    for item in entries:
        if not isinstance(item, dict):
            raise ValueError("Invalid asset manifest entry")
        path = _asset(item.get("path"), item.get("sha256"))
        if path in verified:
            raise ValueError("Duplicate configured asset path")
        verified[path] = item["sha256"]
    return verified


def model_asset_identity(provider, profile):
    """Verify bytes and return the identity to bind a runtime qualification to."""
    if provider not in PROVIDERS or not isinstance(profile, dict):
        raise ValueError("Unsupported provider or invalid administrator profile")
    if profile.get("provider", provider) != provider:
        raise ValueError("Profile provider mismatch")
    for name in ("version", "preprocessing_version"):
        if not isinstance(profile.get(name), str) or not profile[name] or len(profile[name]) > 200:
            raise ModelHeld("pinned_runtime_versions_required")
    binary = _asset(profile.get("executable"), profile.get("executable_sha256"), executable=True)
    models = _asset_list(profile.get("model_files", []), "model file")
    runtime = _asset_list(profile.get("runtime_files", []), "runtime file")
    if provider in MODEL_PROVIDERS and not models:
        raise ModelHeld("pinned_model_assets_required")
    identity = {"provider": provider, "version": profile["version"],
                "preprocessing_version": profile["preprocessing_version"],
                "executable": str(binary), "executable_sha256": profile["executable_sha256"],
                "model_files": [{"path": str(p), "sha256": h} for p, h in sorted(models.items())],
                "runtime_files": [{"path": str(p), "sha256": h} for p, h in sorted(runtime.items())],
                "contract_version": CONTRACT_VERSION,
                "adapter_source_sha256": sha256_file(Path(__file__)),
                "configuration": {key: profile[key] for key in (
                    "cnn_model", "cache", "reuse_directory", "reuse_target_id", "reuse_actual_protein_id", "reuse_receptor_sha256",
                    "maxit_executable", "maxit_sha256", "model_directory", "checkpoint", "model_card", "model_card_sha256",
                    "output_schema", "admitted_transformations") if key in profile}}
    return {**identity, "sha256": _hash(identity)}


def _admission(provider, profile, context):
    identity = model_asset_identity(provider, profile)
    receipts = profile.get("license_receipts")
    if not isinstance(receipts, list) or not receipts or len(receipts) > 100:
        raise ModelHeld("commercial_asset_receipts_required")
    scopes = set()
    for receipt in receipts:
        if not isinstance(receipt, dict):
            raise ValueError("Invalid license receipt")
        path = _asset(receipt.get("path"), receipt.get("sha256"))
        record = _json(path, 256 * 1024)
        if (not isinstance(record, dict) or record.get("license") not in _LICENSES
                or not isinstance(record.get("source_url"), str) or not record["source_url"].startswith("https://")
                or record.get("scope") not in {"code", "weights", "data"}):
            raise ModelHeld("commercial_asset_receipt_unresolved")
        scopes.add(record["scope"])
    if "code" not in scopes or provider in MODEL_PROVIDERS and "weights" not in scopes:
        raise ModelHeld("separate_code_and_weight_receipts_required")
    if context.get("qualification_mode") is not True:
        reference = profile.get("qualification_receipt")
        if not isinstance(reference, dict):
            raise ModelHeld("real_input_runtime_qualification_required")
        path = _asset(reference.get("path"), reference.get("sha256"))
        receipt = _json(path, 256 * 1024)
        if (not isinstance(receipt, dict) or receipt.get("status") != "qualified"
                or receipt.get("provider") != provider or receipt.get("asset_identity_sha256") != identity["sha256"]
                or any(not isinstance(receipt.get(k), str) or not _SHA.fullmatch(receipt[k]) for k in ("input_sha256", "output_sha256"))
                or isinstance(receipt.get("elapsed_seconds"), bool)
                or not isinstance(receipt.get("elapsed_seconds"), (int, float))
                or not math.isfinite(receipt["elapsed_seconds"]) or receipt["elapsed_seconds"] <= 0
                or isinstance(receipt.get("peak_memory_mib"), bool)
                or not isinstance(receipt.get("peak_memory_mib"), (int, float))
                or not math.isfinite(receipt["peak_memory_mib"]) or receipt["peak_memory_mib"] < 0):
            raise ModelHeld("runtime_qualification_does_not_bind_assets")
    return identity


def _inside(path, root, *, maximum=_MAX_INPUT):
    candidate = Path(path)
    if not candidate.is_absolute():
        candidate = root / candidate
    if candidate.is_symlink() or not candidate.is_file() or not candidate.resolve().is_relative_to(root):
        raise ValueError("Artifact must be a regular file inside the owned work directory")
    if candidate.stat().st_size > maximum:
        raise ValueError("Artifact exceeds size bound")
    return candidate.resolve()


def _input(row, name, root):
    path = _inside(row.get(name + "_path", ""), root)
    expected = row.get(name + "_sha256")
    if not isinstance(expected, str) or not _SHA.fullmatch(expected) or sha256_file(path) != expected:
        raise ArtifactIntegrityError("Input artifact checksum mismatch: " + name)
    return path


def _pinned_directory(value, profile, *, required=()):
    if not isinstance(value, str) or not Path(value).is_absolute():
        raise ModelHeld("absolute_offline_asset_directory_required")
    root = Path(value)
    if root.is_symlink() or not root.is_dir():
        raise ModelHeld("offline_asset_directory_missing")
    root = root.resolve()
    pinned = {Path(item["path"]).resolve() for item in profile.get("model_files", [])}
    count = 0
    for path in root.rglob("*"):
        count += 1
        if count > 100000 or path.is_symlink() or path.is_file() and path.resolve() not in pinned:
            raise ModelHeld("offline_directory_contains_unpinned_assets")
    if not count or any((root / relative).resolve() not in pinned for relative in required):
        raise ModelHeld("offline_asset_cache_incomplete")
    return root


def _molecule(smiles):
    from rdkit import Chem
    mol = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
    if mol is None or not 0 < mol.GetNumAtoms() <= 1000:
        raise ValueError("A bounded exact molecular structure is required")
    return mol


def _exact(mol):
    from rdkit import Chem
    from .identity import chemical_identity
    return chemical_identity(Chem.MolToCXSmiles(Chem.RemoveHs(mol)))["exact_structure_id"]


def _pose(path, expected):
    from rdkit import Chem
    mols = list(Chem.SDMolSupplier(str(path), removeHs=False))
    if len(mols) != 1 or mols[0] is None or mols[0].GetNumAtoms() > 1000:
        raise ValueError("Exactly one intact pose required")
    mol = mols[0]
    if (mol.GetNumConformers() != 1 or not mol.GetConformer().Is3D()
            or any(not math.isfinite(float(v)) for xyz in mol.GetConformer().GetPositions() for v in xyz)):
        raise ValueError("An actual finite 3D pose is required")
    if not expected or _exact(mol) != expected:
        raise ValueError("Pose chemical identity differs from the frozen prepared state")
    return mol


def _protein(row, context, root):
    reference_path = _input(row, "sequence_reference", root)
    record = _json(reference_path, 4 * 1024 * 1024)
    sequence = record.get("sequence", {}).get("value")
    identifier = "uniprot:" + str(record.get("primaryAccession", ""))
    if (not record.get("primaryAccession") or identifier != context.get("actual_protein_id")
            or not isinstance(sequence, str) or not re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWYX]{1,10000}", sequence)):
        raise ValueError("Exact frozen protein reference does not match target context")
    return sequence


def _write_binding(root, value):
    payload = {**value, "contract_version": CONTRACT_VERSION}
    (root / _MANIFEST).write_bytes(canonical_json({"payload": payload, "sha256": _hash(payload)}))


def prepare_model_command(provider, row, context, profile, workdir):
    """Build a fixed argv; missing assets/qualification yield an explicit hold.

    qualification_mode is an internal, owned canary operation, never an HTTP
    assertion. Normal execution requires a hash-bound local runtime receipt.
    """
    if provider not in PROVIDERS:
        raise ValueError("Unsupported commercially admitted provider")
    if not isinstance(row, dict) or not isinstance(context, dict):
        raise ValueError("Candidate and verified context must be objects")
    cid = row.get("candidate_id")
    if not isinstance(cid, str) or not cid or len(cid) > 200:
        raise ValueError("Bounded candidate ID required")
    if any(k in row or k in context for k in ("executable", "model_files", "argv", "shell", "command")):
        raise ValueError("Executable/model configuration belongs only to the administrator profile")
    for key in ("target_id", "actual_protein_id"):
        if key in row and row[key] != context.get(key):
            raise ValueError("Candidate target/protein differs from verified model context")
    root = Path(workdir).resolve()
    root.mkdir(parents=True, exist_ok=True)
    if (root / _MANIFEST).exists():
        raise ValueError("Use a fresh owned model work directory")
    try:
        identity = _admission(provider, profile, context)
    except ModelHeld as error:
        return {"status": "held", "provider": provider, "candidate_id": cid, "reason": str(error), "argv": []}
    threads = context.get("cpu_threads", 1)
    seed = context.get("seed", 42)
    if type(threads) is not int or not 1 <= threads <= 8 or type(seed) is not int or not 0 <= seed < 2**31:
        raise ValueError("Reserved CPU threads and bounded seed required")
    env = {"HF_HUB_OFFLINE": "1", "TRANSFORMERS_OFFLINE": "1", "HF_DATASETS_OFFLINE": "1",
           "OMP_NUM_THREADS": str(threads), "MKL_NUM_THREADS": str(threads), "OPENBLAS_NUM_THREADS": str(threads)}
    gpu_uuid = context.get("gpu_uuid")
    if gpu_uuid is not None:
        if not isinstance(gpu_uuid, str) or not re.fullmatch(r"(?:GPU|MIG)-[A-Za-z0-9/-]{1,100}", gpu_uuid):
            raise ValueError("Worker-reserved GPU UUID required")
        env["CUDA_VISIBLE_DEVICES"] = gpu_uuid
    elif provider in {"boltzina", "boltz2", "chai1", "unidock2"}:
        return {"status": "held", "provider": provider, "candidate_id": cid,
                "reason": "worker_gpu_reservation_required", "argv": []}
    else:
        env["CUDA_VISIBLE_DEVICES"] = ""
    argv = [identity["executable"]]
    binding = {"provider": provider, "candidate_id": cid, "target_id": context.get("target_id"),
               "asset_identity_sha256": identity["sha256"], "asset_identity": identity,
               "actual_protein_id": context.get("actual_protein_id"),
               "seed": seed, "inputs": {}, "scientific_validation": "NOT_ESTABLISHED", "contract_version": CONTRACT_VERSION}
    expected_outputs = []
    try:
        if provider in {"gnina", "boltzina", "unidock2"}:
            source = _input(row, "pose_sdf", root)
            mol = _pose(source, row.get("prepared_exact_structure_id"))
            binding["prepared_exact_structure_id"] = _exact(mol)
            binding["inputs"]["pose_sdf"] = {"path": str(source), "sha256": sha256_file(source)}
            receptor_key = "receptor_dms" if provider == "unidock2" else "receptor"
            receptor = _input(row, receptor_key, root)
            binding["inputs"][receptor_key] = {"path": str(receptor), "sha256": sha256_file(receptor)}
            if provider == "gnina":
                model = Path(profile.get("cnn_model", "")).resolve()
                if str(model) not in {m["path"] for m in identity["model_files"]}:
                    raise ModelHeld("external_gnina_checkpoint_not_pinned")
                argv += ["-r", str(receptor), "-l", str(source), "--score_only", "--cnn_model", str(model),
                         "--cpu", str(threads), "--seed", str(seed), "-o", str(root / "scored.sdf")]
                if not gpu_uuid:
                    argv += ["--no_gpu"]
                expected_outputs = ["scored.sdf"]
            elif provider == "boltzina":
                from rdkit import Chem
                cache = _pinned_directory(profile.get("cache"), profile, required=("boltz2_aff.ckpt",))
                _pinned_directory(str(cache / "mols"), profile)
                reuse = _pinned_directory(profile.get("reuse_directory"), profile, required=("processed/manifest.json",))
                if (profile.get("reuse_target_id") != context.get("target_id")
                        or profile.get("reuse_actual_protein_id") != context.get("actual_protein_id")
                        or profile.get("reuse_receptor_sha256") != row.get("receptor_sha256")):
                    raise ModelHeld("boltzina_reuse_target_receptor_provenance_required")
                reuse_manifest = _json(reuse / "processed/manifest.json", 4 * 1024 * 1024)
                records = reuse_manifest.get("records")
                if not isinstance(records, list) or len(records) != 1 or not isinstance(records[0], dict):
                    raise ModelHeld("boltzina_requires_one_frozen_reference_complex")
                maxit = _asset(profile.get("maxit_executable"), profile.get("maxit_sha256"), executable=True)
                if str(maxit) not in {m["path"] for m in identity["runtime_files"]}:
                    raise ModelHeld("boltzina_conversion_binary_not_pinned")
                ligand_tag = "candidate_" + hashlib.sha256(cid.encode()).hexdigest()[:16]
                mol.SetProp("_Name", ligand_tag)
                source = root / "boltzina-pose.sdf"
                with Chem.SDWriter(str(source)) as writer:
                    writer.write(mol)
                binding["ligand_tag"] = ligand_tag
                binding["reuse_directory"] = str(reuse)
                env.update(BOLTZ_CACHE=str(cache), BOLTZINA_MAXIT_PATH=str(maxit))
                argv += ["run", str(source), "--work-dir", str(reuse), "--receptor-pdb", str(receptor),
                         "--skip-docking", "--num-workers", "1", "--vina-cpu", "1", "--batch-size", "1",
                         "--seed", str(seed), "--keep-intermediate-files", "--output-dir", str(root / "boltzina")]
                expected_outputs = ["boltzina/boltzina_results.csv"]
            else:
                import yaml
                if receptor.suffix != ".dms":
                    raise ValueError("UniDock2 requires a separately prepared, pinned DMS receptor")
                center, size = context.get("center"), context.get("size")
                for key, values in (("center", center), ("size", size)):
                    if not isinstance(values, list) or len(values) != 3 or any(type(v) not in (int, float) or not math.isfinite(v) for v in values):
                        raise ValueError("Finite UniDock2 docking box required")
                    if key == "size" and any(not 0 < v <= 100 for v in values):
                        raise ValueError("Bounded UniDock2 box required")
                payload = {"Required": {"receptor": str(receptor), "ligand": str(source), "center": center},
                           "Settings": {"box_size": size, "search_mode": "free"},
                           "Advanced": {"seed": seed, "num_pose": 1},
                           "Preprocessing": {"output_sdf": str(root / "docked.sdf")}}
                (root / "unidock2.yaml").write_text(yaml.safe_dump(payload))
                argv += ["docking", "--config", str(root / "unidock2.yaml"), "--output_sdf", str(root / "docked.sdf")]
                expected_outputs = ["docked.sdf"]
        elif provider in {"boltz2", "chai1"}:
            from .identity import chemical_identity
            smiles = row.get("smiles")
            mol = _molecule(smiles)
            exact = chemical_identity(smiles)["exact_structure_id"]
            if exact != row.get("exact_structure_id"):
                raise ValueError("Ligand identity differs from frozen model input")
            sequence = _protein(row, context, root)
            binding.update(exact_structure_id=exact, smiles=smiles, protein_sequence=sequence,
                           ligand_elements=dict(Counter(a.GetSymbol() for a in mol.GetAtoms() if a.GetAtomicNum() != 1)))
            binding["inputs"]["sequence_reference"] = {"path": str(_input(row, "sequence_reference", root)),
                                                       "sha256": row["sequence_reference_sha256"]}
            if provider == "boltz2":
                import yaml
                cache = _pinned_directory(profile.get("cache"), profile,
                                          required=("mols.tar", "boltz2_conf.ckpt", "boltz2_aff.ckpt"))
                _pinned_directory(str(cache / "mols"), profile)
                msa = _input(row, "msa", root)
                lines = msa.read_text().splitlines()
                query = []
                started = False
                for line in lines:
                    if line.startswith(">"):
                        if started:
                            break
                        started = True
                    elif started:
                        query.append(re.sub(r"[a-z.-]", "", line.strip()))
                if "".join(query) != sequence:
                    raise ValueError("Local MSA query differs from frozen protein reference")
                payload = {"version": 1, "sequences": [{"protein": {"id": "A", "sequence": sequence, "msa": str(msa)}},
                                                        {"ligand": {"id": "B", "smiles": smiles}}],
                           "properties": [{"affinity": {"binder": "B"}}]}
                (root / "input.yaml").write_text(yaml.safe_dump(payload))
                binding["inputs"]["msa"] = {"path": str(msa), "sha256": row["msa_sha256"]}
                env["BOLTZ_CACHE"] = str(cache)
                argv += ["predict", str(root / "input.yaml"), "--out_dir", str(root), "--cache", str(cache),
                         "--checkpoint", str(cache / "boltz2_conf.ckpt"), "--affinity_checkpoint", str(cache / "boltz2_aff.ckpt"),
                         "--model", "boltz2", "--devices", "1", "--seed", str(seed), "--diffusion_samples", "1",
                         "--max_parallel_samples", "1", "--preprocessing_threads", "1", "--num_workers", "1"]
                expected_outputs = ["boltz_results_input/predictions/input/input_model_0.cif",
                                    "boltz_results_input/predictions/input/confidence_input_model_0.json",
                                    "boltz_results_input/predictions/input/affinity_input.json"]
            else:
                required = ["conformers_v1.apkl", "esm/traced_sdpa_esm2_t36_3B_UR50D_fp16.pt"]
                required += ["models_v2/" + name + ".pt" for name in ("feature_embedding", "bond_loss_input_proj", "token_embedder", "trunk", "diffusion_module", "confidence_head")]
                cache = _pinned_directory(profile.get("cache"), profile, required=required)
                (root / "input.fasta").write_text(">protein|name=A\n" + sequence + "\n>ligand|name=B\n" + smiles + "\n")
                env["CHAI_DOWNLOADS_DIR"] = str(cache)
                argv += ["fold", str(root / "input.fasta"), "--output-dir", str(root / "chai"), "--num-diffn-samples", "1",
                         "--num-trunk-samples", "1", "--seed", str(seed), "--fasta-names-as-cif-chains", "--low-memory"]
                expected_outputs = ["chai/pred.model_idx_0.cif", "chai/scores.model_idx_0.npz"]
        elif provider in {"openadmet", "chemeleon"}:
            from .identity import chemical_identity
            smiles = row.get("smiles")
            _molecule(smiles)
            if chemical_identity(smiles)["exact_structure_id"] != row.get("exact_structure_id"):
                raise ValueError("Property model input differs from frozen exact ligand")
            schema = profile.get("output_schema")
            if (not isinstance(schema, list) or not 1 <= len(schema) <= 100
                    or any(not isinstance(t, dict) or set(t) != {"column", "name", "unit"}
                           or any(not isinstance(t[k], str) or not t[k] or len(t[k]) > 200 for k in t) for t in schema)
                    or len({t["column"] for t in schema}) != len(schema)):
                raise ModelHeld("task_model_output_schema_required")
            card = _asset(profile.get("model_card"), profile.get("model_card_sha256"))
            if str(card) not in {m["path"] for m in identity["model_files"]}:
                raise ModelHeld("property_model_card_not_pinned")
            declared = _json(card, 256 * 1024)
            if declared.get("output_schema") != schema:
                raise ModelHeld("output_units_not_bound_to_model_card")
            with (root / "input.csv").open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=["candidate_id", "smiles"])
                writer.writeheader()
                writer.writerow({"candidate_id": cid, "smiles": smiles})
            binding.update(exact_structure_id=row["exact_structure_id"], smiles=smiles, output_schema=schema,
                           model_card_sha256=profile["model_card_sha256"])
            if provider == "openadmet":
                directory = _pinned_directory(profile.get("model_directory"), profile)
                argv += ["predict", "--input-path", str(root / "input.csv"), "--input-col", "smiles", "--model-dir", str(directory),
                         "--output-csv", str(root / "predictions.csv"), "--accelerator", "gpu" if gpu_uuid else "cpu"]
            else:
                if declared.get("checkpoint_kind") != "task_finetuned":
                    raise ModelHeld("descriptor_foundation_weights_are_not_task_predictions")
                checkpoint = Path(profile.get("checkpoint", "")).resolve()
                if str(checkpoint) not in {m["path"] for m in identity["model_files"]} or checkpoint.suffix not in {".pt", ".ckpt"}:
                    raise ModelHeld("task_finetuned_chemeleon_checkpoint_required")
                argv += ["predict", "--test-path", str(root / "input.csv"), "--smiles-columns", "smiles", "--model-paths", str(checkpoint),
                         "--preds-path", str(root / "predictions.csv"), "--num-workers", "0"]
            expected_outputs = ["predictions.csv"]
        else:
            if context.get("post_hit_opt_in") is not True:
                raise ModelHeld("openfe_initial_count_zero_post_hit_opt_in_required")
            transformation = _input(row, "transformation", root)
            admitted = profile.get("admitted_transformations", [])
            if row["transformation_sha256"] not in admitted:
                raise ModelHeld("openfe_transformation_not_admitted_by_operator")
            binding["inputs"]["transformation"] = {"path": str(transformation), "sha256": row["transformation_sha256"]}
            argv += ["quickrun", str(transformation), "-d", str(root / "openfe-work"), "-o", str(root / "openfe-result.json")]
            expected_outputs = ["openfe-result.json"]
    except ModelHeld as error:
        return {"status": "held", "provider": provider, "candidate_id": cid, "reason": str(error), "argv": []}
    for filename in ("input.yaml", "input.fasta", "input.csv", "boltzina-pose.sdf", "unidock2.yaml"):
        path = root / filename
        if path.exists():
            path = _inside(path, root)
            binding["inputs"]["generated_" + filename] = {"path": str(path), "sha256": sha256_file(path)}
    _write_binding(root, binding)
    return {"status": "ready", "provider": provider, "candidate_id": cid, "argv": argv, "env": env,
            "expected_outputs": expected_outputs, "asset_identity": identity, "input_binding_sha256": _hash(binding),
            "offline_assets_required": True, "ranking_weight": 0}


def _binding(provider, root, row):
    envelope = _json(_inside(_MANIFEST, root), 256 * 1024)
    record = envelope.get("payload")
    if not isinstance(record, dict) or envelope.get("sha256") != _hash(record):
        raise ArtifactIntegrityError("Model input binding checksum mismatch")
    if record.get("provider") != provider or record.get("candidate_id") != row.get("candidate_id"):
        raise ValueError("Output belongs to a different provider or candidate")
    for key in ("target_id", "actual_protein_id"):
        if key in row and row[key] != record.get(key):
            raise ValueError("Output target/protein differs from frozen candidate context")
    for key in ("exact_structure_id", "prepared_exact_structure_id"):
        if record.get(key) and record[key] != row.get(key):
            raise ValueError("Output input identity differs from frozen candidate")
    for name, entry in record.get("inputs", {}).items():
        if not name.startswith("generated_") and row.get(name + "_sha256") != entry["sha256"]:
            raise ValueError("Output source artifact differs from frozen candidate input")
        path = _inside(entry["path"], root)
        if sha256_file(path) != entry["sha256"]:
            raise ArtifactIntegrityError("Frozen input changed during model execution")
    assets = record.get("asset_identity", {})
    if assets.get("adapter_source_sha256") != sha256_file(Path(__file__)):
        raise ArtifactIntegrityError("Model adapter source changed during execution")
    _asset(assets.get("executable"), assets.get("executable_sha256"), executable=True)
    for entry in [*assets.get("model_files", []), *assets.get("runtime_files", [])]:
        _asset(entry.get("path"), entry.get("sha256"))
    return record


def _metric(name, value, kind, unit, *, bounded=False):
    if isinstance(value, bool):
        raise ValueError("Boolean model score is not a numeric measurement")
    try:
        number = float(value)
    except (TypeError, ValueError, OverflowError) as error:
        raise ValueError("Nonnumeric model output: " + name) from error
    if not math.isfinite(number) or bounded and not 0 <= number <= 1:
        raise ValueError("Nonfinite or out-of-range model output: " + name)
    return {"name": name, "value": number, "kind": kind, "unit": unit, "calibration_status": "NOT_ESTABLISHED"}


def _affinity_metrics(values):
    result = []
    for name in ("affinity_pred_value", "affinity_pred_value1", "affinity_pred_value2"):
        if name in values:
            result.append(_metric(name, values[name], "affinity_model_output", "log10(IC50_uM)"))
    for name in ("affinity_probability_binary", "affinity_probability_binary1", "affinity_probability_binary2"):
        if name in values:
            result.append(_metric(name, values[name], "binary_model_output", "unitless", bounded=True))
    if not any(item["name"] == "affinity_pred_value" for item in result):
        raise ValueError("Expected affinity model output is missing")
    return result


def _csv_one(path):
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError("Unique CSV output columns required")
        first = next(reader, None)
        if first is None or next(reader, None) is not None or None in first:
            raise ValueError("Exactly one complete candidate output row required")
        return first


def _complex_identity(path, binding):
    from Bio.PDB import MMCIFParser
    from Bio.SeqUtils import seq1
    structure = MMCIFParser(QUIET=True, auth_chains=False).get_structure("model", str(path))
    models = list(structure)
    if len(models) != 1:
        raise ValueError("Exactly one predicted complex model required")
    chains = {chain.id: chain for chain in models[0]}
    if set(chains) != {"A", "B"}:
        raise ValueError("Predicted complex requires exact protein A and ligand B chains")
    atoms = list(structure.get_atoms())
    if not 1 <= len(atoms) <= 100000 or any(not math.isfinite(float(v)) for atom in atoms for v in atom.coord):
        raise ValueError("Predicted complex has invalid coordinates or atom count")
    sequence = "".join(seq1(residue.resname, custom_map={"UNK": "X"}) for residue in chains["A"])
    if sequence != binding["protein_sequence"]:
        raise ValueError("Predicted protein sequence differs from frozen target")
    elements = Counter(atom.element.title() for atom in chains["B"].get_atoms() if atom.element.upper() not in {"H", "D"})
    if dict(elements) != binding["ligand_elements"]:
        raise ValueError("Predicted ligand elemental identity differs from frozen ligand")
    return {"atom_count": len(atoms), "protein_sequence_verified": True, "ligand_elements_verified": True,
            "pose_chemical_graph_validation": "NOT_ESTABLISHED"}


def parse_model_outputs(provider, workdir, row):
    """Require actual artifacts and input identity; exit code alone proves nothing."""
    if provider not in PROVIDERS:
        raise ValueError("Unsupported provider")
    root = Path(workdir).resolve()
    binding = _binding(provider, root, row)
    files = []
    metrics = []
    checks = {}
    if provider in {"gnina", "unidock2"}:
        from rdkit import Chem
        path = _inside("scored.sdf" if provider == "gnina" else "docked.sdf", root)
        mol = _pose(path, binding["prepared_exact_structure_id"])
        if provider == "gnina":
            original = Chem.RemoveHs(_pose(Path(binding["inputs"]["pose_sdf"]["path"]), binding["prepared_exact_structure_id"]))
            scored = Chem.RemoveHs(mol)
            matches = scored.GetSubstructMatches(original, useChirality=True, uniquify=False, maxMatches=1000)
            before, after = original.GetConformer(), scored.GetConformer()
            if not any(all((before.GetAtomPosition(i) - after.GetAtomPosition(j)).Length() <= .002 for i, j in enumerate(match)) for match in matches):
                raise ValueError("GNINA score-only output moved the admitted pose")
            for name, kind, unit, bounded in (("CNNscore", "pose_model_score", "unitless", True),
                                             ("CNNaffinity", "affinity_model_score", "pK_model_output", False)):
                if not mol.HasProp(name):
                    raise ValueError("GNINA output lacks " + name)
                metrics.append(_metric(name, mol.GetProp(name), kind, unit, bounded=bounded))
            if mol.HasProp("minimizedAffinity"):
                metrics.append(_metric("minimizedAffinity", mol.GetProp("minimizedAffinity"), "empirical_energy", "kcal/mol"))
        else:
            # UniDock2 versions use energy/score property names; the admitted
            # output must contain one recognized empirical score, never a CNN.
            for name in ("vina_binding_free_energy", "docking_score", "score", "energy", "Uni-Dock2_score"):
                if mol.HasProp(name):
                    metrics.append(_metric(name, mol.GetProp(name), "empirical_energy", "kcal/mol"))
                    break
            if not metrics:
                raise ValueError("UniDock2 output lacks a recognized empirical docking score")
        checks = {"exact_structure_id": _exact(mol), "pose_identity_verified": True}
        files = [path]
    elif provider == "boltzina":
        path = _inside("boltzina/boltzina_results.csv", root)
        values = _csv_one(path)
        ligand = values.get("ligand_name", "")
        if (Path(ligand).name != binding["ligand_tag"] + ".pdb" or values.get("ligand_idx") != "0"
                or values.get("docking_rank") not in {"0", "1"}):
            raise ValueError("Boltzina output does not belong to the single admitted pose")
        metrics = _affinity_metrics(values)
        # The input SDF supplies exact chemistry and frozen coordinates. The
        # affinity-only branch's CSV carries no independent chemical graph.
        checks = {"prepared_exact_structure_id": binding["prepared_exact_structure_id"],
                  "existing_pose_reused": True, "output_chemical_graph_validation": "NOT_ESTABLISHED"}
        files = [path]
    elif provider in {"boltz2", "chai1"}:
        if provider == "boltz2":
            directory = "boltz_results_input/predictions/input/"
            structure_path = _inside(directory + "input_model_0.cif", root)
            confidence = _inside(directory + "confidence_input_model_0.json", root)
            affinity = _inside(directory + "affinity_input.json", root)
            values = _json(confidence, 4 * 1024 * 1024)
            for name in ("confidence_score", "ptm", "iptm", "ligand_iptm", "complex_plddt"):
                if name in values:
                    metrics.append(_metric(name, values[name], "structure_confidence", "unitless", bounded=True))
            if not metrics:
                raise ValueError("Boltz output lacks confidence metrics")
            metrics.extend(_affinity_metrics(_json(affinity, 4 * 1024 * 1024)))
            files = [structure_path, confidence, affinity]
        else:
            import numpy as np
            structure_path = _inside("chai/pred.model_idx_0.cif", root)
            scores = _inside("chai/scores.model_idx_0.npz", root, maximum=8 * 1024 * 1024)
            with zipfile.ZipFile(scores) as archive:
                info = archive.infolist()
                if len(info) > 100 or sum(item.file_size for item in info) > 8 * 1024 * 1024:
                    raise ValueError("Chai score archive exceeds decompressed size bounds")
            with np.load(scores, allow_pickle=False) as values:
                for name in ("aggregate_score", "ptm", "iptm"):
                    if name in values:
                        array = values[name]
                        if array.size != 1:
                            raise ValueError("Expected scalar Chai confidence metric")
                        metrics.append(_metric(name, array.item(), "structure_confidence", "unitless", bounded=True))
            if not metrics:
                raise ValueError("Chai output lacks actual confidence metrics")
            files = [structure_path, scores]
        checks = _complex_identity(structure_path, binding)
    elif provider in {"openadmet", "chemeleon"}:
        from .identity import chemical_identity
        path = _inside("predictions.csv", root)
        values = _csv_one(path)
        if (values.get("candidate_id") != binding["candidate_id"]
                or chemical_identity(values.get("smiles"))["exact_structure_id"] != binding["exact_structure_id"]):
            raise ValueError("Property output candidate or exact chemical identity differs")
        for task in binding["output_schema"]:
            if task["column"] not in values:
                raise ValueError("Property output lacks declared endpoint: " + task["column"])
            metrics.append(_metric(task["name"], values[task["column"]], "property_model_output", task["unit"]))
        checks = {"exact_structure_id": binding["exact_structure_id"], "model_card_sha256": binding["model_card_sha256"]}
        files = [path]
    else:
        path = _inside("openfe-result.json", root)
        values = _json(path)
        if not isinstance(values, dict) or values.get("estimate") is None or values.get("uncertainty") is None or not isinstance(values.get("protocol_result"), dict):
            raise ValueError("OpenFE simulation result is incomplete")
        # Gufe's JSON decoder reconstitutes units. Do not flatten or guess its
        # versioned quantity representation, or equate one leg with binding dG.
        try:
            from gufe.tokenization import JSON_HANDLER
        except ImportError:
            return {"status": "held", "provider": provider, "candidate_id": binding["candidate_id"],
                    "reason": "openfe_quantity_decoder_unavailable", "outputs": {"metrics": []},
                    "artifacts": [{"name": path.name, "sha256": sha256_file(path), "size_bytes": path.stat().st_size}],
                    "ranking_weight": 0, "scientific_validation": "NOT_ESTABLISHED"}
        quantities = json.loads(path.read_text(), cls=JSON_HANDLER.decoder)
        for name in ("estimate", "uncertainty"):
            quantity = quantities[name].to("kilocalorie_per_mole")
            value = quantity.magnitude
            if name == "uncertainty" and value < 0:
                raise ValueError("OpenFE uncertainty cannot be negative")
            metrics.append(_metric(name, value, "transformation_free_energy", "kcal/mol"))
        checks = {"transformation_sha256": binding["inputs"]["transformation"]["sha256"], "binding_delta_g_established": False}
        files = [path]
    artifacts = [{"name": str(path.relative_to(root)), "sha256": sha256_file(path), "size_bytes": path.stat().st_size} for path in files]
    return {"status": "complete", "provider": provider, "candidate_id": binding["candidate_id"],
            "asset_identity_sha256": binding["asset_identity_sha256"], "input_binding_sha256": _hash(binding),
            "outputs": {"metrics": metrics, "identity_checks": checks}, "artifacts": artifacts,
            "ranking_weight": 0, "scientific_validation": "NOT_ESTABLISHED"}
