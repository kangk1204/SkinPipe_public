"""Optional, separately budgeted scorers; outputs never alter priority_score_v2.

GNINA rescores an admitted pose. Boltz accepts a frozen local input package.
Model confidence, predicted affinity and binary model outputs stay separate.
No provider is enabled or downloaded by this module. Site worker configuration,
not an HTTP request, chooses executable paths and exact model files.
"""
from __future__ import annotations

from copy import deepcopy
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import time
import uuid

from rdkit import Chem

from .artifacts import ArtifactIntegrityError, canonical_json, sha256_file
from .identity import chemical_identity
from .resources import cpu_limit
from .ranking import finite
from ..processes import _stop_group

PROVIDERS = {'gnina': 'pose_rescore', 'boltz2': 'complex_prediction'}
OUTPUTS = {
    'CNNscore': ('pose_model_score', 'unitless'),
    'CNNaffinity': ('affinity_model_score', 'pK_model_output'),
    'minimizedAffinity': ('empirical_energy', 'kcal/mol'),
    'confidence_score': ('structure_confidence', 'unitless'),
    'iptm': ('interface_confidence', 'unitless'),
    'ligand_iptm': ('ligand_interface_confidence', 'unitless'),
    'affinity_pred_value': ('affinity_model_output', 'log10(IC50_uM)'),
    'affinity_probability_binary': ('binary_model_output', 'unitless'),
}


class OptionalCancelled(RuntimeError):
    pass


def _bounded_output(path, root, maximum=64 * 1024 * 1024):
    path, root = Path(path), Path(root)
    if (path.is_symlink() or not path.is_file() or not path.resolve().is_relative_to(root.resolve())
            or path.stat().st_size > maximum):
        raise ValueError('Optional output artifact outside file or size bounds')
    return path


def _stop_process(process):
    # All optional children start in an owned session. A completed direct
    # child can still leave live descendants; cleanup covers that case too.
    _stop_group(process)


def _complex_structure(path):
    from Bio.PDB import MMCIFParser, PDBParser

    parser = MMCIFParser(QUIET=True) if path.suffix == '.cif' else PDBParser(QUIET=True)
    structure = parser.get_structure('prediction', str(path))
    atoms = list(structure.get_atoms())
    chains = {chain.id for chain in structure.get_chains()}
    if (not atoms or len(atoms) > 100000 or not {'A', 'B'} <= chains
            or any(not math.isfinite(float(value)) for atom in atoms for value in atom.coord)):
        raise ValueError('Boltz structure must contain finite protein A and ligand B coordinates')
    return {'file': path.name, 'atom_count': len(atoms), 'chains': sorted(chains)}


def finite_metrics(values):
    if not isinstance(values, dict):
        raise ValueError('Model metrics must be an object')
    metrics = []
    for name, (kind, unit) in OUTPUTS.items():
        value = values.get(name)
        if value is None:
            continue
        if finite(value) is None:
            raise ValueError('Nonfinite or nonnumeric model output: ' + name)
        if unit == 'unitless' and not 0 <= value <= 1:
            raise ValueError('Out-of-range model output: ' + name)
        metrics.append({'name': name, 'kind': kind, 'unit': unit, 'value': value,
                        'calibration_status': 'NOT_ESTABLISHED'})
    return metrics


def parse_outputs(provider, directory):
    """Require actual model artifacts; subprocess success alone is insufficient."""
    root = Path(directory)
    if provider == 'gnina':
        path = root / 'scored.sdf'
        _bounded_output(path, root)
        molecules = list(Chem.SDMolSupplier(str(path), removeHs=False))
        if len(molecules) != 1 or molecules[0] is None:
            raise ValueError('Expected exactly one intact scored pose')
        mol = molecules[0]
        values = {key: float(mol.GetProp(key)) for key in OUTPUTS if mol.HasProp(key)}
        if 'CNNscore' not in values or 'CNNaffinity' not in values:
            raise ValueError('GNINA pose lacks expected CNN outputs')
        return {'metrics': finite_metrics(values), 'files': [path],
                'output_identity': chemical_identity(Chem.MolToSmiles(Chem.RemoveHs(mol)))['exact_structure_id']}
    if provider == 'boltz2':
        files = sorted(root.glob('boltz_results_input/predictions/input/*'))
        if len(files) > 100:
            raise ValueError('Optional model generated too many files')
        for path in files:
            _bounded_output(path, root, 8 * 1024 * 1024 if path.suffix == '.json' else 64 * 1024 * 1024)
        confidences = [p for p in files if p.name.startswith('confidence_') and p.suffix == '.json']
        structures = [p for p in files if p.suffix in {'.cif', '.pdb'}]
        if not confidences or not structures:
            raise ValueError('Boltz requires confidence and actual structure artifacts')
        structure_checks = [_complex_structure(path) for path in structures]
        measurements = []
        for path in files:
            if path.suffix == '.json' and path.name.startswith(('confidence_', 'affinity_')):
                values = json.loads(path.read_text())
                measurements.append({'file': path.name, 'metrics': finite_metrics(values)})
        if not any(record['metrics'] for record in measurements):
            raise ValueError('Boltz output contains no recognized confidence or affinity outputs')
        return {'measurements': measurements, 'files': files,
                'structure_checks': structure_checks,
                'output_identity': None, 'pose_chemical_validation': 'NOT_ESTABLISHED'}
    raise ValueError('Unsupported optional scorer')


class OptionalScorers:
    def __init__(self, root, profiles):
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)
        self.profiles = deepcopy(profiles)

    def status(self):
        return {name: {'configured': name in self.profiles, 'automatic_candidates': 0,
                       'status': 'opt_in_and_runtime_verification_required' if name in self.profiles else 'held',
                       'execution_via_research_api': False,
                       'maximum_opt_in_candidates': 12, 'ranking_weight': 0,
                       'purpose': purpose} for name, purpose in PROVIDERS.items()}

    def _profile(self, provider):
        if provider not in PROVIDERS or provider not in self.profiles:
            raise ValueError('Optional model not configured; ordinary candidate handoff remains available')
        profile = deepcopy(self.profiles[provider])
        if not all(profile.get(k) for k in ('version', 'preprocessing_version', 'executable', 'executable_sha256', 'model_files')):
            raise ValueError('Pinned executable, runtime and model manifest required')
        binary = Path(profile['executable']).resolve()
        if not Path(profile['executable']).is_absolute() or not binary.is_file() or not os.access(binary, os.X_OK):
            raise ValueError('Optional executable must be an absolute executable file')
        if sha256_file(binary) != profile['executable_sha256']:
            raise ArtifactIntegrityError('Optional scorer executable changed')
        if not isinstance(profile['model_files'], list) or len(profile['model_files']) > 100000:
            raise ValueError('Bounded pinned model file list required')
        model_paths = set()
        for item in profile['model_files']:
            if not isinstance(item, dict) or not Path(item.get('path', '')).is_absolute():
                raise ValueError('Model file paths must be absolute administrator configuration')
            if sha256_file(Path(item['path'])) != item['sha256']:
                raise ArtifactIntegrityError('Optional scorer model changed')
            model_paths.add(Path(item['path']).resolve())
        profile['executable'] = str(binary)
        if provider == 'gnina' and profile.get('cnn_model') and Path(profile['cnn_model']).resolve() not in model_paths:
            raise ValueError('GNINA external CNN must be pinned in model_files')
        if provider == 'boltz2':
            cache = Path(profile.get('cache', '')).resolve()
            # The official CLI calls download_boltz2 even when --checkpoint is
            # supplied. All default assets must already exist and be pinned.
            required = [cache / name for name in ('mols.tar', 'boltz2_conf.ckpt', 'boltz2_aff.ckpt')]
            molecules = cache / 'mols'
            if (not Path(profile.get('cache', '')).is_absolute() or not molecules.is_dir()
                    or molecules.is_symlink() or any(path not in model_paths for path in required)):
                raise ValueError('Boltz offline cache is incomplete or unpinned; execution held')
            molecule_files = list(molecules.rglob('*'))
            if not molecule_files or any(path.is_symlink() or path.is_file() and path.resolve() not in model_paths
                                         for path in molecule_files):
                raise ValueError('Boltz molecule cache must be complete and content-pinned')
            profile['cache'] = str(cache)
        return profile

    def execute(self, request, *, fetch_artifact, publish_artifact, cancel_event, progress=lambda _: None):
        request = deepcopy(request)
        request_sha256 = hashlib.sha256(canonical_json(request)).hexdigest()
        provider = request.get('provider')
        profile = self._profile(provider)
        budget = request.get('budget_seconds')
        candidates = request.get('candidates')
        if request.get('opt_in') is not True or finite(budget) is None or not 0 < budget <= 28800:
            raise ValueError('Optional scoring requires explicit opt-in and a separate positive budget <=8h')
        if not isinstance(candidates, list) or not 1 <= len(candidates) <= 12:
            raise ValueError('Optional scoring permits 1..12 explicitly selected candidates')
        if any(not isinstance(row, dict) or not isinstance(row.get('candidate_id'), str)
               or not row['candidate_id'] or len(row['candidate_id']) > 200 for row in candidates):
            raise ValueError('Explicit candidate IDs are required')
        if len({row['candidate_id'] for row in candidates}) != len(candidates):
            raise ValueError('Optional scoring candidate IDs must be unique')
        started = time.monotonic()
        result = {'provider': provider, 'profile_sha256': hashlib.sha256(canonical_json(profile)).hexdigest(),
                  'request_sha256': request_sha256, 'profile': deepcopy(profile),
                  'ranking_weight': 0, 'scientific_validation': 'NOT_ESTABLISHED', 'candidates': []}
        for row in candidates:
            receipt = {'candidate_id': row['candidate_id'], 'inputs': deepcopy(row), 'status': 'held',
                       'request_sha256': request_sha256, 'profile_sha256': result['profile_sha256']}
            remaining = budget - (time.monotonic() - started)
            if remaining <= 0 or cancel_event.is_set():
                receipt['status'] = 'budget_deferred' if remaining <= 0 else 'cancelled'
                result['candidates'].append(receipt)
                progress({'partial_result': deepcopy(result), 'compute_seconds': time.monotonic()-started})
                continue
            work = self.root / uuid.uuid4().hex
            work.mkdir()
            try:
                def materialize(key, filename):
                    digest = row[key]
                    if not isinstance(digest, str) or not re.fullmatch(r'[0-9a-f]{64}', digest):
                        raise ValueError('Optional scorer input requires a SHA256 identifier')
                    data = fetch_artifact(digest)
                    if len(data) > 64 * 1024 * 1024 or hashlib.sha256(data).hexdigest() != digest:
                        raise ArtifactIntegrityError('Optional scorer input hash mismatch')
                    path = work / filename
                    path.write_bytes(data)
                    return path
                if provider == 'gnina':
                    receptor = materialize('receptor_sha256', 'receptor.pdb')
                    ligand = materialize('pose_sdf_sha256', 'ligand.sdf')
                    mols = list(Chem.SDMolSupplier(str(ligand), removeHs=False))
                    if (len(mols) != 1 or mols[0] is None or mols[0].GetNumConformers() != 1
                            or not mols[0].GetConformer().Is3D()
                            or any(not math.isfinite(value) for point in mols[0].GetConformer().GetPositions() for value in point)):
                        raise ValueError('One admitted 3D pose is required')
                    identity = chemical_identity(Chem.MolToSmiles(Chem.RemoveHs(mols[0])))['exact_structure_id']
                    if identity is None or identity != row.get('prepared_exact_structure_id'):
                        raise ValueError('Prepared pose identity differs from frozen candidate binding')
                    cnn = profile.get('cnn', 'default1.0')
                    if not re.fullmatch(r'[A-Za-z0-9_.-]+', cnn):
                        raise ValueError('Invalid configured CNN name')
                    argv = [profile['executable'], '-r', str(receptor), '-l', str(ligand), '--score_only',
                            '--cnn', cnn, '--cpu', str(cpu_limit(os.cpu_count() or 1)), '-o', str(work/'scored.sdf')]
                    if profile.get('cnn_model'):
                        argv += ['--cnn_model', profile['cnn_model']]
                else:
                    # The input is a previously reviewed, hash-frozen package.
                    # It cannot reference files outside this private job folder.
                    import yaml
                    source = json.loads(materialize('input_manifest_sha256', 'manifest.json').read_text())
                    sequence, smiles = source.get('protein_sequence'), source.get('ligand_smiles')
                    if not isinstance(sequence, str) or not re.fullmatch(r'[ACDEFGHIKLMNPQRSTVWYX]{1,10000}', sequence):
                        raise ValueError('A bounded verified protein sequence is required')
                    ligand_identity = chemical_identity(smiles)['exact_structure_id']
                    if (source.get('target_id') != request.get('target_id') or ligand_identity is None
                            or ligand_identity != row.get('exact_structure_id')):
                        raise ValueError('Boltz input target or exact ligand binding differs')
                    reference = json.loads(materialize('sequence_reference_sha256', 'uniprot.json').read_text())
                    protein_id = 'uniprot:' + str(reference.get('primaryAccession', ''))
                    if (not reference.get('primaryAccession') or request.get('actual_protein_id') != protein_id
                            or source.get('actual_protein_id') != protein_id
                            or reference.get('sequence', {}).get('value') != sequence):
                        raise ValueError('Boltz sequence does not match the frozen exact protein reference')
                    msa = materialize('msa_sha256', 'protein.a3m')
                    query, in_query = [], False
                    for line in msa.read_text().splitlines():
                        if line.startswith('>'):
                            if in_query:
                                break
                            in_query = True
                        elif in_query:
                            query.append(re.sub(r'[a-z.-]', '', line.strip()))
                    if ''.join(query) != sequence:
                        raise ValueError('MSA query differs from the frozen protein sequence')
                    payload = {'version': 1, 'sequences': [{'protein': {'id': 'A', 'sequence': sequence, 'msa': str(msa)}},
                                                        {'ligand': {'id': 'B', 'smiles': smiles}}]}
                    if profile.get('affinity_checkpoint'):
                        payload['properties'] = [{'affinity': {'binder': 'B'}}]
                    (work/'input.yaml').write_text(yaml.safe_dump(payload))
                    model_paths = {str(Path(m['path']).resolve()) for m in profile['model_files']}
                    if str(Path(profile['checkpoint']).resolve()) not in model_paths:
                        raise ValueError('Boltz checkpoint must be pinned in model_files')
                    argv = [profile['executable'], 'predict', str(work/'input.yaml'), '--out_dir', str(work),
                            '--checkpoint', profile['checkpoint'], '--cache', profile['cache'], '--devices', '1',
                            '--model', 'boltz2', '--seed', '42', '--preprocessing_threads', '1',
                            '--diffusion_samples', '1', '--max_parallel_samples', '1', '--num_workers', '1']
                    if profile.get('affinity_checkpoint'):
                        if str(Path(profile['affinity_checkpoint']).resolve()) not in model_paths:
                            raise ValueError('Boltz affinity checkpoint must be pinned')
                        argv += ['--affinity_checkpoint', profile['affinity_checkpoint']]
                with (work/'process.log').open('wb') as log:
                    receipt['execution'] = {'argv': list(argv), 'python_version': sys.version,
                                            'configured_gpu_uuid': profile.get('gpu_uuid'),
                                            'model_versions': {'version': profile['version'],
                                                               'preprocessing': profile['preprocessing_version']}}
                    environment = dict(os.environ)
                    if profile.get('gpu_uuid'):
                        environment['CUDA_VISIBLE_DEVICES'] = str(profile['gpu_uuid'])
                    process = subprocess.Popen(argv, cwd=work, stdout=log, stderr=subprocess.STDOUT,
                                               start_new_session=True, env=environment)
                    try:
                        while process.poll() is None:
                            if cancel_event.wait(.1):
                                _stop_process(process)
                                raise OptionalCancelled('Optional model cancelled')
                            if time.monotonic()-started >= budget:
                                _stop_process(process)
                                raise TimeoutError('Optional model budget exhausted')
                        if process.returncode:
                            raise RuntimeError('Optional scorer process failed')
                    finally:
                        _stop_process(process)
                if hashlib.sha256(canonical_json(self._profile(provider))).hexdigest() != result['profile_sha256']:
                    raise ArtifactIntegrityError('Optional runtime changed during execution')
                parsed = parse_outputs(provider, work)
                if provider == 'gnina' and parsed['output_identity'] != row['prepared_exact_structure_id']:
                    raise ValueError('GNINA output changed ligand identity')
                if provider == 'gnina':
                    original = Chem.RemoveHs(mols[0])
                    scored = Chem.RemoveHs(next(iter(Chem.SDMolSupplier(str(work/'scored.sdf'), removeHs=False))))
                    matches = scored.GetSubstructMatches(original, uniquify=False, maxMatches=1000, useChirality=True)
                    before, after = original.GetConformer(), scored.GetConformer()
                    if not any(all((before.GetAtomPosition(i)-after.GetAtomPosition(j)).Length() <= .002
                                   for i, j in enumerate(match)) for match in matches):
                        raise ValueError('GNINA score-only output changed the admitted pose coordinates')
                artifacts = []
                for path in [*parsed.pop('files'), work/'process.log']:
                    _bounded_output(path, work)
                    raw = path.read_bytes()
                    published = publish_artifact(raw)
                    expected = hashlib.sha256(raw).hexdigest()
                    if (not isinstance(published, dict) or published.get('sha256') != expected
                            or published.get('size_bytes', len(raw)) != len(raw)):
                        raise ArtifactIntegrityError('Optional output publication did not confirm exact bytes')
                    artifacts.append({'name': path.name, 'sha256': expected, 'size_bytes': len(raw)})
                receipt.update(status='complete', outputs=parsed, artifacts=artifacts)
            except Exception as error:
                receipt.update(status='cancelled' if isinstance(error, OptionalCancelled) else
                               'budget_deferred' if isinstance(error, TimeoutError) else 'failed',
                               error_type=type(error).__name__, error=str(error))
            result['candidates'].append(receipt)
            progress({'partial_result': deepcopy(result), 'compute_seconds': time.monotonic()-started})
        result['compute_seconds'] = time.monotonic()-started
        return result
