"""Incremental fixed-pocket docking with exact-protocol reuse and checkpoints."""
from contextlib import contextmanager
from copy import deepcopy
import json
from pathlib import Path
import sys
import threading
import time

from .artifacts import atomic_bytes, canonical_json, sha256_file
from .engines import EngineRunner, ComputeBudgetExceeded
from .identity import chemical_identity, digest
from .native_prepare import PreparationUnavailable
from .routine_screening import execution_sources, next_compute_batch, select_representatives, validate_policy


@contextmanager
def _owned_cache_lock(path, stop, *, on_wait=None):
    # The analysis host is Ubuntu (macOS also supports this lock). OS ownership
    # releases it on cancellation or process exit; no stale PID-file takeover.
    import fcntl
    with path.open('a+b') as stream:
        while True:
            stop()
            try:
                fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
                break
            except BlockingIOError:
                if on_wait is not None:
                    on_wait()
                    on_wait = None
                time.sleep(0.05)
        try:
            yield
        finally:
            fcntl.flock(stream, fcntl.LOCK_UN)


def pose_quality(payload):
    from .algorithm_chemistry import execute_chemistry
    result = execute_chemistry('posebusters', payload['rows'], {
        'receptor_path': payload['receptor_path'], 'receptor_sha256': payload['receptor_sha256'],
        'candidate_limit': len(payload['rows']), 'budget_seconds': payload['budget_seconds'],
    }, workdir=Path(payload['workdir']))
    return {'status': 'complete', 'quality': result}


def _cache_key(frozen, engine_request, row, source_sha, execution_sources):
    return digest({'docking_policy': frozen, 'execution_sources': execution_sources,
                   'execution': engine_request, 'observed_protein': source_sha,
                   'target_id': row.get('target_id'),
                   'exact_structure': chemical_identity(row['smiles'])['exact_structure_id']})


def evaluate_routine_complex(rows, descriptor, selected, service, config, request, *, stop, progress):
    from . import research_docking as native
    from ..dock.ligprep import read_preparation_receipt
    rows = deepcopy(rows)
    policy = validate_policy(request['screening_policy'])
    current_sources = execution_sources()
    frozen = request.get('research_docking_policy') or native.freeze_docking_policy(config)
    receipt = {'policy': frozen, 'screening_policy': policy, 'execution_sources': current_sources, 'status': 'held',
        'receptor_valid': False, 'prepared_count': 0, 'selected_count': 0,
        'completed_count': 0, 'qualified_count': 0, 'reused_count': 0,
        'new_docking_count': 0, 'candidate_results': [], 'redocking': {'status': 'not_run'},
        'pocket_box': deepcopy(descriptor.get('pocket_box')),
        'protein_structure_sha256': descriptor.get('protein_structure_sha256')}
    if descriptor.get('docking_hold_reasons'):
        return rows, {**receipt, 'hold_reasons': deepcopy(descriptor['docking_hold_reasons'])}
    if frozen != native.freeze_docking_policy(config):
        return rows, {**receipt, 'hold_reasons': ['configured_docking_cohort_changed_since_submission']}
    if frozen['status'] != 'configured':
        return rows, {**receipt, 'hold_reasons': ['docking_engine_not_configured']}
    work = service.root / 'research-docking' / request['discovery_run_id'] / 'routine'
    work.mkdir(parents=True, exist_ok=True)
    cache = service.root / 'routine-docking-cache'
    cache.mkdir(parents=True, exist_ok=True)
    signature = digest({'rows': [(r['candidate_id'], r.get('smiles')) for r in rows],
        'descriptor': descriptor, 'selected': selected, 'policy': policy, 'docking': frozen,
        'execution_sources': current_sources})
    checkpoint = work / 'checkpoint.json'
    state = {'signature': signature, 'rows': {}, 'engine_rows': {}, 'prepared': {},
             'receipt': {}, 'compute_seconds': 0}
    if checkpoint.exists():
        state = json.loads(checkpoint.read_bytes())
        if state.get('signature') != signature:
            return rows, {**receipt, 'hold_reasons': ['routine_checkpoint_conditions_changed']}
    checkpoint_lock = threading.RLock()
    previous_seconds = state.get('compute_seconds', 0)
    started = time.monotonic()
    settings = config.get('research_docking') or {}
    python = settings.get('preparation_python') or config.get('preparation_python') or sys.executable
    by_id = {r['candidate_id']: r for r in rows}

    def remaining():
        return max(0, policy['budget_seconds'] - previous_seconds - (time.monotonic() - started))

    def save():
        with checkpoint_lock:
            state['compute_seconds'] = previous_seconds + time.monotonic() - started
            state['receipt'] = deepcopy(receipt)
            atomic_bytes(checkpoint, canonical_json(state))

    def restore(cid, stored, reused=False, cache_artifact_sha256=None):
        row = by_id[cid]
        for k in ('docking_execution', 'docking', 'primary_docking', 'routine_pose_qc',
                  'routine_pose_quality', 'routine_quality_software', 'ligand_preparation_sha256', 'routine_reuse'):
            if k in stored:
                row[k] = deepcopy(stored[k])
        if reused:
            row['routine_reuse'] = {'source_candidate_id': stored['candidate_id'],
                                    'source_discovery_run_id': stored.get('discovery_run_id'),
                                    'cached_row_artifact_sha256': cache_artifact_sha256,
                                    'exact_structure_and_protocol': True}
            row['docking_execution']['execution_candidate_id'] = stored['candidate_id']
            row['docking_execution']['candidate_id'] = cid
            if row.get('docking'):
                row['docking'][0]['execution_candidate_id'] = stored['candidate_id']
                row['docking'][0]['candidate_id'] = cid
                row['docking'][0]['observation_source'] = 'stored_exact_protocol_reuse'
                row['primary_docking'] = deepcopy(row['docking'][0])

    def valid_cached(stored, expected):
        try:
            execution = stored.get('docking_execution') or {}
            quality = stored.get('routine_pose_quality') or {}
            if stored.get('routine_pose_qc') not in {'pass', 'failed'}:
                return False
            sdf = execution.get('sdf_artifact', {}).get('sha256')
            inputs = quality.get('metrics', {}).get('inputs', {})
            software = stored.get('routine_quality_software') or {}
            if (not sdf or inputs.get('pose_sha256') != sdf or quality.get('status') != 'complete'
                    or software.get('packages', {}).get('posebusters', {}).get('version') != '0.6.5'
                    or stored.get('target_id') != expected.get('target_id')
                    or chemical_identity(stored.get('smiles'))['exact_structure_id'] != chemical_identity(expected['smiles'])['exact_structure_id']
                    or inputs.get('prepared_exact_structure_id') != chemical_identity(expected['smiles'])['exact_structure_id']
                    or inputs.get('receptor_sha256') != descriptor['protein_structure_sha256']
                    or any(execution.get('protocol', {}).get(k) != v for k, v in engine_request.items())):
                return False
            prepared = json.loads(service.artifacts.get_bytes(stored['ligand_preparation_sha256']))['request']
            if {k: v for k, v in prepared.items() if k != 'smiles'} != preparation_context:
                return False
            service.artifacts.get_bytes(sdf)
            service.artifacts.get_bytes(execution['pose_sha256'])
            for artifact in quality.get('artifacts', []):
                service.artifacts.get_bytes(artifact['sha256'])
        except (OSError, ValueError, KeyError, TypeError):
            return False
        return True

    def prepare(items):
        missing = [r for r in items if r['candidate_id'] not in state['prepared']]
        for offset in range(0, len(missing), 32):
            part = native._runtime_step({'operation': 'ligands', 'workdir': str(work),
                'ligands': [{'candidate_id': r['candidate_id'], 'smiles': r['smiles']}
                            for r in missing[offset:offset + 32]]}, work, python,
                stop=stop, timeout=max(1, min(360, remaining())))
            for item in part['ligands']:
                cid = item['candidate_id']
                if not item.get('path'):
                    receipt.setdefault('preparation_holds', []).append(item)
                    state['prepared'][cid] = None
                    continue
                path = Path(item['path']).resolve()
                if not path.is_relative_to(work.resolve()):
                    raise ValueError('Prepared ligand escaped owned routine directory')
                prep = read_preparation_receipt(path, ligand_id=cid)
                if (prep is None or prep['request']['ph'] is not None
                        or prep['input_smiles'] != by_id[cid]['smiles']
                        or chemical_identity(prep['effective_smiles'])['exact_structure_id']
                        != chemical_identity(by_id[cid]['smiles'])['exact_structure_id']):
                    raise ValueError('Routine preparation differs from the frozen exact input')
                artifact = service.artifacts.put_file(path, media_type='chemical/x-pdbqt')
                by_id[cid]['ligand_preparation_sha256'] = service.artifacts.put_bytes(
                    canonical_json(prep), media_type='application/json')['sha256']
                state['prepared'][cid] = {'candidate_id': cid, 'pdbqt_sha256': artifact['sha256'],
                    'preparation_receipt_sha256': by_id[cid]['ligand_preparation_sha256']}
            save()
        for item in items:
            saved = state['prepared'].get(item['candidate_id'])
            if saved:
                service.artifacts.get_bytes(saved['pdbqt_sha256'])
                by_id[item['candidate_id']]['ligand_preparation_sha256'] = saved['preparation_receipt_sha256']
        return [{k: v for k, v in state['prepared'][r['candidate_id']].items()
                 if k != 'preparation_receipt_sha256'} for r in items if state['prepared'].get(r['candidate_id'])]

    try:
        stop()
        with _owned_cache_lock(cache / 'compute.lock', stop,
                               on_wait=lambda: progress({'stage': 'routine_cache_wait'})):
            progress({'stage': 'routine_receptor_preparation'})
            source, observed = work / 'selected-protein.pdb', work / 'observed-query.sdf'
            atomic_bytes(source, service.artifacts.get_bytes(descriptor['protein_structure_sha256']))
            seed = next(r for r in rows if r.get('role') == 'query_seed' or r.get('redock_reference') is True)
            atomic_bytes(observed, service.artifacts.get_bytes(seed['observed_sdf_sha256']))
            receptor = native._runtime_step({'operation': 'receptor', 'source_path': str(source),
                'output_path': str(work / 'receptor.pdbqt')}, work, python, stop=stop,
                timeout=max(1, min(360, remaining())))
            receptor_artifact = service.artifacts.put_file(receptor['path'], media_type='chemical/x-pdbqt')
            receipt.update(receptor_valid=True, receptor_preparation=receptor,
                           receptor_sha256=receptor_artifact['sha256'])
            card = {'schema': 'skinpipe-observed-complex-receptor-card-v1',
                'artifacts': {'pdbqt': {'present': True, 'sha256': receptor_artifact['sha256']},
                    'cocrystal_sdf': {'present': True, 'sha256': seed['observed_sdf_sha256'], 'path': 'observed-query.sdf'}},
                'receptor': {'cocrystal_sdf': 'observed-query.sdf'},
                'selected_structure_sha256': descriptor['protein_structure_sha256'], 'pocket_box': descriptor['pocket_box']}
            card_sha = service.artifacts.put_bytes(canonical_json(card), media_type='application/json')['sha256']
            receipt['receptor_identity_card_sha256'] = card_sha
            engine_request = {'engine': frozen['engine'], 'receptor_sha256': receptor_artifact['sha256'],
                'center': descriptor['pocket_box']['center'], 'size': descriptor['pocket_box']['size'],
                'seed': frozen['seed'], 'n_poses': frozen['n_poses'], 'exhaustiveness': frozen['exhaustiveness']}
            runner = EngineRunner(work / 'engine', settings['engine_paths'], gpu_uuid=settings.get('gpu_uuid'),
                                  cpu_threads=policy['cpu_threads'])

            def observe_progress(value):
                partial = value.get('partial_result') or {}
                with checkpoint_lock:
                    for row in partial.get('candidates', []):
                        state['engine_rows'][row['candidate_id']] = deepcopy(row)
                    if partial:
                        save()
                progress({**value, 'stage': 'routine_docking', 'target': policy['target_candidates']})

            def dock(ligands):
                missing = [r for r in ligands if r['candidate_id'] not in state['engine_rows']]
                if missing:
                    result = native._execute_engine_fenced(runner, {**engine_request, 'ligands': missing},
                        service.artifacts, stop=stop, progress=observe_progress, remaining_seconds=remaining())
                    for row in result['candidates']:
                        state['engine_rows'][row['candidate_id']] = deepcopy(row)
                    save()
                return [deepcopy(state['engine_rows'][r['candidate_id']]) for r in ligands]

            query = prepare([seed])
            if not query:
                return rows, {**receipt, 'hold_reasons': ['query_redocking_preparation_failed']}
            redocked = dock(query)[0]
            seed['docking_execution'] = deepcopy(redocked)
            seed['docking'] = [native._docking_observation(seed, redocked, descriptor, receptor_artifact['sha256'], card_sha)]
            seed['primary_docking'] = seed['docking'][0]
            pose = work / 'redocked-query.pdbqt'
            if not redocked.get('pose_sha256'):
                return rows, {**receipt, 'hold_reasons': ['query_redocking_engine_failed']}
            atomic_bytes(pose, service.artifacts.get_bytes(redocked['pose_sha256']))
            validation = native._runtime_step({'operation': 'redock_geometry', 'observed_sdf': str(observed),
                'pose_path': str(pose), 'pose_scores': redocked['pose_scores'],
                'output_sdf': str(work / 'redocked-query.sdf')}, work, python, stop=stop,
                timeout=max(1, min(360, remaining())))
            receipt['redocking'] = validation
            if not validation.get('pass'):
                save()
                return rows, {**receipt, 'hold_reasons': ['query_redocking_geometry_failed']}
            candidates = [by_id[c] for c in selected[:policy['max_new_docking']] if c != seed['candidate_id']]
            preparation_request = json.loads(service.artifacts.get_bytes(seed['ligand_preparation_sha256']))['request']
            preparation_context = {k: v for k, v in preparation_request.items() if k != 'smiles'}
            if preparation_context.get('software_identity', {}).get('package_inputs_state') != 'observed':
                raise PreparationUnavailable('Preparation software inputs are not verified')
            execution_context = {**current_sources, 'preparation': preparation_context}
            keys = {r['candidate_id']: _cache_key(frozen, engine_request, r,
                descriptor['protein_structure_sha256'], execution_context) for r in candidates}
            pending = []
            for row in candidates:
                cid = row['candidate_id']
                stored = state['rows'].get(cid)
                if stored and valid_cached(stored, row):
                    restore(cid, stored)
                    receipt['reused_count'] += bool(stored.get('routine_reuse'))
                    continue
                path = cache / (keys[cid] + '.json')
                if path.exists():
                    try:
                        saved = json.loads(path.read_bytes())
                        stored = json.loads(service.artifacts.get_bytes(saved['row_artifact_sha256']))
                    except (OSError, ValueError, KeyError, TypeError):
                        saved, stored = {}, None
                    if saved.get('cache_key') == keys[cid] and stored and valid_cached(stored, row):
                        restore(cid, stored, reused=True, cache_artifact_sha256=saved['row_artifact_sha256'])
                        state['rows'][cid] = deepcopy(row)
                        receipt['reused_count'] += 1
                        continue
                pending.append(row)
            attempts = len([cid for cid in state['engine_rows'] if cid != seed['candidate_id']])
            rate = None
            while pending:
                stop()
                selection = select_representatives(rows, policy)
                admitted = next_compute_batch(policy, qualified=selection['selected_count'], attempted=attempts,
                    pending=len(pending), remaining_seconds=remaining(),
                    measured_seconds_per_candidate=rate, benchmark_context_matches=rate is not None)
                if not admitted['admitted']:
                    receipt['stop_reason'] = admitted['status']
                    break
                batch_started = time.monotonic()
                batch, pending = pending[:admitted['admitted']], pending[admitted['admitted']:]
                ligands = prepare(batch)
                attempts += len(batch)
                results = dock(ligands) if ligands else []
                quality_rows = []
                for result in results:
                    row = by_id[result['candidate_id']]
                    row['docking_execution'] = deepcopy(result)
                    if result.get('status') not in {'complete', 'review_required'}:
                        continue
                    path = work / ('pose-' + digest(row['candidate_id'])[:24] + '.pdbqt')
                    atomic_bytes(path, service.artifacts.get_bytes(result['pose_sha256']))
                    try:
                        geometry = native._runtime_step({'operation': 'candidate_geometry',
                            'pose_path': str(path), 'source_smiles': row['smiles'],
                            'pose_scores': result['pose_scores'], 'output_sdf': str(path.with_suffix('.sdf'))},
                            work, python, stop=stop, timeout=max(1, min(360, remaining())))
                        sdf = service.artifacts.put_file(geometry.pop('top_pose_sdf'), media_type='chemical/x-mdl-sdfile')
                        row['docking_execution'].update(sdf_artifact=sdf, pose_geometry=geometry)
                        quality_rows.append({'candidate_id': row['candidate_id'], 'smiles': row['smiles'],
                            'pose_path': str(path.with_suffix('.sdf')), 'pose_sha256': sdf['sha256'],
                            'prepared_exact_structure_id': chemical_identity(row['smiles'])['exact_structure_id']})
                    except PreparationUnavailable as error:
                        row['routine_pose_qc'] = 'held'
                        row['docking_execution']['sdf_status'] = str(error)[:256]
                    row['docking'] = [native._docking_observation(row, row['docking_execution'], descriptor,
                        receptor_artifact['sha256'], card_sha)]
                    row['primary_docking'] = deepcopy(row['docking'][0])
                if quality_rows:
                    quality = native._runtime_step({'operation': 'routine_quality', 'rows': quality_rows,
                        'receptor_path': str(source), 'receptor_sha256': sha256_file(source),
                        'budget_seconds': max(1, min(360, remaining())), 'workdir': str(work / 'quality')},
                        work, python, stop=stop, timeout=max(1, min(360, remaining())))
                    if quality['quality'].get('software', {}).get('packages', {}).get('posebusters', {}).get('version') != '0.6.5':
                        raise PreparationUnavailable('Routine pose QC requires PoseBusters 0.6.5')
                    for q in quality['quality']['candidates']:
                        row = by_id[q['candidate_id']]
                        for artifact in q.get('artifacts', []):
                            artifact_path = Path(artifact['path']).resolve()
                            if not artifact_path.is_relative_to(work.resolve()):
                                raise ValueError('Quality artifact escaped the owned run directory')
                            saved_artifact = service.artifacts.put_file(artifact_path, media_type='application/json')
                            if saved_artifact['sha256'] != artifact['sha256']:
                                raise ValueError('Quality artifact changed before publication')
                        row['routine_pose_quality'] = q
                        row['routine_quality_software'] = quality['quality']['software']
                        row['routine_pose_qc'] = q.get('metrics', {}).get('pose_check_status', 'held')
                for row in batch:
                    state['rows'][row['candidate_id']] = deepcopy(row)
                    if valid_cached(row, row):
                        artifact = service.artifacts.put_bytes(canonical_json(row), media_type='application/json')
                        atomic_bytes(cache / (keys[row['candidate_id']] + '.json'), canonical_json({
                            'cache_key': keys[row['candidate_id']], 'row_artifact_sha256': artifact['sha256']}))
                receipt['new_docking_count'] = len([cid for cid in state['engine_rows'] if cid != seed['candidate_id']])
                rate = (time.monotonic() - batch_started) / len(batch)
                receipt['measured_seconds_per_candidate'] = rate
                receipt['qualified_count'] = select_representatives(rows, policy)['selected_count']
                save()
                progress({'stage': 'routine_pose_quality', 'processed': attempts,
                          'qualified': receipt['qualified_count'], 'target': policy['target_candidates']})
                if quality_rows and not any(by_id[r['candidate_id']].get('routine_pose_qc') in {'pass', 'failed'} for r in quality_rows):
                    receipt['stop_reason'] = 'pose_quality_runtime_or_input_held'
                    break
            receipt['new_docking_count'] = len([cid for cid in state['engine_rows'] if cid != seed['candidate_id']])
            receipt['status'] = 'complete'
            receipt['planned_candidate_count'] = len(candidates)
            receipt['selected_count'] = sum(bool(r.get('docking_execution')) for r in candidates)
            receipt['prepared_count'] = sum(v is not None for v in state['prepared'].values())
            receipt['candidate_results'] = [deepcopy(r['docking_execution']) for r in candidates if r.get('docking_execution')]
            receipt['completed_count'] = sum(r.get('docking_execution', {}).get('status') == 'complete' for r in candidates)
            receipt['qualified_count'] = select_representatives(rows, policy)['selected_count']
            receipt.setdefault('stop_reason', 'target_reached' if receipt['qualified_count'] >= policy['target_candidates'] else 'eligible_pool_exhausted')
            save()
            return rows, receipt
    except TimeoutError:
        save()
        raise
    except ComputeBudgetExceeded:
        save()
        raise TimeoutError('budget_deferred') from None
    except (PreparationUnavailable, ValueError, OSError, ImportError) as error:
        save()
        return rows, {**receipt, 'hold_reasons': ['routine_docking_execution_held'],
                      'error_type': type(error).__name__, 'diagnostic': str(error)[:512]}
