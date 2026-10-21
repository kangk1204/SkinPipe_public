"""Bounded routine triage and incremental computation, with explicit missingness.

New research inputs use this policy; existing frozen runs keep their settings. Selection uses
the saved within-study priority, not raw docking scores across protocols. Exact
structure identity controls reuse; a scaffold is only a diversity grouping.
"""
from __future__ import annotations

from collections import Counter
from copy import deepcopy
import math
from pathlib import Path
import re

from .identity import chemical_identity, digest
from .ranking import finite

VERSION = 'routine-screening-v1'
CONTROL_ROLES = {'control', 'positive_control', 'negative_control', 'counter_control',
                 'reference', 'reference_control', 'query_seed'}


def execution_sources():
    from .artifacts import sha256_file
    directory = Path(__file__).resolve().parent
    return {name: sha256_file(directory / name) for name in
            ('routine_screening.py', 'routine_docking.py', 'algorithm_chemistry.py', 'engines.py')}


def screening_policy(*, target_candidates=200, compute_budget_seconds=7200):
    if type(target_candidates) is not int or not 1 <= target_candidates <= 1000:
        raise ValueError('target_candidates must be an integer within 1..1000')
    if (type(compute_budget_seconds) not in (int, float)
            or not math.isfinite(compute_budget_seconds) or not 0 < compute_budget_seconds <= 28800):
        raise ValueError('compute budget must be within 8h')
    value = {'version': VERSION, 'target_candidates': target_candidates,
        'budget_seconds': compute_budget_seconds, 'max_new_docking': 1000,
        'canary_candidates': 8, 'batch_candidates': 32, 'cost_reserve': 1.25,
        'score_fraction': 0.8, 'cache': 'exact_structure_model_preprocessing_and_protocol',
        'cpu_threads': 1,
        'missing': 'pending_not_pass', 'model_hold': 'do_not_dispatch_docking',
        'stop': ['qualified_target_reached', 'budget_exhausted', 'eligible_pool_exhausted'],
        'optional_deep_models': False, 'controls_outside_candidate_quota': True,
        'execution_sources': execution_sources(),
        'scope': 'one study and calculation condition; preserve separate observations'}
    return {**value, 'sha256': digest(value)}


def freeze_screening(body=None, *, compute_budget_seconds=7200):
    body = {} if body is None else body
    if not isinstance(body, dict) or set(body) - {'target_candidates'}:
        raise ValueError('screening accepts target_candidates only')
    return screening_policy(target_candidates=body.get('target_candidates', 200),
                            compute_budget_seconds=compute_budget_seconds)


def validate_policy(policy):
    if not isinstance(policy, dict):
        raise ValueError('Routine screening policy differs from the frozen version')
    payload = {key: value for key, value in policy.items() if key != 'sha256'}
    expected = screening_policy(target_candidates=policy.get('target_candidates'),
                                compute_budget_seconds=policy.get('budget_seconds'))
    sources = policy.get('execution_sources')
    if (not isinstance(sources, dict) or set(sources) != set(expected['execution_sources'])
            or any(not isinstance(value, str) or not re.fullmatch(r'[0-9a-f]{64}', value)
                   for value in sources.values()) or policy.get('sha256') != digest(payload)):
        raise ValueError('Routine screening policy integrity check failed')
    # Source provenance is frozen at submission; compatible policy semantics
    # survive updates. Actual execution/cache/checkpoint identities remain current.
    expected['execution_sources'] = deepcopy(sources)
    if payload != {key: value for key, value in expected.items() if key != 'sha256'}:
        raise ValueError('Routine screening policy differs from the frozen version')
    return deepcopy(policy)


def exact_admet_channels(row, identity):
    """Keep distinct models; a name, parent or InChIKey cannot establish reuse."""
    records = row.get('admet') or []
    records = [records] if isinstance(records, dict) else records
    cached = [deepcopy(p) for p in records if isinstance(p, dict)
        and p.get('status') == 'success' and p.get('model_id')
        and p.get('model_manifest_sha256') and p.get('endpoint_definitions_sha256')
        and isinstance(p.get('values'), dict) and p['values']
        and all(finite(v) is not None for v in p['values'].values())
        and p.get('input_smiles')
        and chemical_identity(p['input_smiles'])['exact_structure_id'] == identity['exact_structure_id']]
    extra = row.get('additional_toxicity') or {}
    model = extra.get('model_identity') or {}
    values = {**extra.get('values', {}), **extra.get('physicochemical_values', {})}
    if (extra.get('status') == 'success' and extra.get('exact_structure_id') == identity['exact_structure_id']
            and extra.get('smiles') and chemical_identity(extra['smiles'])['exact_structure_id'] == identity['exact_structure_id']
            and model.get('model_manifest_sha256') and model.get('endpoint_definitions_sha256')
            and values and all(finite(v) is not None for v in values.values())):
        cached.append({'channel': 'additional_model', 'model_id': model.get('model_id'),
                       'model_identity': deepcopy(model), 'values': values})
    return cached


def eligibility(row, *, model_hold=False):
    """Return a reasoned disposition; no observed failures are overwritten."""
    if row.get('role') in CONTROL_ROLES or row['candidate_id'].startswith('ref-'):
        return 'control', ['control_outside_candidate_quota']
    reasons = []
    identity = chemical_identity(row.get('smiles'))
    if identity['structure_status'] != 'valid' or identity.get('parent_ambiguous'):
        reasons.append('unresolved_structure')
    if identity.get('unassigned_stereo') or identity.get('unassigned_bond_stereo'):
        reasons.append('explicit_stereoisomer_required')
    original = row.get('identity') or {}
    if original.get('exact_structure_id') not in (None, identity['exact_structure_id']):
        reasons.append('source_exact_structure_mismatch')
    if model_hold:
        reasons.append('receptor_or_calculation_model_held')
    if row.get('routine_pose_qc') == 'failed':
        reasons.append('pose_QC_failed')
    execution = row.get('docking_execution') or row.get('primary_docking') or (row.get('reuse_evidence') or {}).get('source_observation') or {}
    if execution.get('status') == 'review_required' or execution.get('ranking_eligible') is False:
        reasons.append('docking_score_review_required')
    elif execution.get('status') in {'failed', 'cancelled'}:
        reasons.append('docking_execution_incomplete')
    if row.get('discovery_filter_pass') is False:
        reasons.extend(row.get('discovery_filter_reasons') or ['frozen_search_filter'])
    if row.get('critical_warning') is True or row.get('hazard_severity') in {'critical', 'high'}:
        reasons.append('structural_alert_review_required')
    if row.get('opposing_direction') is True:
        reasons.append('opposing_direction_control')
    if (row.get('admet_search_risk') or {}).get('excluded') is True:
        reasons.append('severe_predicted_admet_review_required')
    if reasons:
        return 'held', sorted(set(reasons))
    missing = []
    if not exact_admet_channels(row, identity):
        missing.append('exact_model_predictions_required')
    if row.get('routine_pose_qc') != 'pass' and row.get('routine_route') != 'ligand_only':
        missing.append('exact_pose_QC_required')
    if finite(row.get('priority_score_v2')) is None:
        missing.append('saved_priority_required')
    if not row.get('structural_admet_context') or not row.get('skin'):
        missing.append('skin_structure_assessment_required')
    return ('pending', missing) if missing else ('eligible', [])


def _order(row):
    return (-row['priority_score_v2'], row['candidate_id'])


def _core(row):
    identity = chemical_identity(row['smiles'])
    if identity.get('scaffold_kind') == 'acyclic':
        from .catalog_hopping import acyclic_topology
        return ('acyclic_topology', (acyclic_topology(row['smiles']) or {}).get('id', 'acyclic'))
    return ('ring_scaffold', identity.get('scaffold_id'))


def select_representatives(rows, policy, *, model_hold=False):
    """80% saved priority, then new scaffold/topology groups, stable backfill.

    Every source row receives a disposition, including missing, held, controls,
    duplicate exact structures and review alternatives. No source is mutated.
    """
    policy = validate_policy(policy)
    if len({r['candidate_id'] for r in rows}) != len(rows):
        raise ValueError('Duplicate candidate IDs require observation-level separation')
    conditions = {(r.get('study'), r.get('condition'), r.get('target_id')) for r in rows}
    if len(conditions) > 1:
        raise ValueError('Select each study/condition/target separately')
    decisions, eligible = {}, []
    for row in rows:
        state, reasons = eligibility(row, model_hold=model_hold)
        decisions[row['candidate_id']] = {'state': state, 'reasons': reasons}
        if state == 'eligible':
            eligible.append(row)
    eligible.sort(key=_order)
    # Controls do not re-enter the list through a catalog alias.
    reserved = {chemical_identity(r.get('smiles'))['exact_structure_id'] for r in rows
                if decisions[r['candidate_id']]['state'] == 'control'}
    exact_seen, unique = set(reserved), []
    for row in eligible:
        exact = chemical_identity(row['smiles'])['exact_structure_id']
        if exact in exact_seen:
            decisions[row['candidate_id']] = {'state': 'duplicate', 'reasons': ['exact_structure_already_represented']}
        else:
            exact_seen.add(exact)
            unique.append(row)
    selected, used, cores = [], set(), set()
    cap = policy['target_candidates']

    def take(row, group):
        if row['candidate_id'] in used or len(selected) >= cap:
            return
        used.add(row['candidate_id'])
        cores.add(_core(row))
        selected.append(row['candidate_id'])
        decisions[row['candidate_id']] = {'state': 'selected', 'reasons': [group]}

    for row in unique[:math.ceil(cap * policy['score_fraction'])]:
        take(row, 'saved_priority')
    for row in unique:
        if _core(row) not in cores:
            take(row, 'scaffold_or_acyclic_topology_diversity')
    for row in unique:
        take(row, 'eligible_priority_backfill')
    for row in unique:
        if row['candidate_id'] not in used:
            decisions[row['candidate_id']] = {'state': 'alternative', 'reasons': ['representative_limit']}
    return {'policy': policy, 'selected_ids': selected, 'selected_count': len(selected),
        'counts': dict(Counter(d['state'] for d in decisions.values())), 'decisions': decisions,
        'shortfall': max(0, cap - len(selected)), 'source_population_preserved': True}


def next_compute_batch(policy, *, qualified, attempted, pending, remaining_seconds,
                       measured_seconds_per_candidate=None, benchmark_context_matches=False,
                       model_hold=False):
    """Admit a bounded next batch; no elapsed-time guess or automatic deep run."""
    policy = validate_policy(policy)
    if any(type(v) is not int or v < 0 for v in (qualified, attempted, pending)):
        raise ValueError('Candidate counts must be nonnegative integers')
    if finite(remaining_seconds) is None or remaining_seconds < 0:
        raise ValueError('Remaining compute budget must be finite and nonnegative')
    if model_hold:
        return {'status': 'model_held', 'admitted': 0}
    if qualified >= policy['target_candidates']:
        return {'status': 'target_reached', 'admitted': 0}
    available = min(pending, max(0, policy['max_new_docking'] - attempted))
    if not available:
        return {'status': 'eligible_pool_exhausted', 'admitted': 0}
    if not remaining_seconds:
        return {'status': 'budget_deferred', 'admitted': 0}
    rate = finite(measured_seconds_per_candidate)
    if not benchmark_context_matches or rate is None or rate <= 0:
        return {'status': 'canary_required', 'admitted': min(available, policy['canary_candidates']),
                'estimated_compute_seconds': None}
    fit = math.floor(remaining_seconds / (rate * policy['cost_reserve']))
    size = min(available, policy['batch_candidates'], fit)
    return {'status': 'ready' if size else 'budget_deferred', 'admitted': size,
            'estimated_compute_seconds': size * rate * policy['cost_reserve']}
