"""Bounded scientific receipt/export regressions; no native model execution."""
from copy import deepcopy
import json
from types import SimpleNamespace

from openpyxl import load_workbook
import pytest

from skinpipe.platform.admet_search import RISK_ENDPOINTS
from skinpipe.platform.artifacts import ArtifactStore
from skinpipe.platform.catalog import Catalog
from skinpipe.platform.exports import RunExports
from skinpipe.platform.identity import chemical_identity
from skinpipe.platform.routine_screening import eligibility, select_representatives, screening_policy


def candidate():
    return {'candidate_id': 'fixture', 'smiles': 'C[C@H](O)F', 'target_id': 'fixture-target',
            'priority_score_v2': 75, 'identity': chemical_identity('C[C@H](O)F'),
            'routine_route': 'ligand_only', 'skin': {'topical_score': 5},
            'structural_admet_context': {'product': None}, 'admet': [{
                'status': 'success', 'model_id': 'admet-ai-fixture',
                'input_smiles': 'C[C@H](O)F', 'model_manifest_sha256': 'a' * 64,
                'endpoint_definitions_sha256': 'b' * 64, 'values': {'AMES': 0}}]}


def additional(row, values=None):
    return {'status': 'success', 'smiles': row['smiles'],
            'exact_structure_id': row['identity']['exact_structure_id'],
            'values': {'AMES': 0} if values is None else values,
            'physicochemical_values': {'logP': -1.2},
            'model_identity': {'model_id': 'additional-fixture',
                               'model_manifest_sha256': 'c' * 64,
                               'endpoint_definitions_sha256': 'd' * 64}}


@pytest.mark.parametrize('endpoint', RISK_ENDPOINTS)
@pytest.mark.parametrize('value', [-0.1, 1.1, True, 10**400, None])
@pytest.mark.parametrize('channel', ['cached', 'additional'])
def test_invalid_classifier_receipt_cannot_fill_representative_quota(endpoint, value, channel):
    row = candidate()
    if channel == 'cached':
        row['admet'][0]['values'] = {endpoint: value}
    else:
        row['admet'] = []
        row['additional_toxicity'] = additional(row, {endpoint: value})
    before = deepcopy(row)
    result = select_representatives([row], screening_policy(target_candidates=1))
    assert result['selected_count'] == 0 and result['shortfall'] == 1
    assert result['decisions']['fixture'] == {
        'state': 'pending', 'reasons': ['exact_model_predictions_required']}
    assert row == before


@pytest.mark.parametrize('value', [0, 1])
@pytest.mark.parametrize('channel', ['cached', 'additional'])
def test_classifier_boundaries_and_negative_regressions_remain_valid(value, channel):
    row = candidate()
    values = {'AMES': value, 'Solubility_AqSolDB': -7.25}
    if channel == 'cached':
        row['admet'][0]['values'] = values
    else:
        row['admet'] = []
        row['additional_toxicity'] = additional(row, values)
    assert eligibility(row) == ('eligible', [])


@pytest.mark.parametrize('channel', ['cached', 'additional'])
@pytest.mark.parametrize('definition', [{'task_type': 'classification'}, {'units': 'probability'},
                                       {'unit': 'probability'}])
def test_declared_classifier_domain_is_honored_beyond_search_risk_endpoints(channel, definition):
    row = candidate()
    record = row['admet'][0]
    if channel == 'additional':
        row['admet'] = []
        row['additional_toxicity'] = record = additional(row)
    record['values'] = {'BBB_Martins': 1.1}
    record['endpoint_definitions'] = [{'id': 'BBB_Martins', **definition}]
    before = deepcopy(row)
    assert eligibility(row) == ('pending', ['exact_model_predictions_required'])
    assert row == before


@pytest.mark.parametrize('definition', [None, {'task_type': 'regression'}])
def test_unknown_or_declared_regression_endpoint_is_not_probability_clamped(definition):
    row = candidate()
    record = row['admet'][0]
    record['values'] = {'synthetic_regression': -5.3}
    if definition:
        record['endpoint_definitions'] = [{'id': 'synthetic_regression', **definition}]
    assert eligibility(row) == ('eligible', [])


@pytest.mark.parametrize('bad', [42, True, 'malformed', [42], {'values': None}])
def test_malformed_cached_containers_remain_pending_without_changing_source(bad):
    row = candidate()
    row['admet'] = bad
    before = deepcopy(row)
    assert eligibility(row) == ('pending', ['exact_model_predictions_required'])
    assert row == before


@pytest.mark.parametrize('bad', [42, True, 'malformed', {'status': 'failed', 'values': None},
                                 {'status': 'success', 'model_identity': 42}])
def test_bad_additional_record_does_not_abort_valid_independent_channel(bad):
    row = candidate()
    row['additional_toxicity'] = bad
    before = deepcopy(row)
    assert eligibility(row) == ('eligible', [])
    assert row == before


def test_physicochemical_values_alone_do_not_establish_admet_predictions():
    row = candidate()
    row['admet'] = []
    row['additional_toxicity'] = additional(row, {})
    assert eligibility(row) == ('pending', ['exact_model_predictions_required'])


@pytest.mark.parametrize('value', [0, 2])
def test_computed_property_cannot_shadow_a_learned_endpoint(value):
    row = candidate()
    row['admet'] = []
    row['additional_toxicity'] = additional(row)
    row['additional_toxicity']['physicochemical_values'] = {'AMES': value}
    before = deepcopy(row)
    assert eligibility(row) == ('pending', ['exact_model_predictions_required'])
    assert row == before


@pytest.mark.parametrize('bad', [42, True, 'malformed', [42]])
def test_malformed_cache_survives_optional_toxicity_stage_without_becoming_a_prediction(bad):
    from skinpipe.platform.research_toxicity import apply_toxicity
    row = candidate()
    row['admet'] = bad
    rows, receipt = apply_toxicity([row], None, {}, stop=lambda: None, progress=lambda _: None)
    assert rows[0]['admet'] == bad
    assert rows[0]['additional_toxicity']['status'] == 'unassessed_model_not_configured'
    assert receipt['skipped_cached_count'] == 0
    assert eligibility(rows[0]) == ('pending', ['exact_model_predictions_required'])


@pytest.mark.parametrize('field,key', [('values', 'AMES'), ('physicochemical_values', 'molecular_weight')])
def test_overflowing_additional_value_does_not_abort_whole_workbook(tmp_path, field, key):
    row = candidate()
    row['additional_toxicity'] = additional(row)
    row['additional_toxicity'][field][key] = 10**400
    catalog = Catalog(tmp_path / 'catalog')
    manifest = catalog.create_snapshot([row])
    before = list(catalog.iter_candidates(manifest['snapshot_id'], target_id='fixture-target'))
    exporter = RunExports(SimpleNamespace(root=tmp_path, catalog=catalog,
                                          artifacts=ArtifactStore(tmp_path / 'cas')))
    output = tmp_path / 'export'
    receipt = exporter.build({'snapshot_id': manifest['snapshot_id'], 'target_id': 'fixture-target',
                              'format': 'xlsx', 'candidate_details': True}, output)
    assert receipt['candidate_count'] == 1
    saved = json.loads((output / receipt['chunks'][0]['file']).read_bytes())[0]
    assert saved['additional_toxicity'] == row['additional_toxicity']
    assert list(catalog.iter_candidates(manifest['snapshot_id'], target_id='fixture-target')) == before
    book = load_workbook(output / receipt['workbooks'][0])
    try:
        cells = [str(cell.value) for cells in book['C000001'] for cell in cells if cell.value]
        assert any('유효 범위·자료형 위반 · 표시 보류' in value for value in cells)
        assert any(row['smiles'] == cell.value for cells in book['후보 검토'] for cell in cells)
    finally:
        book.close()
