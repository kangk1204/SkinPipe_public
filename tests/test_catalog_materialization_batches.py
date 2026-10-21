"""Large retrievals retain exact populations while bounding gateway transfers."""
from copy import deepcopy
import hashlib
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from skinpipe.platform.artifacts import ArtifactStore, canonical_json
from skinpipe.platform.catalog_cancellation import owned_catalog_scope
from skinpipe.platform.discovery import DiscoveryRuns, freeze_filters
from skinpipe.platform.engines import OwnedProcessCleanupError
from skinpipe.platform.identity import digest
from skinpipe.platform.queue import Cancelled
from skinpipe.platform.service import PlatformService, _CatalogTransportError


class Gateway:
    def __init__(self, modify=None):
        self.calls, self.responses, self.cancelled = [], [], []
        self.modify = modify

    def __call__(self, path, *, body=None, response_limit=None):
        if path.endswith('/cancel'):
            self.cancelled.append(path.split('/')[2])
            return {'status': 'cancelled', 'query_id': path.split('/')[2], 'resources_released': True}
        if path != '/materialize':
            raise AssertionError('Unexpected query or retry: ' + path)
        self.calls.append((deepcopy(body), response_limit))
        ids = body['compound_ids']
        if len(ids) > 50:
            raise ValueError('Synthetic bounded transfer exceeded')
        result = {
            'status': 'complete', 'query_id': format(len(self.calls), '032x'),
            'rows': [{'compound': {'compound_id': cid, 'canonical_smiles': 'CCO'},
                      'offers': [], 'admet': []} for cid in ids],
            'source': {'report_sha256': 'b' * 64, 'nested': {'release': 7}},
            'elapsed_seconds': 0.25, 'query_scope': 'full_source_catalog',
            'retained_extension': 'original-contract',
        }
        if self.modify is not None:
            self.modify(result, len(self.calls))
        self.responses.append(result)
        return result


class CatalogMaterializationBatches(unittest.TestCase):
    def setUp(self):
        self.tmp = TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.service = PlatformService.__new__(PlatformService)
        self.service.scout = {'url': 'https://catalog.example.invalid', 'token': 'synthetic'}
        self.service.artifacts = ArtifactStore(Path(self.tmp.name) / 'artifacts')
        self.gateway = Gateway()
        self.service._scout_request_impl = self.gateway

    @staticmethod
    def ids(count):
        return [f'candidate-{i:04d}' for i in range(count)]

    def test_fifty_preserves_the_existing_response_contract(self):
        result = self.service._scout_materialization_chunk(self.ids(50))
        self.assertIs(result, self.gateway.responses[0])
        self.assertNotIn('materialization_chunks', result)
        self.assertEqual(len(self.gateway.calls), 1)
        self.assertEqual(self.gateway.calls[0][1], 64 * 1024**2)

    def test_thousand_preserves_order_source_and_every_cleanup_receipt(self):
        ids = self.ids(1000)
        with owned_catalog_scope(self.service):
            result = self.service._scout_materialization_chunk(ids)
        self.assertEqual([row['compound']['compound_id'] for row in result['rows']], ids)
        self.assertEqual(len(self.gateway.calls), 20)
        self.assertTrue(all(len(body['compound_ids']) == 50 and limit == 64 * 1024**2
                            for body, limit in self.gateway.calls))
        self.assertEqual(result['source'], self.gateway.responses[0]['source'])
        self.assertEqual(result['query_scope'], 'full_source_catalog')
        self.assertEqual(result['retained_extension'], 'original-contract')
        self.assertEqual(result['rows'][:50], self.gateway.responses[0]['rows'])
        chunks = result['materialization_chunks']
        self.assertEqual(sum(row['candidate_count'] for row in chunks), 1000)
        self.assertEqual({row['source_sha256'] for row in chunks},
                         {hashlib.sha256(canonical_json(result['source'])).hexdigest()})
        self.assertEqual(set(self.gateway.cancelled), {row['query_id'] for row in chunks})

    def test_invalid_population_is_rejected_before_any_submission(self):
        for ids in [[], self.ids(1001), self.ids(50) + ['candidate-0000'],
                    [''], [True], ['x' * 257], 'candidate']:
            with self.subTest(ids=str(ids)[:40]), self.assertRaises(ValueError):
                self.service._scout_materialization_chunk(ids)
        self.assertEqual(self.gateway.calls, [])

    def test_changed_source_never_returns_a_mixed_population(self):
        def change(result, batch):
            if batch == 2:
                result['source']['nested']['release'] = 8
        self.gateway.modify = change
        with self.assertRaisesRegex(ValueError, 'source changed'):
            with owned_catalog_scope(self.service):
                self.service._scout_materialization_chunk(self.ids(101))
        self.assertEqual(len(self.gateway.calls), 2)
        self.assertEqual(len(self.gateway.cancelled), 2)

    def test_missing_reordered_and_duplicate_details_still_fail_closed(self):
        for mode in ['missing', 'reordered', 'duplicate']:
            def change(result, batch):
                if batch != 2:
                    return
                if mode == 'missing':
                    result['rows'].pop()
                elif mode == 'reordered':
                    result['rows'].reverse()
                else:
                    result['rows'][1] = deepcopy(result['rows'][0])
            self.gateway = Gateway(change)
            self.service._scout_request_impl = self.gateway
            with self.subTest(mode=mode), self.assertRaisesRegex(ValueError, 'reordered catalog IDs'):
                self.service._scout_materialization_chunk(self.ids(101))
            self.assertEqual(len(self.gateway.calls), 2)

    def test_supply_and_metadata_guards_remain_active_in_later_batches(self):
        for mode in ['offers', 'metadata']:
            def change(result, batch):
                if batch != 2:
                    return
                row = result['rows'][0]
                if mode == 'offers':
                    row['offers_status'] = 'source_unavailable'
                else:
                    row['admet'] = [{'model_manifest': {'fixture': True},
                                     'model_manifest_sha256': '0' * 64}]
            self.gateway = Gateway(change)
            self.service._scout_request_impl = self.gateway
            with self.subTest(mode=mode), self.assertRaises(ValueError):
                self.service._scout_materialization_chunk(self.ids(101))
            self.assertEqual(len(self.gateway.calls), 2)

    def test_optional_offers_contract_is_forwarded_to_every_batch(self):
        def change(result, batch):
            for row in result['rows']:
                row['offers_status'] = 'source_unavailable'
        self.gateway.modify = change
        result = self.service._scout_materialization_chunk(self.ids(101), require_offers=False)
        self.assertEqual(len(result['rows']), 101)
        self.assertTrue(all(row['offers_status'] == 'source_unavailable' for row in result['rows']))

    def test_cancel_at_batch_boundary_prevents_another_post(self):
        stop_calls = [0]
        def stop():
            stop_calls[0] += 1
            # First pre-submit boundary, then first response, then next boundary.
            if stop_calls[0] == 3:
                raise Cancelled('synthetic user cancellation')
        with self.assertRaises(Cancelled):
            with owned_catalog_scope(self.service, stop=stop):
                self.service._scout_materialization_chunk(self.ids(1000))
        self.assertEqual(len(self.gateway.calls), 1)
        self.assertEqual(self.gateway.cancelled, ['0' * 31 + '1'])

    def test_uncertain_later_post_is_never_repeated_or_claimed_clean(self):
        def fail(result, batch):
            if batch == 2:
                raise _CatalogTransportError('synthetic response unavailable')
        self.gateway.modify = fail
        with self.assertRaisesRegex(OwnedProcessCleanupError, 'no verifiable query receipt'):
            with owned_catalog_scope(self.service):
                self.service._scout_materialization_chunk(self.ids(1000))
        self.assertEqual(len(self.gateway.calls), 2)
        self.assertEqual(self.gateway.cancelled, ['0' * 31 + '1'])

    def test_later_size_or_integrity_failure_has_no_retry_or_partial_return(self):
        def fail(result, batch):
            if batch == 2:
                raise ValueError('Global catalog response exceeds its bounded transfer limit')
        self.gateway.modify = fail
        with self.assertRaisesRegex(ValueError, 'bounded transfer limit'):
            self.service._scout_materialization_chunk(self.ids(101))
        self.assertEqual(len(self.gateway.calls), 2)

    def test_each_batch_keeps_the_existing_120_second_deadline(self):
        clock, calls = [0.0], []
        def fetch(path, **kwargs):
            calls.append(path)
            return {'status': 'running', 'query_id': 'a' * 32}
        self.service._scout_request = fetch
        with patch('skinpipe.platform.service.time.monotonic', side_effect=lambda: clock[0]), \
             patch('skinpipe.platform.service.time.sleep', side_effect=lambda _: clock.__setitem__(0, 121.0)):
            with self.assertRaisesRegex(ValueError, '120 second budget'):
                self.service._scout_materialization_chunk(self.ids(101))
        self.assertEqual(calls, ['/materialize', '/queries/' + 'a' * 32])

    def test_discovery_1001_retains_global_top_k_and_original_order(self):
        discovery = DiscoveryRuns.__new__(DiscoveryRuns)
        discovery.service = self.service
        ids = self.ids(1001)
        index = {'rows': [{'compound_id': cid, 'tanimoto_2d': 0.5} for cid in ids]}
        requested = []
        def retrieve(request, references, stop, progress):
            requested.append(request['catalog_limit'])
            return deepcopy(index)
        discovery._reference_index = retrieve
        request = {'input': {'kind': 'compound'}, 'catalog_limit': 1001,
                   'filters': freeze_filters({})}
        rows, receipt = discovery._catalog_rows(request, {'target_id': 'fixture'},
            lambda: None, lambda _: None, {'computational_seeds': [{'smiles': 'CCO'}]})
        self.assertEqual(requested, [1001])
        self.assertEqual([row['candidate_id'] for row in rows], ids)
        self.assertEqual(receipt['similarity_index'], index)
        self.assertEqual(receipt['completeness'], 'bounded_retrieval')
        self.assertEqual([len(body['compound_ids']) for body, _ in self.gateway.calls], [50] * 20 + [1])

    def test_text_query_keeps_its_existing_offset_and_population_contract(self):
        discovery = DiscoveryRuns.__new__(DiscoveryRuns)
        discovery.service = self.service
        ids, queries = self.ids(1001), []
        def query(body, *, owner):
            queries.append(deepcopy(body))
            return {'status': 'complete', 'rows': [{'compound_id': cid}
                for cid in ids[body['offset']:body['offset'] + body['limit']]]}
        self.service._scout_query = query
        rows, receipt = discovery._catalog_rows(
            {'catalog_limit': 1001, 'catalog_query': 'fixture', 'filters': freeze_filters({}),
             'admet_search_policy': {'sha256': 'c' * 64}},
            {'target_id': 'fixture'}, lambda: None, lambda _: None)
        self.assertEqual([(q['offset'], q['limit']) for q in queries], [(0, 1000), (1000, 1)])
        self.assertEqual([row['candidate_id'] for row in rows], ids)
        self.assertEqual(receipt['retrieval_policy'], 'deterministic_catalog_window')

    def test_frozen_snapshot_population_is_not_subdivided_or_truncated(self):
        discovery = DiscoveryRuns.__new__(DiscoveryRuns)
        discovery.service = self.service
        manifest = {'fixture': 'immutable'}
        saved = [{'candidate_id': cid} for cid in self.ids(1001)]
        self.service.catalog = SimpleNamespace(manifest=lambda _: manifest,
            iter_candidates=lambda *args, **kwargs: iter(saved))
        rows, receipt = discovery._catalog_rows(
            {'snapshot_id': 'saved', 'source_manifest_sha256': digest(manifest),
             'source_target_id': 'fixture', 'catalog_limit': 50},
            {}, lambda: None, lambda _: None)
        self.assertEqual(rows, saved)
        self.assertEqual(receipt['scope'], 'frozen source population')
        self.assertEqual(self.gateway.calls, [])


if __name__ == '__main__':
    unittest.main()
