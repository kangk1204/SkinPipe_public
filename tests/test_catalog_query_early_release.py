"""Completed catalog payloads do not retain query lifetimes through chemistry."""
from concurrent.futures import ThreadPoolExecutor
from copy import deepcopy
import threading
import time
from types import SimpleNamespace
import unittest

from skinpipe.platform.catalog_cancellation import (
    catalog_checkpoint, note_catalog_query, note_catalog_transport_failure, owned_catalog_scope,
)
from skinpipe.platform.engines import OwnedProcessCleanupError
from skinpipe.platform.queue import Cancelled
from skinpipe.platform.scout_gateway import Gateway


QUERY = 'a' * 32


def released(query=QUERY):
    return {'query_id': query, 'status': 'cancelled', 'resources_released': True}


class CatalogQueryEarlyReleaseTests(unittest.TestCase):
    def test_complete_payload_released_before_next_step_without_mutation(self):
        calls = []
        service = SimpleNamespace(_scout_request=lambda path, body=None:
                                  calls.append((path, body)) or released())
        payload = {'query_id': QUERY, 'status': 'complete', 'rows': [{'value': 7}], 'source': {'sha': 'fixture'}}
        before = deepcopy(payload)
        with owned_catalog_scope(service):
            note_catalog_query('/materialize', payload)
            self.assertEqual(calls, [('/queries/' + QUERY + '/cancel', {})])
            self.assertEqual(payload, before)
        self.assertEqual(len(calls), 1)

    def test_actual_gateway_old_completed_query_cannot_poison_scope_exit(self):
        gateway = Gateway.__new__(Gateway)
        gateway.executor = ThreadPoolExecutor(max_workers=1)
        self.addCleanup(gateway.executor.shutdown, wait=True)
        gateway.jobs, gateway.query_receipts, gateway.query_stop_events = {}, {}, {}
        gateway.query_release_requested, gateway.released_queries = set(), {}
        gateway.lock = threading.Lock()
        gateway.scout = SimpleNamespace(status=lambda: {'synthetic_fixture': True})
        calls = []
        def call(path, body=None):
            calls.append((path, body))
            return gateway.cancel_query(path.split('/')[2]) if body is not None else gateway.result(path.split('/')[2])
        with owned_catalog_scope(SimpleNamespace(_scout_request=call)):
            first = gateway._submit(lambda: [{'fixture': 'first'}])
            query = first['query_id']
            note_catalog_query('/materialize', first)
            self.assertNotIn(query, gateway.jobs)
            self.assertEqual(gateway.result(query), released(query))
            gateway.released_queries[query] = time.monotonic() - 301
            second = gateway._submit(lambda: [{'fixture': 'second'}])
            note_catalog_query('/materialize', second)
            self.assertNotIn(query, gateway.released_queries)
            self.assertEqual(first['rows'], [{'fixture': 'first'}])
        self.assertEqual(len(calls), 2)
        self.assertEqual(gateway.jobs, {})

    def test_running_query_remains_owned_until_final_cleanup(self):
        calls = []
        service = SimpleNamespace(_scout_request=lambda path, body=None:
                                  calls.append((path, body)) or released())
        with owned_catalog_scope(service):
            note_catalog_query('/similarity', {'query_id': QUERY, 'status': 'running'})
            self.assertEqual(calls, [])
        self.assertEqual(calls, [('/queries/' + QUERY + '/cancel', {})])

    def test_incomplete_release_is_polled_before_forgetting_query(self):
        calls = []
        def call(path, body=None):
            calls.append((path, body))
            return {'query_id': QUERY, 'status': 'stopping', 'resources_released': False} if body is not None else released()
        with owned_catalog_scope(SimpleNamespace(_scout_request=call)):
            note_catalog_query('/exact-structures', {'query_id': QUERY, 'status': 'complete'})
            self.assertEqual(calls, [('/queries/' + QUERY + '/cancel', {}), ('/queries/' + QUERY, None)])
        self.assertEqual(len(calls), 2)

    def test_missing_wrong_id_false_and_wrong_status_never_certify_release(self):
        malformed = [
            {'status': 'cancelled', 'resources_released': True},
            released('b' * 32), {**released(), 'resources_released': False},
            {**released(), 'resources_released': 1}, {**released(), 'status': 'complete'},
        ]
        for receipt in malformed:
            with self.subTest(receipt=receipt):
                service = SimpleNamespace(_scout_request=lambda *args, **kwargs: receipt)
                with self.assertRaises(OwnedProcessCleanupError):
                    with owned_catalog_scope(service, timeout=0):
                        note_catalog_query('/query', {'query_id': QUERY, 'status': 'complete'})
        catalog_checkpoint()

    def test_expired_unknown_query_does_not_become_cleanup_success(self):
        def missing(*args, **kwargs):
            raise KeyError('Query expired; cleanup cannot be confirmed')
        with self.assertRaises(OwnedProcessCleanupError):
            with owned_catalog_scope(SimpleNamespace(_scout_request=missing)):
                note_catalog_query('/materialize', {'query_id': QUERY, 'status': 'complete'})

    def test_unknown_submission_remains_held_after_known_query_released(self):
        calls = []
        service = SimpleNamespace(_scout_request=lambda path, body=None: calls.append(path) or released())
        with self.assertRaisesRegex(OwnedProcessCleanupError, 'no verifiable query receipt'):
            with owned_catalog_scope(service):
                note_catalog_query('/materialize', {'query_id': QUERY, 'status': 'complete'})
                note_catalog_transport_failure('/materialize')
        self.assertEqual(calls, ['/queries/' + QUERY + '/cancel'])

    def test_user_stop_keeps_original_reason_after_completed_query_release(self):
        calls = []
        def stop():
            raise Cancelled('original user stop')
        service = SimpleNamespace(_scout_request=lambda path, body=None: calls.append(path) or released())
        with self.assertRaisesRegex(Cancelled, 'original user stop'):
            with owned_catalog_scope(service, stop=stop):
                note_catalog_query('/diversity', {'query_id': QUERY, 'status': 'complete'})
        self.assertEqual(calls, ['/queries/' + QUERY + '/cancel'])


if __name__ == '__main__':
    unittest.main()
