"""Per-research catalog query lifetime, separate from shared index caches."""
from contextlib import contextmanager
from contextvars import ContextVar
import re
import time

from .engines import OwnedProcessCleanupError

_scope = ContextVar("skinpipe_owned_catalog_queries", default=None)
_SUBMIT = {"/query", "/materialize", "/similarity", "/diversity", "/exact-structures"}


def note_catalog_query(path, result):
    tracked = _scope.get()
    if tracked is None or tracked["releasing"]:
        return
    if path not in _SUBMIT and not re.fullmatch(r"/queries/[0-9a-f]{32}", path):
        return
    query_id = result.get("query_id") if isinstance(result, dict) else None
    if isinstance(query_id, str) and re.fullmatch(r"[0-9a-f]{32}", query_id):
        tracked["queries"].add(query_id)
        if result.get("status") == "complete":
            # The complete payload is already local. Do not retain its remote
            # Future throughout chemistry, past the gateway's result TTL.
            tracked["releasing"] = True
            try:
                _release_queries(tracked["service"], {query_id}, tracked["timeout"])
                tracked["queries"].discard(query_id)
            finally:
                tracked["releasing"] = False
    catalog_checkpoint()


def catalog_checkpoint():
    tracked = _scope.get()
    if tracked is not None and not tracked["releasing"]:
        tracked["stop"]()


def note_catalog_transport_failure(path):
    tracked = _scope.get()
    if tracked is not None and not tracked["releasing"] and path in _SUBMIT:
        # An uncertain POST may have admitted work even when no ID came back.
        # Never infer process/memory release from that missing receipt.
        tracked["uncertain_submission"] = True


@contextmanager
def owned_catalog_scope(service, *, timeout=40, stop=lambda: None):
    tracked = {"queries": set(), "uncertain_submission": False, "releasing": False,
               "stop": stop, "service": service, "timeout": timeout}
    token = _scope.set(tracked)
    try:
        yield
    finally:
        tracked["releasing"] = True
        try:
            _release_queries(service, tracked["queries"], timeout)
            if tracked["uncertain_submission"]:
                raise OwnedProcessCleanupError("Catalog submission has no verifiable query receipt")
        finally:
            _scope.reset(token)


def _release_queries(service, queries, timeout):
    # Request every known stop first, then poll only unfinished owned Futures.
    # Missing/expired queries are not evidence of process or memory release.
    def released(query_id, receipt):
        return (isinstance(receipt, dict) and receipt.get("query_id") == query_id
                and receipt.get("status") == "cancelled"
                and receipt.get("resources_released") is True)

    pending = set()
    failure = None
    deadline = time.monotonic() + timeout
    for query_id in sorted(queries):
        try:
            receipt = service._scout_request("/queries/" + query_id + "/cancel", body={})
        except Exception as error:
            if failure is None:
                failure = ("Owned catalog query cleanup is unavailable", error)
            continue
        if not released(query_id, receipt):
            pending.add(query_id)
    while pending:
        if time.monotonic() >= deadline:
            raise OwnedProcessCleanupError("Owned catalog query is still stopping") from (failure[1] if failure else None)
        time.sleep(0.1)
        for query_id in sorted(pending):
            try:
                receipt = service._scout_request("/queries/" + query_id)
            except Exception as error:
                if failure is None:
                    failure = ("Owned catalog query release could not be verified", error)
                pending.remove(query_id)
                continue
            if released(query_id, receipt):
                pending.remove(query_id)
    # A failed receipt still holds the research, but must not prevent release
    # of its other known Futures on the shared Gateway.
    if failure is not None:
        raise OwnedProcessCleanupError(failure[0]) from failure[1]
