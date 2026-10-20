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
    tracked = {"queries": set(), "uncertain_submission": False, "releasing": False, "stop": stop}
    token = _scope.set(tracked)
    try:
        yield
    finally:
        tracked["releasing"] = True
        try:
            # Release known queries in parallel logically: submit every stop
            # intent first, then poll only those whose actual future still runs.
            pending = set()
            deadline = time.monotonic() + timeout
            for query_id in sorted(tracked["queries"]):
                try:
                    receipt = service._scout_request("/queries/" + query_id + "/cancel", body={})
                except Exception as error:
                    raise OwnedProcessCleanupError("Owned catalog query cleanup is unavailable") from error
                if receipt.get("status") != "cancelled" or receipt.get("resources_released") is not True:
                    pending.add(query_id)
            while pending:
                if time.monotonic() >= deadline:
                    raise OwnedProcessCleanupError("Owned catalog query is still stopping")
                time.sleep(0.1)
                for query_id in list(pending):
                    try:
                        receipt = service._scout_request("/queries/" + query_id)
                    except Exception as error:
                        raise OwnedProcessCleanupError("Owned catalog query release could not be verified") from error
                    if receipt.get("status") == "cancelled" and receipt.get("resources_released") is True:
                        pending.remove(query_id)
            if tracked["uncertain_submission"]:
                raise OwnedProcessCleanupError("Catalog submission has no verifiable query receipt")
        finally:
            _scope.reset(token)
