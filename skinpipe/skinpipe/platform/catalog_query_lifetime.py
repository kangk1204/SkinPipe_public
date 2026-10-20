"""Cooperative cancellation confined to one private catalog Future/thread."""
from contextlib import contextmanager
from contextvars import ContextVar
import threading


class CatalogQueryCancelled(BaseException):
    """Control flow, not a missing evidence component recoverable by a caller."""


_stop = ContextVar("skinpipe_catalog_query_stop", default=None)


def query_checkpoint():
    event = _stop.get()
    if event is not None and event.is_set():
        raise CatalogQueryCancelled("Owned catalog query cancellation requested")


@contextmanager
def query_scope(event):
    token = _stop.set(event)
    try:
        query_checkpoint()
        yield
        query_checkpoint()
    finally:
        _stop.reset(token)


@contextmanager
def interruptible_query(connection):
    """Watch only this query's connection, joining before its owner closes it."""
    query_checkpoint()
    event = _stop.get()
    if event is None:
        yield
        return
    finished = threading.Event()

    def watch():
        while not finished.wait(0.05):
            if event.is_set():
                try:
                    connection.interrupt()
                except Exception:
                    pass  # A failed interrupt never acknowledges completion.
                return

    watcher = threading.Thread(target=watch, name="skinpipe-catalog-query-stop", daemon=True)
    watcher.start()
    try:
        yield
    finally:
        finished.set()
        watcher.join()
