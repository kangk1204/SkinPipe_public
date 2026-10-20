"""Select one exact saved parent while hashing and validating the entire JSON.

Only one undecided or matching candidate row is assembled at a time. A row
with a different ID is still validated and hashed through its end. The full
exploratory pool takes precedence even when it is empty or occurs after the
ranked pool. No compact row substitutes for the original scientific row.
"""
from __future__ import annotations

import hashlib
import gzip
import math
from contextlib import contextmanager
from pathlib import Path

import ijson
from ijson.common import ObjectBuilder

MAX_SOURCE_BYTES = 8 * 1024**3
MAX_RECORD_BYTES = 16 * 1024**2
MAX_DEPTH = 100


class _HashReader:
    def __init__(self, file):
        self.file, self.digest = file, hashlib.sha256()
        self.length = self.last_event = 0

    def read(self, size):
        if size < 0:
            raise ValueError("saved source reads must be bounded")
        raw = self.file.read(min(size, 65536))
        self.length += len(raw)
        if self.length > MAX_SOURCE_BYTES:
            raise ValueError("saved result exceeds the 8 GiB input limit")
        # Bound a token before the parser finishes allocating a huge string.
        if self.length - self.last_event > MAX_RECORD_BYTES + 65536:
            raise ValueError("saved result contains an oversized JSON token")
        self.digest.update(raw)
        return raw


def fingerprint(path):
    stat = Path(path).stat()
    return (stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)


def open_source(path):
    return gzip.open(path, "rb") if str(path).endswith(".gz") else Path(path).open("rb")


def source_sha(path):
    digest = hashlib.sha256()
    length = 0
    with open_source(path) as file:
        for raw in iter(lambda: file.read(1 << 20), b""):
            length += len(raw)
            if length > MAX_SOURCE_BYTES:
                raise ValueError("saved result exceeds the 8 GiB input limit")
            digest.update(raw)
    return digest.hexdigest()


@contextmanager
def _parse_events(reader):
    """Own the parser pipeline, including refusal before the JSON is complete.

    ijson's Python pull adapter closes its pipeline on GeneratorExit but not
    when the reader raises a size-limit exception. Its later finalizer then
    raises an unraisable incomplete-JSON error. Explicit push/close keeps that
    secondary error on the ordinary exception path, with bounded chunks and
    event storage. Primary parse errors and successful-input close errors are
    never suppressed.
    """
    events = ijson.sendable_list()
    parser = ijson.basic_parse_coro(events, use_float=True)
    closed = False

    def close():
        nonlocal closed
        if not closed:
            closed = True
            parser.close()

    def iterate():
        while chunk := reader.read(65536):
            parser.send(chunk)
            yield from events
            events.clear()
        close()  # Validate EOF before accepting the entire source.
        yield from events
        events.clear()

    try:
        yield iterate()
    except BaseException:
        try:
            close()
        except (ijson.JSONError, OverflowError):
            # A reader/consumer already refused this input. Closing an
            # unfinished document must not replace that original refusal.
            pass
        raise
    else:
        close()


def select_parent(path: Path, parent_id: str) -> dict:
    before = fingerprint(path)
    if before[2] > MAX_SOURCE_BYTES:
        raise ValueError("saved result exceeds the 8 GiB input limit")
    pools = {"candidates": [], "discovery_candidates": []}
    present, counts, matches, stack = set(), dict.fromkeys(pools, 0), dict.fromkeys(pools, 0), []
    marker = triage = None
    capture = None
    with open_source(path) as file:
        reader = _HashReader(file)
        try:
            with _parse_events(reader) as events:
                for event, value in events:
                    reader.last_event = reader.length
                    if event == "map_key":
                        frame = stack[-1]
                        if value in frame["keys"]:
                            raise ValueError("duplicate saved-result object key")
                        frame["keys"].add(value)
                        frame["key"] = value
                        path_key = frame["path"]
                    elif event in {"end_map", "end_array"}:
                        path_key = stack[-1]["path"]
                    else:
                        if not stack:
                            path_key = ()
                            if event != "start_map":
                                raise ValueError("saved result must be a JSON object")
                        else:
                            frame = stack[-1]
                            key = frame["key"] if frame["type"] == "map" else frame["index"]
                            path_key = (*frame["path"], key)
                            if frame["type"] == "array":
                                frame["index"] += 1
                        if len(path_key) == 1 and path_key[0] in pools:
                            if event != "start_array":
                                raise ValueError("saved result candidate pool must be a list of objects")
                            present.add(path_key[0])
                        if len(path_key) == 1 and path_key[0] == "marker":
                            marker = value if event == "string" else None
                        is_row = len(path_key) == 2 and path_key[0] in pools and isinstance(path_key[1], int)
                        is_triage = path_key == ("triage",)
                        if is_row or is_triage:
                            if event != "start_map":
                                raise ValueError("saved result candidate and triage records must be objects")
                            capture = {"path": path_key, "builder": ObjectBuilder(), "bytes": 0,
                                       "id_path": (*path_key, "ligand_id") if is_row else None}
                            if is_row:
                                counts[path_key[0]] += 1
                    if event == "number" and isinstance(value, float) and not math.isfinite(value):
                        raise ValueError("nonfinite saved-result value")
                    if capture is not None:
                        capture["bytes"] += len(value.encode("utf-8")) + 8 if isinstance(value, str) else 16
                        if capture["bytes"] > MAX_RECORD_BYTES:
                            raise ValueError("saved result candidate or triage record exceeds the 16 MiB limit")
                        if event == "string" and path_key == capture["id_path"] and value != parent_id:
                            # Only the row's own string ID can rule it out.
                            # Keep byte/depth/key/number checks for the entire
                            # row, but release its now-unneeded object tree.
                            capture["builder"] = None
                        if capture["builder"] is not None:
                            capture["builder"].event(event, value)
                        if event == "end_map" and path_key == capture["path"]:
                            record = capture["builder"].value if capture["builder"] is not None else None
                            if path_key == ("triage",):
                                triage = record
                            elif record is not None and record.get("ligand_id") == parent_id:
                                pool = pools[path_key[0]]
                                matches[path_key[0]] += 1
                                if not pool:
                                    pool.append(record)
                            capture = None
                    if event in {"start_map", "start_array"}:
                        if len(stack) >= MAX_DEPTH:
                            raise ValueError("saved result nesting exceeds the input limit")
                        stack.append({"path": path_key, "type": "map" if event == "start_map" else "array",
                                      "key": None, "keys": set(), "index": 0})
                    elif event in {"end_map", "end_array"}:
                        stack.pop()
        except (ijson.JSONError, OverflowError) as error:
            raise ValueError("invalid or nonfinite saved-result JSON") from error
    if fingerprint(path) != before:
        raise ValueError("saved source changed during parent selection")
    chosen = "discovery_candidates" if "discovery_candidates" in present else "candidates"
    if matches[chosen] != 1:
        raise ValueError("parent ID must identify exactly one candidate in the saved pool")
    node = (triage or {}).get("node")
    if not isinstance(marker, str) or not marker.strip() or not isinstance(node, str) or not node.strip():
        raise ValueError("saved result must name a marker and screening node")
    return {"parent": pools[chosen][0], "marker": marker, "node": node,
            "source_sha256": reader.digest.hexdigest(), "source_bytes": reader.length,
            "source_encoding": "gzip" if str(path).endswith(".gz") else "identity",
            "pool": chosen, "pool_rows": counts[chosen]}
