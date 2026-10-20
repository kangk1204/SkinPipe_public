"""Private read-only Scout gateway; bounded asynchronous catalog queries.

Run only on an explicitly registered Tailscale address. No donor writes and no
arbitrary SQL/path endpoints. Expensive queries have a visible pending state.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, TimeoutError
import hmac
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import json
from pathlib import Path
import re
import threading
import time
from urllib.parse import unquote, urlparse
import uuid

from .scout import ScoutCatalog
from .identity import digest
from .catalog_query_lifetime import CatalogQueryCancelled, query_scope


class Gateway:
    def __init__(self, config):
        self.scout = ScoutCatalog(config["scout_root"], query_timeout=30)
        self.token = config["token"]
        if not isinstance(self.token, str) or len(self.token) < 32:
            raise ValueError("Private gateway token required")
        self.executor = ThreadPoolExecutor(max_workers=2)
        self.jobs = {}
        self.query_receipts = {}
        self.query_release_requested = set()
        self.query_stop_events = {}
        self.released_queries = {}
        self.lock = threading.Lock()
        self.admission_errors = []
        if config.get("published_prediction_cache"):
            try:
                from .scout_predictions import PublishedPredictionCache

                self.scout.published_cache = PublishedPredictionCache(
                    config["published_prediction_cache"]
                )
            except Exception as error:
                self.admission_errors.append(
                    {
                        "component": "published_admet_cache",
                        "error": type(error).__name__,
                    }
                )
        self.similarity = None
        self.filtered_similarity = {}
        if config.get("similarity_index"):
            try:
                from .catalog_index import CatalogSimilarityIndex

                self.similarity = CatalogSimilarityIndex(config["similarity_index"])
            except Exception as error:
                self.admission_errors.append(
                    {"component": "catalog_similarity", "error": type(error).__name__}
                )
        if config.get("admet_search_gate"):
            try:
                from .admet_search import ADMETSearchGate

                self.scout.search_gate = ADMETSearchGate(config["admet_search_gate"])
            except Exception as error:
                self.admission_errors.append(
                    {"component": "admet_search_gate", "error": type(error).__name__}
                )
        for key, root in config.get("filtered_similarity_indices", {}).items():
            if key not in {"risk", "purchase", "risk_purchase"}:
                raise ValueError("Unknown catalog index filter")
            try:
                from .catalog_index import CatalogSimilarityIndex

                self.filtered_similarity[key] = CatalogSimilarityIndex(root)
            except Exception as error:
                self.admission_errors.append(
                    {
                        "component": "catalog_similarity_" + key,
                        "error": type(error).__name__,
                    }
                )
        for item in config.get("prediction_partitions", []):
            try:
                self.scout.admit_prediction_partition(
                    item["file"],
                    model_manifest=item["model_manifest"],
                    endpoint_definitions=item["endpoint_definitions"],
                    model_manifest_sha256=item.get("model_manifest_sha256"),
                )
            except Exception as exc:
                self.admission_errors.append(
                    {"file": Path(item["file"]).name, "error": type(exc).__name__}
                )

    def query(self, body, *, materialize=False):
        if materialize:
            if not isinstance(body, dict) or set(body) != {"compound_ids"}:
                raise ValueError("Only compound_ids accepted for materialization")
            ids = body["compound_ids"]
            if (
                not isinstance(ids, list)
                or not 1 <= len(ids) <= 1000
                or any(not isinstance(v, str) or not v or len(v) > 256 for v in ids)
                or len(set(ids)) != len(ids)
            ):
                raise ValueError(
                    "Materialization requires 1..1000 unique bounded compound IDs"
                )
            return self._submit(self.scout.materialize, ids)
        if not isinstance(body, dict) or set(body) - {
            "q",
            "limit",
            "offset",
            "purchasable_only",
            "hide_severe_admet",
            "admet_policy_sha256",
        }:
            raise ValueError("Unsupported catalog query")
        self._verify_frozen_policy(body)
        return self._submit(
            self.scout.query,
            **{k: v for k, v in body.items() if k != "admet_policy_sha256"},
        )

    def _verify_frozen_policy(self, body):
        expected = body.get("admet_policy_sha256")
        if expected is None:
            return  # Interactive queries use the current, explicitly reported policy.
        if (
            not isinstance(expected, str)
            or len(expected) != 64
            or any(c not in "0123456789abcdef" for c in expected)
        ):
            raise ValueError("Invalid frozen ADMET policy hash")
        if body.get("hide_severe_admet", bool(self.scout.search_gate)):
            actual = (
                self.scout.search_gate.manifest["policy"]["sha256"]
                if self.scout.search_gate
                else None
            )
            if actual != expected:
                raise ValueError("Frozen ADMET policy is unavailable on this gateway")

    def similarity_query(self, body):
        if not self.similarity:
            raise ValueError("Whole-catalog similarity index not configured")
        if not isinstance(body, dict) or set(body) - {
            "reference_smiles",
            "limit",
            "purchasable_only",
            "hide_severe_admet",
            "admet_policy_sha256",
        }:
            raise ValueError("Unsupported similarity query fields")
        self._verify_frozen_policy(body)
        purchase = body.get("purchasable_only", False)
        risk = body.get("hide_severe_admet", bool(self.scout.search_gate))
        if type(purchase) is not bool or type(risk) is not bool:
            raise ValueError("Similarity switches must be boolean")
        index = self.similarity
        if purchase or risk:
            key = (
                "risk_purchase"
                if purchase and risk
                else "purchase"
                if purchase
                else "risk"
            )
            index = self.filtered_similarity.get(key)
            if index is None:
                raise ValueError("Requested filtered similarity index is not ready")
            if index.manifest.get("search_filters") != {
                "hide_severe_admet": risk, "purchasable_only": purchase,
            }:
                raise ValueError("Similarity index does not match requested filters")
        return self._submit(
            self._verified_similarity_rows, index, risk,
            {k: v for k, v in body.items() if k in {"reference_smiles", "limit"}},
        )

    def diversity_query(self, body):
        """Seed-independent library support under the same sealed filter gates."""
        from .catalog_diversity import ALGORITHM, validate_request
        if not isinstance(body, dict) or set(body) - {
            "limit", "seed", "algorithm", "purchasable_only", "hide_severe_admet", "admet_policy_sha256",
        }:
            raise ValueError("Unsupported diversity query fields")
        query = {"limit": body.get("limit", 1000), "seed": body.get("seed", 0),
                 "algorithm": body.get("algorithm", ALGORITHM)}
        validate_request(**query)
        self._verify_frozen_policy(body)
        purchase = body.get("purchasable_only", False)
        risk = body.get("hide_severe_admet", bool(self.scout.search_gate))
        if type(purchase) is not bool or type(risk) is not bool:
            raise ValueError("Diversity switches must be boolean")
        if not self.similarity:
            raise ValueError("Whole-catalog diversity index not configured")
        index = self.similarity
        if purchase or risk:
            key = "risk_purchase" if purchase and risk else "purchase" if purchase else "risk"
            index = self.filtered_similarity.get(key)
            if index is None:
                raise ValueError("Requested filtered diversity index is not ready")
            if index.manifest.get("search_filters") != {
                "hide_severe_admet": risk, "purchasable_only": purchase,
            }:
                raise ValueError("Diversity index does not match requested filters")
            if index.manifest.get("full_index_manifest_sha256") != digest(self.similarity.manifest):
                raise ValueError("Filtered diversity index belongs to a different full index")
        elif index.manifest.get("search_filters", {"hide_severe_admet": False, "purchasable_only": False}) != {
            "hide_severe_admet": False, "purchasable_only": False,
        }:
            raise ValueError("Diversity index does not match requested filters")
        return self._submit(self._verified_diversity_rows, index, risk, query,
                            _receipt={"operation": "diversity", "request_sha256": digest(body)})

    def _verified_diversity_rows(self, index, risk, query):
        from .catalog_query_lifetime import query_checkpoint
        query_checkpoint()
        expected = index.manifest["source_structure_fingerprint"]
        if self.scout.structure_fingerprint() != expected:
            raise ValueError("Catalog changed after diversity indexing; rebuild required")
        if risk and (not self.scout.search_gate or index.manifest.get("search_gate_manifest_sha256")
                     != digest(self.scout.search_gate.manifest)):
            raise ValueError("ADMET eligibility and diversity index differ")
        result = index.diversity(**query)
        query_checkpoint()
        if self.scout.structure_fingerprint() != expected:
            raise ValueError("Catalog changed during diversity selection")
        if risk and (not self.scout.search_gate or index.manifest.get("search_gate_manifest_sha256")
                     != digest(self.scout.search_gate.manifest)):
            raise ValueError("ADMET eligibility changed during diversity selection")
        return result

    def exact_structures(self, body):
        if not isinstance(body, dict) or set(body) != {"smiles"}:
            raise ValueError("Exact lookup accepts only a smiles list")
        values = body["smiles"]
        if (not isinstance(values, list) or not 1 <= len(values) <= 1000
                or any(not isinstance(value, str) or not 1 <= len(value) <= 4096
                       or any(ord(char) < 32 for char in value) for value in values)):
            raise ValueError("Exact lookup requires 1..1000 bounded structure strings")
        return self._submit(self._verified_exact_rows, values, _receipt={
            "operation": "exact_structures", "request_sha256": digest({"smiles": values})})

    def _exact_source_binding(self):
        source = self.scout.status()
        if (source.get("structure_admission") != "reported_verified"
                or not isinstance(source.get("report_sha256"), str)
                or type(source.get("compound_count")) is not int):
            raise ValueError("Verified catalog population is required for exact lookup")
        return {"schema": "scout-exact-source-v1", "report_sha256": source["report_sha256"],
                "compound_count": source["compound_count"], "database_identity": self.scout.database_identity()}

    def _verified_exact_rows(self, values):
        from .local_catalog import lookup_exact_structures
        before = self._exact_source_binding()
        rows = lookup_exact_structures(self.scout, values)
        if self._exact_source_binding() != before:
            raise ValueError("Catalog changed during exact structure lookup")
        return {"rows": rows, "source_binding": before}

    def _verified_similarity_rows(self, index, risk, query):
        # A cold structure seal hashes the DB and all external partitions.
        # Perform it inside the bounded asynchronous job, before lookup. Doing
        # it in the HTTP handler can exceed the client's 25-second timeout.
        if (
            self.scout.structure_fingerprint()
            != index.manifest["source_structure_fingerprint"]
        ):
            raise ValueError(
                "Catalog changed after similarity indexing; rebuild required"
            )
        if risk and (
            not self.scout.search_gate
            or index.manifest.get("search_gate_manifest_sha256")
            != digest(self.scout.search_gate.manifest)
        ):
            raise ValueError("ADMET eligibility and similarity index differ")
        return index.query(**query)

    def _submit(self, function, *args, _receipt=None, **kwargs):
        with self.lock:
            now = time.monotonic()
            self._prune_released_locked(now)
            self.jobs = {
                k: v
                for k, v in self.jobs.items()
                if k in getattr(self, "query_release_requested", set())
                or not v[0].done() or now - v[1] < 300
            }
            self.query_receipts = {k: v for k, v in getattr(self, "query_receipts", {}).items() if k in self.jobs}
            self.query_stop_events = {k: v for k, v in getattr(self, "query_stop_events", {}).items() if k in self.jobs}
            if sum(not f.done() for f, _ in self.jobs.values()) >= 4:
                raise ValueError("Catalog is busy; retry after active queries finish")
            query_id = uuid.uuid4().hex
            event = threading.Event()
            future = self.executor.submit(self._execute_query, event, function, args, kwargs)
            self.jobs[query_id] = (future, now)
            if not hasattr(self, "query_stop_events"):
                self.query_stop_events = {}
            self.query_stop_events[query_id] = event
            if _receipt is not None:
                self.query_receipts[query_id] = dict(_receipt)
        return self.result(query_id, wait=1.5)

    @staticmethod
    def _execute_query(event, function, args, kwargs):
        with query_scope(event):
            return function(*args, **kwargs)

    def _prune_released_locked(self, now=None):
        """Keep only small, bounded cleanup receipts, never query rows."""
        now = time.monotonic() if now is None else now
        released = getattr(self, "released_queries", {})
        self.released_queries = dict(
            sorted(((key, finished) for key, finished in released.items()
                    if now - finished < 300), key=lambda item: item[1])[-512:]
        )

    def _complete_release_locked(self, query_id, future):
        item = self.jobs.get(query_id)
        if (not item or item[0] is not future or not future.done()
                or query_id not in getattr(self, "query_release_requested", set())):
            return
        self.jobs.pop(query_id, None)
        getattr(self, "query_receipts", {}).pop(query_id, None)
        getattr(self, "query_stop_events", {}).pop(query_id, None)
        self.query_release_requested.discard(query_id)
        if not hasattr(self, "released_queries"):
            self.released_queries = {}
        self.released_queries[query_id] = time.monotonic()
        self._prune_released_locked()

    def _finish_release(self, query_id, future):
        # Future callbacks can execute inline when registered on a done Future.
        # Registration and Future.cancel() must therefore happen outside lock.
        with self.lock:
            self._complete_release_locked(query_id, future)

    def _release_receipt_locked(self, query_id):
        self._prune_released_locked()
        if query_id in getattr(self, "query_release_requested", set()):
            item = self.jobs.get(query_id)
            if item:
                self._complete_release_locked(query_id, item[0])
            if query_id in self.query_release_requested:
                return {"status": "stopping", "query_id": query_id,
                        "resources_released": False}
        if query_id in self.released_queries:
            return {"status": "cancelled", "query_id": query_id,
                    "resources_released": True}
        return None

    def cancel_query(self, query_id):
        """Release one caller-owned query without touching shared DB/index state.

        Running native work cannot be interrupted by Future.cancel(). It remains
        admitted until it finishes; only its done callback drops retained rows
        and receipts and acknowledges release. Queued work can be cancelled.
        """
        with self.lock:
            receipt = self._release_receipt_locked(query_id)
            if receipt:
                return receipt
            item = self.jobs.get(query_id)
            if item is None:
                raise KeyError("Query expired; cleanup cannot be confirmed")
            if not hasattr(self, "query_release_requested"):
                self.query_release_requested = set()
            self.query_release_requested.add(query_id)
            future = item[0]
            event = getattr(self, "query_stop_events", {}).get(query_id)
            if event is not None:
                event.set()
        future.cancel()
        future.add_done_callback(lambda done: self._finish_release(query_id, done))
        with self.lock:
            return self._release_receipt_locked(query_id)

    def result(self, query_id, *, wait=0):
        with self.lock:
            released = self._release_receipt_locked(query_id)
            if released:
                return released
            item = self.jobs.get(query_id)
            receipt = getattr(self, "query_receipts", {}).get(query_id, {})
        if not item:
            raise KeyError("Query expired; submit again")
        future, started = item
        try:
            rows = future.result(timeout=wait)
        except TimeoutError:
            if not future.done():
                with self.lock:
                    released = self._release_receipt_locked(query_id)
                    if released:
                        return released
                return {
                    "status": "running",
                    "query_id": query_id,
                    "stage": "catalog_query",
                    "elapsed_seconds": round(time.monotonic() - started, 3),
                    **receipt,
                }
            # A query can raise the futures timeout type too (an alias of the
            # built-in TimeoutError since Python 3.11). Retrieve again without
            # a wait: propagate its real failure or retrieve a raced value.
            try:
                rows = future.result()
            except (Exception, CatalogQueryCancelled):
                with self.lock:
                    released = self._release_receipt_locked(query_id)
                    if released:
                        return released
                raise
        except (Exception, CatalogQueryCancelled):
            with self.lock:
                released = self._release_receipt_locked(query_id)
                if released:
                    return released
            raise
        with self.lock:
            released = self._release_receipt_locked(query_id)
            if released:
                return released
        binding = {}
        if receipt.get("operation") == "exact_structures":
            binding = {"source_binding": rows["source_binding"]}
            rows = rows["rows"]
        return {
            "status": "complete",
            "query_id": query_id,
            "rows": rows,
            "elapsed_seconds": round(time.monotonic() - started, 3),
            "query_scope": "full_source_catalog",
            "source": self.scout.status(),
            **receipt,
            **binding,
        }


def handler_for(gateway):
    class Handler(BaseHTTPRequestHandler):
        def reply(self, data, status=200):
            raw = json.dumps(
                data, ensure_ascii=False, allow_nan=False, default=str
            ).encode()
            self.send_response(status)
            self.send_header("Content-Type", "application/json; charset=utf-8")
            self.send_header("Content-Length", str(len(raw)))
            self.send_header("Cache-Control", "no-store")
            self.end_headers()
            self.wfile.write(raw)

        def authenticated(self):
            if not hmac.compare_digest(
                self.headers.get("Authorization", "").encode("utf-8"),
                ("Bearer " + gateway.token).encode("utf-8"),
            ):
                self.reply({"error": "Unauthorized"}, 401)
                return False
            return True

        def do_GET(self):
            if not self.authenticated():
                return
            path = unquote(urlparse(self.path).path)
            try:
                if path == "/status":
                    return self.reply(
                        {
                            **gateway.scout.status(),
                            "admission_errors": gateway.admission_errors,
                        }
                    )
                if path.startswith("/queries/"):
                    return self.reply(gateway.result(path.split("/")[-1]))
                if path.startswith("/compound/"):
                    result = gateway.scout.get_compound(path.removeprefix("/compound/"))
                    if result is None:
                        raise KeyError("Unknown compound")
                    return self.reply(result)
                raise KeyError("Unknown route")
            except KeyError as exc:
                self.reply({"error": str(exc)}, 404)
            except (ValueError, TypeError) as exc:
                self.reply({"error": str(exc)}, 400)
            except Exception:
                self.reply(
                    {
                        "error": "Read-only catalog query failed; verified target snapshots remain available"
                    },
                    503,
                )

        def do_POST(self):
            if not self.authenticated():
                return
            try:
                length = int(self.headers.get("Content-Length", "0"))
                if (
                    not 0
                    < length
                    <= (
                        300000
                        if self.path == "/materialize"
                        else 1024 * 1024
                        if self.path == "/exact-structures"
                        else 25000
                        if self.path == "/similarity"
                        else 4096
                    )
                ):
                    raise ValueError("Query size exceeds limit")
                body = json.loads(self.rfile.read(length))
                path = unquote(urlparse(self.path).path)
                match = re.fullmatch(r"/queries/([0-9a-f]{32})/cancel", path)
                if match:
                    if body != {}:
                        raise ValueError("Query cancellation requires an empty object")
                    return self.reply(gateway.cancel_query(match[1]))
                if self.path == "/exact-structures":
                    return self.reply(gateway.exact_structures(body))
                if self.path == "/similarity":
                    return self.reply(gateway.similarity_query(body))
                if self.path == "/diversity":
                    return self.reply(gateway.diversity_query(body))
                if self.path not in {"/query", "/materialize"}:
                    raise KeyError("Unknown route")
                return self.reply(
                    gateway.query(body, materialize=self.path == "/materialize")
                )
            except KeyError as exc:
                self.reply({"error": str(exc)}, 404)
            except (ValueError, TypeError) as exc:
                self.reply({"error": str(exc)}, 400)
            except Exception:
                self.reply({"error": "Catalog query failed"}, 503)

        def log_message(self, *_):
            pass  # no private query/credential logging

    return Handler


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    path = Path(args.config)
    if path.stat().st_mode & 0o077:
        raise SystemExit("Configuration must have mode 0600")
    config = json.loads(path.read_text())
    gateway = Gateway(config)
    server = ThreadingHTTPServer(
        (config["bind"], config.get("port", 8944)), handler_for(gateway)
    )
    try:
        server.serve_forever()
    finally:
        server.server_close()
        gateway.executor.shutdown(wait=True, cancel_futures=True)


if __name__ == "__main__":
    main()
