"""Private retained catalog bindings, scoped to one research execution.

Scientific run manifests carry opaque binding IDs and component identities.
Only this operator-owned registry stores installed paths; credentials are never
copied here. A context variable selects an adapter without changing the active
catalog used by other requests or executions.
"""
from __future__ import annotations

from contextlib import contextmanager
from contextvars import ContextVar
from copy import deepcopy
import json
import glob
from pathlib import Path
import re
import threading

from .artifacts import atomic_bytes, canonical_json, sha256_file
from .identity import digest


def source_identity(local, scout):
    scout = scout or {}
    if local is None:
        return {"policy": "catalog-source-binding-v1", "kind": "remote" if scout else "not_configured",
                "endpoint_sha256": digest(scout.get("url")) if scout else None,
                "scope": "remote configuration only; actual source receipts frozen during acquisition" if scout else "No catalog configured at submission"}
    config = local.config
    components = {key: value for key, value in config.items()
                  if key in {"scout_root", "published_prediction_cache", "similarity_index", "admet_search_gate"} and value}
    components.update({"filtered_" + key: value for key, value in config.get("filtered_similarity_indices", {}).items()})
    locations = {name: Path(value).resolve() for name, value in components.items()}
    dependencies, dependency_status = {}, "verified"
    try:
        from .database_bundle import _component_path, _sql_paths
        from .scout import ScoutCatalog
        database = locations["scout_root"] / "data/skinscout.duckdb"
        if not database.is_file():
            raise ValueError("Frozen Scout database unavailable")
        reader = getattr(getattr(local, "gateway", None), "scout", None) or ScoutCatalog(locations["scout_root"])
        # Reuse the bounded Scout connection settings, so concurrent reads of
        # this database do not conflict through a second DuckDB configuration.
        views = reader._query("SELECT sql FROM duckdb_views() WHERE NOT internal ORDER BY schema_name,view_name")
        for row in views:
            for value in _sql_paths(row["sql"]):
                _component_path(value, locations)
                paths = sorted(glob.glob(value, recursive=True)) if glob.has_magic(value) else [value]
                if not paths or len(paths) > 200000:
                    raise ValueError("Parquet views require a nonempty bounded frozen file population")
                for path in paths:
                    parts = Path(_component_path(path, locations)).parts
                    dependencies.setdefault(parts[1], set()).add(str(Path(*parts[2:])))
    except Exception as error:
        # Source opening must remain optional for seed-only assessment. A
        # missing or unreadable dependency inventory is an explicit hold.
        dependency_status = "unavailable_" + type(error).__name__
    frozen = {}
    for name, location in sorted(components.items()):
        root = Path(location).resolve()
        names = ("data/skinscout.duckdb", "reports/audit.json", "reports/db_build.json") if name == "scout_root" else (
            "manifest.json", "predictions.duckdb") if name == "published_prediction_cache" else (
            "manifest.json", "catalog.duckdb") if name == "admet_search_gate" else (
            "manifest.json", "structures.h5", "lookup.duckdb")
        # Only admitted DuckDB view dependencies enter this inventory.
        # Unreferenced raw/quarantined downloads are not this run's source.
        names = sorted(set(names) | dependencies.get(name, set()))
        files = {relative: _file_identity(root / relative) for relative in names}
        frozen[name] = {"location_sha256": digest(str(root)), "files": files}
    # Explicit prediction partitions are also evidence-bearing sources. Small
    # endpoint/model definitions are hashed; large data files use immutable
    # installed-file state plus provider admission on opening.
    for index, partition in enumerate(config.get("prediction_partitions", [])):
        files = {}
        for key in ("file", "model_manifest", "endpoint_definitions"):
            value = partition.get(key)
            if isinstance(value, str):
                path = Path(value).resolve()
                files[key] = {"location_sha256": digest(str(path)), **_file_identity(path)}
            elif value is not None:
                files[key] = {"sha256": digest(value)}
        frozen["prediction_partition_" + str(index)] = {"files": files, "model_manifest_sha256": partition.get("model_manifest_sha256")}
    return {"policy": "catalog-source-binding-v1", "kind": "local", "components": frozen,
            "dependency_inventory_status": dependency_status,
            "scope": "Frozen component locations, manifests and file state; existing provider content-admission checks remain required"}


def _file_identity(path):
    try:
        stat = path.stat()
        facts = {"size_bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns,
                 "ctime_ns": stat.st_ctime_ns, "inode": stat.st_ino, "device": stat.st_dev}
        if path.suffix == ".json" and stat.st_size <= 32 * 1024**2:
            facts["sha256"] = sha256_file(path)
        return facts
    except OSError:
        return {"status": "unavailable"}


def _private_config(config):
    allowed = {"kind", "scout_root", "published_prediction_cache", "similarity_index",
               "filtered_similarity_indices", "admet_search_gate", "prediction_partitions"}
    if not isinstance(config, dict) or set(config) - allowed or config.get("kind") != "local":
        raise ValueError("Only trusted local catalog paths can be retained")
    result = deepcopy(config)
    for key in ("scout_root", "published_prediction_cache", "similarity_index", "admet_search_gate"):
        if result.get(key):
            result[key] = str(Path(result[key]).resolve())
    result["filtered_similarity_indices"] = {key: str(Path(value).resolve()) for key, value in result.get("filtered_similarity_indices", {}).items()}
    partitions = result.get("prediction_partitions", [])
    for partition in partitions:
        if set(partition) - {"file", "model_manifest", "endpoint_definitions", "model_manifest_sha256"}:
            raise ValueError("Only prediction paths and immutable scientific definitions can be retained")
        for key in ("file", "model_manifest", "endpoint_definitions"):
            if isinstance(partition.get(key), str):
                partition[key] = str(Path(partition[key]).resolve())
    return result


class CatalogSources:
    def __init__(self, service, *, adapter_factory=None):
        self.service = service
        self.root = service.root / "catalog-source-bindings"
        self.root.mkdir(mode=0o700, parents=True, exist_ok=True)
        self.root.chmod(0o700)
        self.scope = ContextVar("skinpipe_catalog_" + str(id(self)), default=None)
        self.lock = threading.RLock()
        self.adapters = {}
        if adapter_factory is None:
            from .local_catalog import LocalCatalog
            adapter_factory = LocalCatalog
        self.adapter_factory = adapter_factory

    def runtime(self):
        scoped = self.scope.get()
        if scoped is not None:
            return scoped["adapter"], scoped["scout"]
        return getattr(self.service, "local_catalog", None), getattr(self.service, "scout", {})

    def identity(self):
        scoped = self.scope.get()
        if scoped is not None and scoped.get("unavailable"):
            return {"policy": "catalog-source-binding-v1", "kind": "unavailable", "reason": scoped["unavailable"]}
        return source_identity(*self.runtime())

    def bind(self):
        local, scout = self.runtime()
        identity = source_identity(local, scout)
        record = {"schema": "skinpipe-private-catalog-source-v1", "identity": identity,
                  "config": _private_config(local.config) if local is not None else None}
        binding_id = digest(record)
        path = self.root / (binding_id + ".json")
        with self.lock:
            if not path.exists():
                atomic_bytes(path, canonical_json(record))
                path.chmod(0o600)
            if local is not None:
                self.adapters.setdefault(binding_id, local)
        return {"binding_id": binding_id, "kind": identity["kind"]}, identity

    @contextmanager
    def context(self, binding):
        if binding is None:
            yield
            return
        binding_id = binding.get("binding_id") if isinstance(binding, dict) else None
        if not isinstance(binding_id, str) or not re.fullmatch(r"[0-9a-f]{64}", binding_id):
            raise ValueError("Invalid frozen catalog binding")
        scoped = {"adapter": None, "scout": {}, "unavailable": None}
        try:
            record = json.loads((self.root / (binding_id + ".json")).read_bytes())
            if digest(record) != binding_id:
                raise ValueError("Frozen catalog registry record changed")
            identity = record["identity"]
            if identity["kind"] == "local":
                if identity.get("dependency_inventory_status") != "verified":
                    raise ValueError("Frozen source dependency inventory unavailable")
                config = _private_config(record["config"])
                placeholder = type("CatalogConfig", (), {"config": config})()
                if source_identity(placeholder, {}) != identity:
                    raise ValueError("Frozen installed catalog files changed or are missing")
                with self.lock:
                    adapter = self.adapters.get(binding_id)
                    if adapter is None:
                        adapter = self.adapter_factory(config)
                        self.adapters[binding_id] = adapter
                scoped["adapter"] = adapter
            elif identity["kind"] == "remote":
                # Credentials remain operator-owned and can rotate. An endpoint
                # change cannot silently supply evidence to an old run.
                scout = getattr(self.service, "scout", {})
                if getattr(self.service, "local_catalog", None) is not None or source_identity(None, scout) != identity:
                    raise ValueError("Frozen remote catalog endpoint is no longer configured")
                scoped["scout"] = scout
            elif identity["kind"] != "not_configured":
                raise ValueError("Unknown frozen catalog kind")
        except (ValueError, KeyError, TypeError, OSError, ImportError) as error:
            scoped["unavailable"] = type(error).__name__
        token = self.scope.set(scoped)
        try:
            yield
        finally:
            self.scope.reset(token)

    def close(self):
        with self.lock:
            adapters = list(self.adapters.values()) + [getattr(self.service, "local_catalog", None)]
            self.adapters.clear()
        seen = set()
        for adapter in adapters:
            if adapter is not None and id(adapter) not in seen:
                seen.add(id(adapter))
                adapter.close()
