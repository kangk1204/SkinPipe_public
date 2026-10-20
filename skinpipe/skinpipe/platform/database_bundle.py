"""Separate, versioned local database delivery; donor directories are read-only."""
from __future__ import annotations

import argparse
import glob
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import tarfile
import tempfile

from .artifacts import atomic_bytes, canonical_json, sha256_file
from .identity import digest

SCHEMA = "skinpipe-local-database-v1"
COMPONENTS = {"scout_root", "published_prediction_cache", "similarity_index", "admet_search_gate"}
FILTERED_COMPONENTS = {"filtered_risk", "filtered_purchase", "filtered_risk_purchase"}
ALL_COMPONENTS = COMPONENTS | FILTERED_COMPONENTS
SQL_LITERAL = re.compile(r"'(?:''|[^'])*'")


def _sql_paths(sql, *, positions=False):
    """Return complete SQL string literals; never substitute identifiers/code."""
    code = SQL_LITERAL.sub(lambda value: " " * len(value.group()), sql)
    functions = re.findall(r"\b([A-Za-z_][A-Za-z0-9_]*)\s*\(", code)
    if any((f.lower().startswith("read_") or f.lower().endswith("_scan"))
           and f.lower() not in {"read_parquet", "parquet_scan"} for f in functions):
        raise ValueError("Only explicit local Parquet reader functions may enter a portable database")
    argument_ranges = []
    for reader in re.finditer(r"\b(?:read_parquet|parquet_scan)\s*\(", code, re.I):
        start, depth, cursor = reader.end(), 1, reader.end()
        while cursor < len(code) and depth:
            depth += (code[cursor] == "(") - (code[cursor] == ")")
            cursor += 1
        if depth:
            raise ValueError("Incomplete Parquet reader expression")
        argument_ranges.append((start, cursor - 1))
    literals = [m for m in SQL_LITERAL.finditer(sql)
                if any(start <= m.start() < end for start, end in argument_ranges)]
    paths = []
    for literal in literals:
        value = literal.group()[1:-1].replace("''", "'")
        if re.match(r"(?:https?|s3|file|ftp):", value, re.I):
            raise ValueError("Database views may only read packaged local Parquet")
        if value.startswith("/") or value.lower().endswith(".parquet"):
            path = Path(value)
            if not path.is_absolute() or path.suffix.lower() != ".parquet" or any(c in value for c in "\\\x00"):
                raise ValueError("Database views require explicit local Parquet paths")
            paths.append((literal.start(), value) if positions else value)
    return paths


def _view_dependencies(database, locations):
    mapping = {}
    for _, _, sql in _views(database):
        for value in _sql_paths(sql):
            _component_path(value, locations)
            paths = sorted(glob.glob(value, recursive=True)) if glob.has_magic(value) else [value]
            if not paths or len(paths) > 200000:
                raise ValueError("Parquet view requires a nonempty bounded frozen file population")
            mapping[value] = [_component_path(path, locations) for path in paths]
    return mapping


def _views(database):
    import duckdb
    with duckdb.connect(str(database), read_only=True, config={"autoinstall_known_extensions": False, "autoload_known_extensions": False}) as db:
        return db.execute("SELECT schema_name,view_name,sql FROM duckdb_views() WHERE NOT internal ORDER BY schema_name,view_name").fetchall()


def _component_path(value, locations):
    path = Path(value)
    for name, root in sorted(locations.items(), key=lambda item: -len(str(item[1]))):
        if path.is_relative_to(root):
            relative = path.relative_to(root)
            if name == "scout_root" and (not relative.parts or relative.parts[0] != "data"):
                raise ValueError("Scout views must read admitted data files")
            return "components/" + name + "/" + relative.as_posix()
    raise ValueError("Database view references a file outside the packaged components")


def _entries(config):
    if not isinstance(config, dict) or config.get("kind", "local") != "local" or not config.get("scout_root"):
        raise ValueError("A local Scout configuration is required")
    if set(config) - (COMPONENTS | {"kind", "filtered_similarity_indices"}):
        raise ValueError("Only database component locations may enter a database package")
    if any(not isinstance(config[name], str) or not config[name] for name in COMPONENTS if name in config):
        raise ValueError("Database component locations must be nonempty local paths")
    locations = {name: Path(config[name]).resolve(strict=True) for name in COMPONENTS if config.get(name)}
    filtered = config.get("filtered_similarity_indices", {})
    if not isinstance(filtered, dict):
        raise ValueError("Filtered index locations must be an object")
    for key, value in filtered.items():
        if key not in {"risk", "purchase", "risk_purchase"}:
            raise ValueError("Unknown filtered index")
        if not isinstance(value, str) or not value:
            raise ValueError("Filtered index location must be a nonempty path")
        locations["filtered_" + key] = Path(value).resolve(strict=True)
    scout = locations["scout_root"]
    database = scout / "data/skinscout.duckdb"
    if not database.is_file():
        raise ValueError("Source has no Scout structure database")
    selected = {name: set() for name in locations}
    selected["scout_root"].update({database, scout / "reports/audit.json", scout / "reports/db_build.json"})
    for members in _view_dependencies(database, locations).values():
        for value in members:
            member = PurePosixPath(value)
            selected[member.parts[1]].add(locations[member.parts[1]] / Path(*member.parts[2:]))
    records = {}
    for name, root in sorted(locations.items()):
        if not root.is_dir():
            raise ValueError("Database component must be a directory")
        if name == "scout_root":
            # Ignore unrelated raw downloads, logs and credentials. An admitted
            # source link must never let file membership escape its component.
            for path in (root / "data").rglob("*"):
                if path.is_symlink():
                    raise ValueError("Database package may not contain symbolic links")
        else:
            allowed = ({"manifest.json", "predictions.duckdb"} if name == "published_prediction_cache"
                       else {"manifest.json", "catalog.duckdb"} if name == "admet_search_gate"
                       else {"manifest.json", "structures.h5", "lookup.duckdb", "unsupported-structures.jsonl"})
            selected[name].update(root / filename for filename in allowed if (root / filename).exists())
        for path in sorted(selected[name]):
            if not path.is_file() or path.is_symlink() or not path.resolve().is_relative_to(root):
                raise ValueError("Database package member is missing, symbolic or external")
            if Path(str(path) + ".wal").exists():
                raise ValueError("Database source must be checkpointed before packaging")
            relative = "components/" + name + "/" + path.relative_to(root).as_posix()
            records[relative] = {"source": path, "sha256": sha256_file(path), "size_bytes": path.stat().st_size}
    if "components/scout_root/data/skinscout.duckdb" not in records:
        raise ValueError("Source has no Scout structure database")
    return locations, records


def package_database(config, destination):
    destination = Path(destination)
    if destination.exists():
        raise FileExistsError("Database package destination already exists")
    locations, records = _entries(config)
    from .scout import ScoutCatalog
    source = ScoutCatalog(locations["scout_root"])
    source_seal = source.structure_fingerprint()
    source_identity = source.database_identity()
    _verify_components(locations, source_seal)
    total = sum(row["size_bytes"] for row in records.values())
    if shutil.disk_usage(destination.parent).free < total + 256 * 1024**2:
        raise ValueError("Insufficient space for the separate database archive")
    manifest = {"schema": SCHEMA, "components": sorted(locations),
                "source_locations": {name: str(path) for name, path in locations.items()},
                "files": {name: {k: v for k, v in row.items() if k != "source"} for name, row in records.items()},
                "total_bytes": total, "scope": "local research database; admitted view dependencies only; no credentials or inference",
                "source_structure_fingerprint": source_seal, "source_database_identity": source_identity,
                "view_dependencies": _view_dependencies(source.database, locations)}
    manifest["package_id"] = digest(manifest)
    manifest_raw = canonical_json(manifest)
    import io
    descriptor, temporary_name = tempfile.mkstemp(prefix="." + destination.name + ".building-", dir=destination.parent)
    os.close(descriptor)
    temporary = Path(temporary_name)
    states = {name: (row["source"].stat().st_size, row["source"].stat().st_mtime_ns, row["source"].stat().st_ctime_ns)
              for name, row in records.items()}
    try:
        with tarfile.open(temporary, "w:gz", compresslevel=1, dereference=True) as archive:
            info = tarfile.TarInfo("manifest.json")
            info.size, info.mode = len(manifest_raw), 0o600
            archive.addfile(info, io.BytesIO(manifest_raw))
            for name, row in records.items():
                before = row["source"].stat()
                archive.add(row["source"], arcname=name, recursive=False)
                after = row["source"].stat()
                if ((before.st_size, before.st_mtime_ns, before.st_ctime_ns) !=
                    (after.st_size, after.st_mtime_ns, after.st_ctime_ns)
                    or sha256_file(row["source"]) != row["sha256"]):
                    raise ValueError("Database source changed during packaging")
        if ScoutCatalog(locations["scout_root"]).structure_fingerprint() != source_seal:
            raise ValueError("Source structure fingerprint changed during packaging")
        for name, row in records.items():
            current = row["source"].stat()
            if (current.st_size, current.st_mtime_ns, current.st_ctime_ns) != states[name] or Path(str(row["source"]) + ".wal").exists():
                raise ValueError("Database component changed during packaging")
        os.link(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)
    return {"package_id": manifest["package_id"], "file_count": len(records), "total_bytes": total,
            "archive": str(destination), "archive_sha256": sha256_file(destination)}


def _manifest(archive):
    info = archive.next()
    if info is None or info.name != "manifest.json" or not info.isfile() or info.size > 32 * 1024**2:
        raise ValueError("Database package requires a bounded leading manifest")
    manifest = json.load(archive.extractfile(info))
    if manifest.get("schema") != SCHEMA or manifest.get("package_id") != digest({k: v for k, v in manifest.items() if k != "package_id"}):
        raise ValueError("Database package manifest mismatch")
    files = manifest.get("files")
    if not isinstance(files, dict) or not 1 <= len(files) <= 200000:
        raise ValueError("Invalid database file inventory")
    components, locations = manifest.get("components"), manifest.get("source_locations")
    if (not isinstance(components, list) or any(not isinstance(c, str) or c not in ALL_COMPONENTS for c in components)
            or len(set(components)) != len(components) or "scout_root" not in components
            or not isinstance(locations, dict) or set(locations) != set(components)
            or any(not isinstance(p, str) or not Path(p).is_absolute() or "\x00" in p for p in locations.values())):
        raise ValueError("Invalid database component inventory")
    for name, row in files.items():
        path = PurePosixPath(name)
        if (not isinstance(name, str) or "\\" in name or "\x00" in name or path.as_posix() != name
                or path.is_absolute() or ".." in path.parts or len(path.parts) < 3 or path.parts[0] != "components" or path.parts[1] not in components
                or not isinstance(row, dict) or set(row) != {"size_bytes", "sha256"}
                or type(row.get("size_bytes")) is not int or row["size_bytes"] < 0
                or not isinstance(row.get("sha256"), str) or not re.fullmatch("[0-9a-f]{64}", row["sha256"])):
            raise ValueError("Invalid database package member")
        component, relative = path.parts[1], path.parts[2:]
        if component == "scout_root":
            if relative not in {("data", "skinscout.duckdb"), ("reports", "audit.json"), ("reports", "db_build.json")} and not (relative[0] == "data" and path.suffix == ".parquet"):
                raise ValueError("Unadmitted Scout package file")
        else:
            allowed = ({"manifest.json", "predictions.duckdb"} if component == "published_prediction_cache"
                       else {"manifest.json", "catalog.duckdb"} if component == "admet_search_gate"
                       else {"manifest.json", "structures.h5", "lookup.duckdb", "unsupported-structures.jsonl"})
            if len(relative) != 1 or relative[0] not in allowed:
                raise ValueError("Unadmitted database component file")
    if sum(r["size_bytes"] for r in files.values()) != manifest.get("total_bytes"):
        raise ValueError("Database package size inventory differs")
    return manifest


def _verify_components(locations, source_seal):
    """Use the existing admission checks, including immutable data hashes."""
    from .admet_search import ADMETSearchGate
    from .catalog_index import CatalogSimilarityIndex
    from .scout_predictions import PublishedPredictionCache
    manifests = {}
    if "published_prediction_cache" in locations:
        manifests["published_prediction_cache"] = PublishedPredictionCache(locations["published_prediction_cache"]).manifest
    if "admet_search_gate" in locations:
        manifests["admet_search_gate"] = ADMETSearchGate(locations["admet_search_gate"]).manifest
    for name in ({"similarity_index"} | FILTERED_COMPONENTS) & set(locations):
        manifests[name] = CatalogSimilarityIndex(locations[name]).manifest
    for name, value in manifests.items():
        if name != "published_prediction_cache" and value.get("source_structure_fingerprint") != source_seal:
            raise ValueError("Database component source seal differs: " + name)
    gate = manifests.get("admet_search_gate")
    cache = manifests.get("published_prediction_cache")
    if gate and cache and gate.get("prediction_cache_manifest_sha256") != digest(cache):
        raise ValueError("Search gate and published prediction cache differ")
    full = manifests.get("similarity_index")
    for name in FILTERED_COMPONENTS & set(manifests):
        value = manifests[name]
        if not full or value.get("full_index_manifest_sha256") != digest(full):
            raise ValueError("Filtered and full index source manifests differ")
        if not gate or value.get("search_gate_manifest_sha256") != digest(gate):
            raise ValueError("Filtered index and search gate manifests differ")
        filters = {"hide_severe_admet": name in {"filtered_risk", "filtered_risk_purchase"},
                   "purchasable_only": name in {"filtered_purchase", "filtered_risk_purchase"}}
        if value.get("search_filters") != filters:
            raise ValueError("Filtered index declares a different population")
    return manifests


def _verify_original_seal(root, manifest):
    """Reconstruct the original seal from copied bytes, with donor offline."""
    database = root / "components/scout_root/data/skinscout.duckdb"
    source_root = Path(manifest["source_locations"]["scout_root"])
    dependencies = set()
    locations = {name: Path(value) for name, value in manifest["source_locations"].items()}
    for _, name, sql in _views(database):
        for value in _sql_paths(sql):
            member = _component_path(value, locations)
            members = manifest.get("view_dependencies", {}).get(value, [member])
            if not members or any(name not in manifest["files"] for name in members):
                raise ValueError("Database view dependency is absent from the package")
            if name in {"compounds", "offers", "source_records"}:
                if not Path(value).is_relative_to(source_root / "data"):
                    raise ValueError("Structure dependency is outside its original admitted data")
                if glob.has_magic(value):
                    raise ValueError("Structure-view glob cannot establish the existing Scout structure seal")
                dependencies.add((str(Path(value).relative_to(source_root)), member))
    identity = {"sha256": sha256_file(database), "size_bytes": database.stat().st_size,
                "structure_dependencies": [{"path": relative, "sha256": sha256_file(root / member),
                    "size_bytes": (root / member).stat().st_size} for relative, member in sorted(dependencies)]}
    audit = json.loads((root / "components/scout_root/reports/audit.json").read_bytes())
    build = root / "components/scout_root/reports/db_build.json"
    checks = {key: audit.get("checks", {}).get(key) for key in (
        "missing_or_stale_build_report", "incomplete_admitted_structure_sources", "duplicate_compound_ids",
        "orphan_structure_links", "duplicate_source_records")}
    if any(value != 0 for value in checks.values()) or not json.loads(build.read_bytes()).get("admitted_structure_sources_complete"):
        raise ValueError("Source structure admission is unresolved")
    seal = digest({"build_sha256": sha256_file(build), "structure_checks": checks,
        "compound_count": audit.get("compound_count", json.loads(build.read_bytes()).get("compounds")),
        "database_identity": identity, "seal_version": 2})
    if identity != manifest.get("source_database_identity") or seal != manifest.get("source_structure_fingerprint"):
        raise ValueError("Original structure seal cannot be established from packaged bytes")
    locations = {name: root / "components" / name for name in manifest["components"]}
    _verify_components(locations, seal)
    return identity


def _relocate(staging, final, manifest):
    """Rewrite only packaged Parquet literals; all executable SQL is unchanged."""
    import duckdb
    database = staging / "components/scout_root/data/skinscout.duckdb"
    originals, replacements, changes = _views(database), [], []
    locations = {name: Path(value) for name, value in manifest["source_locations"].items()}
    for schema, name, sql in originals:
        path_values = dict(_sql_paths(sql, positions=True))
        def replace(match):
            value = path_values.get(match.start())
            if value is None:
                return match.group()
            member = _component_path(value, locations)
            members = manifest.get("view_dependencies", {}).get(value, [member])
            if not members or any(name not in manifest["files"] for name in members):
                raise ValueError("Relocation dependency is not in the verified inventory")
            return "'" + str(final / member).replace("'", "''") + "'"
        updated = SQL_LITERAL.sub(replace, sql)
        replacements.append((schema, name, updated))
        if updated != sql:
            changes.append({"schema": schema, "view": name, "source_sql_sha256": digest(sql),
                            "installed_sql_sha256": digest(updated), "operation": "packaged_parquet_path_literals_only"})
    if not changes:
        return []
    with duckdb.connect(str(database), config={"autoinstall_known_extensions": False, "autoload_known_extensions": False}) as db:
        tables = db.execute("SELECT schema_name,table_name,sql FROM duckdb_tables() WHERE NOT internal ORDER BY schema_name,table_name").fetchall()
        for (schema, name, sql), (_, _, updated) in zip(originals, replacements):
            if updated != sql:
                quoted = '"' + schema.replace('"', '""') + '"."' + name.replace('"', '""') + '"'
                db.execute("DROP VIEW " + quoted)
                db.execute("SET schema = '" + schema.replace("'", "''") + "'")
                db.execute(updated)
        if db.execute("SELECT schema_name,table_name,sql FROM duckdb_tables() WHERE NOT internal ORDER BY schema_name,table_name").fetchall() != tables:
            raise ValueError("Relocation changed stored table definitions")
        db.execute("CHECKPOINT")
    installed = _views(database)
    # DuckDB canonicalizes DDL when storing it. The canonical SQL before and
    # after must differ only in those same literal paths, never a SELECT term.
    for (schema, name, actual), (expected_schema, expected_name, expected) in zip(installed, replacements):
        if (schema, name) != (expected_schema, expected_name) or actual != expected:
            raise ValueError("Relocation was not an exact path-only SQL transformation")
    if len(installed) != len(replacements):
        raise ValueError("Relocation changed the view inventory")
    return changes


def _rebind_installed(root, manifest, changes, original_identity):
    """Re-seal preserved populations after a proven path-only DB relocation."""
    from .scout import ScoutCatalog
    old_locations = {name: Path(value) for name, value in manifest["source_locations"].items()}
    for source_path, members in manifest.get("view_dependencies", {}).items():
        pattern = root / _component_path(source_path, old_locations)
        actual = sorted(glob.glob(str(pattern), recursive=True)) if glob.has_magic(str(pattern)) else [str(pattern)]
        if actual != sorted(str(root / name) for name in members):
            raise ValueError("Installed view file population differs from its frozen inventory")
    scout = ScoutCatalog(root / "components/scout_root")
    new_seal, new_identity = scout.structure_fingerprint(), scout.database_identity()
    if new_identity["structure_dependencies"] != original_identity["structure_dependencies"]:
        raise ValueError("Relocation changed structure dependency identities; rebuild required")
    old_seal = manifest["source_structure_fingerprint"]
    if new_seal != old_seal and not changes:
        raise ValueError("Physical source changed without a proven path-only transformation")
    relocation = {"schema": "verified-database-relocation-v1", "source_structure_fingerprint": old_seal,
        "installed_source_structure_fingerprint": new_seal, "source_database_identity": original_identity,
        "installed_database_identity": new_identity, "sql_changes": changes,
        "preserved_payloads": {name: row for name, row in manifest["files"].items()
                               if name != "components/scout_root/data/skinscout.duckdb"},
        "inference_performed": False, "population_reindexed": False}
    relocation_sha = digest(relocation)
    atomic_bytes(root / "relocation.json", canonical_json(relocation))
    locations = {name: root / "components" / name for name in manifest["components"]}
    originals = _verify_components(locations, old_seal)
    rebound = {}
    order = ["admet_search_gate", "similarity_index", *sorted(FILTERED_COMPONENTS)]
    for name in order:
        if name not in originals:
            continue
        original = originals[name]
        updated = {**original, "source_structure_fingerprint": new_seal,
            "relocation": {"receipt_sha256": relocation_sha, "predecessor_manifest_sha256": digest(original),
                           "source_structure_fingerprint": old_seal}}
        if name == "admet_search_gate":
            if original.get("source_database_identity") is not None and original["source_database_identity"] != original_identity:
                raise ValueError("Search gate original database identity differs; rebuild required")
            updated["source_database_identity"] = new_identity
        if name in FILTERED_COMPONENTS:
            updated["search_gate_manifest_sha256"] = digest(rebound["admet_search_gate"])
            updated["full_index_manifest_sha256"] = digest(rebound["similarity_index"])
        atomic_bytes(root / "relocation" / (name + "-predecessor.json"), canonical_json(original))
        atomic_bytes(locations[name] / "manifest.json", canonical_json(updated))
        rebound[name] = updated
    _verify_components(locations, new_seal)
    return {"receipt_sha256": relocation_sha, "source_structure_fingerprint": old_seal,
            "installed_source_structure_fingerprint": new_seal, "rebound_components": sorted(rebound),
            "population_reindexed": False}


def install_database(package, destination):
    package, destination = Path(package), Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    with tarfile.open(package, "r|gz") as archive:
        manifest = _manifest(archive)
        final = destination / ("database-" + manifest["package_id"][:24])
        if final.exists():
            raise FileExistsError("This database version is already installed")
        if shutil.disk_usage(destination).free < manifest["total_bytes"] + 512 * 1024**2:
            raise ValueError("Insufficient database and temporary installation space")
        staging = Path(tempfile.mkdtemp(prefix=".database-install-", dir=destination))
        try:
            seen = set()
            for info in iter(archive.next, None):
                if not info.isfile() or info.name not in manifest["files"] or info.name in seen:
                    raise ValueError("Unexpected, duplicate or non-file database archive member")
                row = manifest["files"][info.name]
                if info.size != row["size_bytes"]:
                    raise ValueError("Database archive member size differs")
                path = staging / info.name
                path.parent.mkdir(parents=True, exist_ok=True)
                with archive.extractfile(info) as source, path.open("xb") as target:
                    shutil.copyfileobj(source, target, length=1024 * 1024)
                if sha256_file(path) != row["sha256"]:
                    raise ValueError("Database archive member checksum differs")
                seen.add(info.name)
            if seen != set(manifest["files"]):
                raise ValueError("Database archive is incomplete")
            original_identity = _verify_original_seal(staging, manifest)
            # Publish before relocation so DuckDB can bind final local paths;
            # the active pointer remains untouched until all checks pass.
            staging.rename(final)
            changes = _relocate(final, final, manifest)
            rebinding = _rebind_installed(final, manifest, changes, original_identity)
            config = {"kind": "local", **{name: str(final / "components" / name) for name in manifest["components"] if name in COMPONENTS}}
            filtered = {name.removeprefix("filtered_"): str(final / "components" / name) for name in manifest["components"] if name.startswith("filtered_")}
            if filtered:
                config["filtered_similarity_indices"] = filtered
            receipt = {"schema": SCHEMA, "package_id": manifest["package_id"], "source_manifest": manifest,
                       "relocation": changes, "installed_database_sha256": sha256_file(final / "components/scout_root/data/skinscout.duckdb"),
                       "config": config, "inference_performed": False, "verified_relocation": rebinding}
            atomic_bytes(final / "installation.json", canonical_json(receipt))
            atomic_bytes(final / "local-catalog.json", canonical_json(config))
            atomic_bytes(destination / "active.json", canonical_json({"version": final.name, "package_id": manifest["package_id"]}))
            return receipt
        except BaseException as error:
            if staging.exists():
                shutil.rmtree(staging)
            elif final.exists():
                atomic_bytes(final / "installation-failed.json", canonical_json({
                    "status": "inactive_rebuild_required", "package_id": manifest["package_id"],
                    "error_type": type(error).__name__, "reason": str(error)[:400],
                    "active_pointer_updated": False, "inference_performed": False}))
            # An unsuccessful version remains separate and never becomes active.
            raise


def register_database(config, destination):
    """Register an existing admitted local DB read-only, without copying it."""
    from .scout import ScoutCatalog
    locations, _ = _entries(config)
    scout = ScoutCatalog(locations["scout_root"])
    seal = scout.structure_fingerprint()
    _verify_components(locations, seal)
    normalized = {"kind": "local", **{name: str(path) for name, path in locations.items() if name in COMPONENTS}}
    filtered = {name.removeprefix("filtered_"): str(path) for name, path in locations.items() if name in FILTERED_COMPONENTS}
    if filtered:
        normalized["filtered_similarity_indices"] = filtered
    record = {"schema": SCHEMA, "registration_mode": "existing_read_only", "config": normalized,
              "source_structure_fingerprint": seal, "database_identity": scout.database_identity(),
              "inference_performed": False, "portable_copy": False}
    identity = digest(record)
    record["package_id"] = identity
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    version = destination / ("registered-" + identity[:24])
    if not version.exists():
        staging = Path(tempfile.mkdtemp(prefix=".database-register-", dir=destination))
        try:
            atomic_bytes(staging / "installation.json", canonical_json(record))
            atomic_bytes(staging / "local-catalog.json", canonical_json(normalized))
            staging.rename(version)
        finally:
            if staging.exists():
                shutil.rmtree(staging)
    elif (version / "installation.json").read_bytes() != canonical_json(record):
        raise ValueError("Existing database registration differs")
    atomic_bytes(destination / "active.json", canonical_json({"version": version.name, "package_id": identity}))
    return record


def database_status(destination, *, verify=False):
    """Report a registration cheaply; explicit verification rechecks file seals."""
    from .scout import ScoutCatalog
    destination = Path(destination)
    pointer = destination / "active.json"
    if not pointer.is_file():
        return {"status": "not_configured", "config": None, "integrity": "not_rechecked"}
    active = json.loads(pointer.read_bytes())
    name, identity = active.get("version"), active.get("package_id")
    if (not isinstance(name, str) or not re.fullmatch(r"(?:database|registered)-[0-9a-f]{24}", name)
            or not isinstance(identity, str) or not re.fullmatch(r"[0-9a-f]{64}", identity)
            or not name.endswith(identity[:24])):
        raise ValueError("Invalid active database pointer")
    version = destination / name
    if version.is_symlink() or not version.is_dir():
        raise ValueError("Active database version is missing or symbolic")
    receipt = json.loads((version / "installation.json").read_bytes())
    config = json.loads((version / "local-catalog.json").read_bytes())
    if receipt.get("package_id") != identity or receipt.get("config") != config:
        raise ValueError("Active database configuration differs from its receipt")
    response = {"status": "registered", "version": name, "package_id": identity, "config": config,
                "integrity": "not_rechecked", "registration_mode": receipt.get("registration_mode", "portable_installed_copy"),
                "scientific_validation": "NOT_ESTABLISHED"}
    if verify:
        scout = ScoutCatalog(config["scout_root"])
        seal = scout.structure_fingerprint()
        expected = receipt.get("source_structure_fingerprint") if name.startswith("registered-") else receipt["verified_relocation"]["installed_source_structure_fingerprint"]
        if seal != expected:
            raise ValueError("Registered database changed; re-registration or rebuild required")
        locations = {key: Path(config[key]) for key in COMPONENTS if config.get(key)}
        locations.update({"filtered_" + key: Path(path) for key, path in config.get("filtered_similarity_indices", {}).items()})
        _verify_components(locations, seal)
        response.update(integrity="file_and_population_seals_verified", source_structure_fingerprint=seal)
    return response


def main(argv=None):
    parser = argparse.ArgumentParser(description="SkinPipe 별도 로컬 DB 패키지·설치")
    sub = parser.add_subparsers(dest="command", required=True)
    package = sub.add_parser("package")
    package.add_argument("--config", required=True)
    package.add_argument("--output", required=True)
    install = sub.add_parser("install")
    install.add_argument("--package", required=True)
    install.add_argument("--destination", required=True)
    register = sub.add_parser("register", help="기존 로컬 DB를 읽기 전용으로 등록")
    register.add_argument("--config", required=True)
    register.add_argument("--destination", required=True)
    status = sub.add_parser("status", help="현재 DB 등록 상태")
    status.add_argument("--destination", required=True)
    status.add_argument("--verify", action="store_true", help="전체 파일·소스 봉인을 다시 확인")
    args = parser.parse_args(argv)
    if args.command == "package":
        result = package_database(json.loads(Path(args.config).read_text()), args.output)
    elif args.command == "install":
        result = install_database(args.package, args.destination)
    elif args.command == "register":
        result = register_database(json.loads(Path(args.config).read_text()), args.destination)
    else:
        result = database_status(args.destination, verify=args.verify)
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
