"""Operational commands; configuration and private catalog data stay outside Git.

Legacy ingestion opens ResultStore read-only. It deliberately never constructs
SkinPipeBridge, whose restart recovery changes the old SQLite job ledger.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import stat
import sys
import tempfile
import time
import uuid

from .artifacts import atomic_bytes, canonical_json, fsync_directory, sha256_file
from .resources import probe
from .service import PlatformService


def _outside_git(path):
    path = Path(path).expanduser().resolve()
    if any((parent / '.git').exists() for parent in [path, *path.parents]):
        raise ValueError('Runtime configuration and private state must be outside Git checkouts')
    return path


def _validate_config(config):
    if not isinstance(config, dict) or not isinstance(config.get('database'), dict) or not config.get('root'):
        raise ValueError('Runtime config requires root and a database connection object')
    config['root'] = str(_outside_git(config['root']))
    tokens = config.get('worker_tokens', {})
    if not isinstance(tokens, dict) or any(not isinstance(key, str) or not isinstance(value, str) or not value for key, value in tokens.items()):
        raise ValueError('worker_tokens must map worker IDs to nonempty tokens')
    return config


def read_config(path):
    path = _outside_git(path)
    if path.stat().st_size > 1024 * 1024:
        raise ValueError('Runtime config exceeds 1 MiB')
    if path.stat().st_mode & (stat.S_IRWXG | stat.S_IRWXO):
        raise ValueError('Runtime config must be private to its owner (chmod 600)')
    return _validate_config(json.loads(path.read_text()))


def create_config(path, *, root, database_file, legacy=None, worker_tokens_file=None):
    path = _outside_git(path)
    root = _outside_git(root)
    database = json.loads(Path(database_file).read_text())
    if not isinstance(database, dict):
        raise ValueError('Database configuration must be an object')
    config = {'schema_version': 1, 'root': str(root), 'database': database,
              'worker_tokens': json.loads(Path(worker_tokens_file).read_text()) if worker_tokens_file else {},
              'legacy': legacy or {},
              'discovery': {'literature': {'enabled': True, 'max_papers': 8,
                                           'timeout': 20, 'review_provider': 'none'}}}
    config = _validate_config(config)
    data = canonical_json(config)
    if len(data) > 1024 * 1024:
        raise ValueError('Runtime config exceeds 1 MiB')
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix='.runtime-config-', dir=path.parent)
    try:
        with os.fdopen(descriptor, 'wb') as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
        # Publish complete private bytes without replacing a concurrent or
        # existing configuration. Invalid input never reserves the final path.
        os.link(temporary, path)
        fsync_directory(path.parent)
    finally:
        os.unlink(temporary)
    return config


def _load_legacy_result_store(directory):
    directory = Path(directory).resolve()
    # Load exactly the configured private bridge implementation. Importing the
    # server module would run optional application imports and is unnecessary.
    modules = {}
    previous = {}
    try:
        for name, filename in (('docking_metadata', 'docking_metadata.py'), ('_skinpipe_legacy_results', 'results.py')):
            path = directory / filename
            if not path.is_file():
                raise ValueError('Legacy ResultStore module missing: ' + filename)
            spec = importlib.util.spec_from_file_location(name, path)
            module = importlib.util.module_from_spec(spec)
            previous[name] = sys.modules.get(name)
            sys.modules[name] = module
            spec.loader.exec_module(module)
            modules[name] = module
        return modules['_skinpipe_legacy_results'].ResultStore
    finally:
        for name, prior in previous.items():
            if prior is None:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = prior


class LegacyResultAdapter:
    """The read-only subset of the legacy bridge needed for source ingestion."""
    def __init__(self, config):
        if not config.get('repo') or not config.get('state_dir'):
            raise ValueError('Legacy configuration requires repo and state_dir')
        self.repo = Path(config['repo']).resolve()
        self.state_dir = Path(config['state_dir']).resolve()
        directory = Path(config.get('bridge_dir') or Path(__file__).resolve().parents[3] / 'web' / 'bridge')
        if not directory.is_dir():
            raise ValueError('Configure legacy.bridge_dir for the installed ResultStore implementation')
        result_store = _load_legacy_result_store(directory)
        self.results = result_store([self.repo / '결과물' / name for name in ('screens', 'screens_staging')]
                                    + [self.state_dir / 'runs', self.state_dir / 'imports'])

    def screen_files(self):
        found = {}
        sources = [(p.parent.name, p) for p in sorted((self.repo / '결과물' / 'screens').glob('*/result.json'))]
        sources += [('site-' + p.parent.parent.name, p) for p in sorted((self.state_dir / 'runs').glob('*/*/result.json'))]
        sources += [('import-' + p.parent.parent.name, p) for p in sorted((self.state_dir / 'imports').glob('*/*/result.json'))]
        for screen_id, path in sources:
            if screen_id in found:
                raise ValueError('Legacy screen ID collision; preserve both sources before assigning explicit IDs')
            found[screen_id] = self.results.safe(path)
        return found

    def screen_file(self, screen_id):
        try:
            return self.screen_files()[screen_id]
        except KeyError as exc:
            raise ValueError('Unknown legacy screen') from exc


def service_from_config(config, *, legacy=False):
    from psycopg.conninfo import make_conninfo

    adapter = LegacyResultAdapter(config.get('legacy', {})) if legacy else None
    return PlatformService(config['root'], make_conninfo(**config['database']), bridge=adapter,
                           worker_tokens=config.get('worker_tokens', {}), scout=config.get('scout'),
                           discovery_config=config.get('discovery'), access_config=config.get('access'))


def legacy_inventory(adapter):
    inventory = []
    for screen_id, path in adapter.screen_files().items():
        digest = sha256_file(path)
        _, data = adapter.results.read(path)
        package = adapter.results.package(path)
        if isinstance(data.get('discovery_candidates'), list):
            ids = [row.get('ligand_id') for row in data['discovery_candidates']]
        else:
            ids = sorted({row['ligand_id'] for key in ('candidates', 'dropped')
                          for row in data.get(key, []) if row.get('ligand_id')})
        if len(ids) != len(set(ids)) or any(not isinstance(value, str) or not value for value in ids):
            raise ValueError('Legacy source has missing or duplicate candidate IDs: ' + screen_id)
        if sha256_file(path) != digest:
            raise ValueError('Legacy source changed during inventory: ' + screen_id)
        inventory.append({'screen_id': screen_id, 'source_path': str(path), 'source_file_sha256': digest,
                          'marker': data.get('marker'),
                          'scientific_source_sha256': package['manifest']['source_sha256'] if package else digest,
                          'candidate_count': len(ids), 'candidate_ids_sha256': hashlib.sha256(canonical_json(sorted(ids))).hexdigest(),
                          'package_manifest_sha256': hashlib.sha256(canonical_json(package['manifest'])).hexdigest() if package else None,
                          'package_files': package['manifest']['files'] if package else None})
    return sorted(inventory, key=lambda row: row['screen_id'])


def import_legacy(service, *, expected_screens=None, expected_targets=None, expected_markers=None, activate=True):
    inventory = legacy_inventory(service.bridge)
    if not inventory:
        raise ValueError('No legacy screens were found; nothing was activated')
    if expected_screens is not None and len(inventory) != expected_screens:
        raise ValueError('Legacy screen inventory does not match the required baseline count')
    markers = {row['marker'] for row in inventory if isinstance(row['marker'], str) and row['marker']}
    if expected_markers is not None and len(markers) != expected_markers:
        raise ValueError('Legacy marker inventory does not match the required baseline count')
    report = {'schema_version': 1, 'kind': 'legacy_source_migration',
              'started_at': datetime.now(timezone.utc).isoformat(), 'inventory': inventory,
              'inventory_sha256': hashlib.sha256(canonical_json(inventory)).hexdigest(),
              'snapshots': [], 'activated': False, 'scientific_validation': 'not_established'}
    receipt_dir = service.root / 'migration-receipts'
    receipt_path = receipt_dir / (datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S') + '-' + uuid.uuid4().hex + '.json')
    atomic_bytes(receipt_path, canonical_json(report))
    try:
        targets = set()
        mapping = {}
        for source in inventory:
            manifest = service.import_legacy(source['screen_id'], activate=False,
                                              expected_source_sha256=source['source_file_sha256'])
            source_provenance = manifest['provenance']
            if manifest['candidate_count'] != source['candidate_count'] or source_provenance['source_candidate_ids_sha256'] != source['candidate_ids_sha256']:
                raise ValueError('Imported candidate IDs or counts differ from the frozen source inventory')
            service.catalog.verify(manifest['snapshot_id'])
            targets.update(manifest['target_counts'])
            mapping[source['screen_id']] = manifest['snapshot_id']
            report['snapshots'].append({'screen_id': source['screen_id'], 'snapshot_id': manifest['snapshot_id'],
                                        'candidate_count': manifest['candidate_count'], 'target_counts': manifest['target_counts']})
            atomic_bytes(receipt_path, canonical_json(report))
        if expected_targets is not None and len(targets) != expected_targets:
            raise ValueError('Imported target inventory does not match the required baseline count')
        if legacy_inventory(service.bridge) != inventory:
            raise ValueError('Legacy sources changed before activation; prior active map is preserved')
        report.update(status='verified', screen_count=len(inventory), target_count=len(targets), marker_count=len(markers),
                      candidate_count=sum(item['candidate_count'] for item in inventory),
                      finished_at=datetime.now(timezone.utc).isoformat())
        if activate:
            report['previous_active_snapshots'] = service._update_active_snapshots(mapping)
            report['activated'] = True
        atomic_bytes(receipt_path, canonical_json(report))
    except Exception as exc:
        report.update(status='failed', error=str(exc), finished_at=datetime.now(timezone.utc).isoformat())
        atomic_bytes(receipt_path, canonical_json(report))
        raise
    return {key: value for key, value in report.items() if key != 'inventory'} | {'receipt_path': str(receipt_path)}


def verify_service(service):
    reports = []
    for snapshot in service.catalog.list_snapshots():
        try:
            reports.append({'snapshot_id': snapshot['snapshot_id'], 'status': 'verified',
                            'verification': service.catalog.verify(snapshot['snapshot_id'])})
        except ValueError as exc:
            reports.append({'snapshot_id': snapshot['snapshot_id'], 'status': 'quarantined', 'reason': str(exc)})
    path = service.root / 'active-snapshots.json'
    active = json.loads(path.read_text()) if path.exists() else {}
    verified = {row['snapshot_id'] for row in reports if row['status'] == 'verified'}
    broken = {screen: snapshot for screen, snapshot in active.items() if snapshot not in verified}
    return {'status': 'verified' if not broken and all(row['status'] == 'verified' for row in reports) else 'needs_attention',
            'snapshots': reports, 'active_screen_count': len(active), 'broken_active_references': broken,
            'scientific_validation': 'not_established'}


def benchmark_catalog(service, snapshot_id, *, target_id=None, runs=50):
    if isinstance(runs, bool) or not isinstance(runs, int) or not 10 <= runs <= 1000:
        raise ValueError('Benchmark runs must be between 10 and 1000')
    manifest = service.catalog.manifest(snapshot_id)
    variants = ({}, {'hide_critical': True, 'hide_opposing': True}, {'measured_only': True})
    for filters in variants:
        service.catalog.query(snapshot_id, target_id=target_id, filters=filters)
    elapsed = []
    for index in range(runs):
        start = time.perf_counter()
        service.catalog.query(snapshot_id, target_id=target_id, filters=variants[index % len(variants)])
        elapsed.append(time.perf_counter() - start)
    ordered = sorted(elapsed)
    p95 = ordered[max(0, math.ceil(len(ordered) * .95) - 1)]
    report = {'kind': 'operational_catalog_query_benchmark', 'snapshot_id': snapshot_id,
              'snapshot_candidate_count': manifest['candidate_count'], 'target_id': target_id,
              'runs': runs, 'warmup_runs': len(variants), 'warm_p95_seconds': p95,
              'warm_median_seconds': ordered[len(ordered) // 2], 'query_variants': list(variants),
              'latencies_seconds': elapsed, 'p95_target_seconds': 2.0, 'target_met': p95 <= 2,
              'scientific_validation': 'not_established',
              'scope': 'Measured local operational query latency; no prediction or biological SOTA claim'}
    destination = service.root / 'benchmarks' / (snapshot_id + '-' + uuid.uuid4().hex + '.json')
    atomic_bytes(destination, canonical_json(report))
    return {**report, 'report_path': str(destination)}


def build_parser():
    parser = argparse.ArgumentParser(prog='skinpipe-platform')
    parser.add_argument('--config', default=os.environ.get('SKINPIPE_PLATFORM_CONFIG'))
    commands = parser.add_subparsers(dest='command', required=True)
    initialize = commands.add_parser('init', help='Initialize PostgreSQL schema and private runtime state')
    initialize.add_argument('--create-config', action='store_true')
    initialize.add_argument('--root')
    initialize.add_argument('--database-config')
    initialize.add_argument('--worker-tokens-file')
    initialize.add_argument('--legacy-repo')
    initialize.add_argument('--legacy-state-dir')
    initialize.add_argument('--legacy-bridge-dir')
    initialize.add_argument('--legacy-cli')
    initialize.add_argument('--allowed-user', default='')
    initialize.add_argument('--origin', action='append', default=[])
    commands.add_parser('probe', help='Read local resource and configured engine availability')
    migrate = commands.add_parser('import-legacy', help='Verify and import every existing legacy screen')
    migrate.add_argument('--expected-screens', type=int)
    migrate.add_argument('--expected-targets', type=int)
    migrate.add_argument('--expected-markers', type=int)
    migrate.add_argument('--no-activate', action='store_true')
    commands.add_parser('verify', help='Verify immutable snapshots and active screen references')
    benchmark = commands.add_parser('benchmark', help='Measure real warm catalog query latency')
    benchmark.add_argument('--snapshot', required=True)
    benchmark.add_argument('--target')
    benchmark.add_argument('--runs', type=int, default=50)
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if not args.config:
        parser.error('--config or SKINPIPE_PLATFORM_CONFIG is required')
    try:
        if args.command == 'init' and args.create_config:
            if not args.root or not args.database_config:
                parser.error('--create-config requires --root and --database-config')
            legacy = {'repo': args.legacy_repo, 'state_dir': args.legacy_state_dir,
                      'bridge_dir': args.legacy_bridge_dir, 'cli': args.legacy_cli,
                      'allowed_user': args.allowed_user, 'origins': args.origin}
            config = create_config(args.config, root=args.root, database_file=args.database_config,
                                   legacy=legacy, worker_tokens_file=args.worker_tokens_file)
        else:
            config = read_config(args.config)
        if args.command == 'probe':
            result = probe(directory=config['root'], engine_paths=config.get('engine_paths', {}))
        else:
            service = service_from_config(config, legacy=args.command == 'import-legacy')
            if args.command == 'init':
                result = {'status': 'initialized', 'operational_store': 'postgresql', 'root': config['root']}
            elif args.command == 'import-legacy':
                service.import_progress = lambda event: print(canonical_json(event).decode(), file=sys.stderr, flush=True)
                result = import_legacy(service, expected_screens=args.expected_screens,
                                       expected_targets=args.expected_targets, expected_markers=args.expected_markers,
                                       activate=not args.no_activate)
            elif args.command == 'verify':
                result = verify_service(service)
            else:
                result = benchmark_catalog(service, args.snapshot, target_id=args.target, runs=args.runs)
        print(canonical_json(result).decode())
        return 0 if result.get('status') not in {'failed', 'needs_attention'} else 1
    except Exception as exc:
        # Connection exceptions can contain DSNs or credentials. Report type and
        # safe validation messages only; never serialize the configuration.
        message = str(exc) if isinstance(exc, (ValueError, FileExistsError, FileNotFoundError)) else type(exc).__name__
        print(canonical_json({'status': 'failed', 'error': message}).decode(), file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
