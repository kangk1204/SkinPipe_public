"""Content identity of files opened by an Open Babel SMILES/pH startup probe.

The executable alone does not identify its loaded libraries or pH rules.
Linux strace observes four fixed SMILES at pH 7.4; actual opened regular files
are hashed. This is a startup-input observation, not a trace of every ligand
or evidence that a physiological protomer has been selected. Unsupported or
failed observation leaves preparation usable, but ineligible for cache reuse.
"""
from __future__ import annotations

import ast
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import threading
from collections import OrderedDict
from typing import Any

_INPUT = 'C\nCN\nCC(=O)O\nNCC(=O)O\n'
_TIMEOUT = 8
_MAX_TRACE = 2 << 20
_OPEN = re.compile(r'openat\((AT_FDCWD|\d+), ("(?:\\.|[^"\\])*"), .*\) = (-?\d+)')
_CACHE: OrderedDict[tuple[str, ...], dict[str, Any]] = OrderedDict()
_LOCK = threading.RLock()


def refresh() -> None:
    with _LOCK:
        _CACHE.clear()


def _fingerprint(path: Path) -> tuple[Any, ...]:
    try:
        s = path.lstat()
        return (s.st_dev, s.st_ino, s.st_mode, s.st_size, s.st_mtime_ns, s.st_ctime_ns)
    except OSError:
        return (None,)


def _watch(paths: set[Path]) -> dict[Path, tuple[Any, ...]]:
    # Missing higher-priority files matter too: e.g. a newly created local
    # phmodel.txt overrides installed rules. Parent/symlink changes are watched.
    watched = {p for path in paths for p in (path, *path.parents)}
    return {p: _fingerprint(p) for p in watched}


def _parse_trace(text: str, cwd: Path) -> tuple[set[Path], set[Path]]:
    opened: set[Path] = set()
    searched: set[Path] = set()
    for line in text.splitlines():
        if 'openat(' not in line:
            continue
        match = _OPEN.search(line)
        if match is None:
            raise ValueError('incomplete_openat_trace')
        base, quoted, result = match.groups()
        value = ast.literal_eval(quoted)
        path = Path(value)
        if not path.is_absolute():
            if base != 'AT_FDCWD':
                raise ValueError('unresolved_relative_openat')
            path = cwd / path
        # Normpath preserves lookup symlinks for freshness monitoring. resolve
        # is used only when hashing content, not to forget a searched location.
        path = Path(os.path.normpath(path))
        searched.add(path)
        if int(result) >= 0 and path.is_file():
            opened.add(path)
    return opened, searched


def _observe(binary: str, observer: str, cwd: Path) -> dict[str, Any]:
    with tempfile.TemporaryDirectory(prefix='skinpipe-obabel-identity-') as directory:
        trace = Path(directory) / 'openat.log'
        result = subprocess.run(
            [observer, '-f', '-s', '4096', '-e', 'trace=openat', '-o', str(trace),
             binary, '-ismi', '-osmi', '-p', '7.4'], input=_INPUT,
            capture_output=True, text=True, timeout=_TIMEOUT, cwd=cwd)
        if result.returncode != 0:
            raise ValueError('startup_probe_failed')
        if len([s for s in result.stdout.splitlines() if s.strip()]) != 4:
            raise ValueError('incomplete_startup_probe')
        if not trace.is_file() or trace.stat().st_size > _MAX_TRACE:
            raise ValueError('missing_or_oversize_trace')
        opened, searched = _parse_trace(trace.read_text(), cwd)
    if not any(p.name == 'phmodel.txt' for p in opened):
        # Compiled-in fallback exists in Open Babel, but is not certified by
        # this external-rule observation. Never infer it from a missing trace.
        raise ValueError('external_ph_rules_unobserved')
    if not any(p.name.startswith('libopenbabel') for p in opened):
        raise ValueError('native_openbabel_unobserved')
    return {'opened': opened, 'watch': _watch(searched | opened | {Path(binary), Path(observer), cwd})}


def identity(binary: str | None) -> dict[str, Any]:
    """Portable content identity; paths/environment remain internal to probing."""
    from .software import executable_identity, file_identity
    implementation = executable_identity(binary)
    result = dict(implementation)
    common = {'schema_version': 1, 'scope': 'four_smiles_ph7.4_startup_opened_files',
              'per_ligand_execution_trace': False}
    if not implementation['available']:
        result['runtime_inputs'] = dict(common, state='unavailable')
        return result
    observer = shutil.which('strace')
    if observer is None:
        result['runtime_inputs'] = dict(common, state='unverified', reason='observer_unavailable')
        return result
    cwd = Path.cwd()
    # Environment values are not exported or persisted, including credentials.
    # Any environment change forces rediscovery; compatibility uses only files.
    environment = hashlib.sha256(json.dumps(sorted(os.environ.items()), separators=(',', ':')).encode()).hexdigest()
    key = (str(binary), str(observer), str(cwd), environment)
    try:
        with _LOCK:
            discovery = _CACHE.get(key)
            if discovery is None or any(_fingerprint(p) != s for p, s in discovery['watch'].items()):
                discovery = _observe(str(binary), observer, cwd)
                _CACHE[key] = discovery
                _CACHE.move_to_end(key)
                while len(_CACHE) > 8:
                    _CACHE.popitem(last=False)
            entries = []
            for path in discovery['opened']:
                digest = file_identity(path)
                if digest is None:
                    raise ValueError('observed_input_unreadable')
                name = path.name
                kind = 'native' if re.search(r'\.so(?:\.|$)', name) else 'runtime_data'
                entries.append({'name': name, 'kind': kind, 'sha256': digest})
            # A list preserves same-named files with differing content. Absolute
            # locations never change compatibility when identical bytes move.
            entries.sort(key=lambda e: (e['name'], e['kind'], e['sha256']))
            manifest = hashlib.sha256(json.dumps(entries, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
            if any(_fingerprint(p) != s for p, s in discovery['watch'].items()):
                _CACHE.pop(key, None)
                raise ValueError('runtime_inputs_changed_during_snapshot')
        result['runtime_inputs'] = dict(common, state='observed', files=entries,
                                       files_sha256=manifest, files_count=len(entries))
    except (OSError, ValueError, SyntaxError, subprocess.TimeoutExpired) as exc:
        reason = str(exc) if type(exc) is ValueError else type(exc).__name__
        result['runtime_inputs'] = dict(common, state='unverified', reason=reason)
    return result
