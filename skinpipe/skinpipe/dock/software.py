"""Version and byte identities for the engines and ligand-preparation code.

Paths are diagnostic metadata; content digests determine cache compatibility.
This is a software identity, not proof of numerical agreement between devices.
"""
from __future__ import annotations

import hashlib
import importlib.metadata
import importlib.util
import json
import os
import platform
import shutil
import copy
import re
import stat
from functools import lru_cache
from pathlib import Path
from typing import Any, Optional


@lru_cache(maxsize=8192)
def _digest(path: str, fingerprint: tuple[int, ...]) -> str:
    # fingerprint participates in lru_cache's key. It invalidates cached bytes
    # after a file-stat change; it is deliberately not part of the byte digest.
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def file_identity(path: Optional[str | Path]) -> Optional[str]:
    if path is None:
        return None
    try:
        resolved = Path(path).resolve(strict=True)
        stat = resolved.stat()
        if not resolved.is_file():
            return None
        return _digest(str(resolved), (stat.st_dev, stat.st_ino, stat.st_size,
                                      stat.st_mtime_ns, stat.st_ctime_ns))
    except OSError:
        return None


def _file_snapshot(path: Path) -> tuple[str, Optional[tuple[int, ...]]]:
    try:
        # stat follows symlinks and binds target device/inode as well as both
        # timestamps. Resolving every ancestor twice for every namespace file
        # adds avoidable latency; content hashing still resolves its real path.
        s = path.stat()
        fingerprint = (s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns, s.st_ctime_ns)
        return str(path.absolute()), fingerprint if stat.S_ISREG(s.st_mode) else None
    except OSError:
        return str(path.absolute()), None


def _directory_snapshot(path: Path) -> Optional[tuple[int, ...]]:
    try:
        s = path.stat()
        return (s.st_dev, s.st_ino, s.st_mtime_ns, s.st_ctime_ns)
    except OSError:
        return None


@lru_cache(maxsize=8)
def _package_snapshot(version: Optional[str], metadata_available: bool,
                      snapshots: tuple[tuple[str, str, Optional[tuple[int, ...]]], ...]) -> dict[str, Any]:
    entries = [(name, file_identity(path)) for name, path, _ in snapshots]
    present = [(name, digest) for name, digest in entries if digest is not None]
    missing = [name for name, digest in entries if digest is None]
    encoded = json.dumps(sorted(entries), separators=(',', ':')).encode()
    return {'version': version,
            'package_sha256': hashlib.sha256(encoded).hexdigest() if present else None,
            'installed_files_hashed': len(present), 'missing_manifest_files': missing,
            'native_files': {name: digest for name, digest in present
                             if re.search(r'\.so(?:\.[0-9]+)*$', name) or
                             name.endswith(('.pyd', '.dll', '.dylib'))},
            'available': metadata_available or bool(present),
            'snapshot_state': 'observed' if present else 'unverified',
            'scope': 'current on-disk package namespace and declared files; not loaded-code memory attestation'}


def _package_identity(distribution: str, prefix: str) -> dict[str, Any]:
    """Fresh namespace/manifest inventory; cache hashes using full file stats.

    A package-name-only cache silently ignores changes until the next screen.
    Conda egg-info manifests can also omit installed compiled modules. Include
    namespace files, and cache only snapshots with matching file fingerprints.
    """
    try:
        dist = importlib.metadata.distribution(distribution)
    except importlib.metadata.PackageNotFoundError:
        dist = None
    entries: dict[str, Path] = {}
    for file in sorted(dist.files or [], key=str) if dist else []:
        name = str(file).replace('\\', '/')
        in_package = name.startswith(prefix + '/') or name.startswith(prefix + '.libs/')
        # Vina's compiled SWIG extension can be installed at site-packages
        # root, outside vina/. It is the executing native code, not optional.
        root_native = (prefix == 'vina' and '/' not in name and
                       name.startswith(('_vina', 'vina.')) and
                       name.endswith(('.so', '.pyd', '.dll', '.dylib')))
        if not (in_package or root_native) or '__pycache__' in name or name.endswith('.pyc'):
            continue
        entries[name] = Path(dist.locate_file(file))
    if prefix == 'vina':
        # Resolve the wrapper before scanning: importing its parent can create
        # ignored bytecode caches during otherwise normal first-use startup.
        try:
            spec = importlib.util.find_spec('vina._vina_wrapper')
            origin = spec.origin if spec else None
        except (ImportError, ValueError):
            origin = None
        name = 'executing_native/' + (Path(origin).name if origin else 'vina._vina_wrapper')
        entries[name] = Path(origin) if origin else Path('/__skinpipe_missing_vina_native__')
    roots: dict[Path, str] = {}
    if dist is not None:
        for folder in (prefix, prefix+'.libs'):
            root = Path(dist.locate_file(folder)).resolve()
            if root.is_dir():
                roots[root] = folder
    try:
        spec = importlib.util.find_spec(prefix)
        locations = getattr(spec, 'submodule_search_locations', None) or []
    except (ImportError, ValueError):
        locations = []
    for location in locations:
        root = Path(location).resolve()
        if root.is_dir() and root not in roots:
            roots[root] = 'executing_namespace/'+prefix
    namespace_paths = set()
    directories = {}
    for root, label in roots.items():
        directories[root] = _directory_snapshot(root)
        for path in sorted(root.rglob('*')):
            if '__pycache__' in path.parts or path.suffix in ('.pyc', '.pyo'):
                continue
            if path.is_file():
                name = label+'/'+path.relative_to(root).as_posix()
                entries[name] = path
                namespace_paths.add(name)
            elif path.is_dir():
                directories[path] = _directory_snapshot(path)
    snapshots = tuple((name, *_file_snapshot(path)) for name, path in sorted(entries.items()))
    version = dist.version if dist else None
    result = copy.deepcopy(_package_snapshot(version, dist is not None, snapshots))
    result['namespace_files_inventoried'] = len(namespace_paths)
    if (any(_file_snapshot(entries[name]) != (path, fingerprint) for name, path, fingerprint in snapshots)
            or any(_directory_snapshot(path) != stamp for path, stamp in directories.items())
            or dist is not None and dist.version != version):
        result.update(snapshot_state='unverified', snapshot_reason='input_changed_during_package_snapshot')
    return result


# Preserve the private cache-reset interface used by existing tests and refresh.
_package_identity.cache_clear = _package_snapshot.cache_clear


def executable_identity(path: Optional[str]) -> dict[str, Any]:
    resolved = shutil.which(path) if path else None
    if not resolved and path and Path(path).is_file():
        resolved = str(Path(path).resolve())
    return {'available': bool(resolved and os.access(resolved, os.X_OK)),
            'sha256': file_identity(resolved)}


def refresh() -> None:
    """Take a fresh installed-package snapshot at the start of a screen."""
    _package_identity.cache_clear()
    from . import protonate, openbabel_inputs
    protonate.obabel_path.cache_clear()
    openbabel_inputs.refresh()


def preparation_identity() -> dict[str, Any]:
    from . import protonate, openbabel_inputs
    try:
        rdkit_version = importlib.metadata.version('rdkit')
    except importlib.metadata.PackageNotFoundError:
        rdkit_version = None
    native = {}
    for module in ('rdkit.rdBase', 'rdkit.Chem.rdDistGeom', 'rdkit.Chem.rdForceFieldHelpers'):
        try:
            spec = importlib.util.find_spec(module)
            native[module] = file_identity(spec.origin if spec else None)
        except (ImportError, ValueError):
            native[module] = None
    directory = Path(__file__).parent
    rdkit_package = _package_identity('rdkit', 'rdkit')
    meeko_package = _package_identity('meeko', 'meeko')
    verified = all(p.get('snapshot_state') == 'observed' for p in (rdkit_package, meeko_package))
    return {'protocol': 'etkdgv3-seed42-smallrings-mmff94s-uff-meeko-v1',
            'rdkit_version': rdkit_version, 'rdkit_selected_native_sha256': native,
            'rdkit_package': rdkit_package, 'meeko': meeko_package,
            'package_inputs_state': 'observed' if verified else 'unverified',
            'openbabel': openbabel_inputs.identity(protonate.obabel_path()),
            'code_sha256': {name: file_identity(directory / name)
                            for name in ('ligprep.py', 'protonate.py', 'software.py', 'openbabel_inputs.py')}}


def receptor_preparation_identity() -> dict[str, Any]:
    """Installed preparation bytes, distinct from code and requested settings."""
    packages = {name: _package_identity(name, prefix) for name, prefix in (
        ('openmm', 'openmm'), ('pdbfixer', 'pdbfixer'), ('rdkit', 'rdkit'),
        ('meeko', 'meeko'), ('numpy', 'numpy'), ('openmmforcefields', 'openmmforcefields'))}
    required = ('openmm', 'pdbfixer', 'rdkit', 'meeko', 'numpy')
    inspected = (*required, *(name for name in ('openmmforcefields',)
                              if packages[name].get('available')))
    verified = all(packages[name].get('snapshot_state') == 'observed'
                   and not packages[name].get('missing_manifest_files') for name in inspected)
    return {'packages': packages, 'package_inputs_state': 'observed' if verified else 'unverified',
            'scope': 'installed package bytes; host BLAS/libm and loaded memory are not attested'}


def engine_identity(engine: str, binary: Optional[str] = None) -> dict[str, Any]:
    if engine == 'vina':
        implementation = dict(_package_identity('vina', 'vina'))
    elif engine == 'unidock':
        from . import unidock
        try:
            path = unidock.executable(binary)
        except unidock.UnidockUnavailable:
            path = None
        implementation = executable_identity(path)
    else:
        raise ValueError(f'unsupported actual engine {engine!r}')
    directory = Path(__file__).parent
    return {'engine': engine, 'implementation': implementation,
            'python': platform.python_version(), 'system': platform.system(),
            'machine': platform.machine(),
            'code_sha256': {name: file_identity(directory / name)
                            for name in ('engine.py', 'poses.py', 'unidock.py')},
            'preparation': preparation_identity()}
