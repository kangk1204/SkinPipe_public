'Validate the selected structure and retain explicit calculation provenance.'
from __future__ import annotations

from contextlib import contextmanager
import ctypes
import errno
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import sys
from typing import Callable
import uuid


def _sync_dir(path: Path) -> None:
    fd = os.open(path, os.O_RDONLY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def _json_atomic(path: Path, data: dict) -> None:
    tmp = path.with_name(path.name + '.tmp')
    with tmp.open('w', encoding='utf-8') as f:
        json.dump(data, f, ensure_ascii=False, indent=2)
        f.write('\n')
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, path)
    _sync_dir(path.parent)


def _record(store: Path, data: dict) -> None:
    history = store / 'transactions'
    history.mkdir(exist_ok=True)
    _json_atomic(history / (data['id'] + '.json'), data)
    _json_atomic(store / 'journal.json', data)


def _sync_tree(root: Path) -> None:
    for base, dirs, files in os.walk(root, topdown=False):
        for name in files:
            p = Path(base) / name
            if p.is_symlink():
                continue
            fd = os.open(p, os.O_RDONLY)
            try:
                os.fsync(fd)
            finally:
                os.close(fd)
        _sync_dir(Path(base))


def _inventory(root: Path) -> list[dict]:
    rows = []
    for base, dirs, files in os.walk(root):
        dirs.sort()
        for name in sorted(files):
            p = Path(base) / name
            if p == root / '_publication.json':
                continue
            if p.is_symlink():
                raise ValueError('publication artifacts must be self-contained files')
            h = hashlib.sha256()
            with p.open('rb') as f:
                for b in iter(lambda: f.read(1024 * 1024), b''):
                    h.update(b)
            rows.append({'path': p.relative_to(root).as_posix(),
                         'bytes': p.stat().st_size, 'sha256': h.hexdigest()})
        if any((Path(base) / d).is_symlink() for d in dirs):
            raise ValueError('publication directories must be self-contained')
    return rows


def verify_generation(root: Path) -> None:
    """Verify the complete artifact set, including additions and removals."""
    card = json.loads((root / '_publication.json').read_text(encoding='utf-8'))
    if card.get('phase') != 'complete' or _inventory(root) != card.get('files'):
        raise ValueError('published generation failed its content manifest')


def _exchange(a: Path, b: Path) -> None:
    """Atomically adopt a nonempty legacy directory without an absent path.

    Unsupported platforms fail before touching the destination. Regular
    generation-to-generation commits use portable POSIX rename instead.
    """
    libc = ctypes.CDLL(None, use_errno=True)
    if sys.platform.startswith('linux') and hasattr(libc, 'renameat2'):
        fn = libc.renameat2
        fn.argtypes = [ctypes.c_int, ctypes.c_char_p, ctypes.c_int,
                       ctypes.c_char_p, ctypes.c_uint]
        rc = fn(-100, os.fsencode(a), -100, os.fsencode(b), 2)
    elif sys.platform == 'darwin' and hasattr(libc, 'renamex_np'):
        fn = libc.renamex_np
        fn.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_uint]
        rc = fn(os.fsencode(a), os.fsencode(b), 2)
    else:
        raise OSError(errno.ENOTSUP, 'atomic legacy-directory adoption unavailable')
    if rc != 0:
        code = ctypes.get_errno()
        raise OSError(code, os.strerror(code))


def _store_for(dst: Path) -> Path:
    return dst.parent / ('.' + dst.name + '.generations')


def _destination(path: Path) -> Path:
    """Canonicalize the parent while retaining the replaceable final alias.

    macOS temp roots and user workspaces can contain parent symlinks. Resolving
    only the parent gives journal paths the same spelling as resolved readers
    without following an existing destination into its active generation.
    """
    absolute = Path(os.path.abspath(path))
    return absolute.parent.resolve() / absolute.name


@contextmanager
def _locked(dst: Path):
    dst.parent.mkdir(parents=True, exist_ok=True)
    with (dst.parent / ('.' + dst.name + '.publish.lock')).open('a') as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(f, fcntl.LOCK_UN)


def _recover_locked(dst: Path, store: Path) -> dict | None:
    journal = store / 'journal.json'
    if not journal.exists():
        return None
    data = json.loads(journal.read_text(encoding='utf-8'))
    if not re.fullmatch(r'[0-9a-f]{32}', str(data.get('id', ''))):
        raise ValueError('invalid publication journal identifier')
    generation = store / data['id']
    temporary = dst.parent / ('.' + dst.name + '.pointer.' + data['id'])
    previous = _destination(Path(data['previous'])) if data.get('previous') else None
    if previous is not None and previous.parent != store:
        raise ValueError('publication journal points outside its generation store')
    if (data.get('phase') == 'rolled_back' and dst.is_symlink()
            and dst.resolve(strict=True) == previous):
        return data
    # The atomic exchange may have committed but the process died before the
    # old directory was moved from the temporary pointer's name to its archive.
    active = dst.is_symlink() and dst.resolve(strict=True) == generation
    if (data.get('phase') == 'committed' and active
            and not (temporary.exists() and not temporary.is_symlink())):
        return data
    if data.get('phase') == 'inactive' and not active:
        return data
    if active and temporary.exists() and not temporary.is_symlink():
        if previous is None or previous.parent != store or previous.exists():
            raise RuntimeError('ambiguous legacy recovery; evidence retained')
        os.replace(temporary, previous)
        _sync_dir(dst.parent)
        _sync_dir(store)
    data['phase'] = 'committed' if active else 'inactive'
    _record(store, data)
    return data


def recover_publication(destination: Path) -> dict | None:
    dst = _destination(destination)
    with _locked(dst):
        store = _store_for(dst)
        if not store.exists():
            return None
        return _recover_locked(dst, store)


@contextmanager
def read_campaign(coverage: Path):
    """Pin a reader, including during initial adoption of a legacy directory.

    The shared lock is held until all files are consumed. This also closes the
    resolve/open race and protects readers whose legacy index names the alias.
    """
    index = Path(os.path.abspath(coverage))
    dst = _destination(index.parent)
    index = dst / index.name
    if not dst.parent.is_dir():
        yield index.resolve()
        return
    lock = dst.parent / ('.' + dst.name + '.publish.lock')
    # An unpacked read-only artifact cannot be published in place by this
    # process. Existing managed lock files can still be locked read-only.
    try:
        f = lock.open('a')
    except PermissionError:
        if lock.exists():
            f = lock.open('r')
        else:
            yield index.resolve()
            return
    with f:
        fcntl.flock(f, fcntl.LOCK_SH)
        try:
            yield index.resolve()
        finally:
            fcntl.flock(f, fcntl.LOCK_UN)


def ensure_mutable_workdir(path: Path) -> None:
    """Refuse an in-place compute/retry/report write into a published artifact."""
    root = Path(path).resolve()
    for p in (root, *root.parents):
        store = p.parent
        if ((p / '_publication.json').is_file()
                or (store.name.startswith('.') and store.name.endswith('.generations')
                    and (store / 'journal.json').is_file())):
            raise ValueError('published generations are immutable; use a new staging workdir')


def publish_generation(destination: Path,
                       build: Callable[[Path, Path | None], dict], *,
                       require_new: bool = False) -> dict:
    """Build, durably validate, and commit one complete generation.

    ``build`` runs under the publication lock and receives the new root and
    previous root. On any failure all old data, new data and the journal remain
    available. A previous generation is never edited or removed here.
    ``require_new`` refuses existing output, including a destination created
    during the build. Installing a hard link to the new symlink is atomic and
    cannot replace an unrelated file, directory or symlink.
    """
    dst = _destination(destination)
    with _locked(dst):
        if require_new and os.path.lexists(dst):
            raise ValueError('new publication destination already exists')
        store = _store_for(dst)
        store.mkdir(exist_ok=True)
        _recover_locked(dst, store)
        if dst.is_symlink():
            current = dst.resolve(strict=True)
            if current.parent != store:
                raise ValueError('destination symlink is not a managed generation')
            if (current / '_publication.json').is_file():
                verify_generation(current)
            previous = current
            legacy = False
        elif dst.exists():
            if not dst.is_dir():
                raise ValueError('destination must be a directory or managed pointer')
            current = dst
            previous = store / ('legacy-' + uuid.uuid4().hex)
            legacy = True
        else:
            current = previous = None
            legacy = False
        ident = uuid.uuid4().hex
        root = store / ident
        pointer = dst.parent / ('.' + dst.name + '.pointer.' + ident)
        state = {'schema_version': 1, 'id': ident, 'phase': 'building',
                 'generation': str(root), 'previous': str(previous) if previous else None,
                 'legacy_adoption': legacy}
        _record(store, state)
        try:
            root.mkdir()
            details = build(root, current)
            _json_atomic(root / '_publication.json', state | {
                'phase': 'complete', 'files': _inventory(root)})
            _sync_tree(root)
            verify_generation(root)
            _sync_dir(store)
            state['phase'] = 'prepared'
            _record(store, state)
            os.symlink(os.path.relpath(root, dst.parent), pointer)
            _sync_dir(dst.parent)
            if legacy:
                _exchange(pointer, dst)
                _sync_dir(dst.parent)
                os.replace(pointer, previous)
                _sync_dir(store)
            elif require_new:
                os.link(pointer, dst, follow_symlinks=False)
                pointer.unlink()
            else:
                os.replace(pointer, dst)
            _sync_dir(dst.parent)
            state['phase'] = 'committed'
            _record(store, state)
            return state | details
        except BaseException as exc:
            state['error'] = type(exc).__name__ + ': ' + str(exc)
            # A failed fsync/journal write must not mask the original failure
            # or remove the last durable checkpoint needed by recovery.
            try:
                _record(store, state)
            except OSError:
                pass
            raise


def rollback_publication(destination: Path) -> dict:
    """Activate the recorded previous generation; retain both generations."""
    dst = _destination(destination)
    with _locked(dst):
        store = _store_for(dst)
        state = _recover_locked(dst, store)
        if state and state['phase'] == 'rolled_back':
            return state
        if not state or state['phase'] != 'committed' or not state.get('previous'):
            raise ValueError('no committed publication with a previous generation')
        prior = Path(state['previous'])
        if prior.parent != store or not prior.is_dir():
            raise ValueError('previous generation is unavailable')
        pointer = dst.parent / ('.' + dst.name + '.rollback.' + uuid.uuid4().hex)
        os.symlink(os.path.relpath(prior, dst.parent), pointer)
        os.replace(pointer, dst)
        _sync_dir(dst.parent)
        state['phase'] = 'rolled_back'
        _record(store, state)
        return state
