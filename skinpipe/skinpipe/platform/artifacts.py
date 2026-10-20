"""Atomic, content addressed artifacts with bounded corruption isolation."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import re
import tempfile


def canonical_json(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def sha256_file(path):
    # Cold catalog seals can span tens of GB. Check the current query/study
    # lifetime between existing blocks without changing hash bytes or order.
    # The independently staged model runner includes only its frozen minimal
    # modules. An active cancellation scope necessarily loaded its helper;
    # inspect that module without importing optional queue/engine dependencies
    # into the model interpreter.
    import sys
    query_checkpoint = getattr(sys.modules.get("skinpipe.platform.catalog_query_lifetime"), "query_checkpoint", lambda: None)
    catalog_checkpoint = getattr(sys.modules.get("skinpipe.platform.catalog_cancellation"), "catalog_checkpoint", lambda: None)
    query_checkpoint()
    catalog_checkpoint()
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            query_checkpoint()
            catalog_checkpoint()
            digest.update(block)
    query_checkpoint()
    catalog_checkpoint()
    return digest.hexdigest()


def atomic_bytes(path, data):
    """Replace one file only after all bytes have reached the filesystem."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=".writing-", dir=path.parent)
    try:
        with os.fdopen(fd, "wb") as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
        fsync_directory(path.parent)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def fsync_directory(path):
    """Persist the published name as well as its file contents on POSIX."""
    fd = os.open(path, os.O_RDONLY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


class ArtifactIntegrityError(ValueError):
    """The requested artifact is unusable; unrelated artifacts remain valid."""


class ArtifactStore:
    def __init__(self, root):
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)

    def _paths(self, digest):
        if not isinstance(digest, str) or not re.fullmatch(r"[0-9a-f]{64}", digest):
            raise ValueError("Expected a SHA256 artifact identifier")
        base = self.root / digest[:2] / digest
        return base.with_suffix(".blob"), base.with_suffix(".json")

    def put_bytes(self, data, *, media_type="application/octet-stream", metadata=None):
        data = bytes(data)
        digest = hashlib.sha256(data).hexdigest()
        blob, manifest = self._paths(digest)
        receipt = {"sha256": digest, "size_bytes": len(data), "media_type": media_type}
        if metadata is not None:
            receipt["metadata"] = metadata
        # Metadata belongs to the reference, not to the content identity. Store
        # only immutable content facts so identical bytes can have many uses.
        if not self._already_stored(digest, len(data)):
            atomic_bytes(blob, data)
            atomic_bytes(manifest, canonical_json({"sha256": digest, "size_bytes": len(data)}))
        return receipt

    def _already_stored(self, digest, size):
        # Existence alone is not evidence. A reused object still undergoes the
        # full content/receipt verification; only redundant durable writes go.
        try:
            return self.verified_path(digest).stat().st_size == size
        except ArtifactIntegrityError:
            return False

    def put_file(self, path, *, media_type="application/octet-stream", metadata=None):
        """Stream large poses/uploads without loading the entire file into RAM."""
        path = Path(path)
        digest = sha256_file(path)
        blob, manifest = self._paths(digest)
        size = path.stat().st_size
        if self._already_stored(digest, size):
            if sha256_file(path) != digest:
                raise ArtifactIntegrityError("Source changed while publishing artifact")
            return {"sha256": digest, "size_bytes": size, "media_type": media_type, "metadata": metadata or {}}
        blob.parent.mkdir(parents=True, exist_ok=True)
        fd, temporary = tempfile.mkstemp(prefix=".writing-", dir=blob.parent)
        try:
            with path.open("rb") as source, os.fdopen(fd, "wb") as out:
                for block in iter(lambda: source.read(1024 * 1024), b""):
                    out.write(block)
                out.flush()
                os.fsync(out.fileno())
            if sha256_file(temporary) != digest:
                raise ArtifactIntegrityError("Source changed while publishing artifact")
            size = os.path.getsize(temporary)
            os.replace(temporary, blob)
            fsync_directory(blob.parent)
            atomic_bytes(manifest, canonical_json({"sha256": digest, "size_bytes": size}))
        finally:
            if os.path.exists(temporary):
                os.unlink(temporary)
        return {"sha256": digest, "size_bytes": size, "media_type": media_type, "metadata": metadata or {}}

    def verified_path(self, digest):
        blob, manifest = self._paths(digest)
        try:
            receipt = json.loads(manifest.read_bytes())
            if (not isinstance(receipt, dict) or receipt.get("sha256") != digest or blob.stat().st_size != receipt.get("size_bytes")
                    or sha256_file(blob) != digest):
                raise ArtifactIntegrityError("Artifact checksum or receipt mismatch: " + digest)
        except (OSError, json.JSONDecodeError) as error:
            raise ArtifactIntegrityError("Artifact missing or unreadable: " + digest) from error
        return blob

    def get_bytes(self, digest):
        path = self.verified_path(digest)
        data = path.read_bytes()
        if hashlib.sha256(data).hexdigest() != digest:
            raise ArtifactIntegrityError("Artifact changed during read: " + digest)
        return data
