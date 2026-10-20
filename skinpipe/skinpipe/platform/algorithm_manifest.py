"""Freeze algorithm inputs as bounded, content-addressed artifacts."""
from __future__ import annotations

import gzip
import hashlib
import io
import zlib

from .artifacts import ArtifactIntegrityError, canonical_json

COMPRESSION_THRESHOLD = 32 * 1024 * 1024
MAX_ENCODED_BYTES = 64 * 1024 * 1024
# Bound decoded and wire sizes separately while retaining original bytes.
MAX_MANIFEST_BYTES = 384 * 1024 * 1024


def freeze_manifest(manifest, artifacts):
    raw = canonical_json(manifest)
    if len(raw) > MAX_MANIFEST_BYTES:
        raise ValueError("Frozen algorithm manifest exceeds384MiB")
    original = artifacts.put_bytes(raw)
    if len(raw) <= COMPRESSION_THRESHOLD:
        return original, None
    encoded = gzip.compress(raw, compresslevel=6, mtime=0)
    if len(encoded) > MAX_ENCODED_BYTES:
        raise ValueError("Compressed algorithm manifest exceeds64MiB")
    blob = artifacts.put_bytes(encoded)
    return original, {"encoding": "gzip", "sha256": blob["sha256"],
                      "size_bytes": len(encoded), "decoded_size_bytes": len(raw),
                      "decoded_sha256": original["sha256"]}


def decode_manifest(source, fetch_verified):
    """Return original bytes after verifying the bounded encoded CAS object."""
    transport = source.get("manifest_transport")
    if transport is None:
        return fetch_verified(source.get("manifest_sha256"))
    if (not isinstance(transport, dict)
            or set(transport) != {"encoding", "sha256", "size_bytes", "decoded_size_bytes", "decoded_sha256"}
            or transport.get("encoding") != "gzip"
            or transport.get("decoded_sha256") != source.get("manifest_sha256")):
        raise ArtifactIntegrityError("Invalid frozen manifest transport binding")
    for key, maximum in (("size_bytes", MAX_ENCODED_BYTES), ("decoded_size_bytes", MAX_MANIFEST_BYTES)):
        value = transport.get(key)
        if type(value) is not int or not 0 < value <= maximum:
            raise ArtifactIntegrityError("Frozen manifest transport size outside bound")
    encoded = fetch_verified(transport["sha256"])
    if len(encoded) != transport["size_bytes"]:
        raise ArtifactIntegrityError("Frozen manifest encoded size mismatch")
    try:
        with gzip.GzipFile(fileobj=io.BytesIO(encoded), mode="rb") as stream:
            raw = stream.read(transport["decoded_size_bytes"] + 1)
    except (OSError, EOFError, zlib.error) as error:
        raise ArtifactIntegrityError("Invalid frozen manifest gzip bytes") from error
    if (len(raw) != transport["decoded_size_bytes"]
            or hashlib.sha256(raw).hexdigest() != source["manifest_sha256"]):
        raise ArtifactIntegrityError("Frozen manifest decoded bytes differ from original identity")
    return raw
