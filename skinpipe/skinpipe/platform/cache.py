"""Exact-identity prediction cache; parent identities never substitute stereochemistry."""
from __future__ import annotations

from contextlib import closing
import hashlib
import json
from pathlib import Path
import sqlite3
from datetime import datetime, timezone

from .artifacts import ArtifactIntegrityError, canonical_json


def prediction_identity(*, exact_structure, model, weights_sha256, preprocessing, **runtime):
    result = {"kind": "prediction", "exact_structure": exact_structure, "model": model,
              "weights_sha256": weights_sha256, "preprocessing": preprocessing, "runtime": runtime}
    _validate(result)
    return result


def docking_identity(*, prepared_ligand_sha256, receptor_sha256, pocket, engine, settings, seed, **runtime):
    result = {"kind": "docking", "prepared_ligand_sha256": prepared_ligand_sha256,
              "receptor_sha256": receptor_sha256, "pocket": pocket, "engine": engine,
              "settings": settings, "seed": seed, "runtime": runtime}
    _validate(result)
    return result


def _validate(identity):
    required = {"prediction": ("exact_structure", "model", "weights_sha256", "preprocessing"),
                "docking": ("prepared_ligand_sha256", "receptor_sha256", "pocket", "engine", "settings", "seed")}
    if not isinstance(identity, dict) or identity.get("kind") not in required:
        raise ValueError("A versioned prediction or docking identity is required")
    for key in required[identity["kind"]]:
        if key not in identity or identity[key] is None or identity[key] == "":
            raise ValueError("Missing cache identity: " + key)
    canonical_json(identity)


class PredictionCache:
    def __init__(self, path):
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with closing(self._connect()) as connection, connection:
            connection.execute("PRAGMA journal_mode=WAL")
            connection.execute("""CREATE TABLE IF NOT EXISTS cache (
                key TEXT PRIMARY KEY, identity TEXT NOT NULL, payload TEXT NOT NULL,
                payload_sha256 TEXT NOT NULL, created_at TEXT NOT NULL)""")

    def _connect(self):
        return sqlite3.connect(self.path, timeout=30)

    @staticmethod
    def key(identity):
        _validate(identity)
        return hashlib.sha256(canonical_json(identity)).hexdigest()

    def put(self, identity, payload):
        key = self.key(identity)
        encoded = canonical_json(payload)
        row = (key, canonical_json(identity).decode(), encoded.decode(), hashlib.sha256(encoded).hexdigest(),
               datetime.now(timezone.utc).isoformat())
        with closing(self._connect()) as connection, connection:
            connection.execute("BEGIN IMMEDIATE")
            existing = connection.execute("SELECT payload_sha256 FROM cache WHERE key=?", (key,)).fetchone()
            if existing and existing[0] != row[3]:
                raise ValueError("Conflicting output for an immutable cache identity; version the computation")
            connection.execute("INSERT OR IGNORE INTO cache VALUES(?,?,?,?,?)", row)
        return key

    def get(self, identity):
        key = self.key(identity)
        with closing(self._connect()) as connection:
            row = connection.execute("SELECT identity,payload,payload_sha256 FROM cache WHERE key=?", (key,)).fetchone()
        if row is None:
            return None
        if row[0] != canonical_json(identity).decode() or hashlib.sha256(row[1].encode()).hexdigest() != row[2]:
            raise ArtifactIntegrityError("Corrupt cache entry: " + key)
        return json.loads(row[1])
