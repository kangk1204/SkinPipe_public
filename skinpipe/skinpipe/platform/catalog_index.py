"""Pinned full-catalog similarity index; build separately from target jobs."""

from __future__ import annotations

import json
import time
from pathlib import Path

from .artifacts import sha256_file, canonical_json, atomic_bytes
from .identity import digest


class CatalogSimilarityIndex:
    def __init__(self, root):
        self.root = Path(root)
        self.manifest = json.loads((self.root / "manifest.json").read_bytes())
        if self.manifest.get("schema") != "catalog-similarity-index-v1":
            raise ValueError("Unsupported catalog similarity index")
        self._state = None
        self._engine = None
        self.verify()

    def verify(self):
        paths = [self.root / name for name in ("structures.h5", "lookup.duckdb")]
        state = tuple(
            (s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns, s.st_ctime_ns)
            for p in paths
            for s in [p.stat()]
        )
        if self._state != state:
            for p in paths:
                if sha256_file(p) != self.manifest["files"][p.name]["sha256"]:
                    raise ValueError("Catalog index checksum mismatch")
            after = tuple(
                (s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns, s.st_ctime_ns)
                for p in paths for s in [p.stat()]
            )
            if state != after:
                raise ValueError("Catalog index changed during verification")
            self._engine = None
            self._state = state

    def query(self, reference_smiles, *, limit=10000):
        import duckdb
        from rdkit import Chem
        from FPSim2 import FPSim2Engine
        from .catalog_query_lifetime import interruptible_query, query_checkpoint
        query_checkpoint()

        if (
            not isinstance(reference_smiles, list)
            or not 1 <= len(reference_smiles) <= 5
            or any(
                not isinstance(s, str) or len(s) > 4096 or Chem.MolFromSmiles(s) is None
                for s in reference_smiles
            )
        ):
            raise ValueError("At most five valid exact reference structures required")
        if (
            isinstance(limit, bool)
            or not isinstance(limit, int)
            or not 1 <= limit <= 10000
        ):
            raise ValueError("Similarity retrieval limit must be 1..10000")
        self.verify()
        query_checkpoint()
        if not self.manifest["record_count"]:
            return {
                "rows": [],
                "index_manifest_sha256": digest(self.manifest),
                "indexed_structures": 0,
                "scope": "whole pinned eligible catalog",
                "search_filters": self.manifest.get("search_filters", {}),
                "search_gate_manifest_sha256": self.manifest.get(
                    "search_gate_manifest_sha256"
                ),
            }
        if self._engine is None:
            # FPSim2's disk path keeps a 17M fingerprint matrix off the API
            # process heap. One worker is used for each bounded search.
            self._engine = FPSim2Engine(
                str(self.root / "structures.h5"), in_memory_fps=False
            )
        scores = {}
        for index, smiles in enumerate(reference_smiles):
            query_checkpoint()
            hits = self._engine.on_disk_top_k(smiles, limit, 0.0, n_workers=1)
            query_checkpoint()
            for hit in hits:
                query_checkpoint()
                mol_id, score = int(hit["mol_id"]), float(hit["coeff"])
                if mol_id not in scores or score > scores[mol_id][0]:
                    scores[mol_id] = (score, index)
        selected = sorted(scores, key=lambda mid: (-scores[mid][0], mid))[:limit]
        if not selected:
            return {
                "rows": [],
                "index_manifest_sha256": digest(self.manifest),
                "scope": "whole pinned indexed catalog",
            }
        with duckdb.connect(
            str(self.root / "lookup.duckdb"), read_only=True
        ) as connection:
            query_checkpoint()
            with interruptible_query(connection):
                try:
                    rows = {
                        r[0]: r[1]
                        for r in connection.execute(
                            "SELECT mol_id,compound_id FROM structures WHERE mol_id IN (SELECT unnest(?))",
                            [selected],
                        ).fetchall()
                    }
                except Exception:
                    query_checkpoint()
                    raise
            query_checkpoint()
        if set(rows) != set(selected):
            raise ValueError("Similarity index-to-structure mapping incomplete")
        self.verify()
        query_checkpoint()
        return {
            "rows": [
                {
                    "compound_id": rows[mid],
                    "tanimoto_2d": scores[mid][0],
                    "reference_index": scores[mid][1],
                }
                for mid in selected
            ],
            "index_manifest_sha256": digest(self.manifest),
            "scope": "whole pinned indexed catalog",
            "indexed_structures": self.manifest["record_count"],
            "threshold": 0.0,
            "source_report_sha256": self.manifest["source_report_sha256"],
            "search_filters": self.manifest.get("search_filters", {}),
            "search_gate_manifest_sha256": self.manifest.get(
                "search_gate_manifest_sha256"
            ),
        }

    def diversity(self, *, limit=1000, seed=0, algorithm="global-stratified-morgan-maxmin128-roundrobin-v1"):
        from .catalog_diversity import query_index
        return query_index(self, limit=limit, seed=seed, algorithm=algorithm)


def build_index(scout_root, output):
    """Single-process offline build. Invalid chemistry is logged, never erased."""
    import duckdb
    from rdkit import Chem, rdBase
    from FPSim2.io import create_db_file
    import pyarrow as pa
    from .scout import ScoutCatalog

    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    scout = ScoutCatalog(scout_root)
    status = scout.status()
    source_fingerprint = scout.structure_fingerprint()
    if (
        status["structure_admission"] != "reported_verified"
        or not status["report_sha256"]
    ):
        raise ValueError("Pinned admitted catalog report required before index build")
    source = duckdb.connect(str(scout.database), read_only=True)
    source.execute("SET threads=1; SET memory_limit='1GB'")
    lookup = duckdb.connect(str(output / "lookup.duckdb"))
    lookup.execute("SET threads=1; SET memory_limit='6GB'; SET preserve_insertion_order=false")
    lookup.execute(
        "CREATE TABLE structures(mol_id BIGINT PRIMARY KEY,compound_id VARCHAR UNIQUE)"
    )
    count, invalid = 0, 0
    started = time.monotonic()
    source_state = (scout.database.stat().st_size, scout.database.stat().st_mtime_ns)

    def records():
        nonlocal count, invalid
        cursor = source.execute(
            "SELECT compound_id,identity_smiles FROM compounds ORDER BY compound_id"
        )
        with (output / "unsupported-structures.jsonl").open("w") as failures:
            while True:
                batch = cursor.fetchmany(32000)
                if not batch:
                    break
                index_rows, valid = [], []
                for compound_id, smiles in batch:
                    mol = Chem.MolFromSmiles(smiles or "")
                    if mol is None:
                        invalid += 1
                        failures.write(
                            json.dumps(
                                {
                                    "compound_id": compound_id,
                                    "reason": "structure_parse_failed",
                                }
                            )
                            + "\n"
                        )
                        continue
                    count += 1
                    index_rows.append((count, compound_id))
                    valid.append((mol, count))
                if index_rows:
                    table = pa.table(
                        {
                            "mol_id": [r[0] for r in index_rows],
                            "compound_id": [r[1] for r in index_rows],
                        }
                    )
                    lookup.register("index_chunk", table)
                    lookup.execute("INSERT INTO structures SELECT * FROM index_chunk")
                    lookup.unregister("index_chunk")
                atomic_bytes(
                    output / "progress.json",
                    canonical_json(
                        {
                            "status": "building",
                            "indexed": count,
                            "unsupported": invalid,
                            "source_count": status["compound_count"],
                            "elapsed_seconds": time.monotonic() - started,
                        }
                    ),
                )
                yield from valid

    try:
        create_db_file(
            records(),
            str(output / "structures.h5"),
            "rdkit",
            "Morgan",
            fp_params={"radius": 2, "fpSize": 2048},
        )
        lookup.execute("CHECKPOINT")
    finally:
        source.close()
        lookup.close()
    scout.status()
    if scout.structure_fingerprint() != source_fingerprint or source_state != (
        scout.database.stat().st_size,
        scout.database.stat().st_mtime_ns,
    ):
        raise ValueError(
            "Source changed during full-catalog index build; index not activated"
        )
    manifest = {
        "schema": "catalog-similarity-index-v1",
        "record_count": count,
        "unsupported_count": invalid,
        "rdkit_version": rdBase.rdkitVersion,
        "source_report_sha256": status["report_sha256"],
        "source_structure_fingerprint": source_fingerprint,
        "elapsed_seconds": time.monotonic() - started,
        "fingerprint": {"type": "Morgan", "radius": 2, "fpSize": 2048},
        "files": {
            p.name: {"sha256": sha256_file(p), "size_bytes": p.stat().st_size}
            for p in output.iterdir()
            if p.is_file() and p.name != "progress.json"
        },
    }
    atomic_bytes(output / "manifest.json", canonical_json(manifest))
    atomic_bytes(
        output / "progress.json",
        canonical_json(
            {
                "status": "complete",
                "indexed": count,
                "unsupported": invalid,
                "elapsed_seconds": time.monotonic() - started,
            }
        ),
    )
    return manifest


def build_filtered_index(
    index_root, gate_root, output, *, hide_severe_admet=True, purchasable_only=True
):
    """Copy already computed fingerprints; never rerun chemistry or predictions."""
    import duckdb
    import numpy as np
    import tables as tb
    from .admet_search import ADMETSearchGate, _attach

    if (
        type(hide_severe_admet) is not bool
        or type(purchasable_only) is not bool
        or not (hide_severe_admet or purchasable_only)
    ):
        raise ValueError("At least one explicit search filter required")
    index, gate = CatalogSimilarityIndex(index_root), ADMETSearchGate(gate_root)
    if (
        index.manifest["source_structure_fingerprint"]
        != gate.manifest["source_structure_fingerprint"]
    ):
        raise ValueError("Eligibility and fingerprint source populations differ")
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    where = []
    if hide_severe_admet:
        where.append("NOT g.severe_admet")
    if purchasable_only:
        where.append("g.purchase_record")
    with duckdb.connect(str(output / "lookup.duckdb")) as conn:
        # A 16M-row ART identity index exceeds the 1GB budget used for a
        # bounded read query. Builds are separately owned and limited to 4GB.
        conn.execute(
            "SET threads=1; SET memory_limit='4GB'; SET preserve_insertion_order=false"
        )
        _attach(conn, index.root / "lookup.duckdb", "original")
        _attach(conn, gate.database, "gate")
        conn.execute(
            "CREATE TABLE structures AS SELECT s.* FROM original.structures s JOIN gate.catalog g USING(compound_id) WHERE "
            + " AND ".join(where)
        )
        conn.execute("CREATE UNIQUE INDEX mol_identity ON structures(mol_id)")
        conn.execute("CREATE UNIQUE INDEX exact_identity ON structures(compound_id)")
        max_id = (
            conn.execute("SELECT max(mol_id) FROM original.structures").fetchone()[0]
            or 0
        )
        mask = np.zeros(max_id + 1, dtype=bool)
        cursor = conn.execute("SELECT mol_id FROM structures")
        while rows := cursor.fetchmany(100000):
            mask[np.fromiter((r[0] for r in rows), dtype=np.int64)] = True
        expected = conn.execute("SELECT count(*) FROM structures").fetchone()[0]
        conn.execute("CHECKPOINT")
    count, bins = 0, []
    with (
        tb.open_file(str(index.root / "structures.h5"), mode="r") as original,
        tb.open_file(str(output / "structures.h5"), mode="w") as filtered,
    ):
        fps = filtered.create_table(
            filtered.root,
            "fps",
            original.root.fps.description,
            filters=original.root.fps.filters,
        )
        params = filtered.create_vlarray(filtered.root, "config", atom=tb.ObjectAtom())
        # FPSim2 0.7.4 stores type/parameters/RDKit/FPSim2 at 0..3;
        # population-specific popcount bins are index 4, not index 3.
        if len(original.root.config) != 5:
            raise ValueError("Unsupported FPSim2 config layout")
        for item in original.root.config[:4]:
            params.append(item)
        for start in range(0, original.root.fps.nrows, 64000):
            chunk = original.root.fps.read(
                start, min(start + 64000, original.root.fps.nrows)
            )
            chunk = chunk[mask[chunk["fp_id"]]]
            if len(chunk):
                fps.append(chunk)
                values, starts = np.unique(chunk["popcnt"], return_index=True)
                for i, value in enumerate(values):
                    first = count + int(starts[i])
                    last = count + (
                        int(starts[i + 1]) if i + 1 < len(starts) else len(chunk)
                    )
                    if bins and bins[-1][0] == int(value):
                        bins[-1] = (int(value), (bins[-1][1][0], last))
                    else:
                        bins.append((int(value), (first, last)))
                count += len(chunk)
            atomic_bytes(
                output / "progress.json",
                canonical_json(
                    {
                        "stage": "copy_filtered_fingerprints",
                        "source_processed": int(
                            min(start + 64000, original.root.fps.nrows)
                        ),
                        "retained": count,
                        "expected": expected,
                        "elapsed_seconds": time.monotonic() - started,
                    }
                ),
            )
        fps.cols.popcnt.create_index(kind="full")
        params.append(bins)
    if count != expected:
        raise ValueError("Filtered fingerprints and identities differ")
    index.verify()
    gate.verify()
    manifest = {
        **index.manifest,
        "record_count": count,
        "full_index_manifest_sha256": digest(index.manifest),
        "search_gate_manifest_sha256": digest(gate.manifest),
        "search_filters": {
            "hide_severe_admet": hide_severe_admet,
            "purchasable_only": purchasable_only,
        },
        "elapsed_seconds": time.monotonic() - started,
        "files": {
            p.name: {"sha256": sha256_file(p), "size_bytes": p.stat().st_size}
            for p in output.iterdir()
            if p.name in {"structures.h5", "lookup.duckdb"}
        },
    }
    atomic_bytes(output / "manifest.json", canonical_json(manifest))
    atomic_bytes(
        output / "progress.json",
        canonical_json(
            {
                "stage": "complete",
                "retained": count,
                "source_indexed": index.manifest["record_count"],
                "elapsed_seconds": time.monotonic() - started,
            }
        ),
    )
    return manifest
