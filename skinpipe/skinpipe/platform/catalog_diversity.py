"""Bounded diversity from a globally spread, sealed fingerprint population.

No source/index writes, arbitrary catalog window or target-activity inference.
The sampled pool, not the entire catalog, is chemically diversified.
"""
from collections import deque

from .identity import digest


ALGORITHM = "global-stratified-morgan-maxmin128-roundrobin-v1"
POOL_LIMIT = 20000
MAXMIN_PIVOTS = 128
STRATA_LIMIT = 128


def selection_policy():
    """Server-only default, frozen only on newly submitted generalized runs."""
    return {"version": "protein-seedless-diversity-v1", "algorithm": ALGORITHM,
            "seed": 0, "pool_limit": POOL_LIMIT, "maxmin_pivots": MAXMIN_PIVOTS}


def selection_summary(receipt, *, selected_count=0, status="complete", reason=None):
    return {"fallback_used": True, "strategy": ALGORITHM, "status": status,
            "source_population": receipt.get("indexed_structures"),
            "source_pool_size": receipt.get("source_pool_size"), "selected_count": selected_count,
            "requested_limit": receipt.get("requested_limit"), "pool_limit": POOL_LIMIT,
            "seed": receipt.get("seed", 0), "index_manifest_sha256": receipt.get("index_manifest_sha256"),
            "scope": "bounded chemically diverse subset from globally stratified pinned index; not exhaustive or target-active",
            "limitations": list(receipt.get("limitations", [])),
            "target_activity": "NOT_ESTABLISHED_NO_REFERENCE", **({"reason": reason} if reason else {})}


def validate_request(*, limit=1000, seed=0, algorithm=ALGORITHM):
    if type(limit) is not int or not 1 <= limit <= 10000:
        raise ValueError("Diversity retrieval limit must be 1..10000")
    if type(seed) is not int or not 0 <= seed < 2**32:
        raise ValueError("Diversity seed must be an unsigned 32-bit integer")
    if algorithm != ALGORITHM:
        raise ValueError("Unsupported diversity algorithm")


def sample_blocks(population, pool_size, seed, manifest_sha256):
    """One seed-selected contiguous block in every nonoverlapping global stratum.

    Blocks bound compressed-HDF5 I/O; random individual reads across all 17M
    rows could otherwise read almost the full fingerprint file.
    """
    if population <= pool_size:
        return [(start, min(1024, population - start)) for start in range(0, population, 1024)]
    strata = min(STRATA_LIMIT, pool_size)
    blocks = []
    for i in range(strata):
        first, last = i * population // strata, (i + 1) * population // strata
        size = pool_size // strata + int(i < pool_size % strata)
        offset = int(digest([ALGORITHM, seed, manifest_sha256, i]), 16) % (last - first - size + 1)
        blocks.append((first + offset, size))
    return blocks


def select_fingerprints(fingerprints, mol_ids, limit, seed, *, checkpoint):
    """MaxMin pivots followed by a balanced round-robin of nearest-pivot groups.

    At most 128 bulk similarity rounds even for a 10,000-candidate budget.
    Distances are Morgan fingerprint diversity, never target similarities.
    """
    import numpy as np
    from rdkit import DataStructs

    if not fingerprints:
        return []
    total = len(fingerprints)
    rounds = min(limit, MAXMIN_PIVOTS, total)
    distances = np.full(total, np.inf)
    nearest = np.zeros(total, dtype=np.int32)
    chosen = np.zeros(total, dtype=bool)
    selected = []
    pivot = int(digest([ALGORITHM, seed, mol_ids]), 16) % total
    for group in range(rounds):
        checkpoint()
        chosen[pivot] = True
        selected.append(pivot)
        observed = 1 - np.asarray(DataStructs.BulkTanimotoSimilarity(fingerprints[pivot], fingerprints))
        closer = observed < distances
        nearest[closer] = group
        distances = np.minimum(distances, observed)
        distances[chosen] = -1
        # np.argmax is stable for ties in the frozen sampled-row order.
        pivot = int(np.argmax(distances))
    groups = [deque() for _ in range(rounds)]
    for i in sorted(range(total), key=lambda i: (-distances[i], mol_ids[i])):
        if not chosen[i]:
            groups[int(nearest[i])].append(i)
    while len(selected) < min(limit, total):
        checkpoint()
        for group in groups:
            if group and len(selected) < limit:
                selected.append(group.popleft())
    return selected


def query_index(index, *, limit=1000, seed=0, algorithm=ALGORITHM):
    import duckdb
    import numpy as np
    import tables
    from rdkit import DataStructs, rdBase
    from .catalog_query_lifetime import query_checkpoint, interruptible_query

    validate_request(limit=limit, seed=seed, algorithm=algorithm)
    query_checkpoint()
    index.verify()
    manifest = index.manifest
    if manifest.get("fingerprint") != {"type": "Morgan", "radius": 2, "fpSize": 2048}:
        raise ValueError("Diversity requires the admitted 2048-bit Morgan index")
    population = manifest.get("record_count")
    if type(population) is not int or population < 0:
        raise ValueError("Invalid diversity index population")
    manifest_sha = digest(manifest)
    pool_size = min(population, max(limit, min(POOL_LIMIT, max(256, limit * 4))))
    blocks = sample_blocks(population, pool_size, seed, manifest_sha)
    fingerprints, mol_ids = [], []
    query_checkpoint()
    with tables.open_file(str(index.root / "structures.h5"), mode="r") as source:
        table = source.root.fps
        words = ["f" + str(i) for i in range(1, 33)]
        if table.nrows != population or set(table.colnames) != {"fp_id", "popcnt", *words}:
            raise ValueError("Fingerprint layout or population differs from the sealed index")
        for start, size in blocks:
            query_checkpoint()
            rows = table.read(start, start + size)
            for row in rows:
                query_checkpoint()
                fp = DataStructs.CreateFromBinaryText(np.asarray([row[name] for name in words], dtype="<u8").tobytes())
                if fp.GetNumOnBits() != int(row["popcnt"]):
                    raise ValueError("Fingerprint popcount differs from admitted row")
                fingerprints.append(fp)
                mol_ids.append(int(row["fp_id"]))
    if len(mol_ids) != pool_size or len(set(mol_ids)) != len(mol_ids):
        raise ValueError("Diversity sampled index identities are incomplete or duplicated")
    selected = select_fingerprints(fingerprints, mol_ids, limit, seed, checkpoint=query_checkpoint)
    ids = [mol_ids[i] for i in selected]
    mapping = {}
    if ids:
        with duckdb.connect(str(index.root / "lookup.duckdb"), read_only=True) as connection:
            query_checkpoint()
            with interruptible_query(connection):
                try:
                    mapped = connection.execute(
                        "SELECT mol_id,compound_id FROM structures WHERE mol_id IN (SELECT unnest(?))", [ids]).fetchall()
                except Exception:
                    query_checkpoint()
                    raise
            mapping = dict(mapped)
        if len(mapped) != len(ids) or set(mapping) != set(ids) or len(set(mapping.values())) != len(ids):
            raise ValueError("Diversity index-to-structure mapping incomplete")
    query_checkpoint()
    index.verify()
    return {"rows": [{"compound_id": mapping[mid], "diversity_rank": rank}
                     for rank, mid in enumerate(ids, 1)],
            "strategy": ALGORITHM, "seed": seed, "requested_limit": limit,
            "pool_limit": POOL_LIMIT, "maxmin_pivots": min(limit, MAXMIN_PIVOTS, pool_size),
            "indexed_structures": population, "source_pool_size": pool_size,
            "sampled_blocks": [{"start_row": first, "count": size} for first, size in blocks],
            "sampled_mol_ids_sha256": digest(mol_ids), "selected_mol_ids_sha256": digest(ids),
            "index_manifest_sha256": manifest_sha,
            "source_report_sha256": manifest["source_report_sha256"],
            "source_structure_fingerprint": manifest["source_structure_fingerprint"],
            "search_filters": manifest.get("search_filters", {}),
            "search_gate_manifest_sha256": manifest.get("search_gate_manifest_sha256"),
            "fingerprint": manifest["fingerprint"], "selection_runtime_rdkit": rdBase.rdkitVersion,
            "scope": "bounded chemically diverse subset from globally stratified pinned index",
            "limitations": ["Only the globally spread bounded pool was diversified; not exhaustive MaxMin of the entire catalog",
                            "Fingerprint diversity does not establish target binding, activity, action direction, skin efficacy or safety"]}
