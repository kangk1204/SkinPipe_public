"""Reference-guided retrieval from a frozen catalog, with per-row provenance.

No BRICS generation, potency transfer or universal similarity exclusion. All
source rows survive, including unsupported chemistry and unfinished alignments.
"""

from __future__ import annotations

from copy import deepcopy
import time

from .identity import chemical_identity, digest
from .ranking import finite


def acyclic_topology(smiles):
    """Compare linker/functional-group graphs while ignoring chain homologs.

    Saturated neutral degree-two carbon runs are contracted. This is a
    separate exploratory descriptor, never a redefinition of Murcko cores.
    """
    from rdkit import Chem

    identity = chemical_identity(smiles)
    if (identity["structure_status"] != "valid" or identity["parent_ambiguous"]
            or identity["scaffold_kind"] != "acyclic"):
        return None
    mol = Chem.RWMol(Chem.MolFromSmiles(identity["parent_smiles"]))
    while True:
        removable = next((a for a in mol.GetAtoms()
            if a.GetAtomicNum() == 6 and a.GetDegree() == 2
            and a.GetFormalCharge() == 0 and a.GetIsotope() == 0
            and a.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED
            and all(b.GetBondType() == Chem.BondType.SINGLE for b in a.GetBonds())), None)
        if removable is None:
            break
        index = removable.GetIdx()
        a, b = [n.GetIdx() for n in removable.GetNeighbors()]
        mol.AddBond(a, b, Chem.BondType.SINGLE)
        mol.RemoveAtom(index)
    value = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
    return {"policy": "acyclic-linker-topology-v1", "smiles": value,
            "id": "acyclic-topology:" + digest(value),
            "limitation": "Contracted saturated carbon runs ignore chain-length homologs; topology novelty is not activity or efficacy evidence"}


HOP_ROUTES = {"reference_catalog_hop", "reference_acyclic_topology_hop"}


class _ShapeFeatures:
    """Request-local feature definitions and fixed reference features only."""

    def __init__(self):
        self.state, self.factory, self.references = None, None, {}

    def features(self, mol, conf, *, reference=False):
        from pathlib import Path
        from rdkit import RDConfig
        from rdkit.Chem import ChemicalFeatures

        path = Path(RDConfig.RDDataDir) / "BaseFeatures.fdef"
        stat = path.stat()
        state = (str(path.resolve()), stat.st_dev, stat.st_ino, stat.st_size,
                 stat.st_mtime_ns, stat.st_ctime_ns)
        if state != self.state:
            self.factory = ChemicalFeatures.BuildFeatureFactory(str(path))
            self.references.clear()
            self.state = state
        key = (id(mol), conf)
        cached = self.references.get(key) if reference else None
        if cached is not None and cached[0] is mol:
            return deepcopy(cached[1])
        values = [{"family": f.GetFamily(), "atoms": list(f.GetAtomIds()),
                   "position": list(f.GetPos(conf))}
                  for f in self.factory.GetFeaturesForMol(mol)]
        if reference:
            self.references[key] = (mol, deepcopy(values))
        return values


def _ensemble(smiles, conformers=10, timeout_seconds=0):
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None or mol.GetNumHeavyAtoms() > 256:
        raise ValueError("unsupported_shape_chemistry")
    if mol.GetStereoGroups():
        raise ValueError("enhanced_stereochemistry_requires_explicit_isomer")
    stereo = chemical_identity(smiles)
    if stereo.get("unassigned_stereo") or stereo.get("unassigned_bond_stereo"):
        raise ValueError("undefined_stereochemistry_requires_explicit_isomer")
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed, params.numThreads, params.maxIterations = 42, 1, 1000
    if timeout_seconds and hasattr(params, "timeout"):
        params.timeout = timeout_seconds
    ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=conformers, params=params))
    if not ids:
        raise ValueError("conformer_generation_failed")
    return mol, ids


def shape_compare(candidate_smiles, references, artifacts, *, conformers=10, stop=None,
                  _feature_cache=None, conformer_timeout_seconds=0):
    """Color-aware shape alignment, on copies of generated comparison structures.

    Conformers are bounded per molecule. One representative reference conformer
    per seed is used; the full candidate ensemble is tested. Saved comparison
    coordinates are never described as docking coordinates.
    """
    from rdkit import Chem
    from rdkit.Chem import rdShapeAlign

    probe, ids = _ensemble(candidate_smiles, conformers, conformer_timeout_seconds)
    best = None
    for reference, mol, ref_ids in references:
        for cid in ids:
            if stop:
                stop()
            copy = Chem.Mol(probe)
            shape, color = rdShapeAlign.AlignMol(
                mol,
                copy,
                refConfId=ref_ids[0],
                probeConfId=cid,
                useColors=True,
                opt_param=0.5,
            )
            key = ((shape + color) / 2, shape, color, reference["candidate_id"], -cid)
            if best is None or key > best[0]:
                best = key, copy, cid, reference, mol, ref_ids[0]
    if best is None:
        raise ValueError("reference_conformers_unavailable")
    key, mol, cid, ref, ref_mol, ref_cid = best
    features = _feature_cache or _ShapeFeatures()
    a = features.features(ref_mol, ref_cid, reference=True)
    b = features.features(mol, cid)
    matches = sorted(
        (sum((x - y) ** 2 for x, y in zip(r["position"], p["position"])) ** 0.5, i, j)
        for i, r in enumerate(a)
        for j, p in enumerate(b)
        if r["family"] == p["family"]
    )
    used_a, used_b, pairs = set(), set(), []
    for distance, i, j in matches:
        if distance > 1.5 or i in used_a or j in used_b:
            continue
        used_a.add(i)
        used_b.add(j)
        pairs.append(
            {
                "family": a[i]["family"],
                "reference_atoms": a[i]["atoms"],
                "candidate_atoms": b[j]["atoms"],
                "distance_angstrom": distance,
            }
        )

    def saved(m, conf):
        raw = (Chem.MolToMolBlock(m, confId=conf) + "\n$$$$\n").encode()
        artifact = artifacts.put_bytes(raw, media_type="chemical/x-mdl-sdfile")
        reread = Chem.MolFromMolBlock(raw.decode().split("$$$$")[0], removeHs=False)
        if reread is None or reread.GetNumAtoms() != m.GetNumAtoms():
            raise ValueError("comparison_atom_mapping_failed")
        return artifact["sha256"]

    return {
        "policy": "catalog-color-shape-v1",
        "reference_id": ref["candidate_id"],
        "shape_score": key[1],
        "pharmacophore_color_score": key[2],
        "retention_score": key[0],
        "candidate_sdf_sha256": saved(mol, cid),
        "reference_sdf_sha256": saved(ref_mol, ref_cid),
        "feature_pairs": pairs,
        "unmatched_reference_features": [v for i, v in enumerate(a) if i not in used_a],
        "unmatched_candidate_features": [v for i, v in enumerate(b) if i not in used_b],
        "atom_index_basis": "zero-based atom indices of these exact saved comparison SDFs including H",
        "coordinates_scope": "generated candidate alignment against the exact observed bound-ligand coordinates; original protein and ligand coordinates unchanged"
        if ref.get("observed_sdf_sha256") else "generated alignment only; original docking poses unchanged",
        "reference_coordinate_source": "observed_bound_ligand" if ref.get("observed_sdf_sha256") else "generated_conformer",
        "reference_observed_sdf_sha256": ref.get("observed_sdf_sha256"),
        "conformers_requested": conformers,
        "candidate_conformers_observed": len(ids),
        "reference_conformers_compared_per_seed": 1,
        "status": "complete",
    }


def catalog_hopping(
    rows,
    reference_set,
    artifacts,
    *,
    shape_limit=2000,
    conformers=10,
    stop=None,
    progress=None,
    shape_seconds=3600,
    conformer_timeout_seconds=0,
):
    from rdkit import Chem, DataStructs
    from rdkit.Chem import rdFingerprintGenerator

    rows = deepcopy(rows)
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    reference_rows = reference_set["population"]
    cores = {r["identity"]["scaffold_id"] for r in reference_rows}
    acyclic_cores = {v["id"] for r in reference_rows
                    if (v := acyclic_topology(r.get("smiles"))) is not None}
    fingerprints, ensembles, reference_conformer_holds = [], [], []
    for ref in reference_set["computational_seeds"]:
        if stop:
            stop()
        m = Chem.MolFromSmiles(ref["smiles"])
        fingerprints.append((ref, generator.GetFingerprint(m)))
        if not shape_limit:
            continue
        try:
            if ref.get("observed_sdf_sha256"):
                from io import BytesIO
                observed = next(Chem.ForwardSDMolSupplier(BytesIO(artifacts.get_bytes(ref["observed_sdf_sha256"])), removeHs=False), None)
                if (observed is None or not observed.GetNumConformers()
                        or chemical_identity(Chem.MolToSmiles(Chem.RemoveHs(observed)))["exact_structure_id"] != ref["identity"]["exact_structure_id"]):
                    raise ValueError("observed_reference_graph_or_coordinates_invalid")
                ensemble, ids = Chem.AddHs(observed, addCoords=True), [observed.GetConformer().GetId()]
            else:
                ensemble, ids = _ensemble(ref["smiles"], conformers, conformer_timeout_seconds)
            ensembles.append((ref, ensemble, ids))
        except ValueError as error:
            reference_conformer_holds.append({"candidate_id": ref["candidate_id"], "reason": str(error)})
    eligible = []
    for row in rows:
        if stop:
            stop()
        previous = {key: row.pop(key) for key in (
            "reference_similarity", "reference_alignment", "acyclic_topology_comparison",
        ) if key in row}
        if previous:
            previous["discovery_origin"] = deepcopy(row.get("discovery_origin"))
            history = row.get("source_reference_comparisons")
            row["source_reference_comparisons"] = [*(history if isinstance(history, list) else []), previous]
        row["discovery_origin"] = {
            "route": "generated_proposal" if row.get("derivation") else "general_catalog",
            "reference_policy_sha256": reference_set["sha256"],
            "activity_inherited": False,
        }
        identity = chemical_identity(row.get("smiles"))
        row["identity"] = identity
        row["novel_scaffold"] = (
            (
                identity["scaffold_kind"] == "murcko"
                and identity["scaffold_id"] not in cores
            )
            if reference_rows
            else None
        )
        topology = acyclic_topology(row.get("smiles"))
        if topology is not None:
            row["acyclic_topology_comparison"] = {
                **topology,
                "different_from_query_topologies": (
                    topology["id"] not in acyclic_cores if acyclic_cores else None
                ),
            }
        if (
            identity["structure_status"] != "valid"
            or identity["parent_ambiguous"]
            or not fingerprints
        ):
            continue
        candidate_fp = generator.GetFingerprint(Chem.MolFromSmiles(row["smiles"]))
        similarities = [
            (DataStructs.TanimotoSimilarity(candidate_fp, fp), ref["candidate_id"])
            for ref, fp in fingerprints
        ]
        similarity, ref_id = max(similarities)
        row["reference_similarity"] = {
            "tanimoto_2d": similarity,
            "reference_id": ref_id,
            "policy": "Morgan-radius2-2048",
            "binding_established": False,
        }
        if (
            row.get("discovery_filter_pass")
            and (row["novel_scaffold"] or
                 row.get("acyclic_topology_comparison", {}).get("different_from_query_topologies") is True)
            and row.get("role", "candidate") == "candidate"
        ):
            eligible.append(row)
    eligible.sort(
        key=lambda r: (
            -r["reference_similarity"]["tanimoto_2d"],
            -(finite(r.get("priority_score_v2")) or 0),
            r["candidate_id"],
        )
    )
    started = time.monotonic()
    completed = 0
    feature_cache = _ShapeFeatures()
    for index, row in enumerate(eligible[:shape_limit]):
        if stop:
            stop()
        elapsed = time.monotonic() - started
        # Small canary first, then conservative admission. Unfinished source
        # rows remain general-catalog fallback candidates under the same filter.
        if index and elapsed / index * 1.25 + elapsed > shape_seconds:
            break
        try:
            comparison = shape_compare(
                row["smiles"], ensembles, artifacts, conformers=conformers, stop=stop,
                _feature_cache=feature_cache,
                conformer_timeout_seconds=conformer_timeout_seconds,
            )
            row["reference_alignment"] = comparison
            row["discovery_origin"].update(
                route="reference_acyclic_topology_hop"
                if row.get("acyclic_topology_comparison") else "reference_catalog_hop",
                reference_id=comparison["reference_id"]
            )
            completed += 1
        except ValueError as error:
            row["reference_alignment"] = {"status": "held", "reason": str(error)}
        if progress:
            progress(
                {
                    "stage": "reference_shape",
                    "processed": index + 1,
                    "complete": completed,
                    "requested": min(shape_limit, len(eligible)),
                }
            )
    return rows, {
        "candidate_population": len(rows),
        "different_scaffold_retrieved": sum(r["novel_scaffold"] is True for r in eligible),
        "different_representation_retrieved": len(eligible),
        "acyclic_topology_hop_retrieved": sum(
            bool(r.get("acyclic_topology_comparison")) for r in eligible),
        "acyclic_topology_policy": "acyclic-linker-topology-v1",
        "conformer_timeout_seconds": conformer_timeout_seconds or None,
        "shape_completed": completed,
        "shape_limit": shape_limit,
        "shape_compute_seconds": time.monotonic() - started,
        "reference_policy_sha256": reference_set["sha256"],
        "reference_conformers_completed": len(ensembles),
        "reference_conformer_holds": reference_conformer_holds,
        "similarity_hard_floor": None,
        "novelty_scope": "Compared only with this frozen query/reference population; no patent, legal novelty or freedom-to-operate determination",
    }


def discovery_order(row):
    hop = row.get("discovery_origin", {}).get("route") in HOP_ROUTES
    return (
        not hop,
        -(finite(row.get("reference_similarity", {}).get("tanimoto_2d")) or 0)
        if row.get("research_search_mode") == "similarity" else 0,
        -(finite(row.get("reference_alignment", {}).get("retention_score")) or 0)
        if hop
        else 0,
        -(finite(row.get("priority_score_v2")) or 0),
        row["candidate_id"],
    )


def select_docking(rows, maximum=1000):
    """Hopping first, same-filter backfill, controls outside the quota."""
    used, selected = set(), []
    controls = {
        r.get("identity", {}).get("parent_id")
        for r in rows
        if r.get("role")
        in {"control", "positive_control", "negative_control", "counter_control", "query_seed"}
    }
    ordered = sorted(rows, key=discovery_order)
    hops = [
        r
        for r in ordered
        if r.get("discovery_origin", {}).get("route") in HOP_ROUTES
    ]
    catalog = [
        r
        for r in ordered
        if r.get("discovery_origin", {}).get("route") not in HOP_ROUTES
    ]
    cores = {}
    for r in catalog:
        core = r.get("identity", {}).get("scaffold_id")
        r["_catalog_diversity_round"] = cores.get(core, 0)
        cores[core] = cores.get(core, 0) + 1
    catalog.sort(key=lambda r: (r["_catalog_diversity_round"], discovery_order(r)))
    for r in catalog:
        r.pop("_catalog_diversity_round", None)
    for row in hops + catalog:
        identity = row["identity"]
        if (
            not row.get("discovery_filter_pass")
            or row.get("role", "candidate") != "candidate"
            or identity["structure_status"] != "valid"
            or identity["parent_ambiguous"]
            or row.get("identity_conflict")
            or row.get("engine_supported") is False
            or identity["parent_id"] in used | controls
        ):
            continue
        selected.append(row["candidate_id"])
        used.add(identity["parent_id"])
        if len(selected) == maximum:
            break
    return selected


def stored_alignment(row, artifacts):
    """Re-map common atoms directly onto the stored comparison SDF atom order."""
    from rdkit import Chem
    from rdkit.Chem import rdFMCS

    comparison = row.get("reference_alignment") or {}
    if comparison.get("status") != "complete":
        raise ValueError("Stored reference alignment unavailable")
    a = artifacts.get_bytes(comparison["reference_sdf_sha256"]).decode()
    b = artifacts.get_bytes(comparison["candidate_sdf_sha256"]).decode()
    ref = Chem.MolFromMolBlock(a.split("$$$$")[0], removeHs=False)
    candidate = Chem.MolFromMolBlock(b.split("$$$$")[0], removeHs=False)
    if ref is None or candidate is None:
        raise ValueError("Stored alignment SDF cannot be parsed")
    if (
        chemical_identity(Chem.MolToSmiles(Chem.RemoveHs(candidate)))[
            "exact_structure_id"
        ]
        != chemical_identity(row["smiles"])["exact_structure_id"]
    ):
        raise ValueError("Stored alignment candidate identity differs from frozen row")
    parameters = rdFMCS.MCSParameters()
    parameters.Timeout = 5
    parameters.AtomTyper = rdFMCS.AtomCompare.CompareElements
    parameters.BondTyper = rdFMCS.BondCompare.CompareOrderExact
    parameters.AtomCompareParameters.MatchChiralTag = True
    parameters.AtomCompareParameters.RingMatchesRingOnly = True
    parameters.BondCompareParameters.RingMatchesRingOnly = True
    parameters.BondCompareParameters.CompleteRingsOnly = True
    result = rdFMCS.FindMCS([Chem.RemoveHs(ref), Chem.RemoveHs(candidate)], parameters)
    query = Chem.MolFromSmarts(result.smartsString) if result.smartsString else None
    ref_atoms = ref.GetSubstructMatch(query, useChirality=True) if query else []
    candidate_atoms = (
        candidate.GetSubstructMatch(query, useChirality=True) if query else []
    )
    pairs = [
        {"reference_atom": i, "candidate_atom": j}
        for i, j in zip(ref_atoms, candidate_atoms)
    ]
    features = []
    for pair in comparison["feature_pairs"]:
        if any(i < 0 or i >= ref.GetNumAtoms() for i in pair["reference_atoms"]) or any(
            i < 0 or i >= candidate.GetNumAtoms() for i in pair["candidate_atoms"]
        ):
            raise ValueError("Stored feature atom index outside exact saved SDF")
        features.append(
            {
                "reference": {
                    "family": pair["family"],
                    "atom_indices": pair["reference_atoms"],
                },
                "candidate": {
                    "family": pair["family"],
                    "atom_indices": pair["candidate_atoms"],
                },
                "distance_angstrom": pair["distance_angstrom"],
            }
        )
    return {
        **comparison,
        "reference_sdf": a,
        "candidate_sdf": b,
        "mcs": {
            "atom_pairs": pairs,
            "partial_timeout": result.canceled,
            "complete_search": not result.canceled,
        },
        "feature_pairs": features,
        "feature_status": "evaluated",
        "limitations": [
            "Generated comparison alignment; no docking coordinates changed",
            "Reference activity is not transferred to the candidate",
        ],
        "atom_index_basis": "zero-based indices in these exact saved SDFs, verified on read",
    }
