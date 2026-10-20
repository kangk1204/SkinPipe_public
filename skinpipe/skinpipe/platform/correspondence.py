"""Inspectable 2D common-substructure and 3D feature correspondence.

Generated conformer alignment illustrates structural hypotheses; it is not a
docking pose or evidence of binding. Original heavy-atom indices are preserved.
"""
from __future__ import annotations

import math
from pathlib import Path
from typing import Any

from rdkit import Chem, RDConfig
from rdkit.Chem import AllChem, ChemicalFeatures, rdFMCS, rdMolAlign


def _molecule(value):
    mol = Chem.Mol(value) if isinstance(value, Chem.Mol) else Chem.MolFromSmiles(value)
    if mol is None or not mol.GetNumAtoms():
        raise ValueError("valid molecular structures are required")
    if len(Chem.GetMolFrags(mol)) != 1:
        raise ValueError("correspondence requires a specified single-component structure")
    return mol


def _conformer(mol, seed):
    if mol.GetNumConformers() and mol.GetConformer().Is3D():
        return Chem.Mol(mol), "supplied_coordinates"
    hydrogenated = Chem.AddHs(mol)
    hydrogenated.RemoveAllConformers()
    params = AllChem.ETKDGv3()
    params.randomSeed = int(seed)
    if AllChem.EmbedMolecule(hydrogenated, params) != 0:
        return None, "embedding_failed"
    optimization = "not_available"
    if AllChem.MMFFHasAllMoleculeParams(hydrogenated):
        optimization = "converged" if AllChem.MMFFOptimizeMolecule(hydrogenated, maxIters=300) == 0 else "iteration_limit"
    return hydrogenated, "generated_ETKDG3_MMFF:" + optimization


def _atoms(mol, count):
    conformer = mol.GetConformer()
    return [{"atom_index": i, "element": mol.GetAtomWithIdx(i).GetSymbol(),
             "position": list(conformer.GetAtomPosition(i))} for i in range(count)]


def compare_structures(reference_smiles, candidate_smiles, *, seed: int = 42,
                       feature_tolerance: float = 2.0, mcs_timeout: int = 5) -> dict[str, Any]:
    """Return atom-level MCS and greedy one-to-one same-family 3D feature pairs.

    RDKit Mol inputs preserve supplied 3D coordinates. SMILES inputs generate
    conformers deterministically. 3D features are returned only after >=3 common
    atoms provide an alignment; symmetry may allow other equally valid mappings.
    """
    if isinstance(seed, bool) or not isinstance(seed, int) or not 0 <= seed < 2**31:
        raise ValueError("seed must be an integer in [0, 2**31)")
    if not isinstance(feature_tolerance, (float, int)) or isinstance(feature_tolerance, bool) or not math.isfinite(feature_tolerance) or feature_tolerance <= 0:
        raise ValueError("feature tolerance must be finite and positive")
    if isinstance(mcs_timeout, bool) or not isinstance(mcs_timeout, int) or mcs_timeout < 1:
        raise ValueError("MCS timeout must be a positive integer")
    ref, cand = _molecule(reference_smiles), _molecule(candidate_smiles)
    mcs = rdFMCS.FindMCS([ref, cand], timeout=mcs_timeout, ringMatchesRingOnly=True,
                        completeRingsOnly=True, matchValences=True, matchChiralTag=True,
                        bondCompare=rdFMCS.BondCompare.CompareOrderExact)
    query = Chem.MolFromSmarts(mcs.smartsString)
    ref_match = ref.GetSubstructMatch(query, useChirality=True) if query is not None else ()
    cand_match = cand.GetSubstructMatch(query, useChirality=True) if query is not None else ()
    out: dict[str, Any] = {"method": "MCS atom mapping + ETKDG3/MMFF aligned feature illustration v1",
        "atom_index_base": 0, "coordinate_units": "angstrom", "seed": seed,
        "mcs": {"smarts": mcs.smartsString, "complete_search": not mcs.canceled,
            "atom_pairs": [{"reference_atom": a, "candidate_atom": b}
                           for a, b in zip(ref_match, cand_match)],
            "reference_fraction": len(ref_match) / ref.GetNumAtoms(),
            "candidate_fraction": len(cand_match) / cand.GetNumAtoms()},
        "feature_pairs": [], "unmatched_reference_features": [],
        "unmatched_candidate_features": [], "feature_tolerance_angstrom": feature_tolerance,
        "feature_matching": "distance-sorted greedy one-to-one same family",
        "feature_status": "not_evaluated", "novelty_established": False,
        "limitations": ["illustrative alignment is not docking or evidence of binding",
                        "symmetry can yield alternative common-substructure mappings"]}
    ref3d, ref_source = _conformer(ref, seed)
    cand3d, cand_source = _conformer(cand, seed)
    out["coordinate_sources"] = {"reference": ref_source, "candidate": cand_source}
    if ref3d is None or cand3d is None:
        out["feature_status"] = "embedding_failed"
        return out
    if len(ref_match) < 3:
        out["feature_status"] = "insufficient_alignment_atoms"
        return out
    rmsd = rdMolAlign.AlignMol(cand3d, ref3d, atomMap=list(zip(cand_match, ref_match)))
    out["common_atom_alignment_rmsd"] = rmsd
    out["reference_sdf"] = Chem.MolToMolBlock(ref3d) + "\n$$$$\n"
    out["candidate_sdf"] = Chem.MolToMolBlock(cand3d) + "\n$$$$\n"
    out["sdf_atom_index_mapping"] = {
        "index_base": 0, "mapping": "input atom index equals SDF atom record index; added H follow input atoms",
        "reference_input_atom_count": ref.GetNumAtoms(), "candidate_input_atom_count": cand.GetNumAtoms()}
    out["reference_atoms"] = _atoms(ref3d, ref.GetNumAtoms())
    out["candidate_atoms"] = _atoms(cand3d, cand.GetNumAtoms())
    factory = ChemicalFeatures.BuildFeatureFactory(str(Path(RDConfig.RDDataDir) / "BaseFeatures.fdef"))
    def features(mol, original_count):
        return [{"feature_index": i, "family": feat.GetFamily(),
                 "atom_indices": [a for a in feat.GetAtomIds() if a < original_count],
                 "position": list(feat.GetPos())}
                for i, feat in enumerate(factory.GetFeaturesForMol(mol))]
    reference_features = features(ref3d, ref.GetNumAtoms())
    candidate_features = features(cand3d, cand.GetNumAtoms())
    possible = []
    for a in reference_features:
        for b in candidate_features:
            if a["family"] != b["family"]:
                continue
            distance = math.dist(a["position"], b["position"])
            if distance <= feature_tolerance:
                possible.append((distance, a["feature_index"], b["feature_index"]))
    used_ref, used_cand = set(), set()
    for distance, a, b in sorted(possible):
        if a in used_ref or b in used_cand:
            continue
        used_ref.add(a)
        used_cand.add(b)
        out["feature_pairs"].append({"reference": reference_features[a],
                                     "candidate": candidate_features[b],
                                     "distance_angstrom": distance})
    out["unmatched_reference_features"] = [f for f in reference_features if f["feature_index"] not in used_ref]
    out["unmatched_candidate_features"] = [f for f in candidate_features if f["feature_index"] not in used_cand]
    out["feature_status"] = "illustrative_correspondence"
    return out
