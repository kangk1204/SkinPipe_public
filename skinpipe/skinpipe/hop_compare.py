"""Exact common-substructure comparison with explicit atom-index provenance."""
from __future__ import annotations
import hashlib
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, rdFMCS, rdFingerprintGenerator
from rdkit.Chem.Draw import rdMolDraw2D

BLUE = (31/255, 110/255, 196/255)
ORANGE = (240/255, 122/255, 20/255)
MAGENTA = (196/255, 46/255, 158/255)
POLICY = "hop-common-structure-v1"


def _svg(molecule, common, common_bonds, changed_color):
    mol = Chem.Mol(molecule)
    AllChem.Compute2DCoords(mol)
    drawer = rdMolDraw2D.MolDraw2DSVG(520, 340)
    common = set(common)
    atom_colors = {i: BLUE if i in common else changed_color for i in range(mol.GetNumAtoms())}
    bond_colors = {b.GetIdx(): BLUE if b.GetIdx() in common_bonds
                   else changed_color for b in mol.GetBonds()}
    drawer.DrawMolecule(mol, highlightAtoms=list(atom_colors), highlightBonds=list(bond_colors),
                        highlightAtomColors=atom_colors, highlightBondColors=bond_colors)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def compare(parent_smiles, candidate_smiles, timeout_s=5):
    parent, candidate = Chem.MolFromSmiles(parent_smiles), Chem.MolFromSmiles(candidate_smiles)
    if parent is None or candidate is None:
        raise ValueError("원본 또는 후보의 분자 구조를 읽을 수 없습니다.")
    if max(parent.GetNumHeavyAtoms(), candidate.GetNumHeavyAtoms()) > 250:
        raise ValueError("구조 비교는 분자당 250 heavy atoms 이하에서 지원합니다.")
    # Drawings and index maps refer to exactly these canonical structures.
    p_smi, c_smi = Chem.MolToSmiles(parent), Chem.MolToSmiles(candidate)
    parent, candidate = Chem.MolFromSmiles(p_smi), Chem.MolFromSmiles(c_smi)
    parameters = rdFMCS.MCSParameters()
    parameters.Timeout = timeout_s
    parameters.AtomTyper = rdFMCS.AtomCompare.CompareElements
    parameters.BondTyper = rdFMCS.BondCompare.CompareOrderExact
    parameters.AtomCompareParameters.MatchValences = True
    parameters.AtomCompareParameters.MatchChiralTag = True
    parameters.AtomCompareParameters.RingMatchesRingOnly = True
    parameters.BondCompareParameters.RingMatchesRingOnly = True
    parameters.BondCompareParameters.CompleteRingsOnly = True
    parameters.BondCompareParameters.MatchStereo = True
    result = rdFMCS.FindMCS([parent, candidate], parameters)
    query = Chem.MolFromSmarts(result.smartsString) if result.smartsString else None
    p_matches = sorted(parent.GetSubstructMatches(query, uniquify=False, useChirality=True,
                                                 maxMatches=1000)) if query else []
    c_matches = sorted(candidate.GetSubstructMatches(query, uniquify=False, useChirality=True,
                                                     maxMatches=1000)) if query else []
    p_atoms = list(p_matches[0]) if p_matches else []
    c_atoms = list(c_matches[0]) if c_matches else []
    def matched_bonds(mol, match):
        return [mol.GetBondBetweenAtoms(match[b.GetBeginAtomIdx()], match[b.GetEndAtomIdx()]).GetIdx()
                for b in query.GetBonds()] if query and match else []
    p_bonds, c_bonds = matched_bonds(parent, p_atoms), matched_bonds(candidate, c_atoms)
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    similarity = DataStructs.TanimotoSimilarity(generator.GetFingerprint(parent),
                                                generator.GetFingerprint(candidate))
    import rdkit
    return {"policy": POLICY, "rdkit_version": rdkit.__version__,
            "status": "PARTIAL_TIMEOUT" if result.canceled else "EVALUATED",
            "timeout_seconds": timeout_s, "parent_smiles": p_smi, "candidate_smiles": c_smi,
            "parent_sha256": hashlib.sha256(p_smi.encode()).hexdigest(),
            "candidate_sha256": hashlib.sha256(c_smi.encode()).hexdigest(),
            "atom_index_basis": "zero-based RDKit atom indices of the canonical SMILES shown here; not docking SDF indices",
            "mcs_smarts": result.smartsString, "common_atoms": result.numAtoms,
            "atom_pairs": [list(pair) for pair in zip(p_atoms, c_atoms)],
            "parent_changed_atoms": [i for i in range(parent.GetNumAtoms()) if i not in p_atoms],
            "candidate_changed_atoms": [i for i in range(candidate.GetNumAtoms()) if i not in c_atoms],
            "mapping_ambiguous": len(p_matches) > 1 or len(c_matches) > 1,
            "mapping_counts_capped_at": 1000,
            "mapping_choice": "lexicographically first chirality-compatible match in each structure",
            "matching_rules": "element, exact bond order, valence, complete rings, ring membership, chirality/stereo",
            "tanimoto_2d": round(similarity, 4),
            "parent_common_bonds": p_bonds, "candidate_common_bonds": c_bonds,
            "parent_svg": _svg(parent, p_atoms, p_bonds, ORANGE),
            "candidate_svg": _svg(candidate, c_atoms, c_bonds, MAGENTA),
            "legend": {"common": "#1f6ec4", "parent_changed": "#f07a14", "candidate_changed": "#c42e9e"},
            "interpretation": "exact common structure only; not pharmacophore equivalence or biological activity",
            "coordinates": "2D depiction; original docking coordinates are unchanged"}
