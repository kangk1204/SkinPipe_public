"""Exact structure identity and explicitly separate parent/preparation relations.

Exact IDs preserve salt, charge, isotope and stereochemistry, including enhanced
stereo. Parent IDs are *only* a selection grouping; never use them as cache keys.
"""
from __future__ import annotations

import hashlib
import json
from copy import deepcopy
from functools import lru_cache
from typing import Any

from rdkit import Chem, rdBase
from rdkit.Chem import Descriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold

IDENTITY_VERSION = "skinpipe-exact-cx-v1"
PARENT_VERSION = "largest-organic-neutral-parent-v1"
STEREO_POLICY = "potential-stereo-heavy-sp2-v2"


def _reversible_protonated_nitrogen(mol, index):
    atom = mol.GetAtomWithIdx(index)
    if atom.GetAtomicNum() != 7 or atom.GetFormalCharge() != 1 or atom.GetTotalNumHs() != 1:
        return False
    # A temporary deprotonated graph is used only to distinguish an inverting
    # tertiary amine from a persistent stereocentre. Docking retains the actual
    # charged graph. Quaternary ammonium (no N-H) is never exempted.
    neutral = Chem.Mol(mol)
    nitrogen = neutral.GetAtomWithIdx(index)
    nitrogen.SetFormalCharge(0)
    nitrogen.SetNumExplicitHs(0)
    nitrogen.SetNoImplicit(True)
    nitrogen.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    try:
        Chem.SanitizeMol(neutral)
        return not any(info.type == Chem.StereoType.Atom_Tetrahedral and info.centeredOn == index
                       for info in Chem.FindPotentialStereo(neutral))
    except (ValueError, RuntimeError):
        return False


def stereochemistry_assessment(mol: Chem.Mol) -> dict[str, Any]:
    """Shared source/effective graph check; never assign an absent isomer.

    RDKit potential stereo recognises actual E/Z-capable bonds; terminal or
    symmetric double bonds are not unspecified stereochemistry. Atom and bond
    indexes are separate namespaces and remain separate exported fields.
    """
    atoms, bonds, reversible = [], [], []
    try:
        potential = Chem.FindPotentialStereo(Chem.Mol(mol))
    except (ValueError, RuntimeError) as error:
        raise ValueError("stereochemistry assessment failed; preparation held") from error
    for info in potential:
        if info.specified == Chem.StereoSpecified.Specified:
            continue
        if info.type == Chem.StereoType.Atom_Tetrahedral:
            if _reversible_protonated_nitrogen(mol, info.centeredOn):
                reversible.append(int(info.centeredOn))
            else:
                atoms.append(int(info.centeredOn))
        elif info.type == Chem.StereoType.Bond_Double:
            bond = mol.GetBondWithIdx(int(info.centeredOn))
            ends = (bond.GetBeginAtom(), bond.GetEndAtom())
            # FindPotentialStereo can also report terminal imines/cumulenes
            # that cannot carry an explicit E/Z descriptor in SMILES. Apply
            # the documented heavy-neighbour criterion before holding input.
            # Keep substituted imines/azo bonds and actual alkene E/Z gates.
            if all(sum(n.GetAtomicNum() > 1 or n.GetIsotope() > 0 for n in atom.GetNeighbors()) >= 2
                   and atom.GetHybridization() == Chem.HybridizationType.SP2
                   for atom in ends):
                bonds.append(int(info.centeredOn))
        else:
            # Unsupported potential stereo is not silently assigned.
            atoms.append(int(info.centeredOn))
    # The new RDKit detector omits a C=C end bearing D/H even though
    # isomeric SMILES can represent it. Consult the isotope-aware legacy bond
    # detector only for that narrow case, on a copy; retain supplied E/Z.
    isotope_ends = {a.GetIdx() for a in mol.GetAtoms()
                    if any(n.GetAtomicNum() == 1 and n.GetIsotope() > 1 for n in a.GetNeighbors())}
    if isotope_ends:
        probe = Chem.Mol(mol)
        Chem.FindPotentialStereoBonds(probe, cleanIt=False)
        for bond in probe.GetBonds():
            original = mol.GetBondWithIdx(bond.GetIdx())
            ends = (bond.GetBeginAtom(), bond.GetEndAtom())
            if (bond.GetBondType() == Chem.BondType.DOUBLE
                    and bond.GetStereo() == Chem.BondStereo.STEREOANY
                    and original.GetStereo() in (Chem.BondStereo.STEREONONE, Chem.BondStereo.STEREOANY)
                    and all(a.GetAtomicNum() == 6 and a.GetHybridization() == Chem.HybridizationType.SP2 for a in ends)
                    and any(a.GetIdx() in isotope_ends for a in ends)):
                bonds.append(bond.GetIdx())
    return {"stereo_policy": STEREO_POLICY, "unassigned_stereo": sorted(atoms),
            "unassigned_bond_stereo": sorted(set(bonds)),
            "reversible_protonated_nitrogen": sorted(reversible)}


def digest(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                    ensure_ascii=True, allow_nan=False).encode()).hexdigest()


def _canonical(mol: Chem.Mol) -> str:
    clean = Chem.Mol(mol)
    for atom in clean.GetAtoms():
        atom.SetAtomMapNum(0)  # Input atom labels are provenance, not chemistry.
    params = Chem.SmilesWriteParams()
    params.canonical = True
    params.doIsomericSmiles = True
    return Chem.MolToCXSmiles(clean, params, Chem.CXSmilesFields.CX_ENHANCEDSTEREO)


def scaffold_identity(smiles: str) -> dict[str, Any]:
    """Murcko cores disregard stereo; acyclic compounds are not novel ring cores."""
    mol = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
    if mol is None or not mol.GetNumAtoms():
        return {"scaffold_id": None, "scaffold_smiles": None, "scaffold_kind": "invalid"}
    core = MurckoScaffold.GetScaffoldForMol(mol)
    if not core.GetNumAtoms():
        return {"scaffold_id": "acyclic", "scaffold_smiles": "", "scaffold_kind": "acyclic"}
    for atom in core.GetAtoms():
        atom.SetIsotope(0)
        atom.SetAtomMapNum(0)
    smi = Chem.MolToSmiles(core, canonical=True, isomericSmiles=False)
    return {"scaffold_id": "murcko:" + digest(smi), "scaffold_smiles": smi,
            "scaffold_kind": "murcko"}


def deferred_identity(smiles, *, reason="identity_timeout"):
    """Keep a traceable identity calculation hold, distinct from invalid chemistry."""
    result = _base_identity(None)
    result.update(original_smiles=smiles, original_structure_id="raw:" + digest(smiles),
                  structure_status="identity_deferred", deferred_reason=reason)
    return result


def candidate_identity(candidate, *, trust_stored=False):
    'Reuse only an explicitly trusted, source/version-bound server derivation.'
    smiles = candidate.get("smiles")
    stored = candidate.get("identity")
    if (trust_stored and isinstance(stored, dict)
            and stored.get("identity_version") == IDENTITY_VERSION
            and stored.get("parent_policy") == PARENT_VERSION
            and stored.get("normalizer_runtime") == rdBase.rdkitVersion
            and stored.get("stereo_policy") == STEREO_POLICY
            and stored.get("original_smiles") == smiles
            and stored.get("original_structure_id") == "raw:" + digest(smiles)
            and stored.get("structure_status") in {"valid", "invalid_structure", "missing_structure", "identity_deferred"}):
        if stored["structure_status"] != "valid" or (
                stored.get("exact_structure_id") == "exact:" + digest([IDENTITY_VERSION, stored.get("identity_smiles")])
                and stored.get("parent_id") == "parent:" + digest([PARENT_VERSION, stored.get("parent_smiles")])
                and all(key in stored for key in ("scaffold_id", "scaffold_kind", "parent_ambiguous"))):
            return deepcopy(stored)
    return chemical_identity(smiles)


def chemical_identity(smiles: str | None, *, prepared_smiles: str | None = None,
                      preparation: dict[str, Any] | None = None) -> dict[str, Any]:
    """Describe exact material, selection parent, and a versioned prepared state.

    A multi-organic-component material is flagged as ambiguous; its largest
    fragment does not establish the active ingredient. Invalid input remains a
    traceable record with null chemical IDs. Prepared IDs require a nonempty
    preparation context (engine/version/settings, supplied by the caller).
    """
    # Repeated compounds occur across many marker/target pools. Cache only the
    # exact input and versioned structural derivation; never share mutable rows.
    if isinstance(smiles, (str, type(None))):
        out = deepcopy(_cached_identity(smiles, IDENTITY_VERSION, PARENT_VERSION, rdBase.rdkitVersion, STEREO_POLICY))
    else:
        out = _base_identity(smiles)
    out["preparation"] = deepcopy(preparation)
    if out["structure_status"] == "valid" and prepared_smiles is not None:
        if not isinstance(preparation, dict) or not preparation:
            raise ValueError("prepared structure requires nonempty preparation context")
        prepared = Chem.MolFromSmiles(prepared_smiles)
        if prepared is None or not prepared.GetNumAtoms():
            raise ValueError("invalid prepared structure")
        prepared_exact = _canonical(prepared)
        out.update(prepared_identity_smiles=prepared_exact,
                   prepared_structure_id="prepared:" + digest(
                       [IDENTITY_VERSION, out["identity_smiles"], prepared_exact, preparation]))
    return out


@lru_cache(maxsize=8192)
def _cached_identity(smiles, identity_version, parent_version, runtime, stereo_policy):
    # Explicit version/runtime arguments invalidate memoized results if policy
    # changes within a long-lived API process. No prepared state enters this cache.
    return _base_identity(smiles)


def _base_identity(smiles):
    out: dict[str, Any] = {"identity_version": IDENTITY_VERSION,
        "stereo_policy": STEREO_POLICY,
        "normalizer_runtime": rdBase.rdkitVersion, "original_smiles": smiles,
        "original_structure_id": "raw:" + digest(smiles), "exact_structure_id": None,
        "canonical_smiles": None, "identity_smiles": None, "parent_id": None,
        "parent_smiles": None, "parent_policy": PARENT_VERSION,
        "prepared_structure_id": None, "prepared_identity_smiles": None,
        "preparation": None, "structure_status": "missing_structure",
        "parent_ambiguous": False, "scaffold_id": None, "scaffold_kind": "invalid"}
    if not isinstance(smiles, str) or not smiles.strip():
        return out
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or not mol.GetNumAtoms():
        out["structure_status"] = "invalid_structure"
        return out
    exact = _canonical(mol)
    fragments = list(Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True))
    organic = [m for m in fragments if any(a.GetAtomicNum() == 6 for a in m.GetAtoms())]
    pool = organic or fragments
    pool.sort(key=lambda m: (-m.GetNumHeavyAtoms(), _canonical(m)))
    parent = rdMolStandardize.Uncharger().uncharge(Chem.Mol(pool[0]))
    parent_smi = _canonical(parent)
    out.update(exact_structure_id="exact:" + digest([IDENTITY_VERSION, exact]),
               canonical_smiles=Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True),
               identity_smiles=exact, parent_id="parent:" + digest([PARENT_VERSION, parent_smi]),
               parent_smiles=parent_smi, parent_ambiguous=len(organic) > 1,
               fragment_count=len(fragments),
               component_identities=sorted(_canonical(m) for m in fragments),
               formal_charge=Chem.GetFormalCharge(mol),
               isotope_atoms=[{"atom_index": a.GetIdx(), "isotope": a.GetIsotope()}
                              for a in mol.GetAtoms() if a.GetIsotope()],
               **stereochemistry_assessment(mol),
               enhanced_stereo=bool(mol.GetStereoGroups()),
               heavy_atoms=mol.GetNumHeavyAtoms(), molecular_weight=Descriptors.MolWt(mol),
               structure_status="valid", **scaffold_identity(parent_smi))
    return out
