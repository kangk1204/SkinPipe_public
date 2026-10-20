"""Versioned, conservative molecular identity. No salt/stereo/isotope collapsing."""

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors
from .common import digest

RDLogger.DisableLog("rdApp.*")
NORMALIZER_VERSION = "rdkit-2025.09.5-exact-isomeric-v1"
SUPPORTED = {1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 35, 53}


def normalize(smiles):
    empty = dict(
        compound_id=None,
        canonical_smiles=None,
        identity_smiles=None,
        inchikey=None,
        parent_smiles=None,
        parent_id=None,
        formula=None,
        mw=None,
        fragments=None,
        formal_charge=None,
        heavy_atoms=None,
        structure_status="missing_structure",
        inference_status="unsupported_missing_structure",
        normalizer=NORMALIZER_VERSION,
    )
    if not smiles or not str(smiles).strip():
        return empty
    try:
        m = Chem.MolFromSmiles(str(smiles).strip())
        if m is None or not m.GetNumAtoms():
            return {
                **empty,
                "structure_status": "invalid_smiles",
                "inference_status": "unsupported_invalid_smiles",
            }
        canonical = Chem.MolToSmiles(m, canonical=True, isomericSmiles=True)
        params = Chem.SmilesWriteParams()
        params.canonical = True
        params.doIsomericSmiles = True
        identity = Chem.MolToCXSmiles(m, params, Chem.CXSmilesFields.CX_ENHANCEDSTEREO)
        fragments = Chem.GetMolFrags(m, asMols=True)
        # A largest-fragment relationship is descriptive only; it never replaces identity.
        ordered = sorted(
            fragments, key=lambda x: (-x.GetNumHeavyAtoms(), Chem.MolToSmiles(x))
        )
        parent = Chem.MolToSmiles(ordered[0], canonical=True, isomericSmiles=True)
        elements = {a.GetAtomicNum() for a in m.GetAtoms()}
        status = "eligible"
        if 0 in elements:
            status = "unsupported_query_or_polymer"
        elif len(fragments) > 1:
            status = "unsupported_multi_component"
        elif not elements <= SUPPORTED:
            status = "unsupported_element"
        elif m.GetNumHeavyAtoms() > 200:
            status = "outside_operational_size_limit"
        elif m.GetStereoGroups():
            status = "unsupported_enhanced_stereo"
        try:
            ik = Chem.MolToInchiKey(m)
        except Exception:
            ik = None
        return dict(
            compound_id=digest(NORMALIZER_VERSION + "|" + identity),
            canonical_smiles=canonical,
            identity_smiles=identity,
            inchikey=ik,
            parent_smiles=parent,
            parent_id=digest(NORMALIZER_VERSION + "|" + parent),
            formula=rdMolDescriptors.CalcMolFormula(m),
            mw=Descriptors.MolWt(m),
            fragments=len(fragments),
            formal_charge=Chem.GetFormalCharge(m),
            heavy_atoms=m.GetNumHeavyAtoms(),
            structure_status="valid",
            inference_status=status,
            normalizer=NORMALIZER_VERSION,
        )
    except Exception:
        return {
            **empty,
            "structure_status": "normalization_error",
            "inference_status": "unsupported_normalization_error",
        }


def normalize_record(row):
    return {**row, **normalize(row.get("raw_smiles"))}
