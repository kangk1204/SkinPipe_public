"""Bounded CMR screening, separating structural features from reference identity.

Reference SMILES were checked through public PubChem CAS queries on 2026-10-05.
They are a small screening reference set, never an exhaustive legal list.
Ignoring stereochemistry and neutralising terminal carboxylates deliberately
provides a conservative connectivity-family warning, not a statute lookup.
"""
from __future__ import annotations

from typing import Any

from rdkit import Chem

from ..types import ADMETProfile
from . import admet_alerts as al
from . import admet_photo as photo

VITAMIN_A_NOTE = (
    "EU Reg. 2024/996 specifies Retinol, Retinyl Acetate and Retinyl Palmitate: "
    "0.05 percent retinol equivalent in body lotion and 0.3 percent in other "
    "leave-on and rinse-off products. These are product concentrations in "
    "retinol equivalents, not limits for every beta-ionone/polyene compound. "
    "Identity, formulation and applicable legal conditions must be checked."
)

# PubChem CIDs 3026, 6782, 8343, 2347 and 444795, respectively. The cited
# Annex II entries motivate screening priority, not a definitive legal result.
_REFERENCES = (
    ("dibutyl_phthalate_reference", "84-74-2",
     "CCCCOC(=O)C1=CC=CC=C1C(=O)OCCCC", "675"),
    ("diisobutyl_phthalate_reference", "84-69-5",
     "CC(C)COC(=O)C1=CC=CC=C1C(=O)OCC(C)C", "1492"),
    ("bis_2_ethylhexyl_phthalate_reference", "117-81-7",
     "CCCCC(CC)COC(=O)C1=CC=CC=C1C(=O)OCC(CC)CCCC", "677"),
    ("benzyl_butyl_phthalate_reference", "85-68-7",
     "CCCCOC(=O)C1=CC=CC=C1C(=O)OCC2=CC=CC=C2", "1152"),
    ("retinoic_acid_reference", "302-79-4",
     "CC1=C(C(CCC1)(C)C)C=CC(=CC=CC(=CC(=O)O)C)C", "375"),
)


def _connectivity_key(mol: Chem.Mol) -> str:
    """Clone before normalising carboxylate salts for a reference-family check."""
    copy = Chem.RemoveHs(Chem.Mol(mol))
    for atom in copy.GetAtoms():
        if atom.GetAtomicNum() != 8 or atom.GetFormalCharge() != -1 or atom.GetDegree() != 1:
            continue
        carbon = atom.GetNeighbors()[0]
        if carbon.GetAtomicNum() != 6:
            continue
        if not any(neighbour.GetAtomicNum() == 8
                   and copy.GetBondBetweenAtoms(carbon.GetIdx(), neighbour.GetIdx()).GetBondType() == Chem.BondType.DOUBLE
                   for neighbour in carbon.GetNeighbors()):
            continue
        atom.SetFormalCharge(0)
        atom.SetNumExplicitHs(1)
        atom.SetNoImplicit(True)
    Chem.SanitizeMol(copy)
    return Chem.MolToSmiles(copy, canonical=True, isomericSmiles=False)


def _compile_references() -> tuple[tuple[str, str, str, str], ...]:
    out = []
    for name, cas, smiles, entry in _REFERENCES:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"unparseable CMR screening reference {name}")
        out.append((name, cas, _connectivity_key(mol), entry))
    return tuple(out)


_COMPILED_REFERENCES = _compile_references()


def assess(mol: Chem.Mol) -> dict[str, Any]:
    """Return screening warnings without claiming toxicity or legal clearance."""
    skip = set() if photo.longest_polyene(mol) >= 4 else {"retinoid"}
    alerts = al.match_alerts(mol, al.CMR_COMPILED, skip=skip)
    # Only candidate families can equal these references. Avoid expensive
    # canonicalisation on the many unrelated structures in a full population.
    if any(row["name"] in {"retinoid", "phthalate_diester"} for row in alerts):
        keys = {_connectivity_key(fragment)
                for fragment in Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)}
        for name, cas, key, entry in _COMPILED_REFERENCES:
            if key not in keys:
                continue
            alerts.append({
                "name": name, "severity": "critical", "n_matches": 1,
                "mechanism": f"component connectivity matches screening reference CAS {cas}; "
                             "stereochemistry, salt form and mixture composition need an identity review",
                "action": f"verify compound identity and EU Cosmetics Regulation Annex II entry {entry}; "
                          "screening connectivity is not a legal classification",
                "basis": "reference_component_connectivity",
                "reference_cas": cas,
                "reference_source": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
                                    + cas + "/property/ConnectivitySMILES,SMILES,IUPACName/JSON",
                "regulatory_reference_source": ("https://eur-lex.europa.eu/eli/reg/2019/831/oj"
                                                if entry == "1492" else
                                                "https://eur-lex.europa.eu/eli/reg/2009/1223/oj"),
            })
    levels = {row["severity"] for row in alerts}
    concern = next((level for level in ADMETProfile.SEVERITY_ORDER if level in levels), "none")
    verdict = (f"CMR screening concern {concern}: " + ", ".join(row["name"] for row in alerts)
               + "; compound-specific toxicity, legal identity and use-level exposure remain to be assessed"
               if alerts else
               "No CMR alert matched the evaluated rules; toxicity and legal classification remain unestablished")
    result = {
        "concern": concern, "alerts": alerts, "verdict": verdict,
        "basis": "structural-feature triage and a limited reference-connectivity screen",
        "legal_classification": "NOT_ESTABLISHED", "experimental_toxicity": "NOT_ESTABLISHED",
        "reference_coverage": {"compounds": len(_REFERENCES), "exhaustive": False,
                               "stereochemistry_checked": False},
        "rules_version": "cmr-specificity-20261005",
    }
    if any(row["name"] == "retinoid" for row in alerts):
        result["vitamin_a_context"] = VITAMIN_A_NOTE
        result["vitamin_a_source"] = "https://eur-lex.europa.eu/eli/reg/2024/996/oj"
    return result
