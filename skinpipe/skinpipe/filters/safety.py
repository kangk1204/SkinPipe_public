'Skin sensitisation alerts for topical candidates.'
from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Optional

from rdkit import Chem, RDLogger
from rdkit.Chem import FilterCatalog

RDLogger.DisableLog("rdApp.*")
log = logging.getLogger(__name__)


@dataclass(frozen=True)
class Alert:
    name: str
    domain: str          # mechanistic reaction domain
    severity: str        # high | moderate | low | context
    smarts: str
    note: str = ""


# Electrophilic mechanistic domains. SMARTS are deliberately conservative:
# a false negative is recoverable in the assay, a false positive that kills a
# whole chemotype before anyone looks at it is not.
_ALERTS: tuple[Alert, ...] = (
    # --- Michael acceptors -------------------------------------------------
    Alert("enone", "Michael acceptor", "high",
          "[$([CX3]=[CX3][CX3](=[OX1])[#6]),$([CX3]=[CX3][CX3H1]=[OX1])]",
          "alpha,beta-unsaturated carbonyl; adds thiol directly"),
    Alert("acrylate_ester", "Michael acceptor", "high",
          "[CH2]=[CH][CX3](=[OX1])[OX2][#6]",
          "acrylates are among the most common occupational skin allergens"),
    Alert("acrylamide", "Michael acceptor", "high",
          "[CH2]=[CH][CX3](=[OX1])[NX3]", ""),
    Alert("beta_substituted_enoate_ester", "Michael acceptor", "context",
          "[$([#6][CX3;H0,H1]=[CX3][CX3](=[OX1])[OX2][#6]),"
          "$([CX3;R]1=[CX3;R][CX3;R](=[OX1])[OX2;R][CX4;R]1)]",
          "beta-substituted conjugated ester, including unsaturated lactones; "
          "reactivity depends on substitution and use conditions. Structural "
          "screening alert, not a sensitisation diagnosis or use prohibition"),
    # The alpha carbon is spelled [CX3]([CX4]), not [CH], so the 2-substituted
    # analogues match. Written [CH] the two alerts above see ethyl acrylate and
    # are blind to ethyl methacrylate, 2-HEMA, HPMA and EGDMA - half the class,
    # and the half the clinic sees most. Requiring a terminal [CH2]= and an
    # sp3 carbon on the alpha position keeps this off unrelated trisubstituted
    # alkenes; the alpha-methylene lactones it does also catch (parthenolide,
    # helenalin) are genuine Michael acceptors and genuine sensitisers.
    Alert("alpha_substituted_acrylate_ester", "Michael acceptor", "high",
          "[CH2]=[CX3]([CX4])[CX3](=[OX1])[OX2][#6]",
          "Alpha-substituted acrylate ester, including alpha-methylene lactones; "
          "the matched group is not necessarily a methacrylate. Graded the same as the parent "
          "acrylate, not below it, although the alpha-methyl makes it the "
          "weaker electrophile - electron-donating and sterically hindering, "
          "so conjugate addition is roughly an order of magnitude slower. "
          "A severity here is a statement about sensitisation risk in use, "
          "and on that evidence the methacrylates are if anything the worse "
          "class: 2-HEMA is among the most frequently positive allergens in "
          "current patch-test series and the nail-product epidemic of "
          "allergic contact dermatitis led the EU to restrict HEMA and "
          "Di-HEMA Trimethylhexyl Dicarbamate to professional nail-product "
          "use (Annex III entries 313 and 314, Reg. (EU) 2020/1682). "
          "These restrictions are substance-specific, not a restriction of "
          "every methacrylate. Grading it below the assay threshold would be wrong "
          "advice about the best-documented member of the class"),
    Alert("methacrylamide", "Michael acceptor", "high",
          "[CH2]=[CX3]([CX4])[CX3](=[OX1])[NX3]",
          "2-substituted acrylamide; see alpha_substituted_acrylate_ester for the grade"),
    Alert("ynoate", "Michael acceptor", "high",
          "[CX2]#[CX2][CX3](=[OX1])[OX2,NX3]",
          "alkyne conjugated to an ester or amide. An ynoate is a Michael "
          "acceptor as much as an enoate is, and adds thiol to give a vinyl "
          "sulfide; methyl 2-octynoate (methyl heptine carbonate) is an EU "
          "Annex III labelled fragrance allergen and IFRA-restricted"),
    Alert("fumarate_maleate_diester", "Michael acceptor", "high",
          "[#6][OX2][CX3](=[OX1])[CX3]=[CX3][CX3](=[OX1])[OX2][#6]",
          "dimethyl fumarate class; potent NRF2 activator and potent sensitiser"),
    Alert("para_quinone", "Michael acceptor", "high",
          "[#6]1(=[OX1])[#6]=,:[#6][#6](=[OX1])[#6]=,:[#6]1",
          "para-quinone pattern; ortho-quinones are covered by the enone and "
          "alpha-dicarbonyl patterns rather than this pattern. Quinones can "
          "haptenate; also the activated form of many prohaptens"),
    Alert("maleimide", "Michael acceptor", "high",
          "[#7;r5]1-[#6;r5](=[OX1])-[#6;r5]=[#6;r5]-[#6;r5]1=[OX1]", ""),
    Alert("vinyl_sulfone", "Michael acceptor", "high",
          "[CH2]=[CH][SX4](=[OX1])(=[OX1])", ""),
    Alert("unsaturated_nitrile", "Michael acceptor", "moderate",
          "[CX3]=[CX3][CX2]#[NX1]", ""),
    Alert("isothiocyanate", "Michael acceptor", "context",
          "[NX2]=[CX2]=[SX1]",
          "sulforaphane class; reactive by design, this is the NRF2 mechanism"),

    # --- Schiff base formers ----------------------------------------------
    Alert("aldehyde", "Schiff base", "moderate",
          "[$([CX3H1](=[OX1])[#6]),$([CX3H2]=[OX1])]",
          "condenses with lysine; cinnamal and many fragrance allergens sit "
          "here. The second alternative is formaldehyde itself, which the "
          "carbon-substituent requirement of the first one excluded - the "
          "archetypal sensitiser, reported at the grade of its class"),
    Alert("alpha_dicarbonyl", "Schiff base", "high",
          "[CX3](=[OX1])[CX3](=[OX1])", ""),

    # --- Acyl transfer agents ---------------------------------------------
    Alert("anhydride", "Acyl transfer", "high",
          "[CX3](=[OX1])[OX2][CX3](=[OX1])", ""),
    Alert("acyl_halide", "Acyl transfer", "high",
          "[CX3](=[OX1])[F,Cl,Br,I]", ""),
    Alert("isocyanate", "Acyl transfer", "high",
          "[NX2]=[CX2]=[OX1]", ""),
    Alert("beta_lactam", "Acyl transfer", "high",
          "C1C(=O)NC1", ""),

    # --- SN2 electrophiles -------------------------------------------------
    Alert("epoxide", "SN2", "high",
          "[OX2r3]1[#6r3][#6r3]1",
          "ring strain drives alkylation of cysteine and lysine"),
    Alert("aziridine", "SN2", "high", "[NX3r3]1[#6r3][#6r3]1", ""),
    Alert("alkyl_halide", "SN2", "moderate",
          "[CX4;H2,H3][Cl,Br,I]", ""),
    # A sultone is by definition the *cyclic* ester of a hydroxy-sulfonic acid,
    # and ring strain is what makes it alkylate. Without `;R` this matched every
    # open-chain alkyl sulfate: sodium lauryl sulfate -- the standard LLNA
    # *negative* control and the most common anionic surfactant in cosmetics --
    # came back "protein-reactive; needs a sensitisation assay before use", and
    # so did SLES. A skin-actives screen that rejects its own commonest
    # surfactant class rejects the chemistry it exists to find.
    Alert("sultone", "SN2", "high",
          "[SX4;R](=[OX1])(=[OX1])[OX2;R][CX4;R]",
          "cyclic sulfonate ester; ring strain drives the alkylation"),
    # The genuine open-chain alkylators the constraint above would otherwise
    # lose. Both oxygens must be esterified -- dimethyl sulfate matches, a
    # lauryl sulfate *salt* does not, because its second oxygen is an anion and
    # the charge kills the SN2.
    # `!R` for the same reason it is on `alkyl_sulfonate_ester` below: a cyclic
    # sulfate (ethylene sulfate) is the ring form of this connectivity and is
    # already counted by `sultone` above, so without it that one compound
    # scored burden 6.0 where every other alert row scores 3.0.
    Alert("dialkyl_sulfate", "SN2", "high",
          "[CX4][OX2][SX4;!R](=[OX1])(=[OX1])[OX2][CX4]",
          "dimethyl/diethyl sulfate class; a salt of a monoester is not this"),
    # The other open-chain alkylator the ring constraint on `sultone` would
    # have dropped: a true sulfonate ester, S bonded straight to carbon.
    # Ethyl methanesulfonate is the textbook SN2 mutagen and it used to be
    # caught here under the wrong name. `!R` keeps it off the sultones above,
    # which are the same connectivity closed into a ring and would otherwise
    # be counted twice in `burden`.
    Alert("alkyl_sulfonate_ester", "SN2", "high",
          "[#6][SX4;!R](=[OX1])(=[OX1])[OX2][CX4]",
          "ethyl methanesulfonate class; direct SN2 alkylation of protein "
          "nucleophiles"),
    Alert("isothiazolinone", "SN2", "high",
          "[#16;r5]1~[#7;r5]~[#6;r5](=[OX1])~[#6;r5]~[#6;r5]~1",
          "MIT/CMIT preservatives; a major cause of cosmetic contact allergy"),

    # --- SNAr --------------------------------------------------------------
    Alert("activated_aryl_halide", "SNAr", "high",
          "[F,Cl,Br,I][c]1[c][c][c]([N+](=O)[O-])[c][c]1",
          "dinitrochlorobenzene class; the experimental positive control"),

    # --- Prohaptens: need metabolic or air oxidation first ------------------
    Alert("hydroquinone", "Prohapten (oxidation)", "high",
          "[OX2H][c]1[c][c][c]([OX2H])[c][c]1",
          "air-oxidises to the quinone; hydroquinone and arbutin's aglycone"),
    Alert("catechol", "Prohapten (oxidation)", "high",
          "[OX2H][c]1[c]([OX2H])[c][c][c][c]1",
          "urushiol class; oxidises to the ortho-quinone"),
    # Every aromatic nitrogen below is spelled with the same exclusion the
    # genotoxicity table uses for the nitrenium alerts: an acylated or
    # sulfonylated aniline has its lone pair delocalised into its own
    # carbonyl and cannot donate into the ring to make a quinoid. Without it,
    # imatinib - an anilide meta to an anilino-pyrimidine - reads as
    # m-phenylenediamine, and every paracetamol-type anilide reads as an
    # aminophenol.
    # Ortho and para are separate patterns, as they already are for the
    # diphenols above. Only the 1,2- and 1,4-isomers can be oxidised to a
    # quinone-imine or a diimine; written para-only, these alerts saw
    # 4-aminophenol and were blind to 2-aminophenol, which is on Korea's
    # 별표 1 (CAS 95-55-6), and to o-phenylenediamine.
    Alert("para_aminophenol", "Prohapten (oxidation)", "high",
          "[NX3;H2,H1;!$(N[CX3]=[OX1,SX1,NX2]);!$(N[SX4])][c]1[c][c][c]([OX2H])[c][c]1",
          "oxidises to the para-quinone-imine, which is the hapten"),
    Alert("ortho_aminophenol", "Prohapten (oxidation)", "high",
          "[NX3;H2,H1;!$(N[CX3]=[OX1,SX1,NX2]);!$(N[SX4])][c]1[c]([OX2H])[c][c][c][c]1",
          "oxidises to the ortho-quinone-imine; o-aminophenol is prohibited "
          "in Korea (별표 1, CAS 95-55-6) and is an EU-restricted hair dye"),
    Alert("para_phenylenediamine", "Prohapten (oxidation)", "high",
          "[NX3;H2,H1;!$(N[CX3]=[OX1,SX1,NX2]);!$(N[SX4])][c]1[c][c][c]([NX3;H2,H1;!$(N[CX3]=[OX1,SX1,NX2]);!$(N[SX4])])[c][c]1",
          "hair-dye allergen; oxidises to Bandrowski's base"),
    Alert("ortho_phenylenediamine", "Prohapten (oxidation)", "high",
          "[NX3;H2,H1;!$(N[CX3]=[OX1,SX1,NX2]);!$(N[SX4])][c]1[c]([NX3;H2,H1;!$(N[CX3]=[OX1,SX1,NX2]);!$(N[SX4])])[c][c][c][c]1",
          "oxidises to the ortho-diimine, the same route as the para isomer"),
    # The meta isomers are graded below the ortho/para pair on purpose, and
    # the reason is drawn rather than asserted: a quinoid needs the two
    # donors 1,2 or 1,4 so the ring can carry two exocyclic double bonds. At
    # 1,3 it cannot, so the direct oxidation route to a hapten is closed. They
    # are not clean - both are regulated hair-dye couplers and LLNA-positive,
    # reacting with the oxidised primary intermediate to give indo dyes - so
    # this is a lower grade, not an exemption.
    Alert("meta_aminophenol", "Prohapten (oxidation)", "moderate",
          "[NX3;H2,H1;!$(N[CX3]=[OX1,SX1,NX2]);!$(N[SX4])][c]1[c][c]([OX2H])[c][c][c]1",
          "1,3-substitution cannot form a quinone-imine directly; regulated "
          "as a hair-dye coupler and sensitising through the coupled dye"),
    Alert("meta_phenylenediamine", "Prohapten (oxidation)", "moderate",
          "[NX3;H2,H1;!$(N[CX3]=[OX1,SX1,NX2]);!$(N[SX4])][c]1[c][c]([NX3;H2,H1;!$(N[CX3]=[OX1,SX1,NX2]);!$(N[SX4])])[c][c][c]1",
          "1,3-substitution cannot form a diimine directly; regulated as a "
          "hair-dye coupler and sensitising through the coupled dye"),
    # --- Defect 3 domains: prehapten/prohapten fragrance chemistry ----------
    Alert("allyl_propenyl_phenol", "Prohapten (oxidation)", "high",
          "[OX2H][c]1[c][c][c]([$([CX4;H2][CX3;H1]=[CX3;H2]),"
          "$([CX3;H1]=[CX3;H1][CX4;H3])])[c][c]1",
          "4-allyl or 4-propenyl phenol; oxidises to the para-quinone "
          "methide, which is the hapten. Eugenol and isoeugenol are the "
          "classic pair and isoeugenol is among the strongest fragrance "
          "sensitisers known; both are EU Annex III labelled allergens. The "
          "alkene is spelled as a propenyl or allyl rather than any [CX3]= "
          "so this stays off the 4-hydroxystilbenes - resveratrol, "
          "oxyresveratrol and tapinarof are marketed actives"),
    Alert("isopropenyl_terpene", "Prehapten (autoxidation)", "moderate",
          "[CH2]=[CX3]([CH3])[CX4]",
          "limonene class; the storage hydroperoxide, not the parent, is the allergen"),
    Alert("tertiary_allylic_alcohol", "Prehapten (autoxidation)", "moderate",
          "[CX4]([OX2H])([CH3])[CX3]=[CX3]",
          "linalool class; autoxidises on air exposure to sensitising hydroperoxides"),





    Alert("allylic_alcohol", "Prehapten (autoxidation)", "moderate",
          "[CX3]=[CX3][CX4;H2][OX2H]",
          "primary allylic alcohol; the module carried only the tertiary "
          "(linalool) case and so missed geraniol, farnesol, cinnamyl and "
          "amylcinnamyl alcohol, all EU Annex III labelled allergens. "
          "Geraniol is the textbook prehapten: on air exposure it autoxidises "
          "to the aldehyde and to allylic hydroperoxides, and the storage "
          "hydroperoxide, not the parent, is the allergen. Read that "
          "literally: this alert is about what the material becomes on the "
          "shelf, and it is not a claim that the parent alcohol haptenates. "
          "Retinol carries the same group and is reported here for the same "
          "reason - it autoxidises readily, which is why it ships in airless "
          "opaque packaging with an antioxidant - but retinol is a "
          "documented irritant and NOT a documented contact allergen, so for "
          "a compound like it this is a packaging and stabilisation question "
          "that `filters.stability` already owns, at a grade that asks for "
          "confirmation at the use level rather than for an assay before use"),
    Alert("peroxide", "Radical initiator (autoxidation product)", "moderate",
          "[OX2][OX2]",
          "peroxide or endoperoxide. This module says twice that the storage "
          "hydroperoxide rather than the parent terpene is the allergen, so a "
          "molecule that already IS one should not read as burden zero here - "
          "ascaridole is the tea-tree oxidation product that sensitises, and "
          "benzoyl peroxide is a documented contact sensitiser. Graded below "
          "the defined electrophiles because radical haptenation is less "
          "direct, and because benzoyl peroxide is a marketed acne active at "
          "2.5-10 percent. `admet_alerts` reports the same group separately "
          "as an irritant; this is the sensitisation half"),
    Alert("formaldehyde_releaser", "Prohapten (formaldehyde release)", "high",
          "[NX3][CX4;H2][OX2H]",
          "N-hydroxymethyl (hemiaminal): DMDM hydantoin, imidazolidinyl and "
          "diazolidinyl urea, quaternium-15. The alert is about what the "
          "molecule liberates, not what it is - the hemiaminal hydrolyses in "
          "the formulation to free formaldehyde, and the formaldehyde is the "
          "hapten. These are among the most widely used cosmetic "
          "preservatives and among the most frequent causes of cosmetic "
          "allergic contact dermatitis; EU and Korea permit them with a "
          "'releases formaldehyde' declaration above 0.001 percent. Bronopol "
          "releases formaldehyde by a different (nitro/bromo decomposition) "
          "route and is deliberately not matched by this pattern"),
)














_MECHANISM_IS_REACTIVITY = {
    "KEAP1": {"Michael acceptor", "SN2"},
    "NFE2L2": {"Michael acceptor", "SN2"},
}

#: What each alert severity contributes to the sensitisation burden. Public,
#: because `score.score_screen` weighs the same alerts and used to carry its
#: own copy -- `{"high": 3.0, "moderate": 1.5}` -- which silently scored `low`
#: and `context` at 0.0 against this table's 0.5. One catalogue, one scale.
SEVERITY_WEIGHT = {"high": 3.0, "moderate": 1.5, "low": 0.5, "context": 0.5}
_SEVERITY_WEIGHT = SEVERITY_WEIGHT                     # pre-existing spelling


def _compiled() -> list[tuple[Alert, Chem.Mol]]:
    out = []
    for a in _ALERTS:
        patt = Chem.MolFromSmarts(a.smarts)
        if patt is None:
            log.error("bad SMARTS for alert %s: %s", a.name, a.smarts)
            continue
        out.append((a, patt))
    return out


_COMPILED = _compiled()
_PAINS = FilterCatalog.FilterCatalog(
    FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)

_CONJUGATED_FREE_ACID = Chem.MolFromSmarts("[CX3]=[CX3][CX3](=[OX1])[$([OX2H1]),$([O-])]")


def coverage_limitations(smiles: str) -> list[str]:
    """Name a catalogue gap without inventing a sensitisation hazard grade.

    Reactivity of a conjugated free acid/carboxylate is not established by its
    ester analogue's severity. Keep the observation distinct from that grade.
    """
    mol = assessable_molecule(smiles)
    if mol is not None and mol.HasSubstructMatch(_CONJUGATED_FREE_ACID):
        return ["conjugated free carboxylic acid/carboxylate: sensitisation severity "
                "is not calibrated in this catalogue; an empty scored alert list "
                "does not establish absence of reactivity. Review substance-specific "
                "evidence and intended-use assays; no new hazard grade is inferred."]
    return []


def assessable_molecule(smiles: str) -> Optional[Chem.Mol]:
    """A missing, invalid or hydrogen-only structure has no structural verdict."""
    if not isinstance(smiles, str) or not smiles.strip():
        return None
    mol = Chem.MolFromSmiles(smiles)
    return mol if mol is not None and mol.GetNumHeavyAtoms() > 0 else None


def sensitisation_alerts(smiles: str, target_gene: Optional[str] = None) -> list[dict]:
    """Structural alerts for skin sensitisation, with mechanistic context.

    `target_gene` lets the caller say which protein is being pursued. When the
    target's pharmacology depends on covalent reactivity, matching alerts are
    *annotated* with that fact and flagged `mechanism_at_target`. Their
    severity is unchanged: the hazard is to skin protein, and the intent
    behind the same chemistry does not lower it.
    """
    mol = assessable_molecule(smiles)
    if mol is None:
        return [{"name": "unassessed", "domain": "unknown", "severity": "unknown",
                 "note": "missing, invalid or no-heavy-atom structure; not assessed",
                 "n_matches": None, "assessment_status": "NOT_EVALUATED"}]

    exempt = _MECHANISM_IS_REACTIVITY.get((target_gene or "").upper(), set())
    hits: list[dict] = []
    for alert, patt in _COMPILED:
        matches = mol.GetSubstructMatches(patt)
        if not matches:
            continue
        severity = alert.severity
        note = alert.note
        mechanism = alert.domain in exempt
        if mechanism:
            note = (f"{note}; at {target_gene} this reactivity is also the "
                    "mechanism of action, so the series is not ruled out by it "
                    "-- but the sensitisation risk is unchanged and has to be "
                    "settled on dose, exposure and an assay").strip("; ")
        hits.append({"name": alert.name, "domain": alert.domain,
                     "severity": severity, "note": note,
                     "mechanism_at_target": mechanism,
                     "n_matches": len(matches)})
    return hits


def pains_alerts(smiles: str) -> list[str]:
    """PAINS substructures: frequent-hitter chemotypes that fake assay activity."""
    mol = assessable_molecule(smiles)
    if mol is None:
        return []
    return [e.GetDescription() for e in _PAINS.GetMatches(mol)]


def sensitisation_burden(alerts: list[dict]) -> Optional[float]:
    """0-10 concern score; None means no assessable structure was supplied."""
    if any(a.get("assessment_status") == "NOT_EVALUATED" for a in alerts):
        return None
    total = sum(_SEVERITY_WEIGHT.get(a["severity"], 0.0) for a in alerts)
    return round(min(10.0, total), 2)


def summarise(smiles: str, target_gene: Optional[str] = None) -> dict:
    if assessable_molecule(smiles) is None:
        return {"alerts": [], "pains": [], "burden": None,
                "assessment_status": "NOT_EVALUATED",
                "reason": "invalid_or_empty_structure",
                "verdict": "unknown sensitisation and PAINS risk; structure not assessable"}
    alerts = sensitisation_alerts(smiles, target_gene)
    pains = pains_alerts(smiles)
    burden = sensitisation_burden(alerts)
    limitations = coverage_limitations(smiles)
    if burden == 0 and not pains:
        verdict = "no electrophilic alert"
    elif burden >= 3.0:
        verdict = "protein-reactive; needs a sensitisation assay (DPRA/KeratinoSens/h-CLAT) before use"
    elif burden > 0:
        verdict = "weak reactivity alert; confirm at intended use level"
    else:
        verdict = "assay-interference risk only"
    if limitations:
        verdict += "; COVERAGE_LIMITATION: " + " ".join(limitations)
    return {"alerts": alerts, "pains": pains, "burden": burden, "verdict": verdict,
            "assessment_status": "EVALUATED", "reason": None,
            "coverage_limitations": limitations}
