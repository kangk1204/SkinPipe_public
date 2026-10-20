"""Structural alert tables for the topical ADMET module.

This file is deliberately only data. The reasoning that turns a match into a
verdict lives in `admet.py`; the chromophore/conjugation engine lives in
`admet_photo.py`; the `StructuralAlert` record and the compile/match helpers
live in `admet_alert_types.py` and are re-exported here, as is the phototoxic
chemotype table from `admet_photo_alerts.py`. Both splits are for file size
alone - import `PHOTO_CHEMOTYPES`, `compile_alerts` or `match_alerts` from
this module exactly as before.

Writing SMARTS that survive RDKit aromaticity perception
--------------------------------------------------------
RDKit re-perceives aromaticity on every parse, and it does not always agree
with the form a chemist would draw. Coumarin comes back with an aromatic
pyranone ring (`O=c1ccc2ccccc2o1`), while anthraquinone keeps its central ring
aliphatic (`O=C1c2ccccc2C(=O)c2ccccc21`) - two carbonyl-containing fused
systems, two different perceptions. An earlier module in this project shipped
an isothiazolinone SMARTS written in Kekule form and silently missed
methylisothiazolinone for exactly this reason.

A second trap sits underneath the first: RDKit's SMARTS `r<n>` primitive tests
the atom's *minimum* ring size, not membership of some ring of size n. Every
fusion carbon of psoralen belongs to both the furan and the benzo ring, so its
minimum ring size is 5 and `[#6;r6]` does not match it. A coumarin pattern
written with `r6` matches coumarin and silently misses every furocoumarin -
the one class this module exists to catch.

The house rule adopted here: specify ring skeletons by element alone (`[#6]`),
let the SMARTS ring-closure bonds enforce ring size, and write bonds as "any"
(`~`) wherever aromatic perception is uncertain. Every pattern in this file
is exercised against named positives and named negatives in
`tests/test_admet.py`; a pattern that fails to compile is dropped loudly at
import time rather than silently never matching.

Alert provenance: genotoxicity after Kazius, McGuire & Bursi (2005) J Med Chem
48:312 (toxicophores from 4337 Ames-tested compounds) and ICH M7(R2);
phototoxic chemotypes from SCCS/1647/22 and ICH S10 (2013); irritation from
functional-group pKa class and amphiphile structure against OECD TG 439; CYP
from mechanism-based-inactivation chemotypes and heme-coordinating azoles.
"""
from __future__ import annotations

from .admet_alert_types import (StructuralAlert, compile_alerts,  # noqa: F401
                                match_alerts)
from .admet_photo_alerts import PHOTO_CHEMOTYPES

__all__ = ["StructuralAlert", "compile_alerts", "match_alerts",
           "PHOTO_CHEMOTYPES", "PHOTO_COMPILED", "GENOTOX_ALERTS",
           "GENOTOX_COMPILED", "IRRITATION_ALERTS", "IRRITATION_COMPILED",
           "LEUKODERMA_ALERTS", "CYP_ALERTS", "CYP_COMPILED",
           "CMR_ALERTS", "CMR_COMPILED"]


# --------------------------------------------------------------------------
# Genotoxicity / Ames mutagenicity (Kazius-Bursi toxicophores + ICH M7)
# --------------------------------------------------------------------------
GENOTOX_ALERTS: tuple[StructuralAlert, ...] = (
    StructuralAlert(
        "aromatic_nitro", "[c][$([NX3](=O)=O),$([NX3+](=O)[O-])]", "high",
        "nitroreductase reduces it to a nitrenium ion that alkylates guanine "
        "N7/C8; the single strongest Kazius toxicophore",
        "Ames (OECD TG 471) required; avoid in a leave-on product"),
    StructuralAlert(
        "primary_aromatic_amine", "[NX3;H2;!$(N[!#6]);!$(N[CX3]=[OX1,SX1,NX2]);"
                                  "!$(N[SX4])][c;$(c1ccccc1)]", "high",
        "N-hydroxylation then O-esterification gives a nitrenium ion (the "
        "2-naphthylamine mechanism). Restricted to carbocyclic aromatic amines: "
        "aminopurines and aminopyrimidines do not follow this route",
        "Ames (OECD TG 471)"),
    StructuralAlert(
        "secondary_aromatic_amine",
        "[NX3;H1;!$(N[!#6]);!$(N[CX3]=[OX1,SX1,NX2]);!$(N[SX4])][c;$(c1ccccc1)]",
        "moderate",
        "the same nitrenium route, but N-substitution slows N-hydroxylation and "
        "the anilino-heteroaryl motif is common in Ames-negative kinase "
        "inhibitors (imatinib among them), so this is weaker evidence than the "
        "primary amine alert",
        "Ames (OECD TG 471) to settle it either way"),
    StructuralAlert(
        "aromatic_azo", "[c][NX2]=[NX2][#6]", "high",
        "azoreductase in skin flora cleaves it to two aromatic amines, which "
        "then follow the nitrenium route",
        "Ames on the reduced cleavage products, not only the parent"),
    StructuralAlert(
        "n_nitroso", "[NX3][NX2]=[OX1]", "critical",
        "nitrosamine; alpha-hydroxylation gives a diazonium that alkylates DNA. "
        "Can form in situ from a secondary amine plus a nitrosating agent",
        "not acceptable at any level; remove the secondary amine or the nitrite"),
    StructuralAlert(
        "epoxide", "[OX2r3]1[#6r3][#6r3]1", "high",
        "ring strain drives direct SN2 alkylation of guanine N7",
        "Ames (OECD TG 471) and an in-vitro micronucleus (OECD TG 487)"),
    StructuralAlert(
        "aziridine", "[NX3r3]1[#6r3][#6r3]1", "high",
        "as the epoxide, and slower to hydrolyse so it survives to reach DNA",
        "Ames (OECD TG 471)"),
    StructuralAlert(
        "alkyl_halide", "[CX4;H1,H2,H3][Cl,Br,I]", "moderate",
        "direct SN2 alkylating agent; potency tracks the leaving group "
        "(I > Br > Cl) and falls off sharply with steric bulk",
        "Ames (OECD TG 471)"),
    StructuralAlert(
        # The tail was `[#6,#8,#7]`, and `[#8]` matches an oxygen whatever its
        # charge -- so the deprotonated oxygen of an alkyl sulfate *salt* read
        # as an ester and sodium lauryl sulfate came back a high-severity
        # genotoxicity alert with an ICH M7 next step. A sulfate monoester
        # anion is not an alkylating agent: the negative charge kills the SN2.
        # `$([OX2][#6])` requires the second oxygen to be esterified, which
        # keeps EMS and dimethyl sulfate and drops the surfactant salts.
        "sulfonate_ester",
        "[CX4][OX2][SX4](=[OX1])(=[OX1])[#6,#7,$([OX2][#6])]", "high",
        "alkyl sulfonate/sulfate ester; ethyl methanesulfonate is the textbook "
        "SN2 mutagen and an ICH M7 cohort-of-concern impurity",
        "ICH M7 control at the TTC (1.5 ug/day) or below"),
    StructuralAlert(
        "azide", "[$(N=[N+]=[N-]),$([N-][N+]#N)]", "high",
        "sodium-azide-type; mutagenic via the metabolite "
        "L-azidoalanine in Ames TA1535",
        "Ames (OECD TG 471)"),
    StructuralAlert(
        "hydrazine", "[NX3][NX3]", "high",
        "oxidises to a diazonium/carbon radical that methylates DNA; "
        "hydrazine and its simple alkyl derivatives are rodent carcinogens",
        "Ames (OECD TG 471); hydrazine is CMR 1B in the EU"),
    StructuralAlert(
        "nitrogen_mustard", "[NX3]([CX4][CX4][Cl,Br,I])[CX4][CX4][Cl,Br,I]",
        "critical",
        "bifunctional alkylator via the aziridinium ion; cross-links DNA",
        "not acceptable in a cosmetic"),
    StructuralAlert(
        "alpha_halo_carbonyl", "[CX3](=[OX1])[CX4;H1,H2][Cl,Br,I]", "moderate",
        "activated alkyl halide; the carbonyl stabilises the SN2 transition state",
        "Ames (OECD TG 471)"),
    StructuralAlert(
        "n_hydroxy_amine", "[NX3;H1,H2][OX2H1]", "moderate",
        "hydroxylamine; already at the oxidation state that esterifies to a "
        "nitrenium ion, so it bypasses the activation step",
        "Ames (OECD TG 471)"),
    StructuralAlert(
        "para_quinone", "[#6X3]1(=[OX1])~[#6X3]~[#6X3]~[#6X3](=[OX1])~[#6X3]~[#6X3]1",
        "moderate",
        "redox-cycles to superoxide and hydroxyl radical, giving oxidative DNA "
        "damage and 8-OHdG rather than adducts; clastogenic more than mutagenic",
        "in-vitro micronucleus (OECD TG 487) reads this better than Ames"),
    StructuralAlert(
        "hydroquinone_catechol",
        "[$([OX2H][c]1[c][c][c]([OX2H])[c][c]1),$([OX2H][c]1[c]([OX2H])[c][c][c][c]1)]",
        "moderate",
        "air- and peroxidase-oxidises to the quinone, which is the actual "
        "clastogen; hydroquinone is Ames-equivocal but reproducibly induces "
        "micronuclei, and is why it is banned from leave-on EU cosmetics",
        "in-vitro micronucleus (OECD TG 487); check the regulatory annex"),
)

# ----------------------------------------------------------------------------
# Chemical leukoderma: tyrosinase-activated melanocyte cytotoxicity.
#
# Why this endpoint exists at all, and why it is graded `context`
# ---------------------------------------------------------------
# Tyrosinase oxidises a 4-substituted phenol to the ortho-quinone, and the
# quinone binds melanocyte protein and kills the cell. The depigmentation is
# permanent. This is not a theoretical pathway: rhododendrol was sold as a
# tyrosinase-inhibiting brightening active and caused leukoderma in 19,605 of
# an estimated 800,000 users, forcing Kanebo's 2013 recall; 4-tert-butylphenol
# is the classic occupational analogue and both are oxidised to ortho-quinones
# recovered as glutathione and cysteine adducts.
#
# The trap for a whitening programme is that this IS the tyrosinase-inhibitor
# chemotype. 4-Butylresorcinol is an approved inhibitor at 0.1-0.3 percent
# with a good clinical record, arbutin and kojic acid are in wide use, and
# this pipeline's own README ranks 4-butylresorcinol third on a tyrosinase
# screen. So the alert is graded `context` on the `resorcinol` precedent
# above: a mechanism worth knowing at the bench, decided on use level and
# clinical evidence, never on the substructure alone.
#
# The free ortho position in the pattern is the monophenolase step, not
# tuning: tyrosinase hydroxylates the carbon next to the phenol OH before it
# oxidises the catechol to the quinone, so a 2,6-blocked phenol is not a
# substrate. That is what separates the alert from tocopherol and BHT.
#
# It earns its place alongside the named-substance check in `regulatory`
# rather than duplicating it. Korea's 별표 1 already prohibits 4-tert-부틸페놀
# (98-54-4), 4-tert-부틸피로카테콜, 4-벤질옥시페놀 (103-16-2, monobenzone) and
# 4-에톡시페놀 (622-62-8) by name, and that list is matched by name and CAS.
# What a name list cannot do is recognise the mechanism in a molecule nobody
# has listed yet - which is exactly what rhododendrol was in 2013.
#
# Deliberately NOT added here: an endocrine-disruption alert. The motifs that
# carry an actual ED classification (4-nonylphenol, 4-tert-octylphenol,
# bisphenol A) are named substances already on 별표 1, while every structural
# generalisation of them collides with ingredients Korea permits at a stated
# limit - parabens, and benzophenone-3 as a notified UV filter to 5 percent.
# A substructure does not carry evidence about EATS-mediated activity the way
# it carries evidence about absorption or electrophilicity.
# ----------------------------------------------------------------------------
LEUKODERMA_ALERTS: tuple[StructuralAlert, ...] = (
    StructuralAlert(
        "leukoderma_4_substituted_phenol",
        "[OX2H][c]1[cH][c][c]([$([CX4]),$([OX2][CX4])])[c][c]1", "context",
        "4-alkyl-, 4-alkoxy- or 4-aralkyloxy-phenol (and the 4-substituted "
        "catechols) with a free ortho position: a tyrosinase substrate, "
        "oxidised in the melanocyte to an "
        "ortho-quinone that binds cellular protein and causes melanocyte "
        "death and permanent depigmentation - the rhododendrol (19,605 "
        "affected users, Kanebo 2013 recall) and 4-tert-butylphenol "
        "mechanism. This is also the tyrosinase-inhibitor chemotype, so the "
        "alert is a mechanism to weigh, not a finding: 4-butylresorcinol is "
        "approved at 0.1-0.3 percent and arbutin and kojic acid are in wide "
        "use. A structural alert is NOT a regulatory classification",
        "a use-level and clinical-evidence question, not a ban: check the "
        "named substance against 별표 1 (4-tert-부틸페놀, 4-벤질옥시페놀 and "
        "4-에톡시페놀 are prohibited by name), then ask for melanocyte "
        "cytotoxicity (B16 or normal human melanocyte, with and without "
        "tyrosinase) and a repeat-insult patch test at the intended level"),
)


# ----------------------------------------------------------------------------
# Skin irritation / corrosion. The acid/base pKa-class table lives in
# admet_irritation, which is the only consumer.
# ----------------------------------------------------------------------------
IRRITATION_ALERTS: tuple[StructuralAlert, ...] = (
    StructuralAlert(
        "peroxide", "[OX2][OX2]", "high",
        "generates radicals on contact with skin metals and keratinocyte "
        "peroxidases; benzoyl peroxide is efficacious and reliably irritant",
        "OECD TG 439 (RhE) and a patch-test protocol"),
    StructuralAlert(
        "aldehyde", "[CX3H1](=[OX1])[#6]", "moderate",
        "condenses with keratin lysine; irritation and sensitisation overlap "
        "here, and `filters.safety` owns the sensitisation half",
        "OECD TG 439"),
    StructuralAlert(
        "acyl_halide_anhydride",
        "[$([CX3](=[OX1])[F,Cl,Br,I]),$([CX3](=[OX1])[OX2][CX3]=[OX1])]",
        "high",
        "hydrolyses on skin to the acid plus a strong acid or a second "
        "equivalent; corrosive rather than merely irritant",
        "OECD TG 431 (corrosion) before TG 439"),
    StructuralAlert(
        "polyethoxylate", "[OX2][CX4;H2][CX4;H2][OX2][CX4;H2][CX4;H2][OX2]", "low",
        "PEG chain; non-ionic surfactants are the mildest class, but the chain "
        "carries a 1,4-dioxane and ethylene-oxide residual question",
        "specify a purified grade; irritation itself is usually low"),
    StructuralAlert(
        "phenol_free", "[OX2H][c]", "low",
        "free phenol; coagulates surface protein at high concentration but is "
        "mild at cosmetic levels, which is why resveratrol, EGCG and quercetin "
        "are all in marketed leave-on products. A note, not a finding",
        "note only; TG 439 if the use level is above a few percent"),
    StructuralAlert(
        "ortho_para_diphenol",
        "[$([OX2H][c]1[c]([OX2H])[c][c][c][c]1),$([OX2H][c]1[c][c][c]([OX2H])[c][c]1)]",
        "moderate",
        "catechol or hydroquinone; air-oxidises on skin to the quinone, which "
        "is both the irritant and the hapten. This is the pair of phenols that "
        "matters, not phenols in general",
        "OECD TG 439, and see filters.safety for the sensitisation reading"),
# Depigmentation is not irritation, but `admet.summarise` has no other local
# skin-toxicity channel, and a finding nobody sees is worse than one filed
# under an imperfect heading. Read the alert name, not the domain name.
) + LEUKODERMA_ALERTS

# --------------------------------------------------------------------------
# CYP inhibition and hERG. Low priority for a topical unless SED is material.
# --------------------------------------------------------------------------
CYP_ALERTS: tuple[StructuralAlert, ...] = (
    StructuralAlert(
        "methylenedioxyphenyl", "[OX2]1[CX4;H2][OX2][c]2[c]1[c][c][c][c]2",
        "high",
        "mechanism-based CYP3A4 inactivation through a carbene-iron complex "
        "(safrole, berberine, piperine)",
        "irreversible; matters only if SED reaches a systemic level"),
    StructuralAlert(
        "terminal_acetylene", "[CX2;H1]#[CX2]", "moderate",
        "oxidised to a ketene that acylates the CYP apoprotein; the "
        "ethinylestradiol mechanism",
        "irreversible CYP inactivation"),
    StructuralAlert(
        "imidazole_triazole",
        "[$([nX2;R1]1[cX3;R1][nX3,nX2;R1][cX3;R1][cX3;R1]1),"
        "$([nX2;R1]1[cX3;R1][nX2;R1][nX2,nX3;R1][cX3;R1]1),"
        "$([nX2;R1]1[nX2;R1][nX2,nX3;R1][cX3;R1][cX3;R1]1)]",
        "moderate",
        "sp2 nitrogen coordinates the heme iron; the azole-antifungal "
        "mechanism, reversible and broad across CYP isoforms",
        "reversible; a real interaction risk only at systemic exposure"),
    StructuralAlert(
        "thiophene_furan", "[$([sX2]),$([oX2;$(o1cccc1)])]", "low",
        "bioactivated to an S-oxide or an epoxide; a metabolism-dependent "
        "toxicity flag, not an inhibition flag", "note only"),
    StructuralAlert(
        "pyridine_nitrogen", "[nX2;r6]", "low",
        "weak heme coordination; niacinamide sits here and is a cosmetic "
        "staple at 5 percent, so this is a note, not a finding", "note only"),
)

# --------------------------------------------------------------------------
# CMR screening features. A motif does not establish a legal classification.
# --------------------------------------------------------------------------
CMR_ALERTS: tuple[StructuralAlert, ...] = (
    StructuralAlert(
        "retinoid", "[CH3][#6]1=[#6]([#6]=[#6][#6]=[#6])[#6]([CH3])([CH3])[#6][#6][#6]1",
        "moderate",
        "beta-ionone/polyene feature shared by some vitamin-A compounds and "
        "carotenoids; it does not identify a retinoid or establish human "
        "teratogenicity. Compound identity, route and exposure are required",
        "review compound-specific reproductive data and applicable use conditions"),
    StructuralAlert(
        "phthalate_diester",
        "[#6][OX2H0][CX3](=[OX1])[c]1[c][c][c][c][c]1[CX3](=[OX1])[OX2H0][#6]",
        "moderate",
        "ortho-phthalate with two carbon-bound ester oxygens; developmental "
        "read-across concern, with compound-dependent evidence. This core "
        "alone does not establish a CMR classification",
        "check compound identity and developmental toxicology; named references are evaluated separately"),
    StructuralAlert(
        "phthalate_monoester",
        "[#6][OX2H0][CX3](=[OX1])[c]1[c][c][c][c][c]1[CX3](=[OX1])[$([OX2H1]),$([O-])]",
        "context",
        "ortho-phthalate monoester, not a diester; parent-specific metabolite "
        "and reproductive evidence must be reviewed",
        "retain for compound-specific review; this alert does not establish safety"),
    StructuralAlert(
        "boron_acid",
        "[#5;!$([#5]~[!#8]);$([#5](~[#8])(~[#8])~[#8])]",
        "critical",
        "boron coordinated only to at least three oxygens, consistent with "
        "boric acid or a borate salt/ester; these groups require a high-priority "
        "identity and restriction check. B-C organoboron does not match this rule",
        "verify membership of EU Cosmetics Regulation Annex II entries 1395/1396; no legal identity is inferred"),
    StructuralAlert(
        "organoboron", "[#5;$([#5]~[#6]);$([#5]~[#8])]", "context",
        "oxygen-containing organoboron with a direct B-C bond; the "
        "boric-acid/borate classification cannot be inherited from this motif",
        "review compound-specific toxicology, exposure and legal identity; an alert is not a safety assessment"),
    StructuralAlert(
        "resorcinol", "[OX2H][c]1[c][c]([OX2H])[c][c][c]1", "context",
        "meta-diphenol; resorcinol itself is a thyroid disruptor and a "
        "sensitiser, yet 4-butylresorcinol is an approved tyrosinase inhibitor "
        "at 0.1-0.3 percent. Treat as a use-level question, not a ban",
        "respect the annex concentration limit"),
)


PHOTO_COMPILED = compile_alerts(PHOTO_CHEMOTYPES)
GENOTOX_COMPILED = compile_alerts(GENOTOX_ALERTS)
IRRITATION_COMPILED = compile_alerts(IRRITATION_ALERTS)
CYP_COMPILED = compile_alerts(CYP_ALERTS)
CMR_COMPILED = compile_alerts(CMR_ALERTS)
