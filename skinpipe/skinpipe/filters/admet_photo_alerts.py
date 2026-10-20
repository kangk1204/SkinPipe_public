"""Phototoxic chemotypes: the table `admet_alerts` re-exports and `admet` uses.

Split out of `admet_alerts` purely for file size; the alert contract, the
`StructuralAlert` record and the aromaticity house rules all live there and
apply here unchanged. Import `PHOTO_CHEMOTYPES` or `PHOTO_COMPILED` from either
module - `admet_alerts` re-exports both so existing importers keep working.

What belongs in this table
--------------------------
Classes for which epsilon above 1000 L/mol/cm in 290-700 nm is a foregone
conclusion and human phototoxicity is documented. A hit here overrides the
conjugation proxy in `admet_photo` in both directions: the proxy is an
inference about absorption, a chemotype is an observation about people.

Furocoumarins are handled in `admet_photo` instead, as a coumarin core fused
to a furan ring, because the linear (psoralen) and angular (angelicin) isomers
share no single tractable SMARTS.
"""
from __future__ import annotations

from .admet_alert_types import StructuralAlert

PHOTO_CHEMOTYPES: tuple[StructuralAlert, ...] = (
    StructuralAlert(
        "benzopyranone", "[#8]1~[#6](=[OX1])~[#6]~[#6]~[#6]2~[#6]~[#6]~[#6]~"
                         "[#6]~[#6]12",
        "high",
        "benzopyran-2-one; coumarin absorbs to 310 nm and the fused ellagic "
        "acid/HHDP lactones of the hydrolysable tannins reach 366 nm, so this "
        "fires on ellagitannins such as punicalagin as well as on coumarins "
        "proper. The absorption call is right in both cases; read the name as "
        "the ring system, not as a claim that the molecule is a coumarin",
        "OECD TG 432 (3T3 NRU) before any sun-exposed use"),
    StructuralAlert(
        "anthraquinone", "[#6]1~[#6]~[#6]~[#6]~[#6]2~[#6]1~[#6](=[OX1])~[#6]1~"
                         "[#6]~[#6]~[#6]~[#6]~[#6]1~[#6]2=[OX1]",
        "high",
        "9,10-anthracenedione; absorbs to ~430 nm and photosensitises "
        "type I/II ROS generation (emodin, hypericin, alizarin)",
        "OECD TG 432; consider a photoallergy (TG 432 + TG 442) package"),
    StructuralAlert(
        "phenothiazine", "[#7]1~[#6]2~[#6]~[#6]~[#6]~[#6]~[#6]2~[#16]~[#6]2~"
                         "[#6]~[#6]~[#6]~[#6]~[#6]12",
        "high",
        "chlorpromazine/promethazine class; the reference human phototoxic and "
        "photoallergic drug chemotype",
        "not usable in a leave-on topical product"),
    StructuralAlert(
        # The fused ring is written [#6,#7], not [#6]. Spelled as six carbons
        # this pattern is a QUINOLONE pattern only, and the 1,8-naphthyridones
        # are exactly the members that matter most: nalidixic acid is the
        # reference bullous-phototoxicity and pseudoporphyria drug, and it
        # dropped through to the proxy at "moderate" while ciprofloxacin, a
        # milder photosensitiser, matched and was graded "high".
        "quinolone_4_one", "[OX1]=[#6]1~[#6](~[#6](=[OX1])[OX2,OX1])~[#6]~[#7]"
                           "~[#6]2~[#6,#7]~[#6,#7]~[#6,#7]~[#6,#7]~[#6]12",
        "high",
        "4-quinolone / 1,8-naphthyridin-4-one 3-carboxylic acid; UVA-driven "
        "defluorination gives an aryl cation (lomefloxacin, ciprofloxacin), "
        "and the unfluorinated naphthyridones photosensitise too - nalidixic "
        "acid caused bullous phototoxicity and pseudoporphyria. The free "
        "3-carboxylic acid on the 4-one is what keeps this off unrelated "
        "fused heteroaromatic ketones such as flavone and chromone",
        "OECD TG 432"),
    StructuralAlert(
        "para_quinone", "[#6X3]1(=[OX1])~[#6X3]~[#6X3]~[#6X3](=[OX1])~[#6X3]~[#6X3]1",
        "high",
        "quinone; absorbs above 290 nm and redox-cycles to superoxide",
        "OECD TG 432 plus a ROS assay (ICH S10 reactive-oxygen-species test)"),
    StructuralAlert(
        "nitroaromatic", "[c][$([NX3](=O)=O),$([NX3+](=O)[O-])]",
        "high",
        "aryl nitro; strong n->pi* absorption above 300 nm, photoreduces to "
        "nitroso and hydroxylamine",
        "OECD TG 432; also a genotoxicity alert"),
    StructuralAlert(
        "azo", "[#6][NX2]=[NX2][#6]",
        "high",
        "azo chromophore; the whole point of the group is visible absorption, "
        "and trans-cis photoisomerisation generates radicals",
        "OECD TG 432; check the azo dye is not a banned aromatic-amine cleaver"),
    StructuralAlert(
        # Named for the scaffold the pattern actually matches. It carries no
        # halogen term, so it fired on unsubstituted salicylanilide too while
        # calling it "halogenated" -- a reader checking the hit against the
        # name found a compound that is not the thing named. The conservative
        # over-call is kept (the parent scaffold is the photoallergophore and a
        # cosmetic screen should raise it), but it is now labelled for what it
        # is, with the documented cases named in the mechanism.
        "salicylanilide",
        "[OX2H][c]1[c][c][c][c][c]1[CX3](=[OX1])[NX3][c]",
        "high",
        "salicylanilide scaffold; the halogenated members (bithionol, "
        "tribromsalan) caused the 1960s photoallergic-dermatitis epidemic from "
        "soaps, and the unhalogenated parent shares the chromophore",
        "banned or restricted in most cosmetic regulations; check the annex"),
    StructuralAlert(
        "halogenated_bisphenol",
        "[$([#17,#35,#53][c]1[c]([#8,#16;R0][c])[c][c][c][c]1),"
        "$([#17,#35,#53][c]1[c][c]([#8,#16;R0][c])[c][c][c]1),"
        "$([#17,#35,#53][c]1[c][c][c]([#8,#16;R0][c])[c][c]1)]",
        "high",
        "halogenated diaryl ether/sulfide. This is the one class the ICH S10 "
        "absorption screen does not catch on epsilon alone: bithionol and "
        "hexachlorophene caused a photoallergic-dermatitis epidemic from "
        "medicated soaps and are banned from cosmetics, and the mechanism is "
        "UVA photodehalogenation to an aryl radical, not a large extinction "
        "coefficient. Triclosan shares the scaffold and is permitted at 0.3 "
        "percent, so read the annex before treating this as a stopper",
        "OECD TG 432 plus a photoallergy read (TG 442 after UV dosing), and "
        "check the regulatory annex for the specific member"),
    StructuralAlert(
        "stilbene", "[c][CX3]=[CX3][c]",
        "moderate",
        "aryl-vinyl-aryl; resveratrol-type, lambda_max ~305-320 nm with "
        "epsilon ~30000, photoisomerises trans->cis",
        "TG 432 at the intended use level; UV-opaque packaging either way"),
    StructuralAlert(
        "aryl_enone", "[c][CX3]=[CX3][CX3]=[OX1]",
        "moderate",
        "chalcone/cinnamoyl; conjugation to the carbonyl pushes the band to "
        "~340 nm (curcumin, licochalcone)",
        "TG 432 at the intended use level"),
    StructuralAlert(
        "aryl_ketone", "[c][CX3](=[OX1])[#6]",
        "moderate",
        "aryl ketone; benzophenone-type n->pi* around 330 nm populates a "
        "triplet that abstracts hydrogen from lipid and DNA",
        "TG 432; benzophenones are separately restricted as UV filters"),
)
