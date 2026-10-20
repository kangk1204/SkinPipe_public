"""Scaffold hopping: generate structurally novel analogues of a known active.

The problem this solves
-----------------------
The rest of this pipeline is good at finding compounds that bind a skin
target, and most of what it finds is already published, already in a product,
or already claimed. A hit list of known actives is a literature review, not a
development programme. What a company needs from a known active is a molecule
that keeps the binding interaction and is nobody else's.

That is a scaffold hop (Schneider et al. 1999, Angew Chem Int Ed 38:2894, who
named it; Bohm, Flohr & Stahl 2004, Drug Discov Today Technol 1:217, who
formalised it): a compound with a different core that presents the same
interaction pattern. Operationally it is a disagreement between two
measurements -- **low 2D similarity, high 3D shape and pharmacophore
similarity** -- and this module exists to measure both and reward the gap.

A molecule that is 2D-similar to the parent is a trivial analogue: cheap to
make, easy to score well, and squarely inside the parent's patent claims. A
molecule that is 3D-dissimilar has probably lost the binding mode. Neither is
worth anything on its own, which is why `hop_quality` is a geometric mean of
the two axes rather than a sum: both have to be true at once.

================================================================================
WHAT THIS MODULE CANNOT DO
================================================================================
These limits are not caveats at the end of a report. They are the reason the
output is a shortlist for a chemist, not a decision.

**A computational scaffold hop is a hypothesis, not an analogue.** Everything
here measures resemblance between two molecules in isolation. There is no
protein in any of these calculations. Scaffold hopping loses potency far more
often than it keeps it -- the whole reason the technique has a literature is
that the success cases are worth publishing -- and neither shape overlap nor
docking nor free-energy calculation can promise that a new core reaches the
same contacts with the same geometry and pays the same desolvation cost. A
high hop_quality means "this is worth making and testing". Only an assay
settles whether it worked, and the base rate says most of these will not.

**Novelty is not patentability.** What this module measures is structural
distance from one parent molecule. Freedom to operate is a different question
entirely: it depends on the claims of every granted and pending patent in
every jurisdiction you care about, including Markush claims broad enough to
cover a compound nobody ever synthesised, and on prior art from anywhere in
the world. A molecule with a novel Murcko scaffold relative to kojic acid can
sit dead centre of someone else's composition-of-matter claim. This module
does not search patents, cannot search patents, and its `scaffold_changed`
flag carries no legal weight whatsoever. Run a real FTO search with a patent
attorney before anyone commits money.

**Generated molecules may be synthetically inaccessible.** Recombination in
particular produces structures that are perfectly valid on paper and have no
route at any price. Every candidate is scored with the Ertl & Schuffenhauer
SA_Score and filtered on it (see `hop_sanity`), but that score is a
similarity-to-known-chemistry heuristic, not a retrosynthesis. A beautiful hop
nobody can make is worth nothing, so the SA score is reported on every row
rather than hidden behind the filter.

**Most of what is generated is junk, and the yield is reported.** Both
generation routes have a high failure rate. `enumerate_hops` returns a census
of exactly what was discarded and why, and `rank_hops` carries the survivor
count. Do not quote a top-5 table from this module without the yield beside
it; the ratio is the honest description of the method. The measured yields on
this project's own seeds are in the table below, and two of the three are zero.

**hop_quality is not a probability of activity.** It deliberately punishes
non-novelty, so a trivial methyl analogue of kojic acid (which very probably
*is* active) scores 1.96, while an unrelated sugar (which certainly is not)
scores 3.57. The score answers "is this a scaffold hop worth making?", not "is
this likely to bind". Read it alongside the pharmacophore and shape columns,
never alone.

================================================================================

The two generation routes, and what each is for
-----------------------------------------------
**BRICS recombination** (Degen et al. 2008, ChemMedChem 3:1503) cuts the seed
at the 16 bond types that medicinal chemists actually break, then rebuilds
new molecules from those pieces -- optionally mixed with a supplied fragment
pool. CAN: rearrange a multi-fragment molecule into genuinely different
topologies, and graft the seed's recognition fragments onto foreign cores.
CANNOT: do anything useful with a small rigid molecule. Kojic acid yields two
BRICS fragments, so there is essentially nothing to recombine; the route is
built for imatinib-shaped molecules, and the validation in the tests shows
exactly that split. It has no notion of chemistry beyond bond-type matching,
which is where most of the junk comes from. Combinatorial growth is bounded by
`max_depth` and by a hard `n` cap on the generator.

**Bioisosteric replacement** (see `hop_bioisostere`) applies a curated,
cited table of classical and non-classical replacements by substructure
match. CAN: make small, defensible, traceable edits with a literature reason
to expect the interaction survives, and it works on small rigid molecules
where recombination does not. CANNOT: change the overall topology much --
these are single-group edits, so they usually produce modest 2D novelty, and
they only fire where a rule matches. A seed with no matching functional group
gets nothing.

The two are complementary and `enumerate_hops` runs both by default.

Measured yield on this project's three seeds
--------------------------------------------
Run with `n=200` and 10 conformers -- both shipped defaults -- and
`max_3d=25`, which is *not* a default: `rank_hops` ships `max_3d=50`. The
3D-scored column below is therefore the cap this particular run was given, not
a property of the method, and a default run would carry twice as many imatinib
candidates into 3D. "Usable" means it survived chemical sanity and the SA
cutoff; "credible" means hop_quality >= 5.5, the band that the decoy
calibration puts above an unrelated molecule.

    seed              raw   usable   3D-scored   credible   best score   decoy baseline
    kojic acid          7    6 (86%)        5        4        8.60        3.6 / 4.5 max
    oxyresveratrol      9    9 (100%)       4        0        4.87        2.3 / 4.0 max
    imatinib          109  109 (100%)      25        0        4.20        2.9 / 5.0 max

Read that table before trusting anything this module produces.

**Kojic acid works, and produces very little.** Four credible hops -- a
thionopyridinone, a deferiprone-like 3-hydroxy-4-pyridinone, a thiopyranone and
a tropolone -- all of which swap the chelating ring while keeping the O,O donor
pair. They are real suggestions. But the seed yielded only seven structures in
total, because BRICS gives two fragments and five bioisostere rules match. The
method's output here is a handful of chemist-obvious chelator swaps, arrived at
rigorously. That is worth something; it is not a screening library.

**Imatinib fails completely, and this is the important result.** 109 structures
were generated and every one passed sanity and SA -- the filters that exist to
catch junk caught nothing, 100%. Not one is a credible hop. The BRICS products
score pharmacophore similarities of 0.10-0.24 against a decoy ceiling of ~0.10:
recombination tore off the amide and the aminopyrimidine hinge binder, which is
the entire binding hypothesis, and kept the piperazine and the aryl rings, which
are the parts that do not matter. Meanwhile all six bioisosteric analogues score
exactly 0.0, because a single-group edit on a 37-heavy-atom molecule leaves
ECFP4 similarity at 0.65-0.91 -- above the same-chemotype threshold, and
correctly so. So for a large ligand with a distributed, multi-point
pharmacophore, recombination is too destructive and single-group replacement is
too conservative, and there is nothing in between here.

**Oxyresveratrol fails too, for a third reason**: only one bioisostere rule
(phenol -> F) matches a polyphenolic stilbene at all, and BRICS's three
fragments recombine into biphenyls that have deleted the vinyl linker holding
the two rings at the right separation.

The honest summary is that this module earns its keep on small chelator-type
actives and does not currently produce credible hops for large or
polyphenolic ligands. If you need hops for an imatinib-class molecule, the
missing ingredient is a pharmacophore constraint during generation -- growing
from the hinge-binding fragment rather than from every fragment equally -- and
that is not implemented here.

Determinism
-----------
Every step that touches randomness is seeded. Conformer embedding uses a fixed
ETKDGv3 seed (`hop_metrics.CONFORMER_SEED`); BRICS reagent scrambling uses the
`random` module, so `_seeded` saves and restores global RNG state around it
rather than leaving the process's RNG perturbed. Identical inputs give
byte-identical output.

Module layout
-------------
This file is the public surface and carries the contract above. The
implementation sits in siblings, kept separate only to respect the project's
file-length limit:

    hop_bioisostere.py   the curated, cited replacement table
    hop_metrics.py       2D / 3D / pharmacophore similarity and scaffolds
    hop_sanity.py        chemical sanity checks and the Ertl SA_Score
    hop_generate.py      fragmentation, BRICS recombination, hop enumeration
    hop_core.py          HopScore and the hop_quality formula
    hop_rank.py          ranking, the cost prescreen, and pipeline glue

Import everything from here; the siblings are implementation detail.

Worked validation on this project's seeds is in `tests/test_hopping.py`.
"""
from __future__ import annotations

import logging

from .hop_bioisostere import Bioisostere, applicable as bioisostere_rules_matching
from .hop_bioisostere import bioisostere_replace, table as bioisostere_table
from .hop_core import (SAME_SCAFFOLD_GATE, T_DISTINCT, T_SAME, W_NOVELTY, W_RETENTION,
                       HopScore, hop_score)
from .hop_generate import BRICS_SEED, enumerate_hops, fragment
from .hop_metrics import (CONFORMER_SEED, N_CONFORMERS_DEFAULT, best_shape_from_smiles,
                          murcko_scaffold,
                          pharmacophore_similarity, shape_tanimoto, tanimoto_2d)
from .hop_rank import rank_hops, to_ligand_records
from .hop_sanity import (SA_CUTOFF_DEFAULT, USING_REFERENCE_SA, filter_generated,
                         is_sane, sa_score)

log = logging.getLogger(__name__)

__all__ = [
    # scaffold analysis
    "murcko_scaffold",
    # generation
    "fragment", "enumerate_hops", "bioisostere_replace", "bioisostere_table",
    "bioisostere_rules_matching", "Bioisostere",
    # scoring
    "HopScore", "hop_score", "rank_hops",
    "tanimoto_2d", "shape_tanimoto", "best_shape_from_smiles", "pharmacophore_similarity",
    # sanity / accessibility
    "sa_score", "is_sane", "filter_generated", "USING_REFERENCE_SA",
    # pipeline glue
    "to_ligand_records",
    # calibration constants, exported so a report can state them
    "T_SAME", "T_DISTINCT", "W_NOVELTY", "W_RETENTION", "SAME_SCAFFOLD_GATE",
    "CONFORMER_SEED", "N_CONFORMERS_DEFAULT", "BRICS_SEED", "SA_CUTOFF_DEFAULT",
]
