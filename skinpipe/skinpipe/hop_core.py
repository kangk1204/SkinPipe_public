"""HopScore and the hop-quality formula -- the measurement core of `hopping`.

Split out of `hopping.py` so that module can carry the honesty contract for
the whole feature without running past the project's file-length limit. Import
from `..hopping`, which re-exports everything here; this module is the
implementation, not the public surface.
"""
from __future__ import annotations

import logging
import math
from dataclasses import asdict, dataclass, field
from typing import Optional

from rdkit import Chem, RDLogger

from .hop_metrics import (N_CONFORMERS_DEFAULT, embed_conformers, murcko_scaffold,
                          pharmacophore_similarity, shape_tanimoto, tanimoto_2d)
from .hop_sanity import sa_score
from .library import standardize

RDLogger.DisableLog("rdApp.*")
log = logging.getLogger(__name__)

# --- hop_quality calibration -------------------------------------------------
# At or above T_SAME, ECFP4 says "same chemotype" and novelty credit is zero.
# At or below T_DISTINCT, the fingerprints have essentially nothing in common
# and novelty credit is full. Between them the credit ramps linearly.
T_SAME, T_DISTINCT = 0.60, 0.25
# Exponents in the weighted geometric mean of novelty and retention. Retention
# counts double because novelty is free -- every random molecule has it -- while
# retention is the scarce property that carries the actual hypothesis. Measured
# consequence: with equal weights, a sugar scores 4.6 against kojic acid purely
# on its automatic novelty; at 1:2 it scores 3.6 and real hops stay above 5.5.
W_NOVELTY, W_RETENTION = 1.0, 2.0
# Sharing the parent's Murcko scaffold is not a hop by definition. The score is
# multiplied down rather than zeroed so the column stays a usable ranking
# signal among same-scaffold analogues, but nothing gated this way can reach
# the top of a mixed list.
SAME_SCAFFOLD_GATE = 0.30


def _clip01(x: float) -> float:
    return max(0.0, min(1.0, x))


@dataclass
class HopScore:
    'The full evidence for one parent/candidate pair. Nothing is hidden.'
    parent_smiles: str
    candidate_smiles: str
    tanimoto_2d: float                       # ECFP4; LOWER is better for novelty
    shape_tanimoto: Optional[float]          # 3D after O3A; HIGHER is better
    pharmacophore_similarity: float          # Gobbi 2D pharm fp; HIGHER is better
    same_murcko: bool
    scaffold_changed: bool
    hop_quality: float                       # 0-10; NaN when unscoreable, see above
    novelty: float = 0.0                     # the 2D half, 0-1, after the ramp
    retention: float = 0.0                   # the 3D/pharm half, 0-1; NaN if no shape
    parent_murcko: str = ""
    candidate_murcko: str = ""
    #: False when either molecule is acyclic and has no Murcko scaffold. Read
    #: it before `scaffold_changed`: that field is false in two different
    #: situations and only one of them is a measurement.
    scaffold_comparable: bool = True
    generic_scaffold_changed: bool = False   # topology change, not just an element swap
    sa_score: Optional[float] = None         # Ertl & Schuffenhauer, 1 easy .. 10 hopeless
    n_conformers: tuple[int, int] = (0, 0)   # (parent, candidate) after RMSD pruning
    shape_measured: bool = True              # False once embedding or MMFF failed
    hop_quality_ceiling: Optional[float] = None   # set only when NOT shape_measured
    notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)


# ------------------------------------------------------------------ scoring

def hop_score(parent_smiles: str, candidate_smiles: str,
              n_confs: int = N_CONFORMERS_DEFAULT,
              parent_mol: Optional[Chem.Mol] = None) -> HopScore:
    """Score one candidate against the parent on all three axes.

    `parent_mol` accepts a pre-embedded conformer ensemble from
    `hop_metrics.embed_conformers`, so `rank_hops` can build the parent's
    conformers once instead of once per candidate. Passing it changes nothing
    about the result, only the cost.

    Deterministic: with a fixed conformer seed the same pair always returns the
    same numbers.

    If the candidate (or the parent) will not embed, `shape_measured` is False
    and `hop_quality` is NaN rather than a pharmacophore-only score, so a
    failed measurement cannot outrank a successful one. See `HopScore` for the
    measurement that forced this.

    >>> s = hop_score("OCC1=CC(=O)C(O)=CO1", "OCC1=CC(=O)C(O)=C(C)O1")
    >>> s.same_murcko                     # a ring methyl is not a hop
    True
    >>> s.hop_quality < 3.0
    True
    """
    parent = standardize(parent_smiles) or parent_smiles
    candidate = standardize(candidate_smiles) or candidate_smiles
    notes: list[str] = []

    parent_scaffold = murcko_scaffold(parent)
    candidate_scaffold = murcko_scaffold(candidate)
    # An acyclic molecule has no Murcko scaffold and returns "". Two empty
    # strings are not a scaffold match, they are two absences, so say so
    # instead of letting `"" == ""` quietly assert sameness.
    # Three states, not two. `CCO` and `CCN` are both acyclic, so both Murcko
    # scaffolds are the empty string -- and `same_murcko = False` reported
    # `scaffold_changed: true` for a comparison that never happened, then the
    # branch below explained it as a heteroatom swap. A boolean consumer reads
    # that as a measured scaffold hop. The note said otherwise and a note is
    # not a field.
    scaffold_comparable = bool(parent_scaffold and candidate_scaffold)
    if not scaffold_comparable:
        # No evidence of a change, which is what "not comparable" licenses --
        # not evidence of one. `SAME_SCAFFOLD_GATE` then declines the pair
        # rather than passing it as a hop.
        same_murcko = False
        scaffold_changed = False
        notes.append("at least one molecule is acyclic and has no Murcko "
                     "scaffold, so the scaffold gate could not be applied: "
                     "`scaffold_changed` is false because nothing was "
                     "compared, not because the scaffold is the same")
    else:
        same_murcko = parent_scaffold == candidate_scaffold
        scaffold_changed = not same_murcko

    generic_changed = (murcko_scaffold(parent, generic=True)
                       != murcko_scaffold(candidate, generic=True))
    # `scaffold_comparable` is redundant here while the branch above keeps
    # `scaffold_changed` false for an incomparable pair, and it stays as the
    # second line: the note is the part a reader believes, and a future edit
    # to that branch must not be able to resurrect it. A mutation removing
    # this guard alone is therefore equivalent, which is the point.
    if scaffold_comparable and scaffold_changed and not generic_changed:
        notes.append("the Murcko scaffold changed but the generic framework did not: "
                     "this is a heteroatom swap, not a topological hop")

    t2d = tanimoto_2d(parent, candidate)
    pharm = pharmacophore_similarity(parent, candidate)

    ensemble_parent = parent_mol if parent_mol is not None else embed_conformers(parent, n_confs)
    ensemble_candidate = embed_conformers(candidate, n_confs)
    shape = shape_tanimoto(ensemble_parent, ensemble_candidate)
    n_conf = (ensemble_parent.GetNumConformers() if ensemble_parent is not None else 0,
              ensemble_candidate.GetNumConformers() if ensemble_candidate is not None else 0)

    pharm_safe = max(0.0, pharm) if pharm == pharm else 0.0        # NaN-safe
    novelty = _clip01((T_SAME - t2d) / (T_SAME - T_DISTINCT)) if t2d == t2d else 0.0
    gate = SAME_SCAFFOLD_GATE if same_murcko else 1.0
    total_weight = W_NOVELTY + W_RETENTION

    def _quality(retention_value: float) -> float:
        return 10.0 * gate * (novelty ** (W_NOVELTY / total_weight)
                              * retention_value ** (W_RETENTION / total_weight))

    shape_measured = shape is not None
    ceiling: Optional[float] = None
    if shape_measured:
        # Geometric mean: shape alone cannot carry a candidate, because its
        # baseline against unrelated molecules is high and size-dependent.
        retention = math.sqrt(max(0.0, shape) * pharm_safe)
        quality = _quality(retention)
    else:
        # Retention is the geometric mean of two terms and one of them is
        # missing, so there is no retention number -- unless the pharmacophore
        # term is itself zero, which drives the product to zero for every
        # possible shape.
        retention = 0.0 if pharm_safe == 0.0 else math.nan
        # The score a perfect overlay (shape = 1.0) would have produced. It is
        # the whole range the missing measurement could have moved within, so
        # when it is zero the measurement could not have mattered.
        best_case = _quality(math.sqrt(pharm_safe))
        if best_case == 0.0:
            certain = [reason for reason, holds in
                       (("novelty credit is already zero", novelty == 0.0),
                        ("the pharmacophore similarity is 0.00", pharm_safe == 0.0))
                       if holds]
            quality = 0.0
            notes.append("3D shape could not be measured (embedding or MMFF failed), "
                         "but the verdict does not depend on it: "
                         + " and ".join(certain) +
                         ", so hop_quality is 0.00 for every possible shape overlap")
        else:
            quality = math.nan
            ceiling = round(best_case, 2)
            notes.append(
                "3D shape could not be measured (embedding or MMFF failed), so this "
                "candidate is UNSCOREABLE: hop_quality is NaN, which is neither a low "
                "score nor a high one. Scoring it on the pharmacophore term alone "
                "would reward the failed measurement -- that is arithmetically the "
                "same as assuming shape equals pharmacophore similarity, and it "
                f"overstates any candidate whose real overlap is below "
                f"{pharm_safe:.2f}. Even a perfect overlay could not take it above "
                f"{ceiling:.2f}; re-run with more conformers or drop the row, but do "
                "not rank it against candidates that were actually measured")

    if same_murcko:
        notes.append("shares the parent's Murcko scaffold, so it is an analogue "
                     "rather than a hop; hop_quality is gated down accordingly")
    if novelty == 0.0 and t2d == t2d:
        notes.append(f"ECFP4 Tanimoto {t2d:.2f} is at or above the {T_SAME} "
                     "same-chemotype threshold: no novelty credit")

    return HopScore(
        parent_smiles=parent, candidate_smiles=candidate,
        tanimoto_2d=t2d, shape_tanimoto=shape, pharmacophore_similarity=pharm,
        same_murcko=same_murcko, scaffold_changed=scaffold_changed,
        scaffold_comparable=scaffold_comparable,
        hop_quality=round(quality, 2), novelty=round(novelty, 4),
        retention=round(retention, 4), parent_murcko=parent_scaffold,
        candidate_murcko=candidate_scaffold, generic_scaffold_changed=generic_changed,
        sa_score=sa_score(candidate), n_conformers=n_conf,
        shape_measured=shape_measured, hop_quality_ceiling=ceiling, notes=notes)
