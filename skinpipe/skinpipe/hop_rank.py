"""Ranking hops, with the skin permeation filter attached.

Split out of `hopping.py`, which re-exports everything here. The cost model
that motivates the prescreen is documented on `_cheap_prescreen`.
"""
from __future__ import annotations

import logging
import math
from typing import Iterable, Optional

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors

from .filters import topical
from .hop_core import T_DISTINCT, T_SAME, W_NOVELTY, W_RETENTION, hop_score
from .hop_metrics import (N_CONFORMERS_DEFAULT, embed_conformers,
                          pharmacophore_similarity, tanimoto_2d)
from .hop_sanity import SA_CUTOFF_DEFAULT, is_sane, sa_score
from .library import standardize
from .types import LigandRecord

RDLogger.DisableLog("rdApp.*")
log = logging.getLogger(__name__)

#: Bucket for a candidate that is a different spelling of one already seen.
#: It is counted in `prescreen_reasons` like any other rejection, inside
#: `n_prescreened_out`, so the accounting still closes:
#: n_input == n_prescreened_out + n_truncated_before_3d + n_scored.
DUPLICATE_REASON = "duplicate of an earlier candidate after canonicalisation"


def _cheap_prescreen(parent: str, candidate: str, max_sa: float) -> tuple[bool, str]:
    """The gates worth applying before paying for conformers and O3A.

    O3A alignment is the dominant cost in this module -- about 0.13 s per
    conformer pair on an imatinib-sized molecule, so ~13 s per candidate at
    ten conformers a side. Anything the 2D metrics can already rule out should
    never reach it. Every rejection here would have scored zero or near-zero
    anyway: `hop_quality` is already zero for a candidate at or above the
    same-chemotype threshold, so skipping its 3D calculation changes the
    ranking not at all and saves most of the runtime.
    """
    ok, reason = is_sane(candidate)
    if not ok:
        return False, f"failed sanity: {reason}"
    score = sa_score(candidate)
    if score is None or not math.isfinite(score) or not 1 <= score <= 10 or score > max_sa:
        return False, f"SA score {score} above cutoff {max_sa}"
    if tanimoto_2d(parent, candidate) >= T_SAME:
        return False, "2D similarity at or above the same-chemotype threshold"
    return True, ""


def _cheap_rank_estimate(parent: str, candidate: str) -> float:
    'hop_quality with the pharmacophore term standing in for retention.'
    novelty = max(0.0, min(1.0, (T_SAME - tanimoto_2d(parent, candidate))
                           / (T_SAME - T_DISTINCT)))
    pharm = pharmacophore_similarity(parent, candidate)
    pharm = max(0.0, pharm) if pharm == pharm else 0.0        # NaN-safe
    return (novelty ** (W_NOVELTY / (W_NOVELTY + W_RETENTION))
            * pharm ** (W_RETENTION / (W_NOVELTY + W_RETENTION)))


def rank_hops(parent_smiles: str, candidates: Iterable[str],
              n_confs: int = N_CONFORMERS_DEFAULT, max_sa: float = SA_CUTOFF_DEFAULT,
              max_3d: int = 50, prescreen: bool = True,
              min_topical_score: float = 0.0) -> pd.DataFrame:
    """Score and rank candidates, with the skin filter attached.

    Returns a DataFrame sorted by `hop_quality` descending, one row per
    candidate that reached 3D scoring, with the full score breakdown plus the
    permeation columns from `filters.topical` -- because a hop that cannot
    cross the stratum corneum is not a topical candidate no matter how elegant
    it is, and the two judgements have to be read together.

    Columns: smiles, hop_quality, tanimoto_2d, shape_tanimoto,
    pharmacophore_similarity, novelty, retention, same_murcko,
    scaffold_changed, generic_scaffold_changed, candidate_murcko, sa_score, mw,
    logp, tpsa, hbd, log_kp, permeation_class, topical_score, n_conf_parent,
    n_conf_candidate, notes.

    `prescreen=True` (the default) skips the expensive 3D step for candidates
    that are already insane, synthetically hopeless, or 2D-similar enough that
    `hop_quality` would be zero regardless. `max_3d` then caps how many of the
    survivors are 3D-scored, ordered by `_cheap_rank_estimate` -- novelty and
    pharmacophore, both free -- so a large candidate list degrades to its most
    promising subset rather than running for an hour. All of it is reported in
    the DataFrame's `.attrs`: `attrs["n_input"]`, `attrs["n_prescreened_out"]`,
    `attrs["n_truncated_before_3d"]`, `attrs["n_scored"]`,
    `attrs["n_duplicate_candidates"]`, `attrs["prescreen_reasons"]`.

    Candidates are identified canonically, not textually. `'Oc1ccccc1'` and
    `'c1ccc(O)cc1'` are one molecule and produce one row; the second is counted
    in `attrs["n_duplicate_candidates"]` and named in
    `attrs["prescreen_reasons"]`. This used to be textual, which meant a caller
    whose candidate list came from two sources got the same compound scored
    twice, ranked twice, and handed on as two `LigandRecord`s -- and, worse,
    the duplicate consumed one of the `max_3d` slots, so a genuinely different
    molecule further down the list was never 3D-scored at all. `n_input` still
    counts what was handed in, so the accounting closes:

        n_input == n_prescreened_out + n_truncated_before_3d + n_scored

    unconditionally. A caller passing already-distinct candidates sees no
    change.

    That invariant used to be quoted as the two-term
    ``n_input == n_prescreened_out + n_scored``, which is only true when
    `max_3d` never binds. Every candidate dropped by the truncation below left
    the run with no bucket at all: on this project's own imatinib figures (109
    generated, 103 surviving the prescreen) and the shipped default of
    `max_3d=50`, that is 53 candidates -- 51% of the survivors -- reported
    nowhere, while `hopping.py` insists the yield is the honest description of
    the method. They are counted in `n_truncated_before_3d` now. The
    `min_topical_score` caveat the old wording carried was never needed
    either: that filter runs on the assembled frame, after `n_scored` is
    fixed, so it changes `len(frame)` and not the census.

    Note that truncation at `max_3d` means the returned table is not the whole
    candidate set ranked; it is the best `max_3d` by a cheap proxy, then ranked
    properly. `attrs["n_truncated_before_3d"]` is how many the proxy alone
    decided about, and a non-zero value is the signal to raise `max_3d` if you
    need exhaustive scoring and can pay for it.

    An empty DataFrame (every candidate rejected) is a legitimate and common
    result. It is returned with the same columns so downstream code does not
    have to special-case it.
    """
    parent = standardize(parent_smiles) or parent_smiles
    # Materialise first: `candidates` may be a generator, and counting it after
    # consuming it gives zero.
    raw = list(candidates)
    n_input = len(raw)
    # `dict.fromkeys` collapses candidates spelled identically. It used to run
    # before anything counted them, so two copies of one string arrived as
    # `n_input=1` with `n_duplicates=0` -- the input shrank and nothing said
    # why. They are counted here instead, in the same bucket as the canonical
    # duplicates found below, so input = dropped + scored holds.
    candidate_list = [c for c in dict.fromkeys(raw) if c]
    n_blank = sum(1 for c in raw if not c)
    n_textual_dupes = len(raw) - len(candidate_list) - n_blank

    survivors: list[str] = []
    reasons: dict[str, int] = {}
    seen: set[str] = set()
    n_duplicates = n_textual_dupes
    if n_textual_dupes:
        reasons[DUPLICATE_REASON] = reasons.get(DUPLICATE_REASON, 0) + n_textual_dupes
    if n_blank:
        reasons["empty input"] = reasons.get("empty input", 0) + n_blank
    for candidate in candidate_list:
        standardised = standardize(candidate)
        if standardised is None or standardised == parent:
            reasons["unparseable or identical to parent"] = \
                reasons.get("unparseable or identical to parent", 0) + 1
            continue
        # Identity is canonical, not textual. `dict.fromkeys` above only
        # collapses candidates spelled the same way, so 'Oc1ccccc1' and
        # 'c1ccc(O)cc1' both reached the 3D stage as separate molecules --
        # two rows, an inflated n_scored, and two of the max_3d conformer
        # slots spent overlaying one phenol against itself. The check sits
        # before the prescreen because a duplicate is a duplicate whatever
        # the prescreen would have said about it, and running the gates twice
        # on the same molecule only pays for the same answer.
        if standardised in seen:
            n_duplicates += 1
            reasons[DUPLICATE_REASON] = reasons.get(DUPLICATE_REASON, 0) + 1
            continue
        seen.add(standardised)
        if prescreen:
            ok, reason = _cheap_prescreen(parent, standardised, max_sa)
            if not ok:
                reasons[reason.split(":")[0]] = reasons.get(reason.split(":")[0], 0) + 1
                continue
        survivors.append(standardised)

    if n_duplicates:
        log.info("%d of %d candidates were other spellings of a molecule already "
                 "in the list; collapsed before 3D scoring", n_duplicates,
                 len(candidate_list))

    survivors.sort(key=lambda s: -_cheap_rank_estimate(parent, s))
    n_survivors = len(survivors)
    if n_survivors > max_3d:
        # Truncate first, then measure, so the count is the difference between
        # two lists that exist rather than an arithmetic guess that a negative
        # or zero `max_3d` would get wrong.
        survivors = survivors[:max(max_3d, 0)]
        log.info("%d candidates passed prescreen; 3D-scoring the %d most promising "
                 "and dropping %d unscored", n_survivors, len(survivors),
                 n_survivors - len(survivors))
    n_truncated = n_survivors - len(survivors)

    parent_ensemble = embed_conformers(parent, n_confs)
    if parent_ensemble is None:
        log.warning("could not embed the parent %s; shape scores will be unavailable",
                    parent[:60])

    rows: list[dict] = []
    for candidate in survivors:
        score = hop_score(parent, candidate, n_confs=n_confs, parent_mol=parent_ensemble)
        mol = Chem.MolFromSmiles(candidate)
        mw, logp = Descriptors.MolWt(mol), Descriptors.MolLogP(mol)
        hbd, tpsa = Descriptors.NumHDonors(mol), Descriptors.TPSA(mol)
        log_kp = topical.potts_guy_log_kp(logp, mw)
        rows.append({
            "smiles": candidate,
            "hop_quality": score.hop_quality,
            "tanimoto_2d": score.tanimoto_2d,
            "shape_tanimoto": score.shape_tanimoto,
            "pharmacophore_similarity": score.pharmacophore_similarity,
            "novelty": score.novelty,
            "retention": score.retention,
            "same_murcko": score.same_murcko,
            "scaffold_comparable": score.scaffold_comparable,
            "scaffold_changed": score.scaffold_changed,
            "generic_scaffold_changed": score.generic_scaffold_changed,
            "parent_murcko": score.parent_murcko,
            "candidate_murcko": score.candidate_murcko,
            "sa_score": score.sa_score,
            "mw": round(mw, 2), "logp": round(logp, 2),
            "tpsa": round(tpsa, 2), "hbd": hbd,
            "log_kp": round(log_kp, 2),
            "permeation_class": topical.permeation_class(log_kp, mw),
            "topical_score": topical.topical_score(mw, logp, hbd, tpsa),
            "n_conf_parent": score.n_conformers[0],
            "n_conf_candidate": score.n_conformers[1],
            "notes": "; ".join(score.notes),
        })

    columns = ["smiles", "hop_quality", "tanimoto_2d", "shape_tanimoto",
               "pharmacophore_similarity", "novelty", "retention", "same_murcko",
               "scaffold_comparable", "scaffold_changed", "generic_scaffold_changed",
               "parent_murcko", "candidate_murcko",
               "sa_score", "mw", "logp", "tpsa", "hbd", "log_kp", "permeation_class",
               "topical_score", "n_conf_parent", "n_conf_candidate", "notes"]
    frame = pd.DataFrame(rows, columns=columns)
    frame["parent_smiles"] = parent
    if min_topical_score > 0.0 and not frame.empty:
        before = len(frame)
        frame = frame[frame["topical_score"] >= min_topical_score]
        log.info("topical filter at %.1f: %d -> %d rows", min_topical_score,
                 before, len(frame))
    if not frame.empty:
        frame = frame.sort_values("hop_quality", ascending=False).reset_index(drop=True)

    # `n_input` is what the caller handed over, not what survived deduplication:
    # the invariant this module documents is input == dropped + scored, and
    # reporting the deduplicated count made it hold by hiding the drop.
    frame.attrs.update({"parent_smiles": parent, "n_input": n_input,
                        "n_prescreened_out": n_input - n_survivors,
                        # Survivors the `max_3d` cap refused to pay for. Not a
                        # rejection -- nothing was measured about them -- which
                        # is why it is its own bucket and not a prescreen reason.
                        "n_truncated_before_3d": n_truncated,
                        "n_scored": len(rows),
                        "n_duplicate_candidates": n_duplicates,
                        "prescreen_reasons": reasons})
    return frame


# ------------------------------------------------------------- pipeline glue

def to_ligand_records(frame: pd.DataFrame, parent_id: str,
                      derivations: Optional[dict[str, str]] = None,
                      prefix: str = "HOP") -> list[LigandRecord]:
    """Turn ranked hops into LigandRecords the rest of the pipeline can consume.

    `parent_id` and `derivation` are written into every record so a generated
    structure can never be mistaken downstream for a compound with real assay
    data behind it. `activity` is deliberately left empty: these molecules have
    no measured activity, and populating it with the parent's would be a
    fabrication.
    """
    records: list[LigandRecord] = []
    for i, row in enumerate(frame.to_dict("records")):
        smiles = row["smiles"]
        records.append(LigandRecord(
            ligand_id=f"{prefix}_{parent_id}_{i:03d}", smiles=smiles,
            name=None, source="scaffold_hop",
            parent_id=parent_id,
            parent_smiles=row.get("parent_smiles") or frame.attrs.get("parent_smiles"),
            hop_metrics={k: row.get(k) for k in ("hop_quality", "tanimoto_2d",
                          "shape_tanimoto", "pharmacophore_similarity")},
            parent_metadata_source="computational scaffold generation",
            derivation=(derivations or {}).get(smiles, "scaffold_hop"),
            props={"mw": row["mw"], "logp": row["logp"], "tpsa": row["tpsa"],
                   "hbd": row["hbd"]},
            flags=["generated: hypothesis only, no assay data"]))
    return records
