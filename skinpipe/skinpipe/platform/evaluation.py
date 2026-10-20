"""Evidence-bounded benchmark metrics and leakage checks for frozen rankings.

Unmeasured candidates are never converted to inactive labels. Rank metrics use
measured rows, while the prospective top-k report retains missing outcome counts.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from copy import deepcopy
from datetime import datetime
import math
import operator
import random
from statistics import mean
from typing import Any, Callable

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator

from .identity import chemical_identity
from .ranking import finite


def _label(row):
    label = row.get("label")
    if row.get("label_measured") is not True or isinstance(label, (bool, np.bool_)):
        return None
    try:
        label = operator.index(label)
    except TypeError:
        return None
    return label if label in {0, 1} else None


def _score(row):
    value = finite(row.get("score", row.get("priority_score_v2")))
    if value is None:
        raise ValueError("benchmark scores must be finite")
    return value


def _average_precision(y, scores):
    positives = sum(y)
    if not positives:
        return None
    order = sorted(range(len(y)), key=lambda i: -scores[i])
    tp = previous_tp = 0
    ap = 0.0
    for rank, i in enumerate(order, 1):
        tp += y[i]
        if rank == len(order) or scores[order[rank]] != scores[i]:
            ap += (tp - previous_tp) / positives * tp / rank
            previous_tp = tp
    return ap


def _bedroc(y, scores, alpha):
    n, positives = len(y), sum(y)
    if not positives or positives == n:
        return None
    order = sorted(range(n), key=lambda i: -scores[i])
    # Subtracting the common 1 cancels in the normalized ratio. expm1 keeps
    # tiny positive alpha informative; exp retains range for large alpha.
    transform = math.expm1 if alpha < 1 else math.exp
    weights = [transform(-alpha * rank / n) for rank in range(1, n + 1)]
    start = 0
    actual = 0.0
    while start < n:
        end = start + 1
        while end < n and scores[order[end]] == scores[order[start]]:
            end += 1
        actual += mean(weights[start:end]) * sum(y[order[i]] for i in range(start, end))
        start = end
    perfect, worst = sum(weights[:positives]), sum(weights[-positives:])
    return min(1.0, max(0.0, (actual - worst) / (perfect - worst)))


def evaluate_ranked_candidates(rows: list[dict[str, Any]], *, k: int = 96,
                               alpha: float = 20.0) -> dict[str, Any]:
    """Evaluate one target; controls are reported separately and never new hits.

    `hit_at_k` is the count of confirmed hits among the frozen top-k. The measured
    hit rate and follow-up coverage are separate. `pr_auc` is threshold-group
    average precision (step PR area), not trapezoidal interpolation.
    """
    if isinstance(k, bool) or not isinstance(k, int) or k < 1:
        raise ValueError("k must be a positive integer")
    if finite(alpha) is None or not 0 < alpha <= 1000:
        raise ValueError("alpha must be finite and in (0, 1000]")
    targets = {r.get("target_id") for r in rows}
    if len(targets) > 1 or None in targets:
        raise ValueError("evaluate one explicit target at a time")
    ids = [r.get("candidate_id") for r in rows]
    if any(not isinstance(cid, str) or not cid for cid in ids) or len(set(ids)) != len(ids):
        raise ValueError("benchmark rows require unique nonempty candidate IDs")
    control_roles = {"control", "positive_control", "negative_control", "counter_control"}
    def is_control(row):
        return row.get("role") in control_roles or row.get("opposing_direction") is True
    controls = [r for r in rows if is_control(r)]
    candidates = [r for r in rows if not is_control(r)]
    # Parent dedup is a benchmark contract; unknown parents remain disclosed.
    parents = [r["parent_id"] for r in candidates if r.get("parent_id")]
    if len(set(parents)) != len(parents):
        raise ValueError("duplicate parents inflate benchmark candidate counts")
    candidates.sort(key=lambda r: (-_score(r), r["candidate_id"]))
    top = candidates[:k]
    observed = [r for r in candidates if _label(r) is not None]
    observed_top = [r for r in top if _label(r) is not None]
    hits = [r for r in observed_top if _label(r) == 1]
    y, scores = [_label(r) for r in observed], [_score(r) for r in observed]
    new_cores = {r["scaffold_id"] for r in hits if r.get("novel_scaffold") is True
                 and r.get("scaffold_kind") == "murcko" and r.get("scaffold_id")}
    positive_count = sum(y)
    ef = None
    if observed and positive_count:
        cutoff = max(1, math.ceil(len(observed) * .01))
        threshold = sorted(scores, reverse=True)[cutoff - 1]
        above = [i for i, s in enumerate(scores) if s > threshold]
        tied = [i for i, s in enumerate(scores) if s == threshold]
        count = sum(y[i] for i in above) + sum(y[i] for i in tied) * (cutoff - len(above)) / len(tied)
        ef = count / cutoff / (positive_count / len(observed))
    procured = [r for r in top if type(r.get("procurement_success")) is bool]
    return {"target_id": next(iter(targets), None), "k": k, "selected_count": len(top),
        "candidate_count": len(candidates), "control_count": len(controls),
        "measured_count": len(observed), "unmeasured_count": len(candidates) - len(observed),
        "invalid_measured_label_count": sum(r.get("label_measured") is True and _label(r) is None for r in candidates),
        "top_k_measured_count": len(observed_top), "hit_at_k": len(hits),
        "hit_rate_among_measured_top_k": len(hits) / len(observed_top) if observed_top else None,
        "top_k_outcome_coverage": len(observed_top) / len(top) if top else None,
        "new_active_scaffolds_at_k": len(new_cores),
        "novelty_unestablished_hit_count": sum(r.get("novel_scaffold") is None for r in hits),
        "pr_auc": _average_precision(y, scores), "pr_auc_method": "threshold-group average precision",
        "bedroc": _bedroc(y, scores, alpha), "bedroc_alpha": alpha,
        "enrichment_factor_1pct": ef,
        "procurement_success_rate": sum(r["procurement_success"] for r in procured) / len(procured) if procured else None,
        "procurement_assessed_count": len(procured),
        "both_measured_classes": 0 < positive_count < len(observed),
        "parent_identity_missing_count": len(candidates) - len(parents),
        "interpretation": "retrospective metrics conditional on measured outcomes; prospective missing outcomes remain explicit"}


def audit_splits(rows: list[dict[str, Any]], *, near_similarity: float = .85,
                 maximum_pairs: int = 2_000_000, check_temporal: bool = False) -> dict[str, Any]:
    """Audit exact/parent/scaffold/study/assay and cross-split near-analogue leakage.

    All acyclic molecules share a conservative leakage group. Pair checks stop at
    the explicit budget and report incomplete; an incomplete check never passes.
    Source IDs must identify studies/assays, not a broad database such as ChEMBL.
    """
    if not 0 < near_similarity <= 1 or not math.isfinite(near_similarity):
        raise ValueError("near similarity must be in (0, 1]")
    if isinstance(maximum_pairs, bool) or not isinstance(maximum_pairs, int) or maximum_pairs < 0:
        raise ValueError("maximum_pairs must be a nonnegative integer")
    groups = {name: defaultdict(set) for name in ("candidate", "exact", "parent", "scaffold", "study", "assay")}
    molecules = []
    incomplete = []
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    for index, row in enumerate(rows):
        split = row.get("split")
        if split not in {"train", "calibration", "validation", "test"}:
            raise ValueError("each row requires a named train/calibration/validation/test split")
        identity = chemical_identity(row.get("smiles"))
        if identity.get("parent_ambiguous"):
            incomplete.append({"row_index": index, "missing": "unambiguous_material_parent"})
        fields = {"candidate": row.get("candidate_id"), "exact": identity["exact_structure_id"],
                  "parent": identity["parent_id"], "scaffold": identity["scaffold_id"],
                  "study": row.get("study_id"), "assay": row.get("assay_id")}
        for name, value in fields.items():
            if value:
                groups[name][value].add(split)
            elif name in {"candidate", "exact", "study", "assay"}:
                incomplete.append({"row_index": index, "missing": name})
        mol = Chem.MolFromSmiles(identity["parent_smiles"]) if identity["parent_smiles"] else None
        molecules.append((split, generator.GetFingerprint(mol) if mol is not None else None))
    overlap = {name: [{"group": group, "splits": sorted(splits)} for group, splits in mapping.items()
                      if len(splits) > 1] for name, mapping in groups.items()}
    split_sizes = Counter(split for split, fp in molecules if fp is not None)
    n_valid = sum(split_sizes.values())
    total = (n_valid * n_valid - sum(size * size for size in split_sizes.values())) // 2
    checked = 0
    near_pairs = []
    for i, (split, fp) in enumerate(molecules):
        if checked >= maximum_pairs:
            break
        for j in range(i + 1, len(molecules)):
            other_split, other_fp = molecules[j]
            if split == other_split or fp is None or other_fp is None:
                continue
            if checked >= maximum_pairs:
                break
            similarity = DataStructs.TanimotoSimilarity(fp, other_fp)
            checked += 1
            if similarity >= near_similarity:
                near_pairs.append({"first_index": i, "second_index": j, "similarity": similarity})
    temporal = {"checked": check_temporal, "passed": None}
    if check_temporal:
        dates: dict[str, list[datetime]] = defaultdict(list)
        for index, row in enumerate(rows):
            try:
                value = datetime.fromisoformat(row["measured_at"].replace("Z", "+00:00"))
                # A calendar date is sufficient for ordered holdout construction.
                dates[row["split"]].append(value.date())
            except (KeyError, TypeError, ValueError):
                incomplete.append({"row_index": index, "missing": "measured_at"})
        earlier = dates["train"] + dates["calibration"] + dates["validation"]
        temporal["passed"] = bool(earlier and dates["test"] and max(earlier) < min(dates["test"]))
    complete = checked == total and not incomplete
    passed = complete and not any(overlap.values()) and not near_pairs and (not check_temporal or temporal["passed"])
    return {"passed": bool(passed), "complete": complete, "overlap": overlap,
        "near_analogue_pairs": near_pairs, "near_similarity": near_similarity,
        "pairs_checked": checked, "pairs_required": total, "missing_metadata": incomplete,
        "temporal": temporal, "acyclic_policy": "one conservative leakage group"}


def cluster_bootstrap(rows: list[dict[str, Any]], *, metric: str = "hit_rate_among_measured_top_k",
                      repeats: int = 500, seed: int = 42, k: int = 96,
                      evaluator: Callable = evaluate_ranked_candidates) -> dict[str, Any]:
    """Hierarchical target then scaffold bootstrap; macro average across targets.

    Resampled copies get synthetic IDs only within this statistical calculation.
    Identical scaffolds are resampled as units. Intervals are percentile intervals
    conditional on the observed target/scaffold population, not prospective proof.
    """
    if isinstance(repeats, bool) or not isinstance(repeats, int) or repeats < 2:
        raise ValueError("at least two bootstrap repeats are required")
    groups: dict[str, dict[str, list[dict]]] = defaultdict(lambda: defaultdict(list))
    for row in rows:
        if not row.get("target_id") or not row.get("scaffold_id"):
            raise ValueError("target_id and scaffold_id are required for clustered inference")
        groups[row["target_id"]][row["scaffold_id"]].append(row)
    targets = sorted(groups)
    rng = random.Random(seed)
    estimates = []
    base = []
    contributing_targets = []
    for target in targets:
        value = evaluator([r for block in groups[target].values() for r in block], k=k).get(metric)
        if finite(value) is not None:
            base.append(value)
            contributing_targets.append(target)
    for _ in range(repeats):
        values = []
        for target in rng.choices(contributing_targets, k=len(contributing_targets)):
            cores = sorted(groups[target])
            sample = []
            for copy_index, core in enumerate(rng.choices(cores, k=len(cores))):
                for row_index, original in enumerate(groups[target][core]):
                    row = deepcopy(original)
                    row["candidate_id"] = f"bootstrap:{copy_index}:{row_index}:{original['candidate_id']}"
                    if row.get("parent_id"):
                        row["parent_id"] = f"bootstrap:{copy_index}:{row['parent_id']}"
                    sample.append(row)
            value = evaluator(sample, k=k).get(metric)
            if finite(value) is not None:
                values.append(value)
        if values:
            estimates.append(mean(values))
    interval = np.quantile(estimates, [.025, .975]).tolist() if len(estimates) >= 2 else None
    return {"metric": metric, "estimate": mean(base) if base else None,
        "interval_95": interval, "valid_repeats": len(estimates), "requested_repeats": repeats,
        "targets": len(targets), "scaffold_clusters": sum(len(g) for g in groups.values()),
        "contributing_targets": len(contributing_targets),
        "contributing_scaffold_clusters": sum(len(groups[t]) for t in contributing_targets),
        "noncontributing_target_ids": sorted(set(targets) - set(contributing_targets)),
        "seed": seed, "method": "hierarchical target/scaffold percentile bootstrap; macro target mean",
        "inferential_support": "insufficient_clusters" if len(contributing_targets) < 2
            or any(len(groups[t]) < 2 for t in contributing_targets) else "conditional_on_observed_clusters"}
