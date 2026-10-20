"""Measured-data optional learning and independently recomputed CAS evaluation.

All optional numerical/ML imports are lazy. Outputs are assay-specific research
results, never binding probabilities, skin efficacy or automatic promotions.
"""

from __future__ import annotations

from collections import defaultdict
from copy import deepcopy
import csv
import hashlib
import importlib.metadata
import json
import math
from pathlib import Path
import pickle
import random
import tempfile

from .artifacts import ArtifactIntegrityError, canonical_json
from .algorithm_pockets import _read, _trusted_profile

CONTROL_ROLES = {"control", "positive_control", "negative_control", "counter_control"}
CONTEXT = ("target_id", "endpoint", "species", "unit", "assay_type")
MAX_ROWS = 2000


def _held(reason, candidates=()):
    return {
        "status": "held",
        "reason": str(reason),
        "candidates": deepcopy(list(candidates)),
        "candidate_population_changed": False,
        "promoted": False,
        "ranking_weight": 0,
        "scientific_validation": "NOT_ESTABLISHED",
    }


def _finite(value):
    try:
        return (
            not isinstance(value, bool)
            and isinstance(value, (float, int))
            and math.isfinite(value)
        )
    except OverflowError:
        return False


def _control(row):
    return row.get("role") in CONTROL_ROLES or row.get("opposing_direction") is True


def _artifact(sha256, fetch_artifact, *, maximum=32 * 1024 * 1024):
    if not isinstance(sha256, str) or len(sha256) != 64:
        raise ArtifactIntegrityError("A SHA256 artifact reference is required")
    raw = bytes(fetch_artifact(sha256))
    if len(raw) > maximum or hashlib.sha256(raw).hexdigest() != sha256:
        raise ArtifactIntegrityError("Evaluation artifact hash/size mismatch")
    return raw


def _json(sha256, fetch_artifact):
    def reject(value):
        raise ValueError("Nonfinite JSON evaluation input")

    return json.loads(_artifact(sha256, fetch_artifact), parse_constant=reject)


def _publish(payload, publish_artifact):
    raw = canonical_json(payload)
    receipt = publish_artifact(raw)
    if receipt.get("sha256") != hashlib.sha256(raw).hexdigest():
        raise ArtifactIntegrityError("Evaluation publication hash mismatch")
    return receipt


def _identities(rows, *, allow_controls=True):
    from .identity import chemical_identity

    if not isinstance(rows, list) or not rows or len(rows) > MAX_ROWS:
        raise ValueError("A frozen population of 1..2000 rows is required")
    out, ids, parents = [], set(), set()
    for original in rows:
        if not isinstance(original, dict):
            raise ValueError("Every population row must be a record")
        row = deepcopy(original)
        if (
            not isinstance(row.get("candidate_id"), str)
            or not row["candidate_id"]
            or row["candidate_id"] in ids
        ):
            raise ValueError("Candidate IDs must be unique nonempty strings")
        ids.add(row["candidate_id"])
        identity = chemical_identity(row.get("smiles"))
        if not identity["parent_id"] or identity.get("parent_ambiguous"):
            raise ValueError(
                "Valid exact chemical structures are required for evaluation"
            )
        if row.get("exact_structure_id") not in (None, identity["exact_structure_id"]):
            raise ValueError("Claimed exact structure differs from source SMILES")
        if row.get("parent_id") not in (None, identity["parent_id"]):
            raise ValueError("Claimed parent differs from versioned chemical identity")
        row.update(
            parent_id=identity["parent_id"],
            exact_structure_id=identity["exact_structure_id"],
            scaffold_id=identity["scaffold_id"],
            scaffold_kind=identity["scaffold_kind"],
        )
        if not _control(row) or not allow_controls:
            if row["parent_id"] in parents:
                raise ValueError(
                    "Parent duplicates cannot inflate candidate/training counts"
                )
            parents.add(row["parent_id"])
        out.append(row)
    control_parents = {row["parent_id"] for row in out if _control(row)}
    if any(row["parent_id"] in control_parents and not _control(row) for row in out):
        raise ValueError("A known control parent cannot also count as a new candidate")
    return out


def _context(rows):
    values = {tuple(row.get(key) for key in CONTEXT) for row in rows}
    if len(values) != 1 or any(
        not isinstance(value, str) or not value for value in next(iter(values))
    ):
        raise ValueError(
            "One explicit target/species/endpoint/unit/assay-type context is required"
        )
    return dict(zip(CONTEXT, next(iter(values))))


def measured_partitions(records):
    """Return measured partitions and a complete recomputed overlap audit.

    Missing observations and controls remain in `retained`; they never become
    negative outcomes. No automatic splitting or test-guided model selection.
    """
    from .evaluation import audit_splits

    rows = _identities(records)
    context = _context(rows)
    retained, measured = [], []
    for row in rows:
        if (
            _control(row)
            or row.get("status") != "measured"
            or row.get("relation") != "="
            or not _finite(row.get("value"))
        ):
            retained.append(row)
            continue
        row["split"] = row.get("partition", row.get("split"))
        row["study_id"] = row.get("study_id", row.get("source_group"))
        if row["split"] not in {"train", "validation", "calibration", "test"}:
            raise ValueError(
                "Measured rows require explicit train/validation/calibration/test partitions"
            )
        measured.append(row)
    if not measured:
        raise ValueError("No measured uncensored outcomes")
    audit = audit_splits(measured, maximum_pairs=2_000_000)
    if not audit["passed"]:
        raise ValueError(
            "Scaffold/source/assay/parent or near-analogue separation failed: "
            + json.dumps(
                {k: v for k, v in audit.items() if k != "near_analogue_pairs"},
                sort_keys=True,
            )
        )
    groups = {
        name: [row for row in measured if row["split"] == name]
        for name in ("train", "validation", "calibration", "test")
    }
    return {
        "partitions": groups,
        "context": context,
        "audit": audit,
        "retained": retained,
        "all_rows": rows,
    }


def evaluate_comparison(request, *, fetch_artifact, publish_artifact):
    """Recompute evaluation; supplied scores, flags and confidence intervals are rejected.

    request={evaluation_manifest_sha256}. Manifest schema=skinpipe-evaluation-v1
    references snapshot_sha256 (rows), model_binding_sha256,
    baseline_model_binding_sha256, policy_sha256, predictions_sha256 and
    baseline_predictions_sha256. Both model bindings contain
    model_sha256 (actual model bytes), training_rows_sha256 and context. Each
    prediction bundle binds snapshot_sha256, model_sha256, budget_seconds,
    compute_seconds and rows[{candidate_id,exact_structure_id,score,probability?}].
    A frozen policy supplies k, bootstrap_repeats,seed and reference_smiles.
    No result here changes the active model or grants automatic promotion.
    """
    from .evaluation import audit_splits, evaluate_ranked_candidates
    from .identity import chemical_identity

    try:
        if set(request) != {"evaluation_manifest_sha256"}:
            raise ValueError(
                "Evaluation accepts one immutable manifest, not promotion assertions"
            )
        manifest = _json(request["evaluation_manifest_sha256"], fetch_artifact)
        if (
            not isinstance(manifest, dict)
            or manifest.get("schema") != "skinpipe-evaluation-v1"
        ):
            raise ValueError("Unsupported evaluation manifest")
        rows = _identities(_json(manifest["snapshot_sha256"], fetch_artifact))
        for row in rows:
            if row.get("label_measured") is True and (
                row.get("status") != "measured"
                or row.get("relation") != "="
                or type(row.get("label")) is not int
                or row["label"] not in (0, 1)
            ):
                raise ValueError(
                    "Binary evaluation labels require explicit measured uncensored outcomes"
                )
        binding = _json(manifest["model_binding_sha256"], fetch_artifact)
        _artifact(binding["model_sha256"], fetch_artifact, maximum=1024 * 1024 * 1024)
        training = _identities(_json(binding["training_rows_sha256"], fetch_artifact))
        baseline_binding = _json(
            manifest["baseline_model_binding_sha256"], fetch_artifact
        )
        _artifact(
            baseline_binding["model_sha256"], fetch_artifact, maximum=1024 * 1024 * 1024
        )
        baseline_training = _identities(
            _json(baseline_binding["training_rows_sha256"], fetch_artifact)
        )
        context = _context(rows)
        if (
            binding.get("context") != context
            or _context(training) != context
            or baseline_binding.get("context") != context
            or _context(baseline_training) != context
        ):
            raise ValueError("Model, evaluation and training contexts differ")
        policy = _json(manifest["policy_sha256"], fetch_artifact)
        if not isinstance(policy, dict):
            raise ValueError("Evaluation policy must be a record")
        k, repeats, seed = (
            policy.get("k", 96),
            policy.get("bootstrap_repeats", 500),
            policy.get("seed", 42),
        )
        if (
            any(
                isinstance(value, bool) or not isinstance(value, int)
                for value in (k, repeats, seed)
            )
            or not 1 <= k <= 1000
            or not 20 <= repeats <= 2000
        ):
            raise ValueError("Invalid bounded evaluation policy")
        audits = {}
        for name, used_rows in (("model", training), ("baseline", baseline_training)):
            independent = []
            for records, split in ((used_rows, "train"), (rows, "test")):
                for row in records:
                    independent.append(
                        {
                            **row,
                            "split": split,
                            "study_id": row.get("study_id", row.get("source_group")),
                        }
                    )
            audits[name] = audit_splits(
                independent,
                maximum_pairs=2_000_000,
                check_temporal=policy.get("check_temporal", False),
            )
            if not audits[name]["passed"]:
                raise ValueError(name + " independent evaluation overlap audit failed")
        overlap = {"passed": all(a["passed"] for a in audits.values()),
                   "complete": all(a["complete"] for a in audits.values()),
                   "temporal_checked": all(a["temporal"]["checked"] for a in audits.values()),
                   "scope": "structure/source separation; temporal separation only when explicitly requested",
                   **audits}
        references = policy.get("reference_smiles", [])
        ref_scaffolds = set()
        if not isinstance(references, list):
            raise ValueError("Frozen reference structures must be a list")
        for smiles in references:
            identity = chemical_identity(smiles)
            if not identity["exact_structure_id"]:
                raise ValueError("Invalid frozen reference structure")
            if identity["scaffold_kind"] == "murcko":
                ref_scaffolds.add(identity["scaffold_id"])
        bundles, ranked = [], []
        for key in ("predictions_sha256", "baseline_predictions_sha256"):
            bundle = _json(manifest[key], fetch_artifact)
            if not isinstance(bundle, dict):
                raise ValueError("Prediction bundle must be a record")
            if bundle.get("snapshot_sha256") != manifest["snapshot_sha256"]:
                raise ValueError(
                    "Comparators must use the same exact frozen population"
                )
            _artifact(
                bundle["model_sha256"], fetch_artifact, maximum=1024 * 1024 * 1024
            )
            expected_binding = (
                binding if key == "predictions_sha256" else baseline_binding
            )
            if bundle["model_sha256"] != expected_binding["model_sha256"]:
                raise ValueError("Evaluation predictions refer to a different model")
            budget, compute = (
                bundle.get("budget_seconds"),
                bundle.get("compute_seconds"),
            )
            if (
                not _finite(budget)
                or not _finite(compute)
                or not 0 <= compute <= budget
                or budget <= 0
            ):
                raise ValueError(
                    "Measured compute and a positive frozen budget are required"
                )
            supplied = bundle["rows"]
            lookup = {item["candidate_id"]: item for item in supplied}
            if len(lookup) != len(supplied) or set(lookup) != {
                row["candidate_id"] for row in rows
            }:
                raise ValueError(
                    "Every preserved candidate/control needs one comparator prediction"
                )
            assembled = []
            for row in rows:
                prediction = lookup[row["candidate_id"]]
                if prediction.get("exact_structure_id") != row[
                    "exact_structure_id"
                ] or not _finite(prediction.get("score")):
                    raise ValueError("Prediction structure or score mismatch")
                probability = prediction.get("probability")
                if probability is not None and (
                    not _finite(probability) or not 0 <= probability <= 1
                ):
                    raise ValueError("Assay outcome probability is outside [0,1]")
                assembled.append(
                    {
                        **row,
                        "score": prediction["score"],
                        "probability": probability,
                        "novel_scaffold": (
                            row["scaffold_kind"] == "murcko"
                            and row["scaffold_id"] not in ref_scaffolds
                        )
                        if references
                        else None,
                    }
                )
            bundles.append(bundle)
            ranked.append(assembled)
        if bundles[0]["budget_seconds"] != bundles[1]["budget_seconds"]:
            raise ValueError("Comparator budgets differ")

        def metrics(items):
            values = evaluate_ranked_candidates(items, k=k)
            observed = [
                r
                for r in items
                if not _control(r)
                and r.get("label_measured") is True
                and type(r.get("label")) is int
                and r["label"] in (0, 1)
            ]
            values["brier"] = (
                sum((r["probability"] - r["label"]) ** 2 for r in observed)
                / len(observed)
                if observed and all(r.get("probability") is not None for r in observed)
                else None
            )
            return values

        measurements = [metrics(items) for items in ranked]
        # Resample connected scaffold/source clusters as paired observations.
        candidates = [r for r in ranked[0] if not _control(r)]
        parent = list(range(len(candidates)))

        def find(index):
            while parent[index] != index:
                parent[index] = parent[parent[index]]
                index = parent[index]
            return index

        owners = {}
        for index, row in enumerate(candidates):
            for key in (
                ("scaffold", row["scaffold_id"]),
                ("source", row.get("study_id", row.get("source_group"))),
                ("assay", row.get("assay_id")),
            ):
                if key in owners:
                    a, b = find(index), find(owners[key])
                    parent[max(a, b)] = min(a, b)
                owners[key] = index
        clusters = defaultdict(list)
        for index, row in enumerate(candidates):
            clusters[find(index)].append(row["candidate_id"])
        lookup = [{r["candidate_id"]: r for r in items} for items in ranked]
        rng, samples = (
            random.Random(seed),
            {
                name: []
                for name in ("pr_auc", "brier", "hit_at_k", "new_active_scaffolds_at_k")
            },
        )
        keys = sorted(clusters)
        # Two clusters have only three distinct count configurations under a
        # size-two bootstrap. Do not turn that degenerate interval into a
        # minimum-improvement claim. Three is a reporting floor, not power.
        min_clusters = 3
        for _ in range(repeats if len(keys) >= min_clusters else 0):
            paired = [[], []]
            for number, key in enumerate(rng.choices(keys, k=len(keys))):
                for cid in clusters[key]:
                    for arm in range(2):
                        row = deepcopy(lookup[arm][cid])
                        row["candidate_id"] = f"resample:{number}:{cid}"
                        row["parent_id"] = f"resample:{number}:" + row["parent_id"]
                        paired[arm].append(row)
            first, second = (metrics(items) for items in paired)
            for name in samples:
                if _finite(first.get(name)) and _finite(second.get(name)):
                    samples[name].append(first[name] - second[name])
        import numpy as np

        intervals = {
            name: np.quantile(values, [0.025, 0.975]).tolist()
            if len(values) >= 20
            else None
            for name, values in samples.items()
        }
        pr, brier = intervals["pr_auc"], intervals["brier"]
        receipt = {
            "schema": "skinpipe-recomputed-evaluation-v1",
            "status": "evaluated",
            "request": deepcopy(request),
            "model_sha256": binding["model_sha256"],
            "baseline_model_sha256": baseline_binding["model_sha256"],
            "context": context,
            "snapshot_sha256": manifest["snapshot_sha256"],
            "training_rows_sha256": binding["training_rows_sha256"],
            "baseline_training_rows_sha256": baseline_binding["training_rows_sha256"],
            "policy_sha256": manifest["policy_sha256"],
            "overlap_audit": overlap,
            "budget_seconds": bundles[0]["budget_seconds"],
            "compute_seconds": [b["compute_seconds"] for b in bundles],
            "model_metrics": measurements[0],
            "baseline_metrics": measurements[1],
            "paired_difference_ci_95": intervals,
            "confidence_interval_unit": "scaffold_and_source",
            "independent_clusters": len(keys),
            "minimum_independent_clusters": min_clusters,
            "inference_status": "UNDERPOWERED" if len(keys) < min_clusters else "CONDITIONAL_RETROSPECTIVE",
            "cluster_count_is_power_calculation": False,
            "minimum_improvement_criteria_met": bool(
                len(keys) >= min_clusters and pr and brier and pr[0] > 0 and brier[1] < 0
            ),
            "promoted": False,
            "promotion_status": "requires_separate_versioned_decision",
            "prespecification_status": "not_verified_by_content_hashes; check frozen policy and outcome event ordering",
            "scientific_validation": "conditional_retrospective_evaluation_only",
        }
        receipt["artifact"] = _publish(receipt, publish_artifact)
        return receipt
    except (OSError, ValueError, KeyError, TypeError, ImportError) as error:
        return _held(error)


def _features(rows):
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Lipinski
    import numpy as np

    values = []
    for row in rows:
        mol = Chem.MolFromSmiles(row["smiles"])
        values.append(
            [
                Descriptors.MolWt(mol),
                Descriptors.MolLogP(mol),
                Lipinski.NumHDonors(mol),
                Lipinski.NumHAcceptors(mol),
                Descriptors.TPSA(mol),
                Lipinski.NumRotatableBonds(mol),
            ]
        )
    return np.asarray(values, dtype=float)


def mapie_intervals(records, candidates, *, confidence=0.9, seed=42):
    """Actual split-conformal RF regression with isolated calibration/test rows."""
    try:
        from mapie.regression import SplitConformalRegressor
        from sklearn.ensemble import RandomForestRegressor
        import numpy as np

        if not _finite(confidence) or not 0.5 < confidence < 1:
            raise ValueError("Confidence must be strictly between .5 and1")
        data = measured_partitions(records)
        train, calibration, test = (
            data["partitions"][name] for name in ("train", "calibration", "test")
        )
        if (
            len(train) < 5
            or len(calibration) < max(10, math.ceil(1 / (1 - confidence)))
            or not test
        ):
            raise ValueError(
                "Insufficient separate training/calibration/test observations"
            )
        queries = _identities(candidates)
        if _context(queries) != data["context"]:
            raise ValueError(
                "Candidate prediction context differs from measured context"
            )
        estimator = RandomForestRegressor(
            n_estimators=128, min_samples_leaf=2, n_jobs=1, random_state=seed
        )
        model = SplitConformalRegressor(
            estimator, confidence_level=confidence, prefit=False
        )
        model.fit(_features(train), [r["value"] for r in train])
        model.conformalize(_features(calibration), [r["value"] for r in calibration])
        points, intervals = model.predict_interval(_features(queries))
        test_points, test_intervals = model.predict_interval(_features(test))
        intervals = np.asarray(intervals)[:, :, 0]
        checked = np.asarray(test_intervals)[:, :, 0]
        if (
            not np.isfinite(points).all()
            or not np.isfinite(intervals).all()
            or np.any(intervals[:, 0] > intervals[:, 1])
        ):
            raise ValueError("Invalid conformal output")
        result = []
        train_features = _features(train)
        for row, features, point, interval in zip(
            queries, _features(queries), points, intervals
        ):
            result.append(
                {
                    **row,
                    "prediction": float(point),
                    "interval": interval.tolist(),
                    "confidence": confidence,
                    "prediction_unit": data["context"]["unit"],
                    "kind": "assay_regression_prediction",
                    "applicability": "within_descriptor_range"
                    if np.all(features >= train_features.min(0))
                    and np.all(features <= train_features.max(0))
                    else "outside_descriptor_range",
                    "eligible_for_ranking": False,
                    "scientific_validation": "NOT_ESTABLISHED",
                }
            )
        y = np.asarray([r["value"] for r in test])
        return {
            "status": "complete",
            "candidates": result,
            "retained_observations": data["retained"],
            "test_coverage": float(
                np.mean((checked[:, 0] <= y) & (y <= checked[:, 1]))
            ),
            "test_count": len(test),
            "context": data["context"],
            "overlap_audit": data["audit"],
            "versions": {
                "mapie": importlib.metadata.version("mapie"),
                "sklearn": importlib.metadata.version("scikit-learn"),
            },
            "coverage_scope": "marginal exchangeability assumption; not guaranteed under scaffold or assay shift",
            "promoted": False,
            "ranking_weight": 0,
        }
    except (
        ImportError,
        OSError,
        ValueError,
        KeyError,
        TypeError,
        RuntimeError,
    ) as error:
        return _held(error, candidates)


def recommend_qlognehvi(records, candidates, *, objectives, q=12, seed=42):
    """Discrete qLogNEHVI over purchasable, unmeasured, same-target unique parents.

    objectives=[{name,unit,higher_is_better}]; each measured training row provides
    outcomes[name]={value,status:'measured',unit}. Only partition=train is used.
    Acquisition values are log hypervolume utility, never activity probability.
    """
    try:
        from .evaluation import audit_splits
        import numpy as np
        import torch
        from botorch.models import SingleTaskGP
        from botorch.models.transforms.outcome import Standardize
        from botorch.fit import fit_gpytorch_mll
        from gpytorch.mlls import ExactMarginalLogLikelihood
        from botorch.acquisition.multi_objective.logei import (
            qLogNoisyExpectedHypervolumeImprovement,
        )
        from botorch.optim import optimize_acqf_discrete
        from botorch.sampling.normal import SobolQMCNormalSampler

        if isinstance(q, bool) or not isinstance(q, int) or not 1 <= q <= 96:
            raise ValueError("Recommendation q must be1..96")
        if not isinstance(objectives, list) or not 2 <= len(objectives) <= 4:
            raise ValueError("Two to four measured objectives are required")
        names = [o["name"] for o in objectives]
        if len(set(names)) != len(names) or any(
            not o.get("unit") or not isinstance(o.get("higher_is_better"), bool)
            for o in objectives
        ):
            raise ValueError(
                "Distinct explicit objective units and directions required"
            )
        all_training = _identities(records)
        queries = _identities(candidates)
        targets = {row.get("target_id") for row in all_training + queries}
        if len(targets) != 1 or None in targets:
            raise ValueError("Acquisition must stay within one explicit target")
        context = _context(all_training + queries)
        audit_rows = [
            {
                **row,
                "split": row.get("partition", row.get("split")),
                "study_id": row.get("study_id", row.get("source_group")),
            }
            for row in all_training
        ]
        overlap = audit_splits(audit_rows, maximum_pairs=2_000_000)
        if not overlap["passed"]:
            raise ValueError("Multi-objective holdout separation failed")
        train, outcomes = [], []
        for row in all_training:
            if _control(row) or row.get("partition", row.get("split")) != "train":
                continue
            values = [row.get("outcomes", {}).get(name, {}) for name in names]
            if any(
                value.get("status") != "measured"
                or value.get("relation", "=") != "="
                or not _finite(value.get("value"))
                or value.get("unit") != objective["unit"]
                for value, objective in zip(values, objectives)
            ):
                continue
            if not row.get("source_group") or not row.get("assay_id"):
                raise ValueError(
                    "Measured multi-objective source/assay identity is required"
                )
            train.append(row)
            outcomes.append(
                [
                    value["value"] * (1 if objective["higher_is_better"] else -1)
                    for value, objective in zip(values, objectives)
                ]
            )
        if len(train) < 4:
            raise ValueError(
                "At least four measured multi-objective training parents are required"
            )
        excluded = {row["parent_id"] for row in all_training}
        eligible = [
            row
            for row in queries
            if not _control(row)
            and row.get("purchasable") is True
            and row["parent_id"] not in excluded
            and row.get("status") != "measured"
        ]
        if not eligible:
            raise ValueError(
                "No unmeasured purchasable noncontrol parents are eligible"
            )
        eligible.sort(key=lambda row: row["candidate_id"])
        train_x, candidates_x = _features(train), _features(eligible)
        low, scale = train_x.min(0), np.maximum(train_x.max(0) - train_x.min(0), 1.0)
        train_x, candidates_x = (train_x - low) / scale, (candidates_x - low) / scale
        # Descriptor collisions remain visible, but cannot masquerade as distinct
        # discrete acquisition points. Exact parent exclusion is independent.
        distinct, seen = [], set()
        for index, vector in enumerate(candidates_x):
            key = tuple(vector.tolist())
            if key not in seen:
                distinct.append(index)
                seen.add(key)
        eligible = [eligible[index] for index in distinct]
        candidates_x = candidates_x[distinct]
        with torch.random.fork_rng():
            torch.manual_seed(seed)
            tx, ty, cx = [
                torch.tensor(value, dtype=torch.double)
                for value in (train_x, outcomes, candidates_x)
            ]
            model = SingleTaskGP(
                tx, ty, outcome_transform=Standardize(m=len(objectives))
            )
            fit_gpytorch_mll(
                ExactMarginalLogLikelihood(model.likelihood, model),
                optimizer_kwargs={"options": {"maxiter": 40}},
            )
            span = (ty.max(0).values - ty.min(0).values).clamp_min(1.0)
            ref = ty.min(0).values - 0.1 * span
            acquisition = qLogNoisyExpectedHypervolumeImprovement(
                model=model,
                ref_point=ref.tolist(),
                X_baseline=tx,
                sampler=SobolQMCNormalSampler(torch.Size([64]), seed=seed),
                prune_baseline=True,
            )
            chosen, values = optimize_acqf_discrete(
                acquisition,
                q=min(q, len(eligible)),
                choices=cx,
                max_batch_size=64,
                unique=True,
            )
            indices = [
                int(torch.argmin(((cx - vector) ** 2).sum(-1))) for vector in chosen
            ]
        ids = [eligible[index]["candidate_id"] for index in indices]
        return {
            "status": "complete",
            "recommended_ids": ids,
            "candidate_count": len(ids),
            "candidates": deepcopy(candidates),
            "deferred_ids": [
                row["candidate_id"] for row in queries if row["candidate_id"] not in ids
            ],
            "objective_specification": deepcopy(objectives),
            "context": context,
            "overlap_audit": overlap,
            "measured_training_count": len(train),
            "excluded_unmeasured_training_count": len(all_training) - len(train),
            "acquisition": "qLogNEHVI",
            "acquisition_values": values.detach().cpu().reshape(-1).tolist(),
            "value_scope": "log expected hypervolume utility in specified objective directions",
            "versions": {
                "botorch": importlib.metadata.version("botorch"),
                "torch": torch.__version__,
            },
            "promoted": False,
            "ranking_weight": 0,
            "scientific_validation": "NOT_ESTABLISHED",
        }
    except (
        ImportError,
        OSError,
        ValueError,
        KeyError,
        TypeError,
        RuntimeError,
    ) as error:
        return _held(error, candidates)


def train_chemprop(request, profile, *, fetch_artifact, publish_artifact, run_owned):
    """Measured Chemprop CLI training, optional pinned local CheMeleon foundation.

    The owned executor enforces lease/cancellation and a <=1h budget. Dataset
    JSON rows require explicit four-way splits; calibration rows never train or
    select checkpoints. Existing `.pt` foundation paths are administrator assets;
    the automatic-download `chemeleon` name is deliberately not used.
    """
    try:
        if set(request) - {
            "dataset_sha256",
            "epochs",
            "seed",
            "task_type",
            "timeout_seconds",
            "foundation",
        }:
            raise ValueError(
                "Training accepts frozen dataset identity and bounded parameters only"
            )
        selected = _trusted_profile(profile, "chemprop")
        records = _json(request["dataset_sha256"], fetch_artifact)
        data = measured_partitions(records)
        parts = data["partitions"]
        if not all(
            parts[name] for name in ("train", "validation", "calibration", "test")
        ):
            raise ValueError(
                "Separate measured train/validation/calibration/test partitions are required"
            )
        task = request.get("task_type", "regression")
        if task not in {"regression", "classification"}:
            raise ValueError(
                "Only explicit measured regression/binary classification is supported"
            )
        if task == "classification" and any(
            row["value"] not in (0, 1) for rows in parts.values() for row in rows
        ):
            raise ValueError("Classification requires measured binary labels")
        epochs, seed, timeout = (
            request.get("epochs", 20),
            request.get("seed", 42),
            request.get("timeout_seconds", 3600),
        )
        if (
            any(isinstance(v, bool) or not isinstance(v, int) for v in (epochs, seed))
            or not 1 <= epochs <= 100
            or not _finite(timeout)
            or not 0 < timeout <= 3600
        ):
            raise ValueError("Training epochs/seed/one-hour budget invalid")
        with tempfile.TemporaryDirectory(prefix="skinpipe-chemprop-") as temporary:
            work = Path(temporary)
            paths = []
            for name in ("train", "validation", "test"):
                path = work / (name + ".csv")
                with path.open("w", newline="") as stream:
                    writer = csv.writer(stream)
                    writer.writerow(["smiles", "outcome"])
                    writer.writerows(
                        [row["smiles"], row["value"]] for row in parts[name]
                    )
                paths.append(str(path))
            output = work / "model"
            argv = [
                selected["executable"],
                "train",
                "--data-path",
                *paths,
                "--smiles-columns",
                "smiles",
                "--target-columns",
                "outcome",
                "--task-type",
                task,
                "--output-dir",
                str(output),
                "--epochs",
                str(epochs),
                "--warmup-epochs",
                "0",
                "--data-seed",
                str(seed),
                "--pytorch-seed",
                str(seed),
                "--num-workers",
                "0",
                "--accelerator",
                "cpu",
                "--devices",
                "1",
                "--message-hidden-dim",
                "64",
                "--ffn-hidden-dim",
                "64",
                "--depth",
                "2",
                "--batch-size",
                "32",
            ]
            if request.get("foundation") is not None:
                if request["foundation"] != "chemeleon":
                    raise ValueError("Unsupported foundation model")
                foundation = selected.get("foundation_path")
                if foundation not in {
                    item["path"] for item in selected.get("model_files", [])
                }:
                    raise ValueError(
                        "CheMeleon requires a pinned local checkpoint; no automatic downloads"
                    )
                # The official CheMeleon file contains message-passing weights,
                # whereas Chemprop's local --from-foundation expects an MPNN.
                # Convert only the pinned weights-only source; never ask the CLI
                # to auto-download through its named foundation option.
                import torch
                from chemprop import models, nn
                from chemprop.models.utils import save_model

                message = torch.load(foundation, weights_only=True, map_location="cpu")
                if not isinstance(message, dict) or set(message) != {
                    "hyper_parameters",
                    "state_dict",
                }:
                    raise ValueError(
                        "Expected a native CheMeleon message-passing checkpoint"
                    )
                mp = nn.BondMessagePassing(**message["hyper_parameters"])
                mp.load_state_dict(message["state_dict"], strict=True)
                converted = work / "pinned-foundation.pt"
                with torch.random.fork_rng():
                    torch.manual_seed(seed)
                    model = models.MPNN(
                        mp,
                        nn.MeanAggregation(),
                        nn.RegressionFFN(input_dim=mp.output_dim),
                    )
                    save_model(converted, model)
                argv += ["--from-foundation", str(converted)]
            execution = run_owned(argv, cwd=work, timeout_seconds=float(timeout))
            if (
                execution.get("status") != "complete"
                or execution.get("returncode") != 0
            ):
                return {
                    **_held("Owned training process did not complete", records),
                    "execution": execution,
                }
            _trusted_profile(profile, "chemprop")
            files = sorted(output.rglob("*.pt"))
            if not files:
                raise ValueError("Training exited without a model artifact")
            models = []
            for path in files[:8]:
                raw = _read(path, work, maximum=1024 * 1024 * 1024)
                receipt = publish_artifact(raw)
                if receipt.get("sha256") != hashlib.sha256(raw).hexdigest():
                    raise ArtifactIntegrityError("Trained model publication mismatch")
                models.append({"name": path.relative_to(output).as_posix(), **receipt})
            binding = {
                "schema": "skinpipe-model-evaluation-binding-v1",
                "model_sha256": models[0]["sha256"],
                "training_rows_sha256": request["dataset_sha256"],
                "context": data["context"],
                "profile_sha256": hashlib.sha256(canonical_json(selected)).hexdigest(),
                "models": models,
                "calibration_status": "held_independent_calibration_required",
                "promoted": False,
            }
            binding["artifact"] = _publish(binding, publish_artifact)
            return {
                "status": "complete",
                "binding": binding,
                "execution": deepcopy(execution),
                "argv": argv,
                "split_audit": data["audit"],
                "retained_observations": data["retained"],
                "candidate_population_changed": False,
                "promoted": False,
                "ranking_weight": 0,
                "scientific_validation": "NOT_ESTABLISHED",
            }
    except (
        OSError,
        ValueError,
        KeyError,
        TypeError,
        ImportError,
        RuntimeError,
        pickle.UnpicklingError,
        EOFError,
    ) as error:
        return _held(error)
