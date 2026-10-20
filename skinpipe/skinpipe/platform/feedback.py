"""Measured feedback and conservatively admitted fingerprint baselines.

The observational unit is one distinct parent compound in one fully specified
target/endpoint/assay context. Missing, failed and censored values are never
converted to inactive labels. Split or sample-size failures hold model scoring,
not candidate exploration or experimental result registration.
"""
from __future__ import annotations

from collections import Counter
import hashlib
import importlib.metadata
import math
from typing import Any

from .experiments import _parent_key, canonical_json

CONTEXT_FIELDS = ("target_id", "endpoint", "species", "assay_id", "unit")
FINGERPRINT_BITS = 2048
MAX_TRAINING_ROWS = 2000


class FeedbackError(ValueError):
    pass


def measured_feedback(observations: list[dict]) -> dict:
    """Expose original measurements without averaging incompatible contexts."""
    groups: dict[str, dict] = {}
    unmeasured = []
    for observation in observations:
        row = observation.get("normalized", observation)
        if row.get("status") != "measured" or row.get("value") is None:
            unmeasured.append(observation)
            continue
        context = {k: row.get(k) for k in (
            "target_id", "candidate_id", "endpoint", "unit", "assay_id", "conditions",
            "biological_context", "concentration", "concentration_unit",
        )}
        key = hashlib.sha256(canonical_json(context).encode()).hexdigest()
        group = groups.setdefault(key, {"context": context, "observations": [], "measured_count": 0,
                                        "censored_count": 0, "interpretation": "measured_assay_evidence"})
        group["observations"].append(observation)
        group["measured_count"] += 1
        group["censored_count"] += row.get("relation") != "="
    return {"kind": "measured_feedback", "groups": list(groups.values()),
            "unmeasured": unmeasured, "model_update": "not_automatic",
            "limitations": "Assay observations do not establish target-specific mechanism or clinical efficacy."}


def _fingerprint(smiles: str):
    from rdkit import Chem
    from rdkit.Chem import rdFingerprintGenerator, Descriptors, Lipinski
    from rdkit.Chem.Scaffolds import MurckoScaffold
    import numpy as np

    if not isinstance(smiles, str) or not 1 <= len(smiles) <= 4096:
        raise FeedbackError("Training structure exceeds supported input bounds")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or not 1 <= mol.GetNumAtoms() <= 512:
        raise FeedbackError("Invalid training structure")
    fp = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=FINGERPRINT_BITS).GetFingerprint(mol)
    scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=True) or "acyclic"
    descriptors = np.asarray([Descriptors.MolWt(mol), Descriptors.MolLogP(mol),
                              Lipinski.NumHDonors(mol), Lipinski.NumHAcceptors(mol)], dtype=float)
    return np.asarray(fp, dtype=float), fp, scaffold, descriptors


def leakage_groups(records: list[dict], *, near_similarity: float = 0.85) -> tuple[list[str], dict]:
    """Union scaffold, study/source, parent and near-neighbour relationships."""
    from rdkit import DataStructs

    if not 0 < near_similarity <= 1 or len(records) > MAX_TRAINING_ROWS:
        raise FeedbackError("Near-neighbour audit requires a bounded dataset and threshold in (0,1]")
    parents = list(range(len(records)))

    def find(index):
        while parents[index] != index:
            parents[index] = parents[parents[index]]
            index = parents[index]
        return index

    def union(a, b):
        a, b = find(a), find(b)
        parents[max(a, b)] = min(a, b)

    features, fingerprints, scaffolds, descriptors = [], [], [], []
    owners: dict[tuple, int] = {}
    for index, row in enumerate(records):
        if not row.get("source_group") or not row.get("parent_id"):
            raise FeedbackError("Explicit source_group and parent_id are required for leakage control")
        feature, fp, scaffold, descriptor = _fingerprint(row["smiles"])
        features.append(feature)
        fingerprints.append(fp)
        scaffolds.append(scaffold)
        descriptors.append(descriptor)
        for key in (("scaffold", scaffold), ("source", row["source_group"]), ("parent", row["parent_id"])):
            if key in owners:
                union(index, owners[key])
            else:
                owners[key] = index
    pairs = 0
    for index, fp in enumerate(fingerprints):
        for other, similarity in enumerate(DataStructs.BulkTanimotoSimilarity(fp, fingerprints[index + 1:]), index + 1):
            pairs += 1
            if similarity >= near_similarity:
                union(index, other)
    groups = [str(find(index)) for index in range(len(records))]
    return groups, {"groups": len(set(groups)), "near_similarity": near_similarity,
                    "pairs_checked": pairs, "near_analogue_screen_complete": True,
                    "method": "union of exact parent, Murcko scaffold, source group and Morgan neighbours",
                    "acyclic_policy": "one conservative group", "features": features,
                    "fingerprints": fingerprints, "descriptors": descriptors}


def _partitions(records: list[dict], groups: list[str], y, seed: int) -> dict:
    import numpy as np
    from sklearn.model_selection import GroupShuffleSplit

    labels = [row.get("partition") for row in records]
    if any(labels):
        if not all(label in {"train", "calibration", "test"} for label in labels):
            raise FeedbackError("Explicit partitions must cover every record with train/calibration/test")
        assigned: dict[str, str] = {}
        for group, label in zip(groups, labels):
            if group in assigned and assigned[group] != label:
                raise FeedbackError("A scaffold, source or near-neighbour group crosses explicit partitions")
            assigned[group] = label
        return {label: np.asarray([i for i, value in enumerate(labels) if value == label], dtype=int)
                for label in ("train", "calibration", "test")}
    if len(set(groups)) < 5:
        raise FeedbackError("At least five independent groups are required for partitioning")
    indices = np.arange(len(records))
    groups = np.asarray(groups)
    for offset in range(100):
        dev, test = next(GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=seed + offset).split(indices, y, groups))
        train_relative, cal_relative = next(GroupShuffleSplit(n_splits=1, test_size=0.25, random_state=seed + 1000 + offset).split(dev, y[dev], groups[dev]))
        parts = {"train": dev[train_relative], "calibration": dev[cal_relative], "test": test}
        if all(len(set(y[part])) == 2 for part in parts.values()):
            return parts
    raise FeedbackError("Cannot form outcome-supported independent partitions")


def _serialize(estimator, family: str) -> dict:
    if family == "logistic":
        return {"family": family, "coefficients": estimator.coef_[0].tolist(), "intercept": float(estimator.intercept_[0])}
    trees = []
    for fitted in estimator.estimators_:
        tree = fitted.tree_
        values = tree.value[:, 0, :]
        sums = values.sum(axis=1)
        trees.append({"left": tree.children_left.tolist(), "right": tree.children_right.tolist(),
                      "feature": tree.feature.tolist(), "threshold": tree.threshold.tolist(),
                      "value": (values[:, 1] / sums).tolist()})
    return {"family": family, "trees": trees}


def _raw_predict(estimator: dict, x):
    import numpy as np

    if estimator["family"] == "logistic":
        coefficients = np.asarray(estimator["coefficients"])
        if coefficients.shape != (x.shape[1],):
            raise FeedbackError("Model fingerprint dimensions differ")
        logits = x @ coefficients + estimator["intercept"]
        return 1 / (1 + np.exp(-np.clip(logits, -700, 700)))
    trees = estimator.get("trees", [])
    if estimator["family"] != "random_forest" or not 1 <= len(trees) <= 512:
        raise FeedbackError("Unsupported model family or tree count")
    prediction = np.zeros(len(x))
    x = x.astype(np.float32)
    for tree in trees:
        n = len(tree["left"])
        if not 1 <= n <= 100000 or any(len(tree[key]) != n for key in ("right", "feature", "threshold", "value")):
            raise FeedbackError("Invalid serialized tree dimensions")
        for i, vector in enumerate(x):
            node = 0
            for _ in range(n + 1):
                if not 0 <= node < n:
                    raise FeedbackError("Invalid serialized tree index")
                if tree["left"][node] == -1:
                    prediction[i] += tree["value"][node]
                    break
                feature = tree["feature"][node]
                if not 0 <= feature < x.shape[1]:
                    raise FeedbackError("Invalid serialized fingerprint feature")
                node = tree["left"][node] if vector[feature] <= tree["threshold"][node] else tree["right"][node]
            else:
                raise FeedbackError("Cyclic serialized tree")
    return prediction / len(trees)


def train_baseline(records: list[dict], *, family: str = "logistic", seed: int = 42) -> dict:
    """Fit a prespecified family; never select a model using its held test set.

    Labels must be explicitly measured binary outcomes from one context. A
    returned model remains exploratory until separate held evidence passes the
    promotion policy; fitting and calibration alone do not promote it.
    """
    import numpy as np
    from rdkit import rdBase
    from sklearn.linear_model import LogisticRegression
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score

    if family not in {"logistic", "random_forest"}:
        raise FeedbackError("Choose a prespecified logistic or random_forest baseline")
    card: dict[str, Any] = {"status": "held", "promoted": False, "family": family, "seed": seed,
                            "n": len(records), "scope": "assay_specific_exploratory_prediction",
                            "versions": {"rdkit": rdBase.rdkitVersion,
                                         "scikit_learn": importlib.metadata.version("scikit-learn")},
                            "minimum_compounds": 200, "minimum_class_count": 50,
                            "minimum_partition_class_count": 20}
    try:
        if not 200 <= len(records) <= MAX_TRAINING_ROWS:
            raise FeedbackError("Model fitting requires 200 to 2000 unique measured parent compounds")
        contexts = set()
        seen = set()
        seen_structures = set()
        for row in records:
            context = tuple(row.get(key) for key in CONTEXT_FIELDS)
            if not all(context):
                raise FeedbackError("Target, endpoint, species, assay_id and unit must be explicit")
            contexts.add(context)
            if (row.get("status") != "measured" or row.get("relation") != "="
                    or isinstance(row.get("value"), (bool, np.bool_))
                    or row.get("value") not in (0, 1) or row.get("unit") != "binary"):
                raise FeedbackError("Only measured uncensored binary labels are eligible; missing is not inactive")
            if not row.get("parent_id") or row["parent_id"] in seen:
                raise FeedbackError("Provide one prespecified outcome per unique parent compound")
            seen.add(row["parent_id"])
            parent_key = _parent_key(row.get("smiles", ""))
            if row["parent_id"] != "parent:" + parent_key:
                raise FeedbackError("Training parent_id differs from the versioned chemical identity")
            if parent_key in seen_structures:
                raise FeedbackError("Salt or supplier duplicates cannot inflate model sample size")
            seen_structures.add(parent_key)
        if len(contexts) != 1:
            raise FeedbackError("Do not pool different assay, species, endpoint or target contexts")
        y = np.asarray([row["value"] for row in records], dtype=int)
        if min(Counter(y).get(0, 0), Counter(y).get(1, 0)) < 50:
            raise FeedbackError("Both measured classes require at least 50 unique compounds")
        groups, audit = leakage_groups(records)
        x = np.asarray(audit.pop("features"))
        audit.pop("fingerprints")
        descriptors = np.asarray(audit.pop("descriptors"))
        card["group_audit"] = audit
        parts = _partitions(records, groups, y, seed)
        card["partitions"] = {name: {"n": len(part), "groups": sorted({groups[int(i)] for i in part}),
                                    "candidate_ids": [records[int(i)]["candidate_id"] for i in part],
                                    "class_counts": {str(label): int(np.sum(y[part] == label)) for label in (0, 1)}}
                               for name, part in parts.items()}
        if any(min(np.sum(y[part] == label) for label in (0, 1)) < 20 for part in parts.values()):
            raise FeedbackError("Each train/calibration/test partition requires 20 measured examples per class")
        train, calibration, test = (parts[name] for name in ("train", "calibration", "test"))
        estimator = (LogisticRegression(C=1.0, max_iter=1000, random_state=seed) if family == "logistic"
                     else RandomForestClassifier(n_estimators=128, min_samples_leaf=3, max_features="sqrt",
                                                 n_jobs=1, random_state=seed))
        estimator.fit(x[train], y[train])
        raw_cal = np.clip(estimator.predict_proba(x[calibration])[:, 1], 1e-8, 1 - 1e-8)
        calibrator = LogisticRegression(C=1.0, random_state=seed).fit(np.log(raw_cal / (1 - raw_cal)).reshape(-1, 1), y[calibration])
        raw_test = np.clip(estimator.predict_proba(x[test])[:, 1], 1e-8, 1 - 1e-8)
        prediction = calibrator.predict_proba(np.log(raw_test / (1 - raw_test)).reshape(-1, 1))[:, 1]
        baseline = np.full(len(test), float(np.mean(y[train])))
        card.update(status="exploratory", context=dict(zip(CONTEXT_FIELDS, next(iter(contexts)))),
                    dataset_sha256=hashlib.sha256(canonical_json(records).encode()).hexdigest(),
                    metrics={"pr_auc": float(average_precision_score(y[test], prediction)),
                             "roc_auc": float(roc_auc_score(y[test], prediction)),
                             "brier": float(brier_score_loss(y[test], prediction)),
                             "baseline_pr_auc": float(average_precision_score(y[test], baseline)),
                             "baseline_brier": float(brier_score_loss(y[test], baseline))},
                    reason="Separate prespecified held benchmark evidence is required for promotion")
        artifact = {"schema_version": 1, "estimator": _serialize(estimator, family), "context": card["context"],
                    "versions": card["versions"], "dataset_sha256": card["dataset_sha256"],
                    "fingerprint": {"radius": 2, "bits": FINGERPRINT_BITS},
                    "calibration": {"method": "held_calibration_platt", "coefficient": float(calibrator.coef_[0, 0]),
                                    "intercept": float(calibrator.intercept_[0])},
                    "domain": {"training_smiles": [records[int(i)]["smiles"] for i in train],
                               "minimum_tanimoto": 0.4, "descriptor_min": descriptors[train].min(axis=0).tolist(),
                               "descriptor_max": descriptors[train].max(axis=0).tolist()},
                    "promoted": False, "status": "exploratory"}
        artifact["sha256"] = hashlib.sha256(canonical_json(artifact).encode()).hexdigest()
        card["model_sha256"] = artifact["sha256"]
        return {"card": card, "model": artifact}
    except (FeedbackError, ValueError, KeyError, TypeError) as exc:
        card["reason"] = str(exc)
        return {"card": card, "model": None}


def promotion_decision(model_card: dict, evidence: dict | None) -> dict:
    """Assess a separately frozen benchmark receipt, never test metrics alone."""
    reasons = []
    evidence = evidence or {}
    # Historical caller assertions are retained as input evidence, but cannot
    # establish a verified independent evaluation or authorize model promotion.
    reasons.append("Artifact-backed independent evaluation and a versioned promotion receipt are required")
    if model_card.get("status") != "exploratory" or not model_card.get("model_sha256"):
        reasons.append("No fitted auditable model")
    if evidence.get("model_sha256") != model_card.get("model_sha256"):
        reasons.append("Held evidence is not bound to this exact model")
    for field in ("prespecified", "data_integrity_verified", "leakage_audit_passed", "same_budget_comparison"):
        if evidence.get(field) is not True:
            reasons.append(f"Missing verified {field}")
    if evidence.get("context") != model_card.get("context"):
        reasons.append("Held evidence population differs from the model context")
    if not evidence.get("evaluation_snapshot_sha256") or evidence.get("evaluation_snapshot_sha256") == model_card.get("dataset_sha256"):
        reasons.append("An independent frozen evaluation dataset is required")
    if evidence.get("confidence_interval_unit") not in {"target_and_scaffold", "scaffold_and_source"}:
        reasons.append("Uncertainty must account for dependent scaffolds and sources or targets")
    for field, direction in (("pr_auc_improvement_ci_lower", "positive"), ("brier_difference_ci_upper", "negative")):
        value = evidence.get(field)
        if isinstance(value, bool) or not isinstance(value, (float, int)) or not math.isfinite(value) or (value <= 0 if direction == "positive" else value >= 0):
            reasons.append(f"Held improvement criterion not satisfied: {field}")
    return {"promoted": not reasons, "status": "eligible_for_versioned_promotion" if not reasons else "held",
            "reasons": reasons, "model_sha256": model_card.get("model_sha256")}


def predict_baseline(artifact: dict, smiles: list[str]) -> list[dict]:
    """Return clearly held exploratory predictions with applicability diagnostics."""
    import numpy as np
    from rdkit import DataStructs, rdBase

    payload = {k: v for k, v in artifact.items() if k != "sha256"}
    if artifact.get("sha256") != hashlib.sha256(canonical_json(payload).encode()).hexdigest():
        raise FeedbackError("Model artifact integrity mismatch")
    if artifact.get("versions", {}).get("rdkit") != rdBase.rdkitVersion:
        raise FeedbackError("Feature runtime version differs from the model")
    refs = [_fingerprint(s)[1] for s in artifact["domain"]["training_smiles"]]
    results = []
    for value in smiles:
        try:
            features, fp, _, descriptors = _fingerprint(value)
            nearest = max(DataStructs.BulkTanimotoSimilarity(fp, refs), default=0)
            domain = artifact["domain"]
            within = nearest >= domain["minimum_tanimoto"] and bool(np.all(descriptors >= domain["descriptor_min"]) and np.all(descriptors <= domain["descriptor_max"]))
            raw = float(np.clip(_raw_predict(artifact["estimator"], features.reshape(1, -1))[0], 1e-8, 1 - 1e-8))
            cal = artifact["calibration"]
            logit = cal["coefficient"] * math.log(raw / (1 - raw)) + cal["intercept"]
            prediction = 1 / (1 + math.exp(-max(-700, min(700, logit))))
            results.append({"smiles": value, "value": prediction, "kind": "predicted",
                            "status": "held_exploratory", "eligible_for_ranking": False,
                            "domain": "within" if within else "outside", "nearest_training_tanimoto": nearest,
                            "context": artifact["context"], "model_sha256": artifact["sha256"]})
        except (FeedbackError, TypeError, ValueError) as exc:
            results.append({"smiles": value, "value": None, "status": "unsupported", "reason": str(exc),
                            "eligible_for_ranking": False})
    return results
