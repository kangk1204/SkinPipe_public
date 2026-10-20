"""Separate the user's intended use from immutable predictor identity."""

import json


def load_policy(root):
    file = root / "config/admet_policy.json"
    return (
        json.loads(file.read_text())
        if file.exists()
        else {"use": "research_validation", "require_skin_validation": True}
    )


def auxiliary_complete(result, policy):
    return bool(
        policy.get("use") == "auxiliary_metrics"
        and result.get("base_prediction_complete")
        and result.get("software_integrity_passed")
        and result.get("hash_verification_performed")
    )
