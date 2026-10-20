"""Bounded public-literature collection and human-reviewed claim proposals.

Literature retrieval and exact-quote checks never establish scientific validity,
target identity, efficacy, or permission to change a screening policy.
"""
from .collect import collect_evidence
from .review import validate_review

__all__ = ["collect_evidence", "validate_review"]
