"""Candidate purpose, separate from vendor availability and computed rank."""
from __future__ import annotations

from .fetch import mechanism as M

FUNCTIONAL_PRIORITY = 'functional_priority'
EXPLORATORY = 'exploratory'
COUNTER_CONTROL = 'counter_control'
POLICY_VERSION = 1


def compact_goal_support(goal: dict) -> dict:
    """Keep decision fields and identity, referring to the screen's raw sources.

    UniProt/ChEMBL receipts belong to ``target_design`` once per screen. They
    can contain large XML responses and must not be embedded in every candidate
    or CSV cell. This view leaves functional support fields unchanged.
    """
    view = dict(goal or {})
    identity = view.get('target_identity')
    if isinstance(identity, dict):
        view['target_identity'] = {k: v for k, v in identity.items() if k != 'sources'}
        view['target_identity_ref'] = '#/target_design/target_identity'
    return view


def classify(candidate) -> dict:
    """Choose a research purpose without claiming experimental efficacy."""
    ligand = candidate.ligand
    direction = candidate.direction or {}
    call = M.call_on(ligand)
    supported = M.evidence_supported(call, ligand, candidate.target)
    required = M.required_effect(direction.get('required'))
    verdict = (M.NO_REFERENCE if required is None else
               M.verdict(call.effect, required) if supported else M.UNKNOWN)
    goal = ligand.activity.get('goal_support') or {}
    actions = {a for e in call.evidence for a in e.get('action_types', [])}
    if verdict == M.OPPOSING:
        role, reason = COUNTER_CONTROL, 'reviewed direction opposes the requested node action'
    elif verdict != M.COMPATIBLE:
        role, reason = EXPLORATORY, 'direction is unconfirmed, conflicting or lacks a required action'
    elif 'PARTIAL AGONIST' in (call.action_type or '') or 'PARTIAL AGONIST' in actions:
        role, reason = EXPLORATORY, 'partial agonism requires assay context and endogenous-tone review'
    elif not (goal.get('status') == 'supported' and goal.get('endpoint')
              and goal.get('source_urls') and goal.get('node_gene') == candidate.target
              and M.required_effect(goal.get('required_action')) == required):
        role, reason = EXPLORATORY, 'target-to-marker direction and functional endpoint require review'
    else:
        role, reason = FUNCTIONAL_PRIORITY, 'reviewed direction and target-to-marker link support functional testing'
    return {'policy_version': POLICY_VERSION, 'role': role, 'reason': reason,
            'direction_supported': supported, 'direction_verdict': verdict,
            'goal_support': goal, 'experimental_validation': 'NOT_EVALUATED'}
