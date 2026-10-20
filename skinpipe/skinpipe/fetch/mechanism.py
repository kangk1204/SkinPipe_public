'Which way a compound pushes its target: ChEMBL `/mechanism` action_type.'
from __future__ import annotations

import logging
import hashlib
import json
import re
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from functools import lru_cache
from typing import Any, Optional, Sequence

from .cache import cached
from .chembl import CHEMBL, _c, _get

log = logging.getLogger(__name__)

# ----------------------------------------------------------------- vocabulary

#: What the compound does to the target's function.
UP = "up"
DOWN = "down"
UNCLEAR = "unclear"           # a row exists; its action_type carries no direction
CONFLICT = "conflict"         # rows exist for both directions; neither is claimed
NONE = "none"                 # no mechanism row for this (molecule, target)
UNAVAILABLE = "unavailable"   # the lookup did not answer; the direction is not known

#: How that sits against the direction the node has to be pushed.
COMPATIBLE = "compatible"
OPPOSING = "opposing"
UNKNOWN = "unknown"           # nothing recorded, nothing that maps, or no answer
NO_REFERENCE = "no_reference"  # skin_targets.yaml does not say which way this node goes

#: ChEMBL action_types that reduce the target's function.
DOWN_ACTIONS = frozenset({
    "ANTAGONIST", "INHIBITOR", "BLOCKER", "INVERSE AGONIST", "ALLOSTERIC ANTAGONIST",
    "NEGATIVE ALLOSTERIC MODULATOR", "NEGATIVE MODULATOR", "ANTISENSE INHIBITOR",
    "RNAI INHIBITOR", "DEGRADER", "DISRUPTING AGENT", "SEQUESTERING AGENT",
    "INACTIVATOR", "DOWNREGULATOR",
})
#: ChEMBL action_types that increase it. PARTIAL AGONIST is here rather than in
#: a class of its own: against a node that has to be inhibited, a partial
#: agonist is still an agonist, and the assay that would show it behaving as a
#: functional antagonist under high endogenous tone is not in this pipeline.
UP_ACTIONS = frozenset({
    "AGONIST", "PARTIAL AGONIST", "ACTIVATOR", "OPENER", "RELEASING AGENT",
    "POSITIVE ALLOSTERIC MODULATOR", "POSITIVE MODULATOR", "STABILISER",
    "STABILIZER", "UPREGULATOR",
})
#: Everything else with a row -- MODULATOR, SUBSTRATE, BINDING AGENT, OTHER --
#: falls through to UNCLEAR rather than being guessed at.

#: `node_action` values, and the synonyms a marker table tends to use.
_REQUIRED = {
    "inhibit": DOWN, "inhibition": DOWN, "inhibitory": DOWN, "block": DOWN,
    "blockade": DOWN, "suppress": DOWN, "suppression": DOWN, "antagonise": DOWN,
    "antagonize": DOWN, "decrease": DOWN, "down": DOWN,
    "activate": UP, "activation": UP, "agonise": UP, "agonize": UP,
    "induce": UP, "induction": UP, "increase": UP, "stimulate": UP, "up": UP,
}


@dataclass(frozen=True)
class ActionCall:
    """What is known about one compound's direction on one node.

    `source` says where the call came from, because the three sources are not
    equally strong: `chembl_mechanism` is a curated (molecule, target) row,
    `curated_note` is a human sentence in `config/cosmetic_actives.yaml`, and
    `none` / `unavailable` are the two ways of knowing nothing.
    """
    action_type: Optional[str]
    effect: str
    source: str                  # chembl_mechanism | curated_note | none | unavailable
    detail: str = ""
    evidence: tuple[dict, ...] = field(default_factory=tuple)


@lru_cache(maxsize=65536)
def structure_key(smiles: str) -> Optional[str]:
    """The complete standardized structure, including recorded stereochemistry."""
    from ..library import standardize, _key
    parent = standardize(smiles)
    return _key(parent) if parent else None


def evidence_supported(call: ActionCall, record: Any, node_gene: str) -> bool:
    """Require structure-, target- and endpoint-linked evidence for a call."""
    key = structure_key(record.smiles)
    return bool(key and call.effect in (UP, DOWN) and any(
        e.get('status') == 'reviewed'
        and e.get('structure_key') == key
        and e.get('node_gene') == node_gene
        and e.get('effect') == call.effect
        and e.get('endpoint') and e.get('source_url')
        and e.get('identity_basis') in ('molecule_id', 'structural_duplicate', 'reviewed_structure')
        for e in call.evidence))


def _reviewed_call(record: Any, node_gene: Optional[str]) -> Optional[ActionCall]:
    path = Path(__file__).resolve().parents[2] / 'config/direction_evidence.json'
    if not path.is_file() or not node_gene:
        return None
    raw = path.read_bytes()
    data = json.loads(raw)
    if not isinstance(data, dict) or data.get('schema_version') != 1:
        raise ValueError('unsupported direction evidence schema')
    rows = data.get('evidence')
    if not isinstance(rows, list) or any(not isinstance(row, dict) for row in rows):
        raise ValueError('invalid direction evidence records')
    key = structure_key(record.smiles)
    ids = {cid for tier in _candidate_ids(record) for cid in tier}
    hits = []
    for row in rows:
        if (row['node_gene'] == node_gene and row['molecule_id'] in ids
                and key and structure_key(row['smiles']) == key):
            hits.append(dict(row, structure_key=key, status='reviewed',
                             identity_basis='reviewed_structure',
                             ledger_sha256=hashlib.sha256(raw).hexdigest()))
    if not hits:
        return None
    directions = {r['effect'] for r in hits}
    effect = next(iter(directions)) if len(directions) == 1 else CONFLICT
    return ActionCall(' / '.join(sorted({r['action_type'] for r in hits})), effect,
                      'reviewed_assay', '; '.join(r['limitation'] for r in hits),
                      tuple(hits))


def effect_of(action_type: Optional[str]) -> str:
    """Map a ChEMBL action_type onto UP / DOWN, or UNCLEAR when it carries no
    direction. An empty action_type is NONE: the row said nothing."""
    if not action_type:
        return NONE
    a = " ".join(str(action_type).upper().split())
    if a in DOWN_ACTIONS:
        return DOWN
    if a in UP_ACTIONS:
        return UP
    return UNCLEAR


def required_effect(node_action: Optional[str]) -> Optional[str]:
    """The effect a node's `node_action` asks for, or None when it is unset.

    None means "nobody has recorded which way this node goes", and it must not
    be filled in from the marker's own increase/decrease column: for a proxy
    node the two can differ. Lowering a marker by *activating* an upstream
    brake is a normal design, and guessing `inhibit` there would flag exactly
    the right compound as wrong.
    """
    if not node_action:
        return None
    return _REQUIRED.get(str(node_action).strip().lower())


def verdict(effect: str, required: Optional[str]) -> str:
    """Compare one compound's effect with the direction the node needs."""
    if required is None:
        return NO_REFERENCE
    if effect not in (UP, DOWN):
        return UNKNOWN
    return COMPATIBLE if effect == required else OPPOSING


def explain(call: ActionCall, required: Optional[str], node_gene: str = "") -> str:
    """One sentence a reader can act on, for the drop list and the report."""
    node = node_gene or "this node"
    want = {DOWN: "inhibited", UP: "activated"}.get(required or "", "moved")
    v = verdict(call.effect, required)
    if v == OPPOSING:
        return (f"{call.action_type or call.effect} at {node}, but the screen needs "
                f"{node} {want} ({call.source})")
    if v == COMPATIBLE:
        return f"{call.action_type} at {node}; {node} needs to be {want}"
    if v == NO_REFERENCE:
        return f"no node_action recorded for {node}, so direction was not checked"
    if call.effect == CONFLICT:
        return (f"ChEMBL records both directions for this compound at {node} "
                f"({call.action_type}), so the direction was not called")
    if call.effect == UNAVAILABLE:
        return f"the mechanism lookup did not answer, so the direction at {node} is not known"
    if call.effect == UNCLEAR:
        return (f"ChEMBL records '{call.action_type}' at {node}, which carries no "
                "direction")
    return f"no mechanism recorded at {node}; direction unknown"


# -------------------------------------------------------- curated free text

#: English phrase -> the ChEMBL action_type it means. The *effect* is not
#: written down here: it comes back through `effect_of`, so a note and a ChEMBL
#: row that say the same word can never be graded differently. Longest phrases
#: first, so "inverse agonist" is consumed before the "agonist" inside it.
_NOTE_TERMS: tuple[tuple[str, str], ...] = (
    ("inverse agonist", "INVERSE AGONIST"), ("partial agonist", "PARTIAL AGONIST"),
    ("antagonist", "ANTAGONIST"), ("agonist", "AGONIST"),
    ("inhibitor", "INHIBITOR"), ("inhibition", "INHIBITOR"),
    ("suppression", "INHIBITOR"), ("blocker", "BLOCKER"),
    ("blockade", "BLOCKER"), ("activator", "ACTIVATOR"),
    ("activation", "ACTIVATOR"), ("induction", "ACTIVATOR"), ("opener", "OPENER"),
)
_NEGATED = re.compile(r"\b(not|no|non|never|rather than|instead of|without)\b[\s-]*$")
_CLAUSE_NEGATION = re.compile(r"\b(not|no|never|without|rather than|instead of)\b")
#: Words that sit between a protein name and the action noun and are not names.
_QUALIFIERS = frozenset({
    "a", "an", "the", "and", "or", "of", "as", "is", "its", "potent", "selective",
    "direct", "indirect", "weak", "strong", "best", "known", "classic", "natural",
    "human", "dual", "novel", "reversible", "irreversible", "competitive",
    "non-competitive", "allosteric", "partial", "inverse", "marketed", "topical",
})
_CHEMBL_ID = re.compile(r"\bCHEMBL\d+\b", re.IGNORECASE)


def _mentions_node(prefix: str, node_gene: Optional[str]) -> bool:
    'Whether an action noun in a curated note is about `node_gene`.'
    if node_gene is None:
        return True
    tail = prefix.rstrip()
    if not tail or tail[-1] in ";,(:":
        return True
    for token in reversed(re.split(r"[\s]+", tail)):
        word = token.strip("(),;:.").lower()
        if not word or word in _QUALIFIERS:
            continue
        parts = {p for p in re.split(r"[^A-Za-z0-9]+", word.upper()) if p}
        return node_gene.strip().upper() in parts
    return True


def action_from_note(note: Optional[str], node_gene: Optional[str] = None) -> ActionCall:
    'Read a direction out of a curated free-text note, conservatively.'
    text = (note or "").strip()
    if not text:
        return ActionCall(None, NONE, "none")
    low = text.lower()
    taken: list[tuple[int, int]] = []
    found: dict[str, str] = {}          # effect -> action_type
    for phrase, action in _NOTE_TERMS:
        eff = effect_of(action)
        for m in re.finditer(rf"\b{re.escape(phrase)}\b", low):
            if any(s < m.end() and m.start() < e for s, e in taken):
                continue
            taken.append((m.start(), m.end()))
            prefix = low[:m.start()]
            # Negation is clause-scoped, not an arbitrary character window.
            # Words between 'no' and an action do not establish that action.
            clause_prefix = re.split(r"[;.!?]", prefix)[-1]
            if _NEGATED.search(prefix) or _CLAUSE_NEGATION.search(clause_prefix):
                continue
            if not _mentions_node(prefix, node_gene):
                continue
            found.setdefault(eff, action)
    if len(found) == 1:
        eff, action = next(iter(found.items()))
        return ActionCall(action, eff, "curated_note", f"curated note: {text}")
    if len(found) > 1:
        return ActionCall(None, UNCLEAR, "curated_note",
                          f"the curated note names both directions: {text}")
    return ActionCall(None, NONE, "none")


# --------------------------------------------------------------- ChEMBL fetch

def _norm(action: Optional[str]) -> str:
    """One canonical spelling for an action_type, or "" when there is none."""
    return " ".join(str(action or "").upper().split())


def _resolve_actions(actions: Sequence[Optional[str]]) -> tuple[Optional[str], str, list[str]]:
    """Fold every action_type recorded for one compound into a single call.

    Order-invariant by construction: the input is reduced to a sorted *set*
    before anything is decided, so no arrival order is reachable from here.

    Rows that agree collapse to their shared direction -- AGONIST and PARTIAL
    AGONIST are one call, not two. Rows that disagree return CONFLICT and no
    action_type, because there is no honest way to choose between two curated
    rows that say opposite things. A direction-carrying row is not cancelled by
    a direction-less one: MODULATOR alongside AGONIST is still AGONIST, since
    "no direction recorded here" is not evidence of the other direction.

    Returns `(action_type, effect, every action_type seen)`. The third element
    is what lets the caller report *what* the disagreement was.
    """
    seen = sorted({_norm(a) for a in actions if _norm(a)})
    if not seen:
        return None, NONE, []
    directional = {effect_of(a) for a in seen} & {UP, DOWN}
    if len(directional) > 1:
        return None, CONFLICT, seen
    if directional:
        want = next(iter(directional))
        return next(a for a in seen if effect_of(a) == want), want, seen
    return seen[0], UNCLEAR, seen


def _row_actions(row: dict) -> list[str]:
    """Every action_type one table entry stands for.

    `fetch_for_target` writes `action_types` when it folds several rows into a
    single entry; a table built by hand -- a test fixture, an older cached one --
    carries only `action_type`. Both are read here so the two shapes resolve
    identically and a fixture cannot pass a table the live code would reject.
    """
    many = row.get("action_types")
    if many:
        return [_norm(a) for a in many if _norm(a)]
    one = _norm(row.get("action_type"))
    return [one] if one else []


def _merge_rows(rows: Sequence[dict]) -> dict:
    """The one table entry for an id, chosen without reference to row order."""
    action, effect, seen = _resolve_actions([r.get("action_type") for r in rows])
    pool = [r for r in rows if _norm(r.get("action_type")) == action] if action else []
    base = min(pool or list(rows),
               key=lambda r: (str(r.get("mechanism_of_action") or ""),
                              str(r.get("molecule_chembl_id") or ""),
                              str(r.get("parent_molecule_chembl_id") or "")))
    unique = {json.dumps(r, sort_keys=True): dict(r) for r in rows}
    return dict(base) | {"action_type": action, "action_types": seen,
                         "evidence_rows": [unique[k] for k in sorted(unique)],
                         "action_conflict": effect == CONFLICT}


_PAGE = 1000
_MAX_ROWS = 5000
_FIELDS = ("molecule_chembl_id", "parent_molecule_chembl_id", "action_type",
           "mechanism_of_action", "direct_interaction", "mechanism_comment",
           "selectivity_comment", "max_phase")


def fetch_for_target(target_chembl_id: str, cache=None) -> tuple[dict[str, dict], dict]:
    """Mechanism rows for one target, keyed by molecule id, plus a meta dict.

    Never raises. A transport failure returns an empty table with `ok=False`
    and the error, because the caller has to tell "ChEMBL did not answer" from
    "ChEMBL has no mechanism row for this compound", and an exception would
    take the whole screen down for the second -- far more common -- case.

    Rows are indexed under both `molecule_chembl_id` and
    `parent_molecule_chembl_id`: ChEMBL files the mechanism against the parent,
    while an activity row can name a salt, and a lookup by the activity's id
    would otherwise miss a mechanism the database does hold.

    An id that several rows land on gets one merged entry, resolved by
    `_merge_rows` from the *set* of action_types rather than from the first row
    to arrive. Two salts of one parent annotated in opposite directions meet on
    the parent key, and that entry says CONFLICT instead of whichever of the two
    ChEMBL happened to serialise first.

    Cached on the target id. A partial pull is returned but never stored, so an
    hour of 500s cannot become a permanently short mechanism table.
    """
    try:
        res = cached(_c(cache), "mechanism",
                     {"target": target_chembl_id, "max_rows": _MAX_ROWS,
                      "fields": ",".join(_FIELDS), "completion_schema": 2},
                     lambda: _mechanisms_live(target_chembl_id),
                     cacheable=lambda v: bool(v[1].get("complete")))
    except Exception as e:                       # ChemblError, or a transport flake
        log.warning("mechanism lookup for %s failed: %s", target_chembl_id, e)
        return {}, {"target_chembl_id": target_chembl_id, "ok": False,
                    "attempted": True, "n_rows": 0, "complete": False,
                    "error": f"{type(e).__name__}: {e}"}
    rows, meta = res.value
    grouped: dict[str, list[dict]] = {}
    for r in rows:
        for key in (r.get("molecule_chembl_id"), r.get("parent_molecule_chembl_id")):
            if key:
                grouped.setdefault(str(key).upper(), []).append(r)
    table: dict[str, dict] = {k: _merge_rows(v) for k, v in grouped.items()}
    out = dict(meta) | {"target_chembl_id": target_chembl_id, "ok": True,
                        "attempted": True, "error": "",
                        "from_cache": res.from_cache, "stale": res.stale,
                        "cache_age_days": (round(res.age_days, 1)
                                           if res.age_days is not None else None)}
    if res.stale:
        out["stale_note"] = res.note
    return table, out


def _mechanisms_live(target_chembl_id: str) -> tuple[list[dict], dict]:
    out: list[dict] = []
    offset = 0
    total: Optional[int] = None
    counts_consistent = True
    while offset < _MAX_ROWS:
        data = _get(f"{CHEMBL}/mechanism.json",
                    {"target_chembl_id": target_chembl_id,
                     "limit": min(_PAGE, _MAX_ROWS - offset), "offset": offset})
        page_meta = data.get("page_meta")
        raw_total = page_meta.get("total_count") if isinstance(page_meta, dict) else None
        page_total = (raw_total if type(raw_total) is int and raw_total >= 0 else
                      int(raw_total.strip()) if isinstance(raw_total, str)
                      and re.fullmatch(r"[0-9]+", raw_total.strip()) else None)
        if offset == 0:
            total = page_total
        if page_total is None or page_total != total:
            counts_consistent = False
        batch = data.get("mechanisms", [])
        if not batch:
            break
        out += [{k: m.get(k) for k in _FIELDS} for m in batch]
        offset += len(batch)
        if len(batch) < _PAGE:
            break
    # A successful HTTP response does not certify a complete target table.
    # Unknown totals, short pulls and the safety cap can hide a contradictory
    # row, so none can support an affirmative direction or an absence claim.
    complete = bool(counts_consistent and total is not None and len(out) == total)
    log.info("target %s: %d mechanism rows", target_chembl_id, len(out))
    return out, {"n_rows": len(out), "total_count": total, "complete": complete}


# ------------------------------------------------------------------ annotate

def _candidate_ids(record: Any) -> list[list[str]]:
    'The ids that may carry a mechanism row for this compound, by tier.'
    act = record.activity if isinstance(getattr(record, "activity", None), dict) else {}
    cited: list[str] = []
    for note_key in ("potency_note", "note"):
        cited += [m.group(0).upper() for m in _CHEMBL_ID.finditer(str(act.get(note_key) or ""))]
    merged = [str(i).upper() for i in (act.get("duplicate_ids") or []) if i]
    return [[str(record.ligand_id).upper()],
            list(dict.fromkeys(cited)), list(dict.fromkeys(merged))]


def _from_table(record: Any, table: dict[str, dict], *,
                node_gene: Optional[str] = None, target_id: Optional[str] = None,
                query_complete: bool = False) -> Optional[ActionCall]:
    """The ChEMBL call for a record, or None when none of its ids reaches a row.

    Each tier of `_candidate_ids` is resolved as a whole: every action_type its
    ids reach goes through `_resolve_actions`, so agreeing rows reinforce and
    disagreeing rows return CONFLICT. Nothing here depends on the order rows
    arrived in or on which id was tried first within a tier.
    """
    own = str(record.ligand_id).upper()
    tiers = _candidate_ids(record)
    for tier_index in (0, 2, 1):
        tier = tiers[tier_index]
        hits = [(cid, table[cid]) for cid in tier
                if cid in table and _row_actions(table[cid])]
        if not hits:
            continue
        action, effect, seen = _resolve_actions(
            [a for _cid, row in hits for a in _row_actions(row)])
        evidence = tuple({
            'status': ('reviewed' if target_id and tier_index != 1 and query_complete
                       else 'review_required'),
            'query_complete': query_complete,
            'structure_key': structure_key(record.smiles), 'node_gene': node_gene,
            'target_id': target_id, 'effect': effect,
            'endpoint': 'target function annotated by ChEMBL mechanism curation',
            'assay_context': 'not specified by mechanism record',
            'source_url': 'https://www.ebi.ac.uk/chembl/api/data/mechanism.json?target_chembl_id=' + str(target_id or ''),
            'molecule_id': cid,
            'identity_basis': ('molecule_id' if tier_index == 0 else
                               'note_chembl_id' if tier_index == 1 else 'structural_duplicate'),
            'action_types': _row_actions(row),
            'rows': row.get('evidence_rows', [row]),
        } for cid, row in hits)



        via = ", ".join(sorted({cid for cid, _row in hits if cid != own}))
        if effect == CONFLICT:
            why = ("ChEMBL records both directions for this compound at this "
                   f"target ({' / '.join(seen)}), so no direction is claimed")
            return ActionCall(" / ".join(seen), CONFLICT, "chembl_mechanism",
                              f"{why} [via {via}]" if via else why, evidence)
        moa = next((str(row.get("mechanism_of_action") or "") for _cid, row in hits
                    if action in _row_actions(row)), "")
        return ActionCall(action, effect, "chembl_mechanism",
                          f"{moa} [direction via {via}]" if via else moa, evidence)
    return None


def call_for(record: Any, table: dict[str, dict], meta: dict,
             node_gene: Optional[str] = None) -> ActionCall:
    """The best direction call available for one record, strongest source first."""
    from_table = _from_table(record, table, node_gene=node_gene,
                             target_id=meta.get('target_chembl_id'),
                             query_complete=meta.get('complete') is True)
    try:
        reviewed = _reviewed_call(record, node_gene)
    except (OSError, ValueError, KeyError, TypeError) as exc:
        # A broken local evidence ledger cannot certify a direction (nor the
        # absence of a contradictory reviewed call), but must not erase the
        # compound library or all previously acquired measurements.
        return ActionCall(None, UNAVAILABLE, "unavailable",
                          f"reviewed direction evidence unavailable: {type(exc).__name__}: {exc}")
    if reviewed is not None:
        if (from_table is not None and meta.get('complete') is True
                and from_table.effect in (UP, DOWN, CONFLICT)
                and from_table.effect != reviewed.effect):
            return ActionCall(None, CONFLICT, 'evidence_conflict',
                              'reviewed assay and mechanism records disagree',
                              reviewed.evidence + from_table.evidence)
        if from_table is not None and meta.get('complete') is not True:
            return ActionCall(reviewed.action_type, reviewed.effect, reviewed.source,
                              reviewed.detail + '; partial ChEMBL annotations retained for review',
                              reviewed.evidence + from_table.evidence)
        return reviewed
    if from_table is not None:
        return from_table
    act = record.activity if isinstance(getattr(record, "activity", None), dict) else {}
    from_note = action_from_note(act.get("note"), node_gene=node_gene)
    if from_note.effect != NONE:
        return from_note
    if meta.get("ok") and meta.get("complete") is True:
        return ActionCall(None, NONE, "none",
                          "ChEMBL records no mechanism for this compound at this target")
    return ActionCall(None, UNAVAILABLE, "unavailable",
                      meta.get("error") or ("the mechanism table was not completely retrieved"
                                           if meta.get("ok") else
                                           "the mechanism table was not retrieved"))


def call_on(record: Any) -> ActionCall:
    """Rebuild the call `annotate` wrote onto a record.

    A record that never went through `annotate` -- a hand-built LigandRecord in
    a test, an older cached screen -- comes back UNAVAILABLE rather than NONE,
    because "nobody asked" and "ChEMBL has nothing" are the two states this
    module exists to keep apart.
    """
    act = record.activity if isinstance(getattr(record, "activity", None), dict) else {}
    return ActionCall(act.get("action_type"),
                      act.get("action_effect") or UNAVAILABLE,
                      act.get("action_source") or "unavailable",
                      act.get("action_detail") or "",
                      tuple(act.get('direction_evidence') or ()))


def annotate(records: Sequence[Any], target_chembl_id: Optional[str] = None,
             node_gene: Optional[str] = None, cache=None) -> dict:
    """Write the direction fields onto every record and report the coverage.

    Each record gains `activity["action_type"]`, `["action_effect"]`,
    `["action_source"]` and `["action_detail"]`. The returned dict is the audit
    the screen prints: how many compounds carry a direction at all, where those
    calls came from, and -- when the table never arrived -- why. Coverage is the
    number that matters to a reader: on a real ChEMBL node it is single-digit
    percent, so a direction gate is a filter on a small, well-evidenced slice of
    the library and never a claim about the rest.
    """
    if target_chembl_id:
        table, meta = fetch_for_target(target_chembl_id, cache=cache)
    else:
        table, meta = {}, {
            "target_chembl_id": None, "ok": False, "attempted": False, "n_rows": 0,
            "complete": False,
            "error": ("no ChEMBL target id for this node, so the mechanism table "
                      "was never requested")}
    effects: Counter = Counter()
    sources: Counter = Counter()
    n_supported = 0
    ledger_errors = set()
    for r in records:
        call = call_for(r, table, meta, node_gene=node_gene)
        if not isinstance(getattr(r, "activity", None), dict):
            r.activity = {}
        r.activity["action_type"] = call.action_type
        r.activity["action_effect"] = call.effect
        r.activity["action_source"] = call.source
        r.activity["action_detail"] = call.detail
        r.activity['direction_evidence'] = list(call.evidence)
        if call.effect == UNAVAILABLE and call.detail.startswith("reviewed direction evidence unavailable:"):
            ledger_errors.add(call.detail)
        effects[call.effect] += 1
        sources[call.source] += 1
        n_supported += int(evidence_supported(call, r, node_gene or ''))
    n = len(records)
    known = effects[UP] + effects[DOWN]
    if ledger_errors:
        meta = dict(meta, ok=False, complete=False,
                    error="; ".join(sorted(ledger_errors)),
                    reviewed_evidence_status="UNAVAILABLE")
    return dict(meta) | {
        "node_gene": node_gene, "n_records": n, "n_with_direction": n_supported,
        'n_annotations_with_direction': known,
        "coverage": round(n_supported / n, 4) if n else 0.0,
        # Counted separately from `effects` because it is the one unknown that
        # means ChEMBL answered and contradicted itself. It is not a lookup
        # failure and it must not read as "nothing is recorded".
        "n_conflicting": effects[CONFLICT],
        "effects": dict(effects), "sources": dict(sources)}
