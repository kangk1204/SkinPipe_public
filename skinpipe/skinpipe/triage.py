'Grade each marker by whether a screen against it can actually succeed.'
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Optional

from .types import DruggableNode, Grade, TargetMode, TargetSpec

log = logging.getLogger(__name__)

#: Grade thresholds on the count of potent AND topically plausible ligands.
THRESHOLDS: dict[Grade, int] = {Grade.A: 50, Grade.B: 15, Grade.C: 3}

#: Below this many potent ligands a node cannot support a meaningful screen,
#: however good the biology behind it is.
THIN_LIGAND_COUNT = 15


@dataclass
class TriageResult:
    marker: str
    grade: Grade
    provisional: bool                 # graded on raw counts, not filtered ligands
    screening_node: Optional[str]
    n_considered: int
    structure_ready: bool
    recommendation: str
    warnings: list[str] = field(default_factory=list)
    #: True when the node counts behind this grade were never retrieved, so the
    #: grade is an artefact of missing data rather than a finding about biology.
    counts_unresolved: bool = False
    #: Why `structure_ready` is what it is, for `screening_node` specifically.
    #: Callers about to prepare a receptor should read this, not the bool.
    docking: Optional["DockingVerdict"] = None

    @property
    def screenable(self) -> Optional[bool]:
        return None if self.counts_unresolved else self.grade is not Grade.D


def _own_node(spec: TargetSpec) -> Optional[DruggableNode]:
    own = next((n for n in spec.nodes if n.primary), None)
    if own is None:
        norm = spec.marker.upper().replace("-", "").replace("/", "")
        own = next((n for n in spec.nodes
                    if n.gene.upper() == norm or norm.endswith(n.gene.upper())), None)
    return own


def screening_metadata_unresolved(spec: TargetSpec, node: Optional[str]) -> list[DruggableNode]:
    """Missing lookups that can change the selection, rather than all missing lookups.

    A resolved configured route wins independently of unrelated alternatives.
    A resolved direct target with a ligand series also cannot be displaced.
    Other selections depend on comparative counts, so a failed alternative
    lookup must not be treated as a measured zero to freeze a group identity.
    """
    if spec.mode is TargetMode.ENDPOINT:
        return []
    missing = [n for n in spec.nodes if not n.resolved]
    if not missing:
        return []
    selected = next((n for n in spec.nodes if n.gene == node), None)
    if selected is None or not selected.resolved:
        return missing
    preferred = [n for n in spec.nodes if n.preferred]
    if preferred:
        return [n for n in preferred if not n.resolved]
    own = _own_node(spec)
    if spec.mode is TargetMode.DIRECT and selected is own and own.n_ligands_p6 > 0:
        return []
    return missing


class ScreeningMetadataUnresolved(ValueError):
    def __init__(self, marker: str, nodes: list[DruggableNode]):
        self.unresolved_nodes = [{"gene": n.gene, "error": n.resolve_error or "lookup failed"}
                                 for n in nodes]
        detail = "; ".join(f"{n['gene']}: {n['error']}" for n in self.unresolved_nodes)
        super().__init__(f"screening metadata unresolved for {marker}: {detail}; "
                         "resolve with a fresh Registry and regroup before acquisition")


def require_screening_metadata(spec: TargetSpec, node: Optional[str]) -> None:
    missing = screening_metadata_unresolved(spec, node)
    if missing:
        raise ScreeningMetadataUnresolved(spec.marker, missing)


def best_node(spec: TargetSpec) -> tuple[Optional[str], list[str]]:
    """Pick the node to screen against, and say so when the pick is a homologue.

    Preference order: nodes that have their own ligand series, then the one
    with the most experimental structures, then the most ligands.
    """
    if not spec.nodes:
        return None, []
    with_ligands = [n for n in spec.nodes if n.n_ligands_p6 > 0]
    pool = with_ligands or spec.nodes
    pick = max(pool, key=lambda n: (n.n_pdb, n.n_ligands_p6))
    warnings: list[str] = []







    preferred = [n for n in spec.nodes if n.preferred]
    if preferred:
        best_pref = max(preferred, key=lambda n: (n.n_pdb, n.n_ligands_p6))
        if not best_pref.resolved:
            warnings.append(
                f"{best_pref.gene} is the configured route for {spec.marker}; "
                f"its counts are unresolved ({best_pref.resolve_error}). "
                "Target identity is retained; no automatic subtype substitution.")
            return best_pref.gene, warnings
        if best_pref.gene != pick.gene:
            thin = best_pref.n_ligands_p6 < THIN_LIGAND_COUNT
            rich = pick.n_ligands_p6 >= THRESHOLDS[Grade.A]
            if thin and rich:
                # Biology says one thing, the chemistry cupboard says another.
                # Neither is a reason to overrule the other silently.
                warnings.append(
                    f"{best_pref.gene} is the configured route for {spec.marker} but has "
                    f"only {best_pref.n_ligands_p6} potent ligands, too few for a "
                    f"meaningful screen. {pick.gene} has {pick.n_ligands_p6}. Run both: "
                    f"{best_pref.gene} for mechanism, {pick.gene} for hit rate.")
            else:
                warnings.append(
                    f"{pick.gene} has more chemical matter ({pick.n_ligands_p6} ligands) "
                    f"than {best_pref.gene} ({best_pref.n_ligands_p6}), but {best_pref.gene} "
                    f"is the configured route for {spec.marker}. Screening {best_pref.gene}; "
                    f"widen to {pick.gene} only if it comes up empty.")
        pick = best_pref
        return pick.gene, warnings

    if not with_ligands:
        warnings.append(f"no node has a ligand series; {pick.gene} chosen on structures alone")

    # For a direct-mode marker the marker's own gene IS the target. A sibling
    # with more crystal structures is a homologue, not a substitute, so it only
    # wins when the real target has no chemical matter at all -- and then the
    # substitution is announced.
    own = _own_node(spec)
    if spec.mode is TargetMode.DIRECT and own is not None:
        if own.n_ligands_p6 > 0 and pick.gene != own.gene:
            warnings.append(
                f"{pick.gene} has more structures ({pick.n_pdb}) than {own.gene} "
                f"({own.n_pdb}), but screening stays on {own.gene}, the actual target.")
            pick = own
        elif own.n_ligands_p6 == 0 and pick.gene != own.gene:
            warnings.append(
                f"{own.gene} has no ligand series, so {pick.gene} is used as a "
                f"homology surrogate. Validate every hit on {own.gene} itself.")

        # Staying on the right target can still leave you with no structure to
        # dock into. Human tyrosinase is the standing example: one entry, while
        # its paralogue TYRP1 has a dozen. Say so, because it changes the plan.
        if not own.structure_usable and own.n_pdb:
            warnings.append(
                f"{own.gene} has {own.n_pdb} PDB entr(ies) but none is usable for docking"
                + (f": {own.structure_note}" if own.structure_note else "")
                + ". Treat this target as structure-free.")
        if pick.gene == own.gene:
            richer = [n for n in spec.nodes
                      if n.gene != own.gene and n.n_pdb >= max(3 * own.n_pdb, own.n_pdb + 3)]
            for r in richer:
                warnings.append(
                    f"{own.gene} has only {own.n_pdb} experimental structure(s) while "
                    f"{r.gene} has {r.n_pdb}. Screening stays on {own.gene}; if docking "
                    f"into it proves unreliable, build a homology model on {r.gene} and "
                    f"say so in the report.")
    return pick.gene, warnings


def metal_warning(node) -> Optional[str]:
    """Vina cannot score what holds a chelating ligand in place.

    Redocking tropolone into mushroom tyrosinase (2Y9X) puts the top-ranked
    pose 2.87 A from the crystal. Raising exhaustiveness to 64 *finds* a 0.69 A
    pose, but Vina ranks it 8th of 9, losing by 0.28 kcal/mol, because the
    copper chelation that defines the binding mode contributes nothing to the
    score. The same applies to the zinc-dependent MMPs. Docking such a target
    without saying this produces a hit list that looks fine and is not.
    """
    if not getattr(node, "metalloenzyme", None):
        return None
    return (f"{node.gene} is a {node.metalloenzyme}-dependent metalloenzyme. Vina has "
            f"no metal-coordination term, so it will under-rank exactly the chelating "
            f"ligands that work here. Keep the metal in the receptor, raise "
            f"exhaustiveness, rescore the top poses with MM/GBSA, and never rank on "
            f"the Vina score alone.")


class DockingRefused(RuntimeError):
    """Raised when a caller asks to dock a node that may not be docked.

    Carries the verdict's reason as its message, so a stage that lets it
    escape reports *why* rather than "no receptor".
    """


@dataclass(frozen=True)
class DockingVerdict:
    'The answer to "may I dock this node?", with the reasoning attached.'
    node: Optional[str]
    allowed: bool
    reason: str
    n_pdb: int = 0
    #: The richest sibling node whose structures ARE usable. Recorded only on a
    #: refusal, because that is the only case where it changes what happens next.
    homology_template: Optional[str] = None
    template_n_pdb: int = 0
    template_note: str = ""

    def __bool__(self) -> bool:
        """`if may_dock(...)` reads as permission, so it must mean permission.

        Every dataclass instance is truthy by default, which would make a
        refusal read as a yes at every call site that forgot `.allowed` -- the
        same silent-yes shape as the bug this type exists to close.
        """
        return self.allowed

    def raise_if_refused(self) -> "DockingVerdict":
        """Turn a refusal into an exception; return self when docking is allowed.

        For callers whose failure path is a try/except that records a skipped
        stage: one line replaces a hand-written structure check that can drift
        out of step with this module.
        """
        if not self.allowed:
            raise DockingRefused(self.reason)
        return self


def _template_for(spec: TargetSpec, node_gene: str,
                  min_pdb: int) -> Optional[DruggableNode]:
    """Richest sibling whose structures are usable, or None."""
    return max((n for n in spec.nodes
                if n.gene != node_gene and n.structure_usable and n.n_pdb >= min_pdb),
               key=lambda n: n.n_pdb, default=None)


def may_dock(spec: TargetSpec, node: Optional[str] = None,
             min_pdb: int = 1) -> DockingVerdict:
    'Whether the node that will actually be screened can be docked into.'
    if node is None:
        node, _ = best_node(spec)
    if node is None:
        return DockingVerdict(
            node=None, allowed=False,
            reason=f"no druggable node is registered for {spec.marker}, so there "
                   f"is nothing to dock into.")
    picked = next((n for n in spec.nodes if n.gene == node), None)
    if picked is None:
        raise ValueError(
            f"{node} is not a node of {spec.marker}; its nodes are "
            f"{', '.join(n.gene for n in spec.nodes) or '(none)'}. Ask about a node "
            f"the registry actually carries.")

    def refuse(reason: str) -> DockingVerdict:
        tmpl = _template_for(spec, picked.gene, min_pdb)
        if tmpl is not None:
            note = (f"{tmpl.gene} has {tmpl.n_pdb} usable experimental structure(s) and "
                    f"may serve as a homology-model template for {picked.gene}. That is "
                    f"a template, not a substitute: say so in the report and confirm "
                    f"every hit on {picked.gene} itself.")
        elif len(spec.nodes) > 1:
            note = (f"no sibling node of {picked.gene} carries a usable structure "
                    f"either, so there is no template for a homology model.")
        else:
            note = (f"{picked.gene} is the only node registered, so there is no "
                    f"template for a homology model.")
        return DockingVerdict(node=picked.gene, allowed=False, reason=reason,
                              n_pdb=picked.n_pdb,
                              homology_template=tmpl.gene if tmpl else None,
                              template_n_pdb=tmpl.n_pdb if tmpl else 0,
                              template_note=note)

    # Curated first. `structure_usable=False` is a human statement about what
    # those entries contain, and it stays true whether or not the counts came
    # back this morning.
    if not picked.structure_usable:
        return refuse(
            f"{picked.gene} has {picked.n_pdb} PDB entr(ies) mapped to it, none usable "
            f"for docking"
            + (f": {picked.structure_note}" if picked.structure_note else "")
            + ".")
    if not picked.resolved:
        return refuse(
            f"{picked.gene}'s structure count was never retrieved "
            f"({picked.resolve_error or 'lookup failed'}), so whether it can be docked "
            f"is unknown. This is missing data, not a finding of no structure: re-run "
            f"the lookup rather than treating {picked.gene} as structure-free.")
    if picked.n_pdb < min_pdb:
        return refuse(
            f"{picked.gene} has {picked.n_pdb} experimental structure(s), fewer than the "
            f"{min_pdb} required, so there is nothing to dock into.")
    return DockingVerdict(
        node=picked.gene, allowed=True, n_pdb=picked.n_pdb,
        reason=f"{picked.gene} has {picked.n_pdb} experimental structure(s) and the "
               f"registry records no reason to distrust them.")


def structure_ready(spec: TargetSpec, min_pdb: int = 1,
                    node: Optional[str] = None) -> bool:
    'Whether the node that will be screened has a structure worth docking into.'
    return may_dock(spec, node=node, min_pdb=min_pdb).allowed































FEASIBLE_PCHEMBL = 7.0
FEASIBLE_TOPICAL = 7.0


def grade(spec: TargetSpec, n_feasible: Optional[int] = None,
          node: Optional[str] = None) -> Grade:
    'Grade the marker on the node that will actually be screened.'
    if spec.mode is TargetMode.ENDPOINT or not spec.nodes:
        return Grade.D
    if n_feasible is not None:
        n = n_feasible
    elif node is not None:
        picked = next((nd for nd in spec.nodes if nd.gene == node), None)
        n = picked.n_ligands_p6 if picked else 0
    else:
        n = max((nd.n_ligands_p6 for nd in spec.nodes), default=0)
    for g in (Grade.A, Grade.B, Grade.C):
        if n >= THRESHOLDS[g]:
            return g
    return Grade.D


def assess(spec: TargetSpec, n_feasible: Optional[int] = None) -> TriageResult:
    """Full triage verdict for one marker."""
    node, warnings = best_node(spec)
    g = grade(spec, n_feasible, node=node)
    provisional = n_feasible is None and spec.mode is not TargetMode.ENDPOINT and bool(spec.nodes)
    verdict = may_dock(spec, node=node)
    ready = verdict.allowed
    picked_node = next((nd for nd in spec.nodes if nd.gene == node), None)
    n_considered = (n_feasible if n_feasible is not None
                    else (picked_node.n_ligands_p6 if picked_node else 0))

    # A node whose counts were never retrieved has n_pdb = n_ligands_p6 = 0,
    # which is exactly what a genuinely barren target looks like. Grading the
    # two the same turns a ChEMBL outage into "no viable small-molecule route",
    # reported with full confidence. Say which one this is.
    unresolved = screening_metadata_unresolved(spec, node)
    blind = bool(unresolved) and n_feasible is None

    if spec.mode is TargetMode.ENDPOINT:
        rec = "use as an efficacy endpoint; do not screen"
    elif blind:
        why = next((n.resolve_error for n in unresolved if n.resolve_error), "lookup failed")
        rec = (f"cannot grade {spec.marker}: the node counts were never retrieved "
               f"({why}). This is missing data, not a negative result -- re-run "
               f"once UniProt/ChEMBL/RCSB answer.")
    elif g is Grade.D:
        rec = "no viable small-molecule route; consider peptides, biologics or another node"
    elif not ready:
        # The reason travels with the verdict: "no structure" and "the structures
        # are the wrong protein" are both ligand-based, and only the second one
        # leaves a homology template on the table.
        rec = f"screen ligand-based against {node}; {verdict.reason}"
    elif g is Grade.A:
        rec = f"dock against {node}; chemical matter is plentiful"
    elif g is Grade.B:
        rec = f"dock against {node}; expect to optimise hits for skin delivery"
    else:
        rec = f"dock against {node}; chemical matter is thin, widen the library"

    picked = next((n for n in spec.nodes if n.gene == node), None)
    if picked is not None:
        mw = metal_warning(picked)
        if mw:
            warnings.append(mw)
    if blind:
        warnings.insert(0, (
            f"grade {g.value} is NOT a finding: "
            + ", ".join(f"{n.gene} ({n.resolve_error or 'lookup failed'})"
                        for n in unresolved)
            + ". Zero ligands and zero structures here mean the lookup did not "
              "return, not that the target is barren. Do not act on this row."))
    elif provisional:
        warnings.append(
            "grade is provisional: computed from raw ChEMBL counts. Run the topical "
            "filter and re-grade, because potent does not imply skin-deliverable.")
    if spec.notes:
        warnings.extend(spec.notes)

    return TriageResult(marker=spec.marker, grade=g, provisional=provisional or blind,
                        screening_node=node, n_considered=n_considered,
                        structure_ready=ready, recommendation=rec, warnings=warnings,
                        counts_unresolved=blind, docking=verdict)
