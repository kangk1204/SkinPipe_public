"""End-to-end orchestration: marker in, ranked topical candidates out.

Stage order is deliberate: the cheap, decisive filters run before the expensive
ones, because there is no point spending an hour of MM/GBSA on a molecule that
cannot cross skin. Docking and MM/GBSA are optional and the screen completes
without them, saying which stage it skipped.

A stage that did not run must never read as a stage that ran and found nothing.
Anything a reader would act on -- the buy list, the hazard gate, the library
size -- records that difference in `ScreenResult.gates_unverified`.
"""
from __future__ import annotations

import math

import logging
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Optional

from . import library, score, sourcing, strategy, triage
from .fetch import chembl, mechanism, provenance
from .filters import safety as safety_f
from .filters import stability as stability_f
from .filters import admet as admet_f
from .filters import topical as topical_f
from .registry import Registry
from .types import (SEVERITY_ORDER, Candidate, LigandRecord, SkinAssessment,
                    TargetMode, TargetSpec, evidence_for, hazard_severity_of,
                    normalise_severity_gate)

log = logging.getLogger(__name__)

_NO_SOURCING = ("the vendor/patent lookup %s, so 'not purchasable' below "
                "means 'not checked'")


#: Where a CAS number can be found on a ligand today. `LigandRecord` has no
#: CAS field -- that is the finding, not an oversight in this lookup -- so the
#: only place one can arrive is the free-form `props` dict, written by whatever
#: stage supplied the record. A user SDF and the cosmetic-actives file can both
#: carry one; nothing in the ChEMBL path does. Spelled out here rather than
#: guessed at the call site, and kept narrow: a wrong CAS matches a different
#: substance on a prohibition list, which is worse than matching nothing.
_CAS_PROP_KEYS = ("cas", "cas_number", "cas_rn", "casrn", "CAS")

#: A registry number the *input file* asserted, which nothing here has checked
#: against the structure. Kept apart from the keys above: those are either a
#: curated file's own field or what `sourcing` resolved from the structure.
_CAS_CLAIM_KEY = "cas_supplied"


def _note_receptor_preparation(result, receptor) -> None:
    """Promote structural limitations, including historical cards, to gates.

    Newly discarded rebuilt residues contained no deposited heavy atoms. That
    is unavailable geometry, not evidence that observed atoms were deleted.
    Historical discard notes retain their original, stronger scope verbatim.
    """
    evidence = getattr(receptor, "preparation_evidence", {}) or {}
    limits = []
    if evidence.get("discarded_rebuilt_residue_count", 0):
        limits.append(f"{evidence['discarded_rebuilt_residue_count']} wholly rebuilt residue(s) "
                      "were discarded; unobserved geometry remains unavailable. "
                      "No deposited heavy-atom loss is asserted by this record.")
    gaps = evidence.get("unobserved_gaps") or []
    if gaps:
        numbered = [gap for gap in gaps
                    if gap.get("residue_count_evidence") == "absent_author_numbered_positions_only"]
        sequence = [gap for gap in gaps
                    if gap.get("residue_count_evidence") != "absent_author_numbered_positions_only"]
        if sequence:
            count = sum(gap.get("residue_count", 0) for gap in sequence)
            limits.append(f"{count} wholly unobserved residue(s) remain unmodelled; "
                          "pocket completeness and fragment-end chemistry are not validated.")
        if numbered:
            count = sum(gap.get("residue_count", 0) for gap in numbered)
            limits.append(f"{count} author-numbered position(s) absent across observed backbone fragments; "
                          "missing residue identities and sequence count are not established from numbering. "
                          "Pocket completeness and fragment-end chemistry are not validated.")
    # Ordinary terminal trimming and a SEQRES coverage note do not themselves
    # establish discarded internal geometry. Keep them in the preparation
    # notes; actual gaps/discards above and old explicit loss notes still gate.
    prefixes = ("RECEPTOR RESIDUES DISCARDED", "REBUILT RESIDUES DISCARDED", "left ")
    for note in receptor.notes or []:
        if note.startswith(prefixes) and (not note.startswith("left ") or "unmodelled" in note):
            limits.append(note)
    if limits:
        result.gates_unverified["receptor_residues"] = " ".join(dict.fromkeys(limits))


def _note_redock_outcome(result, rd: dict, pass_mark: float) -> None:
    "Record what the redock control established, in the control's own words."
    # Kept as the record, not only as the sentence below. The four
    # `failure_kind` values drive triage and the prose cannot be queried.
    result.redock = dict(rd)
    rmsd = rd.get("rmsd")
    # Formatted, not interpolated raw. A float straight out of `CalcRMS` prints
    # as 9.837543134170579, and this string is read by a person deciding
    # whether to trust a whole ranking.
    shown = f"{rmsd:.2f}" if isinstance(rmsd, (int, float)) else "n/a"
    kind = rd.get("failure_kind")
    # FAIL is a claim about this control and needs a number behind it. The
    # test is the number's presence, not a list of kinds, because the list
    # grows: `unmatchable` was once the only not-measured value and the engine
    # failing, the engine returning nothing, and there being no crystal ligand
    # at all all printed "FAIL, in-place RMSD n/a A" -- a verdict on a box that
    # had never been compared to anything.
    measured = isinstance(rmsd, (int, float))
    outcome = ("pass" if rd.get("pass")
               else "FAIL" if measured else "NOT MEASURED")
    # What the control was measured against, not only what it measured. A pass on
    # a fragment is not a pass on the ligand, and the RMSD alone cannot say which
    # it was; `reference_*` comes from the extraction through the Receptor.
    ref = ""
    if rd.get("reference_is_fragment"):
        atoms = rd.get("reference_atoms") or (None, None)
        pct = rd.get("reference_completeness")
        ref = (f", against a FRAGMENT reference "
               f"({atoms[0]} of {atoms[1]} component heavy atoms"
               f"{f', {100 * pct:.0f}%' if isinstance(pct, (int, float)) else ''})")
    result.stages_run.append(
        f"redock control ({outcome}, in-place RMSD {shown} A{ref})")
    if rd.get("pass"):
        if rd.get("reference_is_fragment"):
            # A pass has to be qualified rather than trusted: the geometry was
            # confirmed against part of the ligand, so the part that is missing
            # was never tested.
            atoms = rd.get("reference_atoms") or (None, None)
            result.gates_unverified["redock_reference"] = (
                f"the redock control passed, but against a fragment of "
                f"{rd.get('ligand')}: {atoms[0]} of {atoms[1]} heavy atoms from "
                "the chemical component were usable. Only the compared fragment "
                "was reproduced in this run; the rest of the ligand was not tested. "
                "Other ligands and affinity ranking were not validated.")
        return
    why = rd.get("note") or (f"the redock control did not pass "
                             f"(pass mark {pass_mark} A)")
    if kind == "unmatchable":
        result.gates_unverified["redock"] = (
            f"{why}. So the control neither confirms nor refutes the box: "
            "raising exhaustiveness will not help, because no RMSD was "
            "computed. The poses below are unvalidated rather than "
            "discredited, and validating them needs a control ligand whose "
            "graph can be matched.")
    elif not measured:
        # Same consequence as unmatchable -- the box is unestablished -- but a
        # different cause, and the cause is what says whether re-running helps.
        # An engine error clears on a retry; no crystal ligand never will.
        result.gates_unverified["redock"] = (
            f"{why} No RMSD was obtained, so the box and receptor preparation "
            "are unvalidated rather than discredited: the poses below are not "
            "refuted by this control, and nothing in this screen confirms them "
            "either.")
    else:
        result.gates_unverified["redock"] = (
            f"{why} This control does not establish the geometry or affinity "
            "ranking of the other candidates.")


def _note_scoring_unreliable(result: "ScreenResult", receptor, rd: dict) -> None:
    'Say when the Vina *ordering* is the thing that cannot be trusted.'
    if rd.get("failure_kind") != "scoring":
        return
    try:
        from .structure import metal_sites
        metals = metal_sites.receptor_metals(Path(receptor.path))
    except Exception as exc:                                  # noqa: BLE001
        log.debug("could not read metals for %s: %s", receptor.pdb_id, exc)
        return
    if not metals:
        return
    names = ", ".join(sorted({m.resname.strip() for m in metals}))
    result.gates_unverified["vina_ranking"] = (
        f"the redock failed on scoring, not search, and this pocket carries "
        f"{len(metals)} catalytic metal atom(s) ({names}). Vina has no "
        "metal-coordination term, so it under-ranks the chelators that bind "
        "here -- on the control ligand it generated a near-native pose and "
        "ranked it below a wrong one. MM/GBSA would normally correct this and "
        "`metal_guard` refuses it for the same missing metal, so nothing in "
        "this screen orders the list correctly. Read the candidates as a set "
        "that reached the pocket, not as a ranking, and rank them on assay.")


def _ligand_cas(lig: LigandRecord) -> list[str]:
    """Every CAS number this ligand carries, in preference order.

    A list, not one number. `sourcing.pubchem_cas` returns all the registry
    numbers PubChem lists for a structure -- resveratrol has two, niacinamide
    three -- and which one 별표 1 filed the substance under is not knowable
    here, so `eligibility` is handed all of them. `cas_all` is where the
    sourcing stage puts them; the singular keys stay supported because a user
    SDF or the cosmetic-actives file can carry one and nothing else does.
    """
    out: list[str] = []
    many = lig.props.get("cas_all")
    if isinstance(many, (list, tuple)):
        out += [str(v).strip() for v in many if str(v).strip()]
    for key in _CAS_PROP_KEYS:
        value = lig.props.get(key)
        if isinstance(value, str) and value.strip():
            out.append(value.strip())
    # Last, and deliberately included. A number the user's own file asserted is
    # not evidence that it describes this structure, so it must never be what
    # makes a verdict *definitive* -- `_cas_is_only_a_claim` below sees to that.
    # It is still worth looking up: a match against 별표 1 can only add a
    # restriction, and refusing to look was the strictly worse failure.
    claim = lig.props.get(_CAS_CLAIM_KEY)
    if isinstance(claim, str) and claim.strip():
        out.append(claim.strip())
    seen: list[str] = []
    for c in out:
        if c not in seen:
            seen.append(c)
    return seen


def _cas_is_only_a_claim(lig: LigandRecord) -> bool:
    """True when every CAS this ligand has came from the input file.

    The verdict then rests on an identifier nothing verified against the
    structure, so it is reported as not definitive whichever way it goes.
    """
    if not isinstance(lig.props.get(_CAS_CLAIM_KEY), str):
        return False
    many = lig.props.get("cas_all")
    if isinstance(many, (list, tuple)) and any(str(v).strip() for v in many):
        return False
    return not any(isinstance(lig.props.get(k), str) and lig.props[k].strip()
                   for k in _CAS_PROP_KEYS)


def _fmt_below(value: float, threshold: float) -> tuple[str, str]:
    '`value` and `threshold` printed at the least precision that shows the gap.'
    v, t = float(value), float(threshold)
    for places in range(1, 7):
        sv, st = f"{v:.{places}f}", f"{t:.{places}f}"
        if float(sv) < float(st):
            return sv, st
    # Only reachable if the caller formats a value that is not below the
    # threshold; print both in full rather than inventing a comparison.
    return repr(v), repr(t)


def _hazard_concern_names(admet, severity: str) -> list[str]:
    """Which alerts set `severity`, so the hold names the finding.

    Reads the serialised `to_dict()` form the pipeline stores on a candidate
    and the `ADMETProfile` object alike, because both shapes reach this code
    path -- the same asymmetry that kept rule 6 of `CosmeticEligibility` dark
    for every real screen. Systemic exposure is excluded for the reason
    `hazard_severity_of` excludes it: a dose is not a liability.
    """
    if admet is None:
        return []
    concerns = (admet.get("concerns", []) if isinstance(admet, dict)
                else getattr(admet, "concerns", []))
    names = set()
    for c in concerns:
        fields = c if isinstance(c, dict) else {
            "severity": getattr(c, "severity", None),
            "domain": getattr(c, "domain", None),
            "name": getattr(c, "name", None)}
        if (fields.get("severity") == severity
                and fields.get("domain") != "systemic"
                and fields.get("name")):
            names.add(fields["name"])
    return sorted(names)


def _admet_note(admet) -> str:
    """The ADMET stage's own sentence, appended to a hold reason.

    A severity is a grade; this is what the grade was for, and a reader who
    wants to argue with the hold needs the second one too.
    """
    if admet is None:
        return ""
    overall = (admet.get("overall") if isinstance(admet, dict)
               else getattr(admet, "overall", ""))
    return f". ADMET verdict: {overall[:160]}" if overall else ""


@dataclass
class ScreenResult:
    marker: str
    triage: "triage.TriageResult"
    candidates: list[Candidate] = field(default_factory=list)
    all_candidates: list[Candidate] = field(default_factory=list)
    filter_profile: dict = field(default_factory=dict)
    breakdowns: list[score.ScoreBreakdown] = field(default_factory=list)
    provenance_audit: dict = field(default_factory=dict)
    decision: object | None = None          # strategy.StrategyDecision
    dropped: list[dict] = field(default_factory=list)   # what was removed, and why
    held_for_regulatory: list[dict] = field(default_factory=list)  # not barred, not cleared
    #: The direction gate's own audit: which way the node had to be pushed, how
    #: many compounds carried a mechanism at all, and what was removed. Empty
    #: when the screen stopped before the library stage.
    direction: dict = field(default_factory=dict)
    target_design: dict = field(default_factory=dict)
    control_check: dict = field(default_factory=dict)
    stages_run: list[str] = field(default_factory=list)
    stages_skipped: dict[str, str] = field(default_factory=dict)
    #: Gates whose verdict is not evidence. A gate that never ran answers "no"
    #: for every compound, which reads exactly like one that rejected everything.
    gates_unverified: dict[str, str] = field(default_factory=dict)
    #: The pChEMBL floor this screen ran at and where it came from.
    #: Carried on the result because a library size cannot be read
    #: without it: "14 compounds" is a barren target or a strict
    #: threshold, and only the floor separates the two.
    potency_floor: dict[str, object] = field(default_factory=dict)
    #: What `check_library_scale` found before the library was fetched: how many
    #: activities ChEMBL holds at the chosen floor, whether `max_chembl` will
    #: cap them, and a floor that would fit if so. Recorded even when the answer
    #: is "fine", so a reader can tell a pull that was checked and found whole
    #: from one nobody looked at.
    library_scale: dict[str, object] = field(default_factory=dict)
    #: How much of the gated library was actually docked. Recorded whenever
    #: docking ran, including when it ran on everything, so "all of it" and
    #: "nobody looked" stay distinguishable -- the same reason `library_scale`
    #: is written even when the pull was whole.
    dock_scale: dict[str, object] = field(default_factory=dict)








    redock: dict = field(default_factory=dict)
    structure_validation: dict = field(default_factory=dict)
    export_contract: dict = field(default_factory=dict)
    receptor: object | None = None

    @property
    def trustworthy(self) -> bool:
        return not self.gates_unverified

    def top(self, n: int = 20) -> list[Candidate]:
        return sorted([c for c in self.candidates if c.rank],
                      key=lambda c: c.rank)[:n]

    def partition(self, max_hazard: str = "moderate", *, purpose: str = 'functional_priority'
                  ) -> tuple[list[Candidate], list[dict]]:
        """``(orderable, held)`` for these candidates. Computes, never stores.

        The buy list and the held list are one decision read two ways, so they
        are produced together by a function of the candidates alone: calling it
        twice returns the same two lists and changes nothing on the result.
        `held_for_regulatory` used to exist *only* as a side effect of
        `orderable()`, so anything that read the field without calling
        `orderable()` first -- `report.to_json` on a fresh result, a notebook,
        the next stage in a script -- saw an empty list and could not tell
        "nothing was held" from "nobody has looked yet".
        """
        order = SEVERITY_ORDER
        if purpose not in ('functional_priority', 'research_available', 'formulation_review'):
            raise ValueError('unknown candidate purpose: ' + str(purpose))
        # Unlike a *gate* threshold, "none" is a legitimate value here: it is
        # the strictest setting, "only compounds with no hazard alert at all".
        # What is not legitimate is a value this list does not contain -- the
        # old `else len(order)` silently cut above every severity and returned
        # an empty buy list, which reads as a finding about the chemistry.
        want = max_hazard.strip().lower() if isinstance(max_hazard, str) else max_hazard
        if want not in order:
            raise ValueError(
                f"unknown max_hazard {max_hazard!r}; expected one of "
                f"{', '.join(order)}. An unrecognised value would empty the buy "
                "list, which is indistinguishable from a screen that found "
                "nothing worth buying.")
        cut = order.index(want)
        out: list[Candidate] = []
        held: list[dict] = []

        def hold(c: Candidate, stage: str, verdict: str, reason: str) -> None:
            held.append({"ligand_id": c.ligand.ligand_id,
                         "name": c.ligand.name or c.ligand.ligand_id,
                         "stage": stage, "verdict": verdict, "reason": reason})

        for c in sorted([x for x in self.candidates if x.rank], key=lambda x: x.rank):
            if not c.ligand.props.get("purchasable"):
                continue
            haz = hazard_severity_of(c.admet)













            if haz == "unknown":
                hold(c, "hazard", "hazard not assessed",
                     "the ADMET stage produced no verdict for this structure, "
                     "so there is no hazard finding to clear it"
                     + _admet_note(c.admet))
                continue
            if haz in order and order.index(haz) < cut:
                names = _hazard_concern_names(c.admet, haz)
                hold(c, "hazard", f"hazard {haz}",
                     f"worst non-systemic hazard is '{haz}', above the "
                     f"'{want}' ceiling this buy list was cut at"
                     + (f" ({', '.join(names)})" if names else "")
                     + _admet_note(c.admet))
                continue
            # Direction. `_direction_gate` removes an opposing compound only
            # when it is running in a dropping mode; with the gate flagging
            # rather than dropping -- the default when no `node_action` is
            # recorded, and the documented `direction_gate=None` setting -- an
            # opposing compound survives into the candidate list carrying its
            # own "this pushes the node backwards" verdict, and this list
            # ignored the verdict and ranked it to buy. A compound ChEMBL says
            # opposes the node is not a weak candidate, it is the wrong one,
            # and no vendor or regulatory clearance redeems it.
            direction = c.direction or {}
            if purpose == 'functional_priority' and direction.get("verdict") == mechanism.OPPOSING:
                hold(c, "direction", "direction opposing",
                     direction.get("explanation")
                     or f"pushes {c.target} the way this marker must not go")
                continue
            # `usable` is eligible or restricted. Excluding only the two
            # barred verdicts let `data_required` and `unknown` through, so a
            # compound whose regulatory review had not finished appeared on a
            # list captioned "can be legally formulated". Buying a research
            # reagent and formulating a cosmetic ingredient are different
            # decisions and this list makes the second claim.
            verdict = (c.eligibility or {}).get("verdict")
            if verdict not in ("eligible", "restricted"):
                hold(c, "eligibility", verdict or "not assessed",
                     (c.eligibility or {}).get("reason", ""))
                continue
            # The CAS is the key 별표 1·2 has real coverage on -- 632 of 656 rows
            # against 12 structures -- so a compound whose CAS lookup failed has
            # not been compared against the prohibition list on the one key that
            # would have found it. This list claims the compound can be legally
            # formulated; that claim is not available for a compound nobody
            # could look up. `False` only, never a missing key: absent means the
            # sourcing stage did not run, which the screen records once.
            if (c.eligibility or {}).get("cas_checked") is False:
                hold(c, "eligibility", "CAS not checked",
                     "PubChem did not answer for this compound's registry "
                     "numbers, so the prohibition list was never consulted on "
                     "the key it has coverage on. Re-run before ordering.")
                continue
            if purpose == 'functional_priority':
                from .candidate_use import classify, FUNCTIONAL_PRIORITY
                use = classify(c)
                if use['role'] != FUNCTIONAL_PRIORITY:
                    hold(c, 'direction', 'functional priority unestablished', use['reason'])
                    continue
            out.append(c)
        return out, held

    def regulatory_holds(self, max_hazard: str = "moderate") -> list[dict]:
        """Candidates neither barred nor cleared, computed on demand."""
        return self.partition(max_hazard=max_hazard)[1]

    def orderable(self, n: Optional[int] = None,
                  max_hazard: str = "moderate") -> list[Candidate]:
        'Available candidates supported for functional validation of the goal.'
        out, held = self.partition(max_hazard=max_hazard)
        self.held_for_regulatory = held
        return out[:n] if n else out

    def research_available(self, max_hazard: str = 'moderate') -> list[Candidate]:
        """Research sourcing shortlist. Unassessed hazards remain unassessed."""
        max_hazard = max_hazard.strip().lower() if isinstance(max_hazard, str) else max_hazard
        if max_hazard not in SEVERITY_ORDER:
            raise ValueError(f"unknown max_hazard {max_hazard!r}")
        ceiling = SEVERITY_ORDER.index(max_hazard)
        return [c for c in sorted(self.candidates, key=lambda c: c.rank or float("inf"))
                if c.ligand.props.get("purchasable") is True
                and (hazard_severity_of(c.admet) == "unknown"
                     or SEVERITY_ORDER.index(hazard_severity_of(c.admet)) >= ceiling)]

    def research_candidates(self, purpose: str) -> list[Candidate]:
        from .candidate_use import classify, FUNCTIONAL_PRIORITY, EXPLORATORY, COUNTER_CONTROL
        if purpose not in (FUNCTIONAL_PRIORITY, EXPLORATORY, COUNTER_CONTROL):
            raise ValueError('unknown candidate purpose: ' + str(purpose))
        if purpose == FUNCTIONAL_PRIORITY:
            return self.orderable()
        return [c for c in sorted((c for c in self.candidates if c.rank), key=lambda c: c.rank)
                if classify(c)['role'] == purpose]










DIRECTION_MODES = ("opposing", "strict")


def validate_direction_gate(mode: Optional[str]) -> Optional[str]:
    'The direction mode to run, or None for "flag but do not drop".'
    if mode is None:
        return None
    if not isinstance(mode, str):
        raise ValueError(
            f"direction_gate must be a mode name or None, not "
            f"{type(mode).__name__} ({mode!r}). Expected one of "
            f"{', '.join(DIRECTION_MODES)}, 'none', or None.")
    value = mode.strip().lower()
    if value == "none":
        return None
    if value not in DIRECTION_MODES:
        raise ValueError(
            f"unknown direction_gate {mode!r}; expected one of "
            f"{', '.join(DIRECTION_MODES)}, or 'none'/None to flag opposing "
            "compounds without dropping them. An unrecognised mode would run "
            "the gate in flag-only mode and report the screen as trustworthy, "
            "so it is refused instead.")
    return value


def _required_action(reg: Registry, spec: TargetSpec, node_gene: str) -> Optional[str]:
    'The `node_action` recorded for this node in config/skin_targets.yaml.'
    for n in spec.nodes:
        if n.gene.upper() == node_gene.upper() and getattr(n, "node_action", None):
            return str(n.node_action)
    return None


#: Used when neither the caller nor the registry names a floor. Lower than the
#: 6.5 a drug programme would use, because this pipeline looks for compounds
#: that go into a formulation at percent concentrations rather than ones that
#: have to work at nanomolar. It is a noise floor, not a potency requirement.
DEFAULT_POTENCY_FLOOR = 5.0


def validate_topical_cutoff(value: Any) -> float:
    """A finite topical score cutoff shared by the CLI and Python API.

    Range comparisons alone accept NaN: every compound then fails the gate
    while its reason says that a finite score is less than NaN.
    """
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError("min_topical must be a finite number between 0 and 10")
    if not math.isfinite(value) or not 0 <= value <= 10:
        raise ValueError("min_topical must be a finite number between 0 and 10")
    return float(value)


def _valid_screen_count(value: Any, name: str, *, minimum: int = 0,
                        allow_none: bool = False) -> Optional[int]:
    """Refuse ambiguous budgets before acquisition or slicing a candidate pool."""
    if value is None and allow_none:
        return None
    if type(value) is not int or value < minimum:
        description = "nonnegative" if minimum == 0 else "positive"
        raise ValueError(f"{name} must be a {description} integer"
                         + (" or None" if allow_none else ""))
    return value


def _valid_floor(value: Any, where: str) -> float:
    """A pChEMBL floor, or a refusal that names the offender.

    `float(value)` accepted a YAML typo as a crash and everything else in
    silence: `potency_floor: micromolar` raised deep inside the screen, where
    `run_one` filed it under "unreachable -- says nothing about the target",
    and the retry pass then spent up to three hours waiting for a ChEMBL that
    was healthy. `potency_floor: 99` did not even raise -- the pull returned
    nothing, `audit["n"]` was 0 so no truncation detector fired, and the screen
    proceeded on the curated actives with only a stage-line suffix to show for
    it. A closed-vocabulary field gets a closed-vocabulary error, the way
    `normalise_severity_gate` already does.
    """
    if isinstance(value, bool) or not isinstance(value, (int, float, str)):
        raise ValueError(f"{where} must be a number; got {value!r}")
    try:
        out = float(value)
    except (TypeError, ValueError):
        raise ValueError(f"{where} must be a number; got {value!r}") from None
    if not math.isfinite(out):
        raise ValueError(f"{where} must be finite; got {value!r}")
    if not 0.0 <= out <= 14.0:
        raise ValueError(
            f"{where} must be a pChEMBL value between 0 and 14; got {out:g}. "
            "A floor above ~10 returns nothing and reads as a barren target.")
    return out


#: Floors to probe when a pull would be capped, lowest first. Not a search --
#: `activity_count` is one cached call each, and the point is to hand the caller
#: a floor that fits rather than a number to solve for.
_FLOOR_LADDER: tuple[float, ...] = (5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0)


def check_library_scale(node_gene: str, floor: float, max_chembl: int,
                        reg: Optional[Registry] = None,
                        target_chembl_id: Optional[str] = None,
                        include_without_pchembl: bool = False) -> dict:
    'Would this pull be capped? Asked *before* the library is fetched.'
    out: dict = {"checked": False, "capped": None, "available": None,
                 "max_chembl": max_chembl, "floor": floor}
    effective_max = None if max_chembl == 0 else max_chembl
    out["effective_max_chembl"] = effective_max
    out["activity_scope"] = "all_recorded" if include_without_pchembl else "pchembl_floor"
    try:
        # `Registry.spec()` has already resolved this node -- that is what
        # `DruggableNode.chembl_target` holds -- so re-resolving it here was a
        # second UniProt round trip for an answer in hand, and it made this a new
        # network dependency for every caller of `run_screen`. 25 tests that mock
        # the library fetch came to depend on a live UniProt because of it, and
        # passed anyway: the `except Exception` below turned the failure into
        # "could not check". The caller passes what it resolved.
        tid = target_chembl_id
        if not tid:
            # Deliberately not resolved here. This is an advisory pre-check, and
            # making it *initiate* a UniProt lookup turned an optional warning
            # into a mandatory network call for every caller of `run_screen` --
            # which is how 25 tests that mock the library fetch came to depend on
            # a live UniProt, and pass anyway because the handler below reported
            # "could not check". `Registry.spec()` resolves the node when it can
            # and the caller passes what it got; when it could not, the honest
            # answer is that the cap was not checked.
            out["reason"] = (
                f"{node_gene} has no resolved ChEMBL target, so whether this "
                f"pull will be capped at max_chembl={max_chembl} was not "
                "checked. The library stage reports truncation after the fact.")
            return out
        count_args = {"include_without_pchembl": True} if include_without_pchembl else {}
        available = chembl.activity_count(tid, pchembl_min=floor, **count_args)
    except Exception as e:                                    # noqa: BLE001
        # An outage here must not read as "this pull is fine". The screen goes
        # ahead -- the check is an early warning, not a gate -- but it says it
        # could not look rather than implying it looked and found nothing.
        out["reason"] = (f"could not check ({type(e).__name__}: {e}); whether "
                         "this pull will be capped is unknown")
        return out

    if max_chembl == 0:
        effective_max = available
        out["effective_max_chembl"] = available
    out.update(checked=True, target_chembl_id=tid, available=available,
               capped=available > effective_max)
    if include_without_pchembl:
        out["floor_that_fits"] = None
        out["reason"] = (f"{available} recorded activities including missing pChEMBL and "
                         f"censored endpoints; effective max_chembl={effective_max}. "
                         + ("This is a database-order slice; raise max_chembl to preserve "
                            "the whole query. A potency floor would change its scope."
                            if out["capped"] else "The configured cap covers this query."))
        return out
    if not out["capped"]:
        out["reason"] = (f"{available} activities at pChEMBL >= {floor:g}, "
                         f"within max_chembl={effective_max}")
        return out

    fits = None
    for candidate in _FLOOR_LADDER:
        if candidate <= floor:
            continue
        try:
            n = chembl.activity_count(tid, pchembl_min=candidate)
        except Exception:                                     # noqa: BLE001
            break
        if n <= effective_max:
            fits = (candidate, n)
            break
    out["floor_that_fits"] = fits
    advice = (f"raise the floor to pChEMBL >= {fits[0]:g} ({fits[1]} activities) "
              f"or max_chembl to at least {available}" if fits else
              f"raise max_chembl to at least {available}")
    out["reason"] = (
        f"{node_gene} has {available} activities at pChEMBL >= {floor:g} but "
        f"effective max_chembl={effective_max}, so {available - effective_max} will be dropped "
        f"in ChEMBL's storage order -- not by potency. The library will be an "
        f"arbitrary cross-section of this target's chemistry. Before spending "
        f"the run: {advice}.")
    return out


def _potency_floor(reg: Registry, spec: TargetSpec, node_gene: str,
                   explicit: Optional[float]) -> tuple[float, str]:
    """The pChEMBL floor for this node, and where the number came from.

    Order: an explicit argument, then the node's own `potency_floor`, then
    :data:`DEFAULT_POTENCY_FLOOR`. The source string is returned rather than
    logged because a library size is meaningless without it -- "14 compounds"
    reads as a barren target and as a strict threshold, and only the floor
    tells the two apart.
    """
    if explicit is not None:
        return _valid_floor(explicit, "the pchembl_min argument"), "caller"
    for n in spec.nodes:
        if (n.gene.upper() == node_gene.upper()
                and getattr(n, "potency_floor", None) is not None):
            return (_valid_floor(n.potency_floor,
                                 f"potency_floor on node {node_gene}"),
                    f"registry ({node_gene})")
    # No `reg._raw` fallback here. `Registry.spec` populates
    # `DruggableNode.potency_floor`, so a second lookup into the private raw
    # YAML could only ever agree -- or drift, which is exactly what
    # `_required_action`'s docstring says got that field into trouble. Worse,
    # it masked the primary path: with `Registry.spec` reverted so the field
    # was never populated, the resolver test still passed on the fallback.
    return DEFAULT_POTENCY_FLOOR, "default"


def _direction_gate(result: ScreenResult, ligands: list[LigandRecord], node: str,
                    node_action: Optional[str], mech: dict,
                    mode: Optional[str]) -> list[LigandRecord]:
    'Remove the compounds whose known mechanism pushes the node backwards.'
    required = mechanism.required_effect(node_action)
    from .candidate_use import compact_goal_support
    # Refused here as well as in `run_screen`, because this is a module-level
    # entry point in its own right: a caller that reaches it directly gets the
    # same answer rather than a quietly disarmed gate.
    mode = validate_direction_gate(mode)
    # With no recorded direction there is nothing to compare against, so no
    # mode can run -- "strict" would otherwise drop the entire library for
    # failing a test that was never defined.
    active = mode if required is not None else None
    counts: Counter = Counter()
    kept: list[LigandRecord] = []
    for ligand_record in ligands:
        call = mechanism.call_on(ligand_record)
        supported = mechanism.evidence_supported(call, ligand_record, node)
        v = (mechanism.verdict(call.effect, required) if supported else
             mechanism.NO_REFERENCE if required is None else mechanism.UNKNOWN)
        why = mechanism.explain(call, required, node)
        if not supported and call.effect in (mechanism.UP, mechanism.DOWN):
            why = 'direction annotation requires structure/target/endpoint review: ' + why
        ligand_record.activity["direction_verdict"] = v
        ligand_record.activity['direction_supported'] = supported
        ligand_record.activity['goal_support'] = compact_goal_support(result.target_design)
        ligand_record.activity["direction_required"] = node_action or ""
        ligand_record.activity["direction_explanation"] = why
        counts[v] += 1
        drop = active is not None and (
            v == mechanism.OPPOSING if active == "opposing"
            else v != mechanism.COMPATIBLE)
        if not drop:
            kept.append(ligand_record)
            continue
        ligand_record.flags.append(f"direction {v}: {why}")
        # Retain excluded candidates with their direction evidence.
        result.dropped.append({
            "ligand_id": ligand_record.ligand_id, "name": ligand_record.name or ligand_record.ligand_id,
            "stage": "direction", "severity": v, "reason": why, "smiles": ligand_record.smiles,
            'direction_evidence': list(call.evidence),
            'candidate_use': 'counter_control' if v == mechanism.OPPOSING else 'exploratory'})

    n = len(ligands)
    lookup_ok = bool(mech.get("ok"))
    result.direction = {
        "node": node, "node_action": node_action or None, "required": required,
        "mode": active, "n_screened": n, "n_dropped": n - len(kept),
        "opposing": counts[mechanism.OPPOSING],
        "compatible": counts[mechanism.COMPATIBLE],
        "unknown": counts[mechanism.UNKNOWN],
        "no_reference": counts[mechanism.NO_REFERENCE],
        "n_with_direction": mech.get("n_with_direction", 0),
        "coverage": mech.get("coverage", 0.0),
        "evidence_sources": mech.get("sources", {}),
        "effects": mech.get("effects", {}),
        "lookup": {"ok": lookup_ok, "attempted": bool(mech.get("attempted")),
                   "n_rows": mech.get("n_rows", 0), "stale": bool(mech.get("stale")),
                   "error": mech.get("error", "")},
    }

    if required is None:
        # `node_action: unknown` is a recorded judgement -- the sign is
        # contested -- and an absent key is an unfilled row. Both leave the
        # gate with nothing to compare against, and neither is "fine either
        # way", so the message names which one happened.
        why = (f"records node_action: {node_action} for {node}, so the sign is not "
               "settled" if node_action else f"records no node_action for {node}")
        msg = (f"config/skin_targets.yaml {why}. No compound in this screen was "
               f"checked for direction: {counts[mechanism.NO_REFERENCE]} candidates "
               f"carry a potency against {node} and an unstated relationship to what "
               "the marker needs.")
        result.direction["verdict"] = msg
        result.stages_skipped["direction"] = msg
        # Keeping these compounds is the documented exploration policy and is
        # not what is recorded here. What is recorded is that the direction
        # gate could not be *defined*: with no required sign there is no test,
        # so every candidate below is unscreened for direction and the run is
        # not entitled to `trustworthy`. A gate that never ran clears nothing,
        # and this branch is the one way a screen reached a reader claiming
        # otherwise -- `gates_unverified` was empty, so the report printed no
        # "These results are not complete" section at all.
        kept_phrase = ("the one candidate kept below is unscreened for it"
                       if len(kept) == 1 else
                       f"all {len(kept)} candidates kept below are unscreened for it")
        result.gates_unverified["direction"] = (
            f"the direction gate could not be defined. config/skin_targets.yaml "
            f"{why}. No candidate in this screen was checked against the direction "
            f"{result.marker or 'this marker'} needs, and {kept_phrase}. Fill in "
            f"node_action for {node} and re-run before reading this list as "
            "direction-screened.")
        return kept

    want = "inhibited" if required == mechanism.DOWN else "activated"
    line = (f"direction gate ({node} must be {want}): {counts[mechanism.OPPOSING]} "
            f"opposing {'dropped' if active else 'flagged'}, "
            f"{counts[mechanism.COMPATIBLE]} compatible, "
            f"{counts[mechanism.UNKNOWN]} of {n} with no usable mechanism")
    if not lookup_ok:
        line += f" -- ChEMBL mechanism table unavailable ({mech.get('error', '')[:80]})"
    result.stages_run.append(line)
    result.direction["verdict"] = line






    n_checked = counts[mechanism.COMPATIBLE] + counts[mechanism.OPPOSING]
    short = f"{n_checked} of {n} candidates against {node} still carry a direction, " \
            f"from curated notes alone" if n_checked else \
            f"no candidate was checked against the {node} direction this marker needs"
    if lookup_ok and mech.get("complete") is False:
        result.gates_unverified["direction"] = (
            f"the ChEMBL mechanism table for {node} came back short "
            f"({mech.get('n_rows')} of {mech.get('total_count')} rows), so a compound "
            "that opposes may have been missed")
    elif not lookup_ok and mech.get("attempted"):
        result.gates_unverified["direction"] = (
            f"the ChEMBL mechanism lookup for {node} failed "
            f"({mech.get('error', '')}); {short}. Re-run when the API is back.")
    elif not lookup_ok:
        # No ChEMBL target id reached `annotate`, so the table was never
        # requested at all. Silent here, this reads as a clean sweep.
        result.gates_unverified["direction"] = (
            f"the ChEMBL mechanism table for {node} was never requested "
            f"({mech.get('error', '')}); {short}")


















    if n and n_checked * 2 < n:
        pct = 100.0 * n_checked / n
        # What the visible minority did, because that is node-specific rather
        # than boilerplate: PDGFRB has to go *up* and all 6 of its annotated
        # ligands go down. It is not evidence about the unchecked rows -- 6 of
        # 721 cannot be -- but it is the only direction evidence this node has,
        # and a reader deciding what to order should see which way it points.
        agree = counts[mechanism.COMPATIBLE]
        against = counts[mechanism.OPPOSING]
        seen = (f"Of the {n_checked} that could be read, {agree} push the way "
                f"this marker needs and {against} push against it. "
                if n_checked else "")
        result.gates_unverified["direction_coverage"] = (
            f"the direction gate checked {n_checked} of {n} candidates "
            f"({pct:.1f}%) against the direction {result.marker or 'this marker'} "
            f"needs; {n - n_checked} carry no usable mechanism record. {seen}"
            f"So \"{against} opposing\" counts only what could be "
            f"seen, and is not a statement that the rest agree. Direction for "
            f"those rows is unestablished, not confirmed -- an activity assay "
            f"is what settles it, not this ranking.")
    return kept


def _emptied_by(result: ScreenResult, fallback: str) -> str:
    """Name the stage that actually emptied the screen, not the one it stopped at.

    Every filter fell through to one shared `if not kept` after the sourcing
    block, which reported `no compound survived the sourcing filters` whatever
    had removed the last compound -- including runs with `check_sourcing=False`,
    where the sourcing filters did not execute at all. The detailed `dropped`
    list said `topical` and the one-line summary said `sourcing`, so the two
    halves of the same report named different causes and the summary is the
    half a reader acts on.

    `result.dropped` is appended in stage order by `_filter_out` and the gates,
    so its last entry is the stage that took the final compound.
    """
    if not result.dropped:
        return fallback
    stage = str(result.dropped[-1].get("stage") or "").strip()
    if not stage:
        return fallback
    n = sum(1 for d in result.dropped if d.get("stage") == stage)
    return (f"no compound survived the {stage} stage, which dropped {n} of "
            f"the {len(result.dropped)} compounds this screen removed")


def _filter_out(result: ScreenResult, ligands: list[LigandRecord], *,
                keep, stage: str, severity: str, reason) -> list[LigandRecord]:
    'Apply a filter and name every compound it removes.'
    out: list[LigandRecord] = []
    for ligand_record in ligands:
        if keep(ligand_record):
            out.append(ligand_record)
            continue
        why = reason(ligand_record) if callable(reason) else str(reason)
        ligand_record.flags.append(f"{stage} {severity}: {why}")
        result.dropped.append({
            "ligand_id": ligand_record.ligand_id, "name": ligand_record.name or ligand_record.ligand_id,
            "stage": stage, "severity": severity, "reason": why,
            "smiles": ligand_record.smiles})
    return out


def _direction_of(lig: LigandRecord) -> Optional[dict]:
    """The direction call carried on a surviving candidate, or None if the gate
    never ran. A candidate that survived is not a candidate that was cleared:
    `verdict` is "unknown" for most of them and says so."""
    a = lig.activity
    if "direction_verdict" not in a:
        return None
    return {"verdict": a["direction_verdict"], "action_type": a.get("action_type"),
            "effect": a.get("action_effect"), "source": a.get("action_source"),
            "required": a.get("direction_required") or None,
            'supported': bool(a.get('direction_supported')),
            'evidence': a.get('direction_evidence') or [],
            "explanation": a.get("direction_explanation", "")}


def assess_skin(ligand: LigandRecord, node_gene: Optional[str] = None) -> SkinAssessment:
    'Run the three skin filters and merge them into one assessment.'
    a = topical_f.assess(ligand)
    sens = safety_f.summarise(ligand.smiles, target_gene=node_gene)
    stab = stability_f.summarise(ligand.smiles)
    a.sensitisation_status = sens["assessment_status"]
    a.stability_status = stab["assessment_status"]
    a.sensitiser_alerts = [f"{x['severity']}:{x['name']}" for x in sens["alerts"]]
    a.sensitiser_details = [dict(x) for x in sens["alerts"]]
    a.stability_alerts = [f"{x['mode']}:{x['name']}" for x in stab["alerts"]]
    a.colour = stab["colour"]
    a.n_pains = len(sens["pains"]) if sens["assessment_status"] == "EVALUATED" else None
    if sens["pains"]:
        a.sensitiser_alerts.append(f"low:PAINS({a.n_pains})")
    notes = [a.verdict]
    if sens["verdict"] != "no electrophilic alert":
        notes.append(sens["verdict"])
    if stab["verdict"] != "no stability flag":
        notes.append(stab["verdict"])
    a.verdict = "; ".join(notes)
    return a


def _prepare_docking_receptor(prep_mod, spec: TargetSpec, node: str,
                              pdb_id: Optional[str], workdir: Path,
                              ligand_code: Optional[str] = None):
    'Resolve the entry to dock into and prepare it as the *intended* protein.'
    identity: dict[str, str] = {}
    if pdb_id is None:
        triage.may_dock(spec, node=node).raise_if_refused()
        uniprot = next((n.uniprot for n in spec.nodes if n.gene == node), None)
        if not uniprot:
            raise RuntimeError(f"no UniProt accession registered for node {node}")
        pid = (prep_mod.best_structure_for_uniprot(uniprot) or {}).get("pdb_id")
        identity = {"uniprot": uniprot, "gene": node}
    else:
        pid = pdb_id
        selected = next(n for n in spec.nodes if n.gene == node)
        if selected.goal_support.get("target_identity") is not None:
            if not selected.structure_usable:
                triage.may_dock(spec, node=node).raise_if_refused()
            # A supplied structure must still represent the requested protein.
            identity = {"uniprot": selected.uniprot, "gene": node}
    if not pid:
        raise RuntimeError(f"no usable experimental structure for {node}")
    if ligand_code:
        identity["ligand_code"] = ligand_code
    # Bound each preparation in an isolated process. The same budget is passed
    # to exact-accession alternatives, and explicit structures never substitute.
    identity["preparation_timeout"] = 1800
    alternative_preparer = getattr(prep_mod, "prepare_target_receptor", None)
    try:
        return pid, prep_mod.prepare_receptor(pid, workdir, **identity)
    except (ValueError, RuntimeError) as exc:
        if pdb_id is not None or ligand_code or not callable(alternative_preparer):
            raise
        receptor = alternative_preparer(pid, workdir, initial_failure=exc, **identity)
        return receptor.pdb_id, receptor


def _attach_gbsa_records(result: ScreenResult, posed: list[Candidate], gdf) -> None:
    """Map one batch frame onto the candidates it was computed for.

    A failed pose occupies a row with nothing in its numeric columns, so the
    frame carries ``NaN`` there. ``int(row.get("n_frames") or 1)`` raised on
    it -- ``NaN`` is truthy and ``int(NaN)`` is a ``ValueError`` -- and that
    exception escaped the per-candidate loop, so one unparameterisable ligand
    cost the whole batch its rescoring while leaving the candidates mapped
    before it already carrying results. A failed pose is a failed pose: it is
    recorded against its own candidate, and every other pose still delivers.

    The conversion is :func:`skinpipe.mmgbsa.records.to_result`, the same one
    ``mmgbsa_single`` uses, so the two APIs cannot drift apart again.
    """
    from .mmgbsa.records import to_result

    # Keyed on (ligand, pose), not the ligand alone. `to_frame` sorts by
    # `dG_bind` ascending -- best first -- so a ligand with several scored
    # poses had every earlier row overwritten by the next and the dict kept the
    # *worst* one. `run_screen` poses each candidate once so it cannot reach
    # that today, but `mmgbsa_batch` accepts many poses per ligand and this
    # function is what reads its output back.
    rows = gdf.to_dict("records")
    lookup = {(r["ligand_id"], int(r["pose_index"])): r for r in rows}
    # Older batch frames used pose_index=-1 for a dead process. That row does
    # not identify an individual pose; apply it only to missing exact rows for
    # the same ligand, never overwrite successful or individual failed rows.
    process_failures = {r["ligand_id"]: r for r in rows if int(r["pose_index"]) < 0}
    scored: list[str] = []
    failed: list[str] = []
    for c in posed:
        row = lookup.get((c.ligand.ligand_id, int(c.pose.pose_index)))
        if row is None:
            process_failure = process_failures.get(c.ligand.ligand_id)
            row = {**(process_failure or {}), "ligand_id": c.ligand.ligand_id,
                   "pose_index": int(c.pose.pose_index), "converged": False,
                   "dG_bind": None,
                   "notes": (process_failure.get("notes") if process_failure else None)
                            or "FAILED: batch returned no MM/GBSA record for this requested pose"}
        c.gbsa = to_result(row)
        if c.gbsa.dG_bind is None or not c.gbsa.converged:
            reason = (c.gbsa.notes or "no dG_bind returned").split(";")[0]
            failed.append(f"{c.ligand.ligand_id} ({reason.strip()[:120]})")
        else:
            scored.append(c.ligand.ligand_id)

    result.stages_run.append(
        f"MM/GBSA rescoring (top {len(posed)}: {len(scored)} scored, "
        f"{len(failed)} failed)")
    if failed:
        # Not a gate, so not `gates_unverified` -- but a pose the reader will
        # otherwise see ranked on docking alone with no explanation of why.
        result.stages_skipped["mmgbsa_failed_poses"] = (
            f"{len(failed)} of {len(posed)} rescored poses returned no MM/GBSA "
            f"energy and are ranked on docking alone: "
            + "; ".join(failed[:8])
            + (" ..." if len(failed) > 8 else ""))
        log.warning("MM/GBSA failed for %d of %d poses: %s",
                    len(failed), len(posed), "; ".join(failed[:8]))


def _dock_pool(candidates: list[Candidate], dock_top: Optional[int],
               result: "ScreenResult") -> list[Candidate]:
    'The candidates to dock, and a note on the report when that is not all.'
    dock_top = _valid_screen_count(dock_top, "dock_top", allow_none=True)
    if not dock_top or len(candidates) <= dock_top:
        for c in candidates:
            c.dock_attempted = True
        result.dock_scale = {"candidates": len(candidates),
                             "docked": len(candidates), "capped": False,
                             "reason": f"every one of the {len(candidates)} "
                                       "candidates that passed the gates was "
                                       "docked"}
        return candidates

    if getattr(result, "filter_profile", {}).get("dock_selection") == "balanced":
        return _balanced_dock_pool(candidates, dock_top, result)

    def _potency(c: Candidate) -> float:
        value = score._finite(c.ligand.activity.get("pchembl"))








        return value if value is not None else float("-inf")

    ordered = sorted(candidates, key=_potency, reverse=True)
    pool = ordered[:dock_top]
    # Say which side of the cut each candidate fell on, so `score` can tell
    # "docked and came back empty" from "never offered to the engine".
    chosen = {id(c) for c in pool}
    for c in candidates:
        c.dock_attempted = id(c) in chosen
    unmeasured = sum(1 for c in pool if _potency(c) == float("-inf"))
    cut = _potency(pool[-1])
    why = (
        f"{len(candidates)} candidates passed the gates and the {dock_top} with "
        f"the highest measured pChEMBL were docked (floor {cut:.2f})"
        + (f"; {unmeasured} of those carry no pChEMBL at all and were kept "
           "because unmeasured is not inactive" if unmeasured else "")
        + f". The other {len(candidates) - dock_top} have no pose, so they carry "
        "no docking or MM/GBSA evidence and cannot enter the structure-scored "
        "ranking. This is a potency-ordered slice of the gated library, not the "
        "whole of it -- raise --dock-top, or pass 0 to dock everything.")
    result.dock_scale = {"candidates": len(candidates), "docked": len(pool),
                         "capped": True, "pchembl_floor": cut,
                         "without_pchembl": unmeasured, "reason": why}
    result.stages_run.append(f"docking pool -- {why}")
    log.warning("docking %d of %d candidates (highest measured pChEMBL first)",
                len(pool), len(candidates))
    return pool


def _balanced_dock_pool(candidates: list[Candidate], limit: int,
                        result: ScreenResult) -> list[Candidate]:
    """Reserve exploration slots instead of treating missing potency as inactive."""
    from .discovery import priorities, scaffold
    values = priorities(candidates)
    def measured(candidate):
        value = candidate.ligand.activity.get("pchembl")
        return isinstance(value, (int, float)) and score._finite(value) is not None
    def priority(candidate):
        return (-values[candidate.ligand.ligand_id]["discovery_score"], candidate.ligand.ligand_id)
    known = sorted((c for c in candidates if measured(c)),
                   key=lambda c: (-c.ligand.activity["pchembl"], c.ligand.ligand_id))
    unknown = sorted((c for c in candidates if not measured(c)), key=priority)
    controls = sorted((c for c in candidates if c.ligand.source == "cosmetic_active"), key=priority)
    pool, chosen = [], set()
    def add(items, count):
        count = max(0, min(count, limit - len(pool)))
        added = 0
        for candidate in items:
            if added >= count:
                break
            cid = candidate.ligand.ligand_id
            if cid not in chosen:
                pool.append(candidate)
                chosen.add(cid)
                added += 1
    add(controls, limit)
    unknown_target = min(len(unknown), math.ceil(limit * .30))
    add(unknown, unknown_target - sum(not measured(c) for c in pool))
    known_target = min(len(known), math.ceil(limit * .50))
    add(known, known_target - sum(measured(c) for c in pool))
    # Round-robin chemical cores for the remaining slots, deterministically.
    from collections import defaultdict
    buckets = defaultdict(list)
    for candidate in sorted(candidates, key=priority):
        if candidate.ligand.ligand_id not in chosen:
            buckets[scaffold(candidate.ligand.smiles)].append(candidate)
    while len(pool) < limit and buckets:
        for core in sorted(list(buckets), key=lambda key: priority(buckets[key][0])):
            add(buckets[core][:1], 1)
            buckets[core].pop(0)
            if not buckets[core]:
                del buckets[core]
            if len(pool) == limit:
                break
    for candidate in candidates:
        candidate.dock_attempted = candidate.ligand.ligand_id in chosen
    n_unknown = sum(not measured(c) for c in pool)
    reason = (f"balanced-v1 docked {len(pool)}/{len(candidates)}: {n_unknown} unmeasured; "
              "curated controls reserved first, then 30% unmeasured and 50% measured "
              "targets, with remaining capacity assigned across chemical cores. "
              "Unfilled quotas are reallocated; unevaluated is not inactive.")
    result.dock_scale = {"policy": "balanced-v1", "candidates": len(candidates),
                         "docked": len(pool), "capped": True, "without_pchembl": n_unknown,
                         "known_fraction_target": .50, "unmeasured_fraction_target": .30,
                         "selected_ids": [c.ligand.ligand_id for c in pool], "reason": reason}
    result.stages_run.append("docking pool -- " + reason)
    return pool
















SA_MAKEABLE = 6.0


def _note_synthesis(result, kept) -> None:
    'Score how hard the unpurchasable compounds would be to make.'
    from .hop_sanity import sa_score










    answered, unchecked = [], []
    for lig in kept:
        if lig.props.get("purchasable"):
            continue
        status = str(lig.props.get("sourcing_status") or "unknown")
        (answered if status in ("known_not_sold", "novel") else unchecked).append(lig)
    if not answered and not unchecked:
        return

    n_makeable = 0
    # Scored for both groups -- an unchecked compound is as makeable as a
    # confirmed-absent one -- but the availability state is kept apart.
    for lig in answered + unchecked:
        sa = sa_score(lig.smiles)
        lig.props["sa_score"] = sa
        if sa is not None and sa <= SA_MAKEABLE:
            n_makeable += 1
    parts = [f"synthesis: of {len(kept)} compounds, {len(answered)} were "
             "checked and no vendor sells them"]
    if unchecked:
        parts.append(f"{len(unchecked)} could not be checked at all, so for "
                     "those 'not purchasable' means 'not asked'")
    parts.append("all of them stay in the ranking, which is ordered on the "
                 "science and not on what is in stock")
    parts.append(
        f"{n_makeable} of the {len(answered) + len(unchecked)} score "
        f"SA <= {SA_MAKEABLE:g} -- a rough triage, not a verdict: on this "
        "install paclitaxel scores 5.38 and erythromycin 6.09, so the line "
        "falls between two of the hardest syntheses ever completed. Read "
        "`sa_score` as complexity, and ask a chemist about a route")
    result.stages_run.append(". ".join(parts))


def _note_composite_scale(result) -> None:
    """One ranked list holding rows scored on different term sets is not one
    ranking.

    Module level so the branch itself is tested: as an inline `if` past the
    docking stage no offline test reached it, and replacing the condition with
    `if False` passed every test that existed.
    """
    if not result.breakdowns:
        return
    summary = score.structure_evidence_summary(result.breakdowns)
    if not summary.get("incomparable"):
        return
    result.gates_unverified["composite_scale"] = (
        f"{summary['n_unevaluated']} of {summary['n']} candidates were never "
        "docked. Their composite is a mean over fewer terms than a docked "
        "row's, so the ranked column mixes two scales: an undocked compound is "
        "not competing with a docked one, it is being scored on a different "
        "question. Read "
        "ranks within a `structure_evidence` group; raising --dock-top to "
        "cover every gated candidate is what removes this.")


def _note_docking_absent(result, reason: str) -> None:
    'Docking was asked for and did not happen. Say it where trust is read.'
    result.stages_skipped["docking"] = reason
    result.gates_unverified["docking"] = (
        f"{reason}. Docking was requested and did not run, so no candidate in "
        "this screen carries a pose, a docking score or an MM/GBSA number, and "
        "the structure term is absent from every composite. This is a gated "
        "screen wearing a docked screen's filenames; re-run it before reading "
        "the ranking as structure-informed.")


class _MMGBSARefused(Exception):
    """The rescore was refused for a stated reason, not lost to an error.

    Carried as an exception only so it can leave the middle of the batch
    block; the reason is already on the result by the time it is raised, and
    the handler must not overwrite it with a `type(e).__name__` string.
    """


def _note_provenance(result, node: str, audit: dict, ligands) -> None:
    """Say what provenance QC actually did to this library.

    Module level so the test can call the real thing. The previous
    wording lived inline in `run_screen`, could only be reached through a
    network screen, and was therefore never asserted on.
    """






























    contamination = str(audit.get("warning") or "").strip()
    frac = audit.get("trustworthy_fraction")
    counts = audit.get("counts") or {}
    n_raw = sum(counts.values())
    n_kept = audit.get("n_trustworthy", counts.get("human_cell", 0) + counts.get("human_protein", 0))

    mix: dict[str, int] = {}
    human_target_hosts = 0
    for lig in ligands:
        act = lig.activity or {}
        if getattr(lig, "source", "") == "cosmetic_active":






            cls = {"nonhuman_cell": "nonhuman_cell",
                   "animal": "animal"}.get(str(act.get("evidence") or ""),
                                           "curated_active")
        else:
            cls = str(act.get("provenance") or "unrecorded")
            if (act.get("provenance_classifier_version") == provenance.CLASSIFIER_VERSION
                    and act.get("human_target_activity_eligible")
                    and cls == "nonhuman_cell"):
                human_target_hosts += 1
                cls = "human_target_nonhuman_system"
        mix[cls] = mix.get(cls, 0) + 1
    if human_target_hosts:
        result.gates_unverified["provenance_context"] = (
            f"{human_target_hosts} compound(s) for {node} have human-target activity "
            "in a nonhuman expression system. Target species and cell species are "
            "separate; these measurements do not verify human skin function or binding.")
    survived_nonhuman = {k: v for k, v in mix.items()
                         if k in ("surrogate", "nonhuman_cell", "animal",
                                  "unknown_origin", "unspecified") and v}

    if survived_nonhuman:




        what = ", ".join(f"{k} {v}" for k, v in
                         sorted(survived_nonhuman.items(), key=lambda kv: -kv[1]))
        result.gates_unverified["provenance"] = (
            f"{sum(survived_nonhuman.values())} of the {len(ligands)} "
            f"compound(s) ranked for {node} rest on a measurement that is not "
            f"confirmed human ({what}). Provenance QC did not remove them, so "
            "this ranking mixes evidence about the human protein with evidence "
            "about something else; `evidence` names the class per row and "
            "`score.EVIDENCE_FACTOR` discounts it, but the pChEMBL column does "
            "not. Of the activities fetched, "
            f"{n_kept} of {n_raw} were confirmed human.")
    elif contamination:
        pct = f"{(1 - frac) * 100:.0f}%" if isinstance(frac, (int, float)) else "some"
        why = ", ".join(f"{k} {v}" for k, v in
                        sorted(counts.items(), key=lambda kv: -kv[1])
                        if k not in ("human_cell", "human_protein") and v)
        result.gates_unverified["provenance"] = (
            f"provenance QC discarded {pct} of the activities fetched for "
            f"{node} -- {n_raw - n_kept} of {n_raw}"
            + (f" ({why})" if why else "") + ". The ranking below rests on the "
            f"{n_kept} that are confirmed human measurements, so it is not "
            "contaminated; what is unverified is whether that surviving "
            "fraction is representative."
            # Only say this where it is true. `unspecified` means the assay
            # record was blank, which is not evidence that the row was
            # non-human -- it is the absence of evidence either way, and on
            # HIF1A it is 4072 of 4515 rows. A node that dropped nothing for
            # that reason should not be handed the sentence.
            + (f" {counts['unspecified']} of them were dropped as "
               "`unspecified` -- a blank assay record, not a non-human one -- "
               "and any of those could have changed the order."
               if counts.get("unspecified") else ""))

    # The survivors are all human but not all on one scale: a cell IC50 and a
    # biochemical Ki answer different questions. `score.EVIDENCE_FACTOR`
    # already discounts the biochemical class (0.80 against 0.90), so the
    # composite ordering carries a correction -- the raw `pchembl` column a
    # reader quotes does not. Stated rather than raised as a gate, because a
    # correction does exist and a gate that fires on every mixed library would
    # carry no information.
    ev_mix = {k: v for k, v in mix.items() if k in ("human_cell", "human_protein")}
    if len(ev_mix) > 1:
        result.stages_run.append(
            "evidence mix: "
            + ", ".join(f"{v} {k}" for k, v in sorted(ev_mix.items()))
            + " -- pChEMBL from cell and biochemical assays are not on one "
              "scale; the composite applies an evidence discount, the pChEMBL "
              "column does not")


def run_screen(marker: str, workdir: Path, *,
               registry: Optional[Registry] = None,
               expected_node: Optional[str] = None,
               extra_files: Optional[Iterable[Path]] = None,
               include_actives: bool = True,
               pchembl_min: Optional[float] = None,
               max_chembl: int = 3000,
               trustworthy_only: bool = True,
               min_topical: float = 4.0,
               direction_gate: Optional[str] = None,
               admet_gate: Optional[str] = "critical",
               applied_percent: float = 1.0,
               product: str = "body_lotion",
               check_sourcing: bool = True,
               purchasable_only: bool = False,
               exclude_patented: bool = False,
               dock: bool = True,




               final_top: Optional[int] = 300,






               dock_top: Optional[int] = None,








               engine: str = "vina",
               #: Where the `unidock` binary is, when it is not on PATH. It
               #: ships as a separate conda package with a CUDA build, so on
               #: most machines it lives in its own environment.
               unidock_bin: Optional[str] = None,
               gbsa: bool = False,
               gbsa_top_n: int = 20,
               pdb_id: Optional[str] = None,





               ligand_code: Optional[str] = None,
               exhaustiveness: int = 16,
               #: OpenMM platform for the MM/GBSA rescoring: "CUDA", "OpenCL",
               #: "CPU", "Reference". None leaves `make_platform`'s default,
               #: which is CPU. Changing it can materially change energies,
               #: so a table must not mix two platforms.
               platform: Optional[str] = None,
               weights: Optional[dict[str, float]] = None,
               dock_selection: str = "balanced",
               gbsa_cpu_profile: str = "auto",
               gbsa_threads: Optional[int] = None,
               gbsa_charge_method: str = "am1bcc-fixed",
               gbsa_ligand_forcefield: str = "openff-2.2.0") -> ScreenResult:
    # Both gates are closed vocabularies, and both used to fall through to a
    # default that ran a *different* screen than the caller asked for. They are
    # settled here, before a directory is made or a byte is fetched: a screen
    # that cannot be configured should cost nothing and write nothing.
    # `admet_gate="none"` is the CLI's documented "runs no ADMET gate at all"
    # and resolves to None; it used to resolve to index 5 of SEVERITY_ORDER and
    # delete the entire library.
    score.validate_weights(weights)
    min_topical = validate_topical_cutoff(min_topical)
    _valid_screen_count(max_chembl, "max_chembl")
    _valid_screen_count(dock_top, "dock_top", allow_none=True)
    _valid_screen_count(final_top, "final_top", allow_none=True)
    _valid_screen_count(gbsa_top_n, "gbsa_top_n", minimum=1)
    _valid_screen_count(exhaustiveness, "exhaustiveness", minimum=1)
    if pchembl_min is not None:
        pchembl_min = _valid_floor(pchembl_min, "the pchembl_min argument")
    from .filters.admet_exposure import PRODUCT_USE, validate_applied_percent
    if isinstance(applied_percent, bool):
        raise ValueError("applied_percent must be a finite number between 0 and 100")
    applied_percent = validate_applied_percent(applied_percent)
    if not isinstance(product, str) or product not in PRODUCT_USE:
        raise ValueError(f"product must be one of {sorted(PRODUCT_USE)}")
    if engine not in ("vina", "unidock", "hybrid"):
        raise ValueError("engine must be vina, unidock or hybrid")
    if dock_selection not in ("balanced", "potency"):
        raise ValueError("dock_selection must be balanced or potency")
    admet_gate = normalise_severity_gate(admet_gate, field="admet_gate")
    direction_gate = validate_direction_gate(direction_gate)
    if expected_node is not None and (not isinstance(expected_node, str) or not expected_node.strip()):
        raise ValueError("expected_node must be a non-empty gene name or None")
    # Validate the request before any fetch or output, even when rescoring is
    # disabled: invalid settings must not quietly configure a different run.
    from .mmgbsa.engine import cpu_threads
    cpu_threads(gbsa_threads, gbsa_cpu_profile)
    from .mmgbsa.ligand import validate_screen_protocol
    gbsa_charge_method, gbsa_ligand_forcefield = validate_screen_protocol(
        gbsa_charge_method, gbsa_ligand_forcefield)
    workdir = Path(workdir)
    from .publication import ensure_mutable_workdir
    ensure_mutable_workdir(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    reg = registry or Registry()
    spec: TargetSpec = reg.spec(marker)
    tri = triage.assess(spec)
    triage.require_screening_metadata(spec, tri.screening_node)
    if expected_node is not None and (tri.screening_node is None
            or tri.screening_node.upper() != expected_node.strip().upper()):
        raise ValueError(f"planned screening node {expected_node!r} differs from "
                         f"resolved node {tri.screening_node!r}; regroup before acquisition")
    result = ScreenResult(marker=spec.marker, triage=tri)
    result.decision = strategy.choose_track(spec, tri)

    if spec.mode is TargetMode.ENDPOINT:
        result.stages_skipped["all"] = Registry.explain(spec)
        return result

    if tri.counts_unresolved:
        result.gates_unverified["triage"] = tri.recommendation
    node = tri.screening_node
    if not node:
        result.stages_skipped["all"] = "no druggable node registered for this marker"
        return result

    selected = next(n for n in spec.nodes if n.gene == node)
    result.target_design = dict(selected.goal_support) | {
        'marker': spec.marker, 'node_gene': node,
        'required_action': selected.node_action,
        'marker_direction': spec.direction.value,
        'status': selected.goal_support.get('status', 'review_required')}















    floor, floor_why = _potency_floor(reg, spec, node, pchembl_min)
    # Kept under a second name because the topical-drop loop below
    # rebinds `floor` to a *formatted string* for its reason text, and a
    # later numeric comparison against it then fails on str vs float.
    library_floor = floor
    result.potency_floor = {"pchembl_min": floor, "source": floor_why}
    activity_label = "all recorded endpoints" if floor == 0 else f"pChEMBL >= {floor:g}"

    # Asked before the fetch, because the alternative is finding out after the
    # docking. It informs and does not abort: a capped library is still a
    # library, and which trade-off to take is the caller's.
    scale = check_library_scale(
        node, floor, max_chembl, reg,
        include_without_pchembl=floor == 0,
        target_chembl_id=next((n.chembl_target for n in spec.nodes
                               if n.gene == node and n.chembl_target), None))
    result.library_scale = scale
    # Zero retains its full-query meaning even when the advisory count fails
    # or changes between count and acquisition. The activity client measures
    # its own live coverage; an earlier count must never become a hidden cap.
    if max_chembl == 0:
        available = scale.get("available")
        scale["auto_max_chembl"] = 0
        scale["acquisition_budget_policy"] = "all_live_pages"
        if isinstance(available, int) and not isinstance(available, bool) and available >= 0:
            log.info("library scale: advisory count is %d; retrieving the full live query", available)
        else:
            result.gates_unverified["library_scale"] = (
                "the advisory activity count could not be measured; --max-chembl 0 "
                "still requests all live pages. Coverage must be established by the "
                "acquisition receipt. " + str(scale.get("reason", "")))
            log.warning("library scale: count unavailable; full live acquisition remains uncapped")
    elif scale.get("capped"):
        result.gates_unverified["library_scale"] = scale["reason"]
        log.warning("library scale: %s", scale["reason"])
    elif not scale.get("checked"):
        log.info("library scale: %s", scale.get("reason", "not checked"))

    # ---- library ---------------------------------------------------------
    identity = selected.goal_support.get("target_identity")
    identity_args = {"target_identity": identity} if identity is not None else {}
    ligands, audit = library.assemble(
        node_gene=node, marker=spec.marker, extra_files=extra_files,
        include_actives=include_actives, pchembl_min=floor,
        max_chembl=max_chembl, trustworthy_only=trustworthy_only,
        preserve_all_queried=True, include_without_pchembl=floor == 0, **identity_args)
    audit["pchembl_min"] = floor
    audit["pchembl_min_source"] = floor_why
    result.provenance_audit = audit
    result.stages_run.append(
        f"library ({len(ligands)} compounds; {activity_label})")
    # A screen that quietly proceeds on a fraction of its library reports a
    # confident answer about nothing. The ratio test alone cannot see this: a
    # truncated activity pull shrinks its own denominator, so the proportion
    # stays healthy while the data vanishes. Ask the servers what they had.
    short: list[str] = []
    fetch = audit.get("fetch") or {}
    mols = audit.get("molecules") or {}
    if fetch.get("assay_context_complete") is False:
        result.gates_unverified["assay_context"] = (
            f"Assay context is incomplete: {fetch.get('n_unresolved_assays', 0)} unresolved assay(s), "
            f"{len(fetch.get('stale_assay_ids') or [])} stale assay context(s), and "
            f"{fetch.get('n_activities_without_assay_id', 0)} activity row(s) without an assay identity. "
            "Recorded measurements are retained; failed context acquisition does not establish absent human activity. "
            "Retry when assay metadata is available before treating this provenance assessment as complete.")
    partial_inputs = [receipt for receipt in audit.get("input_files", [])
                      if receipt.get("complete") is not True]
    if partial_inputs:
        result.gates_unverified["input_files"] = (
            f"{len(partial_inputs)} supplied file(s) were not imported completely: "
            + "; ".join(f"{receipt.get('path')}: {receipt.get('accepted', 0)} accepted, "
                        f"{len(receipt.get('rejected_rows') or [])} rejected, "
                        f"status {receipt.get('status', 'UNRECORDED')}"
                        for receipt in partial_inputs)
            + ". Inspect provenance_audit.input_files; original supplied files were not changed.")
    if fetch.get("complete") is False:
        short.append(f"ChEMBL offered {fetch.get('total_count')} activities and only "
                     f"{fetch.get('retrieved')} were retrieved")
    # `complete` cannot see this one. `fetch_activities` sets
    # `expected = min(total, max_records)`, so a pull cut off by the cap reports
    # complete=True -- the denominator moved with the numerator. Only
    # `capped_at_max_records` says otherwise, and nothing here read it: a sweep
    # wrote IL-31 as a finished screen assembled from 3000 of 14919 activities
    # with `stages_skipped` empty. ChEMBL returns activities in database order,
    # not by potency, so the slice is not even the top of the range.
    if fetch.get("capped_at_max_records"):
        short.append(
            f"ChEMBL holds {fetch.get('total_count')} activities ({activity_label}) "
            f"and only the first {fetch.get('retrieved')} were pulled "
            "(max_chembl); this library is a database-order slice, not the "
            "target's chemistry. Raise max_chembl to cover this query")
    if mols.get("complete") is False:
        short.append(f"{mols.get('n_missing')} of {mols.get('requested')} molecule "
                     "structures were never returned")










    n_mols = mols.get("retrieved") or 0
    if n_mols and len(ligands) < 0.1 * n_mols:
        short.append(f"only {len(ligands)} of {n_mols} retrieved structures "
                     "survived into the library")
    # A cached answer served because the live fetch failed is not a fresh one.
    if fetch.get("stale") or mols.get("stale"):
        result.gates_unverified["chembl_cache"] = (
            f"ChEMBL did not answer; this library was served from a cache entry up "
            f"to {fetch.get('cache_age_days') or 0:.0f} days old. Re-run when the "
            "API is back before treating it as current.")
    if short:
        result.stages_skipped["library_truncated"] = (
            "; ".join(short) + ". Re-run before trusting this ranking.")
        result.gates_unverified["library"] = result.stages_skipped["library_truncated"]
        log.warning("library looks truncated: %s", "; ".join(short))
























    warning = str(audit.get("library_unavailable") or "").strip()
    if warning:
        n_other = len(ligands)
        if warning[-1:] not in ".!?":
            warning += "."
        detail = (f"{warning} The ChEMBL half of the library for {node} is "
                  f"therefore absent, and the {n_other} compound(s) screened "
                  "below come from the curated cosmetic actives and any "
                  "supplied files alone -- not from a search that came back "
                  "empty.")
        result.stages_skipped["library"] = detail
        result.gates_unverified["library"] = (
            "; also: ".join(x for x in
                            (result.gates_unverified.get("library"), detail) if x))
        log.warning("library not retrieved for %s: %s", node, warning)
    if not ligands:
        _note_provenance(result, node, audit, ligands)
        result.stages_skipped["screen"] = "library is empty after provenance QC"
        return result

    # Preserve the entire acquired structure pool before any presentation gate.
    # Assessments remain available when a researcher disables a gate in the UI.
    assessments = {ligand_record.ligand_id: assess_skin(ligand_record, node) for ligand_record in ligands}
    prepared_admet = {}
    admet_profiles = {}
    for ligand in ligands:
        try:
            profile = admet_f.summarise(ligand.smiles, applied_percent=applied_percent,
                                       product=product, ligand_id=ligand.ligand_id,
                                       target_gene=node)
            prepared_admet[ligand.ligand_id] = profile
            admet_profiles[ligand.ligand_id] = profile.to_dict()
        except Exception as exc:
            ligand.flags.append(f"ADMET not assessed: {type(exc).__name__}")
    result.all_candidates = [Candidate(ligand=ligand_record, target=node, skin=assessments[ligand_record.ligand_id],
                                      admet=admet_profiles.get(ligand_record.ligand_id), evidence=evidence_for(ligand_record))
                             for ligand_record in ligands]
    result.filter_profile = {"min_topical": min_topical, "admet_gate": admet_gate,
                             "direction_gate": direction_gate, "product": product,
                             "pchembl_min": floor, "display_limit": 1000,
                             "population": "all acquired valid structures before gates",
                             "provenance_gate": "human_target_eligible" if trustworthy_only else "none",
                             "activity_scope": "all_recorded" if floor == 0 else "pchembl_floor",
                             "dock_selection": dock_selection}

    # A preserved organic salt/mixture must not be scored as though all of its
    # fragments formed one active parent. Keep it in the discovery population,
    # and require an explicit structure choice before computational selection.
    single_component = []
    parent_review = []
    for ligand in ligands:
        identity_note = ligand.props.get("structure_identity")
        if not isinstance(identity_note, dict):
            identity_note = library.structure_identity(ligand.smiles, ligand.smiles,
                                                        source_id=ligand.ligand_id)
            ligand.props["structure_identity"] = identity_note
        if identity_note.get("requires_parent_review") is True:
            parent_review.append(ligand.ligand_id)
            result.dropped.append({"ligand_id": ligand.ligand_id,
                "name": ligand.name or ligand.ligand_id, "stage": "structure_identity",
                "severity": "unverified", "smiles": ligand.smiles,
                "reason": "Unresolved multi-component structure: an active parent was not selected. "
                          "Exact source structure and parent identifiers remain recorded; explicit parent review is required before ranking or docking."})
        else:
            single_component.append(ligand)
    if parent_review:
        result.gates_unverified["structure_identity"] = (
            f"{len(parent_review)} multi-component input(s) require explicit parent review; "
            "they remain in the discovery pool and are excluded from computational selection.")
        audit["structure_identity_holds"] = parent_review
    ligands = single_component
    if not ligands:
        result.stages_skipped["screen"] = "no input has a resolved single-component structure; explicit parent review required"
        return result

    # QC controls the computational ranking, never the saved discovery pool.
    # A merged structure qualifies if any of its recorded rows qualifies;
    # curated controls and user libraries keep their distinct evidence basis.
    before_provenance = len(ligands)
    eligible = []
    for ligand in ligands:
        activity = ligand.activity
        measured = activity.get("measurements") or []
        allowed = (not trustworthy_only or ligand.source not in ("chembl", "chembl_local")
                   or not activity.get("provenance_gate_deferred")
                   or any(row.get("provenance_trustworthy") is True for row in measured))
        if activity.get("provenance_gate_deferred"):
            activity["provenance_policy_eligible"] = allowed
        if allowed:
            eligible.append(ligand)
        else:
            result.dropped.append({"ligand_id": ligand.ligand_id,
                "name": ligand.name or ligand.ligand_id, "stage": "provenance",
                "severity": "unverified", "smiles": ligand.smiles,
                "reason": "No recorded activity meets the configured human-target provenance policy; "
                          "structure and measurements remain in the discovery population."})
    ligands = eligible
    result.stages_run.append(f"provenance filter ({len(ligands)}/{before_provenance} compounds pass)")
    _note_provenance(result, node, audit, ligands)
    if not ligands:
        result.stages_skipped["screen"] = "no compound passed provenance QC; all acquired structures remain saved"
        return result

    # ---- direction -------------------------------------------------------
    # Before anything expensive: a compound that moves the node the wrong way
    # is not a weak candidate, it is the wrong candidate, and no permeation or
    # docking number redeems it.
    ligands = _direction_gate(result, ligands, node,
                              _required_action(reg, spec, node),
                              (audit.get("mechanism") or {}), direction_gate)
    if not ligands:
        result.stages_skipped["screen"] = "every compound opposed the required direction"
        return result

    # ---- skin filters ----------------------------------------------------
    kept: list[LigandRecord] = []
    for ligand_record in ligands:
        a = assessments[ligand_record.ligand_id]
        if a.topical_score >= min_topical or ligand_record.source == "cosmetic_active":
            kept.append(ligand_record)
            continue







        shown, floor = _fmt_below(a.topical_score, min_topical)
        result.dropped.append({
            "ligand_id": ligand_record.ligand_id, "name": ligand_record.name or ligand_record.ligand_id,
            "stage": "topical", "severity": a.permeation_class,
            "reason": (f"topical score {shown} < {floor} "
                       f"(permeation {a.permeation_class}"
                       + (f", log Kp {a.log_kp:.2f}" if a.log_kp is not None else "")
                       + ")"),
            # The raw values as well as the sentence. The report had to parse
                # the permeation class back out of the prose to group these,
                # and it had to repair "topical score 6.0 < 6.0" on the
                # rendering side because the rounded score was the only number
                # it was given. A consumer should not have to read a sentence to
                # recover a number the producer already held.
                "topical_score": a.topical_score,
                "topical_floor": min_topical,
                "permeation": a.permeation_class,
                "log_kp": a.log_kp,
                "smiles": ligand_record.smiles})
    result.stages_run.append(
        f"skin filter ({len(kept)}/{len(ligands)} pass topical >= {min_topical})")

    # re-grade now that we know how many are actually deliverable
    n_feasible = sum(1 for ligand_record in kept
                     if (ligand_record.activity.get("pchembl") or 0) >= triage.FEASIBLE_PCHEMBL
                     and assessments[ligand_record.ligand_id].topical_score >= triage.FEASIBLE_TOPICAL)
    result.triage = triage.assess(spec, n_feasible=n_feasible)
    result.decision = strategy.choose_track(spec, result.triage)



















    if library_floor >= triage.FEASIBLE_PCHEMBL:
        result.gates_unverified["grade_comparability"] = (
            f"this library was fetched at pChEMBL >= {library_floor:g}, at or above the "
            f"{triage.FEASIBLE_PCHEMBL:g} the grade counts against, so every "
            f"compound in it clears the potency half by construction and "
            f"n_feasible={n_feasible} reflects the topical gate alone. The "
            "grade stands for this target and must not be ranked against a "
            "target screened from a lower floor.")

    # Stop where the screen actually ended. Letting an empty `kept` fall
    # through to the sourcing block meant the summary named sourcing for a
    # wipeout the skin filter caused -- and sourcing then "ran" on nothing.
    if not kept:
        result.stages_skipped["screen"] = _emptied_by(
            result, f"no compound passed the skin filters "
                    f"(topical >= {min_topical})")
        return result

    # ---- sourcing ---------------------------------------------------------
    # Runs before docking, because a compound nobody sells cannot be assayed and
    # there is no point spending CPU posing it.
    sourcing_data = False           # did vendor/patent annotation actually land?
    n_buy = n_pat = 0








    n_sourced = 0
    if check_sourcing:
        try:
            sourcing.annotate(kept, progress=len(kept) > 50)
            n_sourced = len(kept)
            n_buy = sum(1 for ligand_record in kept if ligand_record.props.get("purchasable"))
            n_pat = sum(1 for ligand_record in kept if ligand_record.props.get("in_patents"))
            # A UniChem/PubChem timeout on one compound leaves status "unknown"
            # and purchasable False, which is indistinguishable from a vendor
            # search that came back empty.
            n_unk = sum(1 for ligand_record in kept
                        if ligand_record.props.get("sourcing_status", "unknown") == "unknown")









            result.stages_run.append(
                f"sourcing ({n_buy}/{len(kept)} purchasable, {n_pat} in patent literature"
                + (f", {n_unk} unresolved)" if n_unk else ")"))
            _note_synthesis(result, kept)
            if n_unk:
                result.gates_unverified["sourcing"] = (
                    f"{n_unk} of {len(kept)} vendor lookups did not resolve; for those "
                    "compounds 'not purchasable' means 'not checked'")
            sourcing_data = True
        except Exception as e:
            result.stages_skipped["sourcing"] = f"{type(e).__name__}: {e}"
            result.gates_unverified["sourcing"] = _NO_SOURCING % f"failed: {e}"
            log.warning("sourcing skipped: %s", e)
    else:
        result.stages_skipped["sourcing"] = "disabled by the caller (check_sourcing=False)"
        result.gates_unverified["sourcing"] = _NO_SOURCING % "was not run"

    # The two filters live outside the lookup, because they are a separate
    # decision from it and used to be buried inside its `try`: with
    # `check_sourcing=False` -- or with `sourcing.annotate` raising --
    # `--purchasable-only` and `--exclude-patented` were silently not applied
    # and no stage line said so, so `screen MARKER --no-sourcing
    # --purchasable-only` printed a full candidate list captioned as though
    # both flags had been honoured. The lookup's skip and the filter's skip are
    # different facts and are recorded separately.
    if sourcing_data:
        if purchasable_only:
            kept = _filter_out(
                result, kept, keep=lambda ligand_record: bool(ligand_record.props.get("purchasable")),
                stage="sourcing", severity="not purchasable",
                reason=lambda ligand_record: ("no vendor found by the sourcing lookup "
                                  f"(status {ligand_record.props.get('sourcing_status', 'unknown')}); "
                                  "removed by purchasable_only"))
            result.stages_run.append(f"purchasable-only filter ({len(kept)} remain)")
        if exclude_patented:
            kept = _filter_out(
                result, kept, keep=lambda ligand_record: not ligand_record.props.get("in_patents"),
                stage="sourcing", severity="patented",
                reason=lambda ligand_record: ("appears in the patent literature "
                                  f"({ligand_record.props.get('n_patents') or 'n/a'} records); "
                                  "removed by exclude_patented"))
            result.stages_run.append(f"patent-free filter ({len(kept)} remain)")
        denom = max(n_sourced, 1)
        result.decision = strategy.choose_track(
            spec, result.triage,
            patent_fraction=n_pat / denom,
            purchasable_fraction=n_buy / denom)
    elif purchasable_only or exclude_patented:
        asked = " and ".join(
            flag for flag, on in (("purchasable_only", purchasable_only),
                                  ("exclude_patented", exclude_patented)) if on)
        why = result.stages_skipped.get("sourcing", "the sourcing stage did not run")
        msg = (f"{asked} was requested, but no vendor or patent data was "
               f"retrieved ({why}), so the filter was not applied. Every "
               f"compound below is unfiltered on it -- this is not a screen in "
               "which all of them are purchasable or patent-free.")
        result.stages_skipped["sourcing_filters"] = msg
        result.gates_unverified["sourcing_filters"] = msg
        log.warning("sourcing filters not applied: %s", msg)
    if not kept:
        result.stages_skipped["screen"] = _emptied_by(
            result, "no compound survived the sourcing filters")
        return result

    # ---- ADMET and skin toxicity -----------------------------------------
    # Cheap, structural, and decisive: a furocoumarin is not going into a
    # leave-on product whatever it docks at. Runs before the expensive stages.
    if admet_gate is not None:
        order = SEVERITY_ORDER
        # `admet_gate` came through `normalise_severity_gate`, so it is one of
        # GATE_SEVERITIES here and the index is real. The `else 0` this
        # replaced is what turned a misspelling into a different gate and
        # "none" into a gate that dropped the whole library.
        cut = order.index(admet_gate)
        dropped: list[str] = []
        unassessed: list[str] = []
        n_entered = len(kept)
        for ligand_record in kept:
            try:
                # The screen knows which protein it is about; the ADMET layer
                # cannot judge a target-activated hazard without being told.
                # A 4-substituted phenol is a footnote in most screens and the
                # central risk in a tyrosinase one, because tyrosinase is the
                # enzyme that makes the melanocyte-killing quinone. The rule
                # existed and the call was not passing the target, so it never
                # fired in a real run.
                prof = prepared_admet[ligand_record.ligand_id]
            except Exception as e:                    # never let one molecule stop a screen
                log.warning("ADMET failed for %s: %s", ligand_record.ligand_id, e)
                unassessed.append(ligand_record.ligand_id)
                ligand_record.flags.append(f"ADMET not assessed: {type(e).__name__}")
                continue
            if not prof.assessed:                     # e.g. the structure would not parse
                unassessed.append(ligand_record.ligand_id)
                ligand_record.flags.append(f"ADMET not assessed: {prof.overall[:80]}")
            admet_profiles[ligand_record.ligand_id] = prof.to_dict()
            sev = hazard_severity_of(prof)
            if sev in order and order.index(sev) <= cut:
                dropped.append(ligand_record.ligand_id)
                ligand_record.flags.append(f"ADMET {sev}: {prof.overall[:120]}")





                reasons = sorted({c.name for c in prof.concerns
                                  if c.severity == sev and c.domain != "systemic"})
                result.dropped.append({
                    "ligand_id": ligand_record.ligand_id, "name": ligand_record.name or ligand_record.ligand_id,
                    "stage": "admet", "severity": sev,
                    "reason": ", ".join(reasons) or prof.overall[:90],
                    "smiles": ligand_record.smiles})
        if dropped:
            kept = [ligand_record for ligand_record in kept if ligand_record.ligand_id not in dropped]
        result.stages_run.append(
            f"ADMET gate at '{admet_gate}' ({len(dropped)} dropped, {len(kept)} remain"
            + (f", {len(unassessed)} NOT ASSESSED)" if unassessed else ")"))
        if unassessed:
            # Without this the run reads "0 dropped", indistinguishable from a
            # clean library, and an unevaluated furocoumarin reaches the buy list.
            result.stages_skipped["admet_incomplete"] = (
                f"{len(unassessed)} of {n_entered} compounds could not be "
                f"assessed, so the gate never saw them: "
                f"{', '.join(unassessed[:8])}{' ...' if len(unassessed) > 8 else ''}")
            result.gates_unverified["admet"] = result.stages_skipped["admet_incomplete"]
    else:
        msg = ("the ADMET gate was disabled (admet_gate=None or 'none'); no "
               "compound in this screen was excluded on hazard grounds and "
               "hazard profiles remain available as annotations")
        result.stages_skipped["admet_gate"] = msg
        result.stages_run.append(f"ADMET annotation ({len(admet_profiles)} profiles; gate disabled)")
    if not kept:
        result.stages_skipped["screen"] = _emptied_by(
            result, "every compound failed the ADMET gate")
        return result

    population_by_id = {c.ligand.ligand_id: c for c in result.all_candidates}
    candidates = [population_by_id[ligand_record.ligand_id] for ligand_record in kept]
    for candidate in candidates:
        candidate.direction = _direction_of(candidate.ligand)

    # ---- structure-based stages (optional) -------------------------------
    # Set the moment `dock_many` is actually dispatched. `stages_skipped["docking"]`
    # cannot answer this: it carries both "the docking modules are unavailable",
    # where nothing was attempted, and "dock_many raised", where everything was.
    # Scoring the second as `not_run` excuses the structure term for a screen
    # whose docking blew up, which is the failure this flag exists to prevent.
    dock_dispatched = False
    if dock:
        try:
            from .structure import prepare as prep_mod
            from .dock import vina as vina_mod
        except ImportError as e:
            _note_docking_absent(
                result, f"docking modules unavailable: {e}")
        else:
            try:












                pid, receptor = _prepare_docking_receptor(
                    prep_mod, spec, node, pdb_id, workdir / "receptor",
                    ligand_code=ligand_code)
                result.receptor = receptor

                # `vina.py` calls this "the only cheap check that the receptor
                # preparation and the box are sane", and until now nothing
                # outside the test suite ran it: a box placed on the wrong
                # chain or off the pocket entirely produced a full ranking with
                # no gate objecting. It runs here, before the library is posed,
                # because a bad box makes every pose after it worthless and the
                # control costs one docking job to find that out.
                #
                # It informs, it does not abort: on a metalloenzyme a failing
                # redock can mean the scoring function, not the geometry (see
                # `redock_validation`'s own docstring), and that is a judgement
                # for the reader. A control that could not run at all is
                # likewise recorded rather than raised.
                try:
                    rd = vina_mod.redock_validation(
                        receptor, workdir / "redock",
                        exhaustiveness=exhaustiveness,
                        engine=engine, unidock_bin=unidock_bin)
                except Exception as e:                       # noqa: BLE001
                    result.gates_unverified["redock"] = (
                        f"the redock control could not run ({type(e).__name__}: "
                        f"{e}), so nothing here confirms the box sits on the "
                        "pocket")
                else:
                    _note_redock_outcome(result, rd, vina_mod.REDOCK_RMSD_PASS)
                    _note_scoring_unreliable(result, receptor, rd)







                amb = next((n for n in (receptor.notes or [])
                            if n.startswith("BOX LIGAND AMBIGUOUS")), None)
                if amb:
                    result.gates_unverified["box_ligand"] = amb




                split = next((n for n in (receptor.notes or [])
                              if n.startswith(("BOX POCKET SPLIT", "BOX POCKET UNASSESSED"))), None)
                if split:
                    result.gates_unverified["box_pocket"] = split




                _note_receptor_preparation(result, receptor)

                # A pin that resolved to nothing. The box then came from
                # burial or from the chain centroid -- whichever the next tier
                # is -- and the caller asked for neither.
                missing = next((n for n in (receptor.notes or [])
                                if n.startswith("BOX LIGAND MISSING")), None)
                if missing:
                    result.gates_unverified["box_ligand_pin"] = missing

                dock_dispatched = True
                dock_pool = _dock_pool(candidates, dock_top, result)
                df = vina_mod.dock_many(receptor, [c.ligand for c in dock_pool],
                                        workdir / "dock",
                                        exhaustiveness=exhaustiveness,
                                        engine=engine,
                                        unidock_bin=unidock_bin)
                result.dock_scale["engine"] = engine
                execution_by_ligand = df.attrs.get("execution_by_ligand", {})
                result.dock_scale["actual_engine_counts"] = {
                    backend: sum(trace.get("actual_engine") == backend and trace.get("state") == "completed"
                                 for trace in execution_by_ligand.values())
                    for backend in ("vina", "unidock")}
                result.dock_scale["execution_receipt_count"] = len(execution_by_ligand)
                best = {}
                for row in df.to_dict("records"):
                    lid = row["ligand_id"]
                    if lid not in best or row["affinity"] < best[lid]["affinity"]:
                        best[lid] = row
                from .types import DockPose
                for c in candidates:
                    c.docking_execution = execution_by_ligand.get(c.ligand.ligand_id)
                    r = best.get(c.ligand.ligand_id)
                    if r:
                        c.pose = DockPose(ligand_id=r["ligand_id"], pose_index=r.get("pose_index", 0),
                                          affinity=r["affinity"], sdf_path=r.get("sdf_path"),
                                          ligand_efficiency=r.get("ligand_efficiency"))
                result.stages_run.append(f"docking into {pid} ({len(best)} ligands posed)")
            except Exception as e:
                _note_docking_absent(result, f"{type(e).__name__}: {e}")
                log.warning("docking skipped: %s", e)




    from .dock.vina import REDOCK_RMSD_PASS
    from .structure_policy import weights_for_run
    caller_set_structure = "structure" in (weights or {})
    scoring_weights, result.structure_validation = weights_for_run(
        weights, explicit=caller_set_structure, docked=dock_dispatched,
        redock=result.redock or {}, threshold=REDOCK_RMSD_PASS)
    if dock_dispatched:
        policy = result.structure_validation
        result.gates_unverified["structure_weight"] = policy["reason"]
        result.stages_run.append(
            f"structure weight {policy['applied_structure_weight']:g} for {node} "
            f"({policy['application']})")
        log.info("%s", policy["reason"])

    if gbsa and result.receptor is not None:
        try:
            from .mmgbsa import openmm_gbsa as gb
        except ImportError as e:
            result.stages_skipped["mmgbsa"] = f"MM/GBSA module unavailable: {e}"
        else:





            score.score_screen(candidates, weights=scoring_weights, use_gbsa=False,
                               dock_attempted=dock_dispatched)
            posed = sorted([c for c in candidates if c.pose and c.pose.sdf_path],
                           key=lambda c: (-c.composite, c.pose.affinity))[:gbsa_top_n]
            try:
                poses = [(c.ligand.ligand_id, c.pose.pose_index, Path(c.pose.sdf_path))
                         for c in posed]
                # Docking kept the structural cofactors; `Receptor.path` has
                # none of them and `mmgbsa.system.prepare_receptor` removes
                # heterogens regardless. Ask now, once, before the hours: a
                # rescore against a pocket missing a wall the poses were docked
                # against is a different quantity, not a better estimate of the
                # same one.
                # Name the platform that will actually be used. `--platform
                # CUDA` was accepted, registered, and then failed on every
                # pose; the report said "20 of 20 rescored poses returned no
                # MM/GBSA energy" and never mentioned the platform.
                if platform:
                    from .mmgbsa.engine import make_platform, platform_error
                    why = platform_error(platform)
                    if why:
                        actual, _ = make_platform(platform)
                        result.stages_run.append(
                            f"MM/GBSA platform: {platform} is registered but "
                            f"will not run here ({why}); scored on "
                            f"{actual.getName()} instead. A historical report "
                            "described a CUDA/OpenCL difference of up to 13.3 "
                            "kcal/mol; that comparison has not been independently "
                            "reproduced. The scorer's resolution does not validate "
                            "platform equivalence or repeat stability. This table is not comparable "
                            f"with one scored on {platform}.")
                        log.warning("MM/GBSA platform %s unusable (%s); using %s",
                                    platform, why, actual.getName())

                from .mmgbsa import cofactor_guard
                refusals = cofactor_guard.check_poses(
                    Path(result.receptor.path), [p for _, _, p in poses])
                if refusals:
                    reason = " ".join(refusals)
                    result.stages_skipped["mmgbsa"] = reason
                    result.gates_unverified["mmgbsa_receptor"] = reason
                    log.warning("MM/GBSA refused: %s", reason)
                    raise _MMGBSARefused(reason)
                # `platform` reaches `make_platform` through mmgbsa_batch's
                # **kwargs. It was never passed, so the CLI had no way to ask
                # for a GPU however capable the host was: every run landed on
                # CPU by default. None keeps that default.
                gdf = gb.mmgbsa_batch(Path(result.receptor.path), poses,
                                      workdir / "mmgbsa", platform=platform,
                                      cpu_profile=gbsa_cpu_profile,
                                      threads=gbsa_threads,
                                      charge_method=gbsa_charge_method,
                                      ligand_forcefield=gbsa_ligand_forcefield)
                _attach_gbsa_records(result, posed, gdf)
            except _MMGBSARefused:
                pass                 # already recorded, with its own reason
            except Exception as e:
                result.stages_skipped["mmgbsa"] = f"{type(e).__name__}: {e}"
                log.warning("MM/GBSA skipped: %s", e)
    elif gbsa:
        result.stages_skipped["mmgbsa"] = "no receptor; docking did not run"

    # ---- cosmetic eligibility --------------------------------------------
    try:
        from .eligibility import CosmeticEligibility
        from .regulatory import Regulations
        regs = Regulations()
        elig = CosmeticEligibility(regulations=regs)
        n_cas = 0
        for c in candidates:



















            cas = _ligand_cas(c.ligand)
            n_cas += 1 if cas else 0
            v = elig.assess(c.ligand.smiles, name=c.ligand.name,
                            admet_profile=c.admet, cas=cas)
            # Whether *this* compound's registry numbers were actually looked
            # up. `sourcing.annotate` writes it; nothing read it, so two
            # candidates differing only in that fact produced byte-identical
            # verdicts -- measured: both `eligible / basis=none /
            # definitive=False / usable=True`, same caveats. One of them had
            # PubChem answer and carry no listed CAS; the other had the lookup
            # fail. The module's standing caveat says a key that was not used
            # makes "not on the list" mean "never compared", and this is the
            # per-candidate half of that sentence: the screen knows which
            # compounds it applies to and was throwing the answer away.
            #
            # Absent key means the sourcing stage never ran at all, which is a
            # property of the screen rather than of the compound and is already
            # recorded once in `stages_skipped["eligibility_cas_lookup"]`.
            if cas and _cas_is_only_a_claim(c.ligand):
                c.ligand.props["cas_verification_status"] = "supplied_unverified"
            elif cas:
                c.ligand.props["cas_verification_status"] = "resolved"
            checked = c.ligand.props.get("cas_checked")
            c.eligibility = {"verdict": v.verdict.value, "basis": v.basis.value,
                             "definitive": v.definitive, "reason": v.reason,
                             "limit": v.limit_text, "caveats": v.caveats,
                             "cas_checked": checked, "cas_used": list(cas),
                             "cas_verification": c.ligand.props.get(
                                 "cas_verification_status")}
            if c.ligand.props.get("cas_verification_status") == "supplied_unverified":
                c.eligibility["definitive"] = False
                c.eligibility["caveats"] = list(c.eligibility["caveats"]) + [
                    "the only registry number for this compound came from the "
                    "input file and has not been checked against its structure, "
                    "so this verdict is what that number says, not what this "
                    "compound is"]
            if checked is False:
                c.eligibility["definitive"] = False
                c.eligibility["caveats"] = list(v.caveats) + [
                    "이 화합물의 CAS 조회가 실패했습니다(PubChem 무응답). 별표 "
                    "1·2에서 실제 대조 범위가 가장 넓은 키가 CAS이므로, 위 "
                    "'목록에 없음'은 이 화합물에 대해 확인된 것이 아니라 "
                    "조회되지 않은 것입니다. 재실행 후 판정하십시오."]
            if not v.usable:
                c.ligand.flags.append(f"cosmetic {v.verdict.value}: {v.reason[:100]}")
        # Counted once for the screen. Holding each unchecked compound off the buy
        # list keeps that list safe, which was the P1; it does not stop
        # `trustworthy` reading true on a screen whose prohibition check was
        # partial. `trustworthy` is defined as "no unverified gates", so the fact
        # has to arrive as a gate rather than as a per-row field.
        n_cas_failed = sum(1 for c in candidates
                           if (c.eligibility or {}).get("cas_checked") is False)
        if n_cas_failed:
            result.gates_unverified["eligibility_cas"] = (
                f"{n_cas_failed} of {len(candidates)} candidates had their CAS "
                "lookup fail, so 별표 1·2 was never consulted on the key it has "
                f"coverage on (CAS for 632 of 656 rows against 12 structures). "
                "Those compounds are held off the orderable list rather than "
                "cleared, and the screen's regulatory check is partial until a "
                "re-run resolves them.")
        if candidates and not n_cas:
            # Declarative, not conditional. `LigandRecord` has no CAS field
            # and no fetch stage writes `props["cas"]`, so on a screen built
            # from ChEMBL this route does not run -- and a caveat phrased as
            # "한글 고시명 또는 CAS로 조회하지 않았다면" invites the reader to
            # check something that never happens. What did not happen is
            # recorded here, once, with the count.
            result.stages_skipped["eligibility_cas_lookup"] = (
                f"0 of {len(candidates)} candidates carried a CAS number, so "
                "the CAS index of 별표 1·2 -- the key that list has real "
                "coverage on -- was never consulted in this screen. The "
                "sourcing stage is what supplies registry numbers (PubChem "
                "synonyms, check-digit verified), so this is what a screen run "
                "with check_sourcing=False looks like, or one where PubChem "
                "answered for nothing. Every 'not on the list' verdict below "
                "rests on structure and name matching alone.")
        n_bad = sum(1 for c in candidates
                    if (c.eligibility or {}).get("verdict") in ("prohibited", "drug_only"))
        # A gate fails in two ways and only one of them raises. The constructor
        # is happy to carry on with no list at all -- a path that does not
        # exist, a YAML that will not parse, a file that parses to nothing --
        # and every compound then comes back "not on the list", which is
        # indistinguishable from a list that cleared it. `list_usable` is the
        # only thing that tells those apart, and nothing used to ask.
        if not elig.list_usable:
            why = elig.list_caveat() or f"list status {elig.list_status.value}"
            result.stages_skipped["eligibility_list_unavailable"] = why
            result.gates_unverified["eligibility"] = (
                f"the prohibited-ingredient list was never loaded "
                f"({elig.list_status.value}: {elig.list_detail}); no candidate in "
                "this screen has been checked against it, so every verdict below "
                "rests on the element screen and structural alerts alone. " + why)
            result.stages_run.append(
                f"cosmetic eligibility (NOT VERIFIED: list {elig.list_status.value}; "
                f"{n_bad} barred by element/drug evidence alone)")
            log.warning("prohibited-ingredient list unusable (%s): %s",
                        elig.list_status.value, elig.list_detail)
        else:
            result.stages_run.append(f"cosmetic eligibility ({n_bad} barred)")
        if not elig.mfds_coverage.get("present"):
            result.gates_unverified["eligibility_mfds"] = (
                "Official MFDS API restriction data were unavailable; the curated "
                "statutory list and structural assessment were retained, but official-only "
                "restriction coverage is unverified.")
    except Exception as e:
        result.stages_skipped["eligibility"] = f"{type(e).__name__}: {e}"
        # Recording this only as a skipped stage let "could not assess" be read
        # as "nothing bars it": trustworthy stayed true and the buy list still
        # shipped. A gate that did not run is unverified, and the report has to
        # lead with that.
        # `gates_unverified` is keyed by stage, not a list: the `.append` this
        # replaced raised AttributeError inside the handler and took the whole
        # screen down with it, and the test guarding the line only grepped the
        # source for the word.
        result.gates_unverified["eligibility"] = (
            f"cosmetic eligibility never ran ({type(e).__name__}); no candidate "
            "in this screen has been checked against the prohibited-ingredient "
            "list or its notified status")
        log.warning("eligibility skipped: %s", e)

    # ---- score -----------------------------------------------------------
    result.candidates = candidates

    result.breakdowns = score.score_screen(
        candidates, weights=scoring_weights, dock_attempted=dock_dispatched)
    result.control_check = score.control_check(candidates)
    # After `result.breakdowns`, and it has to be: `_note_composite_scale`
    # reads them and returns early when the list is empty. It used to be called
    # 138 lines earlier, before scoring had run, so the gate it raises never
    # fired in a real screen -- the whole of commit 90d596b was dead on the
    # only path that matters. The unit test passed because it called the helper
    # directly with a populated result, and the call-site test asserted through
    # the AST that the call *exists*. Existence is not ordering.
    _note_composite_scale(result)
    result.stages_run.append("scoring")
    # Filled in before the result leaves, so nothing has to call `orderable()`
    # first to make the field true. It was a side effect of that call alone, so
    # `report.to_json` on a fresh result serialised `"held_for_regulatory": []`
    # -- "nothing was held" -- for a screen that had held a prohibited
    # ingredient back from the buy list.
    result.held_for_regulatory = result.regulatory_holds()
    return result
