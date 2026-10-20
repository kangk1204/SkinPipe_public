"""Turn a ScreenResult into something a formulator can read and act on."""
from __future__ import annotations

import collections
import csv
from dataclasses import asdict
import json
import os
import re
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING, Optional

from .types import SEVERITY_ORDER, concerns_of, hazard_severity_of
from .candidate_use import classify, compact_goal_support, FUNCTIONAL_PRIORITY, EXPLORATORY, COUNTER_CONTROL
from .fetch.provenance import CONTEXT_FIELDS
from .score import colour_of, _electrophilic_alerts, skin_filter_status
from .measurements import activity_status
from .validation_contract import candidate_contract, aggregate

if TYPE_CHECKING:
    from .pipeline import ScreenResult

_COLUMNS = ["rank", "name", "ligand_id", "source", "composite", "pchembl",
            "activity_status", "observation_status_counts",
            "skin_assessed", "sensitisation_status", "stability_status", "n_pains", "score_terms", "weights_used", "weighted_terms",
            "weight_denominator", "evidence_factor", "scoring_policy", "composite_rounding_policy",
            "composite_without_evidence_factor", "evidence_factor_contribution", "evidence_diagnostic_policy", "structure_tier_policy",
            "structure_missing_policy", "docking_tie_tolerance_kcal", "rank_within_group", "rank_group_size", "rank_interpretation",
            "measurements", "measurement_groups", "measurement_policy",
            "provenance_policy_eligible", "provenance_gate_deferred",
            "primary_measurement_group", "primary_measurement_status", "primary_measurement_limitation", "n_exact_measurements", "n_potency_quality_excluded", "exact_potency_scope", "potency_aggregation", "pchembl_range", "pchembl_range_scope", "pchembl_range_all_contexts",
            "curated_claims", "curated_claim_policy", "curated_claim_status",
            "natural_product", "natural_product_sources", "original_structures",
            "dock_affinity", "docking_execution", "mmgbsa_dG", "mmgbsa_result", "structure_evidence", "mw", "logp", "log_kp",
            "permeation", "topical_score", "sensitiser_alerts", "sensitiser_details",
            "stability_alerts", "colour", "rank_group", "evidence",
            # What the provenance call rests on, and whether a locally rescued
            # structure came from a verified index. Both were computed and
            # then dropped before the row was written.
            "provenance", "provenance_basis",
            "structure_index_verified", "structure_index_status",
            "annotation_gap",







            "potency_disputed",






            "direction_verdict", "direction_required", "direction_effect",
            "direction_source", "direction_explanation", "direction_evidence",
            "candidate_use", "candidate_use_reason", "direction_supported",
            "goal_support", "experimental_validation", "validation_status", "dock_pose_index",
            "purchasable", "sa_score", "vendor_id", "n_vendors", "vendors",
            "in_patents", "sourcing_status",







            "ip_note", "dossier_note",







            "cas", "cas_all", "cas_source", "cas_supplied",
            "cas_verification",
            "admet_worst", "admet_hazard", "admet_verdict",
            "cosmetic_verdict", "cosmetic_basis", "cosmetic_limit",
            "cosmetic_caveats",
            "caveats", "smiles", *CONTEXT_FIELDS]

#: How the Markdown names each `dropped[*]["stage"]`. A drop announced under
#: the wrong gate sends the reader to fix the wrong thing: a VDR agonist
#: removed for pushing the node backwards was reported as an ADMET failure, and
#: an ADMET profile is not what is wrong with it.
_STAGE_LABELS = {
    "direction": "the direction gate",
    "admet": "the ADMET gate",
    "eligibility": "the cosmetic-eligibility gate",
    "sourcing": "the sourcing gate",
    "topical": "the skin-permeation filter",
}

#: What the `direction` column means. It is printed above the first table that
#: carries the column -- the buy list, when there is one -- because a reader who
#: meets `unknown` without it reads it as a pass, and it is not one.
_DIRECTION_LEGEND = (
    "> `direction` is this compound's own mechanism verdict against "
    "the node: `compatible` pushes it the way the marker needs, "
    "`opposing` pushes it back, `unknown` means nothing is recorded "
    "-- which is the common case and is not a pass -- and "
    "`no_reference` means the node's required direction is unset.")

#: Per-compound drop rows printed in the Markdown before the list is summarised
#: instead. `result.json` carries every one of them either way.
_MAX_DROPPED_ROWS = 25

#: Reason groups printed in the roll-up before the rest is summarised instead.
_MAX_REASON_GROUPS = 8

#: Numbers inside a free-text drop reason. Two drops whose reasons differ only
#: here are the same finding about two different compounds.
_NUMBER = re.compile(r"-?\d+(?:\.\d+)?")

#: A numeric comparison inside a drop reason: value, operator, threshold.
_COMPARISON = re.compile(
    r"(?<![\w.])(-?\d+(?:\.\d+)?)\s*(<=|>=|<|>)\s*(-?\d+(?:\.\d+)?)")

#: The permeation class as the skin filter writes it into its drop reason.
_PERMEATION = re.compile(r"permeation ([A-Za-z_]+)")


def _stage_label(stage: Optional[str]) -> str:
    return _STAGE_LABELS.get(stage or "", f"the {stage} stage" if stage
                             else "an unrecorded stage")


def _cell(value) -> str:
    """One Markdown table cell. A drop reason is free text from a mechanism
    row, and an unescaped pipe in it silently shifts every column after it."""
    return str(value if value is not None else "").replace(
        "|", r"\|").replace("\n", " ").strip()


def _cell_or_dash(value) -> str:
    'One Markdown cell, with an em dash where there is no value.'
    return _cell(value) or "—"


def _ip_cell(props) -> str:
    """The IP cell, three-valued, because two of the three are not the same fact.

    `in_patents` is True (SureChEMBL has a record), False (the lookup ran and
    found none), or absent -- no lookup ever ran. A generated scaffold hop is
    the third case: `hop_rank.to_ligand_records` builds its `props` from mw,
    logp, tpsa and hbd, and nothing puts a patent key there. Rendering that as
    `clear` is exactly the confusion `_cell_or_dash` above exists to prevent,
    on the one column where being wrong costs a freedom-to-operate opinion --
    and `hopping` itself warns that a novel Murcko scaffold can still sit
    inside someone else's composition-of-matter claim.
    """
    flag = props.get("in_patents")
    if flag is None:
        return "not checked"
    return "patented" if flag else "none found"


def _ip_unchecked(candidates) -> bool:
    """Does any row carry no patent lookup at all?"""
    return any(c.ligand.props.get("in_patents") is None for c in candidates)


def _held_section(result: "ScreenResult") -> list[str]:
    'The compounds that ranked but did not reach the buy list.'
    held = _regulatory_holds(result)
    if not held:
        return []
    by_stage = collections.Counter((h.get("stage") or "") for h in held)
    per_stage = ", ".join(f"{n} by {_stage_label(s)}"
                          for s, n in by_stage.most_common())
    noun = "compound" if len(held) == 1 else "compounds"
    L = ["## Held off the buy list", "",
         f"{len(held)} ranked {noun} cleared every gate and still did not reach "
         f"`orderable.csv`: {per_stage}. The `buy` column in the tables above is "
         "vendor availability alone and does not carry these holds.", "",
         "| compound | stage | verdict | reason |", "|---|---|---|---|"]
    for h in held[:_MAX_DROPPED_ROWS]:
        L.append(f"| {_cell(h.get('name') or h.get('ligand_id', '?'))} "
                 f"| {_cell(_stage_label(h.get('stage')))} "
                 f"| {_cell(h.get('verdict'))} "
                 f"| {_cell(_readable_reason(h.get('reason')) or 'unrecorded')} |")
    if len(held) > _MAX_DROPPED_ROWS:
        L.append(f"| *+{len(held) - _MAX_DROPPED_ROWS} more* | | | "
                 "*all of them are in `result.json`* |")
    return L + [""]


def _readable_reason(reason) -> str:
    'A drop reason with its rounding artefacts named instead of printed.'
    def repair(m: "re.Match[str]") -> str:
        lhs, op, rhs = m.group(1), m.group(2), m.group(3)
        # `6.0 <= 6.0` is true as written and needs nothing. Only a *strict*
        # comparison between two numbers that print the same is a contradiction.
        if op not in ("<", ">") or float(lhs) != float(rhs):
            return m.group(0)
        return f"{'below' if op == '<' else 'above'} {rhs}, rounds to {lhs}"

    return _COMPARISON.sub(repair, str(reason or ""))


def _permeation_class(d: dict) -> str:
    """The permeation class behind a drop record, however it is carried."""
    v = d.get("permeation") or d.get("permeation_class")
    if v:
        return str(v)
    m = _PERMEATION.search(str(d.get("reason") or ""))
    return m.group(1) if m else ""


def _reason_group_key(d: dict) -> tuple[str, str, str, str]:
    'What makes two drops the same finding.'
    return (str(d.get("stage") or ""), str(d.get("severity") or ""),
            _permeation_class(d), _NUMBER.sub("#", str(d.get("reason") or "")))


def _group_label(reasons: list[str]) -> str:
    """One reason line for a group, with the numbers that vary elided.

    The numbers the group agrees on are the ones worth printing -- the floor
    that was applied is the same for everybody -- and the ones that differ are
    per-compound and belong in the table above, not in a roll-up speaking for
    fifty-two compounds at once.
    """
    first = str(reasons[0] or "")
    found = [_NUMBER.findall(str(r or "")) for r in reasons]
    if not found or not found[0]:
        return first
    keep = [col[0] if len(set(col)) == 1 else "\u2026" for col in zip(*found)]
    parts = _NUMBER.split(first)
    return "".join(lit + num for lit, num in zip(parts, keep)) + parts[-1]


def _marker_notes(result: "ScreenResult") -> list[str]:
    "The registry's note about the target, as against a warning about this run."
    t = getattr(result, "triage", None)
    for obj, attr in ((result, "marker_notes"), (t, "marker_notes"),
                      (t, "notes"), (getattr(result, "spec", None), "notes")):
        v = getattr(obj, attr, None) if obj is not None else None
        if isinstance(v, (list, tuple)):
            return [str(x) for x in v]
    marker = str(getattr(result, "marker", "") or "")
    if not marker:
        return []
    try:
        from .registry import Registry
        return [str(n) for n in Registry().spec(marker, resolve=False).notes]
    except Exception:          # no registry, unknown marker, unreadable YAML
        return []


def _split_warnings(result: "ScreenResult") -> tuple[list[str], list[str]]:
    """`(marker notes, real holds)` out of the one list triage hands over."""
    warnings = list(getattr(getattr(result, "triage", None), "warnings", None) or [])
    if not warnings:           # nothing to classify; do not go to the registry
        return [], []
    notes = _marker_notes(result)
    return ([w for w in warnings if w in notes],
            [w for w in warnings if w not in notes])


def _tints(c) -> bool:
    """True when the structural colour model flags possible tinting."""
    v = str(getattr(getattr(c, "skin", None), "colour", "") or "").strip().lower()
    return bool(v) and not v.startswith(("low", "unknown", "none"))


def _hazard_domain(c) -> str:
    """The domain driving this compound's hazard severity, if one does."""
    h = hazard_severity_of(getattr(c, "admet", None))
    return next((cc.domain for cc in concerns_of(getattr(c, "admet", None))
                 if cc.severity == h and cc.domain != "systemic"), "")


def _hazard_cell(c) -> str:
    'The `ADMET hazard` cell: the severity, and what is driving it.'
    h = hazard_severity_of(getattr(c, "admet", None))
    if h == "unknown":
        return "not assessed"
    dom = _hazard_domain(c)
    return f"{h}: {dom}" if dom else h


def _admet_summary(ranked: list) -> list[str]:
    'What the hazard and colour columns say about the list that survived.'
    if not ranked:
        return []
    sev = collections.Counter(hazard_severity_of(getattr(c, "admet", None))
                              for c in ranked)
    ladder = list(SEVERITY_ORDER) + ["unknown"]
    ordered = [lv for lv in ladder if sev.get(lv)] + \
              [lv for lv in sorted(sev) if lv not in ladder]
    noun = "candidate" if len(ranked) == 1 else "candidates"
    out = [f"The `ADMET hazard` column over {len(ranked)} ranked {noun} reads "
           + ", ".join(f"{sev[lv]} {lv}" for lv in ordered) +
           ". A gate set above a severity does not clear what it leaves "
           "behind, so every count here is a property of a compound that is "
           "still on the list."]
    serious = [c for c in ranked
               if hazard_severity_of(getattr(c, "admet", None)) in ("critical", "high")]
    if serious:
        dom = collections.Counter(_hazard_domain(c) or "unattributed" for c in serious)
        out.append(f"The {len(serious)} at critical or high are driven by "
                   + ", ".join(f"{d} ({n})" for d, n in dom.most_common()) + ".")
    n_tint = sum(1 for c in ranked if _tints(c))
    if n_tint:
        out.append(f"{n_tint} of {len(ranked)} carry a predicted colour risk that may "
                   f"tint a white base; verify the intended formulation.")
    return [" ".join(out), ""]


def _cosmetic_cell(eligibility: Optional[dict]) -> str:
    'The regulatory cell of the buy list.'
    e = eligibility or {}
    verdict = str(e.get("verdict") or "")
    basis = str(e.get("basis") or "")
    if verdict in ("eligible", "unknown") and basis in ("", "none", "not_assessed"):
        return "not_listed"
    return verdict


def _decision_json(d) -> Optional[dict]:
    """`decision` as a JSON object, or None when the screen produced none.

    The whole Discovery route section -- the track, why it was picked, the IP
    outlook, the steps, and the warnings that qualify them -- had no
    machine-readable counterpart at all: `to_json` never serialised it.
    """
    if d is None:
        return None
    track = getattr(d, "track", None)
    return {"track": getattr(track, "value", track),
            "node": getattr(d, "node", None),
            "rationale": getattr(d, "rationale", ""),
            "ip_outlook": getattr(d, "ip_outlook", ""),
            "next_steps": list(getattr(d, "next_steps", None) or []),
            "warnings": list(getattr(d, "warnings", None) or []),
            "n_ligands": getattr(d, "n_ligands", None),
            "n_structures": getattr(d, "n_structures", None),
            "patent_fraction": getattr(d, "patent_fraction", None),
            "purchasable_fraction": getattr(d, "purchasable_fraction", None)}


def _regulatory_holds(result: "ScreenResult") -> list[dict]:
    """The candidates held out of the buy list -- computed, not remembered.

    `held_for_regulatory` used to be filled in as a side effect of
    `ScreenResult.orderable()`, so the *same* result serialised
    `"held_for_regulatory": []` when `to_json` ran first and the real list when
    the Markdown or the CSV had run before it. `result.json` is the
    machine-readable artefact: an automation that writes only the JSON read
    "no candidate is in regulatory limbo", which is the one claim the field
    exists to refute. Ask the result to compute the list whenever it can, and
    read the stored field only for a result that cannot.
    """
    compute = getattr(result, "regulatory_holds", None)
    if callable(compute):
        return list(compute())
    return list(getattr(result, "held_for_regulatory", None) or [])


def _write_atomically(path: Path, text: str) -> Path:
    """Write via a temp file in the same directory, then rename.

    A screen that dies partway through emitting a report leaves a file that is
    syntactically fine and silently short: a truncated CSV is a shorter hit
    list, not an error, and `result.json` half-written is not valid JSON at all.
    The rename is atomic on POSIX, so a reader sees the old file or the new one.
    """
    path = Path(path)
    from .publication import ensure_mutable_workdir
    ensure_mutable_workdir(path.parent)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=str(path.parent), prefix=f".{path.name}.", suffix=".part")
    try:
        with os.fdopen(fd, "w", newline="", encoding="utf-8") as fh:
            fh.write(text)
        os.replace(tmp, path)
    except BaseException:
        Path(tmp).unlink(missing_ok=True)
        raise
    return path


def _row(c, bd=None) -> dict:
    lig, skin = c.ligand, c.skin
    use = classify(c)
    # `None` means the direction stage never ran for this screen, which is not
    # the same as running and finding nothing recorded ("unknown"). An empty
    # cell would read as the second, so the first is spelled out.
    d = c.direction if isinstance(getattr(c, "direction", None), dict) else None
    direction_verdict = d.get("verdict", "") if d else "not run"
    row = {
        "rank": c.rank, "name": lig.name or "", "ligand_id": lig.ligand_id,
        "source": lig.source, "composite": c.composite if bd is not None or c.rank else None,
        "pchembl": lig.activity.get("pchembl"),
        "activity_status": activity_status(lig.activity),
        "observation_status_counts": lig.activity.get("observation_status_counts"),
        "skin_assessed": all(skin_filter_status(c, domain) == "EVALUATED"
                             for domain in ("sensitisation", "stability")),
        "sensitisation_status": skin_filter_status(c, "sensitisation"),
        "stability_status": skin_filter_status(c, "stability"),
        "n_pains": (skin.n_pains if skin_filter_status(c, "sensitisation") == "EVALUATED"
                    else None),
        "score_terms": {k: getattr(bd, k) for k in
                        ("potency", "structure", "topical", "safety", "stability")} if bd else {},
        "weights_used": getattr(bd, "weights_used", {}) if bd else {},
        "weighted_terms": getattr(bd, "weighted_terms", []) if bd else [],
        "weight_denominator": getattr(bd, "weight_denominator", None) if bd else None,
        "evidence_factor": getattr(bd, "evidence_factor", None) if bd else None,
        "scoring_policy": getattr(bd, "scoring_policy", None) if bd else None,
        "composite_rounding_policy": getattr(bd, "composite_rounding_policy", None) if bd else None,
        **{key: getattr(bd, key, None) if bd else None for key in (
            "composite_without_evidence_factor", "evidence_factor_contribution", "evidence_diagnostic_policy", "structure_tier_policy",
            "structure_missing_policy", "docking_tie_tolerance_kcal", "rank_within_group", "rank_group_size", "rank_interpretation")},
        **{k: lig.activity.get(k) for k in ("measurements", "measurement_groups",
            "provenance_policy_eligible", "provenance_gate_deferred",
            "measurement_policy", "primary_measurement_group", "primary_measurement_status", "primary_measurement_limitation", "n_exact_measurements", "n_potency_quality_excluded", "exact_potency_scope", "potency_aggregation",
            "pchembl_range", "pchembl_range_scope", "pchembl_range_all_contexts", "curated_claims", "curated_claim_policy", "curated_claim_status")},
        "natural_product": lig.natural_product,
        "natural_product_sources": lig.natural_product_sources,
        "original_structures": lig.original_structures,
        "potency_disputed": lig.activity.get("potency_disputed", ""),
        "dock_affinity": c.pose.affinity if c.pose else None,
        "dock_pose_index": getattr(c.pose, "pose_index", None),
        "docking_execution": getattr(c, "docking_execution", None),
        "mmgbsa_dG": c.gbsa.dG_bind if (c.gbsa and c.gbsa.converged) else None,


        "mmgbsa_result": asdict(c.gbsa) if c.gbsa is not None else None,
        # MM/GBSA runs on a subset, so one ranked column mixes rescored and
        # docking-only candidates. The ordering is sound; the mixture was
        # invisible, which let a reader take the whole column as one scale.
        "structure_evidence": getattr(bd, "structure_evidence", "") if bd else "",
        # Which block of the ranking the row sits in. Ranks are comparable
        # inside a block and not across one.
        "rank_group": getattr(bd, "rank_group", "") if bd else "",
        "mw": lig.props.get("mw"), "logp": lig.props.get("logp"),
        "log_kp": skin.log_kp if skin else None,
        "permeation": skin.permeation_class if skin else "",
        "topical_score": skin.topical_score if skin else None,
        "sensitiser_alerts": "|".join(skin.sensitiser_alerts) if skin else "",
        "sensitiser_details": getattr(skin, "sensitiser_details", []) if skin else [],
        "stability_alerts": "|".join(skin.stability_alerts) if skin else "",
        # Not an alert and never counted as one, which is exactly why it used
        # to reach the prose and stop there while the stability column read 10.
        "colour": colour_of(c),
        "evidence": c.evidence.value,


        "provenance": lig.activity.get("provenance", ""),
        "provenance_basis": lig.activity.get("provenance_basis", ""),
        **{k: lig.activity.get(k, "") for k in CONTEXT_FIELDS},





        "structure_index_verified": lig.activity.get("structure_index_verified", ""),
        "structure_index_status": lig.activity.get("structure_index_status", ""),
        "annotation_gap": lig.activity.get("annotation_gap", ""),
        "direction_verdict": direction_verdict,
        "direction_required": (d or {}).get("required") or "",
        "direction_effect": (d or {}).get("effect") or "",
        "direction_source": (d or {}).get("source") or "",
        "direction_explanation": (d or {}).get("explanation") or "",
        'direction_evidence': lig.activity.get('direction_evidence') or [],
        'candidate_use': use['role'], 'candidate_use_reason': use['reason'],
        'direction_supported': use['direction_supported'],
        'goal_support': compact_goal_support(use['goal_support']),
        'experimental_validation': use['experimental_validation'],
        "purchasable": lig.props.get("purchasable", ""),
        # How hard it would be to make, for the rows no vendor stocks. Blank
        # for a purchasable compound (nobody needs a route) and for an
        # unparseable SMILES (which must not read as an easy molecule).
        "sa_score": lig.props.get("sa_score", "") if lig.props.get("sa_score") is not None else "",
        "vendor_id": lig.props.get("vendor_id", "") or "",
        "in_patents": lig.props.get("in_patents", ""),
        "sourcing_status": lig.props.get("sourcing_status", ""),
        "ip_note": lig.props.get("ip_note", "") or "",
        "dossier_note": lig.props.get("dossier_note", "") or "",
        "cas": lig.props.get("cas", "") or "",
        "cas_all": "; ".join(lig.props.get("cas_all") or []),
        "cas_source": lig.props.get("cas_source", "") or "",
        # What the input file claimed, kept in its own column so it is never
        # read as a resolved identifier.
        "cas_supplied": lig.props.get("cas_supplied", "") or "",
        "cas_verification": lig.props.get("cas_verification_status", "") or "",
        # worst_severity includes systemic exposure, which the hazard gate
        # deliberately ignores; printing only the former makes the gate look
        # broken (niacinamide reads "high" and is kept). Both are shown.
        "admet_worst": (c.admet or {}).get("worst_severity", ""),
        "admet_hazard": hazard_severity_of(c.admet),
        "admet_verdict": (c.admet or {}).get("overall", ""),
        "n_vendors": lig.props.get("n_vendors", ""),
        "vendors": "; ".join(lig.props.get("vendors", [])[:3])
                   if isinstance(lig.props.get("vendors"), list) else "",
        "cosmetic_verdict": (c.eligibility or {}).get("verdict", ""),
        "cosmetic_basis": (c.eligibility or {}).get("basis", ""),
        "cosmetic_limit": (c.eligibility or {}).get("limit", "") or "",





        "cosmetic_caveats": " | ".join((c.eligibility or {}).get("caveats") or []),
        # Computed for every candidate since the first version of score.py and
        # never emitted anywhere. These are the sentences that qualify the
        # number in the `composite` column.
        "caveats": " | ".join(getattr(bd, "caveats", []) or []),
        "smiles": lig.smiles,
    }
    row["validation_status"] = candidate_contract(
        dict(row, direction_verdict=use["direction_verdict"]),
        reviewed_direction=use["direction_supported"])
    return row


def _csv_text(rows, bd: dict) -> str:
    import io
    buf = io.StringIO()
    buf.write("\ufeff")
    w = csv.DictWriter(buf, fieldnames=_COLUMNS)
    w.writeheader()
    for c in rows:
        row = _row(c, bd.get(c.ligand.ligand_id))
        w.writerow({k: _csv_cell(v) for k, v in row.items()})
    return buf.getvalue()


def _csv_cell(value):
    """Keep numerical cells numeric and escape textual formula prefixes."""
    if isinstance(value, (dict, list)):
        value = json.dumps(value, sort_keys=True, ensure_ascii=False, allow_nan=False)
    if isinstance(value, str) and value.startswith(("=", "+", "-", "@", "\t", "\r", "\n")):
        value = "'" + value
    return value


def to_csv(result: "ScreenResult", path: Path, top_n: Optional[int] = None) -> Path:
    bd = {b.ligand_id: b for b in result.breakdowns}
    rows = result.top(top_n) if top_n else sorted(
        [c for c in result.candidates if c.rank], key=lambda c: c.rank)
    return _write_atomically(Path(path), _csv_text(rows, bd))


def discovery_to_csv(result: "ScreenResult", path: Path) -> Path:
    """Complete pre-gate population, sorted by the versioned exploratory score."""
    import io
    from .discovery import rows
    records = rows(result, _row)
    buf = io.StringIO()
    buf.write("\ufeff")
    keys = list(dict.fromkeys(k for record in records for k in record))
    writer = csv.DictWriter(buf, fieldnames=keys or ["ligand_id"])
    writer.writeheader()
    for record in records:
        writer.writerow({k: _csv_cell(v) for k, v in record.items()})
    return _write_atomically(path, buf.getvalue())


def candidate_exports(result: "ScreenResult", directory: Path,
                      summary_top: Optional[int] = 300) -> dict:
    """The default CSV is complete; ranked and summary lists have explicit names."""
    if summary_top is not None and (isinstance(summary_top, bool)
            or not isinstance(summary_top, int) or summary_top < 0):
        raise ValueError("summary_top must be a non-negative integer or None")
    directory = Path(directory)
    complete = discovery_to_csv(result, directory / "candidates.csv")
    _write_atomically(directory / "all_candidates.csv", complete.read_bytes().decode("utf-8"))
    to_csv(result, directory / "ranked_candidates.csv")
    to_csv(result, directory / "candidate_summary.csv", top_n=summary_top)
    contract = {"contract": "candidate-exports-v2", "csv_cell_policy": "utf8-bom-formula-escaped-v1",
            "status": "EXPORTED", "complete": "candidates.csv",
            "complete_alias": "all_candidates.csv", "ranked": "ranked_candidates.csv",
            "summary": "candidate_summary.csv", "summary_top": summary_top,
            "export_limit": summary_top or None,
            "n_acquired": len(result.all_candidates or result.candidates),
            "n_ranked": sum(bool(c.rank) for c in result.candidates)}
    contract["n_summary"] = min(contract["n_ranked"], summary_top) if summary_top else contract["n_ranked"]
    result.export_contract = contract
    return contract


def orderable_to_csv(result: "ScreenResult", path: Path,
                     max_hazard: str = "moderate") -> tuple[Path, int]:
    """Available functional-testing priorities; cosmetic permission is unestablished."""
    rows = result.orderable(max_hazard=max_hazard)
    bd = {b.ligand_id: b for b in result.breakdowns}
    return _write_atomically(Path(path), _csv_text(rows, bd)), len(rows)


def research_lists_to_csv(result: "ScreenResult", directory: Path) -> dict[str, int]:
    """Export separate purpose lists and research availability without changing rank."""
    bd = {b.ligand_id: b for b in result.breakdowns}
    counts = {}
    for purpose in (FUNCTIONAL_PRIORITY, EXPLORATORY, COUNTER_CONTROL):
        rows = result.research_candidates(purpose)
        _write_atomically(Path(directory) / (purpose + '.csv'), _csv_text(rows, bd))
        counts[purpose] = len(rows)
    _write_atomically(Path(directory) / 'research_available.csv',
                     _csv_text(result.research_available(), bd))
    controls = [r for r in result.dropped if r.get('stage') == 'direction'
                and r.get('candidate_use') == COUNTER_CONTROL]
    _write_atomically(Path(directory) / 'excluded_counter_controls.json',
                     json.dumps(controls, ensure_ascii=False, indent=2))
    return counts


def _porcelain_z_paths(raw: str) -> list[str]:
    """Paths from `git status --porcelain -z`, renames included.

    Records are NUL-separated. A normal one is `XY <path>`; a rename is
    `RX <new>` followed by a bare record holding the old path, which has no
    status column and must not be sliced like one.
    """
    out: list[str] = []
    records = [r for r in raw.split("\0") if r]
    expect_old = False
    for rec in records:
        if expect_old:
            out.append(rec)
            expect_old = False
            continue
        if len(rec) > 3 and rec[2] == " ":
            out.append(rec[3:])
            expect_old = rec[0] in "RC" or rec[1] in "RC"
        else:                                   # not a shape we recognise
            out.append(rec)
    return out





RESULT_SCHEMA_VERSION: int = 9


def _json_safe(value):
    """Replace non-finite floats with `null`, everywhere in the payload.

    `score_screen` already refuses a NaN pChEMBL -- potency 0.0, a caveat
    saying why, a finite composite. `_row` then reads the *original* activity
    dict again, so the NaN reaches `json.dumps`, which writes the bare token
    `NaN`. That is not JSON: Python's own strict parser rejects it and so does
    `JSON.parse`, so one bad input makes a browser or any standard consumer
    fail to open the whole file rather than lose one cell.

    `null` is the honest encoding of "this number does not exist", and the
    reason it does not is already in that row's caveats. Applied here, at the
    boundary, rather than by making the export raise: refusing to write the
    file trades a partly-unreadable result for no result.
    """
    import math

    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, dict):
        return {k: _json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(v) for v in value]
    return value


def to_json(result: "ScreenResult", path: Path) -> Path:
    bd = {b.ligand_id: b for b in result.breakdowns}
    # The marker's registry note and "this gate did not run" used to share one
    # array. They are different kinds of claim, and a consumer reading the JSON
    # cannot separate them after the fact.
    notes, holds = _split_warnings(result)
    payload = {
        "schema_version": RESULT_SCHEMA_VERSION,
        "marker": result.marker,
        "triage": {"grade": result.triage.grade.value,
                   "provisional": result.triage.provisional,
                   "node": result.triage.screening_node,
                   "recommendation": result.triage.recommendation,
                   "notes": notes, "warnings": holds},




        "decision": _decision_json(getattr(result, "decision", None)),
        "stages_run": result.stages_run, "stages_skipped": result.stages_skipped,
        "gates_unverified": getattr(result, "gates_unverified", {}),
        "trustworthy": getattr(result, "trustworthy", False),
        "trustworthy_scope": "Legacy processing flag; not candidate scientific validation.",
        "provenance_audit": {k: v for k, v in result.provenance_audit.items() if k != "calls"},
        "control_check": result.control_check,






        "direction": result.direction,
        'target_design': result.target_design,
        'candidate_lists': {purpose: [c.ligand.ligand_id for c in result.research_candidates(purpose)]
                            for purpose in (FUNCTIONAL_PRIORITY, EXPLORATORY, COUNTER_CONTROL)},
        'research_available': [c.ligand.ligand_id for c in result.research_available()],
        "dropped": list(result.dropped),
        "held_for_regulatory": _regulatory_holds(result),
        "filter_profile": result.filter_profile,
        "library_scale": result.library_scale,
        "potency_floor": result.potency_floor,
        "candidate_exports": dict(getattr(result, "export_contract", {}) or {
            "status": "NOT_EXPORTED", "complete": "candidates.csv",
            "ranked": "ranked_candidates.csv", "summary": "candidate_summary.csv"}),
        "structure_validation": dict(getattr(result, "structure_validation", {}) or {}),
        "discovery_candidates": __import__("skinpipe.discovery", fromlist=["rows"]).rows(result, _row),
        "candidates": [_row(c, bd.get(c.ligand.ligand_id)) for c in
                       sorted([x for x in result.candidates if x.rank], key=lambda x: x.rank)],
    }
    payload["validation_status"] = aggregate([
        c["validation_status"] for c in payload["discovery_candidates"]])
    return _write_atomically(
        Path(path), json.dumps(_json_safe(payload), indent=2,
                               ensure_ascii=False, default=str,
                               # The last line of defence, and it can only
                               # fire if `_json_safe` missed a type: better a
                               # loud failure here than a file no parser reads.
                               allow_nan=False))


def to_markdown(result: "ScreenResult", path: Optional[Path] = None, top_n: int = 25) -> str:
    t = result.triage
    unverified = getattr(result, "gates_unverified", {})
    L: list[str] = [f"# {result.marker} — topical small-molecule screen", ""]

    L += ["## Verdict", "",
          f"- **Grade {t.grade.value}**"
          + (" *(provisional: graded on raw counts)*" if t.provisional else ""),
          f"- Screening node: **{t.screening_node or 'none'}**",
          f"- {t.recommendation}", ""]





    marker_notes, read_first = _split_warnings(result)
    scientific = aggregate(_row(c)["validation_status"] for c in
        (getattr(result, "all_candidates", None) or result.candidates))
    L += ["## Scientific validation scope", "",
          "Processing flags, control redocking and exploratory scores do not certify these candidate claims.", "",
          "| Claim | Status |", "|---|---|"]
    from .validation_contract import DEFINITIONS
    L += [f"| {key} | {scientific[key]['status']} |" for key in DEFINITIONS]
    L += ["", "All states refer to their recorded scope; project experimental validation requires its own evidence.", ""]
    if marker_notes:
        L += ["**Why this target**", ""] + [f"- {w}" for w in marker_notes] + [""]
    if unverified:
        read_first.append(
            f"{len(unverified)} gate(s) did not run ({', '.join(unverified)}). "
            f"Their silence is missing information, not a finding -- see "
            f"\u201cIncomplete stages\u201d below.")
    prov_warning = str((result.provenance_audit or {}).get("warning") or "")
    if prov_warning:
        # It qualifies every ranking number in the report and used to appear
        # only halfway down, under a table of assay counts.
        read_first.append(prov_warning)
    if read_first:
        L += ["**Read before using these results**", ""]
        L += [f"- {w}" for w in read_first] + [""]

    d = getattr(result, "decision", None)
    if d is not None:
        L += ["## Discovery route", "",
              f"**{getattr(d.track, 'value', d.track)}** — {d.rationale}", ""]
        if d.ip_outlook:
            L += [f"IP: {d.ip_outlook}", ""]





        route_warnings = list(getattr(d, "warnings", None) or [])
        if route_warnings:
            L += ["**Before committing to this route**", ""]
            L += [f"- {w}" for w in route_warnings] + [""]
        if d.next_steps:
            L += [f"{i}. {s_}" for i, s_ in enumerate(d.next_steps, 1)] + [""]

    if unverified:
        # Put this above everything a reader would act on. A gate that did not
        # run rejects nothing and clears nothing, and both readings are wrong.
        L += ["## Incomplete stages", "",
              "The following gates did not run. Their results are "
              "unavailable:", ""]
        L += [f"- **{k}** — {v}" for k, v in unverified.items()] + [""]

    if result.stages_skipped:
        L += ["## Stages skipped", ""]
        L += [f"- **{k}** — {v}" for k, v in result.stages_skipped.items()] + [""]
    if result.stages_run:
        L += ["## Stages run", ""] + [f"- {s}" for s in result.stages_run] + [""]





    L += _admet_summary([c for c in result.candidates if c.rank])

    a = result.provenance_audit
    if a.get("counts"):
        L += ["## Assay provenance", "",
              f"Of {a['n']} activities pulled, "
              f"{a['n_trustworthy']} ({a['trustworthy_fraction']*100:.1f}%) "
              "met the human-target activity provenance criterion. "
              "Target species, expression host, cell species and skin-function "
              "context are recorded separately.", "",
              "| class | n |", "|---|---|"]
        L += [f"| {k} | {v} |" for k, v in sorted(a["counts"].items())]
        L += [""]
        if a.get("classifier_version"):
            L += [f"Classifier v{a['classifier_version']}: "
                  f"{a.get('n_human_target_nonhuman_system', 0)} human-target activities "
                  "in nonhuman cellular systems; "
                  f"{a.get('n_human_skin_function', 0)} activities have an annotated "
                  "human skin cell functional context. These counts describe "
                  "the assay context.", ""]
        if a.get("warning"):
            L += [f"> {a['warning']}", ""]

    cc = result.control_check
    if cc:
        L += ["## Positive controls", "", f"{cc['verdict']}", ""]
        if cc.get("control_ranks"):
            L += ["| control | rank |", "|---|---|"]
            L += [f"| {n} | {r} / {cc.get('n_screened','?')} |" for n, r in cc["control_ranks"]]
            L += [""]

    # Independent of `dropped`: a compound is *held* when it ranked, survived
    # every gate, and still did not reach the buy list. A screen can have holds
    # and no drops, so this cannot live inside the block below.
    L += _held_section(result)

    dropped = getattr(result, "dropped", [])
    if dropped:
        # Every drop carries the stage that made it. Announcing the lot as
        # "dropped by the ADMET gate" sent the reader to fix an ADMET profile
        # on a compound whose ADMET was fine and whose mechanism pushed the
        # node the wrong way -- and a VDR agonist mislabelled that way is one
        # somebody proposes again next quarter with the same pChEMBL attached.
        by_stage = collections.Counter((d.get("stage") or "") for d in dropped)
        per_stage = ", ".join(f"{n} by {_stage_label(s)}"
                              for s, n in by_stage.most_common())
        noun = "compound was" if len(dropped) == 1 else "compounds were"
        L += ["## Removed before ranking", "",
              f"{len(dropped)} {noun} removed before ranking: {per_stage}. They "
              "are listed because a compound that disappears without a reason is "
              "one someone will propose again.", "",
              "| compound | stage | severity | reason |", "|---|---|---|---|"]
        for d in dropped[:_MAX_DROPPED_ROWS]:
            L.append(f"| {_cell(d.get('name') or d.get('ligand_id', '?'))} "
                     f"| {_cell(d.get('stage') or 'unrecorded')} "
                     f"| {_cell(d.get('severity'))} "
                     f"| {_cell(_readable_reason(d.get('reason')) or 'unrecorded')} |")
        if len(dropped) > _MAX_DROPPED_ROWS:
            L.append(f"| *+{len(dropped) - _MAX_DROPPED_ROWS} more* | | | "
                     "*all of them are in `result.json`* |")
        L += [""]
        # `d["reason"]` raised KeyError on a drop record that carried no
        # reason, and a KeyError here costs the whole Markdown report rather
        # than one row -- on a screen whose compounds were removed, which is
        # exactly when the removal section is the part worth reading. The
        # sourcing and skin-permeation gates are the newest producers of drop
        # records; the report renders what they send and does not require it.
        def _reason(d: dict) -> str:
            return str(d.get("reason") or "")






        groups: dict[tuple, list[dict]] = {}
        for d in dropped:
            groups.setdefault(_reason_group_key(d), []).append(d)
        if len(groups) < len(dropped):
            ranked_groups = sorted(groups.values(), key=len, reverse=True)
            shown = ranked_groups[:_MAX_REASON_GROUPS]
            n_shown = sum(len(g) for g in shown)
            L += [f"Grouped by reason ({n_shown} of {len(dropped)} compounds, "
                  f"{len(shown)} of {len(groups)} groups):", "",
                  "| reason | n | severity | stage |", "|---|---:|---|---|"]
            for g in shown:
                label = _readable_reason(_group_label([_reason(d) for d in g]))
                L.append(f"| {_cell(label or 'unrecorded')} | {len(g)} "
                         f"| {_cell(g[0].get('severity'))} "
                         f"| {_cell(g[0].get('stage') or 'unrecorded')} |")
            if len(ranked_groups) > len(shown):
                L.append(f"| *+{len(ranked_groups) - len(shown)} more groups* "
                         f"| {len(dropped) - n_shown} | | "
                         "*all of them are in `result.json`* |")
            L += [""]

    orderable = result.orderable()
    n_ranked = len([c for c in result.candidates if c.rank])
    # "carry no barring cosmetic verdict" is a claim about the eligibility gate.
    # When that gate could not run -- it raised, or its list never loaded -- the
    # sentence certifies a check that did not happen, on the one list in this
    # report that costs money to get wrong.
    blocking = [k for k in ("sourcing", "admet", "eligibility") if k in unverified]
    if blocking:
        # "0 of 40 can be bought" is a statement about the world. When the
        # vendor lookup never ran it is a statement about the pipeline, and the
        # reader cannot tell which from the number alone.
        L += ["## Functional-testing purchase priorities", "",
              f"**Not established.** This list requires the "
              f"{' and '.join(blocking)} gate(s), which did not run, so the "
              f"{len(orderable)} entries below are a lower bound on what is "
              f"available and say nothing about the other "
              f"{n_ranked - len(orderable)} ranked candidates.", ""]
    else:
        L += ["## Functional-testing purchase priorities", "",
              f"{len(orderable)} of {n_ranked} ranked "
              "candidates have reviewed direction and goal links and meet this research sourcing filter. "
              "Review formulation, exposure and local ingredient requirements before planning use.", ""]
    use_counts = collections.Counter(classify(c)['role'] for c in result.candidates if c.rank)
    L += [f"Research purpose: {use_counts[EXPLORATORY]} exploratory and "
          f"{use_counts[COUNTER_CONTROL]} counter-direction controls. "
          "Functional priorities are listed separately. "
          "See the separate exploratory, counter-control and research-availability exports.", ""]
    #: Printed once, above the first table that carries a direction column.
    _legend_written = False
    cos_caveats = list(dict.fromkeys(
        cv for c in result.candidates if c.rank
        for cv in (c.eligibility or {}).get('caveats', []) if cv))
    if cos_caveats:
        L += ["> Formulation review: "
              "The recorded eligibility caveats apply to all research purposes:", ">"]
        L += [f"> - {_cell(cv)}" for cv in cos_caveats] + [""]
    ranked_for_ip = [c for c in result.candidates if c.rank]
    if _ip_unchecked(ranked_for_ip) or any(c.ligand.props.get('ip_note') for c in ranked_for_ip):
        L += ["> `IP` records presence in the patent literature. Per-compound "
              "freedom-to-operate limitations are in `ip_note` and formulation dossier "
              "requirements are in `dossier_note` in the candidate and research-purpose CSVs.", ""]
    if orderable:
        # This is the one table captioned as a buy list and it printed no
        # direction at all: a compound confirmed to push the node the way the
        # marker needs and one nobody has ever checked read identically on the
        # list somebody spends money from. The verdict already reaches
        # candidates.csv and result.json; it belongs here first.
        if not _legend_written and any(c.direction for c in orderable[:20]):
            L += [_DIRECTION_LEGEND, ""]
            _legend_written = True
        L += ["| # | compound | composite | dock | MW | vendors | direction "
              "| IP | cosmetic |",
              "|---|---|---:|---:|---:|---:|---|---|---|"]
        for c in orderable[:20]:




            direction = (c.direction or {}).get("verdict") or (
                "not run" if not c.direction else "—")
            L.append(
                f"| {c.rank} | {_cell_or_dash(c.ligand.name or c.ligand.ligand_id)} "
                f"| {c.composite:.2f} "
                f"| {f'{c.pose.affinity:.1f}' if c.pose else '—'} "
                f"| {_cell_or_dash(c.ligand.props.get('mw'))} "
                f"| {_cell_or_dash(c.ligand.props.get('n_vendors'))} "
                f"| {_cell_or_dash(direction)} "
                f"| {_ip_cell(c.ligand.props)} "
                f"| {_cell_or_dash(_cosmetic_cell(c.eligibility))} |")
        L += [""]














        shown = orderable[:20]
        if _ip_unchecked(shown) or any(c.ligand.props.get("ip_note") for c in shown):
            L += ["> `IP` is presence in the patent literature only. The "
                  "per-compound freedom-to-operate caveat is the `ip_note` "
                  "column of `orderable.csv` and `candidates.csv`, and the "
                  "formulation dossier burden is `dossier_note`.", ""]
        if _ip_unchecked(shown):
            L += ["> `not checked` in the `IP` column means no patent lookup "
                  "ran for that compound -- generated structures carry no "
                  "sourcing annotation. It is not a finding that the chemotype "
                  "is free; a novel scaffold can still read on an existing "
                  "composition-of-matter claim.", ""]

    rows = result.top(top_n)
    if rows:
        # structure_evidence_summary reads ScoreBreakdowns, not Candidates --
        # the evidence level is a property of how a candidate was SCORED, not of
        # the candidate. Passing the wrong type here raised on every screen.
        from . import score as _score
        _mix = _score.structure_evidence_summary(list(result.breakdowns))
        if isinstance(_mix, dict):







            if _mix.get("mixed"):
                _counts = _mix.get("counts") or {}
                _shown = ", ".join(
                    f"{_score.STRUCTURE_EVIDENCE.get(k, (k, ''))[0]} {v}"
                    for k, v in _counts.items() if v)
                _note = str(_mix.get("note") or "").strip().rstrip(".")
                # score.py owns the sentence explaining what the mixture
                # does to the ranking; repeating it here in other words is how
                # the two drift apart. Print it, and only fall back to a
                # generic line when the summary carries none.
                L += ["> This column mixes structure-evidence levels"
                      + (f" ({_shown})" if _shown else "") + "."
                      + (f" {_note}." if _note else
                         " The per-row `structure_evidence` column names the "
                         "method behind each row."), ""]
        elif _mix:
            L += [str(_mix), ""]
        # The direction column is the difference between "this compound pushes
        # the node the way the marker needs" and "nobody has ever checked".
        # Survivors carried the verdict in memory and the table printed neither,
        # so a VDR agonist and a confirmed antagonist read identically here.
        if any(c.direction for c in rows) and not _legend_written:
            L += [_DIRECTION_LEGEND, ""]
            _legend_written = True
        # Which ranking group each row belongs to, and only when there is more
        # than one. The scorer ranks candidates it could evaluate separately
        # from those a budget never reached, so a mixed screen prints rank 1 at
        # composite 4.07 above rank 2 at 5.11 -- correct, and unreadable
        # without the group. The CSV has carried `rank_group` all along; this
        # table did not, so the reader with only the report saw two numbers
        # that contradict each other.
        _groups = {b.ligand_id: getattr(b, "rank_group", "")
                   for b in (result.breakdowns or [])}
        _show_group = len({g for g in _groups.values() if g}) > 1
        _ghead = " rank group |" if _show_group else ""
        _gsep = "---|" if _show_group else ""
        if _show_group:
            L += ["> This screen ranks two groups separately: rows the scorer "
                  "could evaluate, and rows a budget stopped before evaluation. "
                  "Compare composites *within* a group -- across the two they "
                  "are not on one scale, which is why a lower composite can "
                  "sit above a higher one.", ""]
        L += [f"## Top {len(rows)} candidates", "",
              "| # | compound | source |" + _ghead + " composite | pChEMBL | dock | MM/GBSA | MW | log Kp | permeation | direction | buy | IP | ADMET hazard | alerts |",
              "|---|---|---|" + _gsep + "---:|---:|---:|---:|---:|---:|---|---|---|---|---|---|"]
        for c in rows:
            s = c.skin
            direction = (c.direction or {}).get("verdict") or (
                "not run" if not c.direction else "—")







            alerts = "; ".join(
                _electrophilic_alerts(s.sensitiser_alerts if s else [])
                + [x for x in (s.stability_alerts if s else [])][:2]) or (
                    ("none" if all(skin_filter_status(c, domain) == "EVALUATED"
                                   for domain in ("sensitisation", "stability"))
                     else "unknown; not assessed") if s else "—")
            # Same pipe hazard as the buy list: `name`, `source`,
            # `permeation` and the alert strings are free text -- from a user
            # CSV in `from_file`'s case -- and an unescaped pipe in any of them
            # shifts every column to its right.
            _gcell = (f"| {_cell(_groups.get(c.ligand.ligand_id) or '—')} "
                      if _show_group else "")
            L.append(
                f"| {c.rank} | {_cell_or_dash(c.ligand.name or c.ligand.ligand_id)} "
                f"| {_cell_or_dash(c.ligand.source)} "
                + _gcell +
                f"| {c.composite:.2f} | {_cell_or_dash(c.ligand.activity.get('pchembl'))} "
                f"| {f'{c.pose.affinity:.1f}' if c.pose else '—'} "
                f"| {f'{c.gbsa.dG_bind:.1f}' if (c.gbsa and c.gbsa.converged) else '—'} "
                f"| {_cell_or_dash(c.ligand.props.get('mw'))} "
                f"| {_cell_or_dash(s.log_kp if s else None)} "
                f"| {_cell_or_dash(s.permeation_class if s else None)} "
                f"| {_cell(direction)} "
                f"| {'yes' if c.ligand.props.get('purchasable') else 'no'} "
                f"| {_ip_cell(c.ligand.props)} "
                f"| {_cell(_hazard_cell(c))} "
                f"| {_cell_or_dash(alerts)} |")
        L += ["", "The pChEMBL column summarizes the selected recorded endpoint group; "
              "inspect `measurement_groups` for the endpoint and assay context.", ""]

    L += ["---", "",
          "Permeability is Potts & Guy (1992) *Pharm Res* 9:663, "
          "`log Kp = -2.7 + 0.71 logP - 0.0061 MW`, fitted to aqueous vehicles. "
          "Sensitisation alerts follow the electrophilic reaction domains of "
          "Aptula & Roberts (2006) *Chem Res Toxicol* 19:1097. Interpret MM/GBSA "
          "affinity ranking within the recorded protocol and validation scope. "
          "Ranking accuracy and binding free energies require separate measurements.", ""]

    text = "\n".join(L)
    if path:
        _write_atomically(Path(path), text)
    return text


#: Files whose contents decide what a screen means, hashed into the manifest so
#: a result can be tied to the regulation lists and target table that produced it.
_MANIFEST_CONFIGS = ("skin_targets.yaml", "cosmetic_actives.yaml",
                     "mfds_restrictions.yaml", "kr_prohibited_ingredients.yaml",
                     "kr_cosmetic_ingredients.yaml", "kr_functional_cosmetics.yaml",








                     "structure_weight_overrides.json", "direction_evidence.json",
                     "curated_literature.json")


def _source_sha256(repo: Path) -> Optional[str]:
    'A digest of the package source, for when git cannot answer.'
    import hashlib

    pkg = repo / "skinpipe" / "skinpipe"
    if not pkg.is_dir():
        # A wheel has the package directly in site-packages; an editable
        # installation may also expose a package root rather than a checkout.
        pkg = Path(__file__).resolve().parent
        if not (pkg / "__init__.py").is_file():
            return None
    h = hashlib.sha256()
    roots = [(pkg, sorted(pkg.rglob("*.py")))]
    cfg = repo / "skinpipe" / "config"
    if not cfg.is_dir():
        cfg = pkg.parent / "config"
    if cfg.is_dir():
        roots.append((cfg, sorted(f for f in cfg.rglob("*") if f.is_file())))
    for base, files in roots:
        for f in files:
            try:
                body = f.read_bytes()
            except OSError:
                return None
            # The base name is in the key so a `config/` and a package file of
            # the same relative name cannot collide.
            h.update(f"{base.name}/{f.relative_to(base)}".encode("utf-8"))
            h.update(b"\0")
            h.update(hashlib.sha256(body).digest())
    return h.hexdigest()[:16]


def _git_state(repo: Path) -> dict:
    """Commit and dirty flag, or a reason neither could be read.

    A commit alone is not provenance if the tree had uncommitted changes when the
    screen ran, so both travel. Failures are recorded rather than raised: a
    manifest is worth writing from a tarball with no .git in it -- and when
    that happens `source_sha256` still identifies the code that ran.
    """
    import subprocess
    out: dict = {}
    for key, args in (("commit", ["rev-parse", "HEAD"]),
                      ("branch", ["rev-parse", "--abbrev-ref", "HEAD"])):
        try:
            out[key] = subprocess.run(["git", "-C", str(repo), *args],
                                      capture_output=True, text=True, timeout=15,
                                      check=True).stdout.strip()
        except Exception as exc:                              # noqa: BLE001
            out[key] = None
            out.setdefault("git_error", f"{type(exc).__name__}: {exc}")
    try:






        raw = subprocess.run(["git", "-C", str(repo), "status", "--porcelain",
                              "-z"], capture_output=True, text=True,
                             timeout=20, check=True).stdout
        out["dirty"] = bool(raw.strip("\0").strip())
        out["dirty_files"] = _porcelain_z_paths(raw)[:20]
    except Exception as exc:                                  # noqa: BLE001
        out["dirty"] = None
        out.setdefault("git_error", f"{type(exc).__name__}: {exc}")
    # Always, not only as a fallback: on a machine that has git this is the
    # cross-check that the commit and the tree agree, and on one that does not
    # it is the whole of the answer.
    out["source_sha256"] = _source_sha256(repo)
    return out


def to_manifest(result: "ScreenResult", path: Path, *,
                parameters: Optional[dict] = None,
                command: Optional[list[str]] = None) -> Path:
    "Write the run's provenance beside its results."
    import hashlib
    import platform
    import sys
    import time

    pkg_root = Path(__file__).resolve().parent            # skinpipe/skinpipe
    repo = pkg_root.parent.parent                          # the checkout
    config_dir = pkg_root.parent / "config"

    configs: dict[str, Optional[str]] = {}
    for name in _MANIFEST_CONFIGS:
        f = config_dir / name
        try:
            configs[name] = hashlib.sha256(f.read_bytes()).hexdigest()[:16]
        except OSError:
            configs[name] = None                           # absent is recorded too

    versions: dict[str, Optional[str]] = {}
    for mod in ("rdkit", "pandas", "numpy", "openmm", "vina", "meeko"):
        try:
            versions[mod] = __import__(mod).__version__
        except Exception:                                  # noqa: BLE001
            versions[mod] = None                           # optional stack absent

    manifest = {
        "screen": result.marker,
        "generated": time.strftime("%Y-%m-%d %H:%M:%S %Z"),
        "command": command if command is not None else list(sys.argv),
        "code": _git_state(repo),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "packages": versions,
        "config_sha256_16": configs,
        # The parameters the caller actually ran with. Defaults are not enough:
        # the potency floor alone changes the library size by an order of
        # magnitude, which is why `potency_floor` carries its source.
        "parameters": parameters or {},
        "potency_floor": dict(result.potency_floor or {}),
        "library_scale": dict(getattr(result, "library_scale", {}) or {}),
        "target_design": dict(result.target_design or {}),
        "n_acquired_candidates": len(result.all_candidates or []),
        "filter_profile": dict(result.filter_profile or {}),
        "structure_validation": dict(getattr(result, "structure_validation", {}) or {}),
        "candidate_exports": dict(getattr(result, "export_contract", {}) or {}),
        # The verdict a reader has to weigh against all of the above.
        "n_candidates": len(result.candidates or []),
        "gates_unverified": sorted((result.gates_unverified or {}).keys()),
        "stages_skipped": sorted((result.stages_skipped or {}).keys()),
        "trustworthy": not (result.gates_unverified or {}),
        "trustworthy_scope": "Legacy processing flag; not candidate scientific validation.",
    }
    # The same atomic write every other export uses. A direct `write_text` left
    # `{partial` behind when it failed mid-write, and because the manifest is
    # written last a broken one sits beside complete CSV and JSON -- the
    # provenance file is then the one thing that cannot be read.
    from .publication import ensure_mutable_workdir
    ensure_mutable_workdir(path.parent)
    path.parent.mkdir(parents=True, exist_ok=True)
    return _write_atomically(path, json.dumps(_json_safe(manifest), indent=2,
                                              ensure_ascii=False))
