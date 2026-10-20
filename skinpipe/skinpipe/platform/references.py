"""Assay-linked references; no potency floor or activity-count routing gate."""

from __future__ import annotations

from collections import defaultdict
from copy import deepcopy
import re

from .identity import chemical_identity, digest
from .ranking import finite, rank_candidates

POLICY = "measured-references-v6-computational-context"
BAO_FORMAT_SOURCE = "https://raw.githubusercontent.com/BioAssayOntology/BAO/fa8efc176642e59bc19094893d0ace92bdcc425a/vocabulary/bao_vocabulary_format.owl"

# These assay names collide with one-letter substitution notation: receptor
# abbreviations, cell lines, kinase/metabolite labels and a mouse strain. A
# bare label is not construct evidence. Explicit mutant/HGVS annotations and
# structured variant fields still take precedence, including for these tokens.
_AMBIGUOUS_SUBSTITUTION_LABELS = frozenset(
    {"A2A", "H3R", "D2R", "T47D", "N2A", "S6K", "S1P", "G6P", "C3H"}
)
_SUBSTITUTION = re.compile(r"\b[ACDEFGHIKLMNPQRSTVWY][1-9][0-9]{0,4}[ACDEFGHIKLMNPQRSTVWY]\b")
_COMPUTATIONAL_ASSAY = re.compile(
    r"\b(?:molecular\s+docking|docking\s+(?:score|energy|affinity)|"
    r"MM[-\s]?(?:GBSA|PBSA)|in\s+silico|computational(?:ly)?|"
    r"predicted\s+(?:binding|affinity|activity))\b", re.IGNORECASE,
)


def _variant_text_evidence(description):
    tokens = _SUBSTITUTION.findall(description)
    explicit = bool(re.search(
        r"\b(?:mutant|mutants|mutation|mutations|mutated|substitution|substitutions)\b",
        description, re.IGNORECASE,
    ))
    # HGVS protein notation disambiguates even a short token such as p.T47D.
    hgvs = bool(re.search(
        r"\bp\.(?:[ACDEFGHIKLMNPQRSTVWY]|Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val)"
        r"[1-9][0-9]{0,4}(?:[ACDEFGHIKLMNPQRSTVWY]|Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val)\b",
        description,
    ))
    ambiguous = sorted(set(tokens) & _AMBIGUOUS_SUBSTITUTION_LABELS)
    substitutions = sorted(set(tokens) - _AMBIGUOUS_SUBSTITUTION_LABELS)
    return {
        "flag": explicit or hgvs or bool(substitutions),
        "explicit_variant_language": explicit,
        "hgvs_protein_annotation": hgvs,
        "substitution_tokens": substitutions,
        "ambiguous_assay_labels": ambiguous,
    }


def measurement(raw, descriptor, assay):
    """Preserve original records. Binding and biochemical function stay distinct."""
    record = deepcopy(raw)
    target = raw.get("target_chembl_id")
    confidence = assay.get("confidence_score", assay.get("assay_confidence_score"))
    try:
        confidence = int(confidence)
    except (ValueError, TypeError):
        confidence = None
    direct = confidence == (7 if descriptor["kind"] == "complex" else 9)
    organism = assay.get("assay_organism", assay.get("organism")) or raw.get(
        "target_organism"
    )
    same_species = organism == descriptor.get("organism")
    assay_type = assay.get("assay_type", raw.get("assay_type"))
    cell = assay.get("assay_cell_type") or assay.get("assay_tissue")
    # BAO_0000357 is the ontology's single protein format, a descendant of
    # biochemical format. Unspecified/cell-free/whole-cell formats are not
    # promoted simply because the ChEMBL assay type is functional.
    biochemical = (
        assay.get("bao_format") == "BAO_0000357"
        and descriptor.get("kind", "protein") == "protein"
    )
    context = (
        "binding"
        if assay_type == "B"
        else "biochemical_function"
        if assay_type == "F" and biochemical and not cell
        else "cellular_or_unspecified_function"
        if assay_type == "F"
        else "other"
    )
    value = finite(raw.get("standard_value"))
    if value is None:
        # ChEMBL JSON numeric values are decimal strings; parsing is confined
        # to this source adapter, not the generic ranking interface.
        try:
            value = finite(float(raw.get("standard_value")))
        except (ValueError, TypeError, OverflowError):
            value = None
    relation, unit, endpoint = (
        raw.get("standard_relation"),
        raw.get("standard_units"),
        raw.get("standard_type"),
    )
    comment = str(raw.get("activity_comment") or "").strip().lower()
    negative = comment in {
        "inactive",
        "not active",
        "no activity",
        "inconclusive",
        "uncertain",
        "conflicting",
    }
    validity = raw.get("data_validity_comment")
    identity = chemical_identity(raw.get("canonical_smiles"))
    reasons = []
    # A protein accession does not establish the assay's mutant construct.
    # Keep variant observations, but do not seed a canonical-target search or
    # rank them as its direct measurements without explicit construct matching.
    variant_id = assay.get("variant_id", raw.get("assay_variant_id"))
    variant_sequence = assay.get("variant_sequence") or raw.get("variant_sequence")
    description = str(assay.get("description", raw.get("assay_description")) or "")
    computational = bool(_COMPUTATIONAL_ASSAY.search(description))
    if computational:
        reasons.append("computational_assay_not_experimental_reference")
    text_evidence = _variant_text_evidence(description)
    variant_text = text_evidence["flag"]
    if variant_id not in (None, "") or variant_sequence or variant_text:
        reasons.append("assay_variant_context_unverified")
    if target != descriptor.get("chembl_target_id") or not direct:
        reasons.append("target_mapping_not_direct")
    if not same_species:
        reasons.append("species_unverified_or_conflicting")
    if context not in {"binding", "biochemical_function"}:
        reasons.append("assay_context_not_binding_or_verified_biochemical")
    if value is None or not unit or not endpoint or relation != "=":
        reasons.append("measurement_missing_or_censored")
    if (
        endpoint in {"Ki", "Kd", "IC50", "EC50", "AC50", "XC50", "Potency"}
        and value is not None
        and value <= 0
    ):
        reasons.append("nonpositive_concentration_endpoint")
    if validity not in {None, "", "Manually validated"}:
        reasons.append("source_validity_warning")
    if negative:
        reasons.append("negative_conflicting_or_uncertain")
    if identity["structure_status"] != "valid" or identity["parent_ambiguous"]:
        reasons.append("invalid_reference_structure")
    assay_id = str(
        raw.get("assay_chembl_id")
        or assay.get("assay_chembl_id")
        or raw.get("assay_id")
        or ""
    )
    if not assay_id:
        reasons.append("assay_identifier_missing")
    # Endpoint semantics are explicit. An unrecognized endpoint remains a
    # reference observation but has no invented higher/lower ranking direction.
    lower = endpoint in {
        "Ki",
        "Kd",
        "IC50",
        "EC50",
        "AC50",
        "XC50",
        "Potency",
    } and unit in {"nM", "uM", "µM", "mM", "M"}
    record.update(
        target_id=descriptor["target_id"],
        species=descriptor.get("organism"),
        assay_id=assay_id,
        assay_type=context,
        endpoint=endpoint,
        unit=unit,
        value=value,
        relation=relation,
        is_measured=value is not None and not computational,
        observation_scope="computational_prediction" if computational else "assay_record",
        higher_is_better=False if lower else None,
        action="unspecified",
        assay_context=digest(assay),
        original=deepcopy(raw),
        construct_context={
            "variant_id": variant_id,
            "variant_sequence": deepcopy(variant_sequence),
            "variant_text_flag": variant_text,
            "variant_text_evidence": text_evidence,
            "matching_status": "unverified_variant" if "assay_variant_context_unverified" in reasons else "no_explicit_variant_annotation",
        },
        biochemical_format_provenance={
            "bao_format": assay.get("bao_format"),
            "source": BAO_FORMAT_SOURCE,
        }
        if assay.get("bao_format") == "BAO_0000357"
        else None,
        reference_eligible=not reasons,
        reference_hold_reasons=reasons,
        reference_label="computational_prediction"
        if computational
        else "assay_reported_active"
        if comment == "active"
        else "measured_reference",
        ranking_eligible=not reasons and lower,
        identity=identity,
    )
    return record


def freeze_references(records, descriptor, *, maximum=5):
    """Full eligible population is retained; only computational seeds are capped."""
    groups, conflicts = defaultdict(list), set()
    for record in records:
        parent = record.get("identity", {}).get("parent_id")
        if parent:
            groups[parent].append(record)
    for parent, rows in groups.items():
        if any(
            "negative_conflicting_or_uncertain" in r.get("reference_hold_reasons", [])
            for r in rows
        ) and any(r.get("reference_eligible") for r in rows):
            conflicts.add(parent)
    eligible = [
        r
        for r in records
        if r.get("reference_eligible") and r["identity"]["parent_id"] not in conflicts
    ]
    by_exact = defaultdict(list)
    for r in eligible:
        by_exact[r["identity"]["exact_structure_id"]].append(r)
    rows = [
        {
            "candidate_id": exact,
            "target_id": descriptor["target_id"],
            "smiles": rs[0]["canonical_smiles"],
            "measurements": rs,
            "reference_label": "assay_reported_active"
            if any(r["reference_label"] == "assay_reported_active" for r in rs)
            else "measured_reference",
            "identity": rs[0]["identity"],
        }
        for exact, rs in by_exact.items()
    ]
    ranked = rank_candidates(rows)
    ranked.sort(
        key=lambda r: (
            r["reference_label"] != "assay_reported_active",
            -r["priority_score_v2"],
            r["candidate_id"],
        )
    )
    seeds, cores, parents = [], set(), set()
    for diverse in (True, False):
        for row in ranked:
            if len(seeds) >= maximum:
                break
            if row["identity"]["parent_id"] in parents or (
                diverse and row["identity"]["scaffold_id"] in cores
            ):
                continue
            seeds.append(row)
            cores.add(row["identity"]["scaffold_id"])
            parents.add(row["identity"]["parent_id"])
    result = {
        "policy": POLICY,
        "target_id": descriptor["target_id"],
        "population": ranked,
        "computational_seeds": seeds,
        "observations": records,
        "conflicting_parents": sorted(conflicts),
        "active_threshold": None,
        "activity_count_gate": None,
        "interpretation": "Measured references do not prove analogue activity or positive-control status",
    }
    return {**result, "sha256": digest(result)}
