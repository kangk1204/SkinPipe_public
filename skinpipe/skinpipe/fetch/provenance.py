"""Versioned assay provenance; target, host and skin context are separate.

Target organism identifies the database target, not necessarily measured
material. Human-target activity inclusion does not verify binding, direction,
skin efficacy, safety or clinical applicability.
"""
from __future__ import annotations

import logging
import re
from dataclasses import dataclass
from enum import Enum
from typing import Optional

log = logging.getLogger(__name__)
CLASSIFIER_VERSION = 3
CONTEXT_FIELDS = (
    'provenance_classifier_version', 'target_species', 'expression_host',
    'assay_cell_species', 'assay_context', 'skin_relevance', 'target_directness',
    'human_target_activity_eligible', 'human_skin_function_eligible',
    'provenance_conflicts', 'provenance_evidence',
)

class Provenance(str, Enum):
    HUMAN_PROTEIN = 'human_protein'
    HUMAN_CELL = 'human_cell'
    SURROGATE_SPECIES = 'surrogate'
    NONHUMAN_CELL = 'nonhuman_cell'
    UNKNOWN_ORIGIN = 'unknown_origin'
    UNSPECIFIED = 'unspecified'

SURROGATES = {
    r'mushroom|agaricus': ('mushroom', 'TYR'),
    r'clostridium|bacterial collagenase': ('Clostridium', 'MMP1'),
    r'porcine pancreatic elastase|\bppe\b|pancreatic elastase': ('porcine pancreas', 'ELANE'),
    r'bovine testicular|testicular hyaluronidase|bovine hyaluronidase': ('bovine testis', 'HYAL1'),
    r'\bbovine\b|\bcalf\b': ('bovine', ''),
    r'\bporcine\b|\bpig\b|\bswine\b': ('porcine', ''),
    r'\bovine\b|\bsheep\b': ('ovine', ''),
    r'\bequine\b|\bhorse\b': ('equine', ''),
    r'\brabbit\b': ('rabbit', ''),
    r'\byeast\b|saccharomyces': ('yeast', ''),
    r'\be\.? ?coli\b|escherichia': ('E. coli', ''),
    r'streptomyces': ('Streptomyces', ''),
}
NONHUMAN_CELLS = {
    r'\bb-?16': ('B16 (murine melanoma)', 'Mus musculus'),
    r'\bmelan-a\b': ('melan-a (murine melanocyte)', 'Mus musculus'),
    r'\b3t3\b|nih-?3t3': ('3T3 (murine fibroblast)', 'Mus musculus'),
    r'raw ?264': ('RAW264.7 (murine macrophage)', 'Mus musculus'),
    r'\bcho\b|chinese hamster ovary': ('CHO (hamster)', 'Cricetulus griseus'),
    r'\bmdck\b': ('MDCK (canine)', 'Canis lupus familiaris'),
    r'\bhepa ?1c1c7\b': ('Hepa1c1c7 (murine)', 'Mus musculus'),
    r'\brat\b|\brattus norvegicus\b': ('rat', 'Rattus norvegicus'),
    r'\bmurine\b|\bmouse\b|\bmus musculus\b': ('mouse', 'Mus musculus'),
}
HUMAN_MARKERS = (
    r'\bhuman\b', r'\bhomo sapiens\b', r'hek ?293', r'\bhela\b', r'\bhacat\b',
    r'\bmnt-?1\b', r'\bsk-?mel', r'\bnhek\b', r'\bnhdf\b', r'\bepiderm\b',
    r'\bhepg2\b', r'\bthp-?1\b', r'\bjurkat\b', r'\bu937\b', r'\ba375\b',
    r'\bmcf-?7\b', r'\bhuvec\b',
)
_HUMAN_ORGANISMS = frozenset({'homo sapiens', 'human', 'h. sapiens', 'homo sapiens sapiens'})
_HUMAN_RE = re.compile('|'.join(HUMAN_MARKERS), re.I)
_HUMAN_CELL_RE = re.compile('|'.join(HUMAN_MARKERS[2:]), re.I)
_UNKNOWN_RE = re.compile(r'unknown origin|unspecified origin|origin not', re.I)
_CELL_HINT = re.compile(r'\bcells?\b|\bcell.line\b|culture|melanocyte|keratinocyte|fibroblast', re.I)
_SKIN_RE = re.compile(r'\bskin\b|dermal|keratinocyte|melanocyte|hair|\bhacat\b|\bnhek\b|\bnhdf\b|\bmnt-?1\b|\bsk-?mel|\bepiderm\b', re.I)
_EXPRESSION_RE = re.compile(r'express(?:ed|ing|ion)|transfect(?:ed|ion)|recombinant|produced', re.I)
_HOSTS = (
    (r'\be\.? ?coli\b|escherichia coli', 'Escherichia coli', ''),
    (r'\bcho\b|chinese hamster ovary', 'CHO', 'Cricetulus griseus'),
    (r'\byeast\b|saccharomyces', 'yeast', ''),
    (r'hek ?293', 'HEK293', 'Homo sapiens'),
    (r'sf-?9|spodoptera', 'Sf9', 'Spodoptera frugiperda'),
)
_SURROGATE_RES = [(re.compile(p, re.I), v) for p, v in SURROGATES.items()]
_NONHUMAN_CELL_RES = [(re.compile(p, re.I), v) for p, v in NONHUMAN_CELLS.items()]

@dataclass(frozen=True)
class ProvenanceCall:
    provenance: Provenance
    detail: str = ''
    trustworthy: bool = True  # legacy name: human-target activity inclusion
    basis: str = 'description'
    target_species: str = 'unknown'
    expression_host: str = 'unknown'
    assay_cell_species: str = 'not_evaluated'
    assay_context: str = 'unknown'
    skin_relevance: str = 'not_evaluated'
    target_directness: str = 'unknown'
    human_skin_function_eligible: bool = False
    conflicts: tuple[str, ...] = ()
    evidence: tuple[tuple[str, str], ...] = ()
    classifier_version: int = CLASSIFIER_VERSION

    @property
    def label(self) -> str:
        return self.provenance.value

    def record_fields(self) -> dict:
        return {
            'provenance': self.label, 'provenance_detail': self.detail,
            'provenance_trustworthy': self.trustworthy, 'provenance_basis': self.basis,
            'provenance_classifier_version': self.classifier_version,
            'target_species': self.target_species, 'expression_host': self.expression_host,
            'assay_cell_species': self.assay_cell_species, 'assay_context': self.assay_context,
            'skin_relevance': self.skin_relevance, 'target_directness': self.target_directness,
            'human_target_activity_eligible': self.trustworthy,
            'human_skin_function_eligible': self.human_skin_function_eligible,
            'provenance_conflicts': list(self.conflicts),
            'provenance_evidence': [{'field': k, 'value': v} for k, v in self.evidence],
        }


def _human(org: str) -> bool:
    return org.casefold() in _HUMAN_ORGANISMS


def _human_target_pattern(gene: Optional[str]) -> re.Pattern:
    if gene:
        # A human cell elsewhere in the text does not identify the construct.
        between = r'(?:(?!\bcells?\b|\bskin\b|[.;]).){0,70}'
        return re.compile(r'\bhuman\b' + between + r'\b' + re.escape(gene) + r'\b', re.I)
    return re.compile(r'\brecombinant human\b|\bhuman (?:enzyme|protein|receptor)\b', re.I)


def classify(description: Optional[str], assay_organism: Optional[str] = None,
             target_gene: Optional[str] = None, *, target_organism: Optional[str] = None,
             assay_type: Optional[str] = None, assay_cell_type: Optional[str] = None,
             assay_tissue: Optional[str] = None) -> ProvenanceCall:
    """Separate measured material from expression system and functional context.

    Explicit surrogate material takes precedence. A database human target
    alone cannot rescue missing or nonhuman assay provenance. Description
    annotations and assay classifications are not proof of physical binding.
    """
    desc = str(description or '').strip()
    org = ' '.join(str(assay_organism or '').split())
    target_org = ' '.join(str(target_organism or '').split())
    text = ' '.join((desc, str(assay_cell_type or ''), str(assay_tissue or ''))).strip()
    constructs = list(_human_target_pattern(target_gene).finditer(desc))
    explicit = bool(constructs)
    conflicts = []
    evidence = tuple((k, str(v)) for k, v in (
        ('assay_description', description), ('assay_organism', assay_organism),
        ('target_organism', target_organism), ('assay_type', assay_type),
        ('assay_cell_type', assay_cell_type), ('assay_tissue', assay_tissue),
    ) if v is not None and str(v).strip())
    host, host_cell, host_spans = 'unknown', '', []
    for pattern, name, species in _HOSTS:
        for match in re.finditer(pattern, desc, re.I):
            linked = False
            for construct in constructs:
                start = min(match.start(), construct.start())
                end = max(match.end(), construct.end())
                between = desc[construct.end():match.start()] if construct.end() <= match.start() else desc[match.end():construct.start()]
                clause = re.sub(r'\be\. ?coli\b|\bh\. ?sapiens\b', 'organism', between, flags=re.I)
                # A host and a different expressed protein elsewhere in the
                # description must not be assigned to this human construct.
                linked |= (end - start <= 180 and not re.search(r'[.;]', clause)
                           and bool(_EXPRESSION_RE.search(between)))
            if linked:
                if host != 'unknown' and host != name:
                    conflicts.append('multiple_expression_hosts')
                host, host_cell = name, species
                host_spans.append((match.start(), match.end()))
    material_text = desc
    for start, end in sorted(host_spans, reverse=True):
        material_text = material_text[:start] + ' ' + material_text[end:]
    surrogate = next((v for rx, v in _SURROGATE_RES if rx.search(material_text)), None)
    cell = next((v for rx, v in _NONHUMAN_CELL_RES if rx.search(text)), None)
    human_cells = bool(_HUMAN_CELL_RE.search(text)) or bool(_CELL_HINT.search(text) and _HUMAN_RE.search(text))
    cellular = bool(cell or human_cells or _CELL_HINT.search(text))
    cell_species = host_cell if host_cell else cell[1] if cell else 'Homo sapiens' if human_cells else org if cellular and org else 'unknown' if cellular else 'not_applicable'
    if host in ('Escherichia coli', 'yeast'):
        # Expression in microbial cells does not make a purified-protein
        # measurement an assay in human cells.
        cell_species = 'not_applicable'
    context = 'recombinant_cellular' if host_cell else 'recombinant_biochemical' if host != 'unknown' else 'cellular' if cellular else 'biochemical_or_unspecified'
    species = 'unknown'
    kind, trust, detail, basis = Provenance.UNSPECIFIED, False, 'source not stated in the assay description', 'description'
    if not text and not org:
        detail = 'no assay description'
    elif surrogate:
        name, usual_gene = surrogate
        species = org if org and not _human(org) else name
        kind, detail = Provenance.SURROGATE_SPECIES, f'{name} enzyme'
        if target_gene and usual_gene.upper() == target_gene.upper():
            detail += f' commonly substituted for human {usual_gene}'
        if org and _human(org):
            conflicts.append('human_assay_organism_vs_nonhuman_material_description')
    elif org and not _human(org) and not (explicit and host != 'unknown'):
        species = org
        kind, detail, basis = Provenance.SURROGATE_SPECIES, f'assay organism is {org}, not human', 'assay_organism'
    elif _UNKNOWN_RE.search(text):
        kind, detail = Provenance.UNKNOWN_ORIGIN, 'description states the enzyme origin is unknown'
        if org and _human(org):
            conflicts.append('human_assay_organism_vs_unknown_material_description')
    elif explicit and host != 'unknown':
        species, trust = 'Homo sapiens', True
        kind = Provenance.NONHUMAN_CELL if host_cell and host_cell != 'Homo sapiens' else Provenance.HUMAN_CELL if host_cell else Provenance.HUMAN_PROTEIN
        detail, basis = f'explicit human {target_gene or "construct"} in {host}; host is separate from target species', 'description_construct'
        if org and not _human(org):
            aliases = {'Escherichia coli': ('escherichia coli', 'e. coli'),
                       'CHO': ('cricetulus griseus', 'chinese hamster', 'chinese hamster ovary', 'cho'),
                       'yeast': ('saccharomyces cerevisiae', 'yeast'),
                       'Sf9': ('spodoptera frugiperda', 'sf9')}
            if org.casefold() not in aliases.get(host, ()):
                conflicts.append('assay_organism_disagrees_with_expression_host')
                trust = False
    elif cell:
        kind, detail, species = Provenance.NONHUMAN_CELL, cell[0], cell[1]
        if org and _human(org):
            conflicts.append('human_assay_organism_vs_nonhuman_cell_description')
    elif (org and _human(org)) or _HUMAN_RE.search(text):
        kind = Provenance.HUMAN_CELL if cellular else Provenance.HUMAN_PROTEIN
        species, trust = 'Homo sapiens', True
        basis = 'assay_organism' if org and _human(org) else 'description'
        detail = 'human source confirmed by the assay record' if basis == 'assay_organism' else 'human source from the assay description; the assay organism field was not available'
    if target_org and species != 'unknown' and (_human(target_org) != _human(species)):
        conflicts.append('target_record_organism_differs_from_measured_material')
        if trust and not _human(target_org):
            trust = False
    if 'multiple_expression_hosts' in conflicts:
        trust = False
    skin = 'not_evaluated'
    if _SKIN_RE.search(text):
        skin = 'human_skin_cell_context' if cellular and cell_species == 'Homo sapiens' else 'nonhuman_skin_cell_context' if cellular and cell_species not in ('unknown', 'not_applicable') else 'skin_context_unverified'
    atype = str(assay_type or '').upper()
    directness = 'binding_assay_reported' if atype == 'B' else 'target_function_reported' if atype == 'F' and explicit and host != 'unknown' else 'cellular_function_reported' if atype == 'F' and cellular else 'functional_assay_reported' if atype == 'F' else 'unknown'
    skin_function = trust and skin == 'human_skin_cell_context' and atype == 'F' and not conflicts
    return ProvenanceCall(kind, detail, trust, basis, species, host, cell_species,
                          context, skin, directness, skin_function,
                          tuple(dict.fromkeys(conflicts)), evidence)


def reclassify(records: list[dict], target_gene: Optional[str] = None) -> list[dict]:
    """Explicit lossless legacy migration; never mutates or writes raw inputs."""
    out = []
    for record in records:
        call = classify(record.get('assay_description'), record.get('assay_organism'),
                        target_gene, target_organism=record.get('target_organism'),
                        assay_type=record.get('assay_type'),
                        assay_cell_type=record.get('assay_cell_type'),
                        assay_tissue=record.get('assay_tissue'))
        previous = {k: v for k, v in record.items() if k.startswith('provenance') or k in CONTEXT_FIELDS}
        out.append(record | {'previous_provenance': previous} | call.record_fields())
    return out


def _call_for(record: dict, target_gene: Optional[str], desc_key: str,
              organism_key: str) -> ProvenanceCall:
    """Restore a recorded call. Legacy reclassification is an explicit action."""
    stored = record.get('provenance')
    if stored is not None:
        try:
            ev = record.get('provenance_evidence') or []
            return ProvenanceCall(
                Provenance(stored), record.get('provenance_detail', ''),
                bool(record.get('provenance_trustworthy')),
                record.get('provenance_basis') or 'unrecorded',
                record.get('target_species') or 'unknown',
                record.get('expression_host') or 'unknown',
                record.get('assay_cell_species') or 'not_evaluated',
                record.get('assay_context') or 'unknown',
                record.get('skin_relevance') or 'not_evaluated',
                record.get('target_directness') or 'unknown',
                bool(record.get('human_skin_function_eligible')),
                tuple(record.get('provenance_conflicts') or []),
                tuple((e['field'], e['value']) for e in ev),
                record.get('provenance_classifier_version', 0))
        except ValueError:
            log.warning('unrecognised stored provenance label %r; reclassifying', stored)
    return classify(record.get(desc_key), record.get(organism_key), target_gene,
                    target_organism=record.get('target_organism'),
                    assay_type=record.get('assay_type'),
                    assay_cell_type=record.get('assay_cell_type'),
                    assay_tissue=record.get('assay_tissue'))


def audit(records: list[dict], target_gene: Optional[str] = None,
          desc_key: str = 'assay_description', organism_key: str = 'assay_organism') -> dict:
    calls = [_call_for(r, target_gene, desc_key, organism_key) for r in records]
    counts: dict[str, int] = {}
    for c in calls:
        counts[c.label] = counts.get(c.label, 0) + 1
    n_trust = sum(c.trustworthy for c in calls)
    frac = n_trust / len(calls) if calls else 0.0
    warning = '' if not calls or frac >= 0.95 else (
        f'{(1-frac)*100:.0f}% of activities are NOT confirmed human-target material; '
        'review selection and individual potencies. Host and skin function are separate.')
    return {'n': len(calls), 'counts': counts, 'n_trustworthy': n_trust,
            'trustworthy_fraction': round(frac, 3), 'warning': warning, 'calls': calls,
            'classifier_version': CLASSIFIER_VERSION,
            'n_human_skin_function': sum(c.human_skin_function_eligible for c in calls),
            'n_human_target_nonhuman_system': sum(c.trustworthy and c.provenance is Provenance.NONHUMAN_CELL for c in calls),
            'n_conflicted': sum(bool(c.conflicts) for c in calls),
            'classification_versions': sorted({c.classifier_version for c in calls})}


def filter_trustworthy(records: list[dict], target_gene: Optional[str] = None,
                       desc_key: str = 'assay_description', organism_key: str = 'assay_organism',
                       *, purpose: str = 'human_target_activity') -> list[dict]:
    if purpose not in ('human_target_activity', 'human_skin_function'):
        raise ValueError(f'unknown provenance inclusion purpose: {purpose}')
    kept = []
    for record in records:
        call = _call_for(record, target_gene, desc_key, organism_key)
        include = call.trustworthy if purpose == 'human_target_activity' else call.human_skin_function_eligible
        if include:
            kept.append(record | {'_provenance': call.label, '_provenance_detail': call.detail})
    return kept
