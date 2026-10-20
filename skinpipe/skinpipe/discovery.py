"""Versioned exploratory triage, over the complete acquired candidate pool.

The base score prioritises diverse, assessable chemistry for experiments. It is
not an activity predictor. Measured potency and unvalidated docking energies do
not award bonus points; their observations remain separate sortable columns.
"""
from __future__ import annotations
from collections import Counter
from copy import deepcopy
import math
from fractions import Fraction
import re

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from .hop_sanity import sa_assessment
from .candidate_use import classify
from .types import hazard_severity_of

POLICY = "exploration-triage-v2.3"
WEIGHTS = {"chemical_diversity": .20, "synthesis": .10,
           "structural_alerts": .20, "topical": .30, "stability": .20}
NATURAL_PRODUCT_BONUS = 5.0


def natural_product_assessment(row):
    """Credit occurrence evidence, never a ChEMBL flag or a natural parent."""
    sources = row.get("natural_product_sources")
    sources = sources if isinstance(sources, list) else []
    evidence = []
    if not row.get("derivation"):
        for source in sources:
            if not isinstance(source, dict) or source.get("value") is not True:
                continue
            urls = source.get("source_urls")
            urls = [u for u in urls if isinstance(u, str) and u.startswith(("https://", "http://"))] if isinstance(urls, list) else []
            basis = source.get("basis")
            if (source.get("source") == "curated" and isinstance(basis, str)
                    and basis.strip() and basis.lower().strip() != "unreviewed" and urls):
                evidence.append(deepcopy(source))
            elif (source.get("source") == "lotus" and source.get("match") == "exact"
                    and isinstance(source.get("organisms"), list) and source["organisms"]
                    and re.fullmatch(r"[0-9a-f]{64}", str(source.get("archive_sha256", "")))):
                smiles = row.get("smiles")
                mol = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
                if mol is not None and Chem.MolToInchiKey(mol) == source.get("inchikey"):
                    evidence.append(deepcopy(source))
    return {"status": "confirmed" if evidence else "unconfirmed",
            "bonus_points": NATURAL_PRODUCT_BONUS if evidence else 0.0,
            "evidence": evidence,
            "scope": "natural occurrence evidence; not safety or efficacy"}


def scored_terms(terms, natural_product):
    base = float(round(sum((Fraction(str(WEIGHTS[k])) * Fraction(str(v))
                            for k, v in terms.items()), Fraction()), 2))
    return {"discovery_base_score": base,
            "discovery_score": round(min(100, base + natural_product["bonus_points"]), 2),
            "natural_product_bonus": natural_product["bonus_points"],
            "natural_product_bonus_policy": {"points": NATURAL_PRODUCT_BONUS, "cap": 100,
                                              "interpretation": "additive percentage points"},
            "natural_product_assessment": natural_product}


def annotate_natural_occurrences(records, *, archive=None, archive_sha256=None):
    """Attach exact LOTUS occurrence matches; preserve original NP annotations."""
    from .fetch import lotus_local
    keys, cache = {}, {}
    for index, row in enumerate(records):
        smiles = row.get("smiles")
        if not isinstance(smiles, str) or row.get("derivation"):
            continue
        if smiles not in cache:
            mol = Chem.MolFromSmiles(smiles)
            cache[smiles] = Chem.MolToInchiKey(mol) if mol is not None else None
        if cache[smiles]:
            keys[index] = cache[smiles]
    matches = lotus_local.lookup(set(keys.values()), archive=archive)
    output = deepcopy(records)
    for index, key in keys.items():
        hit = matches.get(key)
        if not hit or not hit.get("organisms"):
            continue
        sources = output[index].get("natural_product_sources")
        sources = deepcopy(sources) if isinstance(sources, list) else []
        source = {"source": "lotus", "value": True, "match": "exact", "inchikey": key,
                  "organisms": deepcopy(hit["organisms"]), "archive_sha256": archive_sha256,
                  "basis": "exact full InChIKey occurrence match"}
        if source not in sources:
            sources.append(source)
        output[index]["natural_product_sources"] = sources
    return output


def reweight_records(records):
    """Reweight verified source terms without changing the source observations.

    Missing or malformed terms contribute neutral 50 to the heuristic and are
    disclosed separately. They never become observed measurements. The source
    score and policy are retained before the first policy transition.
    """
    output = []
    for original in records:
        row = deepcopy(original)
        if row.get("discovery_policy") != POLICY:
            row.setdefault("previous_discovery", {key: deepcopy(row.get(key)) for key in (
                "discovery_score", "discovery_policy", "discovery_weights", "discovery_terms")})
        source_terms = row.get("discovery_terms")
        source_terms = source_terms if isinstance(source_terms, dict) else {}
        supplied_missing = row.get("discovery_missing")
        missing = {v for v in supplied_missing if isinstance(v, str) and v in WEIGHTS} if isinstance(supplied_missing, list) else set()
        terms = {}
        for name in WEIGHTS:
            value = source_terms.get(name)
            if (name in missing or isinstance(value, bool) or not isinstance(value, (int, float))
                    or not 0 <= value <= 100 or not math.isfinite(value)):
                terms[name] = 50
                missing.add(name)
            else:
                terms[name] = value
        row.update(scored_terms(terms, natural_product_assessment(row)))
        row.update(discovery_policy=POLICY, discovery_weights=dict(WEIGHTS),
                   discovery_scoring_terms=terms, discovery_missing=sorted(missing),
                   discovery_coverage=round(sum(v for k, v in WEIGHTS.items() if k not in missing), 3),
                   score_interpretation="experimental triage heuristic; not efficacy or binding probability")
        output.append(row)
    return output


def scaffold(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    core = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=True)
    return core or "acyclic:" + Chem.MolToSmiles(mol)


def priorities(candidates):
    from .score import safety_score, stability_score, _electrophilic_alerts, colour_of, _nonnegative_count, _finite, skin_filter_status
    from .filters.safety import SEVERITY_WEIGHT
    cores = {c.ligand.ligand_id: scaffold(c.ligand.smiles) for c in candidates}
    counts = Counter(cores.values())
    out = {}
    for c in candidates:
        lig, skin = c.ligand, c.skin
        missing = []
        sa = lig.props.get("sa_score")
        sa_method = "input property; method unverified"
        sa_fallback_reason = None
        if sa is None:
            assessment = sa_assessment(lig.smiles)
            sa, sa_method = assessment["score"], assessment["method"]
            sa_fallback_reason = assessment["fallback_reason"]
        # Bounds first: math.isfinite can overflow on an arbitrarily large int.
        if isinstance(sa, bool) or not isinstance(sa, (int, float)) or not 1 <= sa <= 10 or not math.isfinite(sa):
            sa = None
        if sa is None:
            missing.append("synthesis")
        core = cores[lig.ligand_id]
        if skin is None:
            missing += ["structural_alerts", "topical", "stability"]
        elif _nonnegative_count(skin.n_pains) is None:
            missing.append("structural_alerts")
        if skin_filter_status(c, "sensitisation") != "EVALUATED" and "structural_alerts" not in missing:
            missing.append("structural_alerts")
        if skin_filter_status(c, "stability") != "EVALUATED" and "stability" not in missing:
            missing.append("stability")
        topical = skin.topical_score if skin else None
        topical_valid = (isinstance(topical, (int, float)) and _finite(topical) is not None
                         and 0 <= topical <= 10)
        if skin is not None and not topical_valid:
            missing.append("topical")
        burden = sum(SEVERITY_WEIGHT.get(a.split(":")[0], 0) for a in
                     _electrophilic_alerts(skin.sensitiser_alerts if skin else []))
        terms = {"chemical_diversity": 100 / math.sqrt(counts[core]) if core else 0,
                 "synthesis": max(0, min(100, (10 - sa) / 9 * 100)) if sa is not None else 50,
                 "structural_alerts": 10 * safety_score(burden, skin.n_pains) if skin else 50,
                 "topical": (10 * topical if topical_valid else 0) if skin else 50,
                 "stability": 10 * stability_score(len(skin.stability_alerts), colour_of(c)) if skin else 50}
        scoring_terms = {k: 50 if k in missing else v for k, v in terms.items()}
        out[lig.ligand_id] = {**scored_terms(scoring_terms, natural_product_assessment({
                                 "smiles": lig.smiles,
                                 "natural_product_sources": lig.natural_product_sources,
                                 "derivation": lig.derivation})),
                             "discovery_terms": terms, "discovery_weights": dict(WEIGHTS),
                             "discovery_scoring_terms": scoring_terms,
                             "discovery_policy": POLICY, "discovery_missing": missing,
                             "discovery_coverage": round(sum(WEIGHTS[k] for k in terms if k not in missing), 3),
                             "scaffold_smiles": core, "scaffold_population": counts[core],
                             "sa_score": sa, "sa_method": sa_method,
                             "sa_fallback_reason": sa_fallback_reason,
                             "sa_status": "EVALUATED" if sa is not None else "NOT_EVALUATED",
                             "score_interpretation": "experimental triage heuristic; not efficacy or binding probability"}
    return out


def rows(result, row_builder):
    population = result.all_candidates or result.candidates
    scores = priorities(population)
    drops = {}
    for drop in result.dropped:
        drops.setdefault(drop["ligand_id"], []).append(drop)
    breakdown = {b.ligand_id: b for b in result.breakdowns}
    output = []
    from .pipeline import _direction_of
    from dataclasses import replace
    for candidate in population:
        direction = candidate.direction if candidate.direction is not None else _direction_of(candidate.ligand)
        view = replace(candidate, direction=direction)
        item = row_builder(view, breakdown.get(candidate.ligand.ligand_id))
        item.update(scores[candidate.ligand.ligand_id])
        item["filter_decisions"] = drops.get(candidate.ligand.ligand_id, [])
        item["pipeline_included"] = bool(candidate.rank)
        item["parent_id"] = candidate.ligand.parent_id
        item["parent_smiles"] = candidate.ligand.parent_smiles
        item["hop_metrics"] = candidate.ligand.hop_metrics
        item["parent_metadata_source"] = candidate.ligand.parent_metadata_source
        item["derivation"] = candidate.ligand.derivation
        item["default_visible"] = (hazard_severity_of(candidate.admet) != "critical"
                                   and classify(view)["role"] != "counter_control")
        item["docking_status"] = "EVALUATED" if candidate.pose else "FAILED" if candidate.dock_attempted else "NOT_EVALUATED"
        output.append(item)
    output.sort(key=lambda r: (-r["discovery_score"], r["ligand_id"]))
    for index, item in enumerate(output, 1):
        item["discovery_rank"] = index
    return output
