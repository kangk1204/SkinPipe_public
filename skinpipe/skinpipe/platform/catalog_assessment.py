"""Cheap structure rules reused from SkinPipe; none is an experimental verdict."""

from dataclasses import asdict
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors
from ..types import LigandRecord, Candidate, ADMETProfile
from ..filters import topical, safety, stability, admet
from ..discovery import priorities


def freeze_skin_context(body=None):
    """Freeze the use scenario; a delivery goal is not a biological claim."""
    from ..filters.admet_exposure import PRODUCT_USE, validate_applied_percent

    supplied = {} if body is None else body
    allowed = {"product", "applied_percent", "goal", "research_goal", "application_site"}
    if not isinstance(supplied, dict) or set(supplied) - allowed:
        raise ValueError("Skin context accepts product, applied_percent and goal")
    product = supplied.get("product")
    if product is not None and (not isinstance(product, str) or product not in PRODUCT_USE):
        raise ValueError("Unknown skin product scenario")
    pct = supplied.get("applied_percent")
    if pct is not None and (isinstance(pct, bool) or not isinstance(pct, (int, float))):
        raise ValueError("applied_percent must be a numeric concentration")
    pct = validate_applied_percent(pct) if pct is not None else None
    goal = supplied.get("goal", "exploratory")
    if not isinstance(goal, str) or goal not in {"surface", "epidermal", "dermal", "exploratory"}:
        raise ValueError("Unknown skin delivery goal")
    research_goal, site = supplied.get("research_goal"), supplied.get("application_site")
    for key, value, maximum in (("research_goal", research_goal, 512), ("application_site", site, 128)):
        if value is not None and (not isinstance(value, str) or len(value) > maximum
                                  or any(ord(c) < 32 for c in value)):
            raise ValueError(key + " must be bounded plain text or null")
    return {
        "product": product,
        "applied_percent": pct,
        "goal": goal,
        "research_goal": research_goal,
        "application_site": site,
        "defaults_applied": sorted(allowed - set(supplied)),
        "policy": "skin-use-scenario-v1",
    }


def structural_profile(smiles, ligand_id, context):
    """Do not manufacture an exposure scenario when formulation is unknown."""
    if context["product"] is not None and context["applied_percent"] is not None:
        return admet.summarise(
            smiles, ligand_id=ligand_id, product=context["product"],
            applied_percent=context["applied_percent"],
        ).to_dict()
    from ..filters import admet_report
    blocks = {
        "phototoxicity": admet.phototoxicity(smiles),
        "skin_irritation": admet.skin_irritation(smiles),
        "genotoxicity": admet.genotoxicity(smiles),
        "cmr": admet.cmr_alerts(smiles),
        "cyp_liability": admet.cyp_liability(smiles, sed_mg_kg_day=None),
    }
    concerns = []
    for domain, key, omit in (("phototoxicity", "phototoxicity", "uv_vis_chromophore"),
                              ("irritation", "skin_irritation", "ph_shift"),
                              ("genotoxicity", "genotoxicity", None),
                              ("cmr", "cmr", None), ("cyp", "cyp_liability", None)):
        admet_report.collect(domain, blocks[key], concerns, omit)
    admet_report.rank(concerns)
    blocks["systemic_exposure"] = {
        "status": "not_assessed", "concern": "unknown",
        "reason": "product_and_applied_percent_required", "sed_mg_kg_day": None,
        "product": context["product"], "applied_percent": context["applied_percent"],
        "verdict": "Product format or applied concentration is unknown; no systemic exposure scenario was computed",
    }
    return ADMETProfile(
        smiles=smiles, ligand_id=ligand_id, **blocks, concerns=concerns,
        verdicts={key: block.get("verdict", "") for key, block in blocks.items()},
        overall="Structure-based prioritisation only; exposure, skin benefit and safety are not established",
    ).to_dict()


def assess_catalog(rows, *, skin_context=None, stop=lambda: None, progress=lambda _: None):
    # Historical callers retain their explicitly established legacy scenario.
    context = skin_context if skin_context is not None else {
        "product": "body_lotion", "applied_percent": 1.0, "goal": "exploratory",
    }
    candidates = []
    by_id = {}
    for index, row in enumerate(rows):
        stop()
        if row.get("discovery_terms") and skin_context is None:
            continue  # Preserve observations/terms of previously evaluated sources.
        mol = Chem.MolFromSmiles(row.get("smiles") or "")
        if mol is None or not mol.GetNumHeavyAtoms():
            continue
        props = {
            "mw": Descriptors.MolWt(mol),
            "logp": Crippen.MolLogP(mol),
            "tpsa": rdMolDescriptors.CalcTPSA(mol),
            "hbd": rdMolDescriptors.CalcNumHBD(mol),
            "hba": rdMolDescriptors.CalcNumHBA(mol),
        }
        ligand = LigandRecord(
            row["candidate_id"],
            row["smiles"],
            props=props,
            natural_product_sources=row.get("natural_product_sources") or [],
        )
        skin = topical.assess(ligand)
        sens = safety.summarise(ligand.smiles)
        stab = stability.summarise(ligand.smiles)
        skin.sensitisation_status = sens["assessment_status"]
        skin.stability_status = stab["assessment_status"]
        skin.sensitiser_alerts = [
            f"{x['severity']}:{x['name']}" for x in sens["alerts"]
        ]
        skin.stability_alerts = [f"{x['mode']}:{x['name']}" for x in stab["alerts"]]
        skin.n_pains = len(sens["pains"])
        skin.colour = stab["colour"]
        profile = structural_profile(ligand.smiles, ligand.ligand_id, context)
        source_critical = (
            row.get("critical_warning") is True
            or row.get("hazard_severity") == "critical"
        )
        row.update(
            mw=props["mw"],
            logp=props["logp"],
            topical_score=skin.topical_score,
            skin=asdict(skin),
            structural_admet=profile,
            hazard_severity="critical"
            if source_critical
            else profile["hazard_severity"],
            critical_warning=profile["hazard_severity"] == "critical"
            or source_critical,
            descriptor_provenance={
                "method": "RDKit exact input structure",
                "runtime": Chem.rdBase.rdkitVersion,
            },
            structural_admet_context={
                "method": "SkinPipe structural heuristics",
                **context,
                "experimental_validation": "NOT_ESTABLISHED",
            },
        )
        row["skin_context"] = dict(context)
        row["skin_goal_assessment"] = {
            "goal": context["goal"],
            "permeation_proxy": skin.permeation_class,
            "depth_selectivity": "NOT_ESTABLISHED",
            "interpretation": (
                "Surface retention and barrier effects require formulation and functional assays; higher permeability is not rewarded"
                if context["goal"] == "surface" else
                "Delivery depth is not inferred from the permeability proxy; compare formulation scenarios and verify experimentally"
            ),
        }
        candidates.append(Candidate(ligand, row["target_id"], skin=skin, admet=profile))
        by_id[row["candidate_id"]] = row
        if index % 100 == 0:
            progress(
                {
                    "stage": "cheap_structure_assessment",
                    "processed": index + 1,
                    "total": len(rows),
                }
            )
    scores = priorities(candidates) if candidates else {}
    for cid, score in scores.items():
        by_id[cid].update(score)
        if skin_context is not None and context["goal"] in {"surface", "exploratory"}:
            # A surface/barrier ingredient need not penetrate skin. Keep the
            # descriptor visible, but do not reward permeability in triage.
            by_id[cid]["discovery_terms"]["topical"] = 50
            by_id[cid]["discovery_missing"] = sorted(set(
                by_id[cid].get("discovery_missing", []) + ["topical"]
            ))
            by_id[cid]["skin_goal_assessment"]["permeability_ranking"] = "neutral"
            from ..discovery import reweight_records
            by_id[cid].update(reweight_records([by_id[cid]])[0])
    return rows
