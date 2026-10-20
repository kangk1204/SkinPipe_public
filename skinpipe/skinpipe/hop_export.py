"""Export a saved-parent analogue run without deleting policy-rejected proposals.

This is a hypothesis-generation workflow, not an assay or a target validation.
The legacy hop generator's survivor list is used only for optional conformer
comparisons. Every distinct parseable non-parent proposal reaches the exploratory
table, and every raw emitted proposal remains in the generation census.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import shutil
import tempfile
from pathlib import Path

from rdkit import Chem

from . import hop_generate, hop_rank, report
from .filters import admet
from .library import _descriptors, standardize
from .pipeline import ScreenResult, assess_skin
from .publication import publish_generation
from .saved_parent import select_parent, source_sha
from .triage import TriageResult
from .types import Candidate, Grade, LigandRecord

POLICY = "saved-parent-hop-export-v1"


def _sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _json(data) -> bytes:
    return (json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False) + "\n").encode()


def _finite(data):
    """Convert optional numeric table values; missing values stay missing."""
    if isinstance(data, dict):
        return {k: _finite(v) for k, v in data.items()}
    if isinstance(data, (list, tuple)):
        return [_finite(v) for v in data]
    if hasattr(data, "item"):
        return _finite(data.item())
    if isinstance(data, float) and not math.isfinite(data):
        return None
    return data


def _budgets(n, max_3d, conformers, max_sa):
    for value, low, high, name in ((n, 1, 1000, "brics_products"),
                                  (max_3d, 0, 200, "max_3d"),
                                  (conformers, 1, 30, "conformers")):
        if isinstance(value, bool) or not isinstance(value, int) or not low <= value <= high:
            raise ValueError(f"{name} must be an integer between {low} and {high}")
    if (isinstance(max_sa, bool) or not isinstance(max_sa, (int, float))
            or not math.isfinite(max_sa) or not 1 <= max_sa <= 10):
        raise ValueError("max_sa must be a finite number between 1 and 10")


def export_hops(source: Path, parent_id: str, destination: Path, *,
                brics_products: int = 200, max_3d: int = 50,
                conformers: int = 10, max_sa: float = 6.0) -> dict:
    """Create a new result directory bound to an explicit saved parent.

    ``brics_products`` bounds BRICS output; the finite bioisostere rule table
    has its own count. ``max_sa`` and the legacy sanity policy govern costly
    3D comparisons only. They do not remove proposals from the exported pool.
    The output is published atomically by the existing generation store. An
    existing destination, including a concurrent publisher's result, is refused.
    """
    _budgets(brics_products, max_3d, conformers, max_sa)
    source, destination = Path(source).resolve(strict=True), Path(destination)
    if os.path.lexists(destination):
        raise ValueError("output already exists; choose a new hop run directory")
    selection = select_parent(source, parent_id)
    parent = selection["parent"]
    smi = parent.get("smiles")
    mol = Chem.MolFromSmiles(smi) if isinstance(smi, str) and smi.strip() else None
    if (mol is None or not 1 <= mol.GetNumHeavyAtoms() <= 250
            or len(Chem.GetMolFrags(mol)) != 1
            or any(a.GetAtomicNum() == 0 or a.GetNumRadicalElectrons() for a in mol.GetAtoms())):
        raise ValueError("parent must be one parseable connected molecule without dummies or radicals")
    parent_smiles = Chem.MolToSmiles(mol)
    if standardize(parent_smiles) != parent_smiles:
        raise ValueError("parent would change during generator neutralisation; use an explicitly standardised saved parent")
    marker, node = selection["marker"], selection["node"]
    source_digest = selection["source_sha256"]
    parent_sha = _sha(_json(parent))
    parameters = {"parent_id": parent_id, "brics_products": brics_products,
                  "max_3d": max_3d, "conformers": conformers, "max_sa": max_sa,
                  "shape_prescreen": False, "min_topical_score": 0,
                  "seed_activity_or_patent_filters": False,
                  "admet_applied_percent": 1.0, "admet_product": "body_lotion"}

    def build(root: Path, previous: Path | None) -> dict:
        if previous is not None:
            raise ValueError("output appeared during generation; existing output preserved")
        # Report writers correctly refuse generation-store paths. Compute in
        # a separate mutable staging folder, then copy a complete snapshot
        # into the immutable publication. Failed staging is kept for diagnosis.
        work = Path(tempfile.mkdtemp(prefix=".hop-export-staging-", dir=root.parent.parent))
        kept, census = hop_generate.enumerate_hops(
            parent_smiles, n=brics_products, max_sa=max_sa, return_census=True)
        if census["generation_parent_smiles"] != parent_smiles:
            raise ValueError("generator changed the saved parent structure")
        proposals = census["proposals"]
        inventory = {}
        for proposal in proposals:
            canonical = proposal["canonical_smiles"]
            if canonical is not None and canonical != parent_smiles:
                inventory.setdefault(canonical, []).append(proposal)

        # rank_hops standardises its inputs. Do not attach a neutralised
        # molecule's shape score to the exact charged proposal shown on the web.
        shape_inputs = [s for s in kept if standardize(s) == s]
        shape_rows, shape_census = {}, {"status": "NOT_EVALUATED", "reason": "max_3d is zero"}
        if max_3d and shape_inputs:
            frame = hop_rank.rank_hops(parent_smiles, shape_inputs, n_confs=conformers,
                                      max_3d=max_3d, max_sa=max_sa, prescreen=False,
                                      min_topical_score=0)
            shape_rows = {r["smiles"]: _finite(r) for r in frame.to_dict("records")}
            shape_census = {"status": "COMPUTED", **_finite(dict(frame.attrs))}
        elif max_3d:
            shape_census = {"status": "NOT_EVALUATED", "reason": "no exact-identity legacy survivors"}

        candidates = []
        for canonical, routes in sorted(inventory.items()):
            ligand_id = "HOP_" + _sha(_json([parent_id, parent_smiles, canonical]))[:32]
            accepted = canonical in kept
            reasons = sorted({p["sanity_policy_reason"] for p in routes if p["sanity_policy_reason"]})
            if not accepted and not reasons:
                reasons = ["sa_reject"]
            shape = shape_rows.get(canonical)
            if shape is None:
                shape_status = "NOT_EVALUATED"
                shape_reason = ("legacy comparison policy: " + ", ".join(reasons) if not accepted
                                else "comparison would alter structure" if canonical not in shape_inputs
                                else "3D comparison disabled" if not max_3d else "outside 3D comparison budget")
            else:
                shape_status = "EVALUATED" if shape.get("shape_tanimoto") is not None else "NOT_EVALUATED"
                shape_reason = None if shape_status == "EVALUATED" else "conformer comparison unavailable"
            metrics = {"comparison_status": shape_status, "comparison_reason": shape_reason,
                       "comparison": shape, "legacy_generation_accepted": accepted,
                       "legacy_policy_reasons": reasons,
                       "source_result_sha256": source_digest, "parent_row_sha256": parent_sha,
                       "proposal_ordinals": [p["ordinal"] for p in routes],
                       "derivations": sorted({p["derivation"] for p in routes}),
                       "interpretation": "computed structural hypothesis; no measured activity inherited"}
            if shape:
                metrics.update({k: shape.get(k) for k in ("hop_quality", "tanimoto_2d",
                                "shape_tanimoto", "pharmacophore_similarity",
                                "scaffold_comparable", "scaffold_changed",
                                "generic_scaffold_changed", "parent_murcko", "candidate_murcko")})
            ligand = LigandRecord(ligand_id=ligand_id, smiles=canonical, source="scaffold_hop",
                                  parent_id=parent_id, parent_smiles=parent_smiles,
                                  derivation=metrics["derivations"][0], hop_metrics=metrics,
                                  parent_metadata_source=f"saved result SHA256 {source_digest}; source_parent.json",
                                  flags=["generated: hypothesis only, no assay data"])
            candidate = Candidate(ligand=ligand, target=node, dock_attempted=False)
            # Policy-limited structures remain in the pool with missing model
            # outputs. No disconnected fragment is silently discarded here.
            if not reasons or reasons == ["sa_reject"]:
                ligand.props.update(_descriptors(Chem.MolFromSmiles(canonical)))
                candidate.skin = assess_skin(ligand, node)
                candidate.admet = admet.summarise(canonical, ligand_id=ligand_id,
                                                target_gene=node, applied_percent=1.0,
                                                product="body_lotion").to_dict()
                metrics["chemical_assessment_status"] = "EVALUATED"
            else:
                metrics["chemical_assessment_status"] = "NOT_EVALUATED"
                metrics["chemical_assessment_reason"] = "outside legacy model policy: " + ", ".join(reasons)
            candidates.append(candidate)

        result = ScreenResult(
            marker=marker,
            triage=TriageResult(marker=marker, grade=Grade.D, provisional=True,
                                screening_node=node, n_considered=len(candidates), structure_ready=False,
                                recommendation="Generated analogues only; target feasibility and functional activity NOT_EVALUATED.",
                                counts_unresolved=True),
            all_candidates=candidates,
            stages_run=["saved_parent_binding", "analogue_generation", "exploratory_chemical_assessment"],
            stages_skipped={"docking": "not requested", "functional_validation": "no generated-compound assay data"},
            gates_unverified={"generated_activity": "Activity is not inherited from the parent.",
                              "target_grade": "No target grade was evaluated; D is a provisional schema placeholder.",
                              "sourcing_and_regulation": "No availability, patent or regulatory clearance established."},
            filter_profile={"generation_policy": "annotated, never deletes exported proposals",
                            "default_web_filters": "reversible critical-hazard and counter-control filters"},
            provenance_audit={"policy": POLICY, "source_result_sha256": source_digest,
                              "source_selection": {k: v for k, v in selection.items() if k not in {"parent", "marker", "node"}},
                              "parent_row_sha256": parent_sha, "parent_id": parent_id,
                              "warning": "Computational hypotheses; parent measurements are only in source_parent.json."})
        report.to_json(result, work / "result.json")
        payload = json.loads((work / "result.json").read_bytes())
        payload["triage"]["grade_status"] = "NOT_EVALUATED"
        payload["hop_generation"] = {"policy": POLICY, "parameters": parameters,
                                     "raw_proposals": len(proposals), "exported_candidates": len(candidates),
                                     "legacy_survivors": len(kept), "shape_census": shape_census,
                                     "generation_file": "generation.json", "parent_file": "source_parent.json"}
        (work / "result.json").write_bytes(_json(payload))
        (work / "generation.json").write_bytes(_json({"policy": POLICY, "census": census,
                                                      "shape_census": shape_census}))
        (work / "source_parent.json").write_bytes(_json(parent))
        report.to_manifest(result, work / "run_manifest.json", parameters=parameters)
        summary = {"policy": POLICY, "source_path": str(source), "source_result_sha256": source_digest,
                   "source_bytes": selection["source_bytes"], "source_pool": selection["pool"],
                   "source_pool_rows": selection["pool_rows"],
                   "parent_row_sha256": parent_sha, "parent_id": parent_id,
                   "parent_smiles": parent_smiles, "parameters": parameters,
                   "raw_proposals": len(proposals), "exported_candidates": len(candidates),
                   "unparseable_proposals": sum(p["canonical_smiles"] is None for p in proposals),
                   "regenerated_parent_proposals": sum(p["regenerated_parent"] for p in proposals),
                   "legacy_survivors": len(kept), "shape_rows": len(shape_rows),
                   "result_sha256": _sha((work / "result.json").read_bytes())}
        if source_sha(source) != source_digest:
            raise ValueError("saved source changed during generation; output not published")
        (work / "hop_run.json").write_bytes(_json(summary))
        shutil.copytree(work, root, dirs_exist_ok=True)
        if source_sha(source) != source_digest:
            raise ValueError("saved source changed before publication; output not published")
        shutil.rmtree(work)
        return summary

    publication = publish_generation(destination, build, require_new=True)
    return {"output": str(destination), "publication": publication,
            **json.loads((destination / "hop_run.json").read_bytes())}
