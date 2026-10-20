'Which discovery route a target deserves.'
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from enum import Enum
from typing import Optional

from .triage import FEASIBLE_PCHEMBL, FEASIBLE_TOPICAL, TriageResult, may_dock
from .types import Grade, TargetSpec, TargetMode

log = logging.getLogger(__name__)


class Track(str, Enum):
    HOP = "scaffold_hop"
    DE_NOVO = "de_novo_docking"
    LIGAND_BASED = "ligand_based"
    NOT_VIABLE = "not_viable"










HOP_LIGAND_FLOOR = 50
#: Below this, there is nothing to hop from and the library has to supply the idea.
DE_NOVO_LIGAND_CEILING = 15
#: Above this share of ligands appearing in patents, treat the space as crowded.
CROWDED_PATENT_FRACTION = 0.7

#: A seed whose Murcko scaffold is this large cannot be hopped by generating
#: analogues from its own fragments. Measured, not assumed: swapping the amide
#: of imatinib for an oxadiazole leaves the scaffold 69% identical by ECFP4,
#: because one ring of a five-ring core is most of nothing. The same swap on
#: kojic acid takes its scaffold similarity to 0.37, because its core IS that
#: one ring. Fragment recombination does reach 2D novelty on large seeds, and
#: loses the binding mode doing it -- the BRICS analogues of imatinib scored
#: 0.31 on 2D and 0.46 on shape, which is a different molecule, not a hop.
HOPPABLE_SCAFFOLD_RINGS = 3
HOPPABLE_SCAFFOLD_HEAVY = 20


@dataclass
class StrategyDecision:
    marker: str
    node: Optional[str]
    track: Track
    rationale: str
    next_steps: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    ip_outlook: str = ""
    n_ligands: int = 0
    n_structures: int = 0
    patent_fraction: Optional[float] = None
    purchasable_fraction: Optional[float] = None

    @property
    def actionable(self) -> bool:
        return self.track is not Track.NOT_VIABLE


def choose_track(spec: TargetSpec, triage: TriageResult,
                 patent_fraction: Optional[float] = None,
                 purchasable_fraction: Optional[float] = None,
                 notified_ingredient: Optional[str] = None) -> StrategyDecision:
    'Pick the discovery route and say why.'
    node = triage.screening_node
    picked = next((n for n in spec.nodes if n.gene == node), None)
    n_lig = picked.n_ligands_p6 if picked else 0
    n_pdb = picked.n_pdb if picked else 0

    d = StrategyDecision(marker=spec.marker, node=node, track=Track.NOT_VIABLE,
                         rationale="", n_ligands=n_lig, n_structures=n_pdb,
                         patent_fraction=patent_fraction,
                         purchasable_fraction=purchasable_fraction)

    if spec.mode is TargetMode.ENDPOINT:
        d.rationale = (f"{spec.marker} is a measurement, not a target. "
                       "It belongs in the efficacy protocol, not the screen.")
        d.next_steps = [f"use {spec.marker} as a clinical endpoint"]
        return d

    if node is None or triage.grade is Grade.D:
        # Four different things grade D, and they need four different next
        # moves. A missing lookup in particular must not read as a verdict on
        # the biology -- "go make biologics" is expensive advice to give
        # because ChEMBL returned a 500 -- and a target whose ligands all fail
        # the skin filter is a delivery problem wearing a scarcity grade.
        if node is None:
            d.rationale = (f"no druggable node is registered for {spec.marker}, so "
                           "there is nothing to screen against.")
            d.next_steps = ["register a node for this marker, or find one in the "
                            "same pathway"]
        elif triage.counts_unresolved:
            d.rationale = (f"{spec.marker} cannot be graded: the counts behind {node} "
                           "were never retrieved. This is missing data, not a finding "
                           "that the target is barren.")
            d.next_steps = ["re-run the triage once UniProt/ChEMBL/RCSB answer",
                            "do not drop this marker on the strength of this row"]
        elif not triage.provisional and triage.n_considered < n_lig:


























            n_deliverable = triage.n_considered
            d.rationale = (f"of {node}'s {n_lig} potent ligand(s), {n_deliverable} "
                           f"clear both pChEMBL >= {FEASIBLE_PCHEMBL:g} and topical "
                           f">= {FEASIBLE_TOPICAL:g}, which is too few to start a "
                           "topical campaign from. That bar is higher than the "
                           "screening gates, so compounds this rejects may still "
                           "appear in the ranked list below. The chemical matter "
                           "exists and the skin barrier rejects it; this is a delivery "
                           "finding, not a barren target.")
            d.next_steps = [
                f"re-read the {node} series for analogues that trade potency for "
                "permeability -- MW, logP and HBD are what the topical filter rejects on",
                "or hop to a smaller, less polar core and re-check delivery before "
                "re-checking potency",
                "consider peptides or biologics if the chemotype cannot be made to cross",
                "or find a different node in the same pathway",
            ]
        else:
            d.rationale = (f"{node} carries {n_lig} potent ligand(s), too few to start a "
                           "small-molecule campaign from. A screen here would be "
                           "starting from nothing.")
            d.next_steps = ["consider peptides or biologics",
                            "or find a different node in the same pathway"]
        # The early exit used to drop these, so the loudest caveat in the run --
        # "grade D is NOT a finding" -- never reached the decision object.
        d.warnings.extend(triage.warnings)
        return d

    # A PDB entry mapped to the accession is not automatically a docking target,
    # so the route is chosen on the same verdict the docking stage will act on.
    # Two implementations of "is there a structure" is how this file and triage
    # came to give opposite instructions for tyrosinase.
    verdict = may_dock(spec, node=node)
    has_structure = verdict.allowed
    if not verdict.allowed and verdict.n_pdb:
        # Only when entries exist: the decision reports n_structures=n_pdb, and
        # a reader who sees a number next to "no structure" needs the reason.
        d.warnings.append(verdict.reason)

    if n_lig >= HOP_LIGAND_FLOOR and has_structure:
        d.track = Track.HOP
        d.rationale = (f"{node} has {n_lig} potent activity records (ChEMBL "
                       f"activity rows, not distinct molecules) and {n_pdb} "
                       "experimental structures. The binders already exist, so the "
                       "scarce thing is not a hit but a chemotype you can own.")
        d.next_steps = [
            f"pick 3-5 potent, topically plausible seeds from the {node} screen",
            "check each seed with seed_is_hoppable(); a large multi-ring scaffold needs "
            "3D shape screening of a purchasable library, not fragment generation",
            "for the small seeds, generate hops and keep low 2D with high 3D similarity",
            "dock the hops into the same receptor and rescore with MM/GBSA",
            "check each survivor for vendor availability and patent presence",
        ]
    elif n_lig <= DE_NOVO_LIGAND_CEILING and has_structure:
        d.track = Track.DE_NOVO
        d.rationale = (f"{node} has only {n_lig} potent activity records but "
                       f"{n_pdb} structures. There is nothing to hop from, so the "
                       "library has to supply the starting point.")
        d.next_steps = [
            f"prepare the best {node} structure and validate the box by redocking",
            "dock a purchasable library, natural-product subsets first",
            "rescore the top poses with MM/GBSA",
            "buy and assay the survivors before any chemistry",
        ]
    elif has_structure:
        d.track = Track.DE_NOVO
        d.rationale = (f"{node} sits between the two routes with {n_lig} potent "
                       "activity records. Screening a purchasable library is the "
                       "cheaper first move; hop later if a series emerges.")
        d.next_steps = [
            "dock a purchasable library and see whether a series appears",
            f"revisit hopping once {node} has a confirmed in-house hit",
        ]
    else:
        d.track = Track.LIGAND_BASED
        d.rationale = (f"{node} has {n_lig} potent activity records and no usable "
                       "experimental structure, so docking is not available. "
                       "Similarity and pharmacophore work still are.")
        d.next_steps = [
            f"build a pharmacophore from the known {node} actives",
            "search purchasable libraries by shape and pharmacophore, not docking",
            (f"if a structure-based pass is wanted, build a homology model on "
             f"{verdict.homology_template} ({verdict.template_n_pdb} experimental "
             f"structures), declare it as a model, and confirm every hit on {node}"
             if verdict.homology_template else
             "consider a homology model only if a close template exists, and say so"),
        ]

    # ---- IP outlook -------------------------------------------------------
    if patent_fraction is None:
        d.ip_outlook = ("patent coverage of this chemotype is unmeasured. Run the "
                        "sourcing stage before assuming the space is open.")
    elif patent_fraction >= CROWDED_PATENT_FRACTION:
        d.ip_outlook = (f"{patent_fraction:.0%} of the compounds screened here appear in the "
                        "patent literature. The space is crowded, which strengthens the case "
                        "for hopping to a distinct chemotype rather than optimising inside it.")
        if d.track is Track.DE_NOVO:
            d.warnings.append(
                "the library screen may rediscover compounds that are already claimed; "
                "carry the patent flag through to the hit list")
    else:
        d.ip_outlook = (f"{patent_fraction:.0%} of the compounds screened here appear in "
                        "patents, so there is room. Presence or absence in SureChEMBL is a "
                        "signal, not a freedom-to-operate opinion.")

    if purchasable_fraction is not None and purchasable_fraction < 0.3:
        d.warnings.append(
            f"only {purchasable_fraction:.0%} of this target's ligands are sold by a "
            "vendor. Expect to commission synthesis, and budget for it before "
            "committing to the target.")

    # ---- Korean functional-cosmetic context -------------------------------
    if notified_ingredient:
        d.warnings.append(
            f"{notified_ingredient} is already a 고시원료 for this claim, so a formulation "
            "using it can go through 보고 rather than 심사. A novel ingredient forfeits "
            "that and carries the full safety dossier. Make that trade deliberately.")

    d.warnings.extend(triage.warnings)
    return d


def seed_is_hoppable(smiles: str) -> tuple[bool, str]:
    """Whether generating analogues from this seed can produce a real hop.

    Returns ``(hoppable, reason)``. A large multi-ring scaffold is not hoppable
    by generation: the only way to change its core is to replace most of it,
    and the methods that do that discard the shape that made it bind. Such a
    target is screened by 3D shape against a purchasable library instead, using
    a known binder's conformer as the query.
    """
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
    from rdkit.Chem.Scaffolds import MurckoScaffold

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "구조를 해석할 수 없습니다"
    core = MurckoScaffold.GetScaffoldForMol(mol)
    rings = rdMolDescriptors.CalcNumRings(core)
    heavy = core.GetNumHeavyAtoms()
    if rings == 0:
        return False, ("고리가 없는 골격이라 Murcko 스캐폴드가 비어 있습니다. "
                       "호핑 대상이 아니라 유사체 최적화 대상입니다.")
    if rings >= HOPPABLE_SCAFFOLD_RINGS or heavy >= HOPPABLE_SCAFFOLD_HEAVY:
        return False, (f"스캐폴드가 고리 {rings}개, 중원자 {heavy}개로 큽니다. 작용기 하나를 "
                       "바꿔도 코어의 대부분이 남아 호핑이 되지 않고, 코어를 통째로 바꾸면 "
                       "결합 형태를 잃습니다. 알려진 결합자의 3차원 형태를 쿼리로 구매 "
                       "가능한 라이브러리를 형태 기반 검색하세요.")
    return True, (f"스캐폴드가 고리 {rings}개, 중원자 {heavy}개로 작아 생물학적 등가체 "
                  "치환 하나가 코어 전체를 바꿉니다. 생성형 호핑이 유효합니다.")


def summarise(decision: StrategyDecision) -> str:
    """One readable block for a report or a terminal."""
    lines = [f"{decision.marker} -> {decision.track.value}", "", decision.rationale, ""]
    if decision.ip_outlook:
        lines += [f"IP: {decision.ip_outlook}", ""]
    if decision.next_steps:
        lines += ["Next:"] + [f"  {i}. {s}" for i, s in enumerate(decision.next_steps, 1)] + [""]
    if decision.warnings:
        lines += ["Warnings:"] + [f"  ! {w}" for w in decision.warnings]
    return "\n".join(lines)
