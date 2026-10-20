"""Shared data contracts for the skin small-molecule discovery pipeline.

Every stage consumes and produces these types. Keep them dependency-light:
only stdlib + typing, so any module can import this without pulling in
RDKit / OpenMM / Vina.
"""
from __future__ import annotations

from dataclasses import dataclass, field, asdict
from enum import Enum
import math
from pathlib import Path
from typing import Any, ClassVar, Optional


class Direction(str, Enum):
    """Desired movement of the biomarker, taken from the user's marker table."""
    INHIBIT = "inhibit"
    ACTIVATE = "activate"
    INCREASE = "increase"
    DECREASE = "decrease"
    SUPPRESS = "suppress"
    UNKNOWN = "unknown"


class TargetMode(str, Enum):
    """How a small molecule can actually reach the biomarker."""
    DIRECT = "direct"        # the marker itself has a bindable pocket
    PROXY = "proxy"          # bind an upstream node instead
    ENDPOINT = "endpoint"    # not a molecular target at all (TEWL, 8-OHdG)


class Grade(str, Enum):
    A = "A"   # ready to screen
    B = "B"   # workable, needs optimisation
    C = "C"   # thin chemical matter
    D = "D"   # not a small-molecule proposition


class Evidence(str, Enum):
    """Evidence ladder for a compound-target claim, strongest first."""
    CLINICAL = "clinical"
    EX_VIVO_HUMAN = "ex_vivo_human"
    HUMAN_CELL = "human_cell"
    NONHUMAN_CELL = "nonhuman_cell"
    BIOCHEMICAL = "biochemical"







    SURROGATE = "surrogate"
    ANIMAL = "animal"
    PREDICTED = "predicted"
    UNKNOWN = "unknown"


@dataclass
class DruggableNode:
    """A protein a small molecule can actually bind."""
    gene: str
    uniprot: Optional[str] = None
    chembl_target: Optional[str] = None
    n_pdb: int = 0
    n_ligands_p6: int = 0
    rationale: str = ""          # why this node serves the marker







    node_action: Optional[str] = None
    goal_support: dict[str, Any] = field(default_factory=dict)






    potency_floor: Optional[float] = None
    potency_floor_why: Optional[str] = None
    primary: bool = False        # this node IS the marker, not a proxy or homologue
    preferred: bool = False      # the configured route for this marker, even if
                                 # another node carries more chemical matter
    structure_usable: bool = True         # False when the PDB entries mapped to this
                                 # accession are not the protein you would dock into
    structure_note: str = ""              # why, when structure_usable is False
    metalloenzyme: Optional[str] = None   # catalytic metal, e.g. "Cu" or "Zn".
                                 # Vina has no metal-coordination term, so its
                                 # scores on chelating ligands are unreliable.
    # `n_pdb` and `n_ligands_p6` default to 0, and 0 is also what a failed
    # UniProt/ChEMBL/RCSB call leaves behind. Those two states grade to opposite
    # conclusions -- "no chemical matter, drop this target" versus "we do not
    # know yet" -- so the difference is recorded rather than inferred.
    resolved: bool = True
    resolve_error: str = ""


@dataclass
class TargetSpec:
    """One row of the biomarker table, resolved to something actionable."""
    marker: str
    category: str
    direction: Direction
    mode: TargetMode
    nodes: list[DruggableNode] = field(default_factory=list)
    grade: Optional[Grade] = None
    n_feasible: int = 0          # potent AND topically plausible ligands
    notes: list[str] = field(default_factory=list)

    @property
    def is_screenable(self) -> bool:
        return self.mode is not TargetMode.ENDPOINT and bool(self.nodes)


@dataclass
class LigandRecord:
    """A candidate compound entering the pipeline."""
    ligand_id: str
    smiles: str
    name: Optional[str] = None
    source: str = "unknown"          # chembl | npass | coconut | user_sdf | inci

















    natural_product: Optional[bool] = None
    inci_name: Optional[str] = None
    max_phase: Optional[float] = None
    props: dict[str, Any] = field(default_factory=dict)      # MW, logP, TPSA...
    activity: dict[str, Any] = field(default_factory=dict)   # pchembl, assay info
    flags: list[str] = field(default_factory=list)           # safety / stability alerts
    # Provenance for computationally generated structures (see hopping.py).
    # A generated analogue has no assay data of its own and must never be
    # confused with a fetched compound that does; these two fields are how a
    # reader tells them apart downstream.
    parent_id: Optional[str] = None      # ligand_id of the compound it was derived from
    derivation: Optional[str] = None     # how, e.g. "bioisostere:acid_to_tetrazole"
    natural_product_sources: list[dict[str, Any]] = field(default_factory=list)
    original_structures: list[dict[str, Any]] = field(default_factory=list)
    parent_smiles: Optional[str] = None
    hop_metrics: dict[str, Any] = field(default_factory=dict)
    parent_metadata_source: Optional[str] = None

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d.update({f"prop_{k}": v for k, v in self.props.items()})
        d.pop("props")
        return d


@dataclass
class Receptor:
    """A docking-ready protein structure with a defined search box."""
    gene: str
    pdb_id: str
    path: Path                       # cleaned + protonated PDB
    pdbqt: Optional[Path] = None     # Vina-ready receptor
    center: tuple[float, float, float] = (0.0, 0.0, 0.0)
    size: tuple[float, float, float] = (22.0, 22.0, 22.0)
    chain: Optional[str] = None
    cocrystal_ligand: Optional[str] = None
    # Crystallographic pose of `cocrystal_ligand`, bond orders restored from the
    # RCSB chemical component. Written by structure.prepare, consumed by
    # dock.vina.redock_validation as the RMSD reference.
    cocrystal_sdf: Optional[Path] = None





    cocrystal_completeness: Optional[float] = None
    cocrystal_is_fragment: bool = False
    cocrystal_atoms: Optional[tuple[int, int]] = None   # (matched, component)
    resolution: Optional[float] = None
    method: Optional[str] = None
    pocket_source: str = "unknown"   # cocrystal | uniprot_binding_site | centroid | manual
    notes: list[str] = field(default_factory=list)
    # Source-bound deposited assembly membership and transforms. Empty for
    # historical/manual receptors whose assembly context was not evaluated.
    assembly: dict[str, Any] = field(default_factory=dict)
    preparation_evidence: dict[str, Any] = field(default_factory=dict)


@dataclass
class DockPose:
    ligand_id: str
    pose_index: int
    affinity: float                  # Vina score, kcal/mol (more negative = better)
    sdf_path: Optional[Path] = None
    rmsd_lb: Optional[float] = None
    rmsd_ub: Optional[float] = None
    ligand_efficiency: Optional[float] = None   # -affinity / heavy atoms


@dataclass
class GBSAResult:
    ligand_id: str
    pose_index: int
    dG_bind: Optional[float] = None  # kcal/mol, MM/GBSA interaction + solvation
    e_complex: Optional[float] = None
    e_receptor: Optional[float] = None
    e_ligand: Optional[float] = None
    e_gas: Optional[float] = None
    e_solv: Optional[float] = None       # polar + nonpolar; e_gas + e_solv == dG_bind
    e_polar: Optional[float] = None      # GB polar solvation alone
    e_nonpolar: Optional[float] = None   # ACE surface-area term alone
    dG_std: Optional[float] = None       # spread across MD snapshots, when averaged
    n_frames: int = 1
    wall_s: Optional[float] = None
    converged: bool = False
    notes: str = ""
    execution: Optional[dict[str, Any]] = None
    parameterization: Optional[dict[str, Any]] = None

    @property
    def resolved(self) -> bool:
        """Legacy indicator of reported within-run SD below 2 kcal/mol.

        Missing SD and single-frame calculations do not establish a sampling
        spread. This flag does not establish repeated-run stability, separation
        from another ligand, validated ranking, or experimental affinity.
        """
        return self.sampling_spread_below_legacy_threshold

    @property
    def sampling_spread_below_legacy_threshold(self) -> bool:
        """Apply the historical 2 kcal/mol SD heuristic to usable sampling."""
        return (self.converged is True and isinstance(self.n_frames, int) and
                not isinstance(self.n_frames, bool) and self.n_frames >= 2 and
                isinstance(self.dG_bind, (int, float)) and not isinstance(self.dG_bind, bool) and
                math.isfinite(self.dG_bind) and isinstance(self.dG_std, (int, float)) and
                not isinstance(self.dG_std, bool) and math.isfinite(self.dG_std) and
                0 <= self.dG_std < 2.0)


@dataclass
class ADMETConcern:
    """One ranked toxicological finding, with mechanism and what to do next.

    Severity ladder, highest first:
        critical  regulatory stopper for a leave-on cosmetic (CMR, corrosive)
        high      needs the reference in-vitro test before a use level is set
        moderate  prioritisation signal; confirm at the intended use level
        low       noted, no action expected
        context   the alert IS the intended mechanism, or is irrelevant here
    """
    domain: str          # phototoxicity | irritation | genotoxicity | systemic | cyp | cmr
    name: str
    severity: str
    mechanism: str = ""
    action: str = ""     # the assay or restriction this points at


@dataclass
class ADMETProfile:
    """Topical ADMET/tox profile for one structure.

    Deliberately NOT an oral-drug profile. The exposure being modelled is a
    leave-on product applied daily, for years, to sun-exposed skin, by a
    population that includes pregnant women and children. Phototoxicity and
    dermal-route systemic exposure therefore outrank CYP and hERG, which are
    only meaningful once the Systemic Exposure Dosage is non-trivial.
    """
    smiles: str
    ligand_id: Optional[str] = None
    phototoxicity: dict[str, Any] = field(default_factory=dict)
    skin_irritation: dict[str, Any] = field(default_factory=dict)
    genotoxicity: dict[str, Any] = field(default_factory=dict)
    systemic_exposure: dict[str, Any] = field(default_factory=dict)
    cyp_liability: dict[str, Any] = field(default_factory=dict)
    cmr: dict[str, Any] = field(default_factory=dict)
    concerns: list[ADMETConcern] = field(default_factory=list)   # ranked, worst first
    verdicts: dict[str, str] = field(default_factory=dict)       # one line per domain
    overall: str = ""                                            # plain-language call

    SEVERITY_ORDER: ClassVar[tuple[str, ...]] = (
        "critical", "high", "moderate", "low", "context")

    @property
    def assessed(self) -> bool:
        """False when no assessment was actually made.

        An unparseable structure returns a profile with no concerns and no
        verdicts, which is indistinguishable from a clean bill of health if you
        only look at the severities. Every real assessment fills one verdict per
        domain, so an empty `verdicts` is the honest signal that nothing ran.
        """
        return bool(self.verdicts)

    @property
    def worst_severity(self) -> str:
        """The highest severity present, or 'none' when nothing fired."""
        for level in self.SEVERITY_ORDER:
            if any(c.severity == level for c in self.concerns):
                return level
        return "none"

    @property
    def hazard_severity(self) -> str:
        """Worst severity among the hazard domains, ignoring systemic exposure.

        This -- not `worst_severity` -- is what gates a hit list, and the two
        disagree routinely. A large Systemic Exposure Dosage is a magnitude, not
        a liability: it says how much gets in, and whether that matters depends
        on a relevant NOAEL the screen does not have. Common cosmetic ingredients
        can exceed a modeled exposure threshold; that alone establishes neither
        hazard nor safety. Reporting `worst_severity`
        next to a gate that used this one is how a reader concludes the gate is
        broken, so both travel together.
        """
        for level in self.SEVERITY_ORDER:
            if any(c.severity == level and c.domain != "systemic"
                   for c in self.concerns):
                return level
        return "none"

    def concerns_in(self, domain: str) -> list[ADMETConcern]:
        return [c for c in self.concerns if c.domain == domain]

    def to_dict(self) -> dict[str, Any]:
        # These are properties, so asdict drops them. Every downstream stage
        # keys on them, so they are written in explicitly.
        return asdict(self) | {"worst_severity": self.worst_severity,
                               "hazard_severity": self.hazard_severity,
                               "assessed": self.assessed}


@dataclass
class SkinAssessment:
    """The skin-specific verdict on one compound. This is the layer that makes
    the pipeline dermatological rather than a generic docking wrapper."""
    ligand_id: str
    log_kp: Optional[float] = None        # Potts-Guy, cm/h
    permeation_class: str = "unknown"     # good | moderate | poor | blocked
    topical_score: float = 0.0            # 0-10
    sensitiser_alerts: list[str] = field(default_factory=list)
    sensitiser_details: list[dict[str, Any]] = field(default_factory=list)
    stability_alerts: list[str] = field(default_factory=list)
    n_pains: Optional[int] = 0            # None when the structure was not assessed
    regulatory_flags: list[str] = field(default_factory=list)
    verdict: str = ""                     # short human-readable call
    # None denotes a legacy assessment with unknown execution provenance.
    # Structure parseability never establishes that these filters actually ran.
    sensitisation_status: Optional[str] = None
    stability_status: Optional[str] = None
    colour: Optional[str] = None


@dataclass
class Candidate:
    """Final integrated row: chemistry + structure + skin fitness."""
    ligand: LigandRecord
    target: str
    pose: Optional[DockPose] = None
    gbsa: Optional[GBSAResult] = None
    skin: Optional[SkinAssessment] = None
    admet: Optional[dict[str, Any]] = None       # ADMETProfile.to_dict()
    eligibility: Optional[dict[str, Any]] = None  # cosmetic-eligibility verdict
    # Which way this compound pushes `target`, and whether that is the way the
    # node has to go: {"verdict", "action_type", "effect", "source", "required",
    # "explanation"} from skinpipe.fetch.mechanism. "unknown" is the common
    # case and is never "compatible"; see that module for why the two are kept
    # apart. None means the direction stage did not run at all.
    direction: Optional[dict[str, Any]] = None
    evidence: Evidence = Evidence.PREDICTED
    #: Was *this* compound offered to the docking engine? `--dock-top` docks a
    #: potency-ordered slice, and the screen-level "docking ran" flag says
    #: nothing about the rest: every candidate outside the budget was scored as
    #: "docking returned no pose for this compound ... it is not excused",
    #: which is a failure report about a run that never happened. None means
    #: nothing set it and the screen-level flag is all there is.
    dock_attempted: Optional[bool] = None
    composite: float = 0.0
    rank: Optional[int] = None
    docking_execution: Optional[dict[str, Any]] = None


SEVERITY_ORDER = ["critical", "high", "moderate", "low", "context", "none"]

#: The severities a gate can actually be set to. `"none"` is in
#: `SEVERITY_ORDER` because it is a legitimate *answer* -- "nothing fired" --
#: but it is never a legitimate *threshold*: a gate is "drop everything at
#: least this bad", and `"none"` sits at the bottom of the ladder, so reading
#: it as a threshold drops every compound in the screen including the clean
#: ones. As a gate value it means what the CLI help says it means: run no gate.
GATE_SEVERITIES = tuple(s for s in SEVERITY_ORDER if s != "none")


def normalise_severity_gate(gate, *, field: str = "admet_gate"):
    'The severity threshold to gate on, or None for "run no gate".'
    if gate is None:
        return None
    if not isinstance(gate, str):
        raise ValueError(
            f"{field} must be a severity name or None, not {type(gate).__name__} "
            f"({gate!r}). Expected one of {', '.join(GATE_SEVERITIES)}, "
            "'none', or None.")
    value = gate.strip().lower()
    if value == "none":
        return None
    if value not in GATE_SEVERITIES:
        raise ValueError(
            f"unknown {field} {gate!r}; expected one of "
            f"{', '.join(GATE_SEVERITIES)}, or 'none'/None to run no gate. An "
            "unrecognised gate cannot exclude a single compound, so it is "
            "refused rather than run as a screen with the gate quietly off.")
    return value


def concerns_of(admet) -> list["ADMETConcern"]:
    """The ranked concerns of a profile or of its serialised form.

    The pipeline stores ``ADMETProfile.to_dict()`` on the candidate, not the
    profile object, and every consumer that reached for ``.concerns`` with
    ``getattr`` got an empty list from a dict and said nothing -- a compound
    with a critical structural alert read exactly like a compound with none.
    Reaching for ``admet["concerns"]`` instead does not fix it either: the
    serialised concerns are plain dicts, so the next line, ``c.severity``,
    raises ``AttributeError``. Both shapes are handled here, once, the way
    :func:`hazard_severity_of` already handles them.

    ``None`` and anything unrecognised come back as an empty list: a profile
    that does not exist has no concerns to report, and the *absence* of an
    assessment is signalled by :func:`hazard_severity_of` returning "unknown",
    not by an empty concern list.
    """
    if admet is None:
        return []
    raw = admet.get("concerns", []) if isinstance(admet, dict) else \
        getattr(admet, "concerns", [])
    out: list[ADMETConcern] = []
    for c in raw or []:
        if isinstance(c, ADMETConcern):
            out.append(c)
        elif isinstance(c, dict):
            out.append(ADMETConcern(
                domain=str(c.get("domain", "")), name=str(c.get("name", "")),
                severity=str(c.get("severity", "")),
                mechanism=str(c.get("mechanism", "")),
                action=str(c.get("action", ""))))
        else:                                  # duck-typed stand-in
            out.append(ADMETConcern(
                domain=str(getattr(c, "domain", "")),
                name=str(getattr(c, "name", "")),
                severity=str(getattr(c, "severity", "")),
                mechanism=str(getattr(c, "mechanism", "")),
                action=str(getattr(c, "action", ""))))
    return out


def hazard_severity_of(admet) -> str:
    """Worst hazard severity in a profile or its serialised form.

    Systemic exposure is excluded deliberately: a large Systemic Exposure Dosage
    is a magnitude, not a liability. It says how much gets in, and whether that
    matters depends on a NOAEL the screen does not have. Common cosmetic
    ingredients can exceed a modeled exposure threshold; exposure magnitude
    alone does not establish hazard or safety without a relevant NOAEL.

    An unassessed profile returns "unknown", never "none". A structure nothing
    could evaluate is not a structure with no hazards.
    """
    if admet is None:
        return "unknown"
    if not isinstance(admet, dict):
        return admet.hazard_severity if admet.assessed else "unknown"
    if not admet.get("assessed", bool(admet.get("verdicts"))):
        return "unknown"
    if admet.get("hazard_severity") in SEVERITY_ORDER:
        return admet["hazard_severity"]
    worst = "none"                            # profile serialised by an older build
    for c in admet.get("concerns", []):
        dom = c.get("domain") if isinstance(c, dict) else getattr(c, "domain", None)
        sev = c.get("severity") if isinstance(c, dict) else getattr(c, "severity", None)
        if dom == "systemic" or sev not in SEVERITY_ORDER:
            continue
        if SEVERITY_ORDER.index(sev) < SEVERITY_ORDER.index(worst):
            worst = sev
    return worst


_EVIDENCE_FROM_PROVENANCE = {
    "human_protein": Evidence.BIOCHEMICAL, "human_cell": Evidence.HUMAN_CELL,
    "nonhuman_cell": Evidence.NONHUMAN_CELL, "surrogate": Evidence.SURROGATE,
    "unknown_origin": Evidence.UNKNOWN, "unspecified": Evidence.UNKNOWN,
}
def evidence_for(lig: LigandRecord) -> Evidence:
    """Where a compound-target claim sits on the evidence ladder."""
    if lig.source == "cosmetic_active":
        # A skin/clinical annotation is not a compound-target assay. When
        # structure de-duplication merges actual observations into a reference,
        # use the selected primary measurement context, never the keeper label.
        primary = lig.activity.get("primary_measurement_group")
        group = next((g for g in lig.activity.get("measurement_groups") or []
                      if g.get("group_id") == primary), None) if primary else None
        if not group:
            return Evidence.UNKNOWN
        return _EVIDENCE_FROM_PROVENANCE.get(group.get("context", {}).get("provenance"),
                                             Evidence.UNKNOWN)
    return _EVIDENCE_FROM_PROVENANCE.get(lig.activity.get("provenance", ""), Evidence.UNKNOWN)
