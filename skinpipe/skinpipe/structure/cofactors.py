'Structural non-protein groups that are part of the pocket, not of the buffer.'
from __future__ import annotations

import logging
import math
import hashlib
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

from .rcsb import HET_BLACKLIST

log = logging.getLogger(__name__)

#: Element -> AutoDock receptor type. Polar atoms take the acceptor forms, which
#: is what `prepare_receptor4.py` writes for a receptor; see the module note on
#: what that leaves out. Vina does not distinguish aromatic from aliphatic
#: carbon (that is an AD4 distinction), so every carbon is `C`.
COFACTOR_AD_TYPE: dict[str, str] = {
    "C": "C", "N": "NA", "O": "OA", "S": "SA", "P": "P",
    "F": "F", "CL": "Cl", "BR": "Br", "I": "I",
}





SUPPORTED_COMPONENTS = frozenset({
    "LFA", "FAD", "FMN", "NAD", "NAP", "NDP",
    "ADP", "ATP", "GDP", "GTP", "SAM", "SAH", "PLP", "TPP", "COA",
})
# Recognised for whole-group detection/MMGBSA refusal, not supported for
# docking: their Fe model is deliberately unavailable.
RECOGNIZED_UNTYPABLE_COMPONENTS = frozenset({"HEM", "HEC", "HEA"})
COFACTOR_POLICY = "source_bound_components_complete_heavy_deposited_H_D_sidecar_v2"
COMPONENT_DECISION_POLICY = "source_bound_structural_component_instance_decisions_v1"
_DECISION_FIELDS = frozenset({"assembly_sha256", "resname", "component_id", "chain", "resnum", "icode",
                              "selected_atom_sha256", "action", "role", "rationale", "evidence_refs"})
_DECISION_ACTIONS = frozenset({"retain_observed_rigid_heavy", "exclude_whole_group"})


def normalize_component_decisions(decisions: Optional[list[dict]]) -> list[dict]:
    """Validate a bounded literal policy; this does not authenticate its evidence.

    No wildcard or component-wide decision exists. Source and atom binding are
    checked against the actual assembled input before a decision can be used.
    Empty/absent decisions retain the default classification policy.
    """
    if decisions is None:
        return []
    if not isinstance(decisions, list) or len(decisions) > 256:
        raise ValueError("structural component decisions require a bounded list")
    output, used = [], set()
    for item in decisions:
        if not isinstance(item, dict) or set(item) != _DECISION_FIELDS:
            raise ValueError("structural component decision has missing/unknown fields")
        for key in ("assembly_sha256", "selected_atom_sha256"):
            if not isinstance(item[key], str) or not re.fullmatch(r"[0-9a-f]{64}", item[key]):
                raise ValueError("structural component decision requires exact SHA256 pins")
        if (not isinstance(item["resname"], str) or not re.fullmatch(r"[A-Z0-9]{1,3}", item["resname"])
                or not isinstance(item["component_id"], str) or not re.fullmatch(r"[A-Z0-9]{1,12}", item["component_id"])
                or not isinstance(item["chain"], str) or not re.fullmatch(r"[A-Za-z0-9 ]", item["chain"])
                or not isinstance(item["resnum"], str) or not re.fullmatch(r"-?[0-9]{1,4}", item["resnum"])
                or not isinstance(item["icode"], str) or not re.fullmatch(r"[A-Za-z0-9]?", item["icode"])):
            raise ValueError("structural component decision requires one exact instance")
        if not isinstance(item["action"], str) or item["action"] not in _DECISION_ACTIONS:
            raise ValueError("unknown structural component action")
        for key, maximum in (("role", 128), ("rationale", 4096)):
            if (not isinstance(item[key], str) or not item[key].strip() or len(item[key]) > maximum
                    or any(ord(c) < 32 for c in item[key])):
                raise ValueError("structural component role/rationale must be bounded literal text")
        refs = item["evidence_refs"]
        if (not isinstance(refs, list) or not 1 <= len(refs) <= 16
                or any(not isinstance(ref, str) or not ref.strip() or len(ref) > 2048
                       or any(ord(c) < 32 for c in ref) for ref in refs)
                or len(set(refs)) != len(refs)):
            raise ValueError("structural component decision needs distinct bounded evidence references")
        key = (item["resname"], item["chain"], item["resnum"], item["icode"])
        if key in used:
            raise ValueError("duplicate/conflicting structural component decision")
        used.add(key)
        output.append({**item, "evidence_refs": list(refs)})
    return sorted(output, key=lambda row: (row["resname"], row["chain"], row["resnum"], row["icode"]))


def component_instance(lines: Iterable[str], *, component_ids: Optional[dict] = None) -> dict:
    """Fingerprint one already-selected raw HETATM group, including deposited H/D.

    The digest binds the exact ordered PDB records (newline-normalized only),
    including atom identity, coordinates, occupancy and selected alternate.
    It does not certify a complete CCD graph, protonation or biological role.
    """
    rows = [line.rstrip("\r\n") for line in lines]
    if not rows or any(not line.startswith("HETATM") or len(line) < 60 for line in rows):
        raise ValueError("nonempty selected raw HETATM group required")
    keys = {_group_key(line) for line in rows}
    names = {(line[12:16], line[16:17]) for line in rows}
    if len(keys) != 1 or len(names) != len(rows):
        raise ValueError("mixed or duplicate structural component atoms")
    for line in rows:
        xyz = [float(line[a:b]) for a, b in ((30, 38), (38, 46), (46, 54))]
        occupancy = float(line[54:60])
        if not all(math.isfinite(x) for x in [*xyz, occupancy]) or not _element_of(line):
            raise ValueError("structural component atom identity/coordinates unavailable")
    resname, chain, resnum, icode = next(iter(keys))
    full = (component_ids or {}).get(resname, resname)
    if not isinstance(full, str) or not re.fullmatch(r"[A-Z0-9]{1,12}", full):
        raise ValueError("structural component full CCD identity is ambiguous")
    return {"resname": resname, "component_id": full, "chain": chain, "resnum": resnum, "icode": icode,
            "selected_atom_sha256": hashlib.sha256(("\n".join(rows) + "\n").encode()).hexdigest()}


def bind_component_decisions(decisions: Optional[list[dict]], atom_records: Iterable[str], *,
                             assembly_sha256: Optional[str], component_ids: Optional[dict] = None) -> list[dict]:
    """Reject stale/missing pins even before receptor-cache lookup."""
    normalized = normalize_component_decisions(decisions)
    if not normalized:
        return []
    if not isinstance(assembly_sha256, str) or not re.fullmatch(r"[0-9a-f]{64}", assembly_sha256):
        raise ValueError("exact assembled input SHA256 required for component decisions")
    from .clean import select_altlocs
    groups = {}
    for line in select_altlocs(atom_records):
        if line.startswith("HETATM"):
            groups.setdefault(_group_key(line), []).append(line)
    for item in normalized:
        key = (item["resname"], item["chain"], item["resnum"], item["icode"])
        if item["assembly_sha256"] != assembly_sha256 or key not in groups:
            raise ValueError("stale/missing structural component source instance")
        actual = component_instance(groups[key], component_ids=component_ids)
        if any(item[k] != actual[k] for k in actual):
            raise ValueError("structural component atom/full-CCD binding differs")
        if item["action"] == "retain_observed_rigid_heavy":
            heavy = [_element_of(line) for line in groups[key] if _element_of(line) not in {"H", "D"}]
            if not heavy or any(element not in COFACTOR_AD_TYPE for element in heavy):
                raise ValueError("retained structural component has unsupported/empty heavy atom model")
    return normalized


class CofactorSelectionError(ValueError):
    """Unknown non-polymer chemistry needs a decision, not silent retention."""

    def __init__(self, evidence: dict):
        self.evidence = evidence
        super().__init__("structural cofactor classification required; preparation held for "
                         + ", ".join(evidence["unclassified_groups"] + evidence["unsupported_atom_groups"]))


def _group_key(line: str) -> tuple[str, str, str, str]:
    return (line[17:20].strip().upper(), line[21], line[22:26].strip(),
            line[26].strip())


def _element_of(line: str) -> str:
    """The element, from columns 77-78 if present, else from the atom name."""
    element = line[76:78].strip().upper()
    if element:
        return element
    name = line[12:16].strip().upper()
    return name[:2] if name[:2] in COFACTOR_AD_TYPE else name[:1]


def select_structural_cofactors(
    atom_records: "Iterable[str]",
    chain: Optional[str] = None,
    exclude: "Optional[tuple]" = None,
    *, polymer_residues: "Optional[set[tuple[str, str, str, str]]]" = None,
    evidence: Optional[dict] = None,
    structural_component_decisions: Optional[list[dict]] = None,
    assembly_sha256: Optional[str] = None,
    component_ids: Optional[dict] = None,
) -> list[str]:
    """HETATM lines that are part of the structure rather than of the buffer.

    `exclude` is the `(resname, chain, resnum, icode)` of the ligand this
    receptor will dock, which must not be in the receptor it is docked into.
    """
    return _select_cofactor_records(atom_records, chain, exclude,
                                   polymer_residues=polymer_residues,
                                   evidence=evidence, require_docking_types=True,
                                   structural_component_decisions=structural_component_decisions,
                                   assembly_sha256=assembly_sha256, component_ids=component_ids)


def _select_cofactor_records(
    atom_records: "Iterable[str]", chain: Optional[str] = None,
    exclude: "Optional[tuple]" = None, *,
    polymer_residues: "Optional[set[tuple[str, str, str, str]]]" = None,
    evidence: Optional[dict] = None, require_docking_types: bool,
    structural_component_decisions: Optional[list[dict]] = None,
    assembly_sha256: Optional[str] = None, component_ids: Optional[dict] = None,
) -> list[str]:
    """Classify components; optionally require their complete docking model.

    Inspection of an existing receptor must see recognised cofactors even when
    the PDBQT writer cannot represent them. Preparation takes the strict path
    above and holds instead of silently emitting a partial group.
    """
    from .clean import is_metal, is_water, select_altlocs
    from . import modified

    selected = select_altlocs(atom_records)
    decisions = bind_component_decisions(structural_component_decisions, selected,
                                        assembly_sha256=assembly_sha256, component_ids=component_ids)
    by_instance = {(d["resname"], d["chain"], d["resnum"], d["icode"]): d for d in decisions}
    native = modified.groups("\n".join(selected), known_phospho=True)
    polymer_residues = set(polymer_residues or ()) | {
        (name, cid, number[:-1] if number[-1:].isalpha() else number,
         number[-1] if number[-1:].isalpha() else "") for cid, number, name in native}
    kept: list[str] = []
    unknown: set[str] = set()
    untypable: set[str] = set()
    eligible = {}
    for line in selected:
        if not line.startswith("HETATM"):
            continue
        if chain is not None and line[21] != chain:
            continue
        resname = line[17:20].strip().upper()
        if is_water(resname) or is_metal(resname) or resname in HET_BLACKLIST:
            continue
        if exclude is not None and _group_key(line) == tuple(exclude):
            continue
        if polymer_residues is not None and _group_key(line) in polymer_residues:
            continue
        eligible.setdefault(_group_key(line), []).append(line)
    if by_instance.keys() - eligible.keys():
        raise ValueError("component decision is outside structural-group scope (selected ligand/polymer/buffer/metal/chain)")
    excluded = []
    for key, lines in eligible.items():
        resname, cid, number, icode = key
        decision = by_instance.get(key)
        if decision and decision["action"] == "exclude_whole_group":
            excluded.append({**component_instance(lines, component_ids=component_ids),
                             "atom_count": len(lines), "action": decision["action"]})
            continue
        if resname not in SUPPORTED_COMPONENTS | RECOGNIZED_UNTYPABLE_COMPONENTS and decision is None:
            unknown.add(f"{resname} {cid}:{number}{icode}")
            continue
        # Deposited H/D stay in the raw sidecar. The declared rigid heavy model
        # intentionally omits them at rendering; it never invents donor typing.
        heavy = [_element_of(line) for line in lines if _element_of(line) not in {"H", "D"}]
        unsupported = set(heavy) - COFACTOR_AD_TYPE.keys()
        if require_docking_types and (unsupported or not heavy):
            for element in unsupported or {"empty heavy model"}:
                untypable.add(f"{resname} {cid}:{number}{icode} ({element})")
            continue
        kept.extend(lines)
    record = {"policy": COFACTOR_POLICY, "unclassified_groups": sorted(unknown),
              "unsupported_atom_groups": sorted(untypable),
              "recognised_but_untypable_components": sorted({key[0] for key in eligible if key[0] in RECOGNIZED_UNTYPABLE_COMPONENTS}),
              "retained_groups": sorted({f"{line[17:20].strip()} {line[21]}:{line[22:27].strip()}"
                                         for line in kept}),
              "pocket_role_established": False,
              "deposited_H_D_records_omitted_from_typing": sum(_element_of(line) in {"H", "D"} for line in kept),
              "rigid_model": "observed heavy atoms, no donor hydrogens, zero partial charges; chemistry not qualified",
              "status": "held" if unknown or untypable else "classified_components"}
    if decisions:
        record.update(structural_component_decisions=decisions, excluded_whole_groups=excluded,
                      decision_policy=COMPONENT_DECISION_POLICY,
                      decisions_sha256=hashlib.sha256(json.dumps(decisions, sort_keys=True, separators=(",", ":")).encode()).hexdigest(),
                      manual_model_status="EXPLICIT_MANUAL_MODEL_NOT_CHEMICALLY_QUALIFIED",
                      native_eligibility="NOT_ESTABLISHED", decision_evidence_authenticated=False)
    if evidence is not None:
        evidence.update(record)
    if unknown or untypable:
        raise CofactorSelectionError(record)
    return kept


def render_cofactor_records(
    cofactor_lines: list[str], start_serial: int = 95000,
) -> tuple[list[str], list[str]]:
    """Render retained groups as rigid PDBQT HETATM records, plus honest notes.

    Every offered heavy atom must be representable; no partial group is emitted.
    Deposited H/D are deliberately omitted under the existing no-donor model,
    while the unmodified source records remain in the coordinate sidecar.
    """
    cofactor_lines = list(cofactor_lines)
    records: list[str] = []
    notes: list[str] = []
    written: dict[str, int] = {}
    hydrogen_count = sum(_element_of(src) in {"H", "D"} for src in cofactor_lines)
    bad = sorted({f"{src[17:20].strip()} {src[21]}:{src[22:27].strip()} ({_element_of(src) or '?'})"
                  for src in cofactor_lines if _element_of(src) not in {*COFACTOR_AD_TYPE, "H", "D"}})
    if bad or (cofactor_lines and hydrogen_count == len(cofactor_lines)):
        raise CofactorSelectionError({"unclassified_groups": [], "unsupported_atom_groups": bad or ["empty heavy model"]})
    serial = start_serial
    for src in cofactor_lines:
        element = _element_of(src)
        if element in {"H", "D"}:
            continue
        ad_type = COFACTOR_AD_TYPE[element]
        resname = src[17:20].strip().upper()
        serial += 1
        written[resname] = written.get(resname, 0) + 1
        name = f"{src[12:16].strip().upper():<4}"[:4]
        records.append(
            f"HETATM{serial:5d} {name} {resname:>3s} {src[21]}{src[22:26]}"
            f"{src[26]}   {src[30:54]}  1.00  0.00    {0.0:+6.3f} {ad_type}"
        )
    if written:
        summary = ", ".join(f"{k} ({v} atoms)" for k, v in sorted(written.items()))
        notes.append(
            f"kept {sum(written.values())} structural non-protein atom(s) in the "
            f"receptor: {summary}. Their role in the selected pocket is not "
            "established by component identity. They carry no hydrogens, so they "
            "accept hydrogen bonds and never donate one.")
    if hydrogen_count:
        notes.append(
            f"intentionally omitted {hydrogen_count} deposited cofactor H/D atom(s) from rigid "
            "heavy-atom typing; raw sidecar records preserved; no donor hydrogens or protonation qualification")
    return records, notes


#: Suffix of the record written beside a receptor PDB naming the structural
#: cofactors that are in the *Vina* PDBQT and not in that PDB.
COFACTOR_SIDECAR_SUFFIX = "_cofactors.pdb"


def sidecar_path(receptor_pdb: Path) -> Path:
    """Where ``receptor_pdb``'s held-out structural cofactors are recorded."""
    path = Path(receptor_pdb)
    return path.with_name(path.stem + COFACTOR_SIDECAR_SUFFIX)


def clear_sidecar(receptor_pdb: Path) -> bool:
    """Delete a stale cofactor record. True when one was removed."""
    path = sidecar_path(receptor_pdb)
    if path.is_file():
        path.unlink()
        return True
    return False


def write_sidecar(receptor_pdb: Path,
                  cofactor_lines: "Iterable[str]") -> Optional[Path]:
    """Record the cofactors that are in the PDBQT and not in ``receptor_pdb``.

    Same contract as :func:`skinpipe.structure.metal_sites.write_sidecar`, and
    for the same reason one level along. `prepare_receptor` writes a
    protein-only PDB and re-attaches structural cofactors to the Vina PDBQT;
    ``Receptor.path`` is the protein-only file, and it is what the pipeline
    hands to MM/GBSA -- which then runs ``removeHeterogens`` and would have
    deleted them anyway. So a pose docked into a pocket one wall of which is a
    cofactor gets rescored against a pocket that no longer has that wall, and
    nothing in either number says so.

    Call it unconditionally: an existing record is deleted when the new one is
    empty, and guarding the call with ``if cofactor_lines:`` leaves a stale file
    that would make the guard refuse a receptor that has no cofactor.
    """
    path = sidecar_path(receptor_pdb)
    lines = [pdb_line.rstrip("\n") for pdb_line in cofactor_lines
             if pdb_line.startswith(("ATOM", "HETATM"))]
    if not lines:
        clear_sidecar(receptor_pdb)
        return None
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        f"REMARK   1 structural cofactors held out of {Path(receptor_pdb).name}\n"
        "REMARK   1 written by skinpipe.structure.cofactors; read by\n"
        "REMARK   1 skinpipe.mmgbsa.cofactor_guard, which refuses to score a\n"
        "REMARK   1 pose whose pocket wall one of them forms.\n"
        + "\n".join(lines) + "\nEND\n"
    )
    return path


@dataclass(frozen=True)
class CofactorAtom:
    resname: str
    chain: str
    resnum: str
    xyz: tuple[float, float, float]

    @property
    def label(self) -> str:
        return f"{self.resname} {self.chain}:{self.resnum}".strip()


def parse_cofactor_atoms(lines: "Iterable[str]") -> list[CofactorAtom]:
    """Coordinate records to values. An unparseable record raises.

    Skipping it would let the guard conclude the pocket is cofactor-free, which
    is the conclusion this path exists to prevent.
    """
    out: list[CofactorAtom] = []
    for line in lines:
        if not line.startswith(("ATOM", "HETATM")):
            continue
        try:
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            if not all(math.isfinite(value) for value in xyz):
                raise ValueError("non-finite cofactor coordinates")
        except (ValueError, IndexError) as exc:
            raise ValueError(
                f"unparseable cofactor coordinate record: {line[:54]!r}") from exc
        out.append(CofactorAtom(resname=line[17:20].strip().upper(),
                                chain=line[21], resnum=line[22:27].strip(),
                                xyz=xyz))
    return out


def receptor_cofactors(receptor_pdb: Path) -> list[CofactorAtom]:
    """Structural cofactors for ``receptor_pdb``: its sidecar and its own body.

    The sidecar is what `structure.prepare` writes when it lifts a cofactor out
    of an entry, and it is the only place this function used to look. That is
    right for a receptor this pipeline prepared and wrong for one handed
    straight to `mmgbsa_single`/`mmgbsa_batch`: a raw holo PDB carries its
    cofactors as HETATM in the body, `[]` came back, the guard had nothing to
    refuse, and PDBFixer's `removeHeterogens()` then deleted them -- rescoring
    the pose in a pocket whose wall had gone.

    Both sources are read and merged. Body records go through
    `select_structural_cofactors` first -- the same selector `structure.prepare`
    uses to build a sidecar -- so waters, metals and `rcsb.HET_BLACKLIST`
    crystallisation additives are excluded and not every HETATM becomes a
    cofactor. Metals are deliberately not here: `mmgbsa.metal_guard` owns them
    and refuses on its own terms.
    A record present in both is counted once, keyed by residue name, chain,
    residue number and coordinates -- `CofactorAtom` carries no atom name, and
    the coordinates are what distinguish two atoms of one residue.
    """
    atoms: list[CofactorAtom] = []
    seen: set[tuple] = set()
    path = sidecar_path(receptor_pdb)
    for lines in ((path.read_text().splitlines() if path.is_file() else []),
                  _body_lines(receptor_pdb)):
        for a in parse_cofactor_atoms(lines):
            key = (a.resname, a.chain, a.resnum, a.xyz)
            if key in seen:
                continue
            seen.add(key)
            atoms.append(a)
    return atoms


def _body_lines(receptor_pdb: Path) -> list[str]:
    """Structural-cofactor HETATM records in the receptor file itself."""
    try:
        text = receptor_pdb.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return []
    from . import modified
    from .clean import atom_lines, select_altlocs

    # Native writers can retain MODRES or omit it for a known phospho residue.
    # The latter still needs a peptide backbone, not merely the component ID.
    selected = select_altlocs(atom_lines(text))
    native_text = "\n".join(modified.declaration_lines(text) + selected)
    polymer = {(name, cid, number[:-1] if number[-1:].isalpha() else number,
                number[-1] if number[-1:].isalpha() else "")
               for cid, number, name in modified.groups(native_text, known_phospho=True)}
    return _select_cofactor_records(selected, polymer_residues=polymer,
                                    require_docking_types=False)
