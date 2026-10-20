"""PDB text surgery, PDBFixer repair, and receptor PDBQT generation.

Why the receptor PDBQT is built in two pieces
---------------------------------------------
Meeko 0.8 replaced the old ``prepare_receptor4.py`` flow with a chemical-template
engine (``Polymer`` + ``ResidueChemTemplates``). That engine only knows about
biopolymer residues: a bare ``HETATM ... CU`` has no template, so
``PDBQTWriterLegacy.write_string_from_polymer`` dies with
``Can't write PDBQT because atom_type is None``.

Silently dropping the metal is not an option - tyrosinase (di-copper) and the
MMPs (catalytic Zn) have pockets that are *defined* by the metal, and removing it
leaves a hole that Vina happily fills with nonsense. So we:

1. strip catalytic metals out of the PDB before handing it to Meeko,
2. let Meeko type and write the protein,
3. append the metals back as explicit PDBQT ``HETATM`` records with AutoDock
   element types (``Zn``, ``Cu``, ...) and formal charges.

Vina 1.2 maps those element types onto its internal ``Met_D`` atom type, which
is exactly what ``prepare_receptor4.py`` used to produce. Verified against
Vina 1.2.7: ``set_receptor`` accepts the resulting file.

A second Meeko caveat: residues with unexpected inter-residue bonds (the
Cys-His thioether crosslink in mushroom tyrosinase, for example) fail template
matching. ``forgive_extra_bonds=True`` then crashes inside the padding builder,
so template matching can report ignored residues. Observed heavy atoms are
checked against the typed output before any file is written; a partial observed
receptor is rejected. Native unusual covalent chemistry needs separate review.
"""
from __future__ import annotations

import logging
import math
import random
from pathlib import Path
from typing import Any, Iterable, Optional

from .metals import render_metal_records
from .rcsb import CATALYTIC_METALS
from .altloc_ties import resolve_occupancy_ties as _resolve_occupancy_ties





























PDBFIXER_SEED: int = 20260929



















PDBFIXER_PLATFORM: str = "Reference"


def validate_ph(value: Any) -> float:
    """Require a finite pH before molecular preparation or output."""
    if isinstance(value, bool):
        raise ValueError("ph must be a finite number")
    try:
        ph = float(value)
    except (TypeError, ValueError, OverflowError):
        raise ValueError("ph must be a finite number") from None
    if not math.isfinite(ph):
        raise ValueError("ph must be a finite number")
    if not 0 <= ph <= 14:
        raise ValueError("ph must be in the range 0 to 14")
    return ph


def preparation_model(platform: Optional[str] = None, cpu_threads: int = 1, *,
                      rebuild_missing_residues: bool = False) -> dict:
    """An explicit protocol; accelerated CPU is never a silent fallback."""
    if type(rebuild_missing_residues) is not bool:
        raise ValueError("rebuild_missing_residues must be boolean")
    if rebuild_missing_residues:
        raise ValueError("loop reconstruction requires a separately validated preparation path; held")
    name = PDBFIXER_PLATFORM if platform is None else platform
    if name not in ("Reference", "CPU"):
        raise ValueError("preparation platform must be Reference or CPU")
    if type(cpu_threads) is not int or not 1 <= cpu_threads <= 16:
        raise ValueError("preparation CPU threads must be an integer from 1 to 16")
    if name == "Reference" and cpu_threads != 1:
        raise ValueError("Reference preparation requires one thread")
    return {"policy": "explicit_preparation_platform_v1", "platform": name,
            "seed": PDBFIXER_SEED,
            "cpu_threads": cpu_threads if name == "CPU" else None,
            "deterministic_forces": True if name == "CPU" else None,
            "missing_residue_policy": "preserve_unobserved_gaps_v1",
            "observed_backbone_gap_policy": "numbering_and_geometry_fragment_guard_v1",
            "rebuild_missing_residues": False,
            "reproducibility_scope": "prepared artifact; same-host repeat must be verified"}


def _fixer_platform(name: Optional[str] = None, cpu_threads: int = 1):
    """Resolve the recorded platform; an unavailable default cannot substitute."""
    requested = name or PDBFIXER_PLATFORM
    try:
        import openmm
        platform = openmm.Platform.getPlatformByName(requested)
        if name == "CPU":
            platform.setPropertyDefaultValue("Threads", str(cpu_threads))
            platform.setPropertyDefaultValue("DeterministicForces", "true")
        return platform
    except Exception as exc:                                  # noqa: BLE001
        # preparation_model records the default as Reference just as it records
        # an explicit request. Refuse substitution: returning None would let
        # PDBFixer choose a platform that differs from the recorded identity.
        raise RuntimeError(f"requested preparation platform {requested} unavailable") from exc

log = logging.getLogger(__name__)

_WATERS = frozenset({"HOH", "DOD", "WAT"})

def is_water(resname: str) -> bool:
    return resname.strip().upper() in _WATERS


def is_metal(resname: str) -> bool:
    return resname.strip().upper() in CATALYTIC_METALS


def atom_lines(pdb_text: str) -> list[str]:
    """ATOM/HETATM records of the first model only."""
    out: list[str] = []
    for line in pdb_text.splitlines():
        if line.startswith("ENDMDL"):
            break
        if line.startswith(("ATOM  ", "HETATM")):
            out.append(line)
    return out


def chain_ids(pdb_text: str) -> list[str]:
    """Polymer chain identifiers, in order of first appearance."""
    seen: list[str] = []
    for line in atom_lines(pdb_text):
        if not line.startswith("ATOM  "):
            continue
        cid = line[21]
        if cid not in seen:
            seen.append(cid)
    return seen


ALTLOC_POLICY = "residue_max_occupancy_fixed_peptide_final_edge_guard_v3"


def select_altlocs(records: Iterable[str], *, evidence: Optional[dict] = None,
                   peptide_context: Optional[str] = None) -> list[str]:
    """Choose one residue conformer by mean heavy-atom occupancy.

    Blank atoms are shared. Exact occupancy ties may resolve a unique
    compatible choice against fixed adjacent peptide anchors in full TER
    context; otherwise ties use A, then lexicographic label. Atoms from
    losing labelled conformers never fill holes in the selected conformer.
    Final evaluable peptide edges touching ties must satisfy the same local
    compatibility bounds or selection holds; this is not geometry validation.
    Selected labels are cleared for downstream readers that accept only A.
    Missing/invalid occupancies in alternate records hold selection rather
    than silently asserting that an arbitrary conformer is the major one.
    """
    rows = list(records)
    groups: dict[tuple[str, str, str], list[int]] = {}
    for i, line in enumerate(rows):
        if line.startswith(("ATOM  ", "HETATM")):
            groups.setdefault((line[21], line[22:27], line[17:20]), []).append(i)
    positions: dict[tuple[str, str], set[str]] = {}
    for chain, number, residue in groups:
        positions.setdefault((chain, number), set()).add(residue)
    for (chain, number), names in positions.items():
        if len(names) > 1 and any(rows[i][16] != " " for name in names
                                  for i in groups[(chain, number, name)]):
            raise ValueError(f"alternate residue identities at {chain}:{number.strip()}; selection held")
    selected: set[int] = set()
    decisions = []
    for key, indices in groups.items():
        labels = sorted({rows[i][16] for i in indices} - {" "})
        if not labels:
            selected.update(indices)
            continue
        occupancy = {}
        for i in indices:
            try:
                value = float(rows[i][54:60])
                if not math.isfinite(value) or not 0 <= value <= 1:
                    raise ValueError
            except ValueError:
                raise ValueError(f"invalid alternate-location occupancy for {key}; selection held") from None
            occupancy[i] = value
        scores = {}
        for label in labels:
            atoms = [i for i in indices if rows[i][16] == label
                     and rows[i][76:78].strip().upper() not in {"H", "D"}]
            if not atoms:
                atoms = [i for i in indices if rows[i][16] == label]
            scores[label] = math.fsum(occupancy[i] for i in atoms) / len(atoms)
        chosen = min(labels, key=lambda label: (-scores[label], label != "A", label))
        if scores[chosen] <= 0:
            raise ValueError(f"alternate conformers have zero occupancy for {key}; selection held")
        # Shared blank atoms take precedence on an exact occupancy tie. This
        # also avoids duplicate identifiers in unusual blank/labelled input.
        by_atom: dict[str, int] = {}
        for i in indices:
            if rows[i][16] not in {" ", chosen}:
                continue
            name = rows[i][12:16]
            if name in by_atom:
                previous = by_atom[name]
                if rows[i][16] == rows[previous][16]:
                    raise ValueError(f"duplicate atom in alternate conformer {key}:{name}; selection held")
                winner = min((previous, i), key=lambda n: (-occupancy[n], rows[n][16] != " "))
                by_atom[name] = winner
            else:
                by_atom[name] = i
        selected.update(by_atom.values())
        decisions.append({"chain": key[0], "residue": key[1].strip(), "name": key[2].strip(),
                          "selected_altloc": chosen, "mean_heavy_occupancy": scores,
                          "discarded_serials": [rows[i][6:11].strip() for i in indices
                                                 if i not in by_atom.values()]})
    selected, tie_evidence = _resolve_occupancy_ties(
        rows, groups, selected, decisions, peptide_context)
    dropped = [row for i, row in enumerate(rows)
               if row.startswith(("ATOM  ", "HETATM")) and i not in selected]
    if evidence is not None:
        evidence.update(policy=ALTLOC_POLICY, decisions=decisions,
                        adjacent_peptide_tie_break=tie_evidence,
                        discarded_record_count=len(dropped),
                        discarded_heavy_atom_record_count=sum(row[76:78].strip().upper() not in {"H", "D"}
                                                              for row in dropped),
                        discarded_metal_record_count=sum(is_metal(row[17:20]) for row in dropped))
    return [row[:16] + " " + row[17:] for i, row in enumerate(rows) if i in selected]


def resolved_altloc_text(pdb_text: str) -> str:
    """The same selection for box/native-chemistry readers; keep record order."""
    source = atom_lines(pdb_text)
    if len({line[6:11] for line in source}) != len(source):
        raise ValueError("duplicate source atom serials; alternate-location selection held")
    selected = {line[6:11]: line for line in select_altlocs(source, peptide_context=pdb_text)}
    output = []
    for line in pdb_text.splitlines():
        if line.startswith("ENDMDL"):
            break
        if line.startswith(("ATOM  ", "HETATM")):
            if line[6:11] in selected:
                output.append(selected[line[6:11]])
        else:
            output.append(line)
    return "\n".join(output) + "\n"


def select_receptor_lines(
    pdb_text: str,
    chain: Optional[str],
    keep_waters: bool = False,
    keep_metals: bool = True,
    *, selection_evidence: Optional[dict] = None,
) -> tuple[list[str], list[str]]:
    """Split a PDB into (protein lines, observed metal-name-list lines).

    Alternate locations use :func:`select_altlocs`, before any atom is dropped.
    """
    protein: list[str] = []
    metals: list[str] = []
    from . import modified
    native = modified.declarations(pdb_text, chain)
    source = [line for line in atom_lines(pdb_text) if chain is None or line[21] == chain]
    alternate_evidence: dict = {}
    for line in select_altlocs(source, evidence=alternate_evidence, peptide_context=pdb_text):
        resname = line[17:20].strip().upper()
        het = line.startswith("HETATM")
        if is_water(resname):
            if keep_waters:
                protein.append(line)
            continue
        if is_metal(resname):
            if keep_metals:
                metals.append(line)
            continue
        if het:
            if (line[21], line[22:27].strip(), resname) in native:
                protein.append(line)
                continue
            continue  # every other HETATM (ligand, buffer, cryo) is discarded
        protein.append(line)
    if selection_evidence is not None:
        selection_evidence.update(altloc=alternate_evidence,
                                  observed_metal_record_count=sum(is_metal(line[17:20]) for line in source),
                                  retained_metal_atom_count=len(metals), keep_metals=keep_metals,
                                  metal_selection_status=("retained" if metals else
                                                          "excluded_by_request" if not keep_metals and any(
                                                              is_metal(line[17:20]) for line in source)
                                                          else "none_observed"))
    return protein, metals


#: A residue needs these to be a real residue rather than a placeholder.
_CORE_BACKBONE = frozenset({"N", "CA", "C"})


def trim_ragged_termini(protein: list[str]) -> tuple[list[str], list[str]]:
    """Drop terminal residues that the crystallographer could not resolve.

    Chain ends fray. A terminal residue is often represented by one or two
    backbone atoms because the rest is disordered, and PDBFixer will happily
    rebuild the remainder from that fragment. 3ERT ends on a proline whose only
    observed atom is N; the rebuilt residue is entirely invented, it contributes
    nothing to any binding site, and it broke receptor typing outright.

    Only termini are trimmed. Missing atoms within an observed internal residue
    may be repaired; a wholly unobserved loop remains a gap under the default
    preparation policy.
    """
    if not protein:
        return protein, []

    order: list[tuple[str, str]] = []
    atoms: dict[tuple[str, str], set[str]] = {}
    for line in protein:
        if is_water(line[17:20].strip()):
            continue  # retained water is not an incomplete terminal residue
        key = (line[21], line[22:27])            # chain, resseq+icode
        if key not in atoms:
            atoms[key] = set()
            order.append(key)
        atoms[key].add(line[12:16].strip().upper())

    complete = {k for k, a in atoms.items() if _CORE_BACKBONE <= a}
    drop: set[tuple[str, str]] = set()
    notes: list[str] = []

    by_chain: dict[str, list[tuple[str, str]]] = {}
    for k in order:
        by_chain.setdefault(k[0], []).append(k)

    for ch, keys in by_chain.items():
        for end, seq in (("N", keys), ("C", list(reversed(keys)))):
            for k in seq:
                if k in complete:
                    break
                drop.add(k)
                missing = ",".join(sorted(_CORE_BACKBONE - atoms[k]))
                notes.append(
                    f"trimmed {end}-terminal residue {ch}:{k[1].strip()} "
                    f"(missing {missing}); rebuilding it from {len(atoms[k])} "
                    f"observed atom(s) would have invented the geometry")

    if not drop:
        return protein, []
    return [pdb_line for pdb_line in protein if is_water(pdb_line[17:20].strip())
            or (pdb_line[21], pdb_line[22:27]) not in drop], notes


def seqres_lines(pdb_text: str, chain: Optional[str]) -> list[str]:
    """SEQRES records for ``chain``; PDBFixer needs them to spot missing loops."""
    out = []
    for line in pdb_text.splitlines():
        if line.startswith("SEQRES") and (chain is None or line[11] == chain):
            out.append(line)
    return out


def write_pdb(lines: Iterable[str], path: Path) -> Path:
    body = "\n".join(lines)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(body + "\nTER\nEND\n")
    return path


def run_pdbfixer(
    in_pdb: Path,
    out_pdb: Path,
    ph: float = 7.4,
    keep_waters: bool = False,
    *,
    platform: Optional[str] = None,
    cpu_threads: int = 1,
    rebuild_missing_residues: bool = False,
) -> list[str]:
    """Rebuild missing heavy atoms/residues and protonate. Returns human notes.

    Entire unobserved residues remain gaps. Repair of missing atoms within
    observed residues is separate; observed heavy coordinates are preserved.
    Explicit loop requests are held until a validated reconstruction path exists.
    """
    ph = validate_ph(ph)
    model = preparation_model(platform, cpu_threads, rebuild_missing_residues=rebuild_missing_residues)
    import openmm.app as app
    import pdbfixer
    from .diagnostics import PreparationJournal

    notes: list[str] = []
    from . import modified
    input_text = Path(in_pdb).read_text()
    journal = PreparationJournal(in_pdb, out_pdb, {**model,
                                                 "ph": float(ph), "keep_waters": keep_waters})
    native = modified.groups(input_text)
    with journal.stage("load_and_register_templates"):
        if platform is None:
            fixer = pdbfixer.PDBFixer(filename=str(in_pdb), platform=_fixer_platform())
        else:
            fixer = pdbfixer.PDBFixer(filename=str(in_pdb), platform=_fixer_platform(platform, cpu_threads))
        if native:
            notes.extend(modified.register_templates(fixer, native))

    with journal.stage("find_missing_residues") as stage:
        fixer.findMissingResidues()
        gaps = _preserve_missing_residue_gaps(fixer, input_text)
        stage.update(gaps=gaps, policy=model["missing_residue_policy"], reconstructed_residue_count=0)
    notes.extend(_missing_residue_gap_notes(gaps))
    if not gaps:
        # PDBFixer can only locate gaps when SEQRES aligns exactly onto the
        # modelled residue-number span. Engineered deletions (e.g. the excised
        # kinase-insert domain in 1T46) break that assumption, and it silently
        # reports no gaps. Say so rather than implying the chain is complete.
        modelled = fixer.topology.getNumResidues()
        declared = sum(len(s.residues) for s in fixer.sequences) or 0
        if not declared:
            # No SEQRES at all, which is the shape `mmcif.cif_to_pdb` produces:
            # it writes the file from `_atom_site` alone. Without it PDBFixer
            # has nothing to align against, so `findMissingResidues` finds
            # nothing AND the branch below cannot fire -- both the detection and
            # the caveat that exists to say "I could not check" go silent
            # together, on exactly the large structures the mmCIF fallback is
            # there to handle, where unmodelled loops are most likely. Measured
            # on 3ERT chain A: with SEQRES, 261 declared against 246 modelled;
            # without it, no note at all.
            notes.append(
                f"no SEQRES in the input, so completeness could not be checked: "
                f"{modelled} residue(s) are modelled and pdbfixer had nothing to "
                "align them against. Unmodelled loops bordering the pocket would "
                "not be reported")
        elif declared > modelled:
            notes.append(
                f"SEQRES declares {declared} residues, {modelled} modelled; pdbfixer "
                "could not align the gaps so none were rebuilt (check for loops "
                "bordering the pocket)"
            )

    with journal.stage("find_nonstandard_residues"):
        fixer.findNonstandardResidues()
    if native:
        fixer.nonstandardResidues = [item for item in fixer.nonstandardResidues
                                    if (item[0].chain.id, item[0].id + item[0].insertionCode.strip(),
                                        item[0].name) not in native]
    if fixer.nonstandardResidues:
        names = ", ".join(sorted({r[0].name for r in fixer.nonstandardResidues}))
        notes.append(f"replaced nonstandard residue(s): {names}")
        conversions = [{"chain_id": residue.chain.id,
                        "residue_id": residue.id + residue.insertionCode.strip(),
                        "source_residue_name": residue.name, "prepared_residue_name": replacement,
                        "chemistry_changed": True}
                       for residue, replacement in fixer.nonstandardResidues]
        for item in conversions:
            notes.append(f"nonstandard residue chemistry changed: {item['chain_id']}:{item['residue_id']} "
                         f"{item['source_residue_name']} -> {item['prepared_residue_name']}; "
                         "retained-atom coordinate preservation does not preserve modified-residue chemistry")
        with journal.stage("replace_nonstandard_residues") as stage:
            stage["residue_conversions"] = conversions
            fixer.replaceNonstandardResidues()

    # A legacy MSE->MET conversion creates standard bonds after the first gap
    # scan. Check the actual resolved residue chemistry before atom repair;
    # retained unusual polymer chemistry with a long C--N edge needs review.
    with journal.stage("verify_resolved_backbone_fragments") as stage:
        additional = _fragment_resolved_backbone(fixer, input_text)
        # The find_missing_residues journal retains this list, so evidence is
        # amended once rather than duplicated across stage.gaps records.
        gaps.extend(additional)
        stage.update(additional_fragment_boundaries=additional,
                     observed_backbone_gap_policy=model["observed_backbone_gap_policy"])
    notes.extend(_missing_residue_gap_notes(additional))

    # Freeze the heavy atoms present before reconstruction, including atoms
    # renamed by a standard-residue conversion. Repair may move new atoms only.
    observed = {_topology_atom_key(a) for a in fixer.topology.atoms()
                if a.element != app.element.hydrogen}
    from openmm import unit
    positions = fixer.positions.value_in_unit(unit.nanometer)
    required_observed = {_topology_atom_key(a): tuple(positions[a.index]) for a in fixer.topology.atoms()
                         if a.element != app.element.hydrogen and
                         (a.residue.name in modified.STANDARD_AMINO_ACIDS
                          or _topology_atom_key(a)[:3] in native
                          or keep_waters and is_water(a.residue.name))}

    with journal.stage("find_missing_atoms"):
        fixer.findMissingAtoms()
    if fixer.missingAtoms:
        notes.append(f"pdbfixer added heavy atoms to {len(fixer.missingAtoms)} residue(s)")
    if fixer.missingTerminals:
        notes.append(f"pdbfixer capped {len(fixer.missingTerminals)} terminal(s)")
    # Both repair steps below are stochastic; see PDBFIXER_SEED. The global
    # RNG is restored afterwards so nothing else in the process inherits a
    # pinned one.
    _rng_state = random.getstate()
    random.seed(PDBFIXER_SEED)
    try:
        with journal.stage("add_missing_atoms"):
            fixer.addMissingAtoms(seed=PDBFIXER_SEED)
        if native:
            notes.extend(modified.deduplicate_repair_bonds(fixer))
            modified.validate_native_bonds(fixer.topology, native)

        notes.extend(_repair_terminal_oxygens(fixer, observed))

        with journal.stage("remove_other_heterogens"):
            if not keep_waters:
                if native:
                    modified.remove_other_heterogens(fixer, native)
                else:
                    fixer.removeHeterogens(keepWater=False)
        with journal.stage("add_missing_hydrogens"):
            if native:
                modified.add_hydrogens(fixer, native, float(ph))
            else:
                fixer.addMissingHydrogens(float(ph))
    finally:
        random.setstate(_rng_state)
    notes.append(f"protonated at pH {ph}")

    with journal.stage("relax_rebuilt_atoms"):
        notes.extend(_relax_rebuilt_atoms(fixer, observed))

    with journal.stage("record_reconstructed_heavy_atoms") as stage:
        reconstructed = reconstructed_heavy_atoms(fixer, observed)
        stage.update(atoms=reconstructed, atom_count=len(reconstructed),
                     coordinate_unit="angstrom", coordinates_observed=False,
                     rotamer_or_pocket_validation="NOT_ESTABLISHED")
    if reconstructed:
        groups = {}
        for atom in reconstructed:
            groups.setdefault((atom['chain_id'], atom['residue_id'], atom['residue_name']), []).append(atom['atom_name'])
        names = '; '.join(f"{chain}:{number} {name} {','.join(atoms)}"
                          for (chain, number, name), atoms in groups.items())
        notes.append(f"reconstructed {len(reconstructed)} unobserved heavy atom(s): {names}; "
                     "positions are modelled, not deposited observations; side-chain rotamers and pocket geometry are unvalidated")

    with journal.stage("verify_observed_coordinates") as stage:
        _verify_observed_coordinates(fixer, required_observed)
        stage.update(observed_heavy_atom_count=len(required_observed), coordinate_unit="nanometer",
                     maximum_allowed_displacement=1e-9)
    notes.append(f"preserved coordinates of all {len(required_observed)} observed heavy atoms; "
                 "whole-residue gaps were not reconstructed")

    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    with journal.stage("write_clean_coordinates"), out_pdb.open("w") as handle:
        if native:
            handle.write("\n".join(modified.declaration_lines(input_text)) + "\n")
        app.PDBFile.writeFile(fixer.topology, fixer.positions, handle, keepIds=True)
    journal.finish(out_pdb)
    return notes


def reconstructed_heavy_atoms(fixer, observed) -> list[dict]:
    """Record actual final added heavy atoms; never move or infer observations."""
    from openmm import unit
    coordinates = fixer.positions.value_in_unit(unit.angstrom)
    atoms = []
    for atom in fixer.topology.atoms():
        key = _topology_atom_key(atom)
        if atom.element is None or atom.element.symbol in {'H', 'D'} or key in observed:
            continue
        atoms.append({'chain_id': key[0], 'residue_id': key[1], 'residue_name': key[2],
                      'atom_name': key[3], 'element': atom.element.symbol,
                      'xyz_angstrom': [float(v) for v in coordinates[atom.index]],
                      'observed': False})
    return atoms


def _topology_atom_key(atom) -> tuple[str, str, str, str]:
    res = atom.residue
    return (res.chain.id, res.id + res.insertionCode.strip(), res.name, atom.name)


def _preserve_missing_residue_gaps(fixer, input_text: str = "") -> list[dict]:
    """Freeze sequence gaps as evidence, then prevent entire-residue synthesis."""
    chains = list(fixer.topology.chains())
    gaps = []
    for (chain_index, index), names in sorted(fixer.missingResidues.items()):
        chain = chains[chain_index]
        residues = list(chain.residues())
        gaps.append({"chain_index": chain_index, "chain_id": chain.id,
                     "residue_index": index, "residue_count": len(names), "residue_names": list(names),
                     "before_residue_id": residues[index-1].id + residues[index-1].insertionCode.strip() if index else None,
                     "after_residue_id": residues[index].id + residues[index].insertionCode.strip() if index < len(residues) else None,
                     "terminal": index in (0, len(residues))})
    internal = {(gap["chain_index"], gap["residue_index"]) for gap in gaps if not gap["terminal"]}
    # mmCIF-to-PDB inputs can have no SEQRES, and engineered constructs may
    # defeat sequence alignment. PDBFixer then reports no missing residues,
    # while OpenMM still bonds adjacent records across their absent loop.
    # Require BOTH a skipped author-number span and non-bonding observed
    # geometry. Numbering alone must not cut a legitimate renumbered bond.
    gaps.extend(_observed_backbone_discontinuities(fixer, internal))
    internal.update((gap["chain_index"], gap["residue_index"])
                    for gap in gaps if not gap["terminal"])
    fixer.missingResidues = {}
    if internal:
        _split_observed_fragments(fixer, internal, input_text)
    return gaps


def _missing_residue_gap_notes(gaps) -> list[str]:
    notes = []
    for gap in gaps:
        if gap.get("residue_count_evidence") == "absent_author_numbered_positions_only":
            notes.append(
                f"split observed backbone fragments: chain {gap['chain_id']} "
                f"{gap['before_residue_id']} -> {gap['after_residue_id']}, "
                f"C--N distance {gap['backbone_distance_angstrom']:.3f} A; "
                f"{gap['residue_count']} author-numbered position(s) absent, "
                "residue identities not established from numbering; no loop rebuilt; "
                "fragment-end hydrogen/terminal chemistry is a preparation assumption; "
                "PDBFixer terminal protonation is used and neutral caps are not implemented")
        else:
            notes.append(f"left {gap['residue_count']} wholly unobserved residue(s) unmodelled: "
                         f"chain {gap['chain_id']} between {gap['before_residue_id']} and {gap['after_residue_id']} "
                         f"({' '.join(gap['residue_names'])}); observed-coordinate gap policy; "
                         "observed fragments are separate topology chains with the original chain ID; "
                         "fragment-end hydrogen/terminal chemistry is a preparation assumption; "
                         "missing pocket geometry may affect interpretation; neutral caps are not implemented")
    return notes


def _fragment_resolved_backbone(fixer, input_text: str = "") -> list[dict]:
    """Recheck bonds after actual nonstandard residue resolution."""
    additional = _observed_backbone_discontinuities(fixer, set(), review_nonstandard=True)
    if additional:
        _split_observed_fragments(fixer,
                                  {(gap["chain_index"], gap["residue_index"]) for gap in additional},
                                  input_text)
    return additional


def _observed_backbone_discontinuities(fixer, declared_boundaries, *, review_nonstandard=False) -> list[dict]:
    """Observe discontinuities without inventing missing sequence or atoms.

    The 2 A bound is a conservative refusal of a standard peptide C--N edge,
    not a loop model. A long edge with consecutive/ambiguous numbering needs
    reviewed topology instead of an automatic cut. Explicit native CONECT
    edges are checked separately by the fragment builder.
    """
    from openmm import unit
    from . import modified
    positions = fixer.positions.value_in_unit(unit.angstrom)
    gaps = []
    for chain in fixer.topology.chains():
        residues = list(chain.residues())
        for index in range(1, len(residues)):
            if (chain.index, index) in declared_boundaries:
                continue
            left, right = residues[index-1], residues[index]
            standard = (left.name in modified.STANDARD_AMINO_ACIDS
                        and right.name in modified.STANDARD_AMINO_ACIDS)
            if not standard and not review_nonstandard:
                continue
            carbon = next((a for a in left.atoms() if a.name == "C"), None)
            nitrogen = next((a for a in right.atoms() if a.name == "N"), None)
            if carbon is None or nitrogen is None:
                continue
            distance = math.dist(positions[carbon.index], positions[nitrogen.index])
            if not math.isfinite(distance):
                raise ValueError("observed backbone distance is non-finite; reviewed topology required")
            if distance <= 2.0:
                continue
            if not standard:
                raise ValueError("long observed native modified-polymer backbone edge; reviewed topology required")
            before = int(left.id) if left.id.lstrip("-").isdigit() else None
            after = int(right.id) if right.id.lstrip("-").isdigit() else None
            if before is None or after is None or after - before <= 1:
                raise ValueError("long observed peptide edge without a numbered gap; reviewed topology required")
            gaps.append({"chain_index": chain.index, "chain_id": chain.id,
                         "residue_index": index, "residue_count": after-before-1,
                         "residue_count_evidence": "absent_author_numbered_positions_only",
                         "residue_names": [], "sequence_identity": "NOT_ESTABLISHED",
                         "before_residue_id": left.id + left.insertionCode.strip(),
                         "after_residue_id": right.id + right.insertionCode.strip(),
                         "terminal": False,
                         "detection": "observed_numbering_and_backbone_distance",
                         "backbone_distance_angstrom": distance})
    return gaps


def _split_observed_fragments(fixer, boundaries, input_text: str) -> None:
    """Prevent template C--N/O3'--P bonds across wholly unobserved gaps.

    OpenMM creates standard bonds between adjacent residue records regardless
    of numbering or distance. Separate topology chains prevent it from making
    those edges again during atom repair. Chain IDs, residue/atom IDs and all
    coordinates stay bound to their source. Other native crosslinks are retained.
    """
    from openmm import app
    original = fixer.topology
    explicit = set()
    for line in input_text.splitlines():
        if line.startswith("CONECT"):
            ids = [line[index:index+5].strip() for index in range(6, len(line), 5) if line[index:index+5].strip()]
            if ids:
                explicit.update(frozenset((ids[0], item)) for item in ids[1:])
    adjacent = set()
    for chain in original.chains():
        residues = list(chain.residues())
        adjacent.update(frozenset((residues[index-1], residues[index])) for c, index in boundaries if c == chain.index)
    rebuilt = app.Topology()
    rebuilt.setPeriodicBoxVectors(original.getPeriodicBoxVectors())
    atoms = {}
    for chain in original.chains():
        target_chain = None
        for index, residue in enumerate(chain.residues()):
            if target_chain is None or (chain.index, index) in boundaries:
                target_chain = rebuilt.addChain(chain.id)
            target_residue = rebuilt.addResidue(residue.name, target_chain, residue.id, residue.insertionCode)
            for atom in residue.atoms():
                copied = rebuilt.addAtom(atom.name, atom.element, target_residue, atom.id)
                if hasattr(atom, "formalCharge"):
                    copied.formalCharge = atom.formalCharge
                atoms[atom] = copied
    for bond in original.bonds():
        first, second = bond
        gap_edge = (frozenset((first.residue, second.residue)) in adjacent
                    and {first.name, second.name} in ({"C", "N"}, {"O3'", "P"}))
        if gap_edge:
            if frozenset((first.id, second.id)) in explicit:
                raise ValueError("explicit native backbone bond spans a declared missing-residue gap; reviewed topology required")
            continue
        rebuilt.addBond(atoms[first], atoms[second], type=bond.type, order=bond.order)
    fixer.topology = rebuilt


def _verify_observed_coordinates(fixer, expected) -> None:
    """Reject moved/missing observed heavy atoms before publishing preparation."""
    from openmm import unit
    positions = fixer.positions.value_in_unit(unit.nanometer)
    actual = {_topology_atom_key(atom): tuple(positions[atom.index]) for atom in fixer.topology.atoms()}
    for key, original in expected.items():
        current = actual.get(key)
        if current is None or any(not math.isfinite(value) or abs(value - before) > 1e-9
                                  for value, before in zip(current, original)):
            raise RuntimeError("observed heavy atom coordinates changed during preparation: " + ":".join(key))


def _heavy_atom_keys(text: str) -> set[tuple[str, str, str, str]]:
    return {(line[21], line[22:27].strip(), line[17:20].strip(), line[12:16].strip())
            for line in atom_lines(text)
            if line[76:78].strip() not in ("H", "D")}


def _terminal_oxygen_position(ca, carbon, oxygen):
    """Place new OXT in the observed CA--C--O plane at a 1.25 A C--OXT bond.

    The direction bisects the angle opposite the observed CA and O, giving
    the terminal carbon a planar three-coordinate geometry. Observed atoms
    are inputs only and are never moved.
    """
    import math

    u = tuple(a - c for a, c in zip(ca, carbon))
    v = tuple(o - c for o, c in zip(oxygen, carbon))
    nu, nv = math.sqrt(sum(x*x for x in u)), math.sqrt(sum(x*x for x in v))
    if not math.isfinite(nu + nv) or min(nu, nv) < 1e-8:
        raise ValueError("degenerate terminal carbon geometry")
    direction = tuple(-(a/nu + b/nv) for a, b in zip(u, v))
    norm = math.sqrt(sum(x*x for x in direction))
    if not math.isfinite(norm) or norm < 1e-8:
        raise ValueError("collinear terminal carbon geometry")
    return tuple(c + 0.125*d/norm for c, d in zip(carbon, direction))


def _repair_terminal_oxygens(fixer, observed: set[tuple[str, str, str, str]]) -> list[str]:
    """Correct a bad newly added OXT before hydrogens are placed."""
    import io
    import numpy as np
    from openmm import app, unit

    text = io.StringIO()
    app.PDBFile.writeFile(fixer.topology, fixer.positions, text, keepIds=True)
    bad = _unexpected_residue_bonds(text.getvalue())
    positions = np.array(fixer.positions.value_in_unit(unit.nanometer), copy=True)
    repaired = []
    for residue in fixer.topology.residues():
        atoms = {a.name: a for a in residue.atoms()}
        if not {"CA", "C", "O", "OXT"} <= atoms.keys():
            continue
        key = _topology_atom_key(atoms["OXT"])
        if key in observed or not any("OXT" in pair for pair in bad.get(key[:3], [])):
            continue
        positions[atoms["OXT"].index] = _terminal_oxygen_position(
            positions[atoms["CA"].index], positions[atoms["C"].index], positions[atoms["O"].index])
        repaired.append(":".join(key))
    if not repaired:
        return []
    fixer.positions = unit.Quantity(positions, unit.nanometer)
    return ["repositioned unresolved rebuilt OXT using the observed CA-C-O plane "
            "and a 1.25 A C-OXT bond before protonation; observed atoms unchanged: "
            + ", ".join(repaired)]


def _relax_rebuilt_atoms(fixer, observed: set[tuple[str, str, str, str]]) -> list[str]:
    """Correct false bonds in rebuilt atoms while pinning observed heavy atoms.

    PDBFixer's soft repair can leave a folded side chain (9J9D GLU284:
    CB--OE1 1.815 A). Removing that residue also removes deposited pocket
    atoms. Amber14 minimization acts only on rebuilt heavy atoms and hydrogens;
    the observed heavy atoms have zero mass and are checked after minimization.
    This is a geometry repair, not evidence of binding or affinity.
    """
    import io
    import openmm
    from openmm import app, unit

    buf = io.StringIO()
    app.PDBFile.writeFile(fixer.topology, fixer.positions, buf, keepIds=True)
    notes = []
    bad = _clashing_rebuilt_residues(buf.getvalue())
    try:
        inter = _unexpected_inter_residue_bonds(buf.getvalue())
    except ImportError as exc:
        inter = []
        notes.append(f"inter-residue connectivity assessment unavailable: {exc}")
    if not bad and not inter:
        return notes
    try:
        from .modified import repair_forcefields
        forcefields = repair_forcefields(fixer.topology)
        system = app.ForceField(*forcefields).createSystem(
            fixer.topology, nonbondedMethod=app.NoCutoff, constraints=None)
        fixed = []
        for atom in fixer.topology.atoms():
            if _topology_atom_key(atom) in observed:
                system.setParticleMass(atom.index, 0.0)
                fixed.append(atom.index)
        platform = fixer.platform or openmm.Platform.getPlatformByName(PDBFIXER_PLATFORM)
        integrator = openmm.VerletIntegrator(0.001)
        context = openmm.Context(system, integrator, platform)
        context.setPositions(fixer.positions)
        before = context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        openmm.LocalEnergyMinimizer.minimize(context, tolerance=10.0, maxIterations=200)
        check = io.StringIO()
        app.PDBFile.writeFile(fixer.topology, context.getState(getPositions=True).getPositions(), check, keepIds=True)
        remaining = _clashing_rebuilt_residues(check.getvalue())
        remaining_inter = _unexpected_inter_residue_bonds(check.getvalue())
        elements = {_topology_atom_key(a): a.element for a in fixer.topology.atoms()}
        heavy_contacts = [pair for pair in remaining_inter
                          if all(elements.get(k) != app.element.hydrogen for k in pair)]
        precision_used = bool(remaining or heavy_contacts)
        if precision_used:
            # A small rebuilt loop can retain an extra O--N contact at the
            # default stopping precision. Tighten precision rather than moving
            # or deleting the deposited backbone (8HDO THR216--rebuilt GLU217).
            openmm.LocalEnergyMinimizer.minimize(context, tolerance=1.0, maxIterations=1000)
        notes.extend(_repair_rebuilt_hydroxyl_hydrogens(fixer.topology, context))
        state = context.getState(getPositions=True)
        after = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        import numpy as np
        if not np.isfinite(after).all() or not np.array_equal(before[fixed], after[fixed]):
            raise RuntimeError("constrained repair changed observed coordinates or produced nonfinite positions")
        fixer.positions = state.getPositions()
        check = io.StringIO()
        app.PDBFile.writeFile(fixer.topology, fixer.positions, check, keepIds=True)
        remaining = _clashing_rebuilt_residues(check.getvalue())
        remaining_inter = _unexpected_inter_residue_bonds(check.getvalue())
        notes.append(
            f"{'+'.join(Path(getattr(p, 'name', p)).name for p in forcefields)} constrained repair of {len(bad)} residue(s) with unexpected "
            f"intra-residue connectivity and {len(inter)} unexpected inter-residue "
            f"bond(s); all {len(fixed)} observed heavy atoms held fixed; "
            f"{len(remaining)} residue(s) and {len(remaining_inter)} inter-residue "
            "bond(s) still unresolved. "
            f"{platform.getName()} platform, tolerance 10 kJ/mol/nm, at most 200 iterations; "
            + ("additional tolerance 1 kJ/mol/nm, at most 1000 iterations; " if precision_used else "") +
            "this does not establish a converged or biologically valid structure.")
        del context, integrator
    except Exception as exc:  # noqa: BLE001
        notes.append(f"constrained rebuilt-atom repair unavailable: {type(exc).__name__}: {exc}; "
                     "observed atoms retained; unresolved geometry must not be silently deleted")
        log.warning("constrained rebuilt-atom repair failed: %s", exc)
    return notes


def _rotated_hydroxyl_positions(positions, oxygen: int, carbon: int, hydrogen: int, angle: float):
    """Rotate one OH hydrogen about its C--O axis, retaining length and angle."""
    import numpy as np

    axis = positions[oxygen] - positions[carbon]
    norm = np.linalg.norm(axis)
    if not np.isfinite(norm) or norm < 1e-8:
        raise ValueError("degenerate hydroxyl rotation axis")
    axis = axis / norm
    vector = positions[hydrogen] - positions[oxygen]
    rotated = (vector*np.cos(angle) + np.cross(axis, vector)*np.sin(angle)
               + axis*np.dot(axis, vector)*(1-np.cos(angle)))
    result = np.array(positions, copy=True)
    result[hydrogen] = positions[oxygen] + rotated
    return result


def _repair_rebuilt_hydroxyl_hydrogens(topology, context) -> list[str]:
    """Choose a clear, lowest-Amber-energy rotamer on a fixed 24-angle OH grid.

    Amber's zero-LJ hydroxyl H can remain within Meeko's C--H covalent cutoff
    after heavy-atom repair. Sampling its C--O torsion moves one hydrogen only;
    O--H length, C--O--H angle and every other position stay unchanged. No new
    inter-residue contacts are allowed. This discrete choice is not evidence
    of energy convergence, protonation correctness or affinity.
    """
    import io
    import numpy as np
    from openmm import app, unit

    def contacts(positions):
        text = io.StringIO()
        app.PDBFile.writeFile(topology, unit.Quantity(positions, unit.nanometer), text, keepIds=True)
        return set(_unexpected_inter_residue_bonds(text.getvalue()))

    positions = context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    initial_contacts = contacts(positions)
    if not initial_contacts:
        return []
    atoms = list(topology.atoms())
    by_key = {_topology_atom_key(a): a for a in atoms}
    neighbors = {a.index: [] for a in atoms}
    for a, b in topology.bonds():
        neighbors[a.index].append(b)
        neighbors[b.index].append(a)
    hydrogens = sorted({k for pair in initial_contacts for k in pair
                        if by_key[k].element == app.element.hydrogen})
    notes = []
    for key in hydrogens:
        h = by_key[key]
        if len(neighbors[h.index]) != 1:
            continue
        oxygen = neighbors[h.index][0]
        if oxygen.element != app.element.oxygen or len(neighbors[oxygen.index]) != 2:
            continue
        carbons = [a for a in neighbors[oxygen.index] if a.element == app.element.carbon]
        if len(carbons) != 1:
            continue
        original_contacts = contacts(positions)
        best = None
        for j in range(24):
            candidate = _rotated_hydroxyl_positions(positions, oxygen.index, carbons[0].index, h.index,
                                                   2*np.pi*j/24)
            candidate_contacts = contacts(candidate)
            if any(key in pair for pair in candidate_contacts) or candidate_contacts - original_contacts:
                continue
            context.setPositions(unit.Quantity(candidate, unit.nanometer))
            energy = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
            if np.isfinite(energy) and (best is None or energy < best[0]):
                best = (energy, j, candidate)
        if best is not None:
            positions = best[2]
            notes.append(f"sampled 24 hydroxyl-H torsions for {':'.join(key)}; selected "
                         f"{best[1]*15} degrees, Amber14 energy {best[0]:.9g} kJ/mol; "
                         "no unexpected inter-residue contact involving this H; O-H length, "
                         "C-O-H angle and all other atom coordinates unchanged; "
                         "discrete geometry repair, not energy convergence")
        context.setPositions(unit.Quantity(positions, unit.nanometer))
    return notes


def _residue_number(key: str) -> Optional[int]:
    _, _, num = key.partition(":")
    num = num.strip()
    return int(num) if num.lstrip("-").isdigit() else None


def _retry_without_crosslinks(
    polymer_cls: Any, text: str, kwargs: dict, polymer: Any, ignored: list[str]
) -> Optional[tuple[Any, list[tuple[str, str]]]]:
    'Retry Meeko typing with non-peptide crosslinks cut.'
    candidates: list[tuple[str, str]] = []
    for pair in getattr(polymer, "bonds", {}) or {}:
        if not any(key in pair for key in ignored):
            continue
        n1, n2 = _residue_number(pair[0]), _residue_number(pair[1])
        if n1 is None or n2 is None or abs(n1 - n2) == 1:
            continue
        candidates.append((pair[0], pair[1]))
    if not candidates:
        return None
    try:
        retried = polymer_cls.from_pdb_string(text, bonds_to_delete=candidates, **kwargs)
    except (RuntimeError, NotImplementedError, ValueError) as exc:
        log.warning("crosslink-deletion retry failed: %s", exc)
        return None
    if len(retried.get_ignored_monomers()) >= len(ignored):
        return None
    return retried, candidates


#: Heavy-atom pairs within one residue that are bonded, so a short distance
#: between them is correct. Anything else at close range is the fault.
_BONDED_WITHIN_RESIDUE = {
    # ("N", "CD") closes the proline ring. Leaving it out made every proline
    # in the structure look like a clash -- 12 of them on 6TYM alone -- and
    # this check would then have stripped the side chain off each one, which
    # is both wrong and worse than the failure it exists to prevent. A
    # residue-geometry test has to know the residues that are rings.
    ("N", "CD"),
    ("N", "CA"), ("CA", "C"), ("C", "O"), ("C", "OXT"), ("CA", "CB"),
    ("CB", "CG"), ("CB", "CG1"), ("CB", "CG2"), ("CB", "OG"), ("CB", "OG1"),
    ("CB", "SG"), ("CG", "CD"), ("CG", "CD1"), ("CG", "CD2"), ("CG", "OD1"),
    ("CG", "OD2"), ("CG", "ND1"), ("CG", "ND2"), ("CG", "SD"), ("CG1", "CD1"),
    ("CD", "NE"), ("CD", "OE1"), ("CD", "OE2"), ("CD", "NE2"), ("CD", "CE"),
    ("SD", "CE"), ("NE", "CZ"), ("CZ", "NH1"), ("CZ", "NH2"), ("CZ", "OH"),
    ("CE", "NZ"), ("CD1", "CE1"), ("CD2", "CE2"), ("CE1", "CZ"), ("CE2", "CZ"),
    ("CD1", "NE1"), ("NE1", "CE2"), ("CD2", "CE3"), ("CE3", "CZ3"),
    ("CZ2", "CH2"), ("CZ3", "CH2"), ("CE2", "CZ2"), ("ND1", "CE1"),
    ("CE1", "NE2"), ("CD2", "NE2"),
}


def _drop_residues(text: str, keys: "set[tuple[str, str, str]]") -> str:
    """Remove every atom of the listed ``(chain, resnum, resname)`` residues."""
    keep = []
    for line in text.splitlines(keepends=True):
        if line.startswith(("ATOM", "HETATM")):
            if (line[21], line[22:27].strip(), line[17:20].strip()) in keys:
                continue
        keep.append(line)
    return "".join(keep)


def _discard_unbuildable(text: str, unresolved: "list[tuple[str, str, str]]"
                         ) -> "tuple[str, list[tuple[str, str, str]]]":
    'Delete residues no atom-level repair can make buildable.'
    if not unresolved:
        return text, []
    text = _drop_residues(text, set(unresolved))
    return text, list(unresolved)


def _clashing_rebuilt_residues(text: str) -> list[tuple[str, str, str]]:
    """Unexpected heavy-atom bonds perceived by Meeko's RDKit connectivity.

    A universal 2 A cutoff incorrectly removed deposited HIS183 of 1HFC:
    its nonbonded CD2--CE1 distance is 1.983 A and RDKit perceives no bond.
    Connectivity alone does not tell whether an atom was observed or rebuilt.
    Residue identifiers include insertion codes to keep adjacent records apart.
    """
    return list(_unexpected_residue_bonds(text))


def _unexpected_residue_bonds(text: str) -> dict[tuple[str, str, str], list[tuple[str, str]]]:
    """Perceived heavy-atom bonds outside the standard residue bond list."""
    from collections import defaultdict
    from rdkit import Chem
    from rdkit.Chem import rdDetermineBonds

    residues: dict[tuple[str, str, str], list[tuple[str, str, tuple[float, float, float]]]]
    residues = defaultdict(list)
    for line in text.splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        name = line[12:16].strip()
        element = line[76:78].strip()
        if not name or element in ("H", "D"):
            continue                      # heavy atoms only
        try:
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            continue
        residues[(line[21], line[22:27].strip(), line[17:20].strip())].append((name, element, xyz))

    bad: dict[tuple[str, str, str], list[tuple[str, str]]] = {}
    for key, atoms in residues.items():
        mol = Chem.RWMol()
        conf = Chem.Conformer(len(atoms))
        for i, (name, element, xyz) in enumerate(atoms):
            mol.AddAtom(Chem.Atom(element or name.lstrip("0123456789")[0]))
            conf.SetAtomPosition(i, xyz)
        mol = mol.GetMol()
        mol.AddConformer(conf)
        rdDetermineBonds.DetermineConnectivity(mol)
        from .modified import expected_heavy_bonds
        native_bonds = expected_heavy_bonds(key[2])
        for bond in mol.GetBonds():
            pair = (atoms[bond.GetBeginAtomIdx()][0], atoms[bond.GetEndAtomIdx()][0])
            if native_bonds is not None and tuple(sorted(pair)) in native_bonds:
                continue
            if pair not in _BONDED_WITHIN_RESIDUE and pair[::-1] not in _BONDED_WITHIN_RESIDUE:
                bad.setdefault(key, []).append(pair)
    return bad


def _unexpected_inter_residue_bonds(text: str) -> list[tuple[tuple[str, str, str, str], tuple[str, str, str, str]]]:
    """Meeko-perceived contacts absent from OpenMM's standard PDB topology.

    Include hydrogen contacts and individual atom pairs: an extra O--N contact
    can share the same residue pair as a legitimate peptide C--N bond. This is
    a geometry diagnostic, not a classification of native unusual crosslinks.
    Construct coordinate-only RDKit molecules without the private Meeko parser
    or sanitization, so malformed intra-residue valence cannot hide contacts.
    """
    import io
    from collections import defaultdict
    from meeko.polymer import find_inter_mols_bonds
    from openmm import app
    from rdkit import Chem

    topology = app.PDBFile(io.StringIO(text)).topology
    expected = {frozenset((_topology_atom_key(a), _topology_atom_key(b)))
                for a, b in topology.bonds()}
    residues = defaultdict(list)
    for line in atom_lines(text):
        name = line[12:16].strip()
        key = (line[21], line[22:27].strip(), line[17:20].strip())
        element = line[76:78].strip() or name.lstrip("0123456789")[0]
        element = "H" if element == "D" else element
        xyz = tuple(float(line[a:b]) for a, b in ((30, 38), (38, 46), (46, 54)))
        residues[key].append((name, element, xyz))
    if not residues:
        return []
    raw = {}
    for key, atoms in residues.items():
        mol = Chem.RWMol()
        conf = Chem.Conformer(len(atoms))
        for i, (_name, element, xyz) in enumerate(atoms):
            mol.AddAtom(Chem.Atom(element))
            conf.SetAtomPosition(i, xyz)
        mol = mol.GetMol()
        mol.AddConformer(conf)
        raw[key] = (mol, key[2])
    unexpected = set()
    for (r1, r2), pairs in find_inter_mols_bonds(raw).items():
        for i, j in pairs:
            a, b = (*r1, residues[r1][i][0]), (*r2, residues[r2][j][0])
            if frozenset((a, b)) not in expected:
                unexpected.add(tuple(sorted((a, b))))
    return sorted(unexpected)


def _strip_side_chains(text: str, residues: list[tuple[str, str, str]]) -> str:
    """Cut the listed residues back to their backbone.

    Deleting the residue outright would leave a gap in the chain that Meeko
    then has to template-match across; keeping the backbone and dropping the
    side chain keeps the polymer intact and loses only the atoms that were
    guessed in the first place. The receptor is rigid during docking, so a
    truncated side chain costs the pocket whatever contacts that side chain
    would have made -- which is why it is reported, not done quietly.
    """
    backbone = {"N", "CA", "C", "O", "OXT"}
    keep: list[str] = []
    targets = set(residues)
    for line in text.splitlines(keepends=True):
        if line.startswith(("ATOM", "HETATM")):
            key = (line[21], line[22:27].strip(), line[17:20].strip())
            if key in targets and line[12:16].strip() not in backbone:
                continue
        keep.append(line)
    return "".join(keep)


def _drop_atoms(text: str, victims: set[tuple[str, str, str, str]]) -> str:
    """Delete named atoms, keyed ``(chain, resnum, resname, atom)``."""
    keep: list[str] = []
    for line in text.splitlines(keepends=True):
        if line.startswith(("ATOM", "HETATM")):
            key = (line[21], line[22:27].strip(), line[17:20].strip(),
                   line[12:16].strip())
            if key in victims:
                continue
        keep.append(line)
    return "".join(keep)


def _restore_terminal_oxygen_names(pdbqt: str, source_pdb: str) -> tuple[str, list[str]]:
    """Restore O/OXT labels by exact input coordinates, leaving physics intact.

    Meeko can permute terminal carboxyl oxygen names during template matching.
    Charges, atom types, coordinates and order remain attached to their current
    atoms; only the PDB name field changes when both positions prove a swap.
    """
    def records(text):
        out = {}
        for line in atom_lines(text):
            name = line[12:16].strip()
            if name not in ("O", "OXT"):
                continue
            key = (line[21], line[22:27].strip(), line[17:20].strip())
            xyz = tuple(float(line[a:b]) for a, b in ((30, 38), (38, 46), (46, 54)))
            atoms = out.setdefault(key, {})
            if name in atoms:
                raise ValueError(f"duplicate terminal oxygen atom {key}:{name}")
            atoms[name] = xyz
        return out

    source, output = records(source_pdb), records(pdbqt)
    swaps = set()
    for key, atoms in output.items():
        expected = source.get(key, {})
        if {"O", "OXT"} <= atoms.keys() and {"O", "OXT"} <= expected.keys():
            if (atoms["O"] == expected["OXT"] and atoms["OXT"] == expected["O"]
                    and expected["O"] != expected["OXT"]):
                swaps.add(key)
    if not swaps:
        return pdbqt, []
    lines = []
    for line in pdbqt.splitlines(keepends=True):
        if line.startswith(("ATOM", "HETATM")):
            key = (line[21], line[22:27].strip(), line[17:20].strip())
            name = line[12:16].strip()
            if key in swaps and name in ("O", "OXT"):
                name = "OXT" if name == "O" else "O"
                line = line[:12] + f"{name:^4s}" + line[16:]
        lines.append(line)
    return "".join(lines), [
        f"restored coordinate-matched terminal O/OXT names in {len(swaps)} residue(s); "
        "atom coordinates, charges and types unchanged: "
        + ", ".join(":".join(k) for k in sorted(swaps))]


def _validate_observed_heavy_atoms(pdbqt: str, input_pdb: str,
                                 observed: set[tuple[str, str, str, str]]) -> None:
    """Refuse a partial or moved receptor before writing any output file.

    Standard-residue conversions can change residue names. Match the observed
    atom's chain, residue identifier and atom name against the input to Meeko;
    bind its output coordinates exactly to that retained input atom.
    """
    observed_names = {(c, n, a) for c, n, _r, a in observed}
    def positions(text):
        result = {}
        for line in atom_lines(text):
            key = (line[21], line[22:27].strip(), line[12:16].strip())
            if key not in observed_names:
                continue
            if key in result:
                raise RuntimeError(f"duplicate observed atom key during receptor typing: {key}")
            result[key] = tuple(float(line[a:b]) for a, b in ((30, 38), (38, 46), (46, 54)))
        return result
    source, typed = positions(input_pdb), positions(pdbqt)
    missing = sorted(source.keys() - typed.keys())
    moved = sorted(k for k in source.keys() & typed.keys() if source[k] != typed[k])
    if missing or moved:
        raise RuntimeError(
            f"Meeko receptor typing lost {len(missing)} or moved {len(moved)} "
            "observed heavy atom(s); PDBQT not written; "
            f"missing={missing[:10]}, moved={moved[:10]}")


def _validate_assembled_receptor_records(records: Iterable[str]) -> None:
    """Reject overlapping protein, metal or cofactor residue identities.

    A native modified residue can be omitted from protein selection, rebuilt
    as its standard parent and then appended again as a structural group.
    Each piece can preserve its own atoms while their combined receptor is
    invalid. Residue names therefore cannot distinguish atom identifiers at
    the same chain and residue position; insertion codes and chains can.
    """
    seen: set[tuple[str, str, str]] = set()
    duplicates: set[tuple[str, str, str]] = set()
    residue_names: dict[tuple[str, str], set[str]] = {}
    for line in records:
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        residue = (line[21], line[22:27].strip())
        key = (*residue, line[12:16].strip())
        if key in seen:
            duplicates.add(key)
        seen.add(key)
        residue_names.setdefault(residue, set()).add(line[17:20].strip())
    conflicts = sorted((key, sorted(names)) for key, names in residue_names.items()
                       if len(names) > 1)
    if duplicates or conflicts:
        raise RuntimeError(
            f"assembled receptor has {len(duplicates)} duplicate atom identifier(s) "
            f"and {len(conflicts)} conflicting residue identity/identities; "
            "PDBQT not written; native modified-residue chemistry requires review; "
            f"duplicates={sorted(duplicates)[:10]}, conflicts={conflicts[:10]}")


def _relieve_backbone_clashes(
    text: str, observed: Optional[set[tuple[str, str, str, str]]] = None,
) -> tuple[str, list[tuple[str, str, str]]]:
    'Cutting the side chain does nothing when both clashing atoms are backbone.'
    still = _clashing_rebuilt_residues(text)
    if not still:
        return text, []
    victims = {(c, n, r, "OXT") for c, n, r in still}
    if observed is not None:
        # PDBFixer can rename a nonstandard residue. A name conversion does
        # not turn its deposited OXT into a newly added atom.
        observed_names = {(c, n, a) for c, n, _r, a in observed}
        victims = {k for k in victims if (k[0], k[1], k[3]) not in observed_names}
    text = _drop_atoms(text, victims)
    return text, _clashing_rebuilt_residues(text)


def _typing_input_without_nonpolar_hydrogens(text: str, templates) -> tuple[str, int]:
    """Let curated templates rebuild C-bound H without changing polar H.

    Meeko determines raw residue connectivity from coordinates before matching
    its chemical templates. A misplaced carbon-bound hydrogen can acquire two
    carbon neighbours and make a valid observed heavy graph fail sanitization.
    Only atom names unanimously C-bound in all applicable curated variants are
    removed from the typing input. Polar, unknown and isotope atoms stay intact;
    the saved protonated PDB used by energy calculations is never changed.
    """
    from . import modified

    library = getattr(templates, "residue_templates", {})
    ambiguous = getattr(templates, "ambiguous", {})
    if not isinstance(library, dict) or not isinstance(ambiguous, dict):
        return text, 0
    policies = {}
    records = [line for line in text.splitlines() if line.startswith(("ATOM  ", "HETATM"))]
    present = {(line[21], line[22:27], line[17:20], line[12:16].strip()): line[76:78].strip()
               for line in records}
    for name in {line[17:20].strip() for line in records}:
        variants = ambiguous.get(name, [name])
        maps = []
        for variant in variants:
            value = library.get(variant)
            if value is None:
                maps = []
                break
            carbon_h = {}
            for atom in value.mol.GetAtoms():
                if atom.GetAtomicNum() != 1 or atom.GetIsotope():
                    continue
                parents = list(atom.GetNeighbors())
                if len(parents) == 1 and parents[0].GetAtomicNum() == 6:
                    carbon_h[value.atom_names[atom.GetIdx()]] = value.atom_names[parents[0].GetIdx()]
            maps.append(carbon_h)
        policies[name] = {hydrogen: parent for hydrogen, parent in (maps[0].items() if maps else [])
                          if all(mapping.get(hydrogen) == parent for mapping in maps)}
    removed = set()
    retained = []
    for line in text.splitlines():
        if line.startswith("CONECT"):
            continue
        if line.startswith(("ATOM  ", "HETATM")) and line[76:78].strip() == "H":
            parent = policies.get(line[17:20].strip(), {}).get(line[12:16].strip())
            if parent and present.get((line[21], line[22:27], line[17:20], parent)) == "C":
                removed.add(int(line[6:11]))
                continue
        retained.append(line)
    if not removed:
        return text, 0
    return "\n".join(retained + modified.retained_conect(text, retained)) + "\n", len(removed)


def write_receptor_pdbqt(
    clean_pdb: Path,
    metal_lines: list[str],
    out_pdbqt: Path,
    cofactor_lines: Optional[list[str]] = None,
    observed_pdb: Optional[Path] = None,
    preparation_evidence: Optional[dict] = None,
) -> list[str]:
    """Meeko-type the protein, re-attach metals and structural groups, write PDBQT.

    `cofactor_lines` are non-protein atoms that belong to the pocket rather than
    to the buffer -- a membrane lipid lining a channel's binding site, most of
    all. They are appended the way metals are, as rigid HETATM records, because
    Meeko's template engine has nothing to say about them. See
    `structure.cofactors` for what happens when they are left out.
    """
    from meeko import (
        MoleculePreparation,
        PDBQTWriterLegacy,
        Polymer,
        ResidueChemTemplates,
    )

    notes: list[str] = []
    text = clean_pdb.read_text()
    kwargs = dict(
        chem_templates=ResidueChemTemplates.create_from_defaults(),
        mk_prep=MoleculePreparation(),
        allow_bad_res=True,
    )
    # Only the selected input can establish which heavy atoms were observed.
    # Without it, conservatively treat every input heavy atom as observed.
    observed = _heavy_atom_keys((observed_pdb or clean_pdb).read_text())
    clashing = _clashing_rebuilt_residues(text)
    if clashing:
        before = _heavy_atom_keys(text)
        text, unresolved = _relieve_backbone_clashes(text, observed)
        removed_oxt = before - _heavy_atom_keys(text)
        if removed_oxt:
            notes.append("removed unresolved rebuilt terminal OXT atom(s) before "
                         "any residue deletion: " + ", ".join(":".join(k) for k in sorted(removed_oxt)))
        observed_residues = {k[:2] for k in observed}
        wholly_rebuilt = [k for k in unresolved if k[:2] not in observed_residues]
        text, discarded = _discard_unbuildable(text, wholly_rebuilt)
        if discarded:
            if preparation_evidence is not None:
                preparation_evidence["discarded_rebuilt_residues"] = [
                    {"chain": c, "residue_id": n, "residue_name": r}
                    for c, n, r in discarded]
                preparation_evidence["discarded_rebuilt_residue_count"] = len(discarded)
                preparation_evidence["observed_heavy_atoms_discarded"] = False
            where = ", ".join(f"{c}:{n} {r}" for c, n, r in discarded[:8])
            notes.append(
                f"REBUILT RESIDUES DISCARDED -- {len(discarded)} unresolved "
                f"residue(s) with no observed heavy atoms in the selected input: {where}. "
                "The corresponding unobserved region remains unavailable.")
        unresolved = _clashing_rebuilt_residues(text)
        if unresolved:
            raise RuntimeError(
                "receptor still has unexpected intra-residue connectivity in "
                + ", ".join(f"{c}:{n} {r}" for c, n, r in unresolved)
                + "; observed heavy atoms retained, PDBQT not written")

    text, removed_hydrogens = _typing_input_without_nonpolar_hydrogens(text, kwargs["chem_templates"])
    if removed_hydrogens:
        notes.append(f"rebuilt {removed_hydrogens} curated carbon-bound hydrogen(s) through Meeko templates "
                     "for receptor typing; all heavy atoms and polar/unknown input hydrogens retained; "
                     "the saved protonated PDB is unchanged")
    polymer = Polymer.from_pdb_string(text, **kwargs)
    ignored = list(polymer.get_ignored_monomers())

    if ignored:
        rescued = _retry_without_crosslinks(Polymer, text, kwargs, polymer, ignored)
        if rescued is not None:
            polymer, deleted = rescued
            notes.append(
                "excluded perceived non-peptide inter-residue bond(s) "
                f"{', '.join('-'.join(p) for p in deleted)} to recover residue(s) "
                "Meeko could not template-match; atom positions retained, but "
                "native covalent chemistry and resulting types remain unvalidated"
            )
            ignored = list(polymer.get_ignored_monomers())

    if ignored:
        notes.append(
            "meeko could not template-match "
            f"{len(ignored)} residue(s), dropped from receptor: {', '.join(ignored[:10])}"
        )
        log.warning("meeko dropped residues %s from %s", ignored, clean_pdb.name)

    written = PDBQTWriterLegacy.write_string_from_polymer(polymer)
    rigid = written[0] if isinstance(written, tuple) else written
    if not rigid or not rigid.strip():
        raise RuntimeError(f"meeko produced an empty receptor PDBQT for {clean_pdb}")

    rigid, name_notes = _restore_terminal_oxygen_names(rigid, text)
    notes.extend(name_notes)
    from .modified import restore_symmetric_names
    rigid, native_name_notes = restore_symmetric_names(rigid, text)
    notes.extend(native_name_notes)
    _validate_observed_heavy_atoms(rigid, text, observed)

    body = [ln for ln in rigid.splitlines() if ln.startswith(("ATOM", "HETATM"))]
    metal_records, metal_notes = render_metal_records(metal_lines)
    body.extend(metal_records)
    notes.extend(metal_notes)

    if cofactor_lines:
        from .cofactors import render_cofactor_records
        cofactor_records, cofactor_notes = render_cofactor_records(cofactor_lines)
        body.extend(cofactor_records)
        notes.extend(cofactor_notes)

    _validate_assembled_receptor_records(body)
    out_pdbqt.parent.mkdir(parents=True, exist_ok=True)
    out_pdbqt.write_text("\n".join(body) + "\n")
    return notes
