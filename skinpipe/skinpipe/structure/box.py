'Search-box placement: the three-tier strategy and its geometry helpers.'
from __future__ import annotations

import json
import logging
import math
import re
from pathlib import Path
from typing import Any, Optional

import numpy as np

from . import clean, mmcif, modified, rcsb

log = logging.getLogger(__name__)

#: Smallest sensible box edge. Vina needs room for the ligand plus translation;
#: a tight box around a small fragment traps the search.
MIN_BOX_EDGE: float = 18.0
#: Beyond this the search is effectively blind and exhaustiveness must go up.
WIDE_BOX_EDGE: float = 35.0
#: Below this molecular weight a HETATM is an ion or a fragment, not a ligand.
MIN_LIGAND_MW: float = 120.0
#: Below this atom count a cocrystal group cannot define a pocket on its own.
MIN_LIGAND_ATOMS: int = 6
#: Heavy-atom separation counted as a ligand-protein contact, Angstrom.
POCKET_CONTACT_A: float = 5.0
#: A ligand genuinely buried in a pocket has tens of protein atoms this close.
#: A lattice contact against a neighbouring chain has a handful, so this is a
#: deliberately permissive floor for "this ligand marks a pocket of that chain".
MIN_POCKET_CONTACTS: int = 10

#: Fewest UniProt site residues that may define a box. Below this the tier
#: falls through to the blind box rather than centring on one or two atoms:
#: three points are the fewest that can define a plane, and a single CA gave
#: a minimum-edge cube still labelled `uniprot_binding_site`.
MIN_UNIPROT_SITE_RESIDUES = 3


def validate_box_padding(value: Any) -> float:
    """Padding is a finite nonnegative distance, never a request to crop atoms."""
    if isinstance(value, bool):
        raise ValueError("box_padding must be a finite nonnegative number")
    try:
        padding = float(value)
    except (TypeError, ValueError, OverflowError):
        raise ValueError("box_padding must be a finite nonnegative number") from None
    if not math.isfinite(padding) or padding < 0:
        raise ValueError("box_padding must be a finite nonnegative number")
    return padding


def validate_explicit_pocket(value: Any) -> dict:
    """An explicit source-frame selection, never an inferred functional pocket."""
    if (not isinstance(value, dict)
            or set(value) != {"schema", "chain", "residues", "citation"}
            or value.get("schema") != "source-residue-pocket-v1"
            or not isinstance(value.get("chain"), str)
            or not re.fullmatch(r"[A-Za-z0-9]", value["chain"])
            or not isinstance(value.get("citation"), str)
            or not re.fullmatch(
                r"(?:PMCID:PMC[0-9]{1,10}|PMID:[0-9]{1,10}|DOI:10\.[0-9]{4,9}/[A-Za-z0-9._;()/:-]{1,120})",
                value["citation"])):
        raise ValueError("Explicit pocket requires a source chain and literature identifier")
    residues = value.get("residues")
    if not isinstance(residues, list) or not 3 <= len(residues) <= 32:
        raise ValueError("Explicit pocket requires 3..32 source residues")
    normalized, seen = [], set()
    for residue in residues:
        if (not isinstance(residue, dict)
                or set(residue) - {"number", "amino_acid", "insertion_code"}
                or type(residue.get("number")) is not int
                or not 1 <= residue["number"] <= 9999
                or not isinstance(residue.get("amino_acid"), str)
                or residue["amino_acid"] not in set("ACDEFGHIKLMNPQRSTVWYOU")
                or not isinstance(residue.get("insertion_code", ""), str)
                or not re.fullmatch(r"[A-Za-z]?", residue.get("insertion_code", ""))):
            raise ValueError("Invalid explicit source residue")
        key = (residue["number"], residue.get("insertion_code", ""))
        if key in seen:
            raise ValueError("Explicit source residues must be unique")
        seen.add(key)
        normalized.append({"number": key[0], "insertion_code": key[1],
                           "amino_acid": residue["amino_acid"]})
    return {**value, "residues": sorted(normalized,
                                       key=lambda r: (r["number"], r["insertion_code"]))}


def explicit_residue_box(raw_text, definition, chain, padding, *, retained_serials=None):
    """Require every selected residue and its identity in the retained source frame."""
    definition = validate_explicit_pocket(definition)
    if definition["chain"] != chain:
        raise ValueError("Explicit pocket chain differs from the verified target chain")
    wanted = {(r["number"], r["insertion_code"]): r["amino_acid"]
              for r in definition["residues"]}
    ca_counts = {key: 0 for key in wanted}
    points = []
    for line in clean.atom_lines(raw_text):
        if line[21] != chain or line[76:78].strip() in {"H", "D"}:
            continue
        if retained_serials is not None and int(line[6:11]) not in retained_serials:
            continue
        key = (int(line[22:26]), line[26:27].strip())
        if key not in wanted:
            continue
        if _THREE_TO_ONE.get(line[17:20].strip()) != wanted[key]:
            raise ValueError("Explicit pocket source amino acid differs")
        if line[12:16].strip() == "CA":
            ca_counts[key] += 1
        point = [float(line[start:start + 8]) for start in (30, 38, 46)]
        if not all(math.isfinite(x) for x in point):
            raise ValueError("Explicit pocket has nonfinite source coordinates")
        points.append(point)
    if any(count != 1 for count in ca_counts.values()):
        raise ValueError("Explicit pocket residue is absent, discarded or ambiguous")
    center, size = _box_from_points(np.asarray(points), validate_box_padding(padding))
    return center, size
















COFACTOR_HINTS: frozenset[str] = frozenset({
    "ADP", "ATP", "AMP", "ANP", "ACP", "AGS", "GTP", "GDP", "GMP", "GNP",
    "UTP", "UDP", "UMP", "CTP", "CDP", "CMP", "TTP", "TDP",
    "NAD", "NAI", "NAP", "NDP", "NAJ", "FAD", "FMN", "FDA",
    "SAM", "SAH", "COA", "ACO", "HEM", "HEC", "HEB", "PLP", "TPP", "TDP",
    "BTN", "MTA", "PPS", "PAP", "P5P",
})







BOX_RUNNER_UP_RATIO = 2.0


def box_ligand_ambiguity(ligands: list[dict[str, Any]], chain: Optional[str],
                         chosen: Optional[dict[str, Any]]
                         ) -> "tuple[Optional[str], Optional[str]]":
    '``(severity, message)`` for the choice of box ligand.'
    if chosen is None:
        return None, None
    pool = [lig for lig in ligands if lig.get("drug_like")
            and (chain is None or lig["chain"] == chain)
            and lig is not chosen]
    if not pool:
        return None, None
    best_other = max(pool, key=lambda d: d.get("contacts_per_atom", 0.0))
    mine = chosen.get("contacts_per_atom", 0.0)
    theirs = best_other.get("contacts_per_atom", 0.0)
    if chosen["resname"].upper() in COFACTOR_HINTS:
        others = ", ".join(sorted({lig["resname"] for lig in pool}))
        return "unverified", (
                f"the box was placed on {chosen['resname']}, a cofactor or "
                f"nucleotide, because it is the most enclosed group in this "
                f"chain ({mine:.2f} contacts per heavy atom). {others} "
                f"also present. For an ATP-competitive programme that is the "
                f"right site; for an allosteric one it is not, and nothing in "
                f"the structure can tell those apart. Pass ligand_code to pin "
                f"the intended ligand. The redock control will not catch this: "
                f"it reproduces the pose of whichever ligand was chosen.")
    if theirs and mine < theirs * BOX_RUNNER_UP_RATIO:
        return "note", (
                f"the box was placed on {chosen['resname']} "
                f"({mine:.2f} contacts per heavy atom) over "
                f"{best_other['resname']} ({theirs:.2f}); they are close enough "
                f"that the choice rests on the enclosure measure alone. Pass "
                f"ligand_code to pin one.")
    return None, None


def pick_ligand(
    ligands: list[dict[str, Any]],
    chain: Optional[str],
    ligand_code: Optional[str],
) -> Optional[dict[str, Any]]:
    candidates = [lig for lig in ligands if lig["drug_like"]]
    if ligand_code:
        wanted = ligand_code.strip().upper()






        candidates = [lig for lig in ligands
                      if lig["resname"] == wanted
                      or str(lig.get("component_id") or "").upper() == wanted]
        if not candidates:
            return None
    if chain is not None:
        in_chain = [lig for lig in candidates if lig["chain"] == chain]
        if in_chain:
            candidates = in_chain
    if not candidates:
        return None








    return max(candidates,
               key=lambda d: (d.get("contacts_per_atom", 0.0),
                              not d.get("lipid_like", False),
                              d["n_atoms"]))


def ligand_touches_chain(raw_text: str, lig: dict[str, Any], chain: str,
                        cutoff: float = POCKET_CONTACT_A,
                        min_contacts: int = MIN_POCKET_CONTACTS) -> int:
    """How many atoms of ``chain`` lie within ``cutoff`` of this ligand copy.

    Used to decide whether a ligand annotated to another chain nevertheless
    marks a pocket of ``chain``. Returns the contact count; compare it against
    ``min_contacts`` rather than testing for any contact at all, so a glancing
    crystal-lattice touch is not mistaken for a binding site.
    """
    lig_pts = _ligand_coords(raw_text, lig)
    if lig_pts.size == 0:
        return 0
    chain_pts = np.asarray(
        [
            [float(pdb_line[30:38]), float(pdb_line[38:46]), float(pdb_line[46:54])]
            for pdb_line in clean.atom_lines(raw_text)
            if pdb_line.startswith("ATOM  ") and pdb_line[21] == chain and pdb_line[16] in (" ", "A")
        ],
        dtype=float,
    )
    if chain_pts.size == 0:
        return 0
    d2 = ((chain_pts[:, None, :] - lig_pts[None, :, :]) ** 2).sum(axis=2)
    return int((d2.min(axis=1) <= cutoff * cutoff).sum())


def _box_from_points(points: np.ndarray, padding: float) -> tuple[tuple, tuple]:
    'Centre and size of a box holding `points` with `padding` on every side.'
    padding = validate_box_padding(padding)
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1:] != (3,) or not len(points) or not np.isfinite(points).all():
        raise ValueError("search box requires non-empty finite 3D coordinates")
    lo, hi = points.min(axis=0), points.max(axis=0)
    center = tuple(float(v) for v in (lo + hi) / 2.0)
    extent = (hi - lo) + 2.0 * float(padding)
    size = tuple(float(max(v, MIN_BOX_EDGE)) for v in extent)
    return center, size


def _chain_ca_by_resnum(pdb_text: str, chain: str) -> dict[tuple[int, str], np.ndarray]:
    out: dict[tuple[int, str], np.ndarray] = {}
    for line in clean.atom_lines(pdb_text):
        if line[21] != chain or line[17:20].strip().upper() not in _THREE_TO_ONE:
            continue
        if line[12:16].strip() != "CA":
            continue
        try:
            resnum = int(line[22:26])
        except ValueError:
            continue
        out[(resnum, line[26:27] or " ")] = np.array(
            [float(line[30:38]), float(line[38:46]), float(line[46:54])], dtype=float
        )
    return out


_THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M", "SEC": "U", "PYL": "O",
}


def _chain_sequence(pdb_text: str, chain: str) -> tuple[str, list[tuple[int, str]]]:
    """Modelled sequence and author residue identities, including insertions."""
    seq: list[str] = []
    nums: list[tuple[int, str]] = []
    seen: set[tuple[int, str]] = set()
    for line in clean.atom_lines(pdb_text):
        if (line[21] != chain or line[12:16].strip() != "CA"
                or line[17:20].strip().upper() not in _THREE_TO_ONE):
            continue
        try:
            resnum = int(line[22:26])
        except ValueError:
            continue
        identity = (resnum, line[26:27] or " ")
        if identity in seen:
            continue
        seen.add(identity)
        seq.append(_THREE_TO_ONE.get(line[17:20].strip().upper(), "X"))
        nums.append(identity)
    return "".join(seq), nums


def _map_uniprot_positions(
    uniprot_seq: str, positions: list[int], pdb_text: str, chain: str
) -> tuple[list[np.ndarray], int]:
    """Map 1-based UniProt positions onto CA coordinates via global alignment."""
    from Bio import Align

    struct_seq, struct_nums = _chain_sequence(pdb_text, chain)
    if not struct_seq or not uniprot_seq or not positions:
        return [], len(positions)

    aligner = Align.PairwiseAligner(
        mode="global", open_gap_score=-11.0, extend_gap_score=-1.0, substitution_matrix=None
    )
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    alignment = aligner.align(uniprot_seq, struct_seq)[0]
    # aligned blocks: [[(u_start, u_end), ...], [(s_start, s_end), ...]]
    uni_blocks, str_blocks = alignment.aligned
    uni_to_struct: dict[int, tuple[int, str]] = {}
    for (u0, u1), (s0, s1) in zip(uni_blocks, str_blocks):
        for offset in range(u1 - u0):
            uni_to_struct[u0 + offset + 1] = struct_nums[s0 + offset]

    ca = _chain_ca_by_resnum(pdb_text, chain)
    coords = [ca[uni_to_struct[p]] for p in positions
              if p in uni_to_struct and uni_to_struct[p] in ca]
    # Both numbers. A site residue in an unmodelled region, or one the alignment
    # places on a gap, is dropped here without a word, and the caller used to
    # report only the survivors -- "box from 4 site residues" when seven were
    # annotated reads as a complete mapping.
    return coords, len(positions)



def _comp_id_map(code: str, out_dir: Path,
                 pdb_path: Optional[Path] = None) -> dict[str, Optional[str]]:
    'The PDB-truncation map `mmcif.cif_to_pdb` wrote for this entry, if any.'
    # Named locations only. An `rglob` two levels up crawls whatever else is
    # under a shared parent and can return a sidecar belonging to a different
    # tree -- caught by two test files that both used `Y.comp_ids.json` and read
    # each other's. The converter writes it beside the PDB it produced, so those
    # are the places to look.
    name = f"{code}.comp_ids.json"
    # `pdb_path` is authoritative: `mmcif.cif_to_pdb` writes the sidecar beside
    # the file it produced, so the only place that always holds is next to it.
    # Guessing from `out_dir` was wrong in the one layout that matters --
    # `prepare_receptor` passes `workdir/receptors/<code>_<chain>` while the
    # sidecar is at `workdir/pdb/`, two levels up and across, so the lookup
    # returned an empty map on every production run and the truncated id was used
    # unresolved. The `out_dir` guesses stay as a fallback for callers that have
    # no path to give.
    candidates = []
    if pdb_path is not None:
        candidates.append(Path(pdb_path).with_suffix(".comp_ids.json"))
        candidates.append(Path(pdb_path).parent / name)
    candidates += [out_dir / name, out_dir / "pdb" / name,
                   out_dir.parent / name, out_dir.parent / "pdb" / name,
                   out_dir.parent.parent / name,
                   out_dir.parent.parent / "pdb" / name]
    seen: set[Path] = set()
    for h in candidates:
        if h in seen:
            continue
        seen.add(h)
        try:
            data = json.loads(h.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        if not isinstance(data, dict):
            continue
        nested = "resolved" in data or "ambiguous" in data
        stamped = (data.get("schema") == mmcif._SIDECAR_SCHEMA
                   and isinstance(data.get("source"), dict))
        if nested and stamped:
            out: dict[str, Optional[str]] = {
                str(k): str(v) for k, v in (data.get("resolved") or {}).items()}
            # An ambiguous prefix maps to None: the reader must be able to tell
            # "no entry" from "more than one candidate", because the second has
            # to refuse rather than fall back to the truncation.
            for k in (data.get("ambiguous") or {}):
                out[str(k)] = None
            return out
        # Everything else is a map this converter did not write, and the trust
        # rule is the same whatever shape it has: no stamp, no answer.
        #
        # Two shapes reach here, and the distinction between them is exactly the
        # trap. A *flat* prefix->id file came from the converter that could not
        # express a collision at all -- one prefix, one winner, the loser
        # silently gone. A *nested* file with `resolved`/`ambiguous` but no stamp
        # looks like it came from a converter that knew about collisions, and it
        # did -- but only collisions between two ids longer than three
        # characters. It still recorded `A1C -> A1C8J` as the unique answer for
        # an entry that also holds a genuine three-character `A1C`, because the
        # short id was never grouped. Accepting the nested shape on the strength
        # of its shape was a hole of exactly the size of the defect above it.
        #
        # `prepare.fetch_pdb` re-derives the map from the mmCIF whenever it is
        # still in the cache, so this branch is reached only when the mmCIF has
        # been deleted and there is nothing left to check the map against. Every
        # key it names is returned as ambiguous, which makes the template fetch
        # refuse and say so on the report, rather than validating a pose against
        # a molecule chosen by a converter whose blind spot cannot be ruled out.
        keys = list((data.get("resolved") or {})) + list((data.get("ambiguous") or {})) \
            if nested else list(data)
        log.warning("%s holds a component-id map with no converter stamp, so "
                    "whether it distinguishes components that share a truncated "
                    "prefix cannot be established (%d entr%s). Treating them as "
                    "ambiguous. Delete the cached PDB and mmCIF for this entry to "
                    "have it re-derived.",
                    h.name, len(keys), "y" if len(keys) == 1 else "ies")
        return {str(k): None for k in keys}
    return {}









MIN_OWN_POCKET_SHARE = 0.5


def pocket_chain_share(raw_text: str, lig: dict[str, Any], chain: str,
                       cutoff: float = POCKET_CONTACT_A, *,
                       retained_chains: Optional[set[str]] = None,
                       retained_serials: Optional[set[int]] = None,
                       ) -> tuple[float, dict[str, int]]:
    """Fraction of the observed polymer pocket wall supplied by retained chains.

    ``retained_serials`` measures the atoms actually selected, including any
    omissions within a retained chain. Without either selection argument this
    measures the primary chain alone. Native
    modified polymer atoms count; loose ligands and cofactors do not. Empty
    counts mean unassessed geometry, including when ligand coordinates are
    absent. That case returns zero and must not be interpreted as a split site.
    """
    lig_pts = _ligand_coords(raw_text, lig)
    if lig_pts.size == 0:
        return 0.0, {}
    counts: dict[str, int] = {}
    retained_count = 0
    native = modified.groups(raw_text, known_phospho=True)
    for line in clean.atom_lines(raw_text):
        if line[16] not in (" ", "A") or line[76:78].strip() in ("H", "D"):
            continue
        key = (line[21], line[22:27].strip(), line[17:20].strip())
        if not line.startswith("ATOM") and key not in native:
            continue
        try:
            pt = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        except ValueError:
            continue
        if np.linalg.norm(lig_pts - pt, axis=1).min() <= cutoff:
            counts[line[21]] = counts.get(line[21], 0) + 1
            if retained_serials is not None and int(line[6:11]) in retained_serials:
                retained_count += 1
    total = sum(counts.values())
    if not total:
        return 0.0, counts
    if retained_serials is not None:
        return retained_count / total, counts
    retained = {chain} if retained_chains is None else retained_chains
    return sum(n for cid, n in counts.items() if cid in retained) / total, counts


def place_box(
    raw_text: str,
    ligands: list[dict[str, Any]],
    chain: str,
    ligand_code: Optional[str],
    box_padding: float,
    entity: Optional[dict[str, Any]],
    out_dir: Path,
    code: str,
    pdb_path: Optional[Path] = None,
    *,
    retained_chains: Optional[set[str]] = None,
    retained_serials: Optional[set[int]] = None,
    explicit_pocket: Optional[dict] = None,
) -> tuple[tuple, tuple, str, list[str], Optional[dict[str, Any]]]:
    """Run the three-tier box strategy. Always returns at least one note."""
    box_padding = validate_box_padding(box_padding)
    if explicit_pocket is not None:
        definition = validate_explicit_pocket(explicit_pocket)
        center, size = explicit_residue_box(raw_text, definition, chain, box_padding,
                                          retained_serials=retained_serials)
        return center, size, "explicit_source_residues", [
            f"box from {len(definition['residues'])} explicitly selected source residues "
            f"in target chain {chain}; citation {definition['citation']}; "
            "source amino-acid identities and retained coordinates checked; "
            "UniProt numbering, pocket accessibility and functional binding require separate assessment"
        ], None
    from .ligand import extract_cocrystal_mol

    notes: list[str] = []

    chosen = pick_ligand(ligands, chain, ligand_code)
    if ligand_code and chosen is None:
        notes.append(
            f"BOX LIGAND MISSING -- {ligand_code} was pinned and no cocrystal "
            f"group in this entry carries that code, under either its "
            f"three-character resName or its full component id. The box falls "
            f"through to the next strategy, which is not the site that was "
            f"asked for: a pin exists because burial picked the wrong group, "
            f"so silently reverting to burial or to the chain centroid answers "
            f"a different question. Check the code against the entry's "
            f"components.")

    # What was taken out of the running, and why. A group excluded in silence
    # leaves a reader asking why the box is where it is, with nothing to answer
    # from; this is the sentence that answers it -- and the one that surfaces a
    # peptidomimetic inhibitor wrongly read as a chain residue, the boundary
    # nobody has an entry to test.
    excluded = [lig for lig in ligands
                if lig.get("polymer_residue")
                and (chain is None or lig["chain"] == chain)]
    if excluded:
        where = ", ".join(f"{lig['resname']} {lig['chain']}{lig['resnum']} "
                          f"({lig['n_atoms']} atoms)" for lig in excluded[:6])
        notes.append(
            f"{len(excluded)} HETATM group(s) carry the backbone triad N/CA/C "
            f"and were read as modified chain residues, not ligands, so they "
            f"were not box candidates: {where}. A peptidomimetic ligand would "
            "look the same; if one of these is the compound you meant to box "
            "on, pass its code as ligand_code.")

    # Falling back across chains ranks by contact with the chain being docked,
    # not by the burial that won in the ligand's own chain. `pick_ligand` has no
    # coordinates and cannot do this; here the raw text is in hand. Without it
    # the fallback picks whichever cross-chain ligand is most enclosed *where it
    # sits*, which says nothing about the pocket being prepared, and the check
    # below then only validates that one guess.
    if chosen is not None and chain is not None and chosen["chain"] != chain \
            and not ligand_code:
        others = [lig for lig in ligands
                  if lig.get("drug_like") and lig["chain"] != chain]
        if len(others) > 1:
            ranked = sorted(
                others,
                key=lambda lig: (ligand_touches_chain(raw_text, lig, chain),
                                 lig.get("contacts_per_atom", 0.0)),
                reverse=True)
            if ranked and ranked[0] is not chosen:
                notes.append(
                    f"no cocrystal ligand in chain {chain}; of the "
                    f"{len(others)} in other chains, {ranked[0]['resname']} "
                    f"{ranked[0]['chain']}{ranked[0]['resnum']} has the most "
                    f"contact with it, not {chosen['resname']} "
                    f"{chosen['chain']}{chosen['resnum']}, which is the most "
                    "enclosed where it sits")
                chosen = ranked[0]
    if chosen is not None and chosen["chain"] != chain:
        contacts = ligand_touches_chain(raw_text, chosen, chain)
        if contacts >= MIN_POCKET_CONTACTS:
            notes.append(
                f"cocrystal {chosen['resname']} is annotated to chain {chosen['chain']} "
                f"but has {contacts} atom contacts with chain {chain}, so it sits at the "
                "interface and does define a pocket of the docked chain")
        else:
            notes.append(
                f"IGNORED cocrystal {chosen['resname']} in chain {chosen['chain']}: chain "
                f"{chain} is the one being docked and only {contacts} of its atoms come "
                f"within {POCKET_CONTACT_A} A of that ligand. Boxing on it would put the "
                "search volume outside the receptor. Falling through to the next tier")
            chosen = None
    if chosen is not None:







        if not ligand_code:
            # Only when nobody pinned one. A caller who named the ligand has
            # already answered the question this warning asks.
            severity, ambiguity = box_ligand_ambiguity(ligands, chain, chosen)
            if severity == "unverified":
                notes.append(f"BOX LIGAND AMBIGUOUS -- {ambiguity}")
            elif severity == "note":
                notes.append(f"box ligand alternatives -- {ambiguity}")
        # Touching the docked chain is not the same as being *in* its pocket.
        share, per_chain = pocket_chain_share(raw_text, chosen, chain,
                                            retained_chains=retained_chains,
                                            retained_serials=retained_serials)
        retained = {chain} if retained_chains is None else retained_chains
        if not per_chain:
            notes.append(
                f"BOX POCKET UNASSESSED -- {chosen['resname']} "
                f"{chosen['chain']}{chosen['resnum']} has no observed polymer "
                f"heavy-atom contacts within {POCKET_CONTACT_A:g} A. Pocket "
                "coverage cannot be inferred from an empty contact set.")
        elif share < MIN_OWN_POCKET_SHARE:
            where = ", ".join(f"chain {c} {n}" for c, n in
                              sorted(per_chain.items(), key=lambda kv: -kv[1]))
            notes.append(
                f"BOX POCKET SPLIT -- {chosen['resname']} "
                f"{chosen['chain']}{chosen['resnum']} is lined by {where} within "
                f"{POCKET_CONTACT_A:g} A; retained chains "
                f"{','.join(sorted(retained)) or '(none)'} supply only "
                f"{share:.0%} of the observed pocket wall. Most of the "
                "observed site is absent from the receptor, so a score "
                "against this file does not validate binding at the full site.")
        elif retained_chains is not None:
            notes.append(
                f"observed pocket coverage -- retained chains "
                f"{','.join(sorted(retained)) or '(none)'} supply {share:.0%} "
                f"of the polymer heavy-atom contacts within {POCKET_CONTACT_A:g} A; "
                "this measures atom coverage, not docking or functional validation")
        coords = _ligand_coords(raw_text, chosen)
        center, size = _box_from_points(coords, box_padding)
        notes.append(
            f"box from cocrystal {chosen['resname']} "
            f"({chosen['n_atoms']} atoms, chain {chosen['chain']}) + {box_padding} A padding"
        )
        sdf_path = out_dir / f"{code}_{chosen['resname']}_crystal.sdf"
        mol = extract_cocrystal_mol(raw_text, chosen, sdf_path,
                                    cache_dir=out_dir.parent,
                                    comp_id_map=_comp_id_map(code, out_dir, pdb_path))
        if mol is None:
            # Two different failures, and only one of them is about this
            # structure. `extract_cocrystal_mol` records `template_fetch_error`
            # when RCSB did not serve the chemical-component template; without
            # that distinction both arrive as "bond orders could not be
            # perceived", which reads as a property of the entry and invites
            # giving up on a ligand that is sitting in the PDB file.
            fetch_error = chosen.get("template_fetch_error")
            if fetch_error:
                notes.append(
                    f"the chemical-component template for {chosen['resname']} "
                    f"could not be fetched from RCSB ({fetch_error[:120]}), so "
                    "bond orders are unknown and redock validation is "
                    "unavailable **for this run**. The ligand is present in "
                    f"{code} with {chosen['n_atoms']} atoms; re-run when RCSB "
                    "answers and the control will run."
                )
            else:
                notes.append(
                    f"could not perceive bond orders for {chosen['resname']}; "
                    "redock validation unavailable"
                )
            sdf_path = None
        record = dict(chosen)
        record["sdf_path"] = sdf_path
        return center, size, "cocrystal", notes, record

    uniprot = (entity or {}).get("uniprot")
    if uniprot:
        lookup_failed = False
        try:
            ann = rcsb.uniprot_site_annotations(uniprot)
            coords, n_asked = _map_uniprot_positions(
                ann["sequence"], ann["positions"], raw_text, chain)
        except rcsb.RCSBError as exc:
            lookup_failed = True
            coords, n_asked = [], 0
            notes.append(f"UniProt lookup for {uniprot} failed: {exc}")
        if len(coords) >= MIN_UNIPROT_SITE_RESIDUES:
            center, size = _box_from_points(np.asarray(coords), box_padding + 3.0)
            notes.append(
                f"no cocrystal ligand; box from {len(coords)} of {n_asked} UniProt "
                f"{uniprot} site residue(s) mapped by sequence alignment"
                + ("; the rest fall in unmodelled regions or on alignment gaps"
                   if len(coords) < n_asked else "")
            )
            return center, size, "uniprot_binding_site", notes, None
        if coords:
            # Below the floor this is not a site, it is a point. One surviving
            # residue gave an 18 A cube centred on a single CA, still labelled
            # `uniprot_binding_site` -- the second-highest-confidence label the
            # module has, and not wide enough to trip the WIDE_BOX_EDGE note.
            # Three is the fewest that can define a plane.
            log.warning("UniProt %s site mapping for chain %s kept only %d of %d "
                        "residues; too few to define a pocket, falling through",
                        uniprot, chain, len(coords), n_asked)
            notes.append(
                f"UniProt {uniprot} site mapping kept only {len(coords)} of "
                f"{n_asked} residue(s), below the {MIN_UNIPROT_SITE_RESIDUES} "
                "needed to define a pocket; fell through to a blind box")
        elif not lookup_failed:
            notes.append(
                f"UniProt {uniprot} has no mappable site annotation for chain {chain}")
    else:
        notes.append("no UniProt accession available for this chain")

    coords = np.asarray(
        [
            [float(pdb_line[30:38]), float(pdb_line[38:46]), float(pdb_line[46:54])]
            for pdb_line in clean.atom_lines(raw_text)
            if pdb_line.startswith("ATOM  ") and pdb_line[21] == chain
        ]
    )
    if not len(coords) or not np.isfinite(coords).all():
        raise ValueError(f"chain {chain!r} has no finite observed polymer coordinates for a blind box")
    center = tuple(float(v) for v in coords.mean(axis=0))
    notes.append(
        "WARNING: blind box at whole-protein centroid - no cocrystal ligand and no "
        "UniProt site annotation. Scores from this receptor are not interpretable "
        "without an independent pocket definition."
    )
    return center, (30.0, 30.0, 30.0), "centroid", notes, None


def _ligand_coords(raw_text: str, lig: dict[str, Any]) -> np.ndarray:
    resnum = str(lig["resnum"])
    # Insertion code, column 27. `list_cocrystal_ligands` groups copies by
    # (resname, chain, resnum, icode) and reports them separately, so leaving it
    # out here merges every copy that shares a residue number into one molecule.
    # Reproduced on a synthetic entry with LIG A 301A and LIG A 301B, three atoms
    # each: the listing sees two 3-atom copies and the extraction returns six
    # atoms for either, spanning 102 A and yielding a 112 A box. A box that size
    # is not a pocket, and the redock reference would be two copies at once.
    icode = str(lig.get("icode") or "").strip()
    pts = [
        [float(pdb_line[30:38]), float(pdb_line[38:46]), float(pdb_line[46:54])]
        for pdb_line in clean.atom_lines(raw_text)
        if pdb_line.startswith("HETATM")
        and pdb_line[17:20].strip().upper() == lig["resname"]
        and pdb_line[21] == lig["chain"]
        and pdb_line[22:26].strip() == resnum
        and pdb_line[26].strip() == icode
        and pdb_line[16] in (" ", "A")
    ]
    return np.asarray(pts, dtype=float)
