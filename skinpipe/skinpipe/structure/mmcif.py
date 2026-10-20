"""mmCIF to PDB conversion for structures RCSB does not serve as PDB.

Why not Biopython's ``PDBIO``
-----------------------------
An entry has no PDB-format file for exactly two reasons: more than 99 999 atoms,
or more than 62 chains. ``PDBIO`` refuses to write either
(``Atom serial number ('100000') exceeds PDB format limit``), so using it as the
mmCIF fallback fails on precisely the structures that need the fallback. Verified
on 9MLH (100 k+ atoms) and 4V6X (80S ribosome, ~180 chains).

This writer emits ``_atom_site`` directly and wraps the atom serial, which is
safe because nothing downstream in this pipeline reads columns 7-11: chain,
residue and ligand selection all key on the chain ID, residue name and residue
number. The chain-count limit is genuinely unrepresentable, so that one raises
with a message naming the cached mmCIF rather than quietly dropping chains.

Only the first model is written. Author numbering (``auth_asym_id`` /
``auth_seq_id``) is preserved, because that is what the literature, UniProt
cross-references and every ``ligand_code`` a user will type refer to.
"""
from __future__ import annotations

import hashlib
import json
import logging
from pathlib import Path
from typing import Any, Optional

log = logging.getLogger(__name__)

_PDB_CHAIN_ALPHABET = (
    "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
)


class StructureTooLargeError(RuntimeError):
    """The structure cannot be represented in PDB format at all."""


def _column(data: dict[str, Any], *names: str) -> list[str]:
    for name in names:
        if name in data:
            value = data[name]
            return list(value) if isinstance(value, list) else [value]
    return []


def _format_atom_name(name: str, element: str) -> str:
    """PDB puts one- and two-letter element names in column 14, not 13."""
    name = name.strip().strip('"')
    if len(name) >= 4:
        return name[:4]
    if len(element.strip()) == 2 and name.upper().startswith(element.strip().upper()):
        return f"{name:<4}"
    return f" {name:<3}"


_SIDECAR_SCHEMA = 2
SIDECAR_SUFFIX = ".comp_ids.json"


def sidecar_path(pdb_path: Path) -> Path:
    """Where the truncation map for ``pdb_path`` lives."""
    return Path(pdb_path).with_suffix(SIDECAR_SUFFIX)


def _source_stamp(cif_path: Path) -> dict[str, Any]:
    """Identify the mmCIF a sidecar was derived from.

    The map is only as trustworthy as the file it came from, and a converted PDB
    outlives the code that wrote it. Content hash rather than mtime: a cache
    copied between machines keeps its contents and loses its timestamps.
    """
    h = hashlib.sha256()
    with open(cif_path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return {"name": Path(cif_path).name, "sha256": h.hexdigest(),
            "size": Path(cif_path).stat().st_size}


def comp_id_translation(comp_ids: list[str]) -> tuple[dict[str, str], dict[str, list[str]]]:
    """Split the entry's component ids into ``(resolved, ambiguous)`` by prefix.

    Every distinct id is grouped, not only the ones longer than three
    characters. Grouping the long ones alone misses the collision that matters
    most: a genuine three-character component and a five-character one that
    share a prefix -- `A1C` (a chlorine compound) alongside `A1C8J`
    (Asivatrep). Both are written as `A1C`, so the PDB file cannot tell them
    apart, yet the long id was the only member of its group and was recorded as
    the unique answer for the prefix. Both ligand copies then resolved to
    `A1C8J`, and the real `A1C` was looked up as a molecule it is not.

    A prefix is resolved only when exactly one distinct component maps to it
    *and* the translation is not the identity, so an ordinary entry of
    three-character ids still produces an empty map.
    """
    by_prefix: dict[str, list[str]] = {}
    for c in dict.fromkeys(comp_ids):
        by_prefix.setdefault(c[:3], []).append(c)
    resolved = {k: v[0] for k, v in by_prefix.items() if len(v) == 1 and v[0] != k}
    ambiguous = {k: sorted(v) for k, v in by_prefix.items() if len(v) > 1}
    return resolved, ambiguous


def write_comp_id_sidecar(cif_path: Path, pdb_path: Path,
                          comp_ids: Optional[list[str]] = None) -> Path:
    """Derive the truncation map from ``cif_path`` and write it beside ``pdb_path``.

    Separate from :func:`cif_to_pdb` because a cached PDB has to be able to get
    a current map without its coordinates being rewritten -- see
    :func:`skinpipe.structure.prepare.fetch_pdb`.
    """
    cif_path = Path(cif_path)
    if comp_ids is None:
        from Bio.PDB.MMCIF2Dict import MMCIF2Dict
        data = MMCIF2Dict(str(cif_path))
        comp_ids = _column(data, "_atom_site.auth_comp_id",
                           "_atom_site.label_comp_id")
    resolved, ambiguous = comp_id_translation(comp_ids)
    out = sidecar_path(pdb_path)
    if resolved:
        log.warning("%s: %d component id(s) do not fit the PDB resName field and "
                    "were truncated (%s); the full ids are in %s",
                    cif_path.name, len(resolved),
                    ", ".join(f"{k}<-{v}" for k, v in list(resolved.items())[:4]),
                    out.name)
    if ambiguous:
        log.error("%s: %d truncated prefix(es) are shared by more than one "
                  "component (%s). The PDB file cannot tell those ligands apart, "
                  "so a chemical-component template cannot be chosen for them and "
                  "the redock control is unavailable. Work from the mmCIF.",
                  cif_path.name, len(ambiguous),
                  "; ".join(f"{k} <- {'/'.join(v)}" for k, v in ambiguous.items()))
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({"schema": _SIDECAR_SCHEMA,
                               "source": _source_stamp(cif_path),
                               "resolved": resolved,
                               "ambiguous": ambiguous},
                              indent=2, sort_keys=True), encoding="utf-8")
    return out


def sidecar_is_current(cif_path: Path, pdb_path: Path) -> bool:
    """Was the map beside ``pdb_path`` written by this converter, from this mmCIF?

    An unstamped or differently-stamped sidecar is not evidence of anything: the
    converter that wrote it may have had the short/long collision above, or the
    grouping that preceded the ambiguity guard entirely. The caller re-derives
    rather than trusting it.
    """
    try:
        data = json.loads(sidecar_path(pdb_path).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return False
    if not isinstance(data, dict) or data.get("schema") != _SIDECAR_SCHEMA:
        return False
    source = data.get("source")
    if not isinstance(source, dict):
        return False
    try:
        return source.get("sha256") == _source_stamp(cif_path)["sha256"]
    except OSError:
        return False


def cif_to_pdb(cif_path: Path, out_path: Path) -> Path:
    """Write the first model of ``cif_path`` as a PDB file at ``out_path``."""
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict

    data = MMCIF2Dict(str(cif_path))
    groups = _column(data, "_atom_site.group_PDB")
    if not groups:
        raise ValueError(f"{cif_path.name}: no _atom_site records")

    n = len(groups)
    atom_names = _column(data, "_atom_site.auth_atom_id", "_atom_site.label_atom_id")
    alt_ids = _column(data, "_atom_site.label_alt_id") or ["."] * n
    comp_ids = _column(data, "_atom_site.auth_comp_id", "_atom_site.label_comp_id")
    asym_ids = _column(data, "_atom_site.auth_asym_id", "_atom_site.label_asym_id")
    seq_ids = _column(data, "_atom_site.auth_seq_id", "_atom_site.label_seq_id")
    ins_codes = _column(data, "_atom_site.pdbx_PDB_ins_code") or ["?"] * n
    xs = _column(data, "_atom_site.Cartn_x")
    ys = _column(data, "_atom_site.Cartn_y")
    zs = _column(data, "_atom_site.Cartn_z")
    occs = _column(data, "_atom_site.occupancy") or ["1.00"] * n
    bfacs = _column(data, "_atom_site.B_iso_or_equiv") or ["0.00"] * n
    elements = _column(data, "_atom_site.type_symbol") or [""] * n
    models = _column(data, "_atom_site.pdbx_PDB_model_num") or ["1"] * n

    first_model = models[0]
    unique_chains: list[str] = []
    for chain_id, model in zip(asym_ids, models):
        if model == first_model and chain_id not in unique_chains:
            unique_chains.append(chain_id)
    if len(unique_chains) > len(_PDB_CHAIN_ALPHABET):
        raise StructureTooLargeError(
            f"{cif_path.name}: {len(unique_chains)} chains exceed the {len(_PDB_CHAIN_ALPHABET)} "
            "the PDB format can hold; work from the mmCIF directly "
            f"({cif_path})"
        )
    # Keep every author chain ID the PDB format can actually hold, and rename
    # only the ones it cannot. Assigning the alphabet by arrival order renamed
    # all of them -- auth ["H","L","A"] became {"H":"A","L":"B","A":"C"} -- so
    # `--chain A` prepared the antibody heavy chain while the docstring above
    # promised author numbering was preserved, and `ligand_code`, UniProt
    # cross-references and anything a user reads off the PDB entry pointed at
    # the wrong polymer. Only a multi-character label is unrepresentable.
    chain_map: dict[str, str] = {}
    taken = {c for c in unique_chains if len(c) == 1 and c in _PDB_CHAIN_ALPHABET}
    spare = (c for c in _PDB_CHAIN_ALPHABET if c not in taken)
    for src in unique_chains:
        if src in taken:
            chain_map[src] = src
        else:
            chain_map[src] = next(spare)
    renamed = {k: v for k, v in chain_map.items() if k != v}
    if renamed:
        log.warning("%s: mmCIF chain IDs %s cannot be written in PDB format and "
                    "were remapped; every other chain keeps its author ID",
                    cif_path.name, renamed)

    lines: list[str] = []
    serial = 0
    for i in range(n):
        if models[i] != first_model:
            continue
        serial += 1
        record = "ATOM  " if groups[i] == "ATOM" else "HETATM"
        alt = " " if alt_ids[i] in (".", "?") else alt_ids[i][:1]
        icode = " " if ins_codes[i] in (".", "?") else ins_codes[i][:1]
        try:
            resnum = int(seq_ids[i])
        except ValueError:  # label_seq_id is '.' for waters and ligands
            resnum = 0
        # `% 10000` silently aliased anything the PDB `resSeq` field cannot
        # hold: an expression-tag residue at -5 came out as 9995 and 10001 as
        # 1, so `_chain_ca_by_resnum` and `_ligand_coords` looked up keys that
        # no longer matched the author numbering this module's docstring
        # promises to preserve. Refuse instead, the way the chain-count limit
        # already does -- an unrepresentable structure is not a renumbered one.
        if not -999 <= resnum <= 9999:
            raise StructureTooLargeError(
                f"{cif_path.name}: auth_seq_id {seq_ids[i]} does not fit the "
                f"PDB resSeq field (-999..9999); work from the mmCIF directly "
                f"({cif_path})")
        try:
            occ, bfac = float(occs[i]), float(bfacs[i])
        except ValueError:
            occ, bfac = 1.0, 0.0
        lines.append(
            f"{record}{serial % 100000:5d} "
            f"{_format_atom_name(atom_names[i], elements[i])}"
            f"{alt}{comp_ids[i][:3]:>3s} {chain_map[asym_ids[i]]}{resnum:4d}{icode}   "
            f"{float(xs[i]):8.3f}{float(ys[i]):8.3f}{float(zs[i]):8.3f}"
            f"{occ:6.2f}{bfac:6.2f}          {elements[i].strip()[:2]:>2s}"
        )


















    out_path.parent.mkdir(parents=True, exist_ok=True)
    write_comp_id_sidecar(cif_path, out_path, comp_ids)
    out_path.write_text("\n".join(lines) + "\nEND\n")
    log.info("converted %s -> %s (%d atoms, %d chains)", cif_path.name, out_path.name, serial, len(unique_chains))
    return out_path
