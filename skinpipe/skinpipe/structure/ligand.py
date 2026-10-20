"""Recover a chemically correct molecule from crystallographic HETATM records.

PDB coordinate files carry no bond orders. ``Chem.MolFromPDBBlock`` therefore
perceives connectivity from interatomic distances and leaves every bond single,
which is useless as an RMSD reference: aromatic rings come out as cyclohexanes
and symmetry-equivalent atoms stop being equivalent.

The fix is to pull the idealised chemical-component SDF from RCSB - which *does*
carry bond orders - and transfer them onto the crystallographic coordinates with
``AssignBondOrdersFromTemplate``. The result has the experimental geometry and
the correct chemistry, so ``GetBestRMS`` can do symmetry-corrected matching
against a docked pose.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Optional

from rdkit import Chem
from rdkit.Chem import AllChem

from . import clean, rcsb

log = logging.getLogger(__name__)


def cocrystal_pdb_block(raw_text: str, lig: dict[str, Any]) -> str:
    """HETATM records belonging to one cocrystal ligand copy."""
    resnum = str(lig["resnum"])
    # Insertion code, column 27. `list_cocrystal_ligands` groups copies by
    # (resname, chain, resnum, icode) and reports them separately, so leaving it
    # out here merges every copy that shares a residue number into one molecule.
    # Reproduced on a synthetic entry with LIG A 301A and LIG A 301B, three atoms
    # each: the listing sees two 3-atom copies and the extraction returns six
    # atoms for either, spanning 102 A and yielding a 112 A box. A box that size
    # is not a pocket, and the redock reference would be two copies at once.
    icode = str(lig.get("icode") or "").strip()
    lines = [
        line
        for line in clean.atom_lines(raw_text)
        if line.startswith("HETATM")
        and line[17:20].strip().upper() == lig["resname"]
        and line[21] == lig["chain"]
        and line[22:26].strip() == resnum
        and line[26].strip() == icode
        and line[16] in (" ", "A")
    ]
    return "\n".join(lines) + "\nEND\n"


def _full_component_id(lig: dict[str, Any],
                        comp_id_map: Optional[dict[str, Any]]) -> Optional[str]:
    """The component's real CCD id, recovering a PDB-truncated one when possible.

    `mmcif.cif_to_pdb` writes `{pdb}.comp_ids.json` mapping each three-character
    truncation back to the full id. Without the map the truncated code is all
    there is, and the element-composition guard downstream is what catches a
    mismatch.
    """
    written = str(lig.get("resname") or "")
    if lig.get("component_id"):
        return str(lig["component_id"])
    if comp_id_map and written in comp_id_map:
        full = comp_id_map[written]
        # None marks a prefix shared by several components: not resolvable.
        return None if full is None else str(full)
    return written


def extract_cocrystal_mol(
    raw_text: str,
    lig: dict[str, Any],
    out_sdf: Optional[Path] = None,
    cache_dir: Optional[Path] = None,
    comp_id_map: Optional[dict[str, Optional[str]]] = None,
) -> Optional[Chem.Mol]:
    """Crystallographic ligand as a bond-order-correct RDKit molecule.

    Returns ``None`` (never raises) when RCSB has no chemical-component template
    or the template cannot be matched onto the observed atoms - a partially
    occupied or truncated ligand is common and must not abort receptor prep.
    """








    resname = _full_component_id(lig, comp_id_map)
    if resname is None:
        # More than one component in this entry truncates to the same three
        # characters, so the PDB file cannot say which one this copy is. Guessing
        # would fetch a template for the wrong molecule, which is the defect the
        # resolution exists to prevent.
        lig["template_ambiguous"] = True
        log.warning("%s: more than one component truncates to this code in %s; "
                    "the ligand cannot be identified from the PDB file and no "
                    "template is fetched", lig["resname"],
                    lig.get("chain") or "?")
        return None
    if resname != lig["resname"]:
        log.info("%s is the truncated form of component %s; using the full id for "
                 "the chemical-component template", lig["resname"], resname)
        lig["component_id"] = resname
    block = cocrystal_pdb_block(raw_text, lig)
    observed = Chem.MolFromPDBBlock(block, sanitize=False, removeHs=False, proximityBonding=True)
    if observed is None:
        log.warning("could not parse HETATM block for %s", lig["resname"])
        return None

    try:
        sdf_text = rcsb.fetch_ligand_sdf_text(resname, cache_dir=cache_dir)
    except rcsb.RCSBError as exc:
        # An outage is not an absence. Every other `return None` here means the
        # structure cannot yield a usable ligand -- a truncated or unmatchable
        # one -- and is permanent for that entry; this one means RCSB did not
        # answer and the same call would succeed on a retry. They reached the
        # caller as the same value, and the caller reports "no crystallographic
        # ligand available", which sends the reader to look for a ligand that
        # is sitting in the PDB file. Observed on 1HFC: PLH is present with 25
        # atoms, the template fetch timed out, and the redock control was
        # recorded as impossible rather than as deferred.
        log.warning("chemical component %s: %s", lig["resname"], exc)
        lig["template_fetch_error"] = str(exc)
        return None
    if not sdf_text:
        lig["template_fetch_error"] = (
            f"RCSB returned an empty chemical-component template for "
            f"{lig['resname']}")
        return None

    template = Chem.MolFromMolBlock(sdf_text, sanitize=True, removeHs=True)
    if template is None:
        log.warning("unparsable ideal SDF for %s", lig["resname"])
        return None

    observed = Chem.RemoveHs(observed, sanitize=False)
    try:
        mol = AllChem.AssignBondOrdersFromTemplate(template, observed)
    except (ValueError, RuntimeError) as exc:







        log.info(
            "bond-order transfer failed for %s (template %d heavy atoms, "
            "observed %d): %s; retrying by CCD atom name",
            lig["resname"], template.GetNumHeavyAtoms(),
            observed.GetNumHeavyAtoms(), exc)
        try:
            cif_text = rcsb.fetch_ligand_cif_text(resname, cache_dir=cache_dir)
        except rcsb.RCSBError as cif_exc:
            lig["template_fetch_error"] = str(cif_exc)
            return None
        mol = (mol_from_ccd_names(raw_text, lig, cif_text) if cif_text else None)
        if mol is None:
            log.warning("could not rebuild %s from CCD atom names either",
                        lig["resname"])
            return None
        lig["built_from"] = "ccd_atom_names"
        if out_sdf is not None:
            out_sdf.parent.mkdir(parents=True, exist_ok=True)
            writer = Chem.SDWriter(str(out_sdf))
            try:
                writer.write(mol)
            finally:
                writer.close()
        return mol

    mol.SetProp("_Name", f"{lig['resname']}_{lig['chain']}{lig['resnum']}")
    try:
        Chem.SanitizeMol(mol)
    except (Chem.KekulizeException, ValueError) as exc:
        log.warning("sanitisation of crystallographic %s failed: %s", lig["resname"], exc)
        return None




    tmpl_heavy = template.GetNumHeavyAtoms()
    matched = mol.GetNumHeavyAtoms()
    lig["ccd_heavy_atoms"] = tmpl_heavy
    lig["matched_heavy_atoms"] = matched
    lig["observed_heavy_atoms"] = matched
    lig["unmatched_atom_names"] = []
    lig["completeness"] = (matched / tmpl_heavy) if tmpl_heavy else None
    lig["is_fragment"] = bool(tmpl_heavy and matched < tmpl_heavy)
    lig.setdefault("built_from", "template_bond_orders")

    if out_sdf is not None:
        out_sdf.parent.mkdir(parents=True, exist_ok=True)
        writer = Chem.SDWriter(str(out_sdf))
        try:
            writer.write(mol)
        finally:
            writer.close()
        log.debug("wrote crystallographic ligand %s", out_sdf)
    return mol


#: CCD bond orders, as `_chem_comp_bond.value_order` spells them.
_CCD_ORDER = {"SING": Chem.BondType.SINGLE, "DOUB": Chem.BondType.DOUBLE,
              "TRIP": Chem.BondType.TRIPLE, "QUAD": Chem.BondType.QUADRUPLE,
              "AROM": Chem.BondType.AROMATIC}


def _ccd_loop(cif: str, prefix: str) -> list[dict[str, str]]:
    """Rows of one `loop_` in a CCD cif, keyed by their column names.

    A deliberately small reader: CCD component files are flat, one-line-per-row
    loops with no multi-line values, and pulling in a full mmCIF parser for two
    tables would be a dependency for nothing.
    """
    lines = cif.splitlines()
    cols: list[str] = []
    rows: list[dict[str, str]] = []
    in_loop = False
    for line in lines:
        s = line.strip()
        if s.startswith(prefix + "."):
            cols.append(s.split(".", 1)[1].split()[0])
            in_loop = True
            continue
        if not in_loop:
            continue
        if not s or s.startswith("#") or s.startswith("_") or s.startswith("loop_"):
            if cols and rows:
                break
            if s.startswith("_") or s.startswith("loop_"):
                in_loop = False
                cols = []
            continue
        parts = _cif_tokens(s)
        if len(parts) < len(cols):
            continue
        rows.append(dict(zip(cols, parts[:len(cols)])))
    return rows


def _cif_tokens(line: str) -> list[str]:
    """Split one mmCIF data line, honouring quoting.

    `str.split()` is wrong here and quietly so. CCD quotes any atom id that
    needs it -- the primed names of sugars and nucleotides, `C1'`, `O2'` -- and
    a bare split keeps the quotes, so `"C1'"` never matches the PDB atom name
    `C1'`. The ligand is then rebuilt from whichever atoms happened to be
    unquoted: one atom and no bonds where there should have been two and one,
    and `redock_validation` compares against a fragment while reporting a
    ligand. Single and double quotes both appear in real files.
    """
    out: list[str] = []
    i, n = 0, len(line)
    while i < n:
        if line[i].isspace():
            i += 1
            continue
        if line[i] in "'\"":
            quote = line[i]
            i += 1
            start = i
            while i < n and line[i] != quote:
                i += 1
            out.append(line[start:i])
            i += 1                      # past the closing quote
        else:
            start = i
            while i < n and not line[i].isspace():
                i += 1
            out.append(line[start:i])
    return out


def mol_from_ccd_names(raw_text: str, lig: dict[str, Any],
                       cif_text: str) -> Optional[Chem.Mol]:
    'Observed ligand built by matching PDB atom names to the CCD definition.'
    atoms = _ccd_loop(cif_text, "_chem_comp_atom")
    bonds = _ccd_loop(cif_text, "_chem_comp_bond")
    if not atoms or not bonds:
        log.warning("CCD cif for %s has no atom or bond loop", lig["resname"])
        return None

    element = {a["atom_id"]: a.get("type_symbol", "").capitalize() for a in atoms}
    charge = {a["atom_id"]: a.get("charge", "0") for a in atoms}









    seen: dict[str, tuple[float, float, float]] = {}
    unmatched: list[str] = []
    n_observed_heavy = 0
    for line in cocrystal_pdb_block(raw_text, lig).splitlines():
        if not line.startswith(("HETATM", "ATOM  ")):
            continue
        name = line[12:16].strip()
        known = element.get(name)
        if known is None:
            # Element from the PDB's own columns 77-78 when the CCD cannot name
            # it, so a hydrogen is still excluded and anything else is counted
            # as the mismatch it is.
            pdb_elem = (line[76:78].strip() or name[:1]).capitalize()
            if pdb_elem != "H":
                unmatched.append(name)
                n_observed_heavy += 1
            continue
        if known == "H":
            continue                      # hydrogens are not part of the metric
        n_observed_heavy += 1
        try:
            seen[name] = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            continue
    if not seen:
        log.warning("no CCD-named heavy atoms observed for %s", lig["resname"])
        return None
    if unmatched:
        log.warning("%s: %d observed heavy atom(s) are absent from the CCD table "
                    "(%s); the reference is built from the matched ones only and "
                    "is not the ligand as deposited", lig["resname"],
                    len(unmatched), ", ".join(unmatched[:6]))

    rw = Chem.RWMol()
    index: dict[str, int] = {}
    conf = Chem.Conformer(len(seen))
    for name, xyz in seen.items():
        # `Chem.Atom` raises on a symbol RDKit does not know -- CCD writes `D`
        # for deuterium -- and this helper is documented to return None rather
        # than raise, because every caller treats None as "no usable ligand".
        try:
            atom = Chem.Atom(element.get(name) or "C")
        except (RuntimeError, ValueError) as exc:
            log.warning("%s: RDKit will not make an atom for element %r (%s); "
                        "no usable crystal reference", lig["resname"],
                        element.get(name), exc)
            return None
        try:
            atom.SetFormalCharge(int(charge.get(name, "0")))
        except ValueError:
            pass
        index[name] = rw.AddAtom(atom)
        conf.SetAtomPosition(index[name], xyz)

    n_dropped = 0
    for b in bonds:
        a1, a2 = b.get("atom_id_1"), b.get("atom_id_2")
        if a1 not in index or a2 not in index:
            n_dropped += 1
            continue
        order = _CCD_ORDER.get((b.get("value_order") or "SING").upper(),
                               Chem.BondType.SINGLE)
        rw.AddBond(index[a1], index[a2], order)

    mol = rw.GetMol()
    mol.AddConformer(conf, assignId=True)
    try:
        Chem.SanitizeMol(mol)
    except (Chem.KekulizeException, ValueError, RuntimeError) as exc:
        log.warning("CCD-built %s did not sanitise: %s", lig["resname"], exc)
        return None
    mol.SetProp("_Name", f"{lig['resname']}_{lig['chain']}{lig.get('resnum') or ''}")

    # How much of the component was actually modelled. Rebuilding from whatever
    # atoms the crystallographer could place is what makes a partially occupied
    # ligand usable at all, and it is also what makes the result a *fragment*
    # when occupancy is poor: a redock measured against 18 of 30 heavy atoms is
    # a weaker statement than one measured against the ligand, and nothing said
    # which had happened. Recorded on the molecule and on the listing dict so
    # `redock_validation` and the report can both read it.
    ccd_heavy = sum(1 for a in atoms
                    if (a.get("type_symbol", "") or "").capitalize() != "H")
    matched_heavy = mol.GetNumAtoms()
    completeness = (matched_heavy / ccd_heavy) if ccd_heavy else None
    mol.SetProp("ccd_heavy_atoms", str(ccd_heavy))
    mol.SetProp("matched_heavy_atoms", str(matched_heavy))
    mol.SetProp("observed_heavy_atoms", str(n_observed_heavy))
    mol.SetProp("unmatched_atom_names", ",".join(unmatched))
    if completeness is not None:
        mol.SetProp("completeness", f"{completeness:.3f}")
    lig["ccd_heavy_atoms"] = ccd_heavy
    lig["matched_heavy_atoms"] = matched_heavy
    lig["observed_heavy_atoms"] = n_observed_heavy
    lig["unmatched_atom_names"] = list(unmatched)
    lig["completeness"] = completeness
    # Either direction makes the reference something other than the deposited
    # ligand: fewer matched atoms than the component has, or observed atoms the
    # component does not name. Both have to reach the same flag, because both
    # make a redock a statement about something else.
    lig["is_fragment"] = bool(
        (completeness is not None and completeness < 1.0) or unmatched)
    if n_dropped:
        log.info("%s: %d CCD bonds dropped because an endpoint was not observed "
                 "(partially occupied ligand)", lig["resname"], n_dropped)
    if lig["is_fragment"]:
        log.warning("%s: rebuilt from %d of %d CCD heavy atoms (%.0f%%), with %d "
                    "observed heavy atom(s) unmatched; a control measured against "
                    "this is a statement about a fragment, not about the whole "
                    "ligand", lig["resname"], matched_heavy, ccd_heavy,
                    100 * (completeness or 0.0), len(unmatched))
    return mol
