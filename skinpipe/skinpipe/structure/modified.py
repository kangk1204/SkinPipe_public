"""Native modified amino acids declared by PDB MODRES, never loose cofactors.

Chemical graphs and hydrogen variants come from Meeko's installed curated
templates. No component download or inference is performed by this module.
Unsupported, incomplete or terminal variants are reported rather than changed
to their standard parent. A template's protonation state is an assumption, not
a site-specific pKa calculation.
"""
from __future__ import annotations

from collections import defaultdict
from functools import lru_cache

STANDARD_AMINO_ACIDS = frozenset({
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
})
PHOSPHORYLATED_AMINO_ACIDS = frozenset({"PTR", "SEP", "TPO"})


def declarations(text: str, chain: str | None = None) -> dict[tuple[str, str, str], str]:
    """Exact modified-residue identifiers and their declared standard parent."""
    result = {}
    for line in text.splitlines():
        if not line.startswith("MODRES"):
            continue
        name, cid, number, parent = line[12:15].strip(), line[16], line[18:23].strip(), line[24:27].strip()
        if chain is not None and cid != chain:
            continue
        if parent not in STANDARD_AMINO_ACIDS or name == parent:
            continue
        key = (cid, number, name)
        if key in result and result[key] != parent:
            raise ValueError(f"conflicting MODRES parent for {key}")
        result[key] = parent
    return result


def declaration_lines(text: str, chain: str | None = None) -> list[str]:
    wanted = declarations(text, chain)
    return [line for line in text.splitlines() if line.startswith("MODRES")
            and (line[16], line[18:23].strip(), line[12:15].strip()) in wanted]


def groups(text: str, *, known_phospho: bool = False) -> dict[tuple[str, str, str], list[str]]:
    wanted = declarations(text)
    result = defaultdict(list)
    for line in text.splitlines():
        if line.startswith("ENDMDL"):
            break
        if line.startswith(("ATOM  ", "HETATM")) and line[16] in (" ", "A"):
            key = (line[21], line[22:27].strip(), line[17:20].strip())
            if key in wanted or (known_phospho and key[2] in PHOSPHORYLATED_AMINO_ACIDS):
                result[key].append(line)
    if known_phospho:
        # External PDB writers may omit MODRES. A known phospho component
        # carrying the peptide backbone must not be silently stripped either.
        result = {key: lines for key, lines in result.items()
                  if key in wanted or {"N", "CA", "C"} <= heavy_names(lines)}
    return dict(result)


@lru_cache(maxsize=1)
def _templates():
    from meeko import ResidueChemTemplates
    return ResidueChemTemplates.create_from_defaults()


def template(name: str, names: set[str]):
    """A curated internal amino-acid template with precisely these heavy atoms."""
    value = _templates().residue_templates.get(name)
    if value is None or set(value.link_labels.values()) != {"N-term", "C-term"}:
        raise ValueError(f"native modified residue {name} has no curated internal polymer template")
    heavy = {value.atom_names[a.GetIdx()] for a in value.mol.GetAtoms() if a.GetAtomicNum() != 1}
    if names != heavy:
        raise ValueError(f"native modified residue {name} does not match complete template; "
                         f"missing={sorted(heavy - names)}, extra={sorted(names - heavy)}")
    return value


def heavy_names(lines: list[str]) -> set[str]:
    return {line[12:16].strip() for line in lines if line[76:78].strip() not in ("H", "D")}


def retained_conect(text: str, retained: list[str]) -> list[str]:
    """Keep explicit original bonds only when both endpoints were selected."""
    serials = {int(line[6:11]) for line in retained if line.startswith(("ATOM  ", "HETATM"))}
    result = []
    for line in text.splitlines():
        if not line.startswith("CONECT"):
            continue
        values = [int(line[i:i+5]) for i in range(6, len(line.rstrip()), 5) if line[i:i+5].strip()]
        if values and values[0] in serials:
            values = [values[0]] + [number for number in values[1:] if number in serials]
            if len(values) > 1:
                result.append("CONECT" + "".join(f"{number:5d}" for number in values))
    return result


def register_templates(fixer, native: dict[tuple[str, str, str], list[str]]) -> list[str]:
    """Use the observed complete heavy geometry as a public PDBFixer template."""
    from openmm import app, unit
    from rdkit import Chem

    notes = []
    for key, lines in native.items():
        value = template(key[2], heavy_names(lines))
        by_name = {line[12:16].strip(): line for line in lines if line[76:78].strip() not in ("H", "D")}
        if len(by_name) != sum(line[76:78].strip() not in ("H", "D") for line in lines):
            raise ValueError(f"duplicate native modified atom identifier in {key}")
        topology = app.Topology()
        residue = topology.addResidue(key[2], topology.addChain(key[0]))
        atoms, positions = {}, []
        for atom in value.mol.GetAtoms():
            if atom.GetAtomicNum() == 1:
                continue
            name = value.atom_names[atom.GetIdx()]
            line = by_name[name]
            if line[76:78].strip().upper() != atom.GetSymbol().upper():
                raise ValueError(f"native modified atom element disagrees with template: {key}:{name}")
            atoms[name] = topology.addAtom(name, app.Element.getBySymbol(atom.GetSymbol()), residue)
            positions.append(tuple(float(line[a:b]) for a, b in ((30,38),(38,46),(46,54))))
        for bond in value.mol.GetBonds():
            a, b = (value.atom_names[i] for i in (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
            if a in atoms and b in atoms:
                topology.addBond(atoms[a], atoms[b])
        fixer.registerTemplate(topology, positions * unit.angstrom)
        notes.append(f"retained native modified polymer {':'.join(key)}; curated Meeko "
                     f"template {key[2]} formal charge {Chem.GetFormalCharge(value.mol):+d}; "
                     "protonation is template-defined, not a site-specific pKa estimate")
    validate_native_bonds(fixer.topology, native)
    return notes


def remove_other_heterogens(fixer, native: dict[tuple[str, str, str], list[str]], *, keep_water: bool = False) -> None:
    from openmm import app
    keep = STANDARD_AMINO_ACIDS | {"N", "UNK", "A", "C", "G", "U", "DA", "DC", "DG", "DT"}
    if keep_water:
        keep |= {"HOH", "WAT", "H2O", "DOD"}
    remove = [residue for residue in fixer.topology.residues()
              if residue.name not in keep
              and (residue.chain.id, residue.id + residue.insertionCode.strip(), residue.name) not in native]
    modeller = app.Modeller(fixer.topology, fixer.positions)
    modeller.delete(remove)
    fixer.topology, fixer.positions = modeller.topology, modeller.positions


def add_hydrogens(fixer, native: dict[tuple[str, str, str], list[str]], ph: float) -> None:
    """Explicit template hydrogens for modified residues; usual pH rules elsewhere."""
    from openmm import app
    variants = []
    for residue in fixer.topology.residues():
        key = (residue.chain.id, residue.id + residue.insertionCode.strip(), residue.name)
        if key not in native:
            variants.append(None)
            continue
        value = template(residue.name, heavy_names(native[key]))
        hydrogens = []
        for atom in value.mol.GetAtoms():
            if atom.GetAtomicNum() != 1:
                continue
            parent, = atom.GetNeighbors()
            hydrogens.append((value.atom_names[atom.GetIdx()], value.atom_names[parent.GetIdx()]))
        variants.append(hydrogens)
    modeller = app.Modeller(fixer.topology, fixer.positions)
    modeller.addHydrogens(pH=float(ph), variants=variants, platform=fixer.platform)
    fixer.topology, fixer.positions = modeller.topology, modeller.positions


def expected_heavy_bonds(name: str) -> set[tuple[str, str]] | None:
    if name in STANDARD_AMINO_ACIDS:
        return None
    value = _templates().residue_templates.get(name)
    if value is None:
        return None
    return {tuple(sorted((value.atom_names[bond.GetBeginAtomIdx()], value.atom_names[bond.GetEndAtomIdx()])))
            for bond in value.mol.GetBonds()
            if all(value.mol.GetAtomWithIdx(i).GetAtomicNum() != 1
                   for i in (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))}


def validate_native_bonds(topology, native: dict[tuple[str, str, str], list[str]]) -> None:
    """Require the complete native graph and both original peptide links."""
    def key(atom):
        residue = atom.residue
        return (residue.chain.id, residue.id + residue.insertionCode.strip(), residue.name)
    bonds = list(topology.bonds())
    for identifier, lines in native.items():
        value = template(identifier[2], heavy_names(lines))
        expected = {tuple(sorted((value.atom_names[b.GetBeginAtomIdx()], value.atom_names[b.GetEndAtomIdx()])))
                    for b in value.mol.GetBonds()
                    if all(value.mol.GetAtomWithIdx(i).GetAtomicNum() != 1
                           for i in (b.GetBeginAtomIdx(), b.GetEndAtomIdx()))}
        actual, links = set(), {"N": [], "C": []}
        for a, b in bonds:
            if a.element.symbol == "H" or b.element.symbol == "H":
                continue
            ka, kb = key(a), key(b)
            if ka == kb == identifier:
                actual.add(tuple(sorted((a.name, b.name))))
            elif ka == identifier or kb == identifier:
                inside, outside = (a, b) if ka == identifier else (b, a)
                if inside.name not in links:
                    raise ValueError(f"unsupported additional covalent link in native modified residue {identifier}")
                links[inside.name].append(outside.name)
        if actual != expected or links != {"N": ["C"], "C": ["N"]}:
            raise ValueError(f"native modified residue {identifier} requires a complete internal graph "
                             f"and both peptide links; missing={sorted(expected - actual)}, "
                             f"extra={sorted(actual - expected)}, external={links}; "
                             "terminal or unrecorded covalent chemistry is unsupported")


def deduplicate_repair_bonds(fixer) -> list[str]:
    """Remove identical edges copied twice by PDBFixer, using public APIs.

    addMissingAtoms copies explicit native CONECT edges, then creates standard
    peptide edges; a native C--standard N edge can be inserted twice. Two
    records would create duplicate bonded force terms. Disagreeing bond
    attributes are not treated as redundant.
    """
    from openmm import app
    original = fixer.topology
    unique, duplicates = {}, 0
    for bond in original.bonds():
        pair = tuple(sorted((bond[0].index, bond[1].index)))
        if pair in unique:
            prior = unique[pair]
            if (prior.type, prior.order) != (bond.type, bond.order):
                raise ValueError(f"conflicting duplicate bond attributes for atom indices {pair}")
            duplicates += 1
        else:
            unique[pair] = bond
    if not duplicates:
        return []
    topology, atoms = app.Topology(), {}
    for chain in original.chains():
        copied_chain = topology.addChain(chain.id)
        for residue in chain.residues():
            copied_residue = topology.addResidue(residue.name, copied_chain, residue.id, residue.insertionCode)
            for atom in residue.atoms():
                atoms[atom.index] = topology.addAtom(atom.name, atom.element, copied_residue,
                                                    atom.id, formalCharge=atom.formalCharge)
    for bond in unique.values():
        topology.addBond(atoms[bond[0].index], atoms[bond[1].index], type=bond.type, order=bond.order)
    topology.setPeriodicBoxVectors(original.getPeriodicBoxVectors())
    assert topology.getNumAtoms() == original.getNumAtoms()
    fixer.topology = topology
    return [f"removed {duplicates} identical duplicate topology bond record(s) introduced "
            "during native residue repair; atom order, identifiers and positions unchanged"]


def amber_parameter_files(names: set[str]) -> tuple[str, ...]:
    """Installed, graph-checked phosaa14SB supplement to Amber14."""
    if not names:
        return ()
    if names - PHOSPHORYLATED_AMINO_ACIDS:
        raise ValueError(f"native modified Amber14 parameters unsupported for {sorted(names)}")
    from importlib.resources import files
    import xml.etree.ElementTree as ET
    from rdkit import Chem

    path = files("openmmforcefields").joinpath("ffxml/amber/phosaa14SB.xml")
    if not path.is_file():
        raise ValueError("installed openmmforcefields has no phosaa14SB.xml")
    root = ET.fromstring(path.read_bytes())
    for name in names:
        residue = next((r for r in root.findall("Residues/Residue") if r.get("name") == name), None)
        value = _templates().residue_templates.get(name)
        if residue is None or value is None:
            raise ValueError(f"native modified parameter template absent for {name}")
        atom_names = {a.get("name") for a in residue.findall("Atom")}
        graph = {tuple(sorted((b.get("atomName1"), b.get("atomName2")))) for b in residue.findall("Bond")}
        expected = {tuple(sorted((value.atom_names[b.GetBeginAtomIdx()], value.atom_names[b.GetEndAtomIdx()])))
                    for b in value.mol.GetBonds()}
        charge = sum(float(a.get("charge")) for a in residue.findall("Atom"))
        if (atom_names != set(value.atom_names) or graph != expected
                or abs(charge - Chem.GetFormalCharge(value.mol)) > 1e-5):
            raise ValueError(f"phosaa14SB and curated Meeko chemistry disagree for {name}")
    return (str(path),)


def parameterization_descriptor(text: str) -> dict | None:
    """Content-bound native receptor model for result provenance and reuse."""
    native = groups(text, known_phospho=True)
    if not native:
        return None
    from hashlib import sha256
    from importlib.metadata import version
    from pathlib import Path
    from rdkit import Chem
    paths = amber_parameter_files({key[2] for key in native})
    residues = []
    for key, lines in sorted(native.items()):
        value = template(key[2], heavy_names(lines))
        residues.append({"chain": key[0], "residue": key[1], "name": key[2],
                         "formal_charge": Chem.GetFormalCharge(value.mol)})
    normalized, base_path, phosphorus_paths = _amber_alias_xml({key[2] for key in native})
    return {"protein_forcefield": "amber14-all.xml", "supplement": "phosaa14SB",
            "openmmforcefields_version": version("openmmforcefields"),
            "files": [{"name": Path(p).name, "sha256": sha256(Path(p).read_bytes()).hexdigest()} for p in paths],
            "base_protein_xml_sha256": sha256(base_path.read_bytes()).hexdigest(),
            "shared_phosphorus_type_sources": [{"name": p.name, "sha256": sha256(p.read_bytes()).hexdigest()}
                                                for p in phosphorus_paths],
            "type_alias_protocol": "phosaa14SB_to_OpenMM_amber14_protein_types_v1",
            "normalized_supplement_xml_sha256": sha256(normalized.encode()).hexdigest(),
            "native_residues": residues,
            "protonation_scope": "curated template state; no site-specific pKa calculation"}


def _amber_alias_xml(names: set[str]):
    """Alias external atom type names to built-in Amber14; preserve numbers.

    phosaa14SB references bare N/C/CX/etc types. OpenMM Amber14 prefixes
    *type names* with protein-, but retains the same bare *classes*. Keep the
    supplement's own new types and every charge and bonded/nonbonded numeric
    parameter unchanged. Refuse an unexpected base naming convention.
    """
    from importlib.resources import files
    from pathlib import Path
    import xml.etree.ElementTree as ET
    original, = amber_parameter_files(names)
    root = ET.parse(original).getroot()
    base_path = Path(str(files("openmm.app").joinpath("data/amber14/protein.ff14SB.xml")))
    base = ET.parse(base_path).getroot()
    base_types = {a.get("name"): a for a in base.findall("AtomTypes/Type")}
    declared = {a.get("name") for a in root.findall("AtomTypes/Type")}
    # Standard protein residues have no P, so protein.ff14SB.xml omits it.
    # Amber14-all already loads the shared Amber P type as RNA-P. Require
    # its element, mass, class and LJ constants to equal official ff14SB P.
    rna_path = Path(str(files("openmm.app").joinpath("data/amber14/RNA.OL3.xml")))
    amber_path = Path(str(files("openmmforcefields").joinpath("ffxml/amber/ff14SB.xml")))
    rna, amber = ET.parse(rna_path).getroot(), ET.parse(amber_path).getroot()
    rna_p = next((a for a in rna.findall("AtomTypes/Type") if a.get("name") == "RNA-P"), None)
    amber_p = next((a for a in amber.findall("AtomTypes/Type") if a.get("name") == "P"), None)
    rna_lj = next((a for a in rna.findall("NonbondedForce/Atom") if a.get("type") == "RNA-P"), None)
    amber_lj = next((a for a in amber.findall("NonbondedForce/Atom") if a.get("type") == "P" or a.get("class") == "P"), None)
    if (any(a is None for a in (rna_p, amber_p, rna_lj, amber_lj))
            or any(rna_p.get(a) != amber_p.get(a) for a in ("element", "mass", "class"))
            or any(float(rna_lj.get(a)) != float(amber_lj.get(a)) for a in ("sigma", "epsilon"))):
        raise ValueError("built-in Amber14 RNA-P is not the official shared Amber ff14SB phosphorus type")
    for element in root.iter():
        for field, value in list(element.attrib.items()):
            if field != "type" and not (field.startswith("type") and field[4:].isdigit()):
                continue
            if not value or value in declared:
                continue
            prefixed = "RNA-P" if value == "P" else "protein-" + value
            definition = rna_p if value == "P" else base_types.get(prefixed)
            if definition is None or definition.get("class") != value:
                raise ValueError(f"phosaa14SB external type {value} has no compatible built-in Amber14 alias")
            element.set(field, prefixed)
    return ET.tostring(root, encoding="unicode"), base_path, (rna_path, amber_path)


def amber_supplements(names: set[str]) -> tuple:
    if not names:
        return ()
    from io import StringIO
    xml, _, _ = _amber_alias_xml(names)
    stream = StringIO(xml)
    stream.name = "phosaa14SB.amber14-types.xml"
    return (stream,)


def repair_forcefields(topology) -> tuple:
    names = {r.name for r in topology.residues()} & PHOSPHORYLATED_AMINO_ACIDS
    return ("amber14-all.xml", *amber_supplements(names))


def restore_symmetric_names(pdbqt: str, input_pdb: str) -> tuple[str, list[str]]:
    """Restore exact observed labels only for a template graph automorphism."""
    from rdkit import Chem

    replacements = {}
    for key, source in groups(input_pdb).items():
        expected = {line[12:16].strip(): line[30:54] for line in source
                    if line[76:78].strip() not in ("H", "D")}
        value = template(key[2], set(expected))
        actual = {line[12:16].strip(): line[30:54] for line in pdbqt.splitlines()
                  if line.startswith(("ATOM  ", "HETATM"))
                  and (line[21], line[22:27].strip(), line[17:20].strip()) == key
                  and line.split()[-1] not in ("H", "HD", "HS")}
        if expected == actual:
            continue
        if expected.keys() != actual.keys() or len(set(expected.values())) != len(expected):
            continue  # The strict observed-atom check will report the failure.
        by_position = {position: name for name, position in expected.items()}
        if set(by_position) != set(actual.values()):
            continue
        permutation = {name: by_position[position] for name, position in actual.items()}
        heavy = Chem.RemoveHs(value.mol)
        names = [value.atom_names[atom.GetIdx()] for atom in value.mol.GetAtoms() if atom.GetAtomicNum() != 1]
        matches = heavy.GetSubstructMatches(heavy, uniquify=False, useChirality=True, maxMatches=1000)
        if not any(all(names[match[i]] == permutation[name] for i, name in enumerate(names)) for match in matches):
            raise RuntimeError(f"native modified residue {key} label permutation is not a chemical graph automorphism")
        for old, new in permutation.items():
            if old != new:
                replacements[(*key, old)] = new
    if not replacements:
        return pdbqt, []
    lines = []
    for line in pdbqt.splitlines(keepends=True):
        if line.startswith(("ATOM  ", "HETATM")):
            key = (line[21], line[22:27].strip(), line[17:20].strip(), line[12:16].strip())
            if key in replacements:
                line = line[:12] + f"{replacements[key]:^4s}" + line[16:]
        lines.append(line)
    return "".join(lines), [f"restored {len(replacements)} native modified atom names by exact "
                            "coordinates and template graph automorphism; coordinates, charges and types unchanged"]
