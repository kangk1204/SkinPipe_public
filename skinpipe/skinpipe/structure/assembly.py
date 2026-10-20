"""Deposited biological assemblies in the original target-chain frame.

An asymmetric unit is not an assembly: a nearby crystallographic neighbour
may be absent from the deposited biological assembly, and a real partner may
need an operator to appear. These helpers resolve the deposited definition,
expand its operators, and preserve the original primary-chain coordinates.
They do not establish physiological relevance, protonation or force-field
support. Unsupported/ambiguous definitions are errors rather than monomers.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import itertools
from pathlib import Path
import re

import numpy as np

from . import mmcif, modified


class AssemblyError(ValueError):
    """A source-bound assembly cannot be represented without ambiguity/loss."""


@dataclass
class AssemblyContext:
    pdb_text: str
    metadata: dict


def operator_products(expression: str, known: set[str], limit: int = 4096) -> list[tuple[str, ...]]:
    """Expand wwPDB unions/ranges/products; tuple order is left to right.

    Application order is right to left, e.g. (1,2)(3,4) applies 3/4 then 1/2.
    The limit prevents an unbounded Cartesian expansion; it never truncates.
    """
    value = re.sub(r"\s+", "", expression)
    if not value:
        raise AssemblyError("empty assembly operator expression")
    if "(" in value or ")" in value:
        groups = re.findall(r"\(([^()]+)\)", value)
        if not groups or "".join("(" + g + ")" for g in groups) != value:
            raise AssemblyError(f"invalid assembly operator expression: {expression}")
    else:
        groups = [value]
    choices = []
    size = 1
    for group in groups:
        identifiers = []
        for item in group.split(","):
            if item in known:
                identifiers.append(item)
            elif re.fullmatch(r"[0-9]+-[0-9]+", item):
                low, high = map(int, item.split("-"))
                if high < low or high - low + 1 > limit:
                    raise AssemblyError(f"invalid or excessive operator range: {item}")
                identifiers.extend(str(i) for i in range(low, high + 1))
            else:
                raise AssemblyError(f"unknown assembly operator: {item}")
        if not identifiers or any(i not in known for i in identifiers):
            raise AssemblyError(f"assembly operator range contains missing operators: {group}")
        identifiers = list(dict.fromkeys(identifiers))
        size *= len(identifiers)
        if size > limit:
            raise AssemblyError(f"assembly operator product exceeds {limit}; no instances truncated")
        choices.append(identifiers)
    return list(itertools.product(*choices))


def _rows(data: dict, category: str, fields: list[str]) -> list[dict]:
    columns = {field: mmcif._column(data, category + "." + field) for field in fields}
    lengths = {len(column) for column in columns.values()}
    if lengths == {0}:
        return []
    if len(lengths) != 1 or 0 in lengths:
        raise AssemblyError(f"missing or unequal {category} columns: {fields}")
    return [dict(zip(fields, values)) for values in zip(*columns.values())]


def _operators(data: dict) -> dict[str, np.ndarray]:
    fields = ["id"] + [f"matrix[{i}][{j}]" for i in range(1, 4) for j in range(1, 4)]
    fields += [f"vector[{i}]" for i in range(1, 4)]
    result = {}
    for row in _rows(data, "_pdbx_struct_oper_list", fields):
        matrix = np.eye(4)
        try:
            matrix[:3, :3] = [[float(row[f"matrix[{i}][{j}]"]) for j in range(1, 4)]
                               for i in range(1, 4)]
            matrix[:3, 3] = [float(row[f"vector[{i}]"]) for i in range(1, 4)]
        except ValueError as exc:
            raise AssemblyError("non-numeric assembly operator") from exc
        rotation = matrix[:3, :3]
        if (not np.isfinite(matrix).all()
                or not np.allclose(rotation.T @ rotation, np.eye(3), atol=1e-4, rtol=0)
                or not np.isclose(np.linalg.det(rotation), 1, atol=1e-4, rtol=0)):
            raise AssemblyError(f"operator {row['id']} is not a finite proper rigid transform")
        if row["id"] in result:
            raise AssemblyError(f"duplicate assembly operator {row['id']}")
        result[row["id"]] = matrix
    if not result:
        raise AssemblyError("no deposited assembly operators")
    return result


def _matrix_key(matrix: np.ndarray) -> tuple:
    return tuple(np.round(matrix, 12).ravel())


def _atom_key(line: str) -> tuple:
    return (line[21], line[22:26].strip(), line[26].strip(),
            line[17:20].strip(), line[12:16].strip(), line[16].strip())


def _atom_sources(data: dict, raw_text: str) -> tuple[dict, dict, dict]:
    """Match CIF label instances to original PDB atoms, not to author chains alone."""
    fields = ["group_PDB", "label_asym_id", "auth_asym_id", "auth_seq_id", "auth_comp_id",
              "auth_atom_id", "label_alt_id", "pdbx_PDB_ins_code", "Cartn_x", "Cartn_y",
              "Cartn_z", "type_symbol", "pdbx_PDB_model_num", "label_entity_id"]
    rows = _rows(data, "_atom_site", fields)
    if not rows:
        raise AssemblyError("no deposited atom records")
    occupancies = mmcif._column(data, "_atom_site.occupancy")
    if occupancies and len(occupancies) != len(rows):
        raise AssemblyError("unequal atom occupancy column")
    for i, row in enumerate(rows):
        row["occupancy"] = occupancies[i] if occupancies else None
    first_model = rows[0]["pdbx_PDB_model_num"]
    rows = [row for row in rows if row["pdbx_PDB_model_num"] == first_model]
    authors = list(dict.fromkeys(row["auth_asym_id"] for row in rows))
    alphabet = mmcif._PDB_CHAIN_ALPHABET
    if len(authors) > len(alphabet):
        raise AssemblyError("source author chains cannot fit in PDB format")
    reserved = {name for name in authors if name in alphabet and len(name) == 1}
    spare = iter(char for char in alphabet if char not in reserved)
    author_map = {name: name if name in reserved else next(spare) for name in authors}
    original = {}
    for line in raw_text.splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        key = _atom_key(line)
        if key in original:
            raise AssemblyError(f"duplicate source PDB atom: {key}")
        original[key] = line
    sources, seen, polymer = {}, set(), {}
    types = dict(zip(mmcif._column(data, "_entity_poly.entity_id"),
                     mmcif._column(data, "_entity_poly.type")))
    for row in rows:
        alt = "" if row["label_alt_id"] in (".", "?") else row["label_alt_id"]
        insertion = "" if row["pdbx_PDB_ins_code"] in (".", "?") else row["pdbx_PDB_ins_code"]
        key = (author_map[row["auth_asym_id"]], row["auth_seq_id"], insertion,
               row["auth_comp_id"][:3], row["auth_atom_id"], alt)
        line = original.get(key)
        if line is None or key in seen:
            raise AssemblyError(f"CIF/PDB atom identity mismatch or alias: {key}")
        seen.add(key)
        source_xyz = np.array([float(row[field]) for field in ("Cartn_x", "Cartn_y", "Cartn_z")])
        raw_xyz = np.array([float(line[a:b]) for a, b in ((30, 38), (38, 46), (46, 54))])
        if (not np.isfinite(source_xyz).all() or not np.isfinite(raw_xyz).all()
                or np.max(np.abs(source_xyz - raw_xyz)) > 0.0011
                or line[76:78].strip().upper() != row["type_symbol"].upper()):
            raise AssemblyError(f"CIF/PDB coordinates or element disagree: {key}")
        if alt:
            try:
                occupancy, pdb_occupancy = float(row["occupancy"]), float(line[54:60])
                if (not np.isfinite(occupancy) or not np.isfinite(pdb_occupancy)
                        or not 0 <= occupancy <= 1 or not 0 <= pdb_occupancy <= 1
                        or abs(occupancy - pdb_occupancy) > .0051):
                    raise ValueError
            except (TypeError, ValueError):
                raise AssemblyError(f"CIF/PDB alternate occupancy missing or disagrees: {key}") from None
        label = row["label_asym_id"]
        sources.setdefault(label, []).append(line)
        if row["label_entity_id"] in types:
            polymer[label] = {"author_chain": key[0], "entity_id": row["label_entity_id"],
                              "polymer_type": types[row["label_entity_id"]]}
    if seen != set(original):
        raise AssemblyError("source PDB contains atoms absent from the first CIF model")
    return sources, polymer, author_map


def build_context(cif_path: Path, raw_text: str, primary_chain: str,
                  assembly_id: str | None = None) -> AssemblyContext:
    """Resolve a deposited assembly and express every member in the target frame.

    A unique assembly containing the target is selected automatically. A sole
    author-supported definition takes precedence over software-only alternatives;
    other ambiguities require an explicit ID. Alternatives remain in metadata.
    A repeated primary-chain instance uses
    the identity operator when present, otherwise the first deposited instance;
    inverse anchoring keeps the original pocket/reference frame unchanged.
    All assembly members are retained, including native modified residues.
    """
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    cif_path = Path(cif_path)
    data = MMCIF2Dict(str(cif_path))
    operators = _operators(data)
    sources, polymer, author_map = _atom_sources(data, raw_text)
    targets = [label for label, row in polymer.items() if row["author_chain"] == primary_chain]
    if len(targets) != 1:
        raise AssemblyError(f"target chain {primary_chain} has {len(targets)} CIF polymer instances")
    primary_label = targets[0]
    assemblies = {}
    for row in _rows(data, "_pdbx_struct_assembly_gen",
                     ["assembly_id", "oper_expression", "asym_id_list"]):
        members = assemblies.setdefault(row["assembly_id"], [])
        labels = row["asym_id_list"].replace(" ", "").split(",")
        if any(label not in sources for label in labels):
            raise AssemblyError("assembly names an absent first-model label chain")
        for product in operator_products(row["oper_expression"], set(operators)):
            transform = np.eye(4)
            for identifier in product:
                transform = transform @ operators[identifier]
            members.extend((label, product, transform) for label in labels)
            if len(members) > 100000:
                raise AssemblyError("assembly exceeds 100000 instances; no members truncated")
    available = [identifier for identifier, members in assemblies.items()
                 if any(label == primary_label for label, _, _ in members)]
    details = dict(zip(mmcif._column(data, "_pdbx_struct_assembly.id"),
                       mmcif._column(data, "_pdbx_struct_assembly.details")))
    author_supported = [identifier for identifier in available
                        if details.get(identifier, "").startswith("author_")]
    selection_policy = "explicit_id"
    if assembly_id is None:
        if len(available) == 1:
            assembly_id = available[0]
            selection_policy = "sole_deposited_target_assembly"
        elif len(author_supported) == 1:
            assembly_id = author_supported[0]
            selection_policy = "sole_author_supported_definition_over_software_only_alternatives"
        else:
            raise AssemblyError(f"target belongs to {len(available)} assemblies {available}; specify assembly_id")
    if assembly_id not in available:
        raise AssemblyError(f"assembly {assembly_id} does not contain target chain {primary_chain}")
    members = assemblies[assembly_id]
    anchors = [member for member in members if member[0] == primary_label]
    anchor = next((member for member in anchors if np.allclose(member[2], np.eye(4), atol=1e-10, rtol=0)),
                  anchors[0])
    inverse = np.linalg.inv(anchor[2])
    # Distinct definitions of the same physical instance must not duplicate atoms.
    instances, used = [], set()
    for label, product, transform in members:
        relative = inverse @ transform
        if np.allclose(relative, np.eye(4), atol=1e-10, rtol=0):
            relative = np.eye(4)
        key = (label, _matrix_key(relative))
        if key in used:
            continue
        used.add(key)
        instances.append((label, product, relative))
    instances.sort(key=lambda member: not (member[0] == primary_label
                                           and np.allclose(member[2], np.eye(4), atol=1e-10, rtol=0)))
    chain_map, reserved = {}, {primary_chain}
    # Associated non-polymers share an author/transform instance with their polymer.
    for label, _, transform in instances:
        author = sources[label][0][21]
        key = (author, _matrix_key(transform))
        if key in chain_map:
            continue
        identity = np.allclose(transform, np.eye(4), atol=1e-10, rtol=0)
        if author == primary_chain and identity:
            chain_map[key] = primary_chain
        elif author not in reserved:
            chain_map[key] = author
            reserved.add(author)
        else:
            free = next((char for char in mmcif._PDB_CHAIN_ALPHABET if char not in reserved), None)
            if free is None:
                raise AssemblyError("assembly chain copies exceed the 62 PDB chain identifiers")
            chain_map[key] = free
            reserved.add(free)
    atom_records, declarations, sequence_records, instance_rows, serial_map = [], [], [], [], {}
    serial = 0
    source_declarations = modified.declarations(raw_text)
    # Converted PDBs may lack MODRES. Recover declarations from the actual CIF.
    for row in _rows(data, "_pdbx_struct_mod_residue", ["auth_asym_id", "auth_comp_id", "auth_seq_id",
                                                         "PDB_ins_code", "parent_comp_id"]):
        insertion = "" if row["PDB_ins_code"] in ("?", ".") else row["PDB_ins_code"]
        key = (author_map[row["auth_asym_id"]], row["auth_seq_id"] + insertion, row["auth_comp_id"])
        if row["parent_comp_id"] in modified.STANDARD_AMINO_ACIDS and key[2] != row["parent_comp_id"]:
            if key in source_declarations and source_declarations[key] != row["parent_comp_id"]:
                raise AssemblyError(f"CIF/PDB modified-residue parents disagree: {key}")
            source_declarations[key] = row["parent_comp_id"]
    for label, product, transform in instances:
        source = sources[label]
        author = source[0][21]
        chain_key = (author, _matrix_key(transform))
        destination = chain_map[chain_key]
        identity = np.allclose(transform, np.eye(4), atol=1e-10, rtol=0)
        for original in source:
            serial += 1
            if serial > 99999:
                raise AssemblyError("assembly exceeds 99999 PDB atom serials; no serial wrapping")
            line = original[:6] + f"{serial:5d}" + original[11:21] + destination + original[22:]
            if not identity:
                point = transform @ np.array([float(original[a:b]) for a, b in ((30, 38), (38, 46), (46, 54))] + [1])
                if np.any(point[:3] < -999.999) or np.any(point[:3] > 9999.999):
                    raise AssemblyError("transformed coordinates exceed PDB fields")
                line = line[:30] + "".join(f"{value:8.3f}" for value in point[:3]) + line[54:]
            atom_records.append(line)
            key = (_matrix_key(transform), int(original[6:11]))
            if key in serial_map:
                raise AssemblyError("source atom serial aliases an assembly connection")
            serial_map[key] = serial
        if label in polymer:
            for (cid, number, name), parent in source_declarations.items():
                if cid == author:
                    if len(name) > 3:
                        raise AssemblyError("native modified component name cannot fit PDB MODRES")
                    insertion = number[-1] if number[-1:].isalpha() else " "
                    resnum = int(number[:-1] if insertion != " " else number)
                    declarations.append(f"MODRES {cif_path.stem[:4]:4s} {name:3s} {destination} {resnum:4d}{insertion} {parent:3s}")
            sequence_records.extend(line[:11] + destination + line[12:] for line in raw_text.splitlines()
                                    if line.startswith("SEQRES") and line[11] == author)
        instance_rows.append({"label_asym_id": label, "author_chain": author,
                              "output_chain": destination, "operators": list(product),
                              "relative_transform": transform.tolist(), "observed_atoms": len(source),
                              "polymer_type": polymer.get(label, {}).get("polymer_type")})
    connections = []
    for transform_key in dict.fromkeys(key[1] for key in chain_map):
        for line in raw_text.splitlines():
            if not line.startswith("CONECT"):
                continue
            numbers = [int(line[i:i+5]) for i in range(6, len(line.rstrip()), 5) if line[i:i+5].strip()]
            if not numbers or (transform_key, numbers[0]) not in serial_map:
                continue
            retained = [serial_map[(transform_key, number)] for number in numbers
                        if (transform_key, number) in serial_map]
            if len(retained) > 1:
                connections.append("CONECT" + "".join(f"{number:5d}" for number in retained))
    # Keep each output chain contiguous even when CIF labels separate a native
    # modified residue from its polymer. Polymer HETATM records are sorted into
    # their author residue position; non-polymers follow their polymer.
    native_keys = {(line[16], line[18:23].strip(), line[12:15].strip()) for line in declarations}
    def record_order(line):
        polymer_atom = line.startswith("ATOM") or (line[21], line[22:27].strip(), line[17:20].strip()) in native_keys
        return (list(chain_map.values()).index(line[21]), not polymer_atom, int(line[22:26]), line[26], int(line[6:11]))
    atom_records.sort(key=record_order)
    metadata = {"schema_version": 1, "status": "DEPOSITED_ASSEMBLY_RESOLVED",
                "assembly_id": assembly_id, "available_target_assemblies": available,
                "selection_policy": selection_policy,
                "alternative_definitions": [{"assembly_id": identifier, "definition": details.get(identifier)}
                                             for identifier in available if identifier != assembly_id],
                "definition": details.get(assembly_id), "source_sha256": hashlib.sha256(cif_path.read_bytes()).hexdigest(),
                "original_pdb_sha256": hashlib.sha256(raw_text.encode()).hexdigest(),
                "primary_chain": primary_chain, "primary_label_asym_id": primary_label,
                "anchor_operators": list(anchor[1]), "instances": instance_rows,
                "observed_atoms": len(atom_records), "explicit_connection_records": len(connections),
                "scope": "Deposited assembly and coordinate preservation; physiological relevance, native chemistry and docking validity require separate checks"}
    text = "\n".join(sequence_records + declarations + atom_records + connections) + "\nEND\n"
    metadata["assembled_pdb_sha256"] = hashlib.sha256(text.encode()).hexdigest()
    return AssemblyContext(text, metadata)
