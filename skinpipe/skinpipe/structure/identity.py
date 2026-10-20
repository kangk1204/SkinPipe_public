"""Deposited, observed sequence segments behind a requested protein identity.

A fusion has several reference proteins on one chain. Entity-level API labels
can omit a component, and an antigen peptide can have a correct accession
without containing the protein's functional domain. Keep those claims separate.
"""
from __future__ import annotations

import hashlib
from pathlib import Path
import re

from . import mmcif, rcsb


class TargetStructureError(ValueError):
    """The construct requires an explicit peptide target role."""


def _rows(data: dict, category: str) -> list[dict]:
    columns = {key.split(".", 1)[1]: mmcif._column(data, key)
               for key in data if key.startswith(category + ".")}
    if not columns:
        return []
    lengths = {len(value) for value in columns.values()}
    if len(lengths) != 1:
        raise ValueError(f"unequal {category} identity columns")
    return [dict(zip(columns, values)) for values in zip(*columns.values())]


def source_entities(cif_path: Path, api_entities: list[dict]) -> list[dict]:
    """Use source references only when their segment has observed polymer atoms.

    The source is subsequently checked against the coordinate file by the
    assembly resolver. Missing source references leave API-only evidence;
    disagreeing source/API accessions are retained as an explicit discrepancy.
    """
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict

    data = MMCIF2Dict(str(cif_path))
    atoms = _rows(data, "_atom_site")
    first_model = atoms[0].get("pdbx_PDB_model_num") if atoms else None
    atoms = [row for row in atoms if row.get("pdbx_PDB_model_num") == first_model]
    authors = list(dict.fromkeys(row.get("auth_asym_id", "") for row in atoms))
    alphabet = mmcif._PDB_CHAIN_ALPHABET
    if len(authors) > len(alphabet):
        raise mmcif.StructureTooLargeError("source identity chains cannot fit in PDB format")
    reserved = {cid for cid in authors if cid in alphabet and len(cid) == 1}
    spare = iter(char for char in alphabet if char not in reserved)
    chain_map = {cid: cid if cid in reserved else next(spare) for cid in authors}
    observed: dict[tuple[str, str], set[int]] = {}
    for row in atoms:
        if row.get("label_alt_id", ".") not in (".", "?", "A"):
            continue
        position = row.get("label_seq_id", ".")
        if not position.isdigit():
            continue
        observed.setdefault((row.get("label_entity_id", ""), row.get("auth_asym_id", "")), set()).add(int(position))
    ref_rows = _rows(data, "_struct_ref")
    references = {row.get("id"): row for row in ref_rows}
    if len(references) != len(ref_rows):
        raise ValueError("duplicate deposited reference identifiers")
    segments = _rows(data, "_struct_ref_seq")
    descriptions = {row.get("id"): row.get("pdbx_description", "")
                    for row in _rows(data, "_entity")}
    source_sha = hashlib.sha256(Path(cif_path).read_bytes()).hexdigest()
    out = []
    for polymer in _rows(data, "_entity_poly"):
        eid = polymer.get("entity_id", "")
        source_chains = [cid for entity, cid in observed if entity == eid]
        chains = [chain_map[cid] for cid in source_chains]
        if not chains:
            continue
        api = next((row for row in api_entities if set(row.get("chains") or []) & set(chains)), {})
        value = dict(api)
        sequence = re.sub(r"\s+", "", polymer.get("pdbx_seq_one_letter_code_can", ""))
        if sequence in ("", ".", "?"):
            sequence = api.get("sequence", "")
        value.update(chains=chains, description=descriptions.get(eid) or api.get("description", ""),
                     sequence=sequence,
                     source_cif_sha256=source_sha, source_entity_id=eid, reference_segments=[])
        for segment in segments:
            ref = references.get(segment.get("ref_id"), {})
            if ref.get("entity_id") != eid or ref.get("db_name", "").upper() not in ("UNP", "UNIPROT"):
                continue
            accession = ref.get("pdbx_db_accession", "").strip().upper()
            if accession in ("", ".", "?"):
                continue
            try:
                lo, hi = int(segment["seq_align_beg"]), int(segment["seq_align_end"])
                db_lo, db_hi = int(segment["db_align_beg"]), int(segment["db_align_end"])
            except (KeyError, ValueError) as exc:
                raise ValueError("invalid deposited reference alignment") from exc
            if lo < 1 or hi < lo or db_lo < 1 or db_hi < db_lo:
                raise ValueError("invalid deposited reference alignment range")
            for cid in segment.get("pdbx_strand_id", "").split(","):
                cid = cid.strip()
                positions = sorted(i for i in observed.get((eid, cid), ()) if lo <= i <= hi)
                if positions:
                    value["reference_segments"].append({"accession": accession, "chain": chain_map[cid],
                        "entity_begin": lo, "entity_end": hi, "reference_begin": db_lo,
                        "reference_end": db_hi, "observed_residues": len(positions)})
        source_acc = list(dict.fromkeys(s["accession"] for s in value["reference_segments"]))
        if source_acc:
            api_acc = set(api.get("uniprots") or []) | ({api["uniprot"]} if api.get("uniprot") else set())
            value.update(uniprots=source_acc, uniprot=source_acc[0],
                         api_accession_discrepancy=sorted(api_acc - set(source_acc)),
                         identity_source="deposited_observed_reference_segments")
        else:
            value["identity_source"] = "api_only_no_observed_source_reference"
        out.append(value)
    return out or api_entities


def bind_gene(entities: list[dict], gene: str | None, uniprot: str | None) -> tuple[list[dict], list[str]]:
    """Resolve absent primary gene labels through exact reviewed human records.

    Only accessions already mapped to observed source chains are queried. No
    synonym, substring or description is accepted as a primary gene match.
    Explicit accessions remain authoritative over a coincident gene label.
    """
    wanted = (gene or "").strip().upper()
    if not wanted or uniprot or any(str(row.get("gene") or "").upper() == wanted
                                  and len(row.get("uniprots") or []) <= 1 for row in entities):
        return entities, []
    checked, out, notes = {}, [], []
    for entity in entities:
        value = dict(entity)
        for accession in entity.get("uniprots") or ([entity["uniprot"]] if entity.get("uniprot") else []):
            if accession not in checked:
                try:
                    checked[accession] = rcsb.uniprot_identity(accession)
                except rcsb.RCSBError as exc:
                    checked[accession] = None
                    notes.append(f"primary gene evidence unavailable for {accession}: {exc}")
            record = checked[accession]
            if record and record["reviewed"] and record["taxon_id"] == 9606 and wanted in record["primary_genes"]:
                value.update(gene=wanted, uniprot=accession, primary_identity=record)
                notes.append(f"primary human gene {wanted} bound to exact UniProt {accession}; source chain mapping retained")
                break
        out.append(value)
    return out, notes


def assess(entity: dict | None, chain: str, accession: str | None, role: str = "protein") -> dict:
    """Require explicit peptide admission for short target constructs.

    A protein accession can annotate a substrate or antigen peptide without
    saying "peptide" in the molecule name. Source mapping measures the requested
    component of a fusion, rather than the length of its larger partner. Domain
    coverage against a full reference protein is not an admission threshold.
    """
    if role not in ("protein", "peptide"):
        raise ValueError("target_role must be protein or peptide")
    entity = entity or {}
    accession = accession or entity.get("uniprot")
    segments = [s for s in entity.get("reference_segments", [])
                if s["chain"] == chain and s["accession"] == accession]
    sequence = re.sub(r"\s+", "", entity.get("sequence") or "")
    length = len(sequence) if sequence not in (".", "?") else 0
    intervals = sorted((s["entity_begin"], s["entity_end"]) for s in segments)
    merged = []
    for start, end in intervals:
        if merged and start <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    target_length = sum(end - start + 1 for start, end in merged) if merged else length
    annotated_peptide = bool(re.search(r"\bpeptide\b", entity.get("description", ""), re.I))
    short_peptide = annotated_peptide and 0 < length <= 30
    short_construct = 0 < target_length <= 30
    result = {"policy": "observed_reference_target_role_v2", "target_role": role,
              "accession": accession, "chain": chain, "construct_sequence_length": length,
              "reference_segments": segments, "annotated_short_peptide": short_peptide,
              "target_construct_residues": target_length,
              "target_length_evidence": "deposited_reference_segments" if merged else "entity_sequence",
              "short_target_construct": short_construct,
              "functional_admissibility": "NOT_ESTABLISHED"}
    if role == "protein" and (short_peptide or short_construct):
        raise TargetStructureError(f"chain {chain} contains a {target_length}-residue peptide-sized construct of "
            f"{accession or 'an unspecified protein'}; it does not establish a functional protein receptor. "
            "Select a structure containing the intended domain, or explicitly request target_role='peptide'.")
    return result
