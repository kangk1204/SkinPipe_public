"""Resolve an explicitly requested human protein, without first-hit selection.

Identity is not evidence of a skin function, desired action, usable pocket or
compound activity. Keep the API responses and their hashes with the decision.
"""
from __future__ import annotations

import copy
import hashlib
import json
import re
from datetime import datetime, timezone
from typing import Any

from . import chembl
from .cache import FetchError, cached

SCHEMA_VERSION = "human-protein-identity-v1"
# Keep each provider's retry budget below the bridge's overall 90-second fence.
# A mapping outage must not erase a UniProt identity already verified above.
REQUEST_TIMEOUT = (5, 10)
REQUEST_ATTEMPTS = 2
_ACCESSION = re.compile(r"(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})\Z")
_GENE = re.compile(r"[A-Z][A-Z0-9-]{0,29}\Z")
_TARGET = re.compile(r"CHEMBL[0-9]{1,12}\Z")


class TargetIdentityError(ValueError):
    def __init__(self, message: str, *, status: str = "INVALID",
                 choices: list[dict[str, Any]] | None = None,
                 sources: list[dict[str, Any]] | None = None):
        super().__init__(message)
        self.status = status
        self.choices = choices or []
        self.sources = sources

    def to_dict(self) -> dict[str, Any]:
        result = {"status": self.status, "error": str(self), "choices": self.choices}
        if self.sources is not None:
            result["sources"] = self.sources
        return result


def validate_query(query: str) -> tuple[str, str]:
    if not isinstance(query, str):
        raise TargetIdentityError("Enter one gene symbol, primary UniProt accession or ChEMBL target ID.")
    query = query.strip().upper()
    if _TARGET.fullmatch(query):
        return query, "chembl_target"
    if _ACCESSION.fullmatch(query):
        return query, "uniprot"
    if _GENE.fullmatch(query) and not query.startswith("CHEMBL"):
        return query, "gene"
    raise TargetIdentityError("Use a single identifier (for example JAK1, P23458 or CHEMBL2835); no search expressions.")


def _request(url: str, params: dict, sources: list[dict], cache) -> dict:
    key = {"schema": SCHEMA_VERSION, "url": url, "params": params}

    def fetch():
        return {"data": chembl._get(url, params, retries=REQUEST_ATTEMPTS,
                                   timeout=REQUEST_TIMEOUT),
                "fetched_at": datetime.now(timezone.utc).isoformat()}

    policy = {"connect_timeout_seconds": REQUEST_TIMEOUT[0],
              "read_timeout_seconds": REQUEST_TIMEOUT[1], "max_attempts": REQUEST_ATTEMPTS}
    try:
        response = cached(chembl._c(cache), "target_identity", key, fetch)
    except FetchError as exc:
        sources.append({"url": url, "params": dict(params), "status": "UNAVAILABLE",
                        "attempted_at": datetime.now(timezone.utc).isoformat(),
                        "request_policy": policy, "error": str(exc)})
        raise
    data = response.value["data"]
    digest = hashlib.sha256(json.dumps(data, sort_keys=True, separators=(",", ":"),
                                     ensure_ascii=False).encode()).hexdigest()
    sources.append({"url": url, "params": dict(params),
                    "fetched_at": response.value["fetched_at"],
                    "from_cache": response.from_cache, "stale": response.stale,
                    "cache_age_days": response.age_days,
                    "request_policy": policy if not response.from_cache else None,
                    "raw_response_sha256": digest, "raw_response": copy.deepcopy(data)})
    return data


def _protein(row: dict) -> dict | None:
    accession = row.get("primaryAccession", "")
    reviewed = row.get("entryType") == "UniProtKB reviewed (Swiss-Prot)"
    if (not _ACCESSION.fullmatch(accession) or not reviewed
            or str((row.get("organism") or {}).get("taxonId")) != "9606"):
        return None
    names = sorted({g.get("geneName", {}).get("value") for g in row.get("genes", [])
                    if g.get("geneName", {}).get("value")})
    gene = names[0] if len(names) == 1 else None
    description = row.get("proteinDescription") or {}
    return {"uniprot": accession, "gene": gene, "gene_names": names,
            "protein_name": (description.get("recommendedName", {}).get("fullName", {}).get("value")),
            "organism": "Homo sapiens", "taxon_id": 9606, "reviewed": True,
            "sequence_length": row.get("sequence", {}).get("length")}


def _uniprot(query: str, kind: str, sources: list[dict], cache,
             primary_gene_only: bool = False) -> dict:
    # The search field includes aliases; enforce primary names on the response.
    gene_field = "gene_exact"
    clause = f"accession:{query}" if kind == "uniprot" else f"{gene_field}:{query}"
    data = _request(f"{chembl.UNIPROT}/search",
                    {"query": clause + " AND organism_id:9606 AND reviewed:true",
                     "format": "json", "size": 100}, sources, cache)
    rows = data.get("results")
    if not isinstance(rows, list):
        raise TargetIdentityError("UniProt returned no verifiable result list.", status="UNAVAILABLE")
    if len(rows) >= 100:
        raise TargetIdentityError("The gene query exceeds the verified result page; enter an exact UniProt accession.", status="AMBIGUOUS")
    choices = {}
    for row in rows:
        protein = _protein(row) if isinstance(row, dict) else None
        if protein is None:
            continue
        if kind == "uniprot" and protein["uniprot"] != query:
            continue
        if kind == "gene":
            gene_names = {n.upper() for n in protein["gene_names"]}
            if not primary_gene_only:
                gene_names.update(s.get("value", "").upper() for g in row.get("genes", [])
                                  for s in g.get("synonyms", []))
            elif gene_names != {query}:
                continue
            if query not in gene_names:
                continue
        choices[protein["uniprot"]] = protein
    ordered = [choices[k] for k in sorted(choices)]
    if not ordered:
        raise TargetIdentityError("No reviewed human protein verified for this identifier. Use a current primary UniProt accession.", status="NOT_FOUND")
    if len(ordered) != 1:
        raise TargetIdentityError("Multiple reviewed human proteins match; select an exact UniProt accession.",
                                  status="AMBIGUOUS", choices=ordered)
    return ordered[0]


def _target(row: dict, accession: str | None = None) -> dict | None:
    tid = row.get("target_chembl_id", "")
    components = row.get("target_components") or []
    if (not _TARGET.fullmatch(tid) or row.get("target_type") != "SINGLE PROTEIN"
            or row.get("organism") != "Homo sapiens"
            or str(row.get("tax_id")) != "9606" or len(components) != 1):
        return None
    mapped = components[0].get("accession", "")
    if not _ACCESSION.fullmatch(mapped) or (accession is not None and mapped != accession):
        return None
    return {"chembl_target": tid, "uniprot": mapped,
            "target_name": row.get("pref_name"), "target_type": "SINGLE PROTEIN"}


def resolve(query: str, *, cache=None, primary_gene_only: bool = False) -> dict:
    """Resolve identity only. Unknown counts and biological goals stay unknown."""
    query, kind = validate_query(query)
    if primary_gene_only:
        # Some real primary genes (P2RY12, B3GAT1) also match UniProt's
        # accession grammar. An explicitly gene-only caller must keep its
        # stated namespace rather than reject or reinterpret the symbol.
        if not _GENE.fullmatch(query) or kind == "chembl_target":
            raise TargetIdentityError("Primary-gene matching requires a gene symbol.")
        kind = "gene"
    sources: list[dict] = []
    target = None
    try:
        if kind == "chembl_target":
            row = _request(f"{chembl.CHEMBL}/target/{query}.json", {}, sources, cache)
            target = _target(row)
            if target is None or target["chembl_target"] != query:
                raise TargetIdentityError("This ChEMBL record is not a verified human single-protein target with one UniProt component.", status="UNSUPPORTED")
            protein = _uniprot(target["uniprot"], "uniprot", sources, cache)
        else:
            try:
                protein = _uniprot(query, kind, sources, cache, primary_gene_only)
            except TargetIdentityError as exc:
                # A successful empty accession lookup is not proof that the
                # same string cannot be a primary gene. Retry only this
                # namespace collision; outages and ambiguous identities keep
                # their original failure and never become biological absence.
                if kind != "uniprot" or exc.status != "NOT_FOUND":
                    raise
                protein = _uniprot(query, "gene", sources, cache, True)
                kind = "gene"
    except FetchError as exc:
        raise TargetIdentityError(f"Protein identity could not be verified: {exc}",
                                  status="UNAVAILABLE", sources=sources) from exc

    mapping_status = "VERIFIED" if target else "NOT_FOUND"
    mapping_note = ""
    if target is None:
        try:
            data = _request(f"{chembl.CHEMBL}/target.json",
                            {"target_components__accession": protein["uniprot"],
                             "target_type": "SINGLE PROTEIN", "limit": 100}, sources, cache)
            rows = data.get("targets")
            if not isinstance(rows, list):
                raise TargetIdentityError("ChEMBL returned no verifiable target list.", status="UNAVAILABLE")
            page = data.get("page_meta") or {}
            if page.get("next") or int(page.get("total_count", len(rows))) > len(rows):
                raise TargetIdentityError("ChEMBL target mappings are incomplete; select an exact target ID.", status="AMBIGUOUS")
            choices = {_t["chembl_target"]: _t for row in rows
                       if (_t := _target(row, protein["uniprot"])) is not None}
            ordered = [choices[k] for k in sorted(choices)]
            if primary_gene_only and any(not isinstance(row, dict) or
                    _target(row, protein["uniprot"]) is None for row in rows):
                raise TargetIdentityError("ChEMBL returned an unexpected component, species or target type; absence has not been established.", status="UNAVAILABLE")
            if len(ordered) > 1:
                raise TargetIdentityError("Multiple human ChEMBL targets map to this accession; select an exact ChEMBL target ID.",
                                          status="AMBIGUOUS", choices=ordered)
            if ordered:
                target = ordered[0]
                mapping_status = "VERIFIED"
        except FetchError as exc:
            mapping_status = "NOT_EVALUATED"
            mapping_note = f"ChEMBL mapping lookup failed: {exc}; user-supplied structures can still be explored."
    return {"schema": SCHEMA_VERSION, "status": "VERIFIED", "query": query,
            "query_kind": kind, "gene_match_policy": "primary_only" if primary_gene_only else "primary_or_synonym", **protein,
            "chembl_target": target["chembl_target"] if target else None,
            "chembl_mapping_status": mapping_status, "chembl_mapping_note": mapping_note,
            "target_type": target["target_type"] if target else None,
            "sources": sources, "biological_goal_status": "NOT_EVALUATED",
            "scope": "Reviewed human single-protein identity; no skin-function or activity claim."}


def main() -> int:
    """Small identity-only entry point, independent of the analysis CLI runtime."""
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("query")
    query = parser.parse_args().query
    try:
        result = resolve(query)
    except TargetIdentityError as exc:
        print(json.dumps(exc.to_dict(), ensure_ascii=False, allow_nan=False))
        return 2
    print(json.dumps(result, ensure_ascii=False, allow_nan=False))
    return 0


def validate_primary_binding(identity: dict, gene: str) -> bool:
    """Recheck a legacy node's exact binding against its retained API bytes.

    This validates identity, not the age or scientific meaning of node counts.
    A lookup that failed to establish a ChEMBL mapping cannot certify zero.
    """
    try:
        if not isinstance(identity, dict) or not _GENE.fullmatch(gene):
            return False
        expected = {"schema": SCHEMA_VERSION, "status": "VERIFIED", "query": gene,
                    "query_kind": "gene", "gene_match_policy": "primary_only",
                    "gene": gene, "gene_names": [gene], "reviewed": True,
                    "organism": "Homo sapiens", "taxon_id": 9606}
        if any(identity.get(k) != v for k, v in expected.items()):
            return False
        sources = identity.get("sources")
        if not isinstance(sources, list) or len(sources) != 2:
            return False
        for source in sources:
            digest = hashlib.sha256(json.dumps(source["raw_response"], sort_keys=True,
                separators=(",", ":"), ensure_ascii=False).encode()).hexdigest()
            if digest != source.get("raw_response_sha256") or not source.get("fetched_at"):
                return False
        up, ch = sources
        if up.get("url") != chembl.UNIPROT + "/search" or up.get("params") != {
                "query": f"gene_exact:{gene} AND organism_id:9606 AND reviewed:true",
                "format": "json", "size": 100}:
            return False
        rows = up["raw_response"].get("results")
        if not isinstance(rows, list) or len(rows) >= 100:
            return False
        proteins = {p["uniprot"]: p for row in rows if isinstance(row, dict)
                    and (p := _protein(row)) is not None and p["gene_names"] == [gene]}
        if len(proteins) != 1 or identity.get("uniprot") not in proteins:
            return False
        protein = proteins[identity["uniprot"]]
        if any(identity.get(k) != v for k, v in protein.items()):
            return False
        if ch.get("url") != chembl.CHEMBL + "/target.json" or ch.get("params") != {
                "target_components__accession": protein["uniprot"],
                "target_type": "SINGLE PROTEIN", "limit": 100}:
            return False
        data = ch["raw_response"]
        rows = data.get("targets")
        page = data.get("page_meta") or {}
        if (not isinstance(rows, list) or page.get("next") or
                int(page.get("total_count", len(rows))) != len(rows)):
            return False
        mapped = [_target(row, protein["uniprot"]) if isinstance(row, dict) else None for row in rows]
        if any(t is None for t in mapped):
            return False
        targets = {t["chembl_target"] for t in mapped}
        if not targets:
            return (identity.get("chembl_mapping_status") == "NOT_FOUND" and
                    identity.get("chembl_target") is None and identity.get("target_type") is None)
        return (len(targets) == 1 and identity.get("chembl_target") in targets and
                identity.get("chembl_mapping_status") == "VERIFIED" and
                identity.get("target_type") == "SINGLE PROTEIN")
    except (KeyError, TypeError, ValueError, AttributeError):
        return False


if __name__ == "__main__":
    raise SystemExit(main())
