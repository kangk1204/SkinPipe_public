"""Thin, timeout-bounded HTTP layer over RCSB and UniProt.

Kept separate from :mod:`skinpipe.structure.prepare` so the structure-preparation
logic can be unit-tested without a network, and so every outbound call lives in
one auditable place.

Endpoints used
--------------
``files.rcsb.org/download/{id}.pdb``
    Legacy PDB coordinates. Absent for large assemblies (ribosomes, virus
    capsids, anything over 62 chains or 99 999 atoms), hence the mmCIF fallback.
``files.rcsb.org/download/{id}.cif``
    mmCIF, always available. Converted to PDB locally with Biopython.
``files.rcsb.org/ligands/download/{CODE}_ideal.sdf``
    Idealised chemical-component coordinates *with bond orders*. Crystallographic
    HETATM records carry no bond orders, so this file is the template used to
    perceive the cocrystal ligand's chemistry.
``search.rcsb.org/rcsbsearch/v2/query``
    Structure search; used to enumerate entries mapped to a UniProt accession.
``data.rcsb.org/graphql``
    Batched entry metadata (resolution, method, bound components) in one round
    trip rather than one REST call per PDB ID.
``rest.uniprot.org/uniprotkb/{acc}.json``
    Binding-site / active-site / metal-binding feature annotations, used as the
    fallback pocket definition for apo structures.
"""
from __future__ import annotations

import logging
import time
from pathlib import Path
from typing import Any, Optional

import requests

log = logging.getLogger(__name__)

#: Coordinate downloads. A big mmCIF legitimately takes tens of seconds.
HTTP_TIMEOUT: float = 60.0

#: Search, GraphQL and UniProt. These answer in about a second when healthy,
#: so a minute of waiting on a hung one buys nothing and costs a minute --
#: with three attempts that is the difference between 63 s and 180 s.
QUERY_TIMEOUT: float = 20.0






RETRIES: int = 3

#: Statuses worth trying again, same set as `fetch.chembl._get`. A 404 is not
#: here on purpose: `fetch_ligand_sdf_text` *depends* on it meaning "RCSB has
#: no ideal SDF for this component", and retrying it would spend three
#: timeouts per missing het code to learn what the first answer already said.
RETRY_STATUS = frozenset({429, 500, 502, 503, 504})

_USER_AGENT = "skinpipe/0.1 (structure preparation)"

_FILES = "https://files.rcsb.org"
_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
_GRAPHQL_URL = "https://data.rcsb.org/graphql"
_UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/{acc}.json"

#: Metals that are genuinely part of a catalytic or structural site and must
#: survive receptor cleaning. Tyrosinase is a type-3 di-copper enzyme; the MMPs
#: carry a catalytic Zn plus a structural Zn and several Ca. Dropping these
#: silently turns a metalloenzyme pocket into a hole and the docking scores
#: become meaningless.
#:
#: The oxidation-state variants are not optional extras. The PDB gives each one
#: its own chemical-component id -- ``CU1`` is copper(I), ``FE2`` is iron(II),
#: ``MN3`` manganese(III), ``3CO`` cobalt(III) -- and a depositor modelling a
#: reduced metal uses them instead of the bare element code. Matching on the
#: element codes alone sends those atoms down the "every other HETATM is a
#: buffer component" branch of
#: :func:`skinpipe.structure.clean.select_receptor_lines`, so the metal is
#: deleted without a word. Over 1600 experimental entries use one of these.
# Historical name-list policy: membership alone does not establish a
# catalytic/biological role; preparation inventories each retained instance.
CATALYTIC_METALS: frozenset[str] = frozenset(
    {"ZN", "CU", "FE", "MG", "MN", "CA", "NI", "CO",
     "CU1", "CU3", "FE2", "MN3", "3CO", "3NI"}
)















_GLYCANS: frozenset[str] = frozenset({
    "NAG", "NDG",            # N-acetylglucosamine, beta / alpha
    "GLC", "BGC",            # glucose, alpha / beta
    "MAN", "BMA",            # mannose, alpha / beta
    "GAL", "GLA",            # galactose, beta / alpha
    "A2G",                   # N-acetylgalactosamine, alpha
    "FUC", "FUL",            # fucose, alpha / beta
    "XYS", "XYP",            # xylose, alpha / beta
    "SIA", "SLB",            # sialic acid (N-acetylneuraminate)
    "RAM",                   # rhamnose
    "FRU",                   # fructose
    "SUC", "TRE",            # sucrose, trehalose -- cryoprotectants
})

HET_BLACKLIST: frozenset[str] = frozenset(
    {
        "HOH", "DOD", "SO4", "GOL", "EDO", "PEG", "PGE", "PG4", "1PE", "ACT",
        "CL", "NA", "K", "BR", "IOD", "DMS", "TRS", "MES", "PO4", "FMT", "ACY",
        "MPD", "CIT", "TLA", "NO3", "CO3", "BME", "IMD", "EPE", "SCN", "AZI",
        "FLC", "OXL", "MLI", "UNL", "UNX",
        # Related additives in the wwPDB CCD: MPD enantiomer, tartrate,
        # malonate, bicarbonate and additional ethylene-glycol oligomers.
        "MRD", "TAR", "MLA", "BCT", "P6G", "2PE",
    }
    | _GLYCANS
) | CATALYTIC_METALS


# Public identity list for reporting excluded observed glycans. Membership is
# not a claim that removing a glycan preserves pocket chemistry.
GLYCANS = _GLYCANS


class RCSBError(RuntimeError):
    """Any non-recoverable failure talking to RCSB or UniProt."""


def _get(url: str, timeout: Optional[float] = None,
         **kwargs: Any) -> requests.Response:
    """GET with a hard timeout, retries, and a greppable failure message.

    Returns the response whatever its status: callers read 404 as "RCSB does
    not have this", and that is an answer, not a failure. Only transport
    errors and :data:`RETRY_STATUS` are tried again.
    """
    last: Optional[str] = None
    for attempt in range(RETRIES):
        try:
            resp = requests.get(
                url, timeout=timeout or HTTP_TIMEOUT,
                headers={"User-Agent": _USER_AGENT}, **kwargs)
            if resp.status_code not in RETRY_STATUS:
                return resp
            last = f"HTTP {resp.status_code}"
        except requests.RequestException as exc:  # DNS, TLS, read timeout
            last = f"{type(exc).__name__}: {exc}"
        if attempt + 1 < RETRIES:
            log.warning("GET %s failed (%s); retry %d/%d",
                        url, last, attempt + 2, RETRIES)
            time.sleep(2 ** attempt)
    raise RCSBError(f"GET {url} failed after {RETRIES} attempts; "
                    f"last failure: {last}")


def _post_json(url: str, payload: dict[str, Any]) -> dict[str, Any]:
    last: Optional[str] = None
    for attempt in range(RETRIES):
        try:
            resp = requests.post(
                url,
                json=payload,
                timeout=QUERY_TIMEOUT,
                headers={"User-Agent": _USER_AGENT,
                         "Content-Type": "application/json"},
            )
            if resp.status_code == 204:  # search API says "no hits" this way
                return {}
            if resp.status_code not in RETRY_STATUS:
                if not resp.ok:
                    raise RCSBError(
                        f"POST {url} -> HTTP {resp.status_code}: {resp.text[:200]}")
                data = resp.json()
                if not isinstance(data, dict):
                    raise RCSBError(f"POST {url}: response is not a JSON object")
                if data.get("errors"):
                    raise RCSBError(f"POST {url}: GraphQL errors: {str(data['errors'])[:200]}")
                if url == _SEARCH_URL and not data:
                    raise RCSBError(f"POST {url}: HTTP 200 returned no search result metadata")
                return data
            last = f"HTTP {resp.status_code}"
        except requests.RequestException as exc:
            last = f"{type(exc).__name__}: {exc}"
        if attempt + 1 < RETRIES:
            log.warning("POST %s failed (%s); retry %d/%d",
                        url, last, attempt + 2, RETRIES)
            time.sleep(2 ** attempt)
    raise RCSBError(f"POST {url} failed after {RETRIES} attempts; "
                    f"last failure: {last}")


def download_coordinates(pdb_id: str, dest_dir: Path) -> tuple[Path, str]:
    """Download ``{pdb_id}`` coordinates, preferring PDB format over mmCIF.

    Returns ``(path, fmt)`` where ``fmt`` is ``"pdb"`` or ``"cif"``. Both formats
    are cached on disk; a non-empty cache hit short-circuits the network.
    """
    code = pdb_id.strip().upper()
    if not code:
        raise ValueError("pdb_id must not be empty")
    dest_dir.mkdir(parents=True, exist_ok=True)

    pdb_path = dest_dir / f"{code}.pdb"
    cif_path = dest_dir / f"{code}.cif"
    if pdb_path.is_file() and pdb_path.stat().st_size > 0:
        return pdb_path, "pdb"
    if cif_path.is_file() and cif_path.stat().st_size > 0:
        return cif_path, "cif"

    resp = _get(f"{_FILES}/download/{code}.pdb")
    if resp.ok and resp.content:
        pdb_path.write_bytes(resp.content)
        log.debug("fetched %s as PDB (%d bytes)", code, len(resp.content))
        return pdb_path, "pdb"

    log.info(
        "no PDB-format file for %s (HTTP %d) - falling back to mmCIF",
        code,
        resp.status_code,
    )
    resp = _get(f"{_FILES}/download/{code}.cif")
    if not resp.ok or not resp.content:
        raise RCSBError(
            f"{code}: neither PDB nor mmCIF available (mmCIF HTTP {resp.status_code})"
        )
    cif_path.write_bytes(resp.content)
    return cif_path, "cif"


def download_mmcif(pdb_id: str, dest_dir: Path) -> Path:
    """Cache the authoritative assembly definition even when a PDB is present.

    Coordinate agreement and assembly/operator validity are checked by the
    assembly resolver. A failed request cannot be reinterpreted as a monomer.
    """
    code = pdb_id.strip().upper()
    if not code or any(char not in "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_" for char in code):
        raise ValueError("invalid PDB identifier")
    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)
    path = dest_dir / f"{code}.cif"
    if path.is_file() and path.stat().st_size:
        return path
    response = _get(f"{_FILES}/download/{code}.cif")
    if not response.ok or not response.content:
        raise RCSBError(f"{code}: assembly mmCIF unavailable (HTTP {response.status_code})")
    temporary = path.with_suffix(".cif.tmp")
    temporary.write_bytes(response.content)
    temporary.replace(path)
    return path


def fetch_ligand_sdf_text(het_code: str, cache_dir: Optional[Path] = None) -> Optional[str]:
    """Idealised SDF for a chemical component, or ``None`` if RCSB has none."""
    code = het_code.strip().upper()
    cached = (cache_dir / f"{code}_ideal.sdf") if cache_dir else None
    if cached is not None and cached.is_file() and cached.stat().st_size > 0:
        return cached.read_text()

    resp = _get(f"{_FILES}/ligands/download/{code}_ideal.sdf")
    if not resp.ok or not resp.text.strip():
        log.warning("no ideal SDF for chemical component %s (HTTP %d)", code, resp.status_code)
        return None
    if cached is not None:
        cached.parent.mkdir(parents=True, exist_ok=True)
        cached.write_text(resp.text)
    return resp.text


def fetch_ligand_cif_text(het_code: str, cache_dir: Optional[Path] = None) -> Optional[str]:
    """CCD component definition for a chemical component, or None.

    Wanted for `_chem_comp_atom.atom_id` and `_chem_comp_bond`: the atom names
    and the bond table. Those let an observed ligand be rebuilt by *name*
    instead of by perceived geometry, which is the only route that survives a
    partially occupied ligand -- see `ligand.mol_from_ccd_names`.
    """
    code = het_code.strip().upper()
    cached = (cache_dir / f"{code}.cif") if cache_dir else None
    if cached is not None and cached.is_file() and cached.stat().st_size > 0:
        return cached.read_text()

    resp = _get(f"{_FILES}/ligands/download/{code}.cif")
    if not resp.ok or not resp.text.strip():
        log.warning("no CCD cif for chemical component %s (HTTP %d)", code, resp.status_code)
        return None
    if cached is not None:
        cached.parent.mkdir(parents=True, exist_ok=True)
        cached.write_text(resp.text)
    return resp.text


class StructureSearchResult(list):
    """List-compatible identifiers with explicit query coverage."""

    def __init__(self, ids: list[str], *, total_count: int, limit: int):
        super().__init__(ids)
        self.source_scope = {"total_count": total_count, "retrieved": len(ids),
                             "max_records": limit, "complete": len(ids) == total_count,
                             "capped_at_max_records": bool(limit and len(ids) < total_count)}


def entries_for_uniprot(uniprot: str, limit: int = 0) -> list[str]:
    """Mapped PDB IDs; zero retrieves all hits, a positive limit is explicit.

    The selector ranks the returned structures itself. A hidden first50 limit
    therefore excluded structures before the docking suitability ranking ran.
    RCSB's documented ``return_all_hits`` retains the whole matched set in one
    response, whose count and unique identifiers are checked before returning.
    """
    if type(limit) is not int or limit < 0 or limit > 10000:
        raise ValueError("structure limit must be an integer from 0 to 10000")
    acc = uniprot.strip().upper()
    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_polymer_entity_container_identifiers"
                        ".reference_sequence_identifiers.database_accession",
                        "operator": "in",
                        "value": [acc],
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_polymer_entity_container_identifiers"
                        ".reference_sequence_identifiers.database_name",
                        "operator": "exact_match",
                        "value": "UniProt",
                    },
                },
            ],
        },
        "return_type": "entry",
        "request_options": {
            "results_content_type": ["experimental", "computational"],
        },
    }
    if limit:
        query["request_options"]["paginate"] = {"start": 0, "rows": limit}
    else:
        query["request_options"]["return_all_hits"] = True
    data = _post_json(_SEARCH_URL, query)
    if not data:  # The HTTP204 path is the authoritative empty query.
        return StructureSearchResult([], total_count=0, limit=limit)
    total = data.get("total_count")
    hits = data.get("result_set", [])
    if type(total) is not int or total < 0 or not isinstance(hits, list):
        raise RCSBError(f"UniProt {acc}: invalid structure search count or result set")
    ids = []
    for hit in hits:
        identifier = hit.get("identifier") if isinstance(hit, dict) else None
        if not isinstance(identifier, str) or not identifier.strip():
            raise RCSBError(f"UniProt {acc}: structure hit has no identifier")
        ids.append(identifier.strip().upper())
    expected = min(total, limit) if limit else total
    if len(ids) != expected or len(set(ids)) != len(ids):
        raise RCSBError(f"UniProt {acc}: incomplete or duplicated structure query "
                        f"({len(ids)} identifiers, expected {expected})")
    return StructureSearchResult(ids, total_count=total, limit=limit)


_ENTRY_GQL = """
query($ids:[String!]!){
  entries(entry_ids:$ids){
    rcsb_id
    struct { title }
    rcsb_entry_info { resolution_combined experimental_method }
    nonpolymer_entities {
      rcsb_nonpolymer_entity_container_identifiers { auth_asym_ids }
      nonpolymer_comp {
        chem_comp { id name formula_weight }
        rcsb_chem_comp_info { atom_count_heavy }
      }
    }
  }
}
"""


def entry_metadata(pdb_ids: list[str]) -> dict[str, dict[str, Any]]:
    """Batched entry metadata keyed by upper-case PDB ID.

    Computed-structure-model identifiers (``AF_...``, ``MA_...``) are not served
    by the ``entries`` endpoint; they are returned with ``method="predicted"``
    and no ligands so the ranker can still see and deprioritise them.
    """
    ids = list(dict.fromkeys(i.upper() for i in pdb_ids))
    experimental = [i for i in ids if len(i) == 4]
    out: dict[str, dict[str, Any]] = {
        i: {
            "pdb_id": i,
            "title": "",
            "resolution": None,
            "method": "predicted",
            "ligands": [],
        }
        for i in ids
        if i not in experimental
    }
    if not experimental:
        return out

    for offset in range(0, len(experimental), 500):
        chunk = experimental[offset:offset + 500]
        data = _post_json(_GRAPHQL_URL, {"query": _ENTRY_GQL, "variables": {"ids": chunk}})
        entries = (data.get("data") or {}).get("entries")
        if data.get("errors") or not isinstance(entries, list):
            raise RCSBError(f"Entry metadata batch {offset // 500 + 1} failed")
        received = [entry["rcsb_id"].upper() for entry in entries
                    if isinstance(entry, dict) and isinstance(entry.get("rcsb_id"), str)]
        if len(entries) != len(chunk) or len(received) != len(chunk) or set(received) != set(chunk):
            raise RCSBError(f"Entry metadata batch {offset // 500 + 1} incomplete or duplicated")
        for entry in entries:
            out[entry["rcsb_id"].upper()] = _entry_metadata_record(entry)
    return out


def _entry_metadata_record(entry):
    """Decode one validated entry without changing ligand or chain semantics."""
    info = entry.get("rcsb_entry_info") or {}
    res = info.get("resolution_combined") or []
    ligands = []
    for npe in entry.get("nonpolymer_entities") or []:
        comp = ((npe or {}).get("nonpolymer_comp") or {}).get("chem_comp") or {}
        if comp.get("id"):

            ids = ((npe or {}).get(
                "rcsb_nonpolymer_entity_container_identifiers") or {})
            chains = [str(c) for c in (ids.get("auth_asym_ids") or [])]
            ligands.append(
                {
                    "het_code": comp["id"].upper(),
                    "name": comp.get("name") or "",
                    "mw": comp.get("formula_weight"),
                    "n_heavy_atoms": (((npe or {}).get("nonpolymer_comp") or {})
                                      .get("rcsb_chem_comp_info") or {}).get("atom_count_heavy"),
                    "chains": chains,
                }
            )
    return {
        "pdb_id": entry["rcsb_id"].upper(),
        "title": (entry.get("struct") or {}).get("title") or "",
        "resolution": float(res[0]) if res else None,
        "method": info.get("experimental_method") or "unknown",
        "ligands": ligands,
    }

_ENTITY_GQL = """
query($ids:[String!]!){
  entries(entry_ids:$ids){
    rcsb_id
    polymer_entities {
      rcsb_polymer_entity_container_identifiers {
        auth_asym_ids
        reference_sequence_identifiers { database_accession database_name }
      }
      rcsb_polymer_entity { pdbx_description }
      rcsb_entity_source_organism { rcsb_gene_name { value } }
      entity_poly { pdbx_seq_one_letter_code_can }
    }
  }
}
"""


def polymer_entities(pdb_id: str) -> list[dict[str, Any]]:
    """Per-entity chain IDs, UniProt accession, gene name and canonical sequence."""
    code = pdb_id.strip().upper()
    data = _post_json(_GRAPHQL_URL, {"query": _ENTITY_GQL, "variables": {"ids": [code]}})
    entries = (data.get("data") or {}).get("entries") or []
    if not entries:
        return []
    out: list[dict[str, Any]] = []
    for ent in entries[0].get("polymer_entities") or []:
        ids = ent.get("rcsb_polymer_entity_container_identifiers") or {}
        refs = ids.get("reference_sequence_identifiers") or []
        uniprots = list(dict.fromkeys(
            str(r["database_accession"]).strip().upper() for r in refs
            if (r.get("database_name") or "").lower() == "uniprot"
            and isinstance(r.get("database_accession"), str)
            and r["database_accession"].strip()))
        genes: list[str] = []
        for org in ent.get("rcsb_entity_source_organism") or []:
            for gn in (org or {}).get("rcsb_gene_name") or []:
                if gn and gn.get("value"):
                    genes.append(gn["value"])
        out.append(
            {
                "chains": [c for c in (ids.get("auth_asym_ids") or [])],
                "uniprot": uniprots[0] if uniprots else None,
                "uniprots": uniprots,
                "gene": genes[0] if genes else None,
                "description": (ent.get("rcsb_polymer_entity") or {}).get("pdbx_description") or "",
                "sequence": (ent.get("entity_poly") or {}).get("pdbx_seq_one_letter_code_can") or "",
            }
        )
    return out


def uniprot_identity(uniprot: str) -> dict[str, Any]:
    """Exact accession identity evidence for otherwise absent gene annotations."""
    import hashlib

    accession = uniprot.strip().upper()
    response = _get(_UNIPROT_URL.format(acc=accession), timeout=QUERY_TIMEOUT)
    if not response.ok:
        raise RCSBError(f"UniProt {accession} -> HTTP {response.status_code}")
    data = response.json()
    if not isinstance(data, dict) or data.get("primaryAccession") != accession:
        raise RCSBError(f"UniProt response does not bind exact accession {accession}")
    primary = [str(row["geneName"]["value"]).strip().upper()
               for row in data.get("genes") or []
               if isinstance(row.get("geneName"), dict) and row["geneName"].get("value")]
    return {"accession": accession, "primary_genes": primary,
            "taxon_id": (data.get("organism") or {}).get("taxonId"),
            "reviewed": data.get("entryType") == "UniProtKB reviewed (Swiss-Prot)",
            "sequence_length": len((data.get("sequence") or {}).get("value") or ""),
            "response_sha256": hashlib.sha256(response.content).hexdigest()}


_SITE_FEATURES = frozenset({"Binding site", "Active site", "Site", "Metal binding"})


def uniprot_site_annotations(uniprot: str) -> dict[str, Any]:
    """Sequence and binding/active/metal-site positions (1-based) for ``uniprot``."""
    acc = uniprot.strip().upper()
    resp = _get(_UNIPROT_URL.format(acc=acc), timeout=QUERY_TIMEOUT)
    if not resp.ok:
        raise RCSBError(f"UniProt {acc} -> HTTP {resp.status_code}")
    data = resp.json()
    positions: list[int] = []
    for feat in data.get("features", []) or []:
        if feat.get("type") not in _SITE_FEATURES:
            continue
        loc = feat.get("location") or {}
        start = (loc.get("start") or {}).get("value")
        end = (loc.get("end") or {}).get("value")
        if start is None or end is None:
            continue
        positions.extend(range(int(start), int(end) + 1))
    return {
        "accession": acc,
        "sequence": ((data.get("sequence") or {}).get("value") or ""),
        "positions": sorted(set(positions)),
    }
