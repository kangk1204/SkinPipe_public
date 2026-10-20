"""Species-aware target resolution. Identity is separate from compute readiness.

Source responses are retained by content hash. Neither search order, skin
marker labels nor client-supplied verification flags establish identity.
"""

from __future__ import annotations

from datetime import datetime, timezone
import re

from .identity import digest
from .artifacts import canonical_json
from .ranking import validate_requested_action

SCHEMA = "target-descriptor-v2"
ACCESSION = re.compile(
    r"(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-[0-9]+)?\Z"
)
CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"
UNIPROT = "https://rest.uniprot.org/uniprotkb"


def now():
    return datetime.now(timezone.utc).isoformat()


def target_request(body):
    if not isinstance(body, dict) or set(body) - {
        "query",
        "kind",
        "taxon_id",
        "desired_action",
        "components",
        "marker",
        "isoform",
    }:
        raise ValueError("Unsupported target descriptor fields")
    value = dict(body)
    query = value.get("query", "")
    if not isinstance(query, str) or not re.fullmatch(
        r"[A-Za-z0-9_.:+ -]{1,160}", query
    ):
        raise ValueError("Enter a bounded target identifier")
    kind = value.get("kind", "protein")
    if kind not in {"protein", "complex", "endpoint"}:
        raise ValueError("kind must be protein, complex or endpoint")
    taxon = value.get("taxon_id", 9606)
    if (
        isinstance(taxon, bool)
        or not isinstance(taxon, int)
        or not 1 <= taxon <= 100000000
    ):
        raise ValueError("taxon_id must be an NCBI taxonomy identifier")
    components = value.get("components", [])
    if (
        not isinstance(components, list)
        or len(components) > 16
        or any(not isinstance(c, str) or not ACCESSION.fullmatch(c) for c in components)
    ):
        raise ValueError(
            "Complex components must be at most 16 exact UniProt accessions"
        )
    if len(set(components)) != len(components):
        raise ValueError("Duplicate complex components")
    if value.get("isoform") is not None and (
        not isinstance(value["isoform"], str)
        or not ACCESSION.fullmatch(value["isoform"])
    ):
        raise ValueError("Explicit isoform accession required")
    return {
        **value,
        "query": query.strip(),
        "kind": kind,
        "taxon_id": taxon,
        "species_defaulted": "taxon_id" not in body,
        "desired_action": validate_requested_action(value.get("desired_action", "any")),
        "components": sorted(components),
    }


class SourceClient:
    """Bounded, provider-only requests with immutable raw response receipts."""

    def __init__(self, artifacts, *, get=None):
        self.artifacts, self.get = artifacts, get
        self.receipts = []

    def request(self, url, params=None):
        import requests
        from .artifacts import canonical_json

        if not (url.startswith(CHEMBL + "/") or url.startswith(UNIPROT + "/")):
            raise ValueError("Only configured public identity providers are allowed")
        params = params or {}
        attempted = now()
        try:
            if self.get:
                data = self.get(url, params)
            else:
                with requests.get(
                    url,
                    params=params,
                    timeout=(5, 20),
                    allow_redirects=False,
                    stream=True,
                ) as response:
                    response.raise_for_status()
                    if response.status_code != 200:
                        raise ValueError(
                            "Provider response was not a complete HTTP 200"
                        )
                    raw = bytearray()
                    for block in response.iter_content(65536):
                        raw.extend(block)
                        if len(raw) > 16 * 1024 * 1024:
                            raise ValueError("Provider page exceeds 16 MiB")
                    import json

                    data = json.loads(
                        raw,
                        parse_constant=lambda _: (_ for _ in ()).throw(
                            ValueError("Nonfinite source JSON")
                        ),
                    )
            receipt = self.artifacts.put_bytes(
                canonical_json(data), media_type="application/json"
            )
            self.receipts.append(
                {
                    "url": url,
                    "params": params,
                    "retrieved_at": attempted,
                    "raw_sha256": receipt["sha256"],
                    "status": "complete",
                }
            )
            return data
        except (requests.RequestException, ValueError, OSError):
            self.receipts.append(
                {
                    "url": url,
                    "params": params,
                    "retrieved_at": attempted,
                    "status": "source_unavailable",
                    "raw_sha256": None,
                }
            )
            raise


def _local_mapping_receipt(client, local, result, *, operation, query, error=None):
    """The pinned DB hash identifies a source; a CAS hash seals this query result."""
    receipt = {'source': 'ChEMBL SQLite', 'version': local.source_version,
        'sha256': local.source_sha256, 'operation': operation, 'query': query,
        'retrieved_at': now(), 'raw_sha256': None,
        'status': 'source_unavailable' if error else 'complete' if result.get('complete') else 'partial_acquisition'}
    if error:
        receipt['error_type'] = type(error).__name__
    else:
        receipt['raw_sha256'] = client.artifacts.put_bytes(canonical_json(result), media_type='application/json')['sha256']
    client.receipts.append(receipt)
    return receipt


class ProviderIdentityConflict(ValueError):
    def __init__(self, requested, response, client):
        super().__init__('Requested accession differs from the provider primary accession')
        aliases = response.get('secondaryAccessions', [])
        self.choice = {'requested_accession': requested, 'accession': response.get('primaryAccession'),
            'secondary_alias_declared': isinstance(aliases, list) and requested in aliases,
            'source_record_sha256': client.artifacts.put_bytes(canonical_json(response))['sha256']}


def _requested_protein(client, accession):
    response = client.request(UNIPROT + '/' + accession + '.json')
    if response.get('primaryAccession') != accession:
        # A declared secondary accession is evidence for a reviewable alias,
        # not permission to silently change a frozen target/component identity.
        raise ProviderIdentityConflict(accession, response, client)
    return response


def _exact_gene_names(protein):
    """Provider-declared primary names and gene synonyms, with exact matching.

    UniProt gene_exact also finds explicit synonyms (for example mouse Car2
    has primary gene Ca2). Locus/ORF labels and substring matches do not resolve
    a gene. Species, complete-page and accession-ambiguity checks stay separate.
    """
    names = set()
    for gene in protein.get("genes", []):
        if not isinstance(gene, dict):
            continue
        labels = [gene.get("geneName")]
        synonyms = gene.get("synonyms", [])
        if isinstance(synonyms, list):
            labels.extend(synonyms)
        for label in labels:
            if isinstance(label, dict) and isinstance(label.get("value"), str):
                names.add(label["value"].upper())
    return names


def resolve_target(body, client, *, local=None):
    """Resolve identities; a configured verified local release is mapping authority.

    UniProt still verifies accession/species, including locally found ChEMBL
    components. Local absence is absence in that pinned release, not all sources.
    """
    request = target_request(body)
    descriptor = {
        "schema": SCHEMA,
        "request": request,
        "kind": request["kind"],
        "taxon_id": request["taxon_id"],
        "desired_action": request["desired_action"],
        "isoform": request.get("isoform"),
        "components": [],
        "status": "selection_required",
        "compute_readiness": "identity_pending",
        "choices": [],
        "verified_at": None,
    }
    if request["kind"] == "endpoint":
        return {**descriptor, "reason": "mechanism_selection_required", "sources": []}
    query = request["query"].upper()
    if ACCESSION.fullmatch(query) and "-" in query:
        return {
            **descriptor,
            "reason": "isoform_sequence_verification_required",
            "sources": [],
        }
    proteins, targets = [], []
    mapping_status = "not_queried"
    local_mapping, mapping_source = None, None
    try:
        if re.fullmatch(r"CHEMBL[0-9]+", query):
            if local:
                try:
                    local_mapping = local.target_record(query)
                    _local_mapping_receipt(client, local, local_mapping, operation='exact_target', query={'target_chembl_id': query})
                    if local_mapping['target'] is not None and local_mapping['complete']:
                        targets = [local_mapping['target']]
                        mapping_status, mapping_source = 'complete', local_mapping['source_provenance']
                except Exception as error:
                    _local_mapping_receipt(client, local, {}, operation='exact_target', query={'target_chembl_id': query}, error=error)
            if not targets:
                targets = [client.request(CHEMBL + "/target/" + query + ".json")]
                mapping_status, mapping_source = "complete", {'source': 'ChEMBL online'}
        elif request["kind"] == "complex":
            if not request["components"]:
                return {
                    **descriptor,
                    "reason": "complex_components_or_exact_chembl_id_required",
                    "sources": [],
                }
            for accession in request["components"]:
                proteins.append(_requested_protein(client, accession))
        else:
            # Exact accession requests include the real organism, so a default
            # human setting cannot hide a nonhuman identifier conflict.
            if ACCESSION.fullmatch(query):
                proteins = [_requested_protein(client, query)]
            else:
                if not re.fullmatch(r"[A-Z][A-Z0-9-]{0,29}", query):
                    raise ValueError("Use an exact gene symbol, UniProt or ChEMBL ID")
                result = client.request(
                    UNIPROT + "/search",
                    {
                        "query": f"gene_exact:{query} AND organism_id:{request['taxon_id']} AND reviewed:true",
                        "format": "json",
                        "size": 100,
                    },
                )
                proteins = result.get("results", [])
                if len(proteins) >= 100:
                    return {
                        **descriptor,
                        "reason": "identity_page_incomplete",
                        "sources": client.receipts,
                    }
                proteins = [
                    p
                    for p in proteins
                    if query in _exact_gene_names(p)
                ]
        if targets:
            target = targets[0]
            if re.fullmatch(r"CHEMBL[0-9]+", query) and target.get('target_chembl_id') != query:
                return {**descriptor, 'reason': 'provider_target_id_conflict', 'choices': targets,
                    'sources': client.receipts, 'target_mapping_status': 'identity_conflict',
                    'target_mapping_source': mapping_source}
            required_type = (
                "PROTEIN COMPLEX" if request["kind"] == "complex" else "SINGLE PROTEIN"
            )
            if target.get("target_type") != required_type or not target.get(
                "target_components"
            ):
                return {
                    **descriptor,
                    "reason": "target_kind_mismatch",
                    "choices": targets,
                    "sources": client.receipts,
                }
            for component in target["target_components"]:
                accession = component.get("accession", "")
                if not ACCESSION.fullmatch(accession):
                    return {
                        **descriptor,
                        "reason": "unresolved_complex_component",
                        "sources": client.receipts,
                    }
                proteins.append(_requested_protein(client, accession))
            if request['components'] and sorted(c['accession'] for c in target['target_components']) != request['components']:
                return {**descriptor, 'reason': 'complex_component_set_conflict', 'choices': targets, 'sources': client.receipts}
            if int(target.get("tax_id") or 0) != request["taxon_id"]:
                return {
                    **descriptor,
                    "reason": "species_conflict",
                    "choices": targets,
                    "sources": client.receipts,
                }
        choices = [
            {
                "accession": p.get("primaryAccession"),
                "taxon_id": p.get("organism", {}).get("taxonId"),
                "organism": p.get("organism", {}).get("scientificName"),
                "gene_names": sorted(
                    {
                        g.get("geneName", {}).get("value")
                        for g in p.get("genes", [])
                        if g.get("geneName", {}).get("value")
                    }
                ),
                "reviewed": p.get("entryType") == "UniProtKB reviewed (Swiss-Prot)",
                "source_record_sha256": client.artifacts.put_bytes(canonical_json(p))[
                    "sha256"
                ],
            }
            for p in proteins
        ]
        if not choices:
            return {
                **descriptor,
                "reason": "identity_not_found",
                "sources": client.receipts,
            }
        if any(p["taxon_id"] != request["taxon_id"] for p in choices):
            return {
                **descriptor,
                "reason": "species_conflict",
                "choices": choices,
                "sources": client.receipts,
            }
        if request["kind"] == "protein" and len(choices) != 1:
            return {
                **descriptor,
                "reason": "ambiguous_identity",
                "choices": choices,
                "sources": client.receipts,
            }
        if any(not ACCESSION.fullmatch(p.get("accession") or "") for p in choices):
            raise ValueError("Provider did not return exact protein identities")
        if request.get("isoform"):
            # UniProt's canonical entry does not itself validate an isoform's
            # sequence. Preserve it, but do not substitute the canonical chain.
            return {
                **descriptor,
                "reason": "isoform_sequence_verification_required",
                "choices": choices,
                "sources": client.receipts,
            }
        if not targets:
            if local:
                try:
                    local_mapping = local.target_mappings([p['accession'] for p in choices], taxon_id=request['taxon_id'], kind=request['kind'])
                    _local_mapping_receipt(client, local, local_mapping, operation='verified_component_set', query=local_mapping['query'])
                    targets = local_mapping['targets']
                    mapping_status = 'complete' if local_mapping['complete'] else 'partial_acquisition'
                    mapping_source = local_mapping['source_provenance']
                    if not local_mapping['complete'] and len(targets) <= 1:
                        targets = []  # unresolved components cannot establish absence or a unique mapping
                    if request['kind'] == 'complex' and not targets and mapping_status == 'complete':
                        # A complete registry lookup is not a complete assay
                        # query for a complex that has no exact target record.
                        mapping_status = 'mapping_absent'
                except Exception as error:
                    _local_mapping_receipt(client, local, {}, operation='verified_component_set',
                        query={'accessions': sorted(p['accession'] for p in choices), 'taxon_id': request['taxon_id'], 'kind': request['kind']}, error=error)
            if mapping_source is None:
                mapping_source = {'source': 'ChEMBL online'}
                mapping_status, targets = _online_target_mapping(request, choices, client)
                if request['kind'] == 'complex' and not targets and mapping_status == 'complete':
                    mapping_status = 'mapping_absent'
        if len(targets) > 1:
            return {
                **descriptor,
                "reason": "ambiguous_chembl_mapping",
                "choices": targets,
                "components": choices,
                "sources": client.receipts,
                "target_mapping_status": mapping_status,
                "target_mapping_source": mapping_source,
                "local_target_mapping": local_mapping,
            }
        if targets and sorted(c.get('accession') or '' for c in targets[0].get('target_components', [])) != sorted(p['accession'] for p in choices):
            return {**descriptor, 'reason': 'target_component_identity_conflict', 'choices': targets,
                'components': choices, 'sources': client.receipts,
                'target_mapping_status': 'identity_conflict', 'target_mapping_source': mapping_source}
        target_id = (
            "complex:" + digest(sorted(p["accession"] for p in choices))[:24]
            if request["kind"] == "complex"
            else "uniprot:" + choices[0]["accession"]
        )
        structures = []
        if request["kind"] == "protein":
            for crossref in proteins[0].get("uniProtKBCrossReferences", []):
                if crossref.get("database") != "PDB" or not re.fullmatch(
                    r"[0-9][A-Z0-9]{3}", crossref.get("id", "")
                ):
                    continue
                properties = {
                    v.get("key"): v.get("value") for v in crossref.get("properties", [])
                }
                try:
                    resolution = float(str(properties.get("Resolution", "")).split()[0])
                    if not 0 < resolution < 100:
                        resolution = None
                except (ValueError, IndexError):
                    resolution = None
                structures.append(
                    {
                        "pdb_id": crossref["id"],
                        "resolution_angstrom": resolution,
                        "method": properties.get("Method"),
                        "chains": properties.get("Chains"),
                        "source_record_sha256": choices[0]["source_record_sha256"],
                    }
                )
            structures.sort(
                key=lambda r: (
                    r["resolution_angstrom"] is None,
                    r["resolution_angstrom"] or 0,
                    r["pdb_id"],
                )
            )
        return {
            **descriptor,
            "status": "verified",
            "target_id": target_id,
            "components": choices,
            "chembl_target_id": targets[0].get("target_chembl_id") if targets else None,
            "target_mapping_status": mapping_status,
            "target_mapping_source": mapping_source,
            "local_target_mapping": local_mapping,
            "organism": choices[0]["organism"],
            "compute_readiness": "receptor_verification_required",
            "assembly_verification_required": request["kind"] == "complex",
            "verified_at": now(),
            "sources": client.receipts,
            "structure_choices": structures,
            "structure_selection_policy": "best recorded resolution, deterministic PDB ID; native identity and pocket validation still required",
            "biological_goal_status": "not_established",
        }
    except ProviderIdentityConflict as error:
        return {**descriptor, 'reason': 'provider_accession_conflict', 'choices': [error.choice],
                'sources': client.receipts, 'target_mapping_status': 'identity_conflict'}
    except Exception as error:
        return {
            **descriptor,
            "status": "source_unavailable",
            "reason": type(error).__name__,
            "sources": client.receipts,
        }


def _online_target_mapping(request, choices, client):
    """Fallback lookup verifies the complete component set, never the first result."""
    required_type = 'SINGLE PROTEIN' if request['kind'] == 'protein' else 'PROTEIN COMPLEX'
    expected = sorted(p['accession'] for p in choices)
    try:
        mapping = client.request(CHEMBL + '/target.json', {
            'target_components__accession': expected[0], 'limit': 100})
        if not isinstance(mapping, dict) or not isinstance(mapping.get('targets'), list):
            raise ValueError('Target response does not contain a target list')
        page = mapping.get('page_meta')
        if (not isinstance(page, dict) or 'next' not in page
            or isinstance(page.get('total_count'), bool) or not isinstance(page.get('total_count'), int)
            or page['total_count'] < len(mapping['targets'])
            or (page['next'] is not None and (not isinstance(page['next'], str) or not page['next']))
            or (page['next'] is None and page['total_count'] != len(mapping['targets']))):
            raise ValueError('Target response pagination is incomplete or inconsistent')
        for target in mapping['targets']:
            if (not isinstance(target, dict) or not re.fullmatch(r'CHEMBL[0-9]+', target.get('target_chembl_id') or '')
                or not isinstance(target.get('target_components'), list)
                or any(not isinstance(c, dict) for c in target['target_components'])):
                raise ValueError('Target response does not contain exact target records')
        targets = [t for t in mapping['targets']
            if t.get('target_type') == required_type
            and int(t.get('tax_id') or 0) == request['taxon_id']
            and sorted(c.get('accession') or '' for c in t.get('target_components', [])) == expected]
        if page['next']:
            return 'partial_acquisition', []
        return 'complete', targets
    except Exception as error:
        # An HTTP 200 raw receipt remains available, but its invalid payload
        # must not certify a complete lookup or the absence of references.
        if client.receipts and client.receipts[-1].get('url') == CHEMBL + '/target.json':
            client.receipts[-1]['validation_status'] = 'source_unavailable'
            client.receipts[-1]['validation_error_type'] = type(error).__name__
        return 'source_unavailable', []
