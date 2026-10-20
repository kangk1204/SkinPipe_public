"""In-process use of the read-only Scout gateway; no network or donor writes."""
from __future__ import annotations

import re
import secrets
import time

from .identity import chemical_identity


def lookup_exact_structures(scout, smiles_list):
    """Indexed exact graph lookup shared by local and private remote adapters.

    The supplied text is only a parameter to the fixed indexed query. Parent,
    similarity and name matches never supply predictions for this operation.
    """
    import duckdb
    from rdkit import Chem
    from .catalog_cancellation import catalog_checkpoint
    if not isinstance(smiles_list, list) or len(smiles_list) > 1000:
        raise ValueError("At most 1000 exact compound structures per lookup")
    results = []
    for smiles in smiles_list:
        catalog_checkpoint()
        identity = chemical_identity(smiles)
        result = {"query_smiles": smiles, "identity": identity, "matches": [], "status": "not_found"}
        if identity["structure_status"] != "valid":
            result["status"] = "invalid_structure"
            results.append(result)
            continue
        def exact_matches(response):
            catalog_checkpoint()
            rows = response.get("rows", []) if isinstance(response, dict) else response
            if len(rows) > 1000:
                raise ValueError("Exact structure lookup exceeds the bounded duplicate limit")
            ids = [row["compound_id"] for row in rows if row.get("compound_id")]
            if not ids:
                return []
            catalog_checkpoint()
            materialized = scout.materialize(ids)
            catalog_checkpoint()
            matches = []
            for record in materialized:
                catalog_checkpoint()
                compound = record.get("compound", {})
                actual = chemical_identity(compound.get("identity_smiles") or compound.get("canonical_smiles"))
                if actual.get("exact_structure_id") == identity["exact_structure_id"]:
                    matches.append(record)
            return matches

        # The catalog indexes InChIKey, but not either SMILES column. A text
        # OR scan over 17M rows can exhaust the query deadline and the research
        # budget. The key retrieves candidates only: salts, stereochemistry,
        # isotopes and tautomers still pass the exact graph check above.
        try:
            key = Chem.MolToInchiKey(Chem.MolFromSmiles(identity["identity_smiles"]))
        except (ValueError, RuntimeError):
            key = None
        if key:
            try:
                indexed = scout._query(
                    "SELECT compound_id FROM compounds WHERE inchikey=? ORDER BY compound_id LIMIT 1001",
                    [key],
                )
            except duckdb.BinderException as error:
                # Older private catalogs may lack this optional column.
                # Timeouts, cancellation and other query failures propagate.
                if "inchikey" not in str(error) or "not found" not in str(error):
                    raise
            else:
                result["matches"] = exact_matches(indexed)
        if not result["matches"]:
            result["matches"] = exact_matches(scout._query(
                "SELECT compound_id FROM compounds WHERE identity_smiles=? OR canonical_smiles=? ORDER BY compound_id LIMIT 1001",
                [identity["identity_smiles"], identity["canonical_smiles"]],
            ))
        if result["matches"]:
            result["status"] = "matched"
        results.append(result)
    return results


class LocalCatalog:
    def __init__(self, config):
        from .scout_gateway import Gateway
        allowed = {"kind", "scout_root", "published_prediction_cache", "similarity_index",
                   "filtered_similarity_indices", "admet_search_gate", "prediction_partitions"}
        if not isinstance(config, dict) or set(config) - allowed or not config.get("scout_root"):
            raise ValueError("Local database requires an installed Scout root and bounded component paths")
        self.config = dict(config)
        self.gateway = Gateway({**config, "token": secrets.token_hex(32)})

    def request(self, path, body=None):
        if path == "/query" and body is not None:
            return self.gateway.query(body)
        if path == "/materialize" and body is not None:
            return self.gateway.query(body, materialize=True)
        if path == "/similarity" and body is not None:
            return self.gateway.similarity_query(body)
        if path == "/diversity" and body is not None:
            return self.gateway.diversity_query(body)
        match = re.fullmatch(r"/queries/([0-9a-f]{32})", path)
        if match and body is None:
            return self.gateway.result(match[1])
        match = re.fullmatch(r"/queries/([0-9a-f]{32})/cancel", path)
        if match and body is not None:
            if body != {}:
                raise ValueError("Query cancellation requires an empty object")
            return self.gateway.cancel_query(match[1])
        if path in {"/health", "/status"} and body is None:
            return {**self.gateway.scout.status(), "mode": "local_read_only", "admission_errors": self.gateway.admission_errors}
        raise ValueError("Unsupported local catalog operation")

    def lookup(self, smiles_list):
        return lookup_exact_structures(self.gateway.scout, smiles_list)

    def close(self):
        self.gateway.executor.shutdown(wait=True, cancel_futures=True)

    def wait(self, response, *, timeout=120):
        deadline = time.monotonic() + timeout
        while response.get("status") == "running":
            if time.monotonic() >= deadline:
                raise ValueError("Local catalog operation timed out")
            time.sleep(0.05)
            response = self.gateway.result(response["query_id"])
        return response
