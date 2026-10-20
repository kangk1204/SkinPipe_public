"""Explicit human-protein exploration independent of the skin marker table."""
from __future__ import annotations

import copy
from pathlib import Path

import yaml

from .fetch import chembl
from .types import Direction, DruggableNode, TargetMode, TargetSpec


class ProteinRegistry:
    def __init__(self, identity: dict, *, action: str = "any"):
        if identity.get("status") != "VERIFIED":
            raise ValueError("A verified protein identity is required.")
        if action not in ("any", "inhibit", "activate"):
            raise ValueError("action must be any, inhibit or activate")
        self.identity = copy.deepcopy(identity)
        self.action = action
        self._spec = None

    def spec(self, marker: str) -> TargetSpec:
        if marker.upper() != self.identity["query"]:
            raise ValueError("The requested identifier differs from the verified identity.")
        if self._spec is not None:
            return copy.deepcopy(self._spec)
        identity = self.identity
        gene = identity.get("gene") or identity["uniprot"]
        errors = []
        counts = {"pdb": "NOT_EVALUATED", "activities_p6": "NOT_EVALUATED"}
        n_pdb = n_ligands = 0
        try:
            n_pdb = chembl.pdb_count(identity["uniprot"])
            counts["pdb"] = "EVALUATED"
        except Exception as exc:  # counts are advisory; identity was verified separately
            errors.append(f"PDB count unavailable: {type(exc).__name__}: {exc}")
        tid = identity.get("chembl_target")
        if tid:
            try:
                n_ligands = chembl.activity_count(tid, pchembl_min=6)
                counts["activities_p6"] = "EVALUATED"
            except Exception as exc:
                errors.append(f"Activity count unavailable: {type(exc).__name__}: {exc}")
        elif identity["chembl_mapping_status"] == "NOT_FOUND":
            counts["activities_p6"] = "NO_TARGET_MAPPING"
        else:
            errors.append(identity.get("chembl_mapping_note") or "ChEMBL mapping not evaluated")
        # Explicit input cannot bypass curated structure/metal limitations.
        config = Path(__file__).resolve().parent.parent / "config/skin_targets.yaml"
        configured = [n for m in yaml.safe_load(config.read_text())["markers"]
                      for n in m.get("nodes", []) if n.get("gene") == gene]
        unusable = [n for n in configured if n.get("structure_usable") is False]
        metals = sorted({n["metalloenzyme"] for n in configured if n.get("metalloenzyme")})
        support = {"status": "user_requested_exploration", "target_identity": identity,
                   "biological_goal_status": "NOT_EVALUATED", "counts": counts,
                   "requested_action": self.action,
                   "note": "Identity verified; desired action is a user request, not validated skin biology."}
        node = DruggableNode(gene=gene, uniprot=identity["uniprot"], chembl_target=tid,
                             n_pdb=n_pdb, n_ligands_p6=n_ligands,
                             node_action=None if self.action == "any" else self.action,
                             goal_support=support, primary=True, preferred=True,
                             potency_floor=0, potency_floor_why="broad protein exploration; no default potency exclusion",
                             structure_usable=not unusable,
                             structure_note="; ".join(sorted({n.get("structure_note", "") for n in unusable})),
                             metalloenzyme=", ".join(metals) or None,
                             resolved=not errors, resolve_error="; ".join(errors))
        self._spec = TargetSpec(marker=marker, category="custom_protein",
                                direction=Direction.UNKNOWN if self.action == "any" else Direction(self.action),
                                mode=TargetMode.DIRECT, nodes=[node],
                                notes=[support["note"]])
        return copy.deepcopy(self._spec)
