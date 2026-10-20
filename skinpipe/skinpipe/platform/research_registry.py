"""Owner-bound, immutable input resolutions shared by the novice workflows."""
from __future__ import annotations

from copy import deepcopy
import base64
import io
import json
import re

from .artifacts import canonical_json


def structure_image(smiles):
    from rdkit import Chem
    from rdkit.Chem import Draw
    molecule = Chem.MolFromSmiles(smiles or "")
    if molecule is None or not molecule.GetNumHeavyAtoms():
        return None
    buffer = io.BytesIO()
    Draw.MolToImage(molecule, size=(360, 220)).save(buffer, format="PNG")
    return "data:image/png;base64," + base64.b64encode(buffer.getvalue()).decode("ascii")


class ResearchInputs:
    def __init__(self, service):
        self.service = service
        with service.queue.connection() as db:
            db.execute("""CREATE TABLE IF NOT EXISTS platform_research_resolutions (
                id TEXT NOT NULL, owner_id TEXT NOT NULL, descriptor_sha256 TEXT NOT NULL,
                created_at TIMESTAMPTZ NOT NULL DEFAULT clock_timestamp(),
                PRIMARY KEY(id,owner_id))""")

    def resolve(self, spec, *, owner="researcher"):
        from .research_inputs import resolve_input
        descriptor = resolve_input(
            spec, artifact_store=self.service.artifacts,
            protein_resolver=self.service.discovery._resolve_target,
            source_get=getattr(self.service.discovery, "research_source_get", None),
            structure_resolver=lambda value: self.service.research_structures.resolve(value, owner=owner),
        )
        receipt = self.service.artifacts.put_bytes(canonical_json(descriptor), media_type="application/json")
        identity = receipt["sha256"]
        with self.service.queue.connection() as db:
            db.execute("INSERT INTO platform_research_resolutions (id,owner_id,descriptor_sha256) VALUES (%s,%s,%s) ON CONFLICT DO NOTHING", (identity, owner, identity))
        result = {**descriptor, "resolution_id": identity}
        seed = descriptor.get("query_seed") or {}
        result["structure_image"] = structure_image(descriptor.get("canonical_smiles") or seed.get("smiles"))
        return result

    def get(self, identity, *, owner="researcher"):
        if not isinstance(identity, str) or not re.fullmatch(r"[0-9a-f]{64}", identity):
            raise ValueError("Invalid research input resolution")
        with self.service.queue.connection() as db:
            row = db.execute("SELECT descriptor_sha256 FROM platform_research_resolutions WHERE id=%s AND owner_id=%s", (identity, owner)).fetchone()
        if row is None:
            raise KeyError("Research input resolution is unavailable for this owner")
        return json.loads(self.service.artifacts.get_bytes(row["descriptor_sha256"]))

    def resolve_for_run(self, body, *, owner="researcher"):
        if body.get("resolution_id") is not None:
            if body.get("input") is not None or body.get("target") is not None:
                raise ValueError("Choose a frozen resolution or an input, not both")
            descriptor = self.get(body["resolution_id"], owner=owner)
        else:
            descriptor = self.resolve(body.get("input"), owner=owner)
            descriptor.pop("structure_image", None)
            descriptor.pop("resolution_id", None)
        if descriptor.get("status") != "verified":
            raise ValueError("먼저 입력 대상과 필요한 선택을 확인해 주세요.")
        return deepcopy(descriptor)

    def image(self, identity, *, owner="researcher"):
        descriptor = self.get(identity, owner=owner)
        seed = descriptor.get("query_seed") or {}
        return {"image_data_url": structure_image(descriptor.get("canonical_smiles") or seed.get("smiles")), "scope": "exact resolved chemical structure; 2D"}

    def structure(self, identity, *, owner="researcher"):
        descriptor = self.get(identity, owner=owner)
        if descriptor.get("kind") != "bound_complex" or descriptor.get("status") != "verified":
            raise ValueError("확인된 결합 구조에서만 관측 자세를 볼 수 있습니다.")
        seed = descriptor["query_seed"]
        def asset(sha, format, label):
            raw = self.service.artifacts.get_bytes(sha)
            if len(raw) > 50 * 1024**2:
                raise ValueError("브라우저 구조 표시 크기를 초과했습니다.")
            return {"data": raw.decode("utf-8"), "format": format, "label": label,
                    "source_file": label + "." + format, "source_sha256": sha}
        return {"receptor": asset(descriptor["protein_structure_sha256"],
                    "mmcif" if descriptor.get("protein_structure_format") == "mmcif" else "pdb", "선택한 단백질"),
                "ligand": {**asset(seed["observed_sdf_sha256"], "sdf", "관측된 기준 리간드"),
                           "kind": "observed_reference"},
                "scope": "observed ligand and selected protein in original coordinate frame; no docking",
                "pocket_box": descriptor.get("pocket_box"), "selection": descriptor.get("selection"),
                "warnings": [{"component": "도킹", "status": "보류", "reason": reason}
                             for reason in descriptor.get("docking_hold_reasons", [])]}
