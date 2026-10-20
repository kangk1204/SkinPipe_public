"""Owner-bound, optional explanations beside immutable scientific calculations.

The deterministic summary is always available. An opted-in ChatGPT request can
select and organize supplied facts; it cannot contribute new values or scores.
"""
from __future__ import annotations

from copy import deepcopy
import heapq
import json
import math
import re

from .artifacts import canonical_json

TASKS = {"input_guidance", "execution_explanation", "result_explanation", "followup"}
SCIENTIFIC_STATUS = "계산 결과는 연구 후보를 좁히는 근거입니다. 효능·안전성·결합은 실험으로 확인해야 합니다."
GUIDANCE = [
    ("input.compound", "화합물은 개별 물질 이름이나 정확한 SMILES로 확인합니다. 여러 분자를 포함하는 이름은 구체적인 종류를 선택해야 합니다."),
    ("input.protein", "단백질은 이름 또는 UniProt 식별자와 종을 확인합니다. 입력을 임의의 대리 표적으로 바꾸지 않습니다."),
    ("input.complex", "결합 구조는 단백질 사슬과 결합 화합물을 직접 선택합니다. 관측 좌표와 정확한 화학적 정체성은 구분합니다."),
    ("interpretation.scope", SCIENTIFIC_STATUS),
]
TASK_TITLES = {"input_guidance": "입력 확인 도움", "execution_explanation": "실행 전 확인",
               "result_explanation": "결과 읽기", "followup": "다음 연구 준비"}
MODE_LABELS = {"similarity": "비슷한 후보 찾기", "scaffold_hopping": "다른 골격의 후보 찾기",
               "assessment": "입력 물질 먼저 평가하기"}
KIND_LABELS = {"compound": "화합물", "protein": "단백질", "bound_complex": "단백질–리간드 결합 구조"}
DELIVERY_LABELS = {"exploratory": "탐색적 검토", "surface": "피부 표면", "epidermal": "표피", "dermal": "진피"}
REASON_MESSAGES = {
    "compound_class_requires_exact_structure": "여러 분자를 포함하는 이름입니다. 정확한 종류나 개별 구조를 선택해 주세요.",
    "protein_family_requires_exact_member": "단백질 가족 이름입니다. 계산할 구성원 하나와 생물종을 선택해 주세요.",
    "invalid_smiles": "구조 문자열을 읽을 수 없습니다. 원래 SMILES를 확인해 주세요.",
    "multiple_organic_components_require_selection": "여러 유기 성분이 있습니다. 계산할 개별 물질을 선택해 주세요.",
    "concrete_structure_required": "반복 단위나 미정 원자가 있습니다. 개별 분자를 나타내는 구조가 필요합니다.",
    "ambiguous_compound_name": "이름에 대응하는 구조가 여러 개입니다. 표시된 구조에서 하나를 선택해 주세요.",
    "ambiguous_protein_name": "이름에 대응하는 단백질이 여러 개입니다. 공식 이름과 생물종을 확인해 주세요.",
    "identity_not_found": "일치하는 대상을 찾지 못했습니다. 다른 이름이나 정확한 구조·식별자로 다시 확인해 주세요.",
    "model_selection_required": "구조 모델이 여러 개입니다. 사용할 모델을 선택해 주세요.",
    "protein_chain_selection_required": "계산에 사용할 단백질 사슬을 선택해 주세요.",
    "bound_ligand_selection_required": "결합 화합물이 여러 개입니다. 출발 화합물 하나를 선택해 주세요.",
    "assembly_selection_required": "단백질 조립체를 선택해 주세요. 선택한 조립체의 좌표를 계산에 사용합니다.",
    "ligand_altloc_selection_required": "리간드에 대체 좌표가 있습니다. 사용할 관측 위치를 선택해 주세요.",
    "ligand_exact_chemistry_required": "좌표만으로 화학 구조를 확정하지 못했습니다. 정확한 SMILES 또는 SDF를 보완해 주세요.",
}
STAGE_GUIDANCE = {
    "input_guidance": [
        ("input.choose", "화합물·단백질·결합 구조 중 가지고 있는 자료에 맞는 시작점을 선택합니다."),
        ("input.identity", "이름 조회와 구조 검증으로 대상을 확인합니다. AI가 이름이나 대리 표적을 확정하지 않습니다."),
    ],
    "execution_explanation": [
        ("execution.order", "확인한 입력 → DB 후보 탐색 → 피부·독성 평가 → 필요한 호핑·도킹 → 결과 저장 순서로 진행합니다."),
        ("execution.unknown_exposure", "제형이나 농도가 미정이면 분자 물성과 구조 경보를 평가하고, 해당 조건이 필요한 전신 노출 계산은 미평가로 둡니다."),
        ("execution.optional_ai", "ChatGPT 설명이 없어도 기본 분석과 Excel 내보내기는 진행됩니다."),
    ],
    "result_explanation": [
        ("result.counts", "후보·출발물질·대조물질은 따로 셉니다. 제외 후보도 전체 표와 내보내기에 보존합니다."),
        ("result.channels", "도킹·호핑·피부 특성·독성은 서로 다른 근거입니다. 미평가를 0점이나 안전으로 해석하지 않습니다."),
    ],
    "followup": [
        ("followup.scope", "후속 연구는 저장된 근거의 빈틈을 확인하는 제안입니다. 계산만으로 효과·안전성·적정 사용 농도를 확정하지 않습니다."),
        ("followup.experiment", "후보와 대조물질의 역할을 구분하고, 실측 활성과 피부 적용 조건에 맞는 실험을 검토합니다."),
    ],
}


def _permissions(value):
    if (not isinstance(value, dict) or set(value) - {"research_data", "public_literature"}
            or any(type(x) is not bool for x in value.values())):
        raise ValueError("Research assistance permissions must be explicit booleans")
    return {"public_literature": True, "research_data": False, **value}


def _execution(value):
    """Use the same scientific validators as submission, without creating a run."""
    from .catalog_assessment import freeze_skin_context
    from .discovery import freeze_catalog_limit
    from .research_generation import freeze_generation
    allowed = {"search_mode", "skin_context", "generation", "compute_budget_seconds", "ai_permissions", "catalog_limit", "screening"}
    if not isinstance(value, dict) or set(value) - allowed:
        raise ValueError("Execution explanation accepts only proposed research conditions")
    mode = value.get("search_mode", "similarity")
    if not isinstance(mode, str) or mode not in MODE_LABELS:
        raise ValueError("Unsupported proposed search mode")
    generation = freeze_generation(value.get("generation"))
    if mode == "assessment" and generation["enabled"]:
        raise ValueError("Generation requires similarity or scaffold_hopping mode")
    budget = _number(value.get("compute_budget_seconds", 7200))
    if budget is None or not 0 < budget <= 28800:
        raise ValueError("Discovery budget must be within 8h actual processing")
    from .routine_screening import freeze_screening
    screening = freeze_screening(value.get("screening"), compute_budget_seconds=budget)
    return {"screening": {"target_candidates": screening["target_candidates"]},
            "search_mode": mode, "skin_context": freeze_skin_context(value.get("skin_context")),
            "generation": generation, "compute_budget_seconds": budget,
            "catalog_limit": freeze_catalog_limit(value.get("catalog_limit")),
            "ai_permissions": _permissions(value.get("ai_permissions", {}))}


def _number(value):
    return value if type(value) in (int, float) and math.isfinite(value) else None


def _messages(values):
    if not isinstance(values, list):
        return []
    return [x for x in values if isinstance(x, str) and x][:30]


def _card(row):
    reasons, concerns, missing = [], [], []
    decision = row.get('routine_disposition') or {}
    if decision.get('state') == 'selected':
        labels = {'saved_priority': '저장된 연구 내 우선순위로 선정했습니다.',
            'scaffold_or_acyclic_topology_diversity': '다른 골격 또는 비고리 연결 구조를 포함하도록 선정했습니다.',
            'eligible_priority_backfill': '평가 기준을 통과한 후보의 우선순위로 선정했습니다.'}
        reasons.extend(labels.get(r, r) for r in decision.get('reasons', []))
    selection = row.get("selection_reason")
    if isinstance(selection, str) and selection:
        reasons.append(selection)
    reasons.extend(_messages(row.get("selection_reasons")))
    similarity = row.get("reference_similarity") or {}
    alignment = row.get("reference_alignment") or {}
    topology = row.get("acyclic_topology_comparison") or {}
    # The producer stores these channels separately. The adapter preserves
    # their exact source fields while presenting one comparison record.
    hop = {**deepcopy(row.get("hop") or {}),
           "reference_similarity": deepcopy(similarity), "reference_alignment": deepcopy(alignment),
           "acyclic_topology_comparison": deepcopy(topology)}
    if "novel_scaffold" in row:
        hop["novel_scaffold"] = row["novel_scaffold"]
    if similarity:
        hop.update(similarity_2d=similarity.get("tanimoto_2d"), reference_id=similarity.get("reference_id"))
    if alignment:
        hop.update(alignment_status=alignment.get("status"), alignment_reason=alignment.get("reason"))
        for key in ("shape", "pharmacophore_color_score", "retention_score"):
            hop.pop(key, None)
        if alignment.get("status") == "complete":
            hop.update(shape=alignment.get("shape_score"),
                       pharmacophore_color_score=alignment.get("pharmacophore_color_score"),
                       retention_score=alignment.get("retention_score"),
                       reference_id=alignment.get("reference_id") or hop.get("reference_id"))
        else:
            missing.append("3차원 형태 비교가 완료되지 않았습니다: " + str(alignment.get("reason") or alignment.get("status") or "미확인"))
    if hop.get("novel_scaffold") is True:
        reasons.append("저장한 참조와 다른 골격 비교 후보입니다. 특허 신규성이나 활성을 입증하지 않습니다.")
    if topology.get("different_from_query_topologies") is True:
        reasons.append("비고리 분자의 연결 위상 표현이 출발구조와 다릅니다. 고리 골격 변화나 활성을 입증하지 않습니다.")
    for key, label in (("similarity_2d", "출발구조와의 2차원 지문 유사도"),
                       ("shape", "완료한 3차원 형태 비교"),
                       ("pharmacophore_color_score", "완료한 약리 특징 색상 비교"),
                       ("retention_score", "완료한 형태·특징 비교의 평균")):
        value = _number(hop.get(key))
        if value is not None and 0 <= value <= 1:
            reasons.append(f"{label}: {value:g}")
    if row.get("derivation"):
        reasons.append("선택적으로 생성한 새 분자 구조입니다. 구매·합성 가능성은 확인되지 않았습니다.")
    if row.get("query_identity_match") is True:
        concerns.append("출발구조와 동일한 분자입니다. 새로운 후보나 활성 대조군으로 해석하지 않습니다.")
    profile = row.get("structural_admet") or {}
    severity = row.get("hazard_severity") or profile.get("hazard_severity")
    if severity and severity != "none":
        concerns.append("구조 기반 경보: " + str(severity))
    concerns.extend(_messages(row.get("discovery_filter_reasons")))
    risk = row.get("admet_search_risk") or {}
    if isinstance(risk, dict):
        concerns.extend(_messages(risk.get("reasons")))
    missing.extend(_messages(row.get("priority_missing_channels")))
    if not row.get("measurements"):
        missing.append("이 연구 조건의 실측 활성 자료가 없습니다.")
    if not row.get("primary_docking"):
        missing.append("검증된 주 도킹 결과가 없습니다.")
    if not row.get("admet"):
        missing.append("기존 DB의 독성·ADMET 예측이 없습니다. 추가 모델 예측은 별도 채널에서 확인하세요.")
    context = row.get("structural_admet_context") or {}
    if context.get("product") is None or context.get("applied_percent") is None:
        missing.append("제품 형태 또는 배합 농도가 미정이므로 전신 노출량을 확정하지 않습니다.")
    return {"candidate_id": row["candidate_id"], "name": row.get("name") or row["candidate_id"],
            "smiles": row.get("smiles"), "selection_reasons": list(dict.fromkeys(reasons)) or ["저장한 탐색 정책의 후보 순서에 따라 표시했습니다."],
            "concerns": list(dict.fromkeys(concerns)), "missing": list(dict.fromkeys(missing)),
            "skin": deepcopy(row.get("skin")), "admet": deepcopy(profile),
            "cached_toxicity": deepcopy(row.get("admet") or []),
            "source_chunk_sha256": row.get("source_chunk_sha256"),
            "additional_toxicity": deepcopy(row.get("additional_toxicity")),
            "primary_docking": deepcopy(row.get("primary_docking")),
            "discovery_origin": deepcopy(row.get("discovery_origin")),
            "priority_score_v2": row.get("priority_score_v2"),
            "filter_pass": row.get("discovery_filter_pass"), "hop": deepcopy(hop)}


class ResearchAssistance:
    def __init__(self, service, *, initialize=True):
        self.service = service
        if initialize:
            with service.queue.connection() as db:
                db.execute("""CREATE TABLE IF NOT EXISTS platform_research_assistance (
                    id TEXT PRIMARY KEY, owner_id TEXT NOT NULL, run_id TEXT, resolution_id TEXT,
                    artifact_sha256 TEXT NOT NULL, created_at TIMESTAMPTZ NOT NULL DEFAULT clock_timestamp())""")

    def _last(self, run_id, owner):
        with self.service.queue.connection() as db:
            row = db.execute("SELECT artifact_sha256 FROM platform_research_assistance WHERE run_id=%s AND owner_id=%s ORDER BY created_at DESC,id DESC LIMIT 1", (run_id, owner)).fetchone()
        if not row:
            return {"status": "NOT_REQUESTED", "summary": None}
        result = json.loads(self.service.artifacts.get_bytes(row["artifact_sha256"]))
        if result.get("owner_id") != owner or result.get("run_id") != run_id:
            raise ValueError("Research explanation owner binding failed")
        return {"status": result["status"], "summary": result.get("summary"),
                "artifact_sha256": row["artifact_sha256"], "human_review_required": True}

    def summary(self, run_id, *, owner="researcher"):
        run = self.service.discovery.get(run_id, owner=owner)
        record = self.service.discovery._record(run_id, owner)
        request = record["request"]
        descriptor = deepcopy(run.get("target_descriptor") or request.get("resolved_input") or {})
        research_input = deepcopy(run.get("research_input") or request.get("input") or request.get("target") or {})
        descriptor.setdefault("input", research_input)
        descriptor.setdefault("name", research_input.get("value") or research_input.get("query") or descriptor.get("target_id"))
        snapshot = run.get("result_snapshot_id") or run.get("snapshot_id")
        result = {"schema": "skinpipe-novice-research-summary-v1", "run_id": run_id,
                  "input": descriptor, "analysis_mode": request.get("search_mode") or run.get("route"),
                  "status": run["status"] if (run.get("cancellation") or {}).get("requested") else run.get("pipeline_status") or run["status"],
                  "lifecycle_state": run.get("lifecycle_state"), "cancellation": deepcopy(run.get("cancellation")),
                  "progress": deepcopy(run.get("progress")), "counts": deepcopy(run.get("counts") or {}),
                  "candidates": [], "seed_assessment": deepcopy(run.get("seed_assessment")),
                  "skin_context": deepcopy(run.get("skin_context") or request.get("skin_context") or {}),
                  "protein_receptor_context": deepcopy(run.get("protein_receptor_context") or request.get("protein_receptor_context")),
                  "bound_complex_docking": deepcopy(run.get("bound_complex_docking")),
                  "catalog_selection": deepcopy(run.get("catalog_selection")),
                  "routine_selection": deepcopy(run.get('routine_selection')),
                  "candidate_preview_limit": 20,
                  "limitations": [], "next_actions": [], "scientific_status": SCIENTIFIC_STATUS,
                  "snapshot_id": snapshot, "target_id": run.get("target_id") or descriptor.get("target_id"),
                  "execution": {key: deepcopy(request.get(key)) for key in
                                ("search_mode", "skin_context", "generation", "compute_budget_seconds", "ai_permissions", "catalog_limit", "screening_policy") if key in request},
                  "ai_permissions": deepcopy(request.get("ai_permissions") or {"public_literature": True, "research_data": False})}
        if not snapshot:
            result["limitations"].append("후보 모집단이 아직 저장되지 않았습니다. 미완료 상태를 후보 0개로 해석하지 않습니다.")
            result["next_actions"].append("실행 기록의 처리 단계와 보류 사유를 확인하고 상태를 갱신해 주세요.")
        else:
            if self.service.catalog.snapshot_owner(snapshot) != owner:
                raise PermissionError("Research snapshot belongs to another owner")
            manifest = self.service.catalog.manifest(snapshot)
            if result["target_id"] not in manifest["target_counts"]:
                raise ValueError("Research result target differs from its saved population")
            best, census, ordinal = [], {"candidate_population": 0, "query_seeds": 0, "controls": 0, "filter_pass": 0}, 0
            for row in self.service.catalog.iter_candidates(snapshot, result["target_id"]):
                role = row.get("role", "candidate")
                if role == "query_seed":
                    census["query_seeds"] += 1
                    if result["seed_assessment"] is None:
                        result["seed_assessment"] = deepcopy(row)
                    continue
                if role != "candidate":
                    census["controls"] += 1
                    continue
                census["candidate_population"] += 1
                if row.get("discovery_filter_pass") is False:
                    continue
                census["filter_pass"] += 1
                if request.get('screening_policy') is not None and row.get('routine_disposition', {}).get('state') != 'selected':
                    continue
                score = _number(row.get("priority_score_v2"))
                ordinal += 1
                item = (score if score is not None else float("-inf"), -ordinal, _card(row))
                if len(best) < 20:
                    heapq.heappush(best, item)
                elif item[:2] > best[0][:2]:
                    heapq.heapreplace(best, item)
            result["counts"].update(census)
            result["candidates"] = [item[2] for item in sorted(best, key=lambda x: (-x[0], -x[1]))]
            provenance = manifest.get("provenance") or {}
            catalog_receipt = provenance.get("catalog_receipt") or {}
            if catalog_receipt.get("completeness") not in {None, "complete", "not_requested"}:
                result["limitations"].append("카탈로그 확보가 제한적이거나 미완료입니다. 전체 카탈로그를 탐색했다고 해석하지 않습니다.")
            result["limitations"].extend(_messages(catalog_receipt.get("limitations")))
            if descriptor.get("kind") == "compound":
                result["limitations"].append("화합물 단독 탐색입니다. 단백질 결합이나 효능의 활성 자료를 새로 만들지 않습니다.")
            if descriptor.get("kind") == "bound_complex":
                result["limitations"].append("관측된 결합 구조에서 시작했습니다. 화학 정체성·포즈 품질·단백질 조립체 검토는 별도로 필요합니다.")
            if request.get("search_mode") == "assessment":
                result["next_actions"].append("입력 화합물의 피부 관련 지표·경보·미평가 항목을 검토하고 필요한 실험과 적용 조건을 정해 주세요.")
            elif census["candidate_population"] == 0:
                result["next_actions"].append("출발구조 평가를 먼저 확인하고 카탈로그 연결 또는 다른 후보 모집단을 검토해 주세요.")
            elif census["filter_pass"] == 0:
                result["next_actions"].append("전체 후보 표에서 제외 이유를 확인하고 필요한 필터 조건을 조정해 주세요.")
            else:
                result["next_actions"].append("후보별 선택 이유·경보·미평가 항목을 검토한 뒤 소규모 실험 후보를 정해 주세요.")
            result["next_actions"].append("전체 후보와 개별 물질 탭을 Excel로 받아 실험 조건과 확인 자료를 기록해 주세요.")
        result["limitations"].append("출발구조·생성 구조·카탈로그 후보·측정 대조군은 서로 다른 역할입니다. 계산의 완료는 생물학적 검증을 뜻하지 않습니다.")
        # Draw only the bounded selected DTOs, never every row scanned for the
        # census. These local 2D images do not modify stored observations or
        # enter the explicitly selected ChatGPT fact packet.
        from .research_registry import structure_image
        from .cached_prediction_views import CachedPredictionView
        prediction_view = CachedPredictionView(self.service.artifacts)
        display_rows = list(result["candidates"])
        if isinstance(result["seed_assessment"], dict):
            display_rows.append(result["seed_assessment"])
        for display in display_rows:
            key = "cached_toxicity" if "cached_toxicity" in display else "admet"
            display[key] = prediction_view.for_row(display, key=key)
            try:
                smiles = display.get("smiles")
                display["structure_image"] = structure_image(smiles) if isinstance(smiles, str) else None
            except (ValueError, TypeError, RuntimeError, OSError, ImportError):
                display["structure_image"] = None
        result["ai"] = self._last(run_id, owner)
        return result

    @staticmethod
    def _packet(task, summary, descriptor, private, execution=None):
        """A bounded stage contract, never a serialization of research records.

        Even an opted-in packet excludes raw structures, coordinates, arbitrary
        user questions and files. Only supplied, locally materialized facts may
        be selected by the language model.
        """
        facts, actions, candidate_ids = [], [], []
        structures = {descriptor.get(k) for k in ("canonical_smiles", "exact_smiles") if isinstance(descriptor.get(k), str)}
        spec = descriptor.get("input") or {}
        if spec.get("format") == "smiles" and isinstance(spec.get("value"), str):
            structures.add(spec["value"])
        if summary:
            structures.update(c["smiles"] for c in summary.get("candidates", []) if isinstance(c.get("smiles"), str))
            seed = summary.get("seed_assessment") or {}
            if isinstance(seed.get("smiles"), str):
                structures.add(seed["smiles"])

        def text_only(value):
            text = str(value)
            # Identity fields can legitimately contain the raw SMILES as name.
            # Replace it before either the local or remote fact is produced.
            for structure in sorted(structures, key=len, reverse=True):
                if structure:
                    text = re.sub(r"(?<![A-Za-z0-9])" + re.escape(structure) + r"(?![A-Za-z0-9])", "[화합물 구조]", text)
            return text[:1900]

        def fact(key, text, **metadata):
            if len(facts) < 128:
                facts.append({"fact_id": key, "text": text_only(text), **metadata})

        def action(key, text):
            if len(actions) < 40:
                actions.append({"action_id": key, "text": text_only(text)})

        for key, text in STAGE_GUIDANCE[task]:
            fact(key, text)
        fact("interpretation.scope", SCIENTIFIC_STATUS)
        if task == "input_guidance":
            for key, text in GUIDANCE[:3]:
                if not private or not descriptor or descriptor.get("kind") == key.split(".")[1].replace("complex", "bound_complex"):
                    fact(key, text)
            action("input.confirm", "확인 화면의 이름·구조·생물종 또는 선택한 리간드가 원래 의도와 맞는지 확인해 주세요.")
        elif task == "execution_explanation":
            action("execution.confirm", "아래 탐색 방식·피부 조건·계산 예산을 확인한 뒤 실행 버튼을 눌러 주세요.")
            action("execution.privacy", "ChatGPT에 전달을 허용한 항목을 확인해 주세요. 원본 구조 파일과 좌표는 설명 요청에 포함하지 않습니다.")
        elif task == "result_explanation":
            action("result.compare", "먼저 검토할 후보에서 선정 이유·경보·미평가 항목을 비교해 주세요.")
            action("result.export", "전체 후보와 제외 후보를 Excel로 내려받아 개별 물질 탭의 근거를 확인해 주세요.")
        else:
            action("followup.controls", "연구 목적에 맞는 양성·음성 대조군과 실측 평가 방법을 정해 주세요. 출발물질을 검증된 대조군으로 간주하지 않습니다.")
            action("followup.conditions", "제형·적용 부위·농도를 확인하고, 실험으로 검증할 범위를 정해 주세요.")
            action("followup.export", "Excel에 후보와 조건별 확인 항목을 기록하고 원자료·문헌 근거를 함께 검토해 주세요.")

        if private and descriptor:
            kind = descriptor.get("kind") or spec.get("kind")
            fact("research.kind", "입력 유형: " + KIND_LABELS.get(kind, "확인 중"))
            name = descriptor.get("name") or spec.get("value") or (descriptor.get("request") or {}).get("query")
            if spec.get("format") == "upload":
                name = "업로드한 결합 구조"
            elif spec.get("format") == "smiles":
                name = "SMILES로 입력한 화합물"
            if name:
                fact("research.input", "확인할 연구 입력: " + str(name)[:500])
            status = descriptor.get("status")
            verified = status == "verified"
            fact("input.status", "대상 확인 상태: " + ("확인됨" if verified else "추가 확인 필요"), value=status)
            if not verified:
                message = REASON_MESSAGES.get(descriptor.get("reason"))
                if not message:
                    message = ("조회 출처에 연결하지 못했습니다. 연결을 확인한 뒤 다시 검색해 주세요."
                               if status == "source_unavailable" else "확인 화면의 필요한 선택이나 보완 자료를 입력해 주세요.")
                fact("input.required", message)
                action("input.resolve", message)
                choices = descriptor.get("choices")
                if isinstance(choices, list) and choices:
                    fact("input.choices", f"저장된 선택지 수: {len(choices)}", value=len(choices))
            if kind in {"protein", "bound_complex"}:
                taxon = descriptor.get("taxon_id")
                fact("input.species", "생물종: " + ("사람 (9606)" if taxon == 9606 else f"확인된 분류 ID {taxon}" if type(taxon) is int else "미확인"))
                direction = descriptor.get("desired_action") or spec.get("desired_action")
                fact("input.direction", "원하는 작용 방향: " + {"activate": "활성화", "inhibit": "억제"}.get(direction, "방향 미정"))
                for index, component in enumerate(descriptor.get("components", [])[:8]):
                    accession = component.get("accession") if isinstance(component, dict) else None
                    if isinstance(accession, str) and re.fullmatch(r"[A-Z0-9-]{1,24}", accession):
                        fact(f"input.protein.{index}", "확인한 단백질 구성원 UniProt ID: " + accession)
            if kind == "compound":
                fact("input.no_proxy", "입력 화합물을 합성·분해에 관여하는 단백질로 자동 치환하지 않습니다.")
                cid = descriptor.get("pubchem_cid")
                if type(cid) is int and cid > 0:
                    fact("input.pubchem", f"확인한 개별 화합물 PubChem CID: {cid}", value=cid)
            if kind == "bound_complex":
                fact("input.observed_pose", "선택한 결합 리간드의 관측 좌표를 비교 기준으로 사용합니다. 후보 도킹은 같은 수용체와 포켓 조건을 유지합니다.")
                selection = descriptor.get("selection") or {}
                for key, label in (("model_id", "구조 모델"), ("assembly_id", "단백질 조립체"), ("ligand_code", "기준 리간드"),
                                   ("ligand_chain", "리간드 사슬"), ("ligand_residue", "리간드 잔기")):
                    if key in selection:
                        fact("input.selection." + key, label + ": " + str(selection[key]))
                for index, reason in enumerate(descriptor.get("docking_hold_reasons", [])[:4]):
                    fact(f"input.docking_hold.{index}", "도킹 단계 보류 사유: " + str(reason))

        if private and task == "execution_explanation":
            proposal = execution or ((summary or {}).get("execution") or {})
            kind = descriptor.get("kind") or spec.get("kind")
            mode = proposal.get("search_mode")
            if mode in MODE_LABELS:
                fact("execution.mode", "선택한 탐색 방식: " + MODE_LABELS[mode])
            else:
                fact("execution.conditions_pending", "실행 조건이 아직 전달되지 않았습니다. 실행 화면의 조건을 선택한 뒤 다시 확인해 주세요.")
            if kind == "compound":
                fact("execution.docking", "단백질을 지정하지 않은 화합물 경로는 도킹 없이 완료하며, 표적 결합과 활성 유지 여부는 미평가로 둡니다.")
            elif kind == "protein":
                fact("execution.docking", "선택한 단백질의 구조와 포켓이 확정되면 도킹합니다. 확정하지 못하면 도킹만 보류하고 확보한 후보의 피부·독성 평가는 제공합니다.")
            elif kind == "bound_complex":
                fact("execution.docking", "같은 수용체·관측 좌표계·포켓·검색 조건으로 재도킹과 후보 도킹을 진행합니다. 지원하지 않는 화학·복합체 조건은 도킹만 보류합니다.")
            if mode == "assessment":
                fact("execution.assessment", "입력 물질을 평가합니다. 새로운 후보 모집단 탐색과 신규 구조 생성은 요청하지 않습니다.")
            context = proposal.get("skin_context") or {}
            fact("execution.skin", "모든 후보에 피부 관련 물성·구조 경보를 평가하고, 정확한 구조와 모델 조건이 일치하는 기존 독성 예측을 재사용합니다. 결측값은 미평가로 유지합니다.")
            if context:
                fact("execution.delivery", "피부 도달 목적: " + DELIVERY_LABELS.get(context.get("goal"), "미정"))
                for field, label in (("research_goal", "연구 목적"), ("application_site", "적용 부위"), ("product", "제형")):
                    fact("condition." + field, label + ": " + str(context.get(field) or "아직 모름"))
                percent = context.get("applied_percent")
                fact("condition.applied_percent", "배합 농도: " + (f"{percent}%" if percent is not None else "아직 모름"), value=percent)
            budget = _number(proposal.get("compute_budget_seconds"))
            if budget is not None:
                fact("execution.budget", f"최대 실제 처리 예산: {budget}초", value=budget)
            limit = proposal.get("catalog_limit")
            if type(limit) is int and mode != "assessment":
                fact("execution.catalog_limit", f"전역 검색 인덱스에서 가져와 평가할 최대 DB 후보: {limit}개. DB 전체 분자를 평가했다는 뜻은 아닙니다.", value=limit)
            screening = proposal.get("screening_policy") or proposal.get("screening") or {}
            goal = screening.get("target_candidates")
            if type(goal) is int and mode != "assessment":
                fact("execution.representatives", f"대표 후보 목표: {goal}개" + ". 구조·피부·독성·포즈 QC를 통과한 후보에서 선정하며 부족한 수와 이유를 기록합니다.", value=goal)
            generation = proposal.get("generation") or {}
            fact("execution.generation", "신규 구조 생성: " + ("선택함. 생성 구조의 구매·합성·활성·독성은 별도로 확인합니다." if generation.get("enabled") else "사용하지 않음"))
            if proposal.get("ai_permissions"):
                allowed = proposal["ai_permissions"].get("research_data") is True
                fact("execution.ai_permissions", "이 실행의 비공개 연구 요약 전달: " + ("허용" if allowed else "허용하지 않음"))

        if private and summary is not None and task in {"result_explanation", "followup"}:
            fact("research.status", "저장한 실행 상태: " + str(summary["status"]), value=summary["status"])
            census_keys = ("candidate_population", "query_seeds", "controls", "filter_pass")
            count_items = sorted(summary["counts"].items(), key=lambda item: census_keys.index(item[0]) if item[0] in census_keys else len(census_keys))
            for key, value in count_items[:16]:
                if _number(value) is not None and isinstance(key, str) and re.fullmatch(r"[A-Za-z0-9_.:-]{1,100}", key):
                    labels = {"candidate_population": "후보 모집단", "query_seeds": "출발물질", "controls": "대조물질", "filter_pass": "조건 통과 후보", "collected": "확보 후보"}
                    fact("count." + key, f"저장한 {labels.get(key, key)}: {value}", value=value)
            if summary.get("analysis_mode") == "assessment":
                fact("result.assessment", "입력 물질 평가 결과입니다. 후보 수가 0이어도 출발물질의 평가와 미평가 항목을 확인할 수 있습니다.")
                seed = summary.get("seed_assessment")
                if isinstance(seed, dict) and seed.get("candidate_id"):
                    card = _card(seed)
                    for channel in ("concerns", "missing"):
                        for index, message in enumerate(card[channel][:4]):
                            fact(f"seed.{channel}.{index}", "출발물질: " + message)
            for index, candidate in enumerate(summary["candidates"][:6]):
                cid = candidate["candidate_id"]
                if not isinstance(cid, str) or not re.fullmatch(r"[A-Za-z0-9_.:-]{1,160}", cid):
                    continue
                candidate_ids.append(cid)
                for channel in ("selection_reasons", "concerns", "missing"):
                    for ordinal, message in enumerate(candidate[channel][:3]):
                        fact(f"candidate.{index}.{channel}.{ordinal}", f"{cid}: {message}", candidate_id=cid)
                if task == "followup":
                    if not candidate.get("primary_docking"):
                        action(f"followup.{index}.binding", f"{cid}: 표적 결합이 필요한 연구라면 적절한 결합·기능 실험 또는 검증 가능한 표적·포켓 자료를 확보해 주세요.")
                    if candidate.get("missing"):
                        action(f"followup.{index}.missing", f"{cid}: 미평가 항목을 먼저 확인하고 필요한 측정 자료를 확보해 주세요.")
            for index, message in enumerate(summary["limitations"][:8]):
                fact(f"limitation.{index}", message)
            for index, message in enumerate(summary["next_actions"][:8]):
                action(f"next.{index}", message)
        return {"schema": "skinpipe-research-explanation-packet-v1", "task": task,
                "research_data_allowed": private, "facts": facts, "candidate_ids": candidate_ids,
                "source_ids": [], "actions": actions}

    @staticmethod
    def _deterministic(packet):
        from ..literature.chatgpt import validate_explanation
        sections = {"input": [], "meaning": [], "concerns": [], "missing": [], "next": []}
        essential = ["research.kind", "research.input", "input.status", "input.required", "research.status",
                     "execution.mode", "result.assessment", "count.candidate_population", "count.query_seeds",
                     "count.controls", "count.filter_pass", "execution.budget"]
        order = {key: index for index, key in enumerate(essential)}
        ranked = sorted(packet["facts"], key=lambda item: (order.get(item["fact_id"], 20 if ".selection_reasons." in item["fact_id"] else 30),))
        for fact in ranked:
            key = fact["fact_id"]
            section = ("concerns" if ".concerns." in key or ".docking_hold." in key else
                       "missing" if ".missing." in key or key.startswith("limitation.") or key == "input.required" else
                       "input" if key.startswith(("input.", "research.input", "research.kind")) else
                       "next" if key.startswith("followup.") else "meaning")
            if len(sections[section]) < 20:
                sections[section].append(key)
        # Preserve the exact status and numbers, but never call a zero candidate
        # population a failed input assessment or imply a pipeline is complete.
        available = {x["fact_id"] for x in packet["facts"]}
        preferred = [key for key in essential if key in available]
        if not preferred:
            preferred = [x["fact_id"] for x in packet["facts"][:3]]
        return validate_explanation({"summary_fact_ids": preferred[:12],
            "sections": [{"section_id": key, "fact_ids": value} for key, value in sections.items() if value],
            "next_action_ids": [x["action_id"] for x in packet["actions"][:12]],
            "human_review_required": True}, packet)

    def assist(self, body, *, owner="researcher"):
        if not isinstance(body, dict) or set(body) - {"run_id", "resolution_id", "task", "permissions", "execution"}:
            raise ValueError("Research assistance accepts a saved run or resolution, task and explicit permissions")
        task = body.get("task", "input_guidance")
        if task not in TASKS:
            raise ValueError("Unsupported research assistance task")
        run_id, resolution_id = body.get("run_id"), body.get("resolution_id")
        if run_id and resolution_id:
            raise ValueError("Choose a run or input resolution")
        if not run_id and not resolution_id and task != "input_guidance":
            raise ValueError("This assistance task requires a saved run or input resolution")
        if task in {"result_explanation", "followup"} and not run_id:
            raise ValueError("Result explanation and followup require a saved research run")
        if "execution" in body and (task != "execution_explanation" or not resolution_id):
            raise ValueError("Proposed execution conditions require an execution explanation of a saved input resolution")
        execution = _execution(body["execution"]) if "execution" in body else None
        permissions = _permissions(body.get("permissions", {"research_data": False}))
        summary, descriptor, frozen = None, {}, {}
        if run_id:
            summary = self.summary(run_id, owner=owner)
            descriptor, frozen = summary["input"], summary["ai_permissions"]
        elif resolution_id:
            descriptor = self.service.research_inputs.get(resolution_id, owner=owner)
        if task == "execution_explanation" and descriptor.get("status") != "verified":
            raise ValueError("Confirm the input identity before requesting an execution explanation")
        private = permissions.get("research_data") is True
        if run_id and private and frozen.get("research_data") is not True:
            raise PermissionError("This run did not authorize sharing private research data")
        packet = self._packet(task, summary, descriptor, private, execution)
        local_packet = self._packet(task, summary, descriptor, True, execution)
        deterministic = self._deterministic(local_packet)
        result = {"schema": "skinpipe-research-assistance-v1", "owner_id": owner, "run_id": run_id,
                  "resolution_id": resolution_id, "task": task, "permissions": deepcopy(permissions),
                  "task_title": TASK_TITLES[task], "status": "DETERMINISTIC", "summary": deterministic["summary"],
                  "explanation": deterministic, "explanation_scope": "LOCAL_SAVED_FACTS",
                  "facts": deepcopy(packet["facts"]), "human_review_required": True,
                  "scientific_validation": "NOT_ESTABLISHED", "scores_modified": False,
                  "research_data_transmitted": False, "research_data_transmission": "NOT_ATTEMPTED"}
        personal = getattr(self.service, "personal_literature", None)
        if private and personal is not None and personal.configured:
            personal._guard(owner)
            from ..literature.chatgpt import ChatGPTError, validate_explanation
            try:
                result.update(research_data_transmitted=None, research_data_transmission="ATTEMPTED_OUTCOME_UNKNOWN")
                with personal._client() as client:
                    generated = client.explain(packet)
                result.update(research_data_transmitted=True, research_data_transmission="COMPLETED_RESPONSE_PENDING_VALIDATION")
                checked = validate_explanation(generated.get("selection"), packet)
                if not checked["summary_facts"] and not any(x["facts"] for x in checked["sections"]) and not checked["next_actions"]:
                    raise ChatGPTError("empty_research_explanation")
                result.update(status="AI_FACTS_ORGANIZED", summary=checked["summary"], explanation=checked,
                              explanation_scope="PERMITTED_FACT_SELECTION",
                              llm_receipt=deepcopy(generated.get("llm_receipt")), research_data_transmitted=True,
                              research_data_transmission="COMPLETED")
            except (ChatGPTError, OSError, ValueError, KeyError, TypeError) as error:
                from .personal_literature import _safe_error
                result.update(status="DETERMINISTIC_AI_UNAVAILABLE", error=_safe_error(error))
                if result["research_data_transmitted"] is True:
                    result["research_data_transmission"] = "COMPLETED_RESPONSE_REJECTED"
        elif not private:
            result["status"] = "DETERMINISTIC_PRIVATE_DATA_WITHHELD"
        else:
            result["status"] = "DETERMINISTIC_AI_NOT_CONFIGURED"
        receipt = self.service.artifacts.put_bytes(canonical_json(result), media_type="application/json")
        identity = receipt["sha256"]
        with self.service.queue.connection() as db:
            db.execute("INSERT INTO platform_research_assistance (id,owner_id,run_id,resolution_id,artifact_sha256) VALUES (%s,%s,%s,%s,%s) ON CONFLICT DO NOTHING",
                       (identity, owner, run_id, resolution_id, identity))
        return {key: value for key, value in result.items() if key != "owner_id"} | {"artifact_sha256": identity}
