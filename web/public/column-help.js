// Numeric direction is an interpretation rule, independent of sort direction.
// pChEMBL definition: https://chembl.gitbook.io/chembl-interface-documentation/frequently-asked-questions/chembl-data-questions
export const columnHelp = {
  discovery_rank: ['점수 순위', '현재 탐색 순위가 있으면 그 순위를, 과거 저장본에서는 당시 종합 점수 순위를 표시합니다.', '숫자가 작을수록 앞선 순위입니다. 해당 분석의 점수와 계산 조건을 함께 확인하세요.'],
  name: ['화합물 / ID', '화합물명과 원자료의 식별자입니다. 같은 이름이라도 염·입체 구조가 다를 수 있습니다.', '후보를 선택해 정확한 구조와 원자료를 확인하세요.'],
  candidate_id: ['후보 ID', '구조·계산·실험 기록을 연결하는 후보 식별자입니다.', '크기나 문자 순서는 성능과 관계없습니다.'],
  discovery_score: ['탐색 점수', '다양성·합성 용이성·구조 경고·피부 적용·안정성을 종합한 탐색 휴리스틱입니다.', '높을수록 탐색 우선순위가 높습니다. 항목별 점수 기여도는 상세에서 확인하세요.'],
  composite: ['기존 종합 점수', '당시 분석 정책으로 계산한 종합 점수입니다. 현재 탐색 점수와 산식이 다를 수 있습니다.', '같은 정책의 결과 안에서 높을수록 우선합니다. 정책이 다른 점수는 직접 비교하지 마세요.'],
  priority_score_v2: ['실험 우선순위', '비교 가능한 실측·구조 근거와 탐색 점수를 종합한 실험 선정 휴리스틱입니다.', '높을수록 우선합니다. 결측은 중립 처리하므로 근거 충족률·불확실성도 함께 보세요.'],
  pchembl: ['pChEMBL', '농도 기반 활성·친화도 값을 몰 단위의 −log10으로 나타냅니다. 7은 100 nM에 해당합니다.', '같은 표적·시험 종류에서 높을수록 강한 potency·affinity 방향입니다. 작용 방향과 실측 여부는 별도 확인하세요.'],
  topical_score: ['피부 적합성', '분자 특성과 피부 적용 관련 경고를 반영한 0–10점 휴리스틱입니다.', '높을수록 분자 물성과 구조 경고를 반영한 이 평가 기준에 적합합니다.'],
  mw: ['MW · 분자량', '분자의 질량을 Da 단위로 표시합니다. 분자 크기를 판단하는 지표입니다.', '피부 적용성·용해도·활성 근거와 함께 판단하세요.'],
  logp: ['logP · 지용성', '옥탄올/물 분배계수의 log10 값입니다. 높을수록 지용성, 낮을수록 친수성입니다.', '용해도·피부 전달·제형을 함께 고려해 적합한 범위를 정하세요.'],
  dock_affinity: ['도킹 / 구조', '기록된 도킹 점수(kcal/mol)와 구조 보기입니다. 원본만 저장된 후보에는 점수가 없을 수 있습니다.', '같은 수용체·포켓·엔진·설정 안에서 비교합니다. 더 낮은 값이 계산상 유리합니다.'],
  primary_docking_score: ['새 도킹 원값', '통합 캠페인에서 검증한 포즈의 도킹 점수(kcal/mol)입니다.', '같은 계산 조건에서 더 낮은 값이 계산상 유리합니다. 조건을 선택한 뒤 비교하세요.'],
  mmgbsa_dG: ['MM/GBSA', '분자역학과 용매 모델로 추정한 결합 에너지(kcal/mol)입니다.', '같은 구조 준비·힘장·용매·계산 조건 안에서 비교합니다. 더 낮은 값이 계산상 유리합니다.'],
  docking_engine: ['기록 엔진', '도킹 점수를 생성한 계산 엔진입니다. 미확인은 실제 엔진 기록을 검증하지 못한 상태입니다.', '같은 엔진과 프로토콜로 계산한 후보끼리 비교하세요.'],
  docking_protocol: ['기록 조건', '수용체·포켓·시드·엔진·검색 설정 등 기록된 도킹 비교 조건입니다.', '같은 계산 조건의 후보끼리 원점수를 비교하세요. 조건별 상세 설정은 근거에서 확인합니다.'],
  admet_hazard: ['ADMET 위험 경보', '흡수·분포·대사·배설·독성 관련 예측과 구조 경고를 요약합니다.', '경보 부담이 낮을수록 탐색에 유리합니다. 경고 없음과 미평가를 구분하며, 안전성 시험은 별도입니다.'],
  direction_verdict: ['작용 방향', '기록된 작용이 연구자가 원하는 억제·활성화 방향과 맞는지 구분합니다.', '목표와 일치하는 후보를 우선 검토합니다. 반대 방향은 대조군으로 활용할 수 있고 도킹만으로 방향을 판정하지 않습니다.'],
  provenance: ['근거 유형', '인체 단백질·세포·대체 표적·수기 주석 등 원자료의 근거 범위입니다.', '같은 숫자여도 근거가 다릅니다. 실제 표적·종·시험·출처를 확인하세요.'],
  purchasable: ['구매 정보', '확인된 공급·판매 기록의 유무입니다. 재고·가격·순도·확인 날짜는 상세에서 확인합니다.', '공급 기록이 있는 후보를 먼저 검토하고, 기록이 없으면 공급사에 문의하세요.'],
  cosmetic_verdict: ['원료 판정', '원료 사용과 관련해 기록된 규제·주석 판정입니다.', '적용 지역·제품·농도와 원자료를 확인하세요.'],
  parent_id: ['호핑 원본', '스카폴드 호핑 후보를 생성하거나 비교할 때 사용한 원본 후보입니다.', '원본의 활성을 새 후보에 물려주지 않습니다. 공통 구조와 변경 부위를 상세에서 확인하세요.'],
  priority_evidence: ['근거 충족 / 불확실성', '충족률은 점수 채널의 평가 범위, 불확실성은 결측·비교군 크기에서 계산한 휴리스틱입니다.', '충족률이 높을수록 평가된 채널이 많습니다. 결측과 비교군 크기를 함께 확인하세요.'],
  computation_status: ['계산 상태', '후보의 계산 완료·대기·실패·보류 상태입니다.', '완료는 계산 결과가 있다는 뜻입니다. 활성·효능 검증 완료와 구분하세요.'],
  detail: ['상세', '원자료·구매 정보·계산 조건·구조와 후보의 선정 근거를 확인합니다.', '점수를 해석하거나 실험 후보를 고정하기 전에 검토하세요.'],
  bookmark: ['관심 후보', '이 브라우저에서 다시 찾아볼 후보를 저장합니다.', '관심 저장은 분석 점수나 실험 선정 결과를 바꾸지 않습니다.']
};

const esc = value => String(value).replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
export function helpButton(key, label = columnHelp[key]?.[0]) {
  if (!columnHelp[key]) return '';
  return `<button type="button" class="column-help" data-column-help="${esc(key)}" aria-label="${esc(label)} 설명" aria-describedby="column-help-popover" aria-expanded="false"><svg viewBox="0 0 24 24" aria-hidden="true" focusable="false"><circle cx="12" cy="12" r="9"/><path d="M12 11v6m0-10v.2"/></svg></button>`;
}

export function installColumnHelp() {
  if (document.getElementById('column-help-popover')) return;
  const tip = document.createElement('div');
  tip.id = 'column-help-popover'; tip.className = 'column-help-popover';
  tip.setAttribute('role', 'tooltip'); tip.hidden = true;
  document.body.append(tip);
  let active = null, pinned = false, hideTimer = null;
  const buttonAt = target => target?.closest?.('[data-column-help]');
  const close = () => {
    clearTimeout(hideTimer);
    active?.setAttribute('aria-expanded', 'false');
    active = null; pinned = false; tip.hidden = true;
  };
  const show = button => {
    const help = columnHelp[button?.dataset.columnHelp];
    if (!help) return;
    clearTimeout(hideTimer);
    if (active !== button) { close(); active = button; }
    tip.innerHTML = `<strong>${esc(help[0])}</strong><p>${esc(help[1])}</p><p class="column-help-direction">${esc(help[2])}</p>`;
    tip.hidden = false; button.setAttribute('aria-expanded', 'true');
    const anchor = button.getBoundingClientRect(), box = tip.getBoundingClientRect(), margin = 12;
    const left = Math.max(margin, Math.min(anchor.left + anchor.width / 2 - box.width / 2, window.innerWidth - box.width - margin));
    const below = anchor.bottom + 10;
    const top = Math.max(margin, Math.min(below + box.height <= window.innerHeight - margin ? below : anchor.top - box.height - 10, window.innerHeight - box.height - margin));
    tip.style.left = left + 'px'; tip.style.top = top + 'px';
  };
  const scheduleClose = () => {
    if (!pinned && document.activeElement !== active) hideTimer = setTimeout(close, 140);
  };
  document.addEventListener('pointerover', event => {
    const button = buttonAt(event.target);
    if (button) show(button);
    else if (tip.contains(event.target)) clearTimeout(hideTimer);
  });
  document.addEventListener('pointerout', event => {
    if ((buttonAt(event.target) === active || tip.contains(event.target)) &&
        !active?.contains(event.relatedTarget) && !tip.contains(event.relatedTarget)) scheduleClose();
  });
  document.addEventListener('focusin', event => {
    const button = buttonAt(event.target);
    if (button) show(button); else close();
  });
  document.addEventListener('click', event => {
    const button = buttonAt(event.target);
    if (button) {
      if (active === button && pinned) close();
      else { show(button); pinned = true; }
    } else if (!tip.contains(event.target)) close();
  });
  document.addEventListener('keydown', event => {
    if (event.key === 'Escape' && active) { close(); event.stopPropagation(); }
  }, true);
  document.addEventListener('scroll', close, true);
  window.addEventListener('resize', close);
  window.addEventListener('hashchange', close);
  new MutationObserver(() => { if (active && !active.isConnected) close(); })
    .observe(document.body, {childList:true, subtree:true});
}
