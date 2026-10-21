import {escapeHtml as e,identityLabel,PRODUCTS,GOALS} from './research-inputs.js';
import {metricRadarMarkup} from './metric-radar.js';
const MODES={similarity:'비슷한 화합물 찾기',scaffold_hopping:'다른 골격 후보 찾기',assessment:'입력 화합물만 평가',protein_target_screening:'단백질 표적 탐색',bound_complex_scaffold_hopping:'결합 화합물에서 출발한 골격 탐색'};
const COUNTS={candidate_population:'전체 후보 화합물',collected:'확보한 카탈로그 행',raw_generated_proposals:'새 구조 생성 제안',query_seeds:'출발구조',comparison_selected:'구조 비교 선택',docking_completed:'후보 도킹 계산 완료',query_redocking_completed:'관측 리간드 재도킹 계산 완료',saved_rows:'전체 저장 행',filter_pass:'필터 통과 후보',population:'확보한 전체 후보',catalog_population:'카탈로그 후보',catalog_candidates:'카탈로그 후보',filtered:'필터를 통과한 후보',filtered_candidates:'필터를 통과한 후보',candidates:'후보 화합물',hops:'골격 비교 후보',hop_candidates:'골격 비교 후보',generated:'생성한 새 구조',generated_candidates:'생성한 새 구조',selected:'후속 계산 선택',docking_selected:'후보 도킹 선택',docked:'후보 도킹 계산 완료',measured_references:'측정 자료가 있는 참조',computational_seeds:'출발구조',controls:'대조군'};
const STATUS={complete:'처리 완료',partial:'일부 처리 완료',queued:'실행 대기',running:'처리 중',held:'추가 확인 필요',needs_data:'추가 확인으로 보류',failed:'처리 실패',cancelled:'취소됨',cancelling:'취소 요청됨 · 종료 확인 중',cleanup_unverified:'작업 자원 정리 확인 필요'};
export function researchRunStatus(status){return STATUS[status]||({budget_deferred:'처리 시간 한도로 보류',resource_wait:'계산 자원을 기다리는 중',preparing:'입력과 계산 자료 준비 중',completed:'처리 완료'})[status]||'처리 상태 확인 필요 · '+String(status||'미기록')}
export function researchRunCounts(counts={}){
 const pairs=[['대표 후보',counts.representative_selected],counts.new_unique_candidates!==undefined?['새 고유 후보',counts.new_unique_candidates]:['후보 행',counts.candidate_population],['출발구조',counts.query_seeds],['새 구조 생성',counts.generated],['후보 도킹 계산 완료',counts.docking_completed],['관측 리간드 재도킹 계산 완료',counts.query_redocking_completed],['저장 관측',counts.saved_rows]];
 return pairs.filter(([,value])=>typeof value==='number'&&Number.isFinite(value)).map(([label,value])=>label+' '+number(value)).join(' · ');
}
export function researchProgressMarkup(progress,{state='running'}={}){
 // Progress is a throttled intermediate record, not a final completion census.
 // Only an actively processing run can describe it as its current calculation.
 if(!['running','preparing'].includes(state)||!progress||typeof progress!=='object'||Array.isArray(progress))return'';
 const stages={identity_resolution:'입력이 무엇인지 확인',query_seed_frozen:'출발 화합물 구조 보존',reference_acquisition:'비교할 참조 자료 확보',literature_evidence:'공개 문헌 근거 확인',catalog_acquisition:'카탈로그 후보 확보',cheap_structure_assessment:'피부 관련 구조 지표와 경보 평가',additional_toxicity:'DB에 예측이 없는 화합물의 추가 독성 계산',reference_shape:'후보의 구조·형태 비교',bound_complex_docking:'선택한 결합 구조로 도킹',docking:'도킹 계산',cpu_fallback:'CPU 도킹 계산',routine_docking:'후보를 나누어 도킹',routine_pose_quality:'도킹 자세 QC와 대표 후보 선정',candidate_population_frozen:'전체 후보와 평가 근거 저장',ligand_preparation:'계산용 화합물 구조 준비',resource_wait:'계산 자원을 기다리는 중',budget_deferred:'처리 시간 한도로 보류',failed:'처리 오류 확인'};
 const stage=progress.stage,label=Object.hasOwn(stages,stage)?stages[stage]:'처리 단계 확인',valid=value=>Number.isSafeInteger(value)&&value>=0;
 // Worker events merge across stages. Read only counters owned by this emitter;
 // batch totals, cumulative attempts and selected model inputs are different units.
 const pair=(done,total)=>!valid(done)||(valid(total)&&done>total)?'':number(done)+'개'+(valid(total)?' / '+number(total)+'개':'');
 // The latest QC result remains historical when the next docking batch starts.
 const recentQc=valid(progress.qualified)&&valid(progress.target)&&progress.target>0&&progress.qualified<=progress.target?' · 최근 QC 후 대표 후보 '+number(progress.qualified)+'개 / 목표 '+number(progress.target)+'개':'';
 let detail='';
 if(['docking','routine_docking','bound_complex_docking','cpu_fallback'].includes(stage)){
  const current=pair(progress.completed,progress.total);
  detail=current?' · 현재 도킹 묶음 '+current+' 처리 완료 · 자세 QC 전':' · 현재 도킹 묶음의 처리 수 확인 필요 · 자세 QC 전';
  detail+=recentQc;
 }else if(stage==='routine_pose_quality'){
  if(valid(progress.processed))detail+=' · 누적 도킹 시도 '+number(progress.processed)+'개';
  detail+=recentQc;
 }else if(stage==='additional_toxicity'){
  detail=(valid(progress.selected)?' · 이번 단계 선택 '+number(progress.selected)+'개':'')+' · 완료 수는 아직 확정하지 않았습니다.';
 }else if(stage==='catalog_acquisition'){
  if(valid(progress.acquired))detail=' · 카탈로그 행 '+number(progress.acquired)+'개 확보';
 }else if(['cheap_structure_assessment','reference_shape','ligand_preparation'].includes(stage)){
  const current=pair(progress.processed,stage==='cheap_structure_assessment'?progress.total:progress.requested);
  if(current)detail=' · 이 단계에서 '+current+' 처리';
 }else if(stage==='candidate_population_frozen'&&valid(progress.catalog_count)){
  detail=' · 확보한 카탈로그 행 '+number(progress.catalog_count)+'개';
 }
 return`<p class="research-muted">현재 단계: ${e(label)}${detail}</p>`;
}
export function routineSelectionMarkup(selection,previewCount=0){
 if(!selection?.policy||!Number.isInteger(selection.selected_count))return'';
 const counts=selection.counts||{},goal=selection.policy.target_candidates;
 const labels={selected:'대표 후보',alternative:'검토 대안',held:'보류',pending:'평가 대기',duplicate:'중복 구조',control:'출발구조·대조물질'};
 return`<section class="notice research-routine-selection" role="status"><h3>대표 후보 ${number(selection.selected_count)}개 · 목표 최대 ${number(goal)}개</h3><p>저장된 평가를 재사용하고, 필요한 계산을 나누어 진행한 뒤 우선순위와 골격 다양성으로 선정합니다.</p><dl>${Object.entries(labels).filter(([key])=>Number.isInteger(counts[key])).map(([key,label])=>`<dt>${e(label)}</dt><dd>${number(counts[key])}</dd>`).join('')}</dl>${selection.shortfall?'<p>선정 기준을 충족한 수만 표시했습니다. 보류·평가 대기 후보의 사유는 전체 표에서 확인하세요.</p>':''}<p>먼저 ${number(previewCount)}개를 보여드립니다. 전체 대표 후보와 대안·제외·보류 자료는 전체 후보 표와 Excel에서 확인할 수 있습니다.</p></section>`;
}
const number=x=>typeof x==='number'&&Number.isFinite(x)?x.toLocaleString('ko-KR',{maximumFractionDigits:3}):'미평가';
const SEVERITY={none:'해당 구조 경보 미기록',low:'낮은 수준의 경보',moderate:'중간 수준의 경보',high:'높은 수준의 경보',critical:'강한 구조 경보',unknown:'미평가'};
const severity=value=>typeof value==='string'&&Object.hasOwn(SEVERITY,value)?SEVERITY[value]:'미평가 · 정의 확인 필요';
function text(value){
 const raw=typeof value==='string'?value:value?.message||value?.label||value?.reason||value?.name||value?.code||'세부 근거 확인 필요';
 const known={measured:'현재 표적·조건의 실측 활성 근거가 부족합니다.',structural:'비교 가능한 조건의 구조·도킹 근거가 부족합니다.',discovery:'탐색 지표가 미평가입니다.',unknown:'미평가'};
 if(Object.hasOwn(known,raw))return known[raw];
 return raw.replace(/^구조 기반 경보: (none|low|moderate|high|critical|unknown)$/,(all,level)=>'구조 기반 경보: '+severity(level));
}
function items(values,empty){return`<ul>${(Array.isArray(values)&&values.length?values:[empty]).map(x=>`<li>${e(text(x))}</li>`).join('')}</ul>`}
const dockingCount=value=>Number.isSafeInteger(value)&&value>=0?value:null;
const reasonCodes=value=>Array.isArray(value)?value.filter(code=>typeof code==='string'&&code):[];
const DOCKING_HOLDS={explicit_pdb_and_observed_ligand_selection_required:'비교할 PDB 구조와 관측된 결합 화합물을 선택해 주세요.',observed_ligand_required_for_fixed_pocket_redocking:'고정한 포켓의 재도킹에 사용할 관측 화합물이 필요합니다.',requested_protein_ph_not_supported_by_as_written_cpu_preparation:'현재 고정 좌표 준비법은 요청한 pH별 단백질 상태를 지원하지 않아 도킹을 보류했습니다.',docking_engine_not_configured:'도킹 계산 엔진이 아직 설정되지 않았습니다.',query_redocking_preparation_failed:'관측 화합물의 재도킹 준비에 실패했습니다.',query_redocking_engine_failed:'관측 화합물의 재도킹 계산에 실패했습니다.',query_redocking_geometry_failed:'관측 결합 자세의 재현 검사를 통과하지 못해 후보 도킹을 보류했습니다.',query_redocking_top_pose_rmsd_above_threshold:'관측 결합 자세를 재현하는 기준을 통과하지 못해 후보 도킹을 보류했습니다.',research_docking_execution_held:'도킹 처리 기록을 추가로 확인해야 합니다.',routine_docking_execution_held:'후보 도킹 실행을 보류했습니다. 상세 처리 기록을 확인하세요.'};
const SELECTION_REASONS={critical_warning:'강한 경보',structural_alert_review_required:'구조 경보 검토 필요',purchase_record_missing:'구매 자료 미기록',severe_predicted_admet:'독성·ADMET 예측 경보',severe_predicted_admet_review_required:'독성·ADMET 예측 경보 검토 필요',explicit_stereoisomer_required:'정확한 입체이성질체 필요',ambiguous_multicomponent_structure:'여러 성분의 구조 확인 필요',unresolved_structure:'구조 미확정',exact_pose_QC_required:'해당 조건의 도킹 자세 QC 필요',pose_QC_failed:'도킹 자세 QC 불통과',receptor_or_calculation_model_held:'수용체 또는 계산 모델 확인 필요',docking_score_review_required:'도킹 점수 검토 필요',docking_execution_incomplete:'도킹 실행 미완료',exact_model_predictions_required:'해당 구조의 모델 예측 필요',source_exact_structure_mismatch:'원본과 정확 구조 불일치'};
export function researchDockingView(summary={}){
 const counts=summary.counts||{},context=summary.protein_receptor_context,receipt=summary.bound_complex_docking;
 const candidateCount=dockingCount(counts.docking_completed),redockingCount=dockingCount(counts.query_redocking_completed);
 const holds=[...new Set([...reasonCodes(context?.hold_reasons),...reasonCodes(receipt?.hold_reasons)])];
 // Only the explicit current ligand-only contract establishes "not requested".
 // A missing receipt, unknown input kind, or zero count alone cannot do so.
 const notRequested=summary.schema==='skinpipe-novice-research-summary-v1'&&summary.input?.kind==='compound'&&['similarity','scaffold_hopping'].includes(summary.analysis_mode)&&context===null&&receipt===null&&candidateCount===0&&redockingCount===0&&counts.docking_selected===0;
 let state='unknown',message='도킹 요청·실행 근거가 충분히 기록되지 않았습니다. 0만으로 미요청이나 실패를 판단하지 않습니다.';
 if(notRequested){state='not_requested';message='화합물만 입력한 연구로 수용체 도킹은 요청하지 않았습니다. 화합물 구조·피부·독성 평가 결과를 확인하세요.'}
 else if(candidateCount>0){state='calculated';message='후보의 도킹 계산 결과가 있습니다. 계산 완료 수에는 동일 조건의 저장 결과 재사용이 포함될 수 있으며, 자세 QC 통과 수와 다릅니다.'}
 else if(receipt?.redocking?.pass===false){state='qc_failed';message='관측 리간드 재도킹이 자세 재현 기준을 통과하지 못했습니다. 후보 도킹 실행 여부와 보류 기록을 아래에서 확인하세요.'}
 else if(holds.length||context?.status==='held'||receipt?.status==='held'){state='held';message='도킹이 보류되었습니다. 기록된 구조 선택·준비·실행 사유를 확인하세요.'}
 else if(receipt?.status==='failed'){state='failed';message='도킹 처리 실패가 기록되었습니다. 실패 단계와 원자료를 확인하세요.'}
 else if(receipt?.stop_reason==='eligible_pool_exhausted'&&receipt?.planned_candidate_count===0){state='no_eligible_candidates';message='이번 도킹 조건에 맞는 추가 후보를 확보하지 못했습니다. 관측 리간드 재도킹과 후보 도킹은 별개입니다.'}
 else if(receipt&&counts.candidate_population===0&&candidateCount===0){state='no_candidates';message='이번 연구에 후보 행이 0개로 기록되어 있습니다. 후보 확보 범위와 도킹 처리 기록을 확인하세요.'}
 else if(receipt?.status==='complete'&&candidateCount===0){state='none_completed';message='후보 도킹 계산 완료는 0개입니다. 처리 완료 상태만으로 후보 도킹 성공을 뜻하지 않습니다.'}
 const redocking=receipt?.redocking;
 const redockingQc=redocking?.status==='complete'?(redocking.pass===true?'관측 리간드 재도킹 QC 통과':redocking.pass===false?'관측 리간드 재도킹 QC 불통과':'관측 리간드 재도킹 QC 미확인'):redocking?.status==='not_run'?'관측 리간드 재도킹 미실행':redocking?.status==='failed'?'관측 리간드 재도킹 계산 실패':'관측 리간드 재도킹 QC 미확인';
 const reasons=new Map(),decisions=summary.routine_selection?.decisions;
 for(const decision of Object.values(decisions&&typeof decisions==='object'?decisions:{})){
  if(!['held','pending'].includes(decision?.state))continue;
  for(const code of new Set(reasonCodes(decision.reasons))){const key=decision.state+':'+code;reasons.set(key,{state:decision.state,code,count:(reasons.get(key)?.count||0)+1})}
 }
 return{state,message,candidateCount,redockingCount,redockingQc,holds,selectionReasons:[...reasons.values()]};
}
export function researchExecutionScopeMarkup(summary={}){
 const execution=summary.execution||{},counts=summary.counts||{},receipt=summary.bound_complex_docking||{},selection=summary.routine_selection||{};
 const count=value=>dockingCount(value)===null?'미기록':number(value)+'개',requestCount=value=>dockingCount(value)>0?count(value):'미기록';
 const seconds=dockingCount(execution.compute_budget_seconds),time=seconds>0?(seconds%3600===0?number(seconds/3600)+'시간 ('+number(seconds)+'초)':seconds%60===0?number(seconds/60)+'분 ('+number(seconds)+'초)':number(seconds)+'초'):'미기록';
 const facts=rows=>'<dl class="research-evaluation-facts">'+rows.map(([label,value])=>'<dt>'+e(label)+'</dt><dd>'+value+'</dd>').join('')+'</dl>';
 const heading='<section class="research-options research-execution-scope"><h3>이번 연구의 요청과 처리 범위</h3>';
 if(summary.analysis_mode==='assessment')return heading+'<p>입력 화합물만 평가한 연구입니다. 카탈로그 후보 탐색·대표 후보 선정·수용체 도킹은 요청하지 않았습니다.</p>'+facts([['전체 연구 최대 처리 시간',time]])+'</section>';
 const notRequested=researchDockingView(summary).state==='not_requested',docking=value=>notRequested?'미요청':count(value);
 const requested=[['공통 카탈로그 요청 한도',requestCount(execution.catalog_limit)],['대표 후보 목표',requestCount(execution.screening_policy?.target_candidates)],['전체 연구 최대 처리 시간',time]];
 const actual=[['확보한 카탈로그 행',count(counts.collected)],['검색 필터 통과 후보 행',count(counts.filter_pass)],['도킹 계획 후보',docking(receipt.planned_candidate_count)],['후보 도킹 계산 완료',docking(counts.docking_completed)],['이번 연구의 신규 도킹 계산',docking(receipt.new_docking_count)],['동일 조건의 저장 도킹 재사용',docking(receipt.reused_count)],['대표 선정 결과',count(selection.selected_count)]];
 const stops={eligible_pool_exhausted:'이번 조건에 맞는 추가 후보가 없어 후보 도킹을 마쳤습니다.',target_reached:'대표 후보 목표에 도달해 후보 도킹을 마쳤습니다.',budget_deferred:'남은 계산 시간 예산으로 다음 도킹 묶음을 시작할 수 없어 후속 계산을 보류했습니다.',model_held:'계산 모델을 확인해야 해 후속 도킹을 보류했습니다.',pose_quality_runtime_or_input_held:'도킹 자세 QC의 실행 환경 또는 입력을 확인해야 해 후속 계산을 보류했습니다.'};
 const shortfall=dockingCount(selection.shortfall),stop=typeof receipt.stop_reason==='string'?'<p>기록된 후보 도킹 종료 사유: '+e(Object.hasOwn(stops,receipt.stop_reason)?stops[receipt.stop_reason]:'처리 코드: '+receipt.stop_reason)+'</p>':'';
 return heading+'<h4>제출 당시 저장된 요청</h4>'+facts(requested)+'<p>요청 한도와 대표 목표는 보장되는 결과 수가 아닙니다. 과거 기록에 값이 없으면 현재 기본값으로 채우지 않습니다.</p><p class="research-muted">카탈로그 요청 한도는 공통 카탈로그에서 새로 가져올 때 적용됩니다. 기존 연구를 후보 모집단으로 선택했다면 저장된 전체 모집단을 보존합니다. 이 요약에 모집단 출처가 없으면 한도의 실제 적용 여부를 단정하지 않습니다.</p><h4>저장된 처리·선정 결과</h4>'+facts(actual)+'<p class="research-muted">검색 필터 통과는 도킹 자세 QC 통과나 대표 선정을 뜻하지 않습니다. 도킹 계획 후보도 QC를 통과한 최종 후보 수가 아닙니다. 후보 도킹 계산 완료에는 동일 조건의 저장 결과 재사용이 포함될 수 있습니다. 관측 리간드 재도킹은 후보 도킹 수에 합치지 않습니다.</p>'+(shortfall>0?'<p>저장된 대표 선정 결과는 목표보다 '+number(shortfall)+'개 적습니다. 보류·평가 대기 사유는 아래 대표 선정 내역과 전체 후보 표에서 확인하세요.</p>'+stop:'')+'</section>';
}
export function dockingContextMarkup(summary){
 if(summary.analysis_mode==='assessment')return'';
 const context=summary.protein_receptor_context,receipt=summary.bound_complex_docking,view=researchDockingView(summary);
 const labels={unknown:'도킹 근거 미확인',not_requested:'도킹 미요청',calculated:'후보 도킹 계산 기록 있음',held:'도킹 보류',qc_failed:'재도킹 QC 불통과',failed:'도킹 처리 실패',no_eligible_candidates:'추가 도킹 대상 0개',no_candidates:'후보 행 0개',none_completed:'후보 도킹 계산 완료 0개'};
 const count=value=>value===null?'미기록':number(value)+'개';
 const stops={eligible_pool_exhausted:'이번 조건에 맞는 추가 후보가 없어 후보 도킹을 마쳤습니다.',target_reached:'대표 후보 목표에 도달해 후보 도킹을 마쳤습니다.',pose_quality_runtime_or_input_held:'후보 도킹 자세 QC의 실행 환경 또는 입력 확인이 필요해 후속 계산을 보류했습니다.'};
 const redocking=receipt?.redocking,rmsd=redocking?.rmsd_angstrom,threshold=redocking?.threshold_angstrom;
 const dimensions=[['도킹 계획에 포함된 후보',receipt?.planned_candidate_count],['이번 연구의 신규 도킹 계산 기록',receipt?.new_docking_count],['동일 조건의 저장 도킹 재사용',receipt?.reused_count]].filter(([,value])=>dockingCount(value)!==null);
 return`<section class="research-options research-docking-summary" data-docking-state="${e(view.state)}"><h3>후보 도킹과 관측 리간드 재도킹</h3><p><strong>${e(labels[view.state])}</strong> · ${e(view.message)}</p><div class="research-counts"><div><strong>${view.state==='not_requested'?'미요청':count(view.candidateCount)}</strong><span>후보 도킹 계산 완료</span></div><div><strong>${view.state==='not_requested'?'미요청':count(view.redockingCount)}</strong><span>관측 리간드 재도킹 계산 완료</span></div></div>${view.state==='not_requested'?'':`<p>${e(view.redockingQc)}${typeof rmsd==='number'&&Number.isFinite(rmsd)&&rmsd>=0?' · RMSD '+rmsd.toLocaleString('ko-KR',{maximumSignificantDigits:8})+' Å':''}${typeof threshold==='number'&&Number.isFinite(threshold)&&threshold>=0?' · 기준 '+threshold.toLocaleString('ko-KR',{maximumSignificantDigits:8})+' Å':''}</p><p>관측 리간드 재도킹은 원래 결합 자세를 재현하는 확인 계산입니다. 후보 도킹 수에 합치지 않습니다. 후보 도킹 계산 완료도 자세 QC 통과나 실험적 결합·효능을 뜻하지 않습니다.</p>`}${dimensions.length?'<dl>'+dimensions.map(([label,value])=>`<dt>${e(label)}</dt><dd>${number(value)}개</dd>`).join('')+'</dl>':''}${typeof receipt?.stop_reason==='string'?'<p>후보 도킹 종료 사유: '+e(Object.hasOwn(stops,receipt.stop_reason)?stops[receipt.stop_reason]:receipt.stop_reason)+'</p>':''}${view.holds.length?items(view.holds.map(code=>Object.hasOwn(DOCKING_HOLDS,code)?DOCKING_HOLDS[code]:'추가 확인이 필요한 처리 코드: '+code),''):''}${view.state!=='not_requested'&&view.selectionReasons.length?'<p>대표 후보 선정에서 기록된 보류·평가 대기 사유입니다. 도킹 대상 선정 사유와는 별도이며, 한 후보에 여러 사유가 있을 수 있습니다.</p>'+items(view.selectionReasons.map(({state,code,count})=>(state==='held'?'보류':'평가 대기')+' · '+(Object.hasOwn(SELECTION_REASONS,code)?SELECTION_REASONS[code]:'처리 코드: '+code)+' · '+number(count)+'개'),''):''}${context?.protonation_scope||receipt?.preparation?.protonation_scope?'<p>단백질은 원래 구조·템플릿에 기록된 상태로 준비합니다. 입력 pH에 따른 단백질 양성자화는 미평가입니다.</p>':''}<details><summary>구조 준비·재도킹·pH 범위의 원자료</summary><pre>${e(JSON.stringify({protein_receptor_context:context,bound_complex_docking:receipt},null,2))}</pre></details></section>`;
}
const TOXICITY_STATUS={success:'계산 완료',partial_nonfinite:'일부 값만 계산됨',failed:'계산 실패',failed_invalid_probability:'출력값 검사 실패',not_requested_cached_prediction:'기존 예측이 있어 추가 모델은 실행하지 않음',unassessed_model_not_configured:'미평가 · 이 연구를 제출할 당시 추가 모델이 등록되지 않았습니다.',not_configured:'미평가 · 이 실행에는 추가 모델이 설정되지 않았습니다.',unsupported_structure:'추가 모델에서 지원하지 않는 구조',unassessed_model_candidate_limit:'미평가 · 이번 실행의 후보 수 한도',unassessed_model_stage_held:'미평가 · 추가 모델 계산 보류',unassessed_model_profile_held:'미평가 · 이 실행의 모델 설정을 확인하지 못해 계산 보류',configured:'설정됨 · 실제 계산 결과 없음',held:'미평가 · 추가 계산 확인 필요'};
const TOXICITY_ENDPOINTS={AMES:'에임스 시험 관련 예측',hERG:'심장 이온채널 관련 예측',DILI:'약물 유발 간손상 관련 예측',ClinTox:'임상 독성 관련 예측',Carcinogens_Lagunin:'발암성 관련 예측',Skin_Reaction:'피부 반응 관련 예측',SkinReaction:'피부 반응 관련 예측',LD50_Zhu:'급성 독성 지표 예측',molecular_weight:'분자량',logP:'분배 계수 logP',hydrogen_bond_acceptors:'수소 결합 받개 수',hydrogen_bond_donors:'수소 결합 주개 수',tpsa:'극성 표면적',QED:'약물 유사성 지표',Lipinski:'Lipinski 규칙 지표',stereo_centers:'입체 중심 수'};
const TOXICITY_CLASSIFIERS=new Set(['AMES','hERG','DILI','ClinTox','Carcinogens_Lagunin','Skin_Reaction','SkinReaction']);
const DERIVED_PROPERTIES=new Set(['molecular_weight','logP','hydrogen_bond_acceptors','hydrogen_bond_donors','tpsa','QED','Lipinski','stereo_centers']);
function toxicityDefinitions(prediction){
 const definitions=new Map(),fields=new Set(['id','task_type','category','name','units']);
 if(!/^[0-9a-f]{64}$/.test(prediction?.endpoint_definitions_sha256||'')||!Array.isArray(prediction?.endpoint_definitions)||prediction.endpoint_definitions.length>200)return definitions;
 for(const item of prediction.endpoint_definitions){
  if(!item||typeof item!=='object'||Array.isArray(item)||Object.keys(item).some(key=>!fields.has(key))||typeof item.id!=='string'||!item.id||item.id.length>200||definitions.has(item.id))return new Map();
  if(Object.values(item).some(value=>value!==null&&(typeof value!=='string'||value.length>500||/[\x00-\x1f\x7f]/.test(value))))return new Map();
  definitions.set(item.id,item);
 }
 return definitions;
}
const IMPORTANT_TOXICITY=['Skin_Reaction','SkinReaction','AMES','DILI','hERG','ClinTox','Carcinogens_Lagunin','LD50_Zhu'];
function toxicityValuesMarkup(values,{derived=false,failed=false,usable=true,definitions=new Map(),strictDefinitions=false,channel='additional',candidateId='',qualityFlagged=new Set()}={}){
 const entries=values&&typeof values==='object'&&!Array.isArray(values)?Object.entries(values).slice(0,200):[];
 if(!entries.length)return'';
 const recognized=key=>{const definition=definitions.get(key);return definition?(derived?definition.category==='Physicochemical':definition.category!=='Physicochemical'&&['classification','regression'].includes(definition.task_type)):!strictDefinitions&&(derived?DERIVED_PROPERTIES:TOXICITY_CLASSIFIERS).has(key)};
 const state=(key,x)=>failed?'failed':!usable?'unassessed':x===null?'missing':typeof x!=='number'||!Number.isFinite(x)?'invalid':!recognized(key)?'undefined':!derived&&(definitions.get(key)?.task_type==='classification'||!definitions.has(key))&&(x<0||x>1)?'out_of_range':qualityFlagged.has(key)?'quality_warning':'valid';
 const value=(key,x)=>({failed:'실패한 출력 · 원자료 확인',unassessed:'현재 상태에서는 평가값으로 사용하지 않음',missing:'값 없음',invalid:'유효한 숫자 아님',undefined:'항목 정의 확인 필요',out_of_range:'출력 범위 오류'})[state(key,x)]||x.toLocaleString('ko-KR',{maximumSignificantDigits:6})+(state(key,x)==='quality_warning'?' · 품질 경고':'');
 const label=key=>{const known=TOXICITY_ENDPOINTS[key]||definitions.get(key)?.name;return known?known+' ('+key+')':key.replaceAll('_',' ')};
 const units=key=>definitions.get(key)?.units?.trim()||(!strictDefinitions&&!derived&&!definitions.has(key)&&TOXICITY_CLASSIFIERS.has(key)?'0–1 출력 척도':'단위 미기록');
 const task=key=>derived?'분자 물성':definitions.get(key)?.task_type==='classification'||!strictDefinitions&&!definitions.has(key)&&TOXICITY_CLASSIFIERS.has(key)?'분류 모델 · 0–1 원출력':definitions.get(key)?.task_type==='regression'?'회귀 모델 원출력':'지표 정의 미확인';
 const table=(items,view='all')=>`<div class="research-value-table" data-output-view="${e(view)}"><table><thead><tr><th>항목</th><th>${derived?'계산값':'모델 원출력'}</th><th>원래 단위·척도</th><th>값의 종류</th></tr></thead><tbody>${items.map(([key,x])=>`<tr data-endpoint="${e(key)}" data-channel="${e(channel)}" data-candidate-id="${e(candidateId)}" data-value-status="${e(state(key,x))}"><th scope="row">${e(label(key))}</th><td>${e(value(key,x))}</td><td>${e(units(key))}</td><td>${e(task(key))}</td></tr>`).join('')}</tbody></table></div>`;
 const valid=entries.filter(([key,x])=>state(key,x)==='valid'),important=IMPORTANT_TOXICITY.flatMap(key=>valid.filter(([id])=>id===key));
 const preview=!derived&&!failed?(important.length?important:valid.slice(0,4)):[];
 const title=derived?(channel==='cached'?'기존 DB에 함께 저장된 분자 물성':'별도로 계산한 분자 물성'):failed?'실패한 계산의 원출력 · 유효한 평가로 사용하지 않음':'모델이 예측한 항목과 값';
 return`${preview.length?'<h5>먼저 살펴볼 모델 원출력</h5>'+table(preview,'key'):''}<details><summary>${e(title)} (${entries.length}개)</summary>${table(entries)}${derived?'<p>이 물성값은 학습된 독성 모델의 예측값과 구분합니다.</p>':'<p>각 항목의 원래 출력 척도입니다. 서로 다른 항목의 수치를 같은 위험도로 비교하지 마세요.</p>'}<p>정의나 범위를 확인하지 못한 수치는 원자료에 보존하며 여기에서는 표시하지 않습니다.</p></details>`;
}
function cachedDefinitions(prediction){
 const hash=value=>typeof value==='string'&&/^[0-9a-f]{64}$/.test(value),refs=prediction?.metadata_refs;
 if(!hash(prediction?.endpoint_definitions_sha256)||!hash(prediction?.model_manifest_sha256)||refs?.endpoint_definitions_sha256?.sha256!==prediction.endpoint_definitions_sha256||refs?.model_manifest_sha256?.sha256!==prediction.model_manifest_sha256||!Array.isArray(prediction.endpoint_definitions)||prediction.endpoint_definitions.length>200)return new Map();
 const definitions=new Map();
 for(const entry of prediction.endpoint_definitions){
  if(!entry||typeof entry!=='object'||Array.isArray(entry)||typeof entry.id!=='string'||!entry.id||entry.id.length>200||definitions.has(entry.id)||!['classification','regression'].includes(entry.task_type))return new Map();
  const selected=Object.fromEntries(['id','task_type','category','name','units'].map(key=>[key,entry[key]]));
  if(Object.values(selected).some(value=>value!==null&&(typeof value!=='string'||value.length>500||/[\x00-\x1f\x7f]/.test(value))))return new Map();
  definitions.set(entry.id,selected);
 }
 return definitions;
}
function validatedRadarMarkup({values,properties,definitions,source,usable,blockedKeys=[]}){
 if(!usable||!definitions.size)return'';
 const allowedProperties=[...definitions].filter(([,d])=>d.category==='Physicochemical').map(([key])=>key);
 const allowedToxicity=[...definitions].filter(([,d])=>d.category!=='Physicochemical'&&d.task_type==='classification').map(([key])=>key);
 return`<div class="metric-radar-grid">${metricRadarMarkup({kind:'properties',values:properties,source,allowedKeys:allowedProperties,blockedKeys})}${metricRadarMarkup({kind:'toxicity',values,source,allowedKeys:allowedToxicity,blockedKeys})}</div>`;
}
export function cachedToxicityMarkup(predictions,{candidateId=''}={}){
 const records=Array.isArray(predictions)?predictions:predictions&&typeof predictions==='object'?[predictions]:[];
 if(!records.length)return'<section class="research-options research-evaluation-channel"><h4>기존 DB 독성·ADMET 예측</h4><p>기존 DB의 독성·ADMET 예측이 없습니다. 추가 모델 예측은 아래의 별도 채널에서 확인하세요.</p></section>';
 return`<section class="research-options research-evaluation-channel"><h4>기존 DB 독성·ADMET 예측</h4>${records.map(prediction=>{
  const definitions=cachedDefinitions(prediction),state=prediction?.status,failed=!['success','complete','partial_nonfinite'].includes(state),model=prediction?.model_manifest||{},values=prediction?.values,qualityFlagged=new Set(Object.values(prediction?.quality_flags&&typeof prediction.quality_flags==='object'?prediction.quality_flags:{}).flatMap(value=>Array.isArray(value)?value:[]).filter(value=>typeof value==='string'));
  const unsupported=typeof state==='string'&&(state.startsWith('unsupported_')||state==='outside_operational_size_limit');
  const unsupportedReason=state==='unsupported_multi_component'?'성분이 둘 이상인 구조':state==='unsupported_element'?'원모델에서 지원하지 않는 원소':state==='unsupported_query_or_polymer'?'질의 구조 또는 고분자 표기':state==='outside_operational_size_limit'?'운영 크기 제한 초과':'원모델의 지원 범위 밖 구조';
  const statusText=unsupported?'기존 DB 모델에서 미지원으로 기록된 구조입니다. 사유: '+unsupportedReason+'. 추가 모델의 평가는 아래 별도 채널에서 확인하세요.':failed?'기존 예측이 실패했거나 처리 상태가 미확인입니다.':state==='partial_nonfinite'?'기존 예측의 일부 값만 유효합니다.':'기존 DB에 저장된 모델 예측 기록입니다.';
  const learned=Object.fromEntries(Object.entries(values&&typeof values==='object'&&!Array.isArray(values)?values:{}).filter(([key])=>definitions.get(key)?.category!=='Physicochemical'));
  const derived=Object.fromEntries(Object.entries(values&&typeof values==='object'&&!Array.isArray(values)?values:{}).filter(([key])=>definitions.get(key)?.category==='Physicochemical'));
  return`<p>${e(statusText)}</p>${model.name||model.version?`<p>원모델: ${e(model.name||'이름 미기록')} · 버전 ${e(model.version||'미기록')}</p>`:''}${qualityFlagged.size?'<p>저장된 출력에 품질 경고가 있습니다. 해당 지표의 원출력을 확인하고 유효한 평가로 단정하지 마세요.</p>':''}${definitions.size?'':'<p>원모델·지표 정의의 연결 근거를 확인하지 못해 수치를 기본 표에 표시하지 않습니다. 원출력은 상세 기록에 보존됩니다.</p>'}${validatedRadarMarkup({values:learned,properties:derived,definitions,source:'기존 DB · '+(prediction.model_id||model.name||'모델 미기록'),usable:!failed,blockedKeys:[...qualityFlagged]})}${toxicityValuesMarkup(learned,{definitions,strictDefinitions:true,failed,channel:'cached',candidateId,qualityFlagged})}${toxicityValuesMarkup(derived,{definitions,strictDefinitions:true,derived:true,failed,channel:'cached',candidateId,qualityFlagged})}<details><summary>기존 DB 예측의 모델·입력·품질 경고·원자료</summary><pre>${e(JSON.stringify(prediction,null,2))}</pre></details>`;
 }).join('')}<p>모델 출력은 6개 유효숫자로 표시합니다. 분류 지표는 0–1 점수이며, 전체 자릿수·지표 정의·단위는 원자료와 Excel에 있습니다.</p></section>`;
}
export function additionalToxicityMarkup(prediction,{candidateId=prediction?.candidate_id||''}={}){
 const model=prediction?.model_identity||{},state=prediction?.status;
 const status=prediction&&Object.hasOwn(TOXICITY_STATUS,state)?TOXICITY_STATUS[state]:prediction?'처리 상태 확인 필요':'추가 모델 결과 미기록';
 const failed=['failed','failed_invalid_probability'].includes(state),usable=['success','partial_nonfinite'].includes(state),definitions=toxicityDefinitions(prediction);
 const version=typeof model.version==='string'&&model.version?`<p>모델: ${e(model.provider==='admet_ai'?'ADMET-AI':model.provider||'이름 미기록')} · 버전 ${e(model.version)}</p>`:'';
 return`<section class="research-options research-evaluation-channel"><h4>추가 CPU 모델의 독성·ADMET 예측</h4><p>${e(status)}</p>${version}${failed?'<p>실패한 원출력은 성공한 평가나 안전성 판단에 사용하지 않습니다.</p>':''}${usable&&!definitions.size?'<p>지표 정의를 확인하지 못해 수치를 기본 표에 표시하지 않습니다. 원출력은 상세 기록에 보존됩니다.</p>':''}${validatedRadarMarkup({values:prediction?.values,properties:prediction?.physicochemical_values,definitions,source:'추가 모델 · '+(model.provider||'모델 미기록')+' · '+(model.version||'버전 미기록'),usable})}${toxicityValuesMarkup(prediction?.values,{failed,usable,definitions,strictDefinitions:true,channel:'additional',candidateId})}${toxicityValuesMarkup(prediction?.physicochemical_values,{derived:true,failed,usable,definitions,strictDefinitions:true,channel:'additional',candidateId})}<p>모델 출력은 6개 유효숫자로 표시합니다. 분류 지표는 0–1 점수이며, 전체 자릿수·지표 정의·단위는 원자료와 Excel에 있습니다.</p><p>추가 모델 결과는 순위·제외 조건을 자동으로 바꾸지 않습니다. 이 화합물에 대한 모델 적용 범위는 확인되지 않았습니다. 구조 경보·기존 DB 예측과 구분해 검토합니다.</p></section>`;
}
function structureImageMarkup(row,label){
 const png=row?.structure_image;
 let inline='';
 if(typeof png==='string'&&png.length<=512*1024&&/^data:image\/png;base64,iVBORw0KGgo[A-Za-z0-9+/]*={0,2}$/.test(png))try{
  const encoded=png.slice('data:image/png;base64,'.length),bytes=atob(encoded),integer=offset=>[...bytes.slice(offset,offset+4)].reduce((sum,value)=>sum*256+value.charCodeAt(0),0);
  if(encoded.length%4===0&&bytes.length>45&&bytes.slice(12,16)==='IHDR'&&bytes.slice(-8,-4)==='IEND'&&integer(16)>0&&integer(16)<=2048&&integer(20)>0&&integer(20)<=2048)inline=png;
 }catch{}
 const image=inline||(typeof row?.image_url==='string'?safeImage(row.image_url):'');
 return image?`<figure class="research-result-molecule"><img src="${e(image)}" alt="${e(label)}의 2차원 분자 구조" loading="lazy"><figcaption>2차원 분자 구조</figcaption></figure>`:'<p class="research-muted">2차원 구조 그림 미기록 · 정확한 SMILES는 식별자·구조·원자료에서 확인하세요.</p>';
}
export function skinEvaluationMarkup(row,profile){
 const skin=row.skin||{},physical=[typeof row.mw==='number'&&Number.isFinite(row.mw)?'분자량 '+number(row.mw)+' Da':null,typeof row.logp==='number'&&Number.isFinite(row.logp)?'logP '+number(row.logp):null].filter(Boolean),channels=[['genotoxicity','유전독성 관련 구조 경보'],['skin_irritation','피부 자극 관련 구조 경보'],['phototoxicity','광독성 관련 구조 경보']];
 return`<section class="research-options research-evaluation-channel"><h4>피부 특성·구조 경보</h4><p>${physical.length?e(physical.join(' · '))+' · ':''}피부 특성 참고점수 ${number(row.topical_score??skin.topical_score)}</p><p>종합 구조 경보: ${e(severity(row.hazard_severity||profile.hazard_severity))}</p><dl class="research-evaluation-facts">${channels.map(([key,label])=>`<dt>${e(label)}</dt><dd>${e(severity(profile[key]?.concern))}</dd>`).join('')}</dl><p>구조와 물성에서 계산한 참고 지표입니다. 미기록·미평가를 0이나 안전으로 해석하지 않습니다. 투과가 높을수록 모든 피부 목적에 유리한 것은 아닙니다.</p></section>`;
}
export function candidateName(candidate,index=null){
 const suffix=Number.isInteger(index)?' '+(index+1):'';
 const name=candidate.name;
 if(typeof name==='string'&&/^gen-[a-f0-9]{32}$/i.test(name))return'생성 후보'+suffix+' · 이름 미기록';
 if(typeof name==='string'&&name&&!/^[a-f0-9]{64}$/i.test(name)&&name!==candidate.smiles)return name;
 if(candidate.discovery_origin?.route==='query_seed')return'출발 화합물';
 return'후보'+suffix+' · 이름 미기록';
}
function originLabel(origin){const route=typeof origin==='string'?origin:origin?.route;return({general_catalog:'DB 후보 탐색',general_catalog_diversity:'일반 DB 다양성 후보',generated_proposal:'선택적으로 생성한 새 구조',query_seed:'입력한 출발 화합물',observed_complex_redock_reference:'관측된 결합 자세의 재도킹 대조',frozen_measured_reference:'실측 자료가 있는 참조 화합물'})[route]||'탐색 경로는 상세 기록에서 확인해 주세요.'}
export function seedAssessmentMarkup(assessment){
 if(!assessment)return'';
 const row=assessment.row||assessment,profile=row.structural_admet||{};
 return`<section class="research-seed-assessment"><h3>출발 화합물의 평가</h3>${structureImageMarkup(row,'출발 화합물')}<p>연구의 출발 구조이며 새로운 후보나 실측 활성 대조군으로 세지 않습니다.</p>${skinEvaluationMarkup(row,profile)}${cachedToxicityMarkup(row.admet,{candidateId:row.candidate_id||''})}${additionalToxicityMarkup(row.additional_toxicity,{candidateId:row.candidate_id||''})}<details><summary>출발 화합물의 식별자·평가 근거·원자료</summary><pre>${e(JSON.stringify(assessment,null,2))}</pre></details></section>`;
}
function primaryCountsMarkup(summary){
 const counts=summary.counts||{},valid=value=>typeof value==='number'&&Number.isFinite(value)&&value>=0;
 if(summary.analysis_mode==='assessment')return`<div class="research-counts">${valid(counts.query_seeds)?`<div><strong>${number(counts.query_seeds)}</strong><span>입력 물질</span></div>`:''}<div><strong>평가 전용</strong><span>새 후보 탐색 요청 없음</span></div></div><p>입력 물질 평가 모드입니다. 새 후보 탐색은 요청하지 않았습니다. 수용체 도킹도 요청하지 않았습니다. 처리 수는 아래에서 확인할 수 있습니다.</p>${!valid(counts.query_seeds)?'<p>처리 수를 집계하고 있습니다. 집계가 끝나면 후보 수를 표시합니다.</p>':''}`;
 const population=counts.candidate_population,passed=counts.filter_pass,excluded=valid(population)&&valid(passed)&&passed<=population?population-passed:null;
 const cards=[['새 고유 후보',counts.new_unique_candidates],['출발 구조',counts.query_seeds],['전체 후보 행',population],['필터 통과 후보 행',passed],['필터 제외 후보 행',excluded],['실측 자료 참조',counts.reference_controls??counts.measured_references],['선택 생성 구조',counts.generated],['전체 저장 관측',counts.saved_rows]];
 const visible=cards.filter(([,value])=>valid(value));
 return`<div class="research-counts">${visible.map(([label,value])=>`<div><strong>${number(value)}</strong><span>${e(label)}</span></div>`).join('')||'<p>처리 수를 집계하고 있습니다. 집계가 끝나면 후보 수를 표시합니다.</p>'}</div><p class="research-muted">새 고유 후보는 출발 구조와 같은 분자를 제외한 집계입니다. 후보 행·출발 구조·실측 참조·전체 저장 관측은 서로 다른 집계이며 합쳐서 새 후보 수로 세지 않습니다. 선택 생성 구조는 새 후보에 포함될 수 있습니다.</p>`;
}
export function catalogSelectionMarkup(selection){
 if(selection?.fallback_used!==true)return'';
 const complete=selection.status==='complete',held=selection.status==='held',valid=value=>typeof value==='number'&&Number.isInteger(value)&&value>=0;
 const title=complete?'계획한 일반 DB 다양성 부분집합 선택을 완료했습니다':held?'일반 DB의 구조 다양성 부분집합 탐색을 보류했습니다':'일반 DB 부분집합의 선정 상태를 확인해 주세요';
 const counts=[['이번 검색의 인덱스 모집단',selection.source_population],['구조 다양성 비교 풀',selection.source_pool_size],['확보한 후보',selection.selected_count],['요청 후보 한도',selection.requested_limit],['비교 풀 한도',selection.pool_limit],['난수 시드',selection.seed]].filter(([,value])=>valid(value));
 return`<section class="notice research-catalog-selection" role="status"><div><h3>${e(title)}</h3><p>이 연구에서는 알려진 활성 화합물 기준 대신 일반 DB에서 화학 구조의 다양성을 고려해 제한된 후보를 고릅니다.${complete&&valid(selection.selected_count)?' 확보한 후보는 '+number(selection.selected_count)+'개입니다.':held?' 이번에는 유효한 후보 부분집합을 확보하지 못했습니다.':''}</p><p><strong>현재 단백질의 활성·결합은 이 후보 선택 단계에서 미평가입니다.</strong> 선정된 부분집합의 피부·독성 평가와 후속 도킹 결과를 확인하세요.</p><details><summary>DB 후보 선정 범위·방법·제한</summary>${counts.length?'<dl>'+counts.map(([label,value])=>'<dt>'+e(label)+'</dt><dd>'+number(value)+'</dd>').join('')+'</dl>':''}${items(selection.limitations,'구체적인 선정 제한은 아래 원자료에서 확인해 주세요.')}<pre>${e(JSON.stringify(selection,null,2))}</pre></details></div></section>`;
}
export function researchSummaryMarkup(summary){
 const input=summary.input||{},context=summary.skin_context||{},candidates=Array.isArray(summary.candidates)?summary.candidates:[];
 const counts=Object.entries(summary.counts||{}).filter(([,v])=>typeof v==='number'&&Number.isFinite(v)),assessment=summary.analysis_mode==='assessment';
 const selection=summary.catalog_selection;
 const limitations=Array.isArray(summary.limitations)?summary.limitations.map(value=>value==='카탈로그 확보가 제한적이거나 미완료입니다. 전체 카탈로그를 탐색했다고 해석하지 않습니다.'?(selection?.fallback_used===true&&selection.status==='complete'?'계획한 일반 DB 다양성 부분집합 선택을 완료했습니다. 전체 DB 전수 평가는 아닙니다.':'카탈로그의 후보·상세 근거 중 일부를 확보하지 못했습니다. 검색 범위와 미확보 항목은 원자료에서 확인하세요.'):value):summary.limitations;
 const mode=input.kind==='protein'&&selection?.fallback_used===true?'단백질 입력 · 일반 DB 후보의 구조·피부 평가':input.kind==='bound_complex'&&!assessment?'결합 구조에서 시작 · 같은 포켓의 골격 비교':input.kind==='protein'?'확인한 단백질의 표적 기반 탐색':MODES[summary.analysis_mode]||summary.analysis_mode||'입력에 따른 탐색';
 const cancellation=summary.cancellation,summaryStatus=cancellation?.requested===true?(cancellation.state==='cleanup_unverified'?'작업 자원 정리 확인 필요':cancellation.state==='stopped'&&cancellation.resources_released===true?'취소 완료 · 작업 자원 정리 확인됨':'취소 요청됨 · 종료·자원 정리 확인 중'):(STATUS[summary.status]||summary.status||'상태 확인 중');
 return`<div class="page-heading"><div><span class="eyebrow">연구 결과</span><h2>${e(identityLabel(input))}</h2><p>${e(mode)} · ${e(summaryStatus)}</p></div></div>
 ${researchExecutionScopeMarkup(summary)}${primaryCountsMarkup(summary)}${dockingContextMarkup(summary)}${routineSelectionMarkup(summary.routine_selection,candidates.length)}${counts.length?`<details><summary>전체 처리 수와 역할별 집계</summary><dl>${counts.map(([k,v])=>`<dt>${e(COUNTS[k]||k)}</dt><dd>${number(v)}</dd>`).join('')}</dl></details>`:''}
 <p>${context.product?`${e(PRODUCTS[context.product]||context.product)} · `:'제품 형태 미지정 · '}${context.applied_percent!==null&&context.applied_percent!==undefined?`배합 ${number(context.applied_percent)}%`:'배합 농도 미지정'} · ${e(GOALS[context.goal]||'탐색 목적')}</p>
 ${catalogSelectionMarkup(selection)}
 ${seedAssessmentMarkup(summary.seed_assessment)}
 ${assessment?'<p>이번 실행은 입력 화합물의 평가만 요청했습니다. 새로운 후보는 탐색하지 않았습니다.</p>':'<h3>먼저 검토할 후보</h3><p>각 후보의 선정 이유·피부 구조 경보·기존 DB 예측·추가 CPU 모델을 구분해 확인해 주세요. 전체 후보와 제외 이유는 아래의 전체 후보 표와 Excel에서 확인할 수 있습니다.</p>'}
 <p class="research-muted">추가 모델 결과는 기존 예측·구조 경보와 구분해 검토하며 순위·제외 조건을 자동으로 바꾸지 않습니다. 모델 출력에 임의 위험 등급이나 안전 확률을 부여하지 않습니다.</p>
 <div class="research-candidate-grid">${candidates.map((c,index)=>`<article class="research-candidate" data-candidate-id="${e(c.candidate_id)}"><h4>${e(candidateName(c,index))}</h4><p class="research-muted">${e(originLabel(c.discovery_origin))}</p>${structureImageMarkup(c,candidateName(c,index))}<strong>선정 이유</strong>${items(c.selection_reasons,'선정 이유가 아직 기록되지 않았습니다.')}<strong>살펴볼 경보</strong>${items(c.concerns,'기록된 경보가 없습니다. 안전성이 확인은 해당 계산 조건에 한정됩니다.')}<strong>아직 모르는 점</strong>${items(c.missing,'추가 원자료와 실험 검토가 필요합니다.')}${skinEvaluationMarkup(c,c.admet||{})}${cachedToxicityMarkup(c.cached_toxicity,{candidateId:c.candidate_id||''})}${additionalToxicityMarkup(c.additional_toxicity,{candidateId:c.candidate_id||''})}<details><summary>식별자·구조·계산 근거</summary><code class="research-smiles">${e(c.smiles||'구조 미기록')}</code><pre>${e(JSON.stringify(c,null,2))}</pre></details></article>`).join('')||(assessment?'':'<div class="notice"><p>현재 표시할 후보가 없습니다. 아래 처리 제한과 다음 할 일을 확인해 주세요.</p></div>')}</div>
 <div class="research-options"><h3>이 결과를 어떻게 해석하나요?</h3><p>${e(summary.scientific_status?.message||summary.scientific_status||'선정 이유·계산 조건·미평가 항목을 읽고 실험 검토 목록을 정하세요.')}</p>${items(limitations,'추가 제한은 세부 근거와 원자료를 확인해 주세요.')}<h3>다음 할 일</h3>${items(summary.next_actions,'선정 이유와 미평가 항목을 검토한 뒤 후보의 시험 조건을 정해 주세요.')}</div>
 ${summary.ai?.summary?`<details><summary>선택적 연구 도우미 설명</summary><p>${e(summary.ai.summary)}</p><p>도우미 상태: ${e(summary.ai.status||'미기록')}</p></details>`:''}<details><summary>결과 요약의 전체 근거</summary><pre>${e(JSON.stringify(summary,null,2))}</pre></details>`;
}

function safeImage(url){try{const resolved=new URL(url,globalThis.location?.href||'https://skinpipe.invalid');return['http:','https:'].includes(resolved.protocol)&&!resolved.username&&!resolved.password?resolved.href:''}catch{return''}}
