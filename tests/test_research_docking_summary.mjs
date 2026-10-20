import assert from 'node:assert/strict';
import test from 'node:test';
import {researchDockingView,dockingContextMarkup,researchSummaryMarkup,researchRunCounts} from '../web/public/research-results.js';

// Synthetic summary contracts only; no saved research inputs or private IDs.
const compound=()=>({
 schema:'skinpipe-novice-research-summary-v1',input:{kind:'compound',name:'Synthetic'},
 analysis_mode:'scaffold_hopping',status:'complete',
 protein_receptor_context:null,bound_complex_docking:null,
 counts:{candidate_population:3,docking_selected:0,docking_completed:0,query_redocking_completed:0}
});
const bound=()=>({
 ...compound(),input:{kind:'bound_complex',name:'Synthetic'},
 counts:{candidate_population:3,docking_selected:0,docking_completed:0,query_redocking_completed:1},
 bound_complex_docking:{status:'complete',planned_candidate_count:0,completed_count:0,new_docking_count:0,reused_count:0,stop_reason:'eligible_pool_exhausted',redocking:{status:'complete',pass:true,rmsd_angstrom:1.99971092152089,threshold_angstrom:2}}
});
const visible=html=>html.replace(/<details>[\s\S]*?<\/details>/g,'');

test('explicit ligand-only summary says not requested while preserving raw zero counts',()=>{
 for(const mode of ['similarity','scaffold_hopping']){
  const summary={...compound(),analysis_mode:mode},before=JSON.stringify(summary);
  assert.equal(researchDockingView(summary).state,'not_requested');
  const html=dockingContextMarkup(summary);
  assert.match(visible(html),/도킹 미요청/);assert.doesNotMatch(visible(html),/<strong>0개<\/strong>/);
  assert.match(researchSummaryMarkup(summary),/&quot;docking_completed&quot;: 0/);
  assert.equal(JSON.stringify(summary),before);
 }
});
test('unknown, incomplete, contradictory, or future contracts never mean not requested',()=>{
 const base=compound();
 for(const summary of [
  {},{counts:{docking_completed:0}},
  {...base,schema:undefined},{...base,schema:'future-schema'},{...base,input:{}},
  {...base,protein_receptor_context:undefined},{...base,bound_complex_docking:undefined},
  {...base,counts:{...base.counts,docking_selected:1}},
  {...base,counts:{...base.counts,query_redocking_completed:1}},
  {...base,counts:{...base.counts,docking_completed:1}}
 ])assert.notEqual(researchDockingView(summary).state,'not_requested');
 assert.match(visible(dockingContextMarkup({})),/근거 미확인/);
 assert.doesNotMatch(visible(dockingContextMarkup({})),/<strong>0개<\/strong>/);
});
test('assessment keeps its dedicated evaluation-only scope',()=>{
 const summary={...compound(),analysis_mode:'assessment'};
 assert.equal(dockingContextMarkup(summary),'');
 assert.match(visible(researchSummaryMarkup(summary)),/수용체 도킹도 요청하지 않았습니다/);
});
test('missing receptor selection has an explicit hold, not a generic docking failure',()=>{
 const summary={...compound(),input:{kind:'protein'},protein_receptor_context:{status:'held',hold_reasons:['explicit_pdb_and_observed_ligand_selection_required']},bound_complex_docking:{status:'held',hold_reasons:['explicit_pdb_and_observed_ligand_selection_required']}};
 assert.equal(researchDockingView(summary).state,'held');
 const html=visible(dockingContextMarkup(summary));
 assert.match(html,/도킹 보류/);assert.match(html,/PDB 구조와 관측된 결합 화합물을 선택/);
 assert.equal((html.match(/PDB 구조와 관측된 결합 화합물을 선택/g)||[]).length,1);
 assert.doesNotMatch(html,/도킹 미요청|도킹 처리 실패/);
});
test('one successful reference redocking cannot become one candidate docking',()=>{
 const summary=bound(),view=researchDockingView(summary),html=visible(dockingContextMarkup(summary));
 assert.equal(view.state,'no_eligible_candidates');assert.equal(view.candidateCount,0);assert.equal(view.redockingCount,1);
 assert.match(html,/<strong>0개<\/strong><span>후보 도킹 계산 완료/);
 assert.match(html,/<strong>1개<\/strong><span>관측 리간드 재도킹 계산 완료/);
 assert.match(html,/관측 리간드 재도킹 QC 통과/);assert.match(html,/도킹 계획에 포함된 후보<\/dt><dd>0개/);
 assert.match(html,/RMSD 1\.9997109 Å · 기준 2 Å/,'Near-threshold RMSD must not round to the threshold');
 assert.match(html,/후보 도킹 수에 합치지 않습니다/);
});
test('a completed reference calculation can fail QC and remains separate from engine failure',()=>{
 const summary=bound();summary.bound_complex_docking.redocking.pass=false;
 summary.bound_complex_docking.hold_reasons=['query_redocking_geometry_failed'];
 const view=researchDockingView(summary),html=visible(dockingContextMarkup(summary));
 assert.equal(view.state,'qc_failed');assert.equal(view.redockingCount,1);
 assert.match(html,/재도킹 QC 불통과/);assert.doesNotMatch(html,/재도킹 QC 통과|재도킹 계산 실패/);
});
test('engine failure and no-candidate records have distinct explanations',()=>{
 const failed={...bound(),bound_complex_docking:{status:'failed',redocking:{status:'failed'}}};
 assert.equal(researchDockingView(failed).state,'failed');
 assert.match(visible(dockingContextMarkup(failed)),/관측 리간드 재도킹 계산 실패/);
 const empty={...bound(),counts:{candidate_population:0,docking_completed:0},bound_complex_docking:{status:'complete'}};
 assert.equal(researchDockingView(empty).state,'no_candidates');
 assert.match(visible(dockingContextMarkup(empty)),/후보 행이 0개로 기록/);
});
test('candidate calculations including reuse are not asserted to pass QC',()=>{
 const summary=bound();
 summary.counts.docking_completed=3;
 Object.assign(summary.bound_complex_docking,{completed_count:3,planned_candidate_count:3,new_docking_count:1,reused_count:2,qualified_count:0,stop_reason:'pose_quality_runtime_or_input_held'});
 summary.routine_selection={decisions:{a:{state:'held',reasons:['pose_QC_failed']},b:{state:'pending',reasons:['exact_pose_QC_required']}}};
 const view=researchDockingView(summary),html=visible(dockingContextMarkup(summary));
 assert.equal(view.state,'calculated');assert.equal(view.candidateCount,3);
 assert.match(html,/저장 결과 재사용이 포함/);assert.match(html,/동일 조건의 저장 도킹 재사용<\/dt><dd>2개/);
 assert.match(html,/자세 QC 불통과 · 1개/);assert.match(html,/평가 대기 · 해당 조건의 도킹 자세 QC 필요 · 1개/);
 assert.match(html,/QC의 실행 환경 또는 입력 확인이 필요/);
 assert.doesNotMatch(html,/후보 도킹 QC 통과|신규 후보 도킹 성공/);
});
test('selection holds and pending reasons stay separate from unrecorded docking selection causes',()=>{
 const summary=bound();
 summary.routine_selection={decisions:{
  a:{state:'held',reasons:['critical_warning','critical_warning','structural_alert_review_required']},
  b:{state:'held',reasons:['critical_warning']},
  c:{state:'pending',reasons:['exact_pose_QC_required']},
  reference:{state:'control',reasons:['critical_warning']},
  duplicate:{state:'duplicate',reasons:['critical_warning']}
 }};
 const view=researchDockingView(summary),html=visible(dockingContextMarkup(summary));
 assert.deepEqual(view.selectionReasons.find(reason=>reason.code==='critical_warning'),{state:'held',code:'critical_warning',count:2});
 assert.match(html,/도킹 대상 선정 사유와는 별도/);assert.match(html,/한 후보에 여러 사유/);
 assert.match(html,/평가 대기 · 해당 조건의 도킹 자세 QC 필요 · 1개/);
 assert.doesNotMatch(html,/대조군.*중복|모든 후보.*보류/,'Do not reconstruct an unrecorded backend gate from downstream pending state');
});
test('missing or invalid counts and unconfirmed QC remain unknown',()=>{
 for(const value of [undefined,null,'1',-1,0.5,NaN,Infinity,Number.MAX_SAFE_INTEGER+1]){
  const summary={...bound(),counts:{docking_completed:value,query_redocking_completed:value}};
  assert.equal(researchDockingView(summary).candidateCount,null);
  assert.equal(researchDockingView(summary).redockingCount,null);
  assert.match(visible(dockingContextMarkup(summary)),/<strong>미기록<\/strong>/);
 }
 for(const redocking of [{status:'complete'}, {pass:true}, {status:'running',pass:true}]){
  const summary=bound();summary.bound_complex_docking.redocking=redocking;
  assert.equal(researchDockingView(summary).redockingQc,'관측 리간드 재도킹 QC 미확인');
 }
});
test('unrecognized evidence is escaped and preserved, not translated into a known cause',()=>{
 const summary=bound(),code='<img src=x onerror=alert(1)>';
 summary.bound_complex_docking.hold_reasons=[code,'__proto__'];
 summary.bound_complex_docking.stop_reason=code;
 summary.routine_selection={decisions:{a:{state:'held',reasons:[code,'constructor']}}};
 const before=JSON.stringify(summary),html=dockingContextMarkup(summary);
 assert.match(html,/&lt;img/);assert.doesNotMatch(html,/<img|\[object Object\]|function Object/);
 assert.match(html,/추가 확인이 필요한 처리 코드/);assert.match(html,/처리 코드: constructor/);
 assert.equal(JSON.stringify(summary),before);
});
test('summary places the docking explanation before long seed assessments and keeps count roles',()=>{
 const summary=bound();summary.seed_assessment={row:{smiles:'CCO'}};
 const html=researchSummaryMarkup(summary);
 assert.ok(html.indexOf('research-docking-summary')<html.indexOf('research-seed-assessment'));
 assert.equal((html.match(/class="research-options research-docking-summary"/g)||[]).length,1);
 assert.doesNotMatch(html,/<span>도킹 완료<\/span>/);
 assert.equal(researchRunCounts({candidate_population:10}),'후보 행 10');
 assert.equal(researchRunCounts({docking_completed:0,query_redocking_completed:1}),'후보 도킹 계산 완료 0 · 관측 리간드 재도킹 계산 완료 1');
});
test('existing pH boundaries remain visible and malformed reason arrays do not crash',()=>{
 const summary={protein_receptor_context:{protonation_scope:'as written',hold_reasons:['requested_protein_ph_not_supported_by_as_written_cpu_preparation']},bound_complex_docking:{hold_reasons:'unexpected-string'},routine_selection:{decisions:{a:null,b:{state:'held',reasons:'unexpected-string'}}}};
 const html=dockingContextMarkup(summary);
 assert.match(html,/pH별 단백질 상태를 지원하지 않아 도킹을 보류/);
 assert.match(html,/입력 pH에 따른 단백질 양성자화는 미평가입니다/);
});
