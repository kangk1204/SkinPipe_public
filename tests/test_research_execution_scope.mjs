import assert from 'node:assert/strict';
import test from 'node:test';
import {researchSummaryMarkup} from '../web/public/research-results.js';

// Synthetic contracts only: no private research IDs, structures, or results.
const summary=()=>({schema:'skinpipe-novice-research-summary-v1',input:{kind:'protein',name:'Synthetic'},analysis_mode:'protein_target_screening',status:'complete',execution:{catalog_limit:50,screening_policy:{target_candidates:3},compute_budget_seconds:1800},counts:{collected:50,filter_pass:44,docking_selected:3,docking_completed:3,query_redocking_completed:1},bound_complex_docking:{status:'complete',planned_candidate_count:3,new_docking_count:1,reused_count:2,stop_reason:'target_reached'},routine_selection:{selected_count:3,shortfall:0}});
const panel=s=>{const html=researchSummaryMarkup(s),match=html.match(/<section class="research-options research-execution-scope">[\s\S]*?<\/section>/);assert.ok(match,'The saved request and actual scope must be visible together');return match[0]};
const value=(html,label)=>{const match=html.match(new RegExp('<dt>'+label+'</dt><dd>(.*?)</dd>'));assert.ok(match,'Missing scope field: '+label);return match[1]};

test('saved 50/3 request is shown with separate collected, planned, completed, reused and representative counts',()=>{
 const s=summary(),before=JSON.stringify(s),html=panel(s);
 assert.equal(value(html,'공통 카탈로그 요청 한도'),'50개');assert.equal(value(html,'대표 후보 목표'),'3개');assert.equal(value(html,'전체 연구 최대 처리 시간'),'30분 (1,800초)');
 for(const [label,count]of [['확보한 카탈로그 행',50],['검색 필터 통과 후보 행',44],['도킹 계획 후보',3],['후보 도킹 계산 완료',3],['이번 연구의 신규 도킹 계산',1],['동일 조건의 저장 도킹 재사용',2],['대표 선정 결과',3]])assert.equal(value(html,label),count+'개');
 assert.match(html,/요청 한도와 대표 목표는 보장되는 결과 수가 아닙니다/);assert.match(html,/검색 필터 통과는 도킹 자세 QC 통과나 대표 선정을 뜻하지 않습니다/);
 assert.match(html,/관측 리간드 재도킹은 후보 도킹 수에 합치지 않습니다/);assert.equal(JSON.stringify(s),before);
});
test('1000/200 and larger saved requests use their own values and exact time, not fixed examples',()=>{
 for(const [limit,target,seconds,time]of [[1000,200,7200,'2시간 (7,200초)'],[2505,80,61,'61초']]){
  const s=summary();s.execution={catalog_limit:limit,screening_policy:{target_candidates:target},compute_budget_seconds:seconds};
  const html=panel(s);assert.equal(value(html,'공통 카탈로그 요청 한도'),limit.toLocaleString('ko-KR')+'개');assert.equal(value(html,'대표 후보 목표'),target+'개');assert.equal(value(html,'전체 연구 최대 처리 시간'),time);
 }
});
test('explicit zero processing counts remain zero while zero request bounds are invalid',()=>{
 const s=summary();s.execution={catalog_limit:0,screening_policy:{target_candidates:0},compute_budget_seconds:0};s.counts={collected:0,filter_pass:0,docking_completed:0};s.bound_complex_docking={planned_candidate_count:0,new_docking_count:0,reused_count:0};s.routine_selection={selected_count:0,shortfall:3};
 const html=panel(s);for(const label of ['공통 카탈로그 요청 한도','대표 후보 목표','전체 연구 최대 처리 시간'])assert.equal(value(html,label),'미기록');
 for(const label of ['확보한 카탈로그 행','검색 필터 통과 후보 행','도킹 계획 후보','후보 도킹 계산 완료','이번 연구의 신규 도킹 계산','동일 조건의 저장 도킹 재사용','대표 선정 결과'])assert.equal(value(html,label),'0개');
});
test('requested 1000/200 with zero acquired rows does not imply 1000 candidates failed QC',()=>{
 const s=summary();s.execution={catalog_limit:1000,screening_policy:{target_candidates:200},compute_budget_seconds:7200};s.counts={collected:0,filter_pass:0,docking_completed:0};s.bound_complex_docking={planned_candidate_count:0,new_docking_count:0,reused_count:0,stop_reason:'eligible_pool_exhausted'};s.routine_selection={selected_count:0,shortfall:200};
 const html=panel(s);assert.equal(value(html,'공통 카탈로그 요청 한도'),'1,000개');assert.equal(value(html,'대표 후보 목표'),'200개');assert.equal(value(html,'확보한 카탈로그 행'),'0개');assert.equal(value(html,'도킹 계획 후보'),'0개');assert.equal(value(html,'대표 선정 결과'),'0개');
 assert.doesNotMatch(html,/모든 후보.*QC|1,000개.*평가 완료|카탈로그 연결 실패|source_unavailable/,'Do not infer an acquisition error or QC exclusion from zero counts alone');
});
test('missing or malformed saved values never become zero or current 200/1000 defaults',()=>{
 for(const invalid of [undefined,null,'50',true,-1,0.5,NaN,Infinity,Number.MAX_SAFE_INTEGER+1,{},[]]){
  const s=summary();s.execution={catalog_limit:invalid,screening_policy:{target_candidates:invalid},compute_budget_seconds:invalid};s.counts={collected:invalid,filter_pass:invalid,docking_completed:invalid,representative_goal:200};s.bound_complex_docking={planned_candidate_count:invalid,new_docking_count:invalid,reused_count:invalid};s.routine_selection={selected_count:invalid,policy:{target_candidates:200},shortfall:invalid};
  const html=panel(s);assert.equal((html.match(/<dd>미기록<\/dd>/g)||[]).length,10);assert.doesNotMatch(html,/<dd>0개|<dd>200개|<dd>1,000개|목표보다/);
 }
 const html=panel({counts:{representative_goal:200}});assert.equal(value(html,'대표 후보 목표'),'미기록');assert.equal(value(html,'공통 카탈로그 요청 한도'),'미기록');
});
test('assessment does not present inactive stored search or docking bounds as requested computation',()=>{
 const s=summary();s.analysis_mode='assessment';
 const html=panel(s);assert.match(html,/카탈로그 후보 탐색·대표 후보 선정·수용체 도킹은 요청하지 않았습니다/);assert.equal(value(html,'전체 연구 최대 처리 시간'),'30분 (1,800초)');assert.doesNotMatch(html,/<dt>공통 카탈로그|<dt>대표 후보|<dt>도킹 계획|<dt>후보 도킹/);
});
test('explicit ligand-only contract says docking not requested; missing evidence remains unknown',()=>{
 const s={...summary(),input:{kind:'compound'},analysis_mode:'similarity',protein_receptor_context:null,bound_complex_docking:null,counts:{collected:50,docking_selected:0,docking_completed:0,query_redocking_completed:0}};
 let html=panel(s);for(const label of ['도킹 계획 후보','후보 도킹 계산 완료','이번 연구의 신규 도킹 계산','동일 조건의 저장 도킹 재사용'])assert.equal(value(html,label),'미요청');
 delete s.protein_receptor_context;s.counts={collected:50};html=panel(s);assert.equal(value(html,'후보 도킹 계산 완료'),'미기록');assert.doesNotMatch(html,/<dd>미요청/);
});
test('catalog limit describes common-catalog scope without asserting it truncated a saved source population',()=>{
 const html=panel(summary());assert.match(html,/공통 카탈로그에서 새로 가져올 때 적용/);assert.match(html,/기존 연구를 후보 모집단으로 선택했다면 저장된 전체 모집단을 보존/);assert.match(html,/이 요약에 모집단 출처가 없으면 한도의 실제 적용 여부를 단정하지 않습니다/);
});
test('a stored shortfall and stop reason are explained without reconstructing absent eligibility reasons',()=>{
 const s=summary();s.routine_selection={selected_count:2,shortfall:1};s.bound_complex_docking.stop_reason='eligible_pool_exhausted';
 let html=panel(s);assert.match(html,/저장된 대표 선정 결과는 목표보다 1개 적습니다/);assert.match(html,/조건에 맞는 추가 후보가 없어/);assert.match(html,/보류·평가 대기 사유는 아래 대표 선정 내역/);
 s.bound_complex_docking.stop_reason='budget_deferred';assert.match(panel(s),/남은 계산 시간 예산으로 다음 도킹 묶음을 시작할 수 없어/);
 s.bound_complex_docking.stop_reason='<img src=x onerror=alert(1)>';html=panel(s);assert.match(html,/&lt;img/);assert.doesNotMatch(html,/<img/);
 delete s.bound_complex_docking.stop_reason;html=panel(s);assert.doesNotMatch(html,/시간 예산|추가 후보가 없어|중복 때문에/);
});
test('scope is a single visible panel above candidate totals and long per-compound details',()=>{
 const html=researchSummaryMarkup(summary());assert.equal((html.match(/class="research-options research-execution-scope"/g)||[]).length,1);assert.ok(html.indexOf('research-execution-scope')<html.indexOf('research-counts'));assert.ok(html.indexOf('research-execution-scope')<html.indexOf('research-candidate-grid'));
});
