import assert from 'node:assert/strict';
import test from 'node:test';
import {researchProgressMarkup} from '../web/public/research-results.js';
import fs from 'node:fs';
import {researchProgressView,applyResearchProgressState,RESEARCH_PROGRESS_STALE_MS} from '../web/public/research-progress.js';
import {renderDiscovery} from '../web/public/discovery.js';

for(const state of ['submitting','checking','queued','preparing','running','cancelling'])test(state+' indicates its own active operation',()=>{
 const view=researchProgressView(state,{checkedAt:100,now:101});
 assert.equal(view.active,true);assert.equal(view.busy,true);
 assert.equal(view.showProgress,['running','preparing'].includes(state));
});
for(const [state,tone]of [['complete','complete'],['completed','complete'],['partial','attention'],['failed','failed'],['cancelled','cancelled'],['cancellation_unverified','attention'],['cleanup_unverified','attention'],['needs_data','attention'],['paused','attention'],['held','attention']])test(state+' has a distinct static terminal presentation',()=>{
 const view=researchProgressView(state,{checkedAt:100,now:100+RESEARCH_PROGRESS_STALE_MS});
 assert.equal(view.tone,tone);assert.equal(view.active,false);assert.equal(view.busy,false);assert.equal(view.showProgress,false);
});
test('connection loss and aged observations cannot animate or present old counts as current',()=>{
 for(const options of [{connection:'unavailable'},{checkedAt:100,now:20100},{checkedAt:NaN,now:100},{checkedAt:100,now:99}]){
  const view=researchProgressView('running',options);assert.equal(view.tone,'unknown');assert.equal(view.active,false);assert.equal(view.showProgress,false);
 }
 assert.equal(researchProgressView('running',{checkedAt:100,now:20099}).active,true);
 assert.equal(researchProgressView('complete',{connection:'unavailable'}).tone,'unknown');
});
test('accessible busy state follows actual presentation and reduced motion is supported',()=>{
 const attrs={},node={setAttribute:(key,value)=>{attrs[key]=value}};
 applyResearchProgressState(node,researchProgressView('running'));assert.equal(attrs['aria-busy'],'true');assert.equal(attrs['data-active'],'true');
 applyResearchProgressState(node,researchProgressView('complete'));assert.equal(attrs['aria-busy'],'false');assert.equal(attrs['data-active'],'false');assert.equal(attrs['data-state'],'complete');
 const css=fs.readFileSync(new URL('../web/public/styles.css',import.meta.url),'utf8');
 assert.match(css,/@media\(prefers-reduced-motion:reduce\)[^\n]+animation:none/);
});
function harness(){
 const nodes=new Map();
 const make=()=>({_html:'',_text:'',value:'',checked:false,disabled:false,hidden:false,dataset:{},attributes:{},get innerHTML(){return this._html},set innerHTML(v){this._html=v;this._text=v.replace(/<[^>]*>/g,'')},get textContent(){return this._text},set textContent(v){this._text=String(v);this._html=String(v)},setAttribute(key,value){this.attributes[key]=value},removeAttribute(key){delete this.attributes[key]},focus(){},scrollIntoView(){}});
 const get=selector=>{if(!nodes.has(selector)){const node=make();node.querySelector=child=>child==='h2'&&node.innerHTML.includes('<h2>')?get(selector+' h2'):null;nodes.set(selector,node)}return nodes.get(selector)};
 const kinds=['compound','protein','bound_complex'].map(value=>Object.assign(make(),{value}));
 const modes=['similarity','scaffold_hopping','assessment'].map(value=>Object.assign(make(),{value,checked:value==='similarity'}));
 const main={innerHTML:'',querySelector(selector){return selector==='input[name="research-mode"]:checked'?modes.find(mode=>mode.checked):get(selector)},querySelectorAll(selector){if(selector==='input[name="research-kind"]')return kinds;if(selector==='input[name="research-mode"]')return modes;if(selector==='[data-view-run]')return[get('#view')];return[]}};
 return{main,get};
}
test('actual discovery renderer transitions running → lost/stale → recovered → completed without new submissions',async()=>{
 const originalTimer=globalThis.setTimeout,originalClear=globalThis.clearTimeout,timers=[];
 globalThis.setTimeout=(callback,delay)=>{const token={callback,delay,cancelled:false,unref(){}};timers.push(token);return token};
 globalThis.clearTimeout=token=>{token.cancelled=true};
 const next=delay=>{const token=timers.findLast(item=>!item.cancelled&&item.delay===delay);assert.ok(token,'Expected live timer '+delay);token.cancelled=true;return token.callback()};
 const fixture=harness(),calls=[];let fail=false,state='running',progressPayload={stage:'catalog_acquisition',acquired:2,processed:99,total:5};
 const record=()=>({run_id:'fixture',status:state,input:{value:'Synthetic public input'},counts:{},candidates:[],progress:progressPayload,...(state==='complete'?{snapshot_id:'synthetic-snapshot'}:{})});
 const api=async(path,options)=>{
  calls.push({path,method:options?.method||'GET'});
  if(path==='/v2/snapshots?view=summary')return{snapshots:[]};
  if(path==='/v2/research-capabilities')return{};
  if(path==='/v2/discovery-runs?view=status')return{runs:[record()]};
  if(path.endsWith('/summary'))return record();
  if(path==='/v2/discovery-runs/fixture?view=status'){if(fail)throw Error('synthetic connection loss');return record()}
  throw Error('Unexpected fixture request '+path);
 };
 try{
  await renderDiscovery(fixture.main,api);fixture.get('#view').dataset.viewRun='0';await fixture.get('#view').onclick();
  const panel=fixture.get('#research-live-status');assert.equal(panel.attributes['data-state'],'running');assert.equal(panel.attributes['data-active'],'true');assert.match(panel.innerHTML,/카탈로그 후보 확보.*카탈로그 행 2개 확보/);
  progressPayload={stage:'routine_docking',completed:9,total:32,processed:8,selected:14,qualified:11,target:200};await next(5000);assert.match(panel.innerHTML,/현재 도킹 묶음 9개 \/ 32개 처리 완료/);
  progressPayload={...progressPayload,completed:24};await next(5000);assert.match(panel.innerHTML,/현재 도킹 묶음 24개 \/ 32개 처리 완료/);assert.doesNotMatch(panel.innerHTML,/이 단계에서 8개|선택 14개/);
  progressPayload={stage:'routine_docking',completed:4,total:32,processed:72,qualified:75,target:200,selected:14};await next(5000);assert.match(panel.innerHTML,/현재 도킹 묶음 4개 \/ 32개 처리 완료 · 자세 QC 전 · 최근 QC 후 대표 후보 75개 \/ 목표 200개/);assert.doesNotMatch(panel.innerHTML,/72개|선택 14개/);
  fail=true;await next(5000);assert.equal(panel.attributes['data-state'],'unknown');assert.equal(panel.attributes['aria-busy'],'false');assert.match(panel.textContent,/마지막으로 확인한 상태/);assert.doesNotMatch(panel.innerHTML,/현재 단계|research-progress-icon/);
  fail=false;await next(10000);assert.equal(panel.attributes['data-active'],'true');await next(RESEARCH_PROGRESS_STALE_MS);assert.equal(panel.attributes['data-active'],'false');assert.match(panel.textContent,/최근 20초/);
  await next(5000);assert.equal(panel.attributes['data-active'],'true');
  state='complete';await next(5000);assert.equal(panel.attributes['data-state'],'complete');assert.equal(panel.attributes['aria-busy'],'false');assert.match(panel.innerHTML,/✓.*연구 완료/);assert.doesNotMatch(panel.innerHTML,/현재 단계|2개 \/ 5개/);assert.match(fixture.get('#research-result-actions').innerHTML,/Excel/);
  assert.equal(timers.filter(token=>!token.cancelled).length,0);assert.ok(calls.every(call=>call.method==='GET'),'Progress never submits or restarts a study');
 }finally{globalThis.setTimeout=originalTimer;globalThis.clearTimeout=originalClear}
});

// Stage payloads are merged by the worker; old fields must not acquire a new meaning.
test('current docking batch advances using completed and never previous QC processed or toxicity selected',()=>{
 for(const stage of ['docking','routine_docking','bound_complex_docking','cpu_fallback'])for(const completed of [0,9,24,32]){
  const progress={stage,completed,total:32,processed:8,selected:14,acquired:1000,qualified:11,target:200},before=JSON.stringify(progress),html=researchProgressMarkup(progress);
  assert.match(html,new RegExp('현재 도킹 묶음 '+completed+'개 / 32개 처리 완료 · 자세 QC 전'));
  assert.match(html,/최근 QC 후 대표 후보 11개 \/ 목표 200개/);
  assert.doesNotMatch(html,/이 단계에서 8개|선택 14개|1,000|QC 통과|도킹 성공|%|완료 예정/);
  assert.equal(JSON.stringify(progress),before);
 }
});
test('pose quality reports cumulative attempts and its own qualified goal without a stale batch denominator',()=>{
 const html=researchProgressMarkup({stage:'routine_pose_quality',processed:40,qualified:30,target:200,total:32,completed:32,selected:14,acquired:1000});
 assert.match(html,/누적 도킹 시도 40개/);assert.match(html,/최근 QC 후 대표 후보 30개 \/ 목표 200개/);
 assert.doesNotMatch(html,/40개 \/ 32개|현재 도킹 묶음|선택 14개|1,000|40개.*계산 완료/);
});
test('additional toxicity selection is not a completion counter and ignores merged earlier counters',()=>{
 const html=researchProgressMarkup({stage:'additional_toxicity',selected:14,processed:1001,total:1001,completed:32,qualified:11,target:200});
 assert.match(html,/이번 단계 선택 14개 · 완료 수는 아직 확정하지 않았습니다/);assert.doesNotMatch(html,/1,001|32개|대표 후보|14개 처리|100%/);
});
test('catalog acquisition owns acquired and never uses a prior processed or total',()=>{
 const html=researchProgressMarkup({stage:'catalog_acquisition',acquired:1000,processed:8,total:32,selected:14,completed:24,qualified:11,target:200});
 assert.match(html,/카탈로그 행 1,000개 확보/);assert.doesNotMatch(html,/8개|32개|14개|24개|대표 후보| \/ /);
});
test('shape comparison and ligand preparation use requested while cheap assessment owns total',()=>{
 for(const stage of ['reference_shape','ligand_preparation']){
  const html=researchProgressMarkup({stage,processed:2,requested:3,total:1001,complete:1,acquired:1000});assert.match(html,/이 단계에서 2개 \/ 3개 처리/);assert.doesNotMatch(html,/1,001|1,000/);
 }
 const html=researchProgressMarkup({stage:'cheap_structure_assessment',processed:101,total:1001,requested:3,selected:14,qualified:11,target:200});
 assert.match(html,/이 단계에서 101개 \/ 1,001개 처리/);assert.doesNotMatch(html,/ \/ 3개|선택 14개|대표 후보/);
});
test('population freeze only reports its catalog count and never old stage completion',()=>{
 const html=researchProgressMarkup({stage:'candidate_population_frozen',catalog_count:1000,processed:40,total:32,completed:32,qualified:30,target:200});
 assert.match(html,/확보한 카탈로그 행 1,000개/);assert.doesNotMatch(html,/40개|32개|대표 후보/);
});
test('counter-free and unknown stages cannot inherit old counts or qualified totals',()=>{
 for(const stage of ['identity_resolution','query_seed_frozen','reference_acquisition','literature_evidence','resource_wait','budget_deferred','failed','future-stage','__proto__','constructor']){
  const html=researchProgressMarkup({stage,processed:8,total:32,completed:24,requested:50,acquired:1000,selected:14,qualified:11,target:200,catalog_count:1000});
  assert.doesNotMatch(html,/\d|\[object Object\]|function Object/);
 }
});
test('invalid or missing counters never fall back to merged numbers; valid zero stays zero',()=>{
 const invalid=[undefined,null,'24',true,-1,0.5,NaN,Infinity,Number.MAX_SAFE_INTEGER+1];
 for(const completed of invalid){const html=researchProgressMarkup({stage:'routine_docking',completed,total:32,processed:8,acquired:1000});assert.doesNotMatch(html,/\d/);assert.match(html,/처리 수 확인 필요/)}
 for(const selected of invalid)assert.doesNotMatch(researchProgressMarkup({stage:'additional_toxicity',selected,processed:8,total:32}),/\d/);
 assert.match(researchProgressMarkup({stage:'routine_docking',completed:0,total:32}),/0개 \/ 32개 처리 완료/);
 assert.match(researchProgressMarkup({stage:'routine_docking',completed:9}),/현재 도킹 묶음 9개 처리 완료/);
 assert.doesNotMatch(researchProgressMarkup({stage:'routine_docking',completed:40,total:32}),/\d/);
 assert.doesNotMatch(researchProgressMarkup({stage:'routine_pose_quality',qualified:201,target:200}),/최근 QC 후 대표 후보/);
 assert.doesNotMatch(researchProgressMarkup({stage:'routine_pose_quality',qualified:0,target:0}),/최근 QC 후 대표 후보/);
 for(const value of [null,[],true,'stage'])assert.equal(researchProgressMarkup(value),'');
});

test('latest QC representatives remain historical while the current docking batch resets and advances',()=>{
 for(const stage of ['docking','routine_docking','bound_complex_docking','cpu_fallback'])for(const completed of [32,0,4,9]){
  const progress={stage,completed,total:32,processed:72,qualified:75,target:200,selected:14,acquired:1000},before=JSON.stringify(progress),html=researchProgressMarkup(progress);
  assert.match(html,new RegExp('현재 도킹 묶음 '+completed+'개 / 32개 처리 완료 · 자세 QC 전 · 최근 QC 후 대표 후보 75개 / 목표 200개'));
  assert.doesNotMatch(html,/72개|14개|1,000|75개.*QC 통과|현재 묶음.*75개.*통과/);
  assert.equal(JSON.stringify(progress),before);
 }
});
test('historical QC figures require a valid recorded pair without inventing missing values',()=>{
 const invalid=[undefined,null,'75',true,-1,0.5,NaN,Infinity,Number.MAX_SAFE_INTEGER+1];
 for(const stage of ['docking','routine_docking','bound_complex_docking','cpu_fallback','routine_pose_quality']){
  for(const qualified of invalid)assert.doesNotMatch(researchProgressMarkup({stage,completed:4,total:32,qualified,target:200,processed:72,selected:14}),/최근 QC 후 대표 후보/);
  for(const target of [...invalid,0,74])assert.doesNotMatch(researchProgressMarkup({stage,completed:4,total:32,qualified:75,target,processed:72,selected:14}),/최근 QC 후 대표 후보/);
  assert.match(researchProgressMarkup({stage,completed:4,total:32,qualified:0,target:200}),/최근 QC 후 대표 후보 0개 \/ 목표 200개/);
  assert.match(researchProgressMarkup({stage,completed:4,total:32,qualified:200,target:200}),/최근 QC 후 대표 후보 200개 \/ 목표 200개/);
 }
});

for(const published of [false,true])test('budget-deferred guidance preserves records and only opens a separate input flow; published='+published,async()=>{
 const h=harness(),calls=[],record={run_id:'budget-fixture',status:'budget_deferred',input:{kind:'protein',value:'Synthetic target'},execution:{compute_budget_seconds:7200},progress:{stage:'budget_deferred',processed:72,qualified:199,target:200},counts:{representative_selected:199},candidates:[],...(published?{snapshot_id:'synthetic-published'}:{})},saved=JSON.stringify(record);
 const api=async(path,options)=>{calls.push({path,method:options?.method||'GET'});if(path==='/v2/snapshots?view=summary')return{snapshots:[]};if(path==='/v2/research-capabilities')return{};if(path==='/v2/discovery-runs?view=status')return{runs:[record]};if(path.endsWith('/summary'))return record;throw Error('Unexpected request '+path)};
 await renderDiscovery(h.main,api);h.get('#view').dataset.viewRun='0';await h.get('#view').onclick();
 const controls=h.get('#research-run-controls').innerHTML;
 assert.match(controls,/저장된 계산·체크포인트와 입력·조건은 보존/);assert.match(controls,/완료된 상태가 아닙니다/);assert.match(controls,/기존 연구의 시간 한도를 늘리거나 이어서 실행하지 않습니다/);
 assert.match(controls,/별도 연구/);assert.match(controls,/2 → 4시간, 최대 8시간/);assert.match(controls,/자동 복원하지 않습니다/);assert.match(controls,/서버가 동일한 계산 조건을 검증한 경우에만 재사용/);assert.match(controls,/보장하지 않습니다/);
 assert.doesNotMatch(controls,/research-resume-run|research-cancel-run/);
 if(published){assert.match(controls,/이미 등록된 결과/);assert.match(h.get('#research-result-actions').innerHTML,/research-result-excel/);assert.doesNotMatch(controls,/아직 등록되지 않아/)}
 else{assert.match(controls,/최종 결과 스냅샷이 아직 등록되지 않아/);assert.match(controls,/체크포인트는 최종 결과가 아닙니다/);assert.doesNotMatch(h.get('#research-result-actions').innerHTML,/research-result-excel|research-result-zip|research-open-population/)}
 for(const[id,value]of Object.entries({'discovery-hours':'2','research-catalog-limit':'1000','research-representatives':'200','discovery-source':'','discovery-pdb':'1ABC'}))h.get('#'+id).value=value;
 h.get('#discovery-critical').checked=true;const callCount=calls.length;h.get('#research-budget-new').onclick();
 assert.equal(h.get('#research-step-1').hidden,false);assert.equal(h.get('#research-step-4').hidden,true);assert.equal(h.get('#research-confirm').disabled,true);
 assert.equal(h.get('#discovery-hours').value,'2','Navigation never increases the old or new form budget automatically');assert.equal(h.get('#research-catalog-limit').value,'1000');assert.equal(h.get('#research-representatives').value,'200');assert.equal(h.get('#discovery-source').value,'');assert.equal(h.get('#discovery-pdb').value,'1ABC');assert.equal(h.get('#discovery-critical').checked,true);
 assert.equal(calls.length,callCount,'Opening input preparation performs no request');assert.ok(calls.every(call=>call.method==='GET'));assert.equal(JSON.stringify(record),saved);
});
test('budget guidance is removed after switching studies and an obsolete button cannot navigate',async()=>{
 const h=harness(),calls=[];let state='budget_deferred';
 const record=()=>({run_id:'budget-switch',status:state,input:{value:'Synthetic'},counts:{},candidates:[],...(state==='complete'?{snapshot_id:'synthetic-published'}:{})});
 await renderDiscovery(h.main,async(path,options)=>{calls.push(options?.method||'GET');if(path==='/v2/snapshots?view=summary')return{snapshots:[]};if(path==='/v2/research-capabilities')return{};if(path==='/v2/discovery-runs?view=status')return{runs:[record()]};if(path.endsWith('/summary'))return record();throw Error('Unexpected request '+path)});
 h.get('#view').dataset.viewRun='0';await h.get('#view').onclick();const obsolete=h.get('#research-budget-new').onclick;
 state='complete';await h.get('#research-result-refresh').onclick();assert.doesNotMatch(h.get('#research-run-controls').innerHTML,/research-budget-new|시간 한도|체크포인트/);assert.match(h.get('#research-result-actions').innerHTML,/research-result-excel/);
 obsolete();assert.equal(h.get('#research-step-4').hidden,false);assert.ok(calls.every(method=>method==='GET'));
});

test('cleanup quarantine is an explicit hold even after a native checkpoint reaches its goal',async()=>{
 const h=harness(),calls=[],record={run_id:'quarantine-fixture',status:'needs_data',progress:{stage:'routine_pose_quality',qualified:200,target:200,resource_quarantine:{active:true,reason:'owned_resource_cleanup_unverified'}},counts:{},candidates:[]},saved=JSON.stringify(record);
 await renderDiscovery(h.main,async(path,options)=>{calls.push(options?.method||'GET');if(path==='/v2/snapshots?view=summary')return{snapshots:[]};if(path==='/v2/research-capabilities')return{};if(path==='/v2/discovery-runs?view=status')return{runs:[record]};if(path.endsWith('/summary'))return record;throw Error('Unexpected fixture request '+path)});
 h.get('#view').dataset.viewRun='0';await h.get('#view').onclick();
 const panel=h.get('#research-live-status'),controls=h.get('#research-run-controls').innerHTML;
 assert.equal(panel.attributes['data-state'],'attention');assert.equal(panel.attributes['aria-busy'],'false');assert.equal(panel.attributes['data-active'],'false');assert.match(panel.innerHTML,/자원 정리 확인 필요/);assert.doesNotMatch(panel.innerHTML,/연구 완료|현재 단계|needs_data/);
 assert.match(controls,/작업 자원 정리를 확인하지 못해 연구를 보류/);assert.match(controls,/최종 결과로 등록되지 않은 체크포인트/);assert.match(controls,/설치 담당자/);assert.match(controls,/자동으로 재실행하지 않습니다/);
 assert.doesNotMatch(controls,/research-resume-run|research-cancel-run|research-budget-new/);assert.doesNotMatch(h.get('#research-result-actions').innerHTML,/research-result-excel|research-result-zip|research-open-population/);
 assert.match(h.get('#research-result-summary').innerHTML.split('<section')[0],/추가 확인으로 보류/);assert.doesNotMatch(h.get('#research-result-summary').innerHTML.split('<section')[0],/needs_data/);assert.ok(calls.every(method=>method==='GET'));assert.equal(JSON.stringify(record),saved);
});
