import assert from 'node:assert/strict';
import test from 'node:test';
import {renderDiscovery} from '../web/public/discovery.js';
import {researchIdentityMarkup} from '../web/public/research-inputs.js';

const deferred=()=>{let resolve,reject;const promise=new Promise((yes,no)=>{resolve=yes;reject=no});return{promise,resolve,reject}};
const ready=(name='A')=>({export_id:'a'.repeat(32),status:'complete',result:{filename:`SkinPipe-${name}-xlsx.zip`}});
const flush=async()=>{for(let i=0;i<8;i++)await Promise.resolve()};

// A small DOM fixture models containment and replacement of action/history
// nodes, so a delayed response cannot pass by reusing an obsolete anchor.
function harness(){
 const nodes=new Map(),history=[],links=[];
 const make=(id='',parent=null)=>({
  id,parent,isConnected:true,_html:'',_text:'',value:'',checked:false,disabled:false,hidden:false,dataset:{},attributes:{},
  get innerHTML(){return this._html},
  set innerHTML(value){this._html=String(value);this._text=String(value).replace(/<[^>]*>/g,'');for(const node of [...nodes.values(),...links])if(node.parent===this)node.remove();if(this.id==='discovery-list')buildHistory(value)},
  get textContent(){return this._text},set textContent(value){this._text=String(value);this._html=String(value)},
  setAttribute(key,value){this.attributes[key]=value},removeAttribute(key){delete this.attributes[key]},
  querySelector(selector){if(selector==='[data-run-export-status]')return this.output;return null},
  closest(selector){for(let node=this;node;node=node.parent)if(selector==='.research-run'&&node.historyRow||selector==='#research-result-actions'&&node.id==='research-result-actions')return node;return null},
  remove(){this.isConnected=false;for(const node of [...nodes.values(),...links])if(node.parent===this)node.remove()},
  focus(){},scrollIntoView(){}
 });
 const get=selector=>{
  if(!nodes.get(selector)?.isConnected){
   const id=selector.startsWith('#')?selector.slice(1):selector;
   const parent=['research-result-excel','research-result-zip','research-result-refresh'].includes(id)?get('#research-result-actions'):null;
   nodes.set(selector,make(id,parent));
  }
  return nodes.get(selector);
 };
 function buildHistory(html){
  for(const row of history){row.remove();row.view.remove();row.export.remove();row.output.remove()}
  history.length=0;
  for(const [,index]of String(html).matchAll(/data-view-run="(\d+)"/g)){
   const row=make();row.historyRow=true;row.view=make('',row);row.view.dataset.viewRun=index;row.export=make('',row);row.export.dataset.exportRun=index;row.export.textContent='Excel 받기';row.output=make('',row);history.push(row);
  }
 }
 const kinds=['compound','protein','bound_complex'].map(value=>Object.assign(make(),{value}));
 const modes=['similarity','scaffold_hopping','assessment'].map(value=>Object.assign(make(),{value,checked:value==='similarity'}));
 const main={innerHTML:'',contains:node=>Boolean(node?.isConnected),querySelector:selector=>selector==='input[name="research-mode"]:checked'?modes.find(mode=>mode.checked):get(selector),querySelectorAll(selector){
  if(selector==='input[name="research-kind"]')return kinds;if(selector==='input[name="research-mode"]')return modes;
  if(selector==='[data-view-run]')return history.map(row=>row.view);if(selector==='[data-export-run]')return history.map(row=>row.export);return[];
 }};
 return{main,get,history,links,make};
}
async function fixture(exportResponse=()=>ready(),{native=true,download,status}={}){
 const h=harness(),calls=[];
 const record=id=>({run_id:id,status:'complete',snapshot_id:'snapshot-'+id,input:{value:id},counts:{},candidates:[]});
 const api=async(path,options)=>{
  const body=options?.body?JSON.parse(options.body):null;calls.push({path,body});
  if(path==='/v2/snapshots?view=summary')return{snapshots:[]};
  if(path==='/v2/research-capabilities')return{};
  if(path==='/v2/discovery-runs')return{runs:['A','B'].map(record)};
  if(path.endsWith('/summary'))return record(path.split('/')[3]);
  if(path==='/v2/exports')return exportResponse(body);
  if(path.endsWith('/download'))return download();
  if(path.startsWith('/v2/exports/'))return status();
  throw Error('Unexpected fixture request '+path);
 };
 if(native)api.downloadExport=({anchor,filename})=>{const link=h.make('',anchor.parent);link.filename=filename;h.links.push(link);return link};
 await renderDiscovery(h.main,api);h.get('#research-step-4').hidden=true;
 return{...h,calls,api,view:async index=>h.history[index].view.onclick(),exportResult:()=>h.get('#research-result-excel').onclick()};
}
function assertNoNotice(h){for(const id of['research-export-status','research-result-export-status']){assert.equal(h.get('#'+id).textContent,'');assert.equal(h.get('#'+id).attributes['aria-busy'],'false')}}

test('opening a verified saved run clears an old lost-submission acknowledgement',async()=>{
 const h=await fixture(),notice=h.get('#research-submit-status');
 notice.textContent='접수 응답을 확인하지 못했습니다. 접수 여부를 먼저 확인해 주세요.';notice.hidden=false;
 await h.view(0);
 assert.equal(notice.hidden,true);assert.equal(notice.textContent,'');
 assert.match(h.get('#research-result-summary').innerHTML,/>A</);
 assert.equal(h.calls.filter(call=>call.body&&call.path==='/v2/discovery-runs').length,0,'Opening the saved run never submits it again');
});

test('new input clears a finished result export notice and ready link',async()=>{
 const h=await fixture();await h.view(0);await h.exportResult();assert.match(h.get('#research-result-export-status').textContent,/SkinPipe-A/);assert.equal(h.links[0].isConnected,true);
 h.get('#research-new').onclick();assertNoNotice(h);assert.equal(h.links[0].isConnected,false);
});
test('opening another saved result clears both notices and the old ready link',async()=>{
 const h=await fixture();await h.view(0);await h.exportResult();await h.view(1);assertNoNotice(h);assert.equal(h.links[0].isConnected,false);assert.match(h.get('#research-result-summary').innerHTML,/>B</);
});
for(const transition of['new-input','different-run','same-run-refresh'])test('late result export cannot mutate '+transition,async()=>{
 const pending=deferred(),h=await fixture(()=>pending.promise);await h.view(0);const button=h.get('#research-result-excel'),request=h.exportResult();
 if(transition==='new-input')h.get('#research-new').onclick();else if(transition==='different-run')await h.view(1);else await h.get('#research-result-refresh').onclick();
 pending.resolve(ready());await request;assertNoNotice(h);assert.equal(h.links.length,0);assert.equal(button.disabled,false);
});
test('late export failure cannot put an error on another result',async()=>{
 const pending=deferred(),h=await fixture(()=>pending.promise);await h.view(0);const request=h.exportResult();await h.view(1);pending.reject(Error('Old export failure'));await request;
 assertNoNotice(h);assert.equal(h.get('.platform-errors').textContent,'');
});
test('late completed-export metadata cannot attach a download to a different result',async()=>{
 const pending=deferred(),h=await fixture(()=>({export_id:'a'.repeat(32),status:'complete'}),{status:()=>pending.promise});await h.view(0);const request=h.exportResult();await flush();await h.view(1);pending.resolve(ready());await request;
 assertNoNotice(h);assert.equal(h.links.length,0);
});
test('export polling stops when its response arrives after a result switch',async()=>{
 const originalTimer=globalThis.setTimeout,wait=deferred(),pending=deferred();let polls=0;
 try{
  globalThis.setTimeout=(callback,delay,...args)=>delay===1500?wait.promise.then(callback):originalTimer(callback,delay,...args);
  const h=await fixture(()=>({export_id:'a'.repeat(32),status:'queued'}),{status:()=>{polls++;return pending.promise}});await h.view(0);const request=h.exportResult();await flush();wait.resolve();await flush();assert.equal(polls,1);
  await h.view(1);pending.resolve({export_id:'a'.repeat(32),status:'running'});await request;assertNoNotice(h);assert.equal(polls,1);assert.equal(h.links.length,0);
 }finally{globalThis.setTimeout=originalTimer}
});
test('queued exports distinguish execution waiting from running file preparation',async()=>{
 const originalTimer=globalThis.setTimeout,waits=[];let polls=0;
 try{
  globalThis.setTimeout=(callback,delay,...args)=>delay===1500?(waits.push(callback),0):originalTimer(callback,delay,...args);
  const h=await fixture(()=>({export_id:'a'.repeat(32),status:'queued'}),{status:()=>++polls===1?{export_id:'a'.repeat(32),status:'running'}:ready()});
  await h.view(0);const button=h.get('#research-result-excel'),request=h.exportResult();await flush();
  assert.equal(button.textContent,'내보내기 실행 대기…');
  for(const id of['research-export-status','research-result-export-status']){
   const output=h.get('#'+id);assert.equal(output.textContent,'내보내기 실행을 기다리고 있습니다. 서버의 앞선 작업이 끝나면 준비를 시작합니다.');assert.equal(output.attributes['aria-busy'],'true');
  }
  assert.equal(waits.length,1);waits.shift()();await flush();
  assert.equal(button.textContent,'내보내기 준비 중…');assert.equal(h.get('#research-result-export-status').textContent,'서버에서 내보내기 파일을 준비하고 있습니다.');assert.equal(h.get('#research-result-export-status').attributes['aria-busy'],'true');
  assert.equal(waits.length,1);waits.shift()();await request;
  assert.equal(h.get('#research-result-export-status').dataset.state,'ready');assert.equal(h.get('#research-result-export-status').attributes['aria-busy'],'false');assert.equal(button.disabled,false);assert.equal(polls,2);assert.equal(h.links.length,1);
 }finally{globalThis.setTimeout=originalTimer}
});
test('a stale download response is cancelled without handing bytes to the browser',async()=>{
 const pending=deferred();let cancelled=0;
 const h=await fixture(()=>ready(),{native:false,download:()=>pending.promise});await h.view(0);const request=h.exportResult();await flush();assert.ok(h.calls.some(call=>call.path.endsWith('/download')));
 await h.view(1);pending.resolve({body:{cancel:async()=>{cancelled++}}});await request;assertNoNotice(h);assert.equal(cancelled,1);
});
test('history export can finish in its own row after switching results',async()=>{
 const pending=deferred(),h=await fixture(()=>pending.promise);await h.view(0);const row=h.history[0],request=row.export.onclick();await h.view(1);pending.resolve(ready());await request;
 assertNoNotice(h);assert.match(row.output.textContent,/SkinPipe-A/);assert.equal(h.links.length,1);assert.equal(h.links[0].parent,row);
 h.get('#research-new').onclick();assert.equal(h.links[0].isConnected,true,'A history download remains available in its own row');
});
test('history refresh invalidates detached export anchors',async()=>{
 const pending=deferred(),h=await fixture(()=>pending.promise),request=h.history[0].export.onclick();await h.get('#discovery-refresh').onclick();pending.resolve(ready());await request;
 assert.equal(h.links.length,0);
});
test('a history export never annotates the displayed result',async()=>{
 const h=await fixture();await h.view(1);await h.history[0].export.onclick();assert.equal(h.get('#research-result-export-status').textContent,'');assert.match(h.history[0].output.textContent,/SkinPipe-A/);
});
test('older concurrent history completion cannot replace the latest global notice',async()=>{
 const a=deferred(),b=deferred(),h=await fixture(body=>body.run_id==='A'?a.promise:b.promise);
 const first=h.history[0].export.onclick(),second=h.history[1].export.onclick();b.resolve(ready('B'));await second;a.resolve(ready('A'));await first;
 assert.match(h.get('#research-export-status').textContent,/SkinPipe-B/);assert.match(h.history[0].output.textContent,/SkinPipe-A/);assert.match(h.history[1].output.textContent,/SkinPipe-B/);
});
test('changing input while a save reservation is open aborts its writer before POST',async()=>{
 const originalPicker=globalThis.showSaveFilePicker,picker=deferred();let aborted=0;
 try{
  globalThis.showSaveFilePicker=()=>picker.promise;const h=await fixture(()=>ready(),{native:false});await h.view(0);const request=h.exportResult();await flush();h.get('#research-new').onclick();
  picker.resolve({createWritable:async()=>({abort:async()=>{aborted++}})});await request;
  assert.equal(aborted,1);assert.equal(h.calls.some(call=>call.path==='/v2/exports'),false);assertNoNotice(h);
 }finally{globalThis.showSaveFilePicker=originalPicker}
});
test('generation control names only the BRICS limit, with a separate substitution explanation',async()=>{
 const h=await fixture();assert.match(h.main.innerHTML,/BRICS 재조합 상한<input id="research-products"/);assert.match(h.main.innerHTML,/작용기 치환 구조는 별도로 추가되어 총 생성 수가 BRICS 재조합 상한보다 많을 수 있습니다/);
});
test('uploaded assembly labels do not claim a public source and preserve raw provenance',()=>{
 for(const [assembly_id,label]of[['1','구조에 기록된 조립체 1'],['asymmetric_unit','구조의 비대칭 단위'],['as_uploaded','업로드한 좌표 그대로']]){
  const result={kind:'bound_complex',status:'verified',input:{kind:'bound_complex',format:'upload',value:'synthetic-upload'},selection:{assembly_id,protein_chains:['A']},source_receipts:[{source:'upload',sha256:'b'.repeat(64)}]},before=JSON.stringify(result);
  const html=researchIdentityMarkup(result);assert.ok(html.includes(label));assert.doesNotMatch(html,/공개 구조의/);assert.equal(JSON.stringify(result),before);assert.ok(html.includes('synthetic-upload'));
 }
});
