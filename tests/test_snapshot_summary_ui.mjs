import assert from 'node:assert/strict';
import test from 'node:test';
import {renderPlatform} from '../web/public/platform.js';
import {renderDiscovery} from '../web/public/discovery.js';

// Synthetic compact API records deliberately omit candidate payloads and files.
const compact=()=>({
 snapshot_id:'synthetic-summary',status:'verified',candidate_count:52,
 created_at:'2026-10-20T00:00:00Z',target_counts:{'uniprot:SYNTHETIC':52},
 provenance:{marker:'Synthetic <target>',target_id:'uniprot:SYNTHETIC',search_mode:'similarity',
  primary_docking_protocols:[{target_id:'uniprot:SYNTHETIC',protocol_id:'synthetic-protocol',engine:'vina',seed:17,candidate_count:3}]}
});
function harness(){
 const nodes=new Map(),make=()=>({_html:'',textContent:'',value:'',checked:false,disabled:false,hidden:false,dataset:{},attributes:{},
  get innerHTML(){return this._html},set innerHTML(value){this._html=value},
  setAttribute(key,value){this.attributes[key]=value},removeAttribute(key){delete this.attributes[key]},
  addEventListener(){},focus(){},scrollIntoView(){}});
 const get=key=>{if(!nodes.has(key))nodes.set(key,make());return nodes.get(key)};
 const kinds=['compound','protein','bound_complex'].map(value=>Object.assign(make(),{value}));
 const modes=['similarity','scaffold_hopping','assessment'].map(value=>Object.assign(make(),{value,checked:value==='similarity'}));
 const select=selector=>selector==='input[name="research-mode"]:checked'?modes.find(mode=>mode.checked):get(selector);
 const all=selector=>selector==='input[name="research-kind"]'?kinds:selector==='input[name="research-mode"]'?modes:[];
 const main={innerHTML:'',querySelector:select,querySelectorAll:all},body=get('#platform-body');
 body.querySelector=select;body.querySelectorAll=all;
 return{main,body,get};
}
async function platform({snapshots=[compact()],selection={snapshot:'synthetic-summary',target:'uniprot:SYNTHETIC',protocol:'synthetic-protocol'},namespace='legacy',route='research',runs=[]}={}){
 const previous=globalThis.sessionStorage;
 globalThis.sessionStorage={getItem:()=>null,removeItem(){}};
 const h=harness(),calls=[],api=async(path,options)=>{
  const body=options?.body?JSON.parse(options.body):null;calls.push({path,body});
  if(path==='/v2/snapshots?view=summary')return{snapshots};
  if(path==='/v2/discovery-runs')return{runs};
  if(path==='/v2/catalog-status')return{compound_count:100,prediction_admitted_partitions:0};
  if(path==='/v2/candidate-queries')return{snapshot_id:body.snapshot_id,target_id:body.target_id,rows:[],next_offset:null,
   counts:{population:52,filtered:0,displayed:0,population_calculated:3,reasons:{critical:0,opposing:0},missing:{purchase:0,pchembl:0}}};
  throw Error('Unexpected initial UI request '+path);
 };
 try{await renderPlatform({main:h.main,api,selection,namespace,route})}
 finally{if(previous===undefined)delete globalThis.sessionStorage;else globalThis.sessionStorage=previous}
 return{...h,calls};
}
test('candidate comparison uses compact metadata with the exact target and recorded protocol',async()=>{
 const h=await platform(),query=h.calls.find(call=>call.path==='/v2/candidate-queries');
 assert.deepEqual(h.calls.map(call=>call.path),['/v2/snapshots?view=summary','/v2/candidate-queries']);
 assert.equal(query.body.snapshot_id,'synthetic-summary');assert.equal(query.body.target_id,'uniprot:SYNTHETIC');
 assert.equal(query.body.filters.primary_docking_protocol,'synthetic-protocol');
 assert.match(h.body.innerHTML,/vina.*17.*3개/);assert.match(h.body.innerHTML,/Synthetic &lt;target&gt;/);
 assert.match(h.get('#platform-query-status').textContent,/조회 완료/);
});
test('compact metadata does not promote held or unknown snapshots to verified',async()=>{
 for(const status of ['held','tampered','unknown',undefined]){
  const denied={...compact(),status,reason:'Synthetic integrity hold'};
  const h=await platform({snapshots:[denied],selection:null});
  assert.equal(h.calls.some(call=>call.path==='/v2/candidate-queries'),false);
  assert.match(h.body.innerHTML,/Synthetic integrity hold/);
  assert.doesNotMatch(h.body.innerHTML,/id="platform-table"/);
 }
});
test('compact selection rejects unknown snapshot, target and protocol without querying candidates',async()=>{
 for(const selection of [
  {snapshot:'not-owned',target:'uniprot:SYNTHETIC',protocol:null},
  {snapshot:'synthetic-summary',target:'uniprot:OTHER',protocol:null},
  {snapshot:'synthetic-summary',target:'uniprot:SYNTHETIC',protocol:'invented'}
 ]){
  const h=await platform({selection});
  assert.equal(h.calls.some(call=>call.path==='/v2/candidate-queries'),false);
  assert.match(h.body.innerHTML,/명시한/);
 }
});
test('generalized comparison retains stored study name, count and exact result link',async()=>{
 const h=await platform({namespace:'generalized',runs:[{run_id:'synthetic-run',snapshot_id:'synthetic-summary',
  input_descriptor:{kind:'compound',name:'Synthetic <compound>',input:{kind:'compound',format:'name',value:'Synthetic compound'}},
  search_mode:'similarity',created_at:'2026-10-20T00:00:00Z'}]});
 assert.deepEqual(h.calls.map(call=>call.path),['/v2/snapshots?view=summary','/v2/discovery-runs','/v2/candidate-queries']);
 assert.match(h.body.innerHTML,/Synthetic &lt;compound&gt;.*저장 관측 52개/);
 assert.match(h.get('#platform-study-context').innerHTML,/#discovery\?run=synthetic-run/);
});
test('catalog target choices use only verified compact snapshots and do not start calculations',async()=>{
 const h=await platform({route:'catalog',snapshots:[compact(),{...compact(),snapshot_id:'held-summary',status:'held',provenance:{marker:'Hidden held marker'}}]});
 assert.deepEqual(h.calls.map(call=>call.path),['/v2/catalog-status','/v2/snapshots?view=summary']);
 assert.match(h.body.innerHTML,/Synthetic &lt;target&gt;.*uniprot:SYNTHETIC/);
 assert.doesNotMatch(h.body.innerHTML,/held-summary|Hidden held marker/);
});
test('discovery source choices use compact verified target counts without a full snapshot request',async()=>{
 const h=harness(),calls=[],snapshots=[compact(),{...compact(),snapshot_id:'held-summary',status:'held'}];
 await renderDiscovery(h.main,async path=>{
  calls.push(path);
  if(path==='/v2/snapshots?view=summary')return{snapshots};
  if(path==='/v2/discovery-runs')return{runs:[]};
  if(path==='/v2/research-capabilities')return{};
  throw Error('Unexpected discovery request '+path);
 });
 assert.deepEqual(calls,['/v2/snapshots?view=summary','/v2/discovery-runs','/v2/research-capabilities']);
 assert.match(h.main.innerHTML,/value="synthetic-summary\|uniprot:SYNTHETIC">Synthetic &lt;target&gt; · 52개/);
 assert.doesNotMatch(h.main.innerHTML,/held-summary/);
 assert.equal(typeof h.get('#research-input-form').onsubmit,'function');
});
