import assert from 'node:assert/strict';
import test from 'node:test';
import fs from 'node:fs';
import {researchProgressView,applyResearchProgressState,RESEARCH_PROGRESS_STALE_MS} from '../web/public/research-progress.js';
import {renderDiscovery} from '../web/public/discovery.js';

for(const state of ['submitting','checking','queued','preparing','running','cancelling'])test(state+' indicates its own active operation',()=>{
 const view=researchProgressView(state,{checkedAt:100,now:101});
 assert.equal(view.active,true);assert.equal(view.busy,true);
 assert.equal(view.showProgress,['running','preparing'].includes(state));
});
for(const [state,tone]of [['complete','complete'],['completed','complete'],['partial','attention'],['failed','failed'],['cancelled','cancelled'],['cancellation_unverified','attention'],['cleanup_unverified','attention'],['paused','attention'],['held','attention']])test(state+' has a distinct static terminal presentation',()=>{
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
 const fixture=harness(),calls=[];let fail=false,state='running';
 const record=()=>({run_id:'fixture',status:state,input:{value:'Synthetic public input'},counts:{},candidates:[],progress:{stage:'catalog_acquisition',processed:2,total:5},...(state==='complete'?{snapshot_id:'synthetic-snapshot'}:{})});
 const api=async(path,options)=>{
  calls.push({path,method:options?.method||'GET'});
  if(path==='/v2/snapshots')return{snapshots:[]};
  if(path==='/v2/research-capabilities')return{};
  if(path==='/v2/discovery-runs')return{runs:[record()]};
  if(path.endsWith('/summary'))return record();
  if(path==='/v2/discovery-runs/fixture'){if(fail)throw Error('synthetic connection loss');return record()}
  throw Error('Unexpected fixture request '+path);
 };
 try{
  await renderDiscovery(fixture.main,api);fixture.get('#view').dataset.viewRun='0';await fixture.get('#view').onclick();
  const panel=fixture.get('#research-live-status');assert.equal(panel.attributes['data-state'],'running');assert.equal(panel.attributes['data-active'],'true');assert.match(panel.innerHTML,/카탈로그 후보 확보.*2개 \/ 5개/);
  fail=true;await next(5000);assert.equal(panel.attributes['data-state'],'unknown');assert.equal(panel.attributes['aria-busy'],'false');assert.match(panel.textContent,/마지막으로 확인한 상태/);assert.doesNotMatch(panel.innerHTML,/현재 단계|research-progress-icon/);
  fail=false;await next(10000);assert.equal(panel.attributes['data-active'],'true');await next(RESEARCH_PROGRESS_STALE_MS);assert.equal(panel.attributes['data-active'],'false');assert.match(panel.textContent,/최근 20초/);
  await next(5000);assert.equal(panel.attributes['data-active'],'true');
  state='complete';await next(5000);assert.equal(panel.attributes['data-state'],'complete');assert.equal(panel.attributes['aria-busy'],'false');assert.match(panel.innerHTML,/✓.*연구 완료/);assert.doesNotMatch(panel.innerHTML,/현재 단계|2개 \/ 5개/);assert.match(fixture.get('#research-result-actions').innerHTML,/Excel/);
  assert.equal(timers.filter(token=>!token.cancelled).length,0);assert.ok(calls.every(call=>call.method==='GET'),'Progress never submits or restarts a study');
 }finally{globalThis.setTimeout=originalTimer;globalThis.clearTimeout=originalClear}
});
