import assert from 'node:assert/strict';
import test from 'node:test';
import {showPlatformDocking} from '../web/public/structure-viewer.js';

async function renderedSeed(fields){
 const previousDocument=globalThis.document,nodes=new Map(),calls=[];
 const get=selector=>{
  if(!nodes.has(selector))nodes.set(selector,{
   innerHTML:'',textContent:'',hidden:false,disabled:false,open:false,dataset:{},
   classList:{add(){},remove(){}},setAttribute(){},removeAttribute(){},
   showModal(){this.open=true},close(){this.open=false;this.onclose?.()}
  });
  return nodes.get(selector);
 };
 globalThis.document={querySelector:get,querySelectorAll:()=>[],activeElement:null};
 let release;
 const payload=new Promise(resolve=>{release=resolve});
 const row={engine:'vina',stage:'broad_docking',score:-5.25,artifact_valid:true,pose_sha256:'a'.repeat(64),...fields};
 const candidate={candidate_id:'synthetic-candidate',target_id:'uniprot:SYNTHETIC',docking:[row]};
 const original=structuredClone(candidate);
 try{
  showPlatformDocking({
   title:'Synthetic pose',candidate,context:{snapshot:'synthetic-snapshot',target:'uniprot:SYNTHETIC'},
   api:(path,options)=>{calls.push({path,options});return payload}
  });
  const dialog=get('#structure-dialog'),html=dialog.innerHTML;
  assert.equal(dialog.open,true);
  assert.equal(calls.length,1);
  assert.equal(calls[0].path,'/v2/candidates/synthetic-candidate/docking?snapshot_id=synthetic-snapshot&target_id=uniprot%3ASYNTHETIC&observation_index=0');
  const methods=JSON.parse(get('#structure-method-record').textContent);
  assert.deepEqual(methods.selected_entry,JSON.parse(JSON.stringify({kind:'accepted',pose_index:1,protocol:fields.protocol})));
  assert.deepEqual(candidate,original,'Opening a pose must not modify the stored seed or protocol');
  return html.match(/<option[^>]*>([^<]+)<\/option>/)?.[1];
 }finally{
  get('#structure-dialog').close();
  release({});
  // Finish the actual display task after closing, before restoring the document.
  await Promise.resolve();await Promise.resolve();await Promise.resolve();
  if(previousDocument===undefined)delete globalThis.document;else globalThis.document=previousDocument;
 }
}
test('platform pose label prefers legacy source_seed, including zero',async()=>{
 for(const seed of [11,0])assert.match(await renderedSeed({source_seed:seed,seed:22,protocol:{seed:42}}),new RegExp('시드 '+seed+' ·'));
});
test('platform pose label retains top-level seed before protocol seed, including zero',async()=>{
 for(const seed of [22,0])assert.match(await renderedSeed({source_seed:null,seed,protocol:{seed:42}}),new RegExp('시드 '+seed+' ·'));
});
test('platform pose label uses the recorded routine protocol seed, including zero',async()=>{
 for(const seed of [42,0])assert.match(await renderedSeed({source_seed:null,seed:null,protocol:{seed}}),new RegExp('시드 '+seed+' ·'));
});
test('platform pose label keeps missing seeds unrecorded without inferring a default',async()=>{
 for(const fields of [{},{protocol:{}},{source_seed:null,seed:null,protocol:{seed:null}}]){
  assert.match(await renderedSeed(fields),/시드 미기록 ·/);
 }
});
