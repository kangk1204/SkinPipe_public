import assert from 'node:assert/strict';
import test from 'node:test';
import {requestJson,requestTimeout} from '../web/public/api-request.js';
import {createLocalResearchApi,localRequestTimeout} from '../web/public/local-research-app.js';

// Advance a deterministic clock without network requests or wall-clock sleeps.
// Native AbortControllers exercise the production wrapper's AbortSignal.any path.
function clock(){
 let now=0,serial=0;
 const timers=new Map(),budgets=[];
 const later=(delay,fn)=>{const id=++serial;timers.set(id,{at:now+delay,fn});return id};
 const signalFactory=ms=>{
  budgets.push(ms);const controller=new AbortController();
  later(ms,()=>controller.abort(new DOMException('deadline','TimeoutError')));
  return controller.signal;
 };
 const advance=async ms=>{
  const target=now+ms;
  while(true){
   const next=[...timers].filter(([,timer])=>timer.at<=target).sort((a,b)=>a[1].at-b[1].at||a[0]-b[0])[0];
   if(!next)break;
   timers.delete(next[0]);now=next[1].at;next[1].fn();
   await Promise.resolve();
  }
  now=target;
  for(let i=0;i<6;i++)await Promise.resolve();
 };
 return{later,signalFactory,advance,budgets};
}
function server(clock,{delay=77000,loss=null}={}){
 const calls=[],accepted=[];
 const fetcher=(url,request)=>new Promise((resolve,reject)=>{
  calls.push({url,request});
  const abort=()=>reject(request.signal.reason);
  if(request.signal?.aborted){abort();return}
  request.signal?.addEventListener('abort',abort,{once:true});
  clock.later(delay,()=>{
   // The server can accept after a client abort; no second POST is safe.
   accepted.push({run_id:'synthetic-run'});
   request.signal?.removeEventListener('abort',abort);
   if(loss==='headers'){reject(new TypeError('accepted; response lost'));return}
   resolve({ok:true,json:async()=>{
    if(loss==='body')throw new TypeError('accepted; response body lost');
    return{run_id:'synthetic-run',status:'queued'};
   }});
  });
 });
 return{fetcher,calls,accepted};
}
const path='/v2/discovery-runs',post={method:'POST',body:'{"synthetic":true}'};

test('only the exact discovery collection path receives the bounded 600-second deadline',()=>{
 assert.equal(requestTimeout(path),600000);
 assert.equal(localRequestTimeout(path),600000);
 for(const unchanged of ['/v2/health','/v2/setup/status','/v2/discovery-runs/','/v2/discovery-runs-other','/v2/discovery-runs/synthetic-run','/v2/discovery-runs/synthetic-run/summary','/v2/discovery-runs/synthetic-run/cancel']){
  assert.equal(requestTimeout(unchanged),20000,unchanged);
  assert.equal(localRequestTimeout(unchanged),20000,unchanged);
 }
 assert.equal(requestTimeout('/v2/research-inputs/resolve'),600000);
 assert.equal(localRequestTimeout('/v2/research-inputs/resolve'),120000);
 assert.equal(requestTimeout('/v2/exports/synthetic/download'),1800000);
});
for(const delay of [77000,240000])test('production requestJson accepts a simulated '+delay+'ms receptor-preparation response once',async()=>{
 const time=clock(),remote=server(time,{delay}),signal=time.signalFactory(requestTimeout(path));
 let settled=false;
 const pending=requestJson(remote.fetcher,path,path,{...post,signal}).then(value=>{settled=true;return value});
 await time.advance(20000);assert.equal(settled,false);assert.equal(signal.aborted,false);
 await time.advance(delay-20000);
 assert.deepEqual(await pending,{run_id:'synthetic-run',status:'queued'});
 assert.equal(remote.calls.length,1);assert.equal(remote.accepted.length,1);
});
for(const delay of [77000,240000])test('createLocalResearchApi inherits the real deadline for a '+delay+'ms response',async()=>{
 const time=clock(),remote=server(time,{delay});
 const api=createLocalResearchApi({fetcher:remote.fetcher,signalFactory:time.signalFactory,document:null});
 const pending=api(path,post);
 assert.deepEqual(time.budgets,[600000]);
 await time.advance(20000);assert.equal(remote.calls[0].request.signal.aborted,false);
 await time.advance(delay-20000);assert.equal((await pending).run_id,'synthetic-run');
 assert.equal(remote.calls.length,1);assert.equal(remote.calls[0].request.credentials,'same-origin');
 assert.equal(remote.calls[0].request.mode,'same-origin');
});
test('the original 20-second deadline reproduces a lost acceptance before the server finishes',async()=>{
 const time=clock(),remote=server(time);
 const pending=requestJson(remote.fetcher,path,path,{...post,signal:time.signalFactory(20000)});
 const rejected=assert.rejects(pending,{name:'TimeoutError'});
 await time.advance(20000);await rejected;
 await time.advance(57000);assert.equal(remote.accepted.length,1);assert.equal(remote.calls.length,1);
});
test('the new deadline still aborts at 600 seconds and never resubmits a late acceptance',async()=>{
 const time=clock(),remote=server(time,{delay:600001});
 const api=createLocalResearchApi({fetcher:remote.fetcher,signalFactory:time.signalFactory,document:null});
 let settled=false;
 const pending=api(path,post).finally(()=>{settled=true});
 const rejected=assert.rejects(pending,{name:'TimeoutError'});
 await time.advance(599999);assert.equal(settled,false);
 await time.advance(1);await rejected;
 assert.equal(remote.calls.length,1);assert.equal(remote.accepted.length,0);
 await time.advance(1);assert.equal(remote.accepted.length,1);assert.equal(remote.calls.length,1);
});
for(const loss of ['headers','body'])test('lost acceptance '+loss+' never causes an automatic duplicate POST',async()=>{
 const time=clock(),remote=server(time,{loss});
 const api=createLocalResearchApi({fetcher:remote.fetcher,signalFactory:time.signalFactory,document:null});
 const rejected=assert.rejects(api(path,post),TypeError);
 await time.advance(77000);await rejected;
 await time.advance(600000);assert.equal(remote.accepted.length,1);assert.equal(remote.calls.length,1);
});
test('caller cancellation stays effective and cannot automatically duplicate an accepted run',async()=>{
 const time=clock(),remote=server(time),caller=new AbortController();
 const api=createLocalResearchApi({fetcher:remote.fetcher,signalFactory:time.signalFactory,document:null});
 const rejected=assert.rejects(api(path,{...post,signal:caller.signal}),{name:'AbortError'});
 await time.advance(10000);caller.abort();await rejected;
 assert.equal(remote.calls[0].request.signal.aborted,true);
 await time.advance(67000);assert.equal(remote.accepted.length,1);assert.equal(remote.calls.length,1);
});
