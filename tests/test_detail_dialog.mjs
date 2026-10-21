// Real shipped HTML and browser handlers; all API records are synthetic.
import assert from 'node:assert/strict';
import {createServer} from 'node:http';
import {readFile} from 'node:fs/promises';
import {existsSync} from 'node:fs';
import {createRequire} from 'node:module';
import {fileURLToPath} from 'node:url';
import path from 'node:path';
const {chromium}=createRequire(import.meta.url)('playwright');
const root=['../web/public/','../public/'].map(name=>fileURLToPath(new URL(name,import.meta.url))).find(name=>existsSync(path.join(name,'local-research.html')));
assert.ok(root,'shipped public assets must exist');
const server=createServer(async(req,res)=>{
 try{
  const name=decodeURIComponent(new URL(req.url,'http://localhost').pathname);
  const file=path.resolve(root,'.'+(name==='/'?'/local-research.html':name));
  if(!file.startsWith(root))throw Error('outside fixture root');
  const body=await readFile(file);
  res.setHeader('Content-Type',({'.html':'text/html','.js':'text/javascript','.css':'text/css'})[path.extname(file)]||'application/octet-stream');
  res.end(body);
 }catch{res.writeHead(404);res.end()}
});
await new Promise(resolve=>server.listen(0,'127.0.0.1',resolve));
const base='http://127.0.0.1:'+server.address().port;
let browser;
try{
 browser=await chromium.launch({headless:true,...(process.env.SKINPIPE_CHROMIUM_EXECUTABLE?{executablePath:process.env.SKINPIPE_CHROMIUM_EXECUTABLE}:{})});
 const page=await browser.newPage({viewport:{width:1440,height:1080}}),errors=[];
 page.on('pageerror',error=>errors.push(error.message));
 const rows=[
  {candidate_id:'synthetic-seed',name:'Synthetic input',role:'query_seed'},
  {candidate_id:'synthetic-candidate',name:'Synthetic candidate',role:'candidate'},
  {candidate_id:'synthetic-reference',name:'Synthetic reference',role:'control',control_kind:'measured_reference'},
  {candidate_id:'synthetic-redock',name:'Synthetic redock',role:'control',redock_reference:true}
 ].map(row=>({...row,target_id:'synthetic-target',smiles:'CCO',priority_score_v2:60.269,
  priority_score_observed_channels:84.23,priority_imputed_weight:0.7,priority_coverage:0.3,
  topical_score:9.18,discovery_terms:{topical:50},skin_goal_assessment:{permeability_ranking:'neutral'},
  offers:[{product_url:'https://example.org/purchase/synthetic',supplier:'Synthetic supplier'}],
  admet:{model_id:'admet-ai-synthetic',status:'success',values:{AMES:0.2}}}));
 const counts={population:4,filtered:4,displayed:4,population_calculated:0,
  reasons:{critical:0,opposing:0},missing:{purchase:4,pchembl:4}};
 await page.context().route('**/*',async route=>{
  const url=new URL(route.request().url());
  if(url.href==='https://example.org/purchase/synthetic')return route.fulfill({contentType:'text/html',body:'Synthetic supplier destination'});
  if(url.origin!==base)return route.abort();
  if(!url.pathname.startsWith('/v2/'))return route.continue();
  let value;
  if(url.pathname==='/v2/snapshots')value={snapshots:[{snapshot_id:'synthetic-snapshot',status:'verified',
   candidate_count:4,target_counts:{'synthetic-target':4},provenance:{target_id:'synthetic-target'}}]};
  else if(url.pathname==='/v2/discovery-runs')value={runs:[]};
  else if(url.pathname==='/v2/candidate-queries')value={snapshot_id:'synthetic-snapshot',target_id:'synthetic-target',rows,counts,next_offset:null};
  else if(url.pathname.startsWith('/v2/candidates/'))value=rows.find(row=>row.candidate_id===url.pathname.split('/').at(-1));
  else if(url.pathname==='/v2/catalog-status')value={compound_count:1};
  else if(url.pathname==='/v2/catalog-queries')value={status:'complete',elapsed_seconds:0,rows:[{compound_id:'synthetic-catalog',identity_smiles:'CCO'}]};
  else if(url.pathname==='/v2/catalog-candidates/synthetic-catalog')value={compound:{compound_id:'synthetic-catalog'},offers:[]};
  else throw Error('Unexpected API request '+url.pathname);
  return route.fulfill({status:200,contentType:'application/json',body:JSON.stringify(value)});
 });
 await page.goto(base+'/#research');
 await page.locator('[data-detail="0"]').waitFor();
 for(const [index,label] of ['출발구조','탐색 후보','측정 자료 참조','관측 리간드 재도킹 대조'].entries()){
  assert.ok((await page.locator('#platform-table tbody tr').nth(index).innerText()).includes(label));
 }
 assert.ok((await page.locator('#platform-table').innerText()).includes('비중 70%'));
 assert.ok((await page.locator('#platform-table').innerText()).includes('연구 목적에 따라 중립 50'));
 await page.locator('[data-detail="0"]').click();
 await page.locator('#detail-dialog[open]').waitFor();
 assert.equal(await page.locator('#detail-title').innerText(),'Synthetic input');
 assert.ok((await page.locator('#detail-content').innerText()).includes('84.23'));
 const purchase=page.getByRole('link',{name:'구매처 홈페이지 열기 ↗'});
 assert.equal(await purchase.getAttribute('href'),'https://example.org/purchase/synthetic');
 assert.equal(await purchase.getAttribute('target'),'_blank');
 const destination=page.waitForEvent('popup');
 await purchase.click();
 const popup=await destination;
 await popup.waitForLoadState();
 assert.equal(popup.url(),'https://example.org/purchase/synthetic');
 await popup.close();
 if(process.env.SKINPIPE_BROWSER_EVIDENCE)await page.screenshot({path:path.join(process.env.SKINPIPE_BROWSER_EVIDENCE,'candidate-detail.png')});
 await page.keyboard.press('Escape');
 assert.equal(await page.locator('#detail-dialog').evaluate(node=>node.open),false);
 await page.locator('[data-detail="1"]').click();
 await page.locator('#detail-dialog[open]').waitFor();
 assert.equal(await page.locator('#detail-title').innerText(),'Synthetic candidate');
 await page.locator('#platform-close').click();
 await page.evaluate(async()=>{
  const {renderPlatform}=await import('/platform.js');
  const {createLocalResearchApi}=await import('/local-research-app.js');
  await renderPlatform({main:document.getElementById('local-research-main'),api:createLocalResearchApi(),route:'catalog'});
 });
 await page.locator('#platform-catalog-search button[type="submit"]').click();
 await page.locator('[data-catalog-detail="0"]').click();
 await page.locator('#detail-dialog[open]').waitFor();
 assert.equal(await page.locator('#detail-title').innerText(),'구조·공급·실험 근거');
 if(process.env.SKINPIPE_BROWSER_EVIDENCE)await page.screenshot({path:path.join(process.env.SKINPIPE_BROWSER_EVIDENCE,'catalog-detail.png')});
 await page.locator('#catalog-detail-close').click();
 assert.equal(await page.locator('#detail-dialog').count(),1);
 assert.deepEqual(errors,[]);
 console.log('PASS: shipped HTML candidate/catalog dialogs, roles, score scope, close/reopen; no browser errors');
}finally{
 await browser?.close();
 await new Promise(resolve=>server.close(resolve));
}
