let manifest,files={},generation=0,ready=false;
const primary=r=>r.primary_docking?.score??null;
const escText=r=>[r.name,r.candidate_id,r.smiles].map(x=>x??'').join(' ').toLowerCase();
const sha=async file=>[...new Uint8Array(await crypto.subtle.digest('SHA-256',await file.arrayBuffer()))].map(v=>v.toString(16).padStart(2,'0')).join('');
const field=(r,k)=>k==='primary_score'?primary(r):r[k];
onmessage=async({data})=>{let operationGeneration=generation;try{
 if(data.type==='init'){
  files=data.files;manifest=data.manifest;generation++;operationGeneration=generation;ready=false;const ticket=generation,sourceFiles=files,source=manifest;
  for(const chunk of source.chunks){if(!sourceFiles[chunk.file]||await sha(sourceFiles[chunk.file])!==chunk.sha256)throw Error('후보 묶음 누락 또는 해시 불일치: '+chunk.file)}
  if(ticket!==generation)return;ready=true;postMessage({type:'ready',count:source.candidate_count});return;
 }
 if(data.type==='query'){
  if(!ready)return;const ticket=generation,source=manifest,sourceFiles=files;
  const candidates=[],q=(data.q||'').toLowerCase();let total=0,matched=0;
  const compare=(a,b)=>{const x=field(a,data.sort),y=field(b,data.sort);if(x==null&&y!=null)return 1;if(y==null&&x!=null)return -1;let n=typeof x==='number'&&typeof y==='number'?x-y:String(x??'').localeCompare(String(y??''));return (data.descending?-n:n)||a.candidate_id.localeCompare(b.candidate_id)};
  for(const chunk of source.chunks){
   const rows=JSON.parse(await sourceFiles[chunk.file].text());if(ticket!==generation)return;if(rows.length!==chunk.count)throw Error('후보 묶음 개수 불일치');
   total+=rows.length;
   for(const r of rows){if(q&&!escText(r).includes(q))continue;if(data.purchase&&!(r.purchasable===true||r.procurement?.catalog_record))continue;if(data.critical&&(r.critical_warning===true||r.hazard_severity==='critical'))continue;if(data.admetRisk&&r.admet_search_risk?.excluded===true)continue;if(data.opposing&&r.opposing_direction===true)continue;matched++;candidates.push({...r,_chunk:chunk.file})}
   candidates.sort(compare);if(candidates.length>1000)candidates.length=1000;
  }
  if(total!==source.candidate_count)throw Error('전체 저장 후보 개수가 매니페스트와 다릅니다.');
  postMessage({type:'results',requestId:data.requestId,total,matched,displayed:Math.min(matched,1000),rows:candidates.slice(data.offset,data.offset+50)});return;
 }
}catch(e){if(operationGeneration!==generation)return;postMessage({type:'error',error:String(e.message||e)})}};
