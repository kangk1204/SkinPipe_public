import {literatureMarkup} from './literature-evidence.js';
import {pendingSubmission,acknowledgeSubmission} from './submission-state.js';

const BASE='/v2/literature',PREFERRED='gpt-6-luna',EFFORTS=['none','low','medium','high','xhigh','max'];
const MAX_SOURCE=4*1024*1024,JOB=/^[a-f0-9]{32}$/,HASH=/^[a-f0-9]{64}$/;
const record=value=>value!==null&&typeof value==='object'&&!Array.isArray(value);
const esc=value=>String(value??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const text=(value,fallback='미기록')=>typeof value==='string'&&value.trim()?value.slice(0,3000):fallback;
const rows=value=>Array.isArray(value)?value.filter(record):[];
const pending=job=>job?.status==='queued'||job?.status==='running';
export function literatureModelChoice(value,previous=''){
 const models=new Map();for(const row of record(value?.error)||value?.status==='HOLD'?[]:rows(value?.models).slice(0,1000)){const id=row.id??row.slug;if(typeof id==='string'&&/^[A-Za-z0-9][A-Za-z0-9._:-]{0,127}$/.test(id)&&!models.has(id))models.set(id,row)}
 const declared=typeof value?.default_model==='string'&&models.has(value.default_model)?value.default_model:'';
 const id=models.has(previous)?previous:models.has(PREFERRED)?PREFERRED:declared;
 const row=models.get(id),supported=Array.isArray(row?.supported_reasoning_efforts)?EFFORTS.filter(e=>row.supported_reasoning_efforts.includes(e)):EFFORTS;
 const effort=id===PREFERRED&&supported.includes('max')?'max':supported.includes(row?.default_reasoning_effort)?row.default_reasoning_effort:'';
 return{id,effort,supported,models,basis:id===PREFERRED?'preferred':id&&id===declared?'provider_default':id?'manual':'manual_required'};
}
const jobID=job=>record(job)&&JOB.test(job.review_id)&&(!job.job_id||job.job_id===job.review_id)?job.review_id:null;
const pathOK=value=>typeof value==='string'&&value.length<=256&&!/[\\\u0000-\u001f\u007f]/.test(value)&&value.split('/').every(part=>part!==''&&part!=='.'&&part!=='..');
const labels={queued:'대기 중',running:'문헌 수집·검토 중',complete:'작업 종료 · 아래 근거 판정 확인',failed:'작업 실패 · 자동 재시도 없음',cancelled:'작업 취소됨',budget_deferred:'시간 한도로 보류 · 자동 재개 없음'};
const statusLabel=value=>typeof value==='string'&&Object.hasOwn(labels,value)?labels[value]:'상태 미확인 · 원본 확인 필요';
function errorText(error){
 const code=error?.code;
 if(['usage_exhausted','usage_limit_exceeded'].includes(code))return '이 앱에서 사용할 수 있는 ChatGPT 사용량 한도에 도달했습니다. 수집 문헌은 보존되며 AI 검토는 보류합니다. 상단 사용량을 확인하세요.';
 if(['rate_limited','rate_limit_exceeded'].includes(code)||error?.status===429)return 'ChatGPT 요청이 제한됐습니다. 자동 재호출 없이 보류합니다. 상단 사용량을 확인한 뒤 직접 다시 실행하세요.';
 if(error?.status===401||error?.status===403)return '이 문헌 검토 기록에 접근할 권한이 없습니다. 분석 서버 연결과 개인 계정 설정을 확인하세요.';
 if(error?.status===400)return '입력한 이름·종·문헌 수를 확인해 주세요. 화합물은 분자군이나 혼합물 대신 정확한 개별 물질 이름이 필요하며, 단백질은 유전자명으로 입력해 주세요. 구조·실험 결과는 보내지 않습니다.';
 if(error?.status===404)return '이 서버에서 개인 문헌 검토를 제공하지 않거나 요청한 기록을 찾을 수 없습니다.';
 return '요청 결과를 확인하지 못했습니다. 저장 목록을 새로고침해 접수 상태를 먼저 확인하세요. 자동 재접수는 하지 않습니다.';
}
function modelErrorText(error){return ['usage_exhausted','usage_limit_exceeded','rate_limited','rate_limit_exceeded'].includes(error?.code)||[401,403,429].includes(error?.status)?errorText(error):'모델 목록을 가져오지 못했습니다. 모델 목록 새로고침으로 다시 조회하세요. 이 조회는 AI 요약을 실행하지 않습니다.'}
function dateText(value){
 if(typeof value==='number'&&Number.isFinite(value)&&value>=0&&value<253402300799)return new Date(value*1000).toISOString();
 return typeof value==='string'&&value.length<=80?value:'시점 미기록';
}
function seconds(value){return typeof value==='number'&&Number.isFinite(value)&&value>=0?value.toFixed(1)+'초':'미기록'}
function safeStorage(storage){
 const memory=new Map();
 return {getItem(key){try{return storage?.getItem(key)??memory.get(key)??null}catch{return memory.get(key)??null}},
  setItem(key,value){memory.set(key,value);try{storage?.setItem(key,value)}catch{}},
  removeItem(key){memory.delete(key);try{storage?.removeItem(key)}catch{}}};
}
export function personalLiteratureLink(query){return typeof query==='string'&&query.trim()&&query.length<=160?'#literature?'+new URLSearchParams({query:query.trim(),source:'discovery'}):'#literature'}
export function personalLiteratureMarkup({prefill=null,publicResearch=false,legacyHistory=true}={}){
 const query=!publicResearch&&typeof prefill?.query==='string'?prefill.query.slice(0,160):'';
 return `<div class="page-heading"><div><h1>${publicResearch?'공개 논문으로 피부 관련성 살펴보기':'논문 근거 검토'}</h1><p>공개 문헌을 수집하고 본인 ChatGPT 플랜으로 검토합니다. gpt-6-luna · max를 우선 선택하며 실행 전 모델을 확인합니다.</p><p>${publicResearch?(legacyHistory?'<a class="inline-link" href="#literature?source=legacy">기존 단백질 문헌 기록 보기 →</a>':''):'<a class="inline-link" href="#literature?source=generalized">화합물·단백질의 새 공개 문헌 검토 →</a>'}</p></div></div><section class="platform-panel"><h2>내 ChatGPT 플랜 연결</h2><p id="literature-account" role="status">개인 연결 설정을 확인하고 있습니다.</p><p>${publicResearch?'이 화면은 새 연구 서버의 연결과 실제 모델 목록을 확인합니다. 기존 단백질 문헌 기록과 연결 설정이 다를 수 있습니다. 플랜 잔여 한도는 조회하지 않으며 실제 응답 가능 여부는 실행 결과에서 확인합니다.':'상단의 ChatGPT 사용량은 이 실행 환경의 기록입니다. 플랜 잔여 한도는 조회되지 않으며 본인 사용량 설정에서 확인하세요. Site 로그인과 ChatGPT 플랜 사용 동의는 별개입니다.'}</p><button id="literature-refresh" class="button" type="button">연결·저장 목록 새로고침</button></section>${prefill?.fromDiscovery&&!publicResearch?'<div class="notice" role="note">기존 분석의 표적을 가져온 새 문헌 검토 초안입니다. 이번 실행은 문헌 근거 기록만 새로 만들며 기존 입력·후보·순위·계산 결과를 변경하지 않습니다.</div>':''}<form id="literature-form" class="platform-panel"><h2>새 문헌 근거 기록</h2><div class="platform-toolbar"><label><span id="literature-query-label">${publicResearch?'공개 화합물 이름':'표적 유전자·식별자'}</span><input id="literature-query" required maxlength="${publicResearch?512:160}" placeholder="${publicResearch?'예: Nicotinamide':'예: GSK3B'}" value="${esc(query)}" autocomplete="off"></label><label>${publicResearch?'무엇의 문헌을 살펴볼까요?':'표적 종류'}<select id="literature-kind">${publicResearch?'<option value="compound">화합물 이름</option><option value="protein">단백질 유전자명</option>':'<option value="protein">단백질</option><option value="complex">복합체 · 단일 분자 표적 선택 필요</option><option value="endpoint">효능 지표 · 분자 표적 선택 필요</option>'}</select></label>${publicResearch?'<label id="literature-topic-label">관심 주제 (영문·선택)<input id="literature-topic" maxlength="128" placeholder="예: skin barrier" autocomplete="off"></label>':''}<label>종 (NCBI taxon ID)<input id="literature-taxon" type="number" min="1" max="100000000" value="9606" required></label><label id="literature-action-label">원하는 작용<select id="literature-action"><option value="any">지정하지 않음</option><option value="inhibit">억제</option><option value="activate">활성화</option></select></label><label>최대 문헌 수<input id="literature-papers" type="number" min="1" max="8" value="2" required></label></div><div class="platform-toolbar"><label>검토 모델<select id="literature-model"><option value="">실제 모델 목록 확인 후 선택</option></select></label><label>추론 강도<select id="literature-effort">${EFFORTS.map(value=>`<option value="${value}"${value==='max'?' selected':''}>${value}</option>`).join('')}</select></label><button id="literature-models" class="button" type="button">모델 목록 새로고침</button></div><p id="literature-model-status" role="status">Luna·max 우선 → 제공자가 명시한 계정 기본 모델 → 목록에서 직접 선택. 선택 후 실행 버튼을 눌러야 호출합니다.</p>${publicResearch?'<p id="literature-public-scope" class="platform-note">직접 입력한 공개 이름·관심 주제와 공개 논문만 전달합니다. 기존 연구의 이름·SMILES·구조 파일·좌표·후보 결과를 자동으로 가져오지 않습니다. 화합물 이름은 공식 PubChem에서 하나의 분자로 확인한 뒤 검색합니다.</p>':''}<label><input id="literature-consent" type="checkbox" required> ${publicResearch?'직접 입력한 공개 이름·관심 주제와 수집한 공개 문헌을 내 ChatGPT 플랜으로 보내 검토합니다.':'입력한 표적 식별자와 수집한 공개 문헌을 내 ChatGPT 플랜으로 보내 검토합니다.'}</label><p>후보 구조·실험 데이터는 이 요청에 포함하지 않습니다. 문헌 검토는 도킹·후보 재스크리닝을 실행하지 않습니다. 제한된 검색이며 사람의 논문·실험 조건 검토와 실제 실험이 필요합니다.</p><button id="literature-submit" class="button primary" type="submit" disabled>공개 문헌 수집·AI 검토 실행</button><p id="literature-message" role="status" aria-live="polite"></p></form><section class="platform-panel"><h2>저장한 문헌 검토 작업</h2><p>접수된 작업은 최대 900초의 실행 예산을 사용합니다. AI 응답 한 번의 한도는 최대 300초입니다. 대기·실패·보류 기록을 자동 재개하지 않습니다.</p><div id="literature-list"></div></section><section id="literature-review" aria-live="polite"></section>`;
}
export function boundLiteratureSources(dossier){
 if(!record(dossier)||dossier.schema!=='skinpipe-literature-evidence-v1')return [];
 const closure=new Map(),conflicts=new Set();
 for(const item of [...rows(dossier.literature_artifacts),...rows(dossier.artifact_paths)]){
  const sha=item.literature_artifact_sha256??item.sha256;
  if(!pathOK(item.path)||!HASH.test(sha)||!Number.isSafeInteger(item.size_bytes)||item.size_bytes<0||item.size_bytes>MAX_SOURCE)continue;
  const previous=closure.get(item.path);if(previous&&(previous.sha256!==sha||previous.size_bytes!==item.size_bytes))conflicts.add(item.path);
  closure.set(item.path,{path:item.path,sha256:sha,size_bytes:item.size_bytes});
 }
 return rows(dossier.sources).flatMap(source=>{const receipt=source.artifacts?.text,entry=record(receipt)?closure.get(receipt.path):null;
  return entry&&!conflicts.has(entry.path)&&entry.sha256===receipt.sha256&&entry.size_bytes===receipt.size_bytes?[{...entry,title:text(source.title,text(source.source_id,'저장 문헌'))}]:[];});
}
export async function verifyLiteratureSource(dto,expected,cryptoImpl=globalThis.crypto){
 if(!record(dto)||!record(expected)||!pathOK(dto.path)||dto.path!==expected.path||!HASH.test(dto.sha256)||dto.sha256!==expected.sha256||dto.size_bytes!==expected.size_bytes||typeof dto.text!=='string')throw Error('원문 경로·해시·크기 기록이 일치하지 않습니다.');
 const bytes=new TextEncoder().encode(dto.text);
 if(bytes.byteLength>MAX_SOURCE||bytes.byteLength!==dto.size_bytes||!cryptoImpl?.subtle)throw Error('원문 크기 또는 해시 검증 환경을 확인할 수 없습니다.');
 const actual=Array.from(new Uint8Array(await cryptoImpl.subtle.digest('SHA-256',bytes)),n=>n.toString(16).padStart(2,'0')).join('');
 if(actual!==expected.sha256)throw Error('원문 SHA256이 저장 기록과 다릅니다.');
 return {text:dto.text,path:dto.path,sha256:actual,size_bytes:bytes.byteLength};
}
export function prepareLiteratureJSON(dossier,id,{anchor,document=globalThis.document,urlAPI=globalThis.URL}={}){
 if(!JOB.test(id)||!record(dossier)||dossier.schema!=='skinpipe-literature-evidence-v1'||!anchor?.insertAdjacentElement)throw Error('문헌 파일의 기록과 저장 위치를 확인하세요.');
 const raw=JSON.stringify(dossier,null,2);if(new TextEncoder().encode(raw).length>MAX_SOURCE)throw Error('JSON 표시·다운로드 한도를 초과했습니다.');
 const link=document.createElement('a'),url=urlAPI.createObjectURL(new Blob([raw],{type:'application/json'}));
 try{link.href=url;link.download='SkinPipe-literature-'+id+'.json';link.id='literature-ready-json';link.className='button';link.textContent='준비된 문헌 JSON 저장';anchor.insertAdjacentElement('afterend',link);link.focus();link.scrollIntoView?.({block:'nearest',behavior:'instant'});}catch(error){urlAPI.revokeObjectURL(url);link.remove();throw error;}
 let disposed=false;return{dispose(){if(disposed)return;disposed=true;urlAPI.revokeObjectURL(url);link.remove()}};
}
export function createPersonalLiteratureController({main,api,isCurrent=()=>true,offline=false,prefill=null,publicResearch=false,legacyHistory=true,onUsageChange=async()=>{},storage=null,uuid=()=>crypto.randomUUID(),schedule=(fn,delay)=>setTimeout(fn,delay),cancel=timer=>clearTimeout(timer),cryptoImpl=globalThis.crypto,saveDossier=prepareLiteratureJSON,maxPolls=180}={}){
 if(!main||typeof api!=='function')throw TypeError('Personal literature container and API are required');
 main.innerHTML=personalLiteratureMarkup({prefill,publicResearch,legacyHistory});
 const $=id=>main.querySelector('#'+id),store=safeStorage(storage),pendingKey=publicResearch?'skinpipe-public-research-literature-pending':'skinpipe-personal-literature-pending';
 let stopped=false,timer=null,selected=null,catalog=new Set(),modelDetails=null,connection=null,submitting=false,polls=0,detailVersion=0,sourceVersion=0,listVersion=0,refreshing=false,modelsLoading=false;
 let preparedDownload=null;const clearDownload=()=>{preparedDownload?.dispose?.();preparedDownload=null};
 const current=()=>!stopped&&isCurrent();
 const configured=()=>connection?.configured===true&&['personal_self_hosted','local_single_user'].includes(connection.mode);
 const allowed=()=>configured()&&connection?.connected===true&&connection?.plan_enabled===true;
 function message(value){if(current())$('literature-message').textContent=value}
 function controls(){if(!current())return;$('literature-submit').disabled=!allowed()||submitting||refreshing||modelsLoading||pending(selected)||!catalog.has($('literature-model').value)||!EFFORTS.includes($('literature-effort').value);$('literature-models').disabled=!allowed()||submitting||refreshing||modelsLoading;$('literature-refresh').disabled=refreshing}
 function modelControls(keep=false){const choice=literatureModelChoice(modelDetails,$('literature-model').value),previous=$('literature-effort').value;$('literature-effort').innerHTML='<option value="">추론 강도를 선택하세요</option>'+choice.supported.map(e=>`<option value="${e}">${e}</option>`).join('');$('literature-effort').value=keep&&choice.supported.includes(previous)?previous:choice.effort;controls()}
 function modelMessage(){const id=$('literature-model').value,effort=$('literature-effort').value;$('literature-model-status').textContent=id?'선택 모델 '+id+' · 추론 '+(effort||'선택 필요')+'. 실행 버튼을 눌러야 본인 플랜으로 요청합니다.':'목록에서 모델과 지원 추론 강도를 선택하세요. 아직 AI 요청을 보내지 않았습니다.';controls()}
 async function usage(){try{await onUsageChange()}catch{}}
 async function loadModels(){
  if(!current()||!allowed()||modelsLoading)return;
  modelsLoading=true;const previous=$('literature-model').value;catalog=new Set();$('literature-model-status').textContent='연결한 계정의 모델 목록을 불러오고 있습니다… AI 요약은 실행하지 않습니다.';controls();
  try{const value=await api(BASE+'/models');if(!current())return;
   $('literature-model').innerHTML='<option value="">모델을 선택하세요</option>';
   modelDetails=value;const choice=literatureModelChoice(value,previous);catalog=new Set(choice.models.keys());
   $('literature-model').innerHTML+=[...catalog].map(id=>`<option value="${esc(id)}">${esc(id)}</option>`).join('');
   $('literature-model').value=choice.id;modelControls(!!previous&&previous===choice.id);
   $('literature-model-status').textContent=record(value?.error)?modelErrorText(value.error):choice.basis==='preferred'?'gpt-6-luna 우선 선택 · 지원되는 추론 강도를 표시합니다. 실행 전 모델·추론·플랜 사용을 확인하세요.':choice.basis==='provider_default'?'Luna가 없어 제공자가 명시한 계정 기본 모델 '+choice.id+'을 선택했습니다. 실행 전 확인하세요.':catalog.size?'gpt-6-luna와 제공자 기본 모델을 확인하지 못했습니다. 목록에서 직접 선택하세요. 임의 모델은 자동 선택하지 않습니다.':'모델 목록을 확인하지 못했습니다. 자동 대체·AI 요청은 실행하지 않습니다.';
  }catch(error){if(current()){$('literature-model').innerHTML='<option value="">모델 목록 확인 필요</option>';$('literature-model').value='';$('literature-model-status').textContent=modelErrorText(error)}}
  finally{modelsLoading=false;controls()}
 }
 function cancelPoll(){cancel(timer);timer=null}
 function nextPoll(id,version){
  cancelPoll();if(!current()||!pending(selected))return;
  if(++polls>maxPolls){message('자동 상태 확인 대기가 끝났습니다. 저장 목록을 새로고침해 실제 상태를 확인하세요. 작업을 재접수하거나 재개하지 않았습니다.');return}
  timer=schedule(()=>readReview(id,version),5000);
 }
 function showReview(job){
  if(!current())return;const id=jobID(job);if(!id)throw Error('Invalid review identity');selected=job;clearDownload();
  const dossier=record(job.result)?job.result.dossier:null,bound=boundLiteratureSources(dossier);
  $('literature-review').innerHTML=`<article class="platform-panel"><h2>${esc(text(job.request?.query,'문헌 검토 기록'))}</h2><p><strong>${esc(statusLabel(job.status))}</strong> · 실제 처리 ${seconds(job.compute_seconds)}</p><p>기록 시점: ${esc(dateText(job.created_at))} · 모델 ${esc(text(job.request?.model))} · 추론 ${esc(text(job.request?.reasoning_effort))}</p>${job.status==='running'?`<p role="status">${job.progress?.stage==='literature_review'?'ChatGPT가 수집 문헌을 요약하고 있습니다.':'공개 문헌과 검토 자료를 준비하고 있습니다.'}${Number.isInteger(job.progress?.review_timeout_seconds)?' AI 응답 한도 '+esc(job.progress.review_timeout_seconds)+'초.':''} 목록 새로고침은 같은 작업의 상태만 확인합니다.</p>`:''}${job.status==='failed'?'<p>실행이 실패했습니다. 자동 재시도·유료 대체·작업 재개는 실행하지 않습니다.</p>':''}${dossier?literatureMarkup(dossier):'<p>문헌 근거 기록이 아직 제공되지 않았습니다. 작업 상태와 과학 판정은 별개입니다.</p>'}${dossier?'<button id="literature-save" class="button" type="button">근거 JSON 저장</button><p>JSON은 응답의 근거 기록이며 원본 문헌 바이트와 구분됩니다.</p>':''}${bound.length?`<h3>저장한 문헌 발췌문</h3>${bound.map((source,index)=>`<button class="button" type="button" data-literature-source="${index}">발췌문 보기 · ${esc(source.title)}</button>`).join('')}`:''}<pre id="literature-source" role="status" tabindex="0" aria-label="저장 문헌 본문과 검증 상태" style="white-space:pre-wrap;overflow-wrap:anywhere;max-height:28rem;overflow:auto"></pre></article>`;
  for(const button of main.querySelectorAll('[data-literature-source]'))button.onclick=()=>viewSource(id,bound[Number(button.dataset.literatureSource)]);
  const save=$('literature-save');if(save)save.onclick=async()=>{if(save.disabled)return;save.disabled=true;const version=detailVersion;try{const value=await api(BASE+'/reviews/'+id+'/dossier');if(!current()||version!==detailVersion)return;if(!record(value)||value.schema!=='skinpipe-literature-evidence-v1')throw Error('Invalid dossier');clearDownload();preparedDownload=saveDossier(value,id,{anchor:save});message('문헌 JSON이 준비되었습니다. 결과의 저장 링크를 눌러 받으세요.')}catch{message('근거 JSON을 확인하거나 저장하지 못했습니다. 실제 근거 기록을 다시 확인하세요.')}finally{if(current()&&version===detailVersion)save.disabled=false}};
  controls();
 }
 async function viewSource(id,expected){
  if(!current()||jobID(selected)!==id||!expected)return;const version=detailVersion,sourceTicket=++sourceVersion,panel=$('literature-source');panel.textContent='저장 원문과 해시를 확인하고 있습니다.';
  try{const dto=await api(BASE+'/reviews/'+id+'/source?path='+encodeURIComponent(expected.path));const source=await verifyLiteratureSource(dto,expected,cryptoImpl);if(!current()||version!==detailVersion||sourceTicket!==sourceVersion)return;panel.textContent='저장 경로·크기·SHA256 일치 확인: '+source.path+'\n'+source.text;}
  catch{if(current()&&version===detailVersion&&sourceTicket===sourceVersion)panel.textContent='원문 경로·크기·SHA256을 확인하지 못해 본문을 표시하지 않습니다. 저장 기록을 확인하세요.'}
 }
 async function readReview(id,version=detailVersion){
  if(!current()||!JOB.test(id))return;
  try{const job=await api(BASE+'/reviews/'+id);if(!current()||version!==detailVersion)return;if(jobID(job)!==id)throw Error('Mismatched review identity');const wasPending=pending(selected);showReview(job);nextPoll(id,version);if(wasPending&&!pending(job)){const refreshed=await Promise.allSettled([usage(),listReviews({openLatest:false})]);if(current()&&version===detailVersion&&refreshed[1].status==='rejected')message('작업 종료는 확인했으나 저장 목록 갱신을 확인하지 못했습니다. 목록만 새로고침하세요. 자동 재접수·재개는 하지 않았습니다.');}}
  catch(error){if(current()&&version===detailVersion){cancelPoll();message(errorText(error)+' 자동 상태 확인을 멈췄습니다. 저장 목록에서 직접 다시 확인하세요.')}}
 }
 async function openReview(id){if(!JOB.test(id)||!current())return;cancelPoll();polls=0;const version=++detailVersion;await readReview(id,version)}
 async function listReviews({openLatest=true}={}){
  if(!current()||!configured())return;
  const ticket=++listVersion,out=await api(BASE+'/reviews');if(!current()||ticket!==listVersion)return;if(!record(out)||!Array.isArray(out.reviews))throw Error('Invalid review list');
  const jobs=rows(out.reviews).filter(job=>jobID(job)).slice(0,100);
  $('literature-list').innerHTML=jobs.length?jobs.map((job,index)=>`<article><button class="button" type="button" data-literature-review="${index}">${esc(text(job.request?.query,'표적 미기록'))} · ${esc(statusLabel(job.status))}</button><p>${esc(dateText(job.created_at))}</p></article>`).join(''):'<p>저장한 문헌 검토 작업이 없습니다.</p>';
  for(const button of main.querySelectorAll('[data-literature-review]'))button.onclick=()=>openReview(jobID(jobs[Number(button.dataset.literatureReview)]));
  if(openLatest&&!selected&&jobs.length)await openReview(jobID(jobs[0]));
 }
 async function refresh(){
  if(!current()||refreshing)return;refreshing=true;controls();
  try{
   if(offline){$('literature-account').textContent='오프라인 저장본에서는 개인 ChatGPT 문헌 검토를 실행하지 않습니다. 분석 서버를 명시적으로 연결해 주세요.';return}
   const value=await api(BASE+'/status');if(!current())return;connection=record(value)?value:null;
   $('literature-account').textContent=!configured()?'이 서버에는 개인 ChatGPT 플랜 문헌 검토가 설정되지 않았습니다. 공유 Site 로그인만으로 플랜이 연결되지는 않습니다.':connection.connected!==true?'개인 ChatGPT 계정 연결을 확인하지 못했습니다. 저장한 문헌 기록은 아래에서 조회할 수 있습니다.':connection.plan_enabled!==true?'개인 계정 연결 기록은 있으나 플랜 사용 동의가 확인되지 않았습니다.':'본인 ChatGPT 플랜 연결·사용 동의 기록 확인 · 현재 응답 가능 여부와 잔여 한도는 미확인입니다.';
   if(configured())await Promise.all([usage(),loadModels(),listReviews()]);
  }catch(error){connection=null;if(current()){$('literature-account').textContent=errorText(error);message(errorText(error))}}
  finally{refreshing=false;controls()}
 }
 async function submit(event){
  event?.preventDefault();if(!current()||submitting||pending(selected))return;
  if(!allowed()||!catalog.has($('literature-model').value)){message('개인 계정·플랜 사용 동의와 실제 모델 목록을 확인하고 모델을 선택하세요.');return}
  if(!$('literature-consent').checked){message('본인 ChatGPT 플랜으로 공개 문헌을 검토할지 확인란에 표시해 주세요.');return}
  const body={query:$('literature-query').value.trim(),kind:$('literature-kind').value,taxon_id:Number($('literature-taxon').value),action:publicResearch&&$('literature-kind').value==='compound'?'any':$('literature-action').value,max_papers:Number($('literature-papers').value),model:$('literature-model').value,reasoning_effort:$('literature-effort').value};
  if(publicResearch&&body.kind==='compound'&&$('literature-topic').value.trim())body.public_topic=$('literature-topic').value.trim();
  const compound=publicResearch&&body.kind==='compound',validQuery=compound?/^[\p{L}\p{N}][\p{L}\p{N} ,().'’+\-]{0,511}$/u.test(body.query):/^[A-Za-z][A-Za-z0-9-]{0,159}$/.test(body.query),validTopic=!body.public_topic||/^[A-Za-z][A-Za-z0-9 -]{0,127}$/.test(body.public_topic)&&body.public_topic.split(/\s+/).length<=12;
  if(!validQuery||!validTopic||!(publicResearch?['protein','compound']:['protein','complex','endpoint']).includes(body.kind)||!Number.isSafeInteger(body.taxon_id)||body.taxon_id<1||body.taxon_id>100000000||!['any','inhibit','activate'].includes(body.action)||!Number.isInteger(body.max_papers)||body.max_papers<1||body.max_papers>8||(!EFFORTS.includes(body.reasoning_effort)||!literatureModelChoice(modelDetails,body.model).supported.includes(body.reasoning_effort))){message(publicResearch?'공개 화합물 이름 또는 단백질 식별자와 문헌 수를 확인하세요. SMILES·비공개 구조를 입력하지 마세요. 관심 주제는 영문 단어로 입력하세요.':'표적 식별자·종·문헌 수·추론 강도를 입력 범위에 맞게 확인하세요.');return}
  let receipt;try{receipt=pendingSubmission(store,pendingKey,body,uuid)}catch{message('접수 식별자를 생성하지 못했습니다. 새로고침한 뒤 다시 확인하세요. AI 요청은 보내지 않았습니다.');return}
  let accepted=false;submitting=true;controls();message('새 문헌 검토를 한 번 접수합니다. 응답이 늦어도 자동 재접수하지 않습니다.');
  try{const job=await api(BASE+'/reviews',{method:'POST',body:JSON.stringify({...body,idempotency_key:receipt.id})});if(!current())return;if(!jobID(job))throw Error('Invalid accepted job');acknowledgeSubmission(store,pendingKey,receipt.id);accepted=true;cancelPoll();polls=0;++detailVersion;showReview(job);nextPoll(jobID(job),detailVersion);message('문헌 검토가 접수됐습니다. 대기 상태는 읽기 조회로 확인하며 새 계산·재접수는 하지 않습니다.');await Promise.all([usage(),listReviews({openLatest:false})]);}
  catch(error){if(current())message(accepted?'문헌 검토는 접수됐으나 저장 목록 갱신을 확인하지 못했습니다. 아래 작업 기록을 확인하고 목록만 새로고침하세요. 자동 재접수하지 않습니다.':errorText(error)+' 동일한 입력을 직접 다시 접수하면 보존한 접수 키를 사용합니다.')}
  finally{submitting=false;controls()}
 }
 function kindControls(){if(!publicResearch)return;const compound=$('literature-kind').value==='compound';$('literature-query-label').textContent=compound?'공개 화합물 이름':'공개 단백질 유전자명';$('literature-query').placeholder=compound?'예: Nicotinamide':'예: TYR';$('literature-topic-label').hidden=!compound;$('literature-action-label').hidden=compound;if(compound)$('literature-action').value='any';$('literature-consent').checked=false;}
 if(publicResearch){$('literature-kind').onchange=kindControls;kindControls();}
 $('literature-form').onsubmit=submit;$('literature-model').onchange=()=>{modelControls();modelMessage()};$('literature-effort').onchange=modelMessage;$('literature-models').onclick=loadModels;$('literature-refresh').onclick=refresh;
 controls();
 return {start:refresh,refresh,submit,openReview,viewSource,stop(){stopped=true;++detailVersion;cancelPoll();clearDownload()}};
}
export async function renderPersonalLiterature(options){const controller=createPersonalLiteratureController(options);await controller.start();return controller}
