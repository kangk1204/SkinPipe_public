import {renderDiscovery} from './discovery.js';
import {renderPersonalLiterature} from './personal-literature.js';
import {renderPlatform} from './platform.js';
import {requestJson,requestTimeout,requestErrorText} from './api-request.js';
import {escapeHtml as e} from './research-inputs.js';

export function localRequestTimeout(path){
 if(path==='/v2/research-inputs/resolve')return 120000;
 return requestTimeout(path);
}
export function createLocalResearchApi({fetcher=globalThis.fetch,signalFactory=ms=>AbortSignal.timeout(ms),document=globalThis.document}={}){
 const api=async function(path,options={}){
  if(typeof path!=='string'||!/^\/v2\/[A-Za-z0-9_/?=&.%:-]+$/.test(path)||path.includes('..'))throw Error('SkinPipe 내부 API 경로가 필요합니다.');
  const request={...options,credentials:'same-origin',mode:'same-origin',redirect:'error',headers:{'Content-Type':'application/json','X-SkinPipe-Local':'1',...(options.headers||{})}};
  const timeout=signalFactory(localRequestTimeout(path));request.signal=options.signal&&typeof AbortSignal.any==='function'?AbortSignal.any([timeout,options.signal]):options.signal||timeout;
  if(options.rawResponse){delete request.rawResponse;const response=await fetcher(path,request);if(!response.ok){let message;try{message=(await response.json()).error}catch{}throw Error(message||'다운로드 응답을 받지 못했습니다.')}return response}
  return requestJson(fetcher,path,path,request);
 };
 // Same-origin HttpOnly cookies allow the browser's native streaming download
 // without buffering a large ZIP in JavaScript or requiring a save-picker API.
 if(document?.createElement)api.downloadExport=({exportId,filename,anchor})=>{
  if(typeof exportId!=='string'||!/^[0-9a-f]{32}$/.test(exportId))throw Error('확인된 내보내기 작업이 필요합니다.');
  const link=document.createElement('a');link.href='/v2/exports/'+exportId+'/download';link.download=filename||'SkinPipe-results.zip';
  document.getElementById?.('local-ready-export')?.remove();
  link.id='local-ready-export';link.className='button primary';link.textContent='준비된 Excel 묶음 ZIP 저장';
  if(anchor?.insertAdjacentElement)anchor.insertAdjacentElement('afterend',link);else document.body.append(link);
  link.focus?.();return link;
 };
 return api;
}
function statusText(value){
 if(typeof value==='boolean')return value?'사용 가능':'아직 준비되지 않음';
 if(typeof value==='string')return({ready:'기본 구성 준비됨',not_configured:'아직 설정되지 않음',unavailable:'현재 사용할 수 없음',pending:'연결 확인 중',connected:'연결 기록 있음',idle:'연결 전'})[value]||value;
 if(value&&typeof value==='object'){
  if(typeof value.connected==='boolean')return value.connected?'연결 기록 있음':value.configured||value.connect_path?'선택 사항 · 연결 전':'선택 사항 · 미설정';
  if(typeof value.compound_count==='number')return value.compound_count.toLocaleString('ko-KR')+'개 구조 · '+(value.structure_admission==='reported_verified'?'구조 검사 기록 있음':'구조 검사 상태 확인 필요');
  if(typeof value.engine==='string')return value.engine+' 경로 등록됨 · 실제 실행 여부는 처리 기록에서 확인';
  return value.message||statusText(value.status||value.reason)||(value.configured===false?'아직 설정되지 않음':'상세 상태 확인');
 }
 return'상태 미확인';
}
export function setupStatusMarkup(status){
 const fields=[['기본 실행',status.status??status.ready],['연구 데이터베이스',status.database||status.postgres],['화합물 카탈로그',status.catalog||status.local_catalog],['구조 계산',status.docking||status.engines],['선택적 ChatGPT 연결',status.chatgpt||status.literature]];
 return`<div class="page-heading"><div><h1>설치 상태</h1><p>실제로 사용할 수 있는 구성과 아직 필요한 준비를 확인합니다.</p></div></div><section class="platform-panel"><ul class="local-research-status-list">${fields.map(([label,value])=>`<li><strong>${e(label)}</strong> · ${e(statusText(value))}</li>`).join('')}</ul>${(Array.isArray(status.next_actions)?status.next_actions:[]).map(x=>`<p>${e(typeof x==='string'?x:x.message||x.action||'추가 준비 사항을 확인해 주세요.')}</p>`).join('')}<p>ChatGPT 연결은 선택 사항입니다. 기본 입력 확인·탐색·피부 평가·Excel 출력은 해당 계산 구성의 준비 상태에 따라 진행합니다.</p><details><summary>설치 검사 근거</summary><pre>${e(JSON.stringify(status,null,2))}</pre></details><a class="button primary" href="#new-research">연구 화면으로</a></section>`;
}
export function additionalToxicityStatusMarkup(profile){
 const value=profile?.additional_toxicity;
 const state={available:'설정됨 · 실제 실행 여부는 결과에서 확인',unconfigured:'미설정',invalid:'설정 확인 필요'}[value?.state]||'준비 상태 확인 불가';
 const version=value?.state==='available'&&typeof value.model_version==='string'&&value.model_version.trim()?` (모델 ${value.model_version.trim().slice(0,200)})`:'';
 return`<p class="research-muted" role="status" aria-label="실행 환경 요약">추가 독성 예측: ${e(state+version)}</p>`;
}
export function savedPopulationNoticeMarkup(){
 return'<section class="local-research-connection" role="note"><h2>전체 저장 행과 새 후보를 구분해 주세요</h2><p>이 표의 전체 모집단은 출발구조·측정 참조·관측 구조의 재도킹 대조를 포함한 저장 행 수입니다. <strong>새 후보 수는 연구 결과 요약의 역할별 집계를 기준으로 확인해 주세요.</strong></p><p>출발구조는 새로운 후보나 측정된 활성 대조군이 아닙니다. 표의 우선순위·탐색 점수는 계산용 휴리스틱이며, 활성·결합·피부 효능을 확정하지 않습니다. 실험 배치를 만들기 전에 개별 물질의 역할과 미평가 항목을 확인해 주세요.</p><a class="button" href="#new-research">연구 요약과 실행 기록으로</a></section>';
}
export async function startLocalResearch({document=globalThis.document,window=globalThis.window,api=createLocalResearchApi(),storage=globalThis.sessionStorage}={}){
 const main=document.getElementById('local-research-main'),banner=document.getElementById('local-research-banner');
 let revision=0,literature=null,setup=null,connectionTimer=null,connectionPolls=0,authorizationUrl=null;
 const current=version=>version===revision;
 async function status(){setup=await api('/v2/setup/status');return setup}
 async function connection(version){
  let result={};try{result=await api('/v2/literature/status')}catch{}if(!current(version))return;
  const capability=setup?.chatgpt||{},loginPath=capability.connect_path;
  const ready=typeof loginPath==='string'&&loginPath==='/v2/setup/chatgpt/connect';
  const configured=result.configured===true||capability.configured===true||ready,connected=result.connected===true||capability.connected===true;
  const pending=capability.signin_status==='pending'||capability.login_in_progress===true;
  if(!pending){authorizationUrl=null;connectionPolls=0}
  const errorCode=typeof capability.error==='string'?capability.error:null;
  main.innerHTML=`<div class="page-heading"><div><h1>ChatGPT 연결 (선택)</h1><p>기본 계산은 연결 없이도 진행합니다. 연결하면 공개 문헌 검토와 동의한 연구 데이터 설명을 추가로 사용할 수 있습니다.</p></div></div><section class="local-research-connection"><h2>${connected?'연결 기록이 있습니다':pending?'브라우저에서 ChatGPT 연결을 확인해 주세요':configured?'ChatGPT 계정을 연결해 주세요':'아직 선택적 연결이 설정되지 않았습니다'}</h2><p>${connected?'현재 계정의 플랜 사용 동의와 실제 응답 가능 여부를 확인한 뒤 문헌 검토를 실행하세요.':configured?'Continue with ChatGPT를 누르면 브라우저의 새 탭에서 본인 계정으로 연결에 동의합니다. 비밀번호나 API 토큰을 이 화면에 입력하지 않습니다.':'현재 구성에서 ChatGPT 도우미는 선택적 미설정 상태입니다. 연구를 먼저 진행하고 설치 상태에서 연결 준비 방법을 확인할 수 있습니다.'}</p><div id="local-chatgpt-connect">${ready&&!connected?`<button class="button primary" id="local-chatgpt-login" ${pending?'disabled':''}>Continue with ChatGPT</button>`:''}${pending?'<button class="button" id="local-chatgpt-cancel">이번 연결 취소</button>':''}${connected?'<button class="button" id="local-chatgpt-logout">ChatGPT 연결 해제</button>':''}</div><div class="form-actions"><a class="button primary" href="#new-research">기본 연구 계속하기</a>${connected?'<a class="button" href="#literature">공개 문헌 검토</a>':''}<button class="button" id="local-connection-refresh">연결 상태 갱신</button></div><p id="local-chatgpt-message" role="status" aria-live="polite">${pending?'연결 대기 중입니다. 동의가 끝나면 이 화면이 상태를 다시 확인합니다.':errorCode?'연결을 완료하지 못했습니다. 다시 연결하거나 기본 연구를 계속할 수 있습니다. ('+e(errorCode)+')':''}</p>${pending&&authorizationUrl?`<p><a href="${e(authorizationUrl)}" target="_blank" rel="noopener noreferrer">새 탭이 열리지 않았다면 여기에서 ChatGPT 연결 열기 ↗</a></p>`:''}</section>`;
  document.getElementById('local-connection-refresh').onclick=()=>route();
  if(ready&&!connected)document.getElementById('local-chatgpt-login').onclick=async()=>{
   const button=document.getElementById('local-chatgpt-login');if(button.disabled)return;button.disabled=true;
   let popup=null;try{popup=window.open('about:blank','_blank');if(popup)popup.opener=null}catch{}
   document.getElementById('local-chatgpt-message').textContent='브라우저 연결 화면을 준비하고 있습니다…';
   try{
    const response=await api(loginPath,{method:'POST',body:'{}'});
    if(response.authorization_url){
     const url=new URL(response.authorization_url);if(url.protocol!=='https:'||url.hostname!=='auth.openai.com'||url.username||url.password||url.port||url.hash)throw Error('허용된 ChatGPT 로그인 주소를 확인하지 못했습니다.');
     authorizationUrl=url.href;if(popup)popup.location.href=url.href;
    }else if(popup){popup.close();popup=null}
    if(current(version))await route();
   }catch(error){if(popup)try{popup.close()}catch{}if(current(version)){document.getElementById('local-chatgpt-message').textContent=requestErrorText(error)+' 연결 상태를 갱신해 접수 여부를 먼저 확인해 주세요.';button.disabled=false}}
  };
  if(pending)document.getElementById('local-chatgpt-cancel').onclick=async()=>{
   const button=document.getElementById('local-chatgpt-cancel');button.disabled=true;
   try{await api('/v2/setup/chatgpt/cancel',{method:'POST',body:'{}'});authorizationUrl=null;if(current(version))await route()}catch(error){if(current(version)){document.getElementById('local-chatgpt-message').textContent=requestErrorText(error);button.disabled=false}}
  };
  if(connected)document.getElementById('local-chatgpt-logout').onclick=async()=>{
   const button=document.getElementById('local-chatgpt-logout');button.disabled=true;
   try{await api('/v2/setup/chatgpt/logout',{method:'POST',body:'{}'});authorizationUrl=null;if(current(version))await route()}catch(error){if(current(version)){document.getElementById('local-chatgpt-message').textContent=requestErrorText(error);button.disabled=false}}
  };
  if(pending&&connectionPolls++<400)connectionTimer=setTimeout(()=>{if(current(version))route()},1500);
 }
 async function route(){
  const version=++revision;if(connectionTimer){clearTimeout(connectionTimer);connectionTimer=null}literature?.stop();literature=null;banner.textContent='';
  try{
   const hash=window.location.hash.split('?')[0];
   if(hash==='#setup'){main.innerHTML='<p>설치 구성과 실행 준비 상태를 확인하고 있습니다…</p>';const value=await status();if(current(version))main.innerHTML=setupStatusMarkup(value);return}
   if(hash==='#connection'){try{await status()}catch{}await connection(version);return}
   if(hash==='#literature'){
    const parameters=new URLSearchParams(window.location.hash.split('?')[1]||'');
    const controller=await renderPersonalLiterature({main,api,storage,publicResearch:true,legacyHistory:false,isCurrent:()=>current(version),prefill:{query:parameters.get('query')||'',fromDiscovery:parameters.get('source')==='discovery'}});if(current(version))literature=controller;else controller.stop();return;
   }
   if(hash==='#research'){await renderPlatform({main,api,route:'research',screens:[]});if(current(version))main.insertAdjacentHTML('afterbegin',savedPopulationNoticeMarkup());return}
   await renderDiscovery(main,api,'researcher',{isCurrent:()=>current(version)});
   if(!current(version))return;
   let profile;try{profile=await api('/v2/local/profile')}catch{}
   if(current(version))main.insertAdjacentHTML('afterbegin',additionalToxicityStatusMarkup(profile));
  }catch(error){if(current(version)){banner.textContent=requestErrorText(error);main.innerHTML='<section class="platform-panel"><h1>화면을 준비하지 못했습니다</h1><p>설치 상태에서 필요한 연결·계산 구성을 확인한 뒤 다시 시도해 주세요. 기존 결과 파일은 그대로 보존됩니다.</p><a class="button" href="#setup">설치 상태 확인</a></section>'}}
 }
 document.getElementById('local-setup-status').onclick=()=>{if(window.location.hash==='#setup')return route();window.location.hash='setup'};
 window.addEventListener('hashchange',route);
 await route();
 return{refresh:route,stop(){revision++;if(connectionTimer)clearTimeout(connectionTimer);literature?.stop();window.removeEventListener('hashchange',route)}};
}
if(globalThis.document?.getElementById('local-research-main'))startLocalResearch();
