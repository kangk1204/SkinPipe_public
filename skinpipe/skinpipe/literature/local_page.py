"""No external scripts, browser credentials, or browser persistence."""

PAGE = r'''<!doctype html><html lang="ko"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>SkinPipe · 문헌 근거</title>
<style>body{font:16px/1.65 system-ui,sans-serif;max-width:1000px;margin:40px auto;padding:0 22px;color:#183e3b;background:#f7faf9}section,article{background:white;border:1px solid #d8e6e0;border-radius:12px;padding:22px;margin:18px 0}h1,h2,h3{line-height:1.3}label{display:inline-block;margin:10px 20px 10px 0}input,select,button{font:inherit;padding:9px;border:1px solid #a7beb4;border-radius:6px}button{cursor:pointer;background:#edf5f1}button:disabled{cursor:wait;opacity:.6}input[type=checkbox]{width:18px;height:18px}input[type=text]{min-width:230px}a{color:#156856}pre{white-space:pre-wrap;word-break:break-word;background:#f4f7f6;padding:12px}blockquote{border-left:3px solid #87ab9a;padding-left:14px;margin-left:0}small{color:#526d65}.notice{border-left:4px solid #b7892a;padding-left:15px}[role=status]{white-space:pre-wrap}ul{padding-left:24px}button,a{margin-right:8px}#report:empty{display:none}</style>
<style>#usage{border:2px solid #658e7c;background:#f0f7f3}#usage>p>strong{font-size:1.4rem}#usage-detail>strong{display:block;font-size:1.2rem;background:#fff1d6;padding:12px;margin:12px 0;color:#65420b}#usage a,#usage button{display:inline-block;margin-top:8px}</style>
<h1>SkinPipe · 문헌 근거</h1><p>새 표적의 공개 문헌을 수집하고, 선택한 ChatGPT 모델로 근거를 검토합니다. 저장된 원본 분석과 후보 순위는 유지합니다.</p>
<section id="usage" aria-label="ChatGPT 사용량"><h2>ChatGPT 사용량</h2><p><strong>잔여 한도 조회 불가</strong></p><p>플랜의 남은 비율·토큰·초기화 시각은 이 연결에서 조회되지 않습니다. 이 앱의 기록 합계를 플랜 잔여량으로 환산하지 않습니다.</p><div id="usage-detail" role="status" aria-live="polite">사용량 기록 확인 중</div><button id="usage-refresh" type="button">사용량 기록 새로고침</button><a href="https://chatgpt.com/settings/usage" target="_blank" rel="noopener noreferrer">ChatGPT 사용량 설정 열기 ↗</a><a href="#literature-guide">연결·사용 방법 보기</a></section>
<section><h2>ChatGPT 계정 연결</h2><p>이 PC의 개인용 연결입니다. OpenAI 화면에서 로그인하고 플랜 사용에 동의합니다. AI 사용량은 본인 ChatGPT 플랜 한도에 포함됩니다. 계정 연결과 과학적 검증은 별개입니다.</p>
<button id="login">Continue with ChatGPT</button><button id="login-other" type="button">다른 계정으로 로그인</button><button id="login-cancel" type="button" hidden disabled>로그인 취소</button><button id="refresh">연결 상태 확인</button><button id="logout">연결 해제</button><a href="https://chatgpt.com/settings/usage" target="_blank" rel="noopener noreferrer">ChatGPT 사용량 설정 열기</a><p id="account" role="status">연결 상태 확인 중</p><p>다른 계정으로 로그인하면 OpenAI 화면에서 새 계정과 플랜 권한을 선택합니다. 새 연결 검증이 끝날 때까지 현재 연결을 유지합니다. 사용량은 활성 계정별로 조회하며, 계정 변경 후 모델을 다시 조회합니다.</p>
<label>문헌 검토 모델 <select id="model"><option value="">모델 조회 후 선택 · 선호 gpt-6-luna</option></select></label><button id="models">사용 가능한 모델 조회</button><label>추론 강도 <select id="reasoning"><option value="none">none</option><option value="low">low</option><option value="medium">medium</option><option value="high">high</option><option value="xhigh">xhigh</option><option value="max" selected>max</option></select></label><small>선호 설정은 gpt-6-luna · max입니다. Luna가 없으면 제공자가 명시한 기본 모델을 표시합니다. 기본 모델을 알 수 없으면 목록에서 직접 선택합니다. 실패 후 자동 전환·API 과금 전환은 하지 않습니다. 실제 요청이 완료되어야 모델·추론 강도 접근이 확인됩니다.</small></section>
<section><h2>새 표적 문헌 수집·검토</h2><form id="collect"><label>표적 식별자 <input id="query" type="text" placeholder="예: GSK3B" maxlength="160" required></label><label>종 ID <input id="taxon" type="number" value="9606" min="1" required></label><label>표적 종류 <select id="kind"><option value="protein">단백질</option><option value="endpoint">실험 지표</option></select></label><label>원하는 작용 <select id="action"><option value="any">미지정</option><option value="inhibit">억제</option><option value="activate">활성화</option></select></label><label>최대 문헌 수 <input id="count" type="number" value="8" min="1" max="8" required></label><label><input id="review" type="checkbox"> ChatGPT로 근거 검토 (공개 문헌 텍스트만 전송)</label><p>체크하지 않으면 문헌 수집과 출처·원문 대조 검사만 실행합니다. 문헌 검색은 제한된 범위이며 전체 문헌 검토를 뜻하지 않습니다. 검토를 선택하면 입력한 표적 식별자와 수집한 공개 문헌 텍스트를 전송합니다. 후보 구조와 실험 데이터는 이 검토 요청에 포함하지 않습니다.</p><button type="submit">새 근거 기록 만들기</button></form><p id="message" role="status"></p></section>
<section><h2>저장한 근거</h2><button id="reports">목록 새로고침</button><div id="list"></div></section><section id="report"></section>
<p id="literature-guide" class="notice">먼저 Continue with ChatGPT를 눌러 이 PC의 계정을 연결하고 플랜 사용에 동의한 뒤 모델을 선택하세요. 공개 문헌만 수집하려면 ChatGPT 검토를 체크하지 않아도 됩니다. 공유 사이트 로그인과 이 PC의 ChatGPT 플랜 연결은 별개이며 Codex 계정을 재사용하지 않습니다. 자동 검토는 사람의 문헌·표적·실험 조건 검토를 대신하지 않습니다. 도킹·활성 방향·효능·프로젝트 실험 검증은 해당 입력과 실제 검증 영수증이 있어야 평가합니다. 이 화면에서 실험 검증은 미평가로 유지합니다. 공유 웹사이트의 계정 로그인과 권한 부여는 별도 운영 설정이 필요합니다.</p>
<script nonce="__CSP_NONCE__">
const $=id=>document.getElementById(id),csrf='__NONCE__';
const record=value=>value!==null&&typeof value==='object'&&!Array.isArray(value);
const text=(value,fallback='미기록')=>typeof value==='string'&&value.trim()?value:fallback;
const limited=(value,max=12000)=>{const out=text(value);return out.length>max?out.slice(0,max)+'\n[표시 한도 초과 · 전체는 근거 JSON을 확인하세요.]':out};
const rows=value=>Array.isArray(value)?value.filter(record):[];
const preferredModel='gpt-6-luna',reasoningEfforts=['none','low','medium','high','xhigh','max'];
let modelMetadata=new Map();
function modelReasoning(){
 const id=$('model').value,row=modelMetadata.get(id),supported=Array.isArray(row?.supported_reasoning_efforts)?reasoningEfforts.filter(e=>row.supported_reasoning_efforts.includes(e)):reasoningEfforts;
 $('reasoning').replaceChildren(node('option','추론 강도를 선택하세요'));$('reasoning').firstChild.value='';
 for(const e of supported)node('option',e,$('reasoning')).value=e;
 $('reasoning').value=id===preferredModel&&supported.includes('max')?'max':supported.includes(row?.default_reasoning_effort)?row.default_reasoning_effort:'';
}
const statuses={HOLD_USAGE_LIMIT:'이 앱의 ChatGPT 사용량 한도 도달 · AI 검토 보류',HOLD_RATE_LIMIT:'일시적 요청 속도 제한 · AI 검토 보류',NOT_RUN:'문헌 수집 미실행',HOLD:'문헌 수집·검토 보류',HOLD_UNAVAILABLE:'수집 또는 검토 연결 보류',HOLD_REVIEW_UNAVAILABLE:'문헌 검토 연결 실패 · 제안 보류',HOLD_PARTIAL_RETRIEVAL:'문헌 수집 일부 실패 · 검토 보류',HOLD_UNSUPPORTED_TARGET:'지원 범위 밖 표적 · 문헌 수집 보류',HOLD_NO_SOURCES:'제한된 검색에서 문헌을 확보하지 못해 보류',COLLECTED_HUMAN_REVIEW_REQUIRED:'문헌 수집됨 · 사람 검토 필요',HOLD_INVALID_REVIEW:'제안 형식·인용 검사 미충족 · 보류',REVIEW_PROPOSED_HUMAN_REVIEW_REQUIRED:'LLM 해석 제안 있음 · 사람 검토 필요'};
const coverages={ABSTRACT_ONLY:'초록만 확보',ABSTRACT_EXCERPT:'초록 일부만 확보',OA_FULL_TEXT:'공개 원문 확보',OA_FULL_TEXT_EXCERPT:'공개 원문 일부만 확보',METADATA_ONLY:'메타데이터만 확보 · 본문 없음'};
const label=(map,key,fallback)=>typeof key==='string'&&Object.hasOwn(map,key)?map[key]:fallback;
function jsonText(value){try{return limited(JSON.stringify(value??null,null,2),30000)}catch{return '기록 형식 미확인'}}
async function request(path,body){
 const response=await fetch(path,body===undefined?{cache:'no-store'}:{method:'POST',headers:{'Content-Type':'application/json','X-SkinPipe-CSRF':csrf},body:JSON.stringify(body)});
 let value;try{value=await response.json()}catch{throw Error('서버 응답을 읽을 수 없습니다. 연결 상태를 확인해 주세요.')}
 if(!response.ok){const error=Error(record(value)?text(value.error,'요청이 완료되지 않았습니다.'):'요청이 완료되지 않았습니다.');const code=record(value)?value.code||value.error_code:null;if(typeof code==='string'&&/^[A-Za-z][A-Za-z0-9_]{0,63}$/.test(code))error.code=code;error.http_status=response.status;if(record(value)){if(typeof value.retry_after_seconds==='number'&&Number.isFinite(value.retry_after_seconds)&&value.retry_after_seconds>=0&&value.retry_after_seconds<=31536000)error.retry_after_seconds=value.retry_after_seconds;if(typeof value.reset_at==='string'&&validResetAt(value.reset_at))error.reset_at=value.reset_at;if(typeof value.usage_exhausted==='boolean')error.usage_exhausted=value.usage_exhausted}throw error}return value;
}
function node(tag,value,parent){const el=document.createElement(tag);el.textContent=String(value??'');if(parent)parent.append(el);return el}
let taskBusy=false,operationVersion=0,loginVersion=0,loginPending=false,localLoginPending=false,cancelPending=false;
function controls(){
 const locked=taskBusy||loginPending;
 document.querySelectorAll('button').forEach(button=>button.disabled=locked);
 const cancel=$('login-cancel');cancel.hidden=!loginPending;cancel.disabled=!loginPending||cancelPending;
}
function busy(value){taskBusy=value;controls()}
function validResetAt(value){return typeof value==='string'&&value.length<=120&&/^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(?:\.\d+)?Z$/.test(value)&&Number.isFinite(Date.parse(value))}
function quotaTiming(error){
 const messages=[];
 if(record(error)&&typeof error.retry_after_seconds==='number'&&Number.isFinite(error.retry_after_seconds)&&error.retry_after_seconds>=0&&error.retry_after_seconds<=31536000)messages.push('제공된 재시도 대기: '+error.retry_after_seconds+'초');
 if(record(error)&&validResetAt(error.reset_at))messages.push('제공된 한도 갱신 시각: '+error.reset_at);
 return messages;
}
function quotaKind(value){return ['usage_exhausted','usage_limit_exceeded','subscription_sharing_usage_limit_exceeded'].includes(value?.code)||value?.usage_exhausted===true?'usage':['rate_limit_exceeded','rate_limited'].includes(value?.code)?'rate':null}
const quotaExplanation='수집한 문헌·출처 기록은 보존하며 AI 검토만 보류합니다. 자동 유료 대체 경로로 전환하거나 반복 호출하지 않습니다. 사용 가능해진 뒤 직접 다시 시도하세요.';
function errorMessage(error,fallback='요청이 완료되지 않았습니다.'){
 const code=typeof error?.code==='string'&&/^[A-Za-z][A-Za-z0-9_]{0,63}$/.test(error.code)?'['+error.code+'] ':'';
 const kind=quotaKind(error);
 if(kind)return code+(kind==='usage'?'이 앱에서 사용할 수 있는 ChatGPT 사용량 한도에 도달해 AI 검토를 보류했습니다. ':'일시적 요청 속도 제한 · AI 검토 보류. ')+quotaExplanation+(quotaTiming(error).length?' '+quotaTiming(error).join(' · '):'');
 return code+text(error?.message,fallback);
}
function showQuota(dossier,panel){
 const error=record(dossier.review)&&record(dossier.review.error)?dossier.review.error:null;
 const kind=dossier.status==='HOLD_USAGE_LIMIT'?'usage':dossier.status==='HOLD_RATE_LIMIT'?'rate':quotaKind(error);
 if(!kind)return;
 node('h3',kind==='usage'?'이 앱의 ChatGPT 사용량 한도 도달 · AI 검토 보류':'일시적 요청 속도 제한 · AI 검토 보류',panel);
 if(quotaKind(error))node('p','오류 코드: '+error.code,panel);
 if(kind==='usage')node('p','이 앱에서 사용할 수 있는 ChatGPT 사용량 한도에 도달해 AI 검토를 보류했습니다.',panel);
 node('p',quotaExplanation,panel);for(const message of quotaTiming(error))node('p',message,panel);
}
let usagePending=null;
const usageCount=value=>Number.isSafeInteger(value)&&value>=0?value.toLocaleString('ko-KR'):'미기록';
function showUsage(value,errorStatus=null){
 const panel=$('usage-detail');panel.replaceChildren();
 if(!record(value)||value.schema!=='skinpipe-chatgpt-usage-v1'){
  node('p',errorStatus===401?'사용량 기록에 접근하려면 SkinPipe 연결·권한을 확인하세요.':errorStatus===404?'이 서버는 사용량 기록 조회를 아직 제공하지 않습니다.':errorStatus?'사용량 기록을 가져오지 못했습니다. 직접 새로고침하세요.':'연결된 사용량 기록이 없습니다.',panel);return;
 }
 node('p',value.account?.connected===true?'이 PC의 SkinPipe 계정 연결 기록이 있습니다. 실제 모델 접근·잔여 한도는 미확인입니다.':'이 PC의 ChatGPT 계정이 연결되지 않았습니다.',panel);
 const error=record(value.last_limit_error)?value.last_limit_error:null,kind=quotaKind(error);
 if(kind){const badge=node('strong',kind==='usage'?'최근 요청: 이 앱의 사용량 한도 도달 · AI 검토 대기':'최근 요청: 일시적 요청 속도 제한 · AI 검토 대기',panel);badge.className='notice';node('p','현재 회복 여부는 미확인입니다. '+quotaExplanation,panel);for(const message of quotaTiming(error))node('p',message,panel);if(validResetAt(error.observed_at))node('p','제한이 기록된 시각: '+error.observed_at,panel)}
 const app=record(value.app_usage)?value.app_usage:null;
 if(value.scope!=='current_registration_this_installation'||app?.accounting_status!=='available'){node('p','이 설치의 실제 AI 요청·토큰 사용량 기록을 아직 확인할 수 없습니다. 미기록 값을 0으로 간주하지 않습니다.',panel);return}
 node('p','이 앱 연결·이 설치에서 기록한 요청만 집계합니다. 기록 이전 요청, 다른 기기·앱·ChatGPT의 사용량은 포함하지 않습니다.',panel);
 node('p','시작한 AI 요청 '+usageCount(app.requests_started)+' · 응답 완료 '+usageCount(app.responses_completed)+' · 실패한 요청 '+usageCount(app.failed_requests)+' · 미종결 요청 '+usageCount(app.pending_requests),panel);
 node('p','응답 완료는 전송 완료 기록이며 문헌 검토 승인이나 과학적 검증을 뜻하지 않습니다.',panel);
 node('p','제공자가 보고한 토큰 합계: 입력 '+usageCount(app.input_tokens)+' · 출력 '+usageCount(app.output_tokens)+' · 전체 '+usageCount(app.total_tokens),panel);
 node('p','토큰 사용량 미제공 요청: '+usageCount(app.usage_missing_requests)+' · 기록 범위: '+(app.coverage==='COMPLETE_FOR_RECORDED_REQUESTS'?'기록된 요청 범위에서 완전':app.coverage==='PARTIAL'?'일부 기록만 있음':'확인 불가')+' · 미기록이 있으면 알려진 값만의 합입니다.',panel);
 if(validResetAt(app.tracking_since))node('p','기록 시작: '+app.tracking_since,panel);
}
async function refreshUsage(){
 if(usagePending)return usagePending;
 usagePending=Promise.resolve().then(()=>request('/api/usage')).then(value=>{showUsage(value);return value},error=>{showUsage(null,error?.http_status||'unavailable');return null}).finally(()=>{usagePending=null});
 return usagePending;
}
$('usage-refresh').onclick=()=>operation(refreshUsage,'사용량 기록을 확인합니다. AI 검토 요청은 보내지 않습니다.');

async function operation(fn,message='처리 중입니다.'){
 const ticket=++operationVersion;busy(true);$('message').textContent=message;
 try{const result=await fn();if(ticket===operationVersion)$('message').textContent=loginPending?'로그인이 진행 중입니다. 로그인 창을 닫았다면 로그인 취소를 누른 뒤 다시 시도하세요.':typeof result==='string'?result:'요청 처리가 끝났습니다. 아래 검토 상태와 필요한 후속 확인을 읽어 주세요.'}
 catch(error){if(ticket===operationVersion)$('message').textContent=errorMessage(error)}
 finally{if(ticket===operationVersion)busy(false)}
}
async function status(){
 const ticket=loginVersion,s=await request('/api/status');if(!record(s))throw Error('계정 상태 응답 형식 미확인');
 $('account').textContent=s.connected===true?'연결됨 · '+text(s.email,'본인 계정')+' · 플랜 사용 '+(s.plan_enabled===true?'동의 확인':'권한 없음'):'연결되지 않음';
 if(ticket===loginVersion){loginPending=localLoginPending||s.login_in_progress===true;controls();if(loginPending)$('account').textContent+=' · 로그인 진행 중 (창을 닫았다면 로그인 취소)'}
 await refreshUsage();
 return s;
}
async function startLogin(addAccount=false){
 if(loginPending||taskBusy)return;
 const ticket=++loginVersion;++operationVersion;localLoginPending=true;loginPending=true;cancelPending=false;controls();
 $('message').textContent='OpenAI 로그인 창에서 계정과 플랜 사용 권한을 확인해 주세요. 창을 닫았다면 로그인 취소를 눌러 다시 시도할 수 있습니다.';
 try{
  const connection=await request('/api/login',addAccount?{add_account:true}:{});if(ticket!==loginVersion)return;
  localLoginPending=false;const current=await status();
  if(addAccount&&connection?.connected===true&&current.connected===true){$('model').replaceChildren(node('option','계정 변경됨 · 사용 가능한 모델을 다시 조회하세요'));$('model').firstChild.value='';$('model').value=''}
  if(ticket===loginVersion)$('message').textContent=loginPending?'로그인이 아직 진행 중입니다. 창을 닫았다면 로그인 취소를 누르세요.':'로그인 요청이 끝났습니다. 연결 상태와 플랜 사용 동의를 확인하세요.';
 }catch(error){
  if(ticket!==loginVersion)return;localLoginPending=false;loginPending=false;controls();
  try{await status()}catch{}
  if(ticket===loginVersion)$('message').textContent=errorMessage(error)+' Continue with ChatGPT를 눌러 다시 시도하세요. 로그인 진행 중으로 남아 있으면 로그인 취소를 먼저 누르세요.';
 }finally{if(ticket===loginVersion){localLoginPending=false;controls()}}
}
async function cancelLogin(){
 if(!loginPending||cancelPending)return;
 const ticket=loginVersion;cancelPending=true;controls();$('message').textContent='로그인 취소 요청 중입니다…';
 try{
  const result=await request('/api/login-cancel',{});if(ticket!==loginVersion)return;
  if(record(result)&&result.cancelled===true){
   ++loginVersion;localLoginPending=false;loginPending=false;cancelPending=false;taskBusy=false;controls();
   $('message').textContent='[LOGIN_CANCELLED] 로그인 요청을 취소했습니다. 연결 상태를 확인하고 Continue with ChatGPT를 눌러 다시 시도하세요.';
   $('account').textContent='로그인 요청 취소됨 · 연결 상태 확인 중';try{await status()}catch{$('account').textContent='로그인 요청 취소됨 · 연결 상태 확인 필요'}
  }else if(record(result)&&result.login_in_progress===false){
   const finished=++loginVersion;localLoginPending=false;loginPending=false;cancelPending=false;taskBusy=false;controls();
   let current=null;try{current=await status()}catch{$('account').textContent='로그인 요청 종료됨 · 연결 상태 확인 필요'}
   if(finished===loginVersion)$('message').textContent=current?.connected===true?'[LOGIN_NOT_IN_PROGRESS] 진행 중인 로그인 요청이 없습니다. 현재 계정 연결은 유지됩니다. 다른 로그인이 필요하면 Continue with ChatGPT를 누르세요.':'[LOGIN_NOT_IN_PROGRESS] 로그인 요청이 이미 끝났습니다. 연결 상태를 확인하고 Continue with ChatGPT를 눌러 다시 시도하세요.';
  }else{
   loginPending=true;$('message').textContent='[LOGIN_CANCELLATION_PENDING] 로그인 종료를 기다리고 있습니다. 잠시 뒤 로그인 취소를 다시 누르거나 연결 상태를 확인하세요.';
  }
 }catch(error){if(ticket===loginVersion)$('message').textContent=errorMessage(error)+' 로그인 취소를 다시 시도하세요.'}
 finally{if(ticket===loginVersion){cancelPending=false;controls()}}
}
$('refresh').onclick=()=>operation(status);
$('model').onchange=modelReasoning;
$('login').onclick=()=>startLogin(false);
$('login-other').onclick=()=>startLogin(true);
$('login-cancel').onclick=cancelLogin;
$('logout').onclick=()=>operation(async()=>{const out=await request('/api/logout',{});await status();if(record(out)&&out.remote_revocation_confirmed===false)throw Error('로컬 연결은 해제했습니다. 원격 권한 취소는 확인되지 않았으므로 ChatGPT 설정에서도 연결을 확인해 주세요.')});
$('models').onclick=()=>operation(async()=>{
 const out=await request('/api/models');$('model').replaceChildren(node('option','모델을 선택하세요'));$('model').firstChild.value='';$('model').value='';
 const models=Array.isArray(out)?out:Array.isArray(out?.models)?out.models:Array.isArray(out?.models?.data)?out.models.data:Array.isArray(out?.data)?out.data:[];
 let count=0,preferredAvailable=false;modelMetadata=new Map();
 for(const model of models){const id=typeof model==='string'?model:record(model)?text(model.id,text(model.slug,'')):'';if(!id||id.length>200||modelMetadata.has(id))continue;modelMetadata.set(id,record(model)?model:{});const option=node('option',id,$('model'));option.value=id;count++;if(id===preferredModel)preferredAvailable=true}
 const declared=record(out)&&typeof out.default_model==='string'&&models.some(m=>record(m)&&(m.id??m.slug)===out.default_model)?out.default_model:'';
 if(preferredAvailable)$('model').value=preferredModel;else if(declared)$('model').value=declared;
 modelReasoning();
 return count?'모델 '+count+'개 표시 · '+(preferredAvailable?'선호 모델 gpt-6-luna 선택됨.':declared?'Luna가 없어 제공자가 명시한 기본 모델 '+declared+' 선택됨. 실행 전 모델·추론을 확인하세요.':'선호 모델 gpt-6-luna와 제공자 기본 모델을 확인하지 못했습니다. 목록에서 직접 선택하세요. 임의 모델은 자동 선택하지 않습니다.')+' 실제 요청 완료 전까지 접근 권한은 미확인입니다.':'모델 목록을 확보하지 못했습니다. 연결·권한과 제공자의 상태를 확인하세요. 사용할 모델이 없다고 확정한 것은 아닙니다.';
});
function literatureURL(value){
 if(typeof value!=='string'||value.length>4096||/[\u0000-\u0020\u007f\\]/.test(value))return null;
 try{const u=new URL(value);if(u.protocol!=='https:'||u.username||u.password||u.port||u.hash)return null;const host=u.hostname,path=u.pathname;
  if(host==='doi.org'&&/^\/10\.\d{4,9}\/\S+$/i.test(path))return u.href;
  if(host==='pubmed.ncbi.nlm.nih.gov'&&/^\/\d+\/$/.test(path)&&!u.search)return u.href;
  if(host==='pmc.ncbi.nlm.nih.gov'&&/^\/articles\/PMC\d+\/?$/.test(path)&&!u.search)return u.href;
  if(host==='www.ncbi.nlm.nih.gov'&&/^\/pmc\/articles\/PMC\d+\/?$/.test(path)&&!u.search)return u.href;
  if(['europepmc.org','www.europepmc.org'].includes(host)&&/^\/articles?\/(?:MED\/\d+|PMC\/PMC\d+)\/?$/i.test(path)&&!u.search)return u.href;
  if(host==='www.ebi.ac.uk'&&/^\/europepmc\/webservices\/rest\/(?:search|PMC\d+\/fullTextXML)$/.test(path))return u.href;
 }catch{}return null;
}
function sourceLink(parent,url,value){const safe=literatureURL(url);if(safe){const a=node('a',limited(value,1500),parent);a.href=safe;a.target='_blank';a.rel='noopener noreferrer';return}node('span',limited(value,1500)+' · 문헌 링크 미확인',parent)}
function reviewRequired(parent){node('p','사람 검토 필요 · LLM 제안이나 인용문 포함 검사는 과학적 타당성·효능 검증을 대신하지 않습니다.',parent);node('p','실제 프로젝트 실험 검증: NOT_EVALUATED',parent)}
function showReport(d,id){
 const panel=$('report');panel.replaceChildren();
 if(!record(d)||d.schema!=='skinpipe-literature-evidence-v1'){node('h2','문헌 근거 · 형식 미확인',panel);node('p','지원하는 기록 계약을 확인할 수 없어 검증 미완료로 표시합니다.',panel);reviewRequired(panel);return}
 node('h2','분석 근거 · '+text(d.target?.query,text(d.target?.request?.query,text(id))),panel);
 node('p','기록 상태: '+label(statuses,d.status,'상태 미확인 · 검증 미완료'),panel);
 if(typeof d.reason==='string'&&d.reason.trim())node('p','사유: '+limited(d.reason,6000),panel);
 node('p','조회일: '+text(d.retrieved_at,text(d.created_at,'기록 JSON 확인')),panel);reviewRequired(panel);showQuota(d,panel);
 if(d.human_review_required!==true)node('p','사람 검토 필수 표식 미확인 · 승인된 근거로 해석할 수 없습니다.',panel);
 if(typeof id==='string'&&/^[0-9a-f]{32}$/.test(id)){const a=node('a','근거 JSON 저장',panel);a.href='/api/report?id='+encodeURIComponent(id);a.download='SkinPipe-literature-'+id+'.json'}
 node('h3','원 출처와 읽은 범위',panel);
 const sources=rows(d.sources),claims=rows(d.claims),ids=new Map();
 for(const source of sources){const key=text(source.source_id,'');if(key)ids.set(key,(ids.get(key)||0)+1)}
 for(const source of sources.slice(0,100)){
  const article=node('article','',panel);sourceLink(article,source.url,text(source.title,text(source.source_id)));
  node('p','출처 ID: '+text(source.source_id)+' · PMID '+text(source.pmid)+' · DOI '+text(source.doi),article);
  node('p','읽은 범위: '+label(coverages,source.coverage,'범위 미확인')+' · '+text(source.coverage)+(source.review_text_truncated===true?' · 검토용 본문이 잘려 있음':''),article);
  node('p',(source.retracted===true?'철회 표시 있음':source.retracted===false?'수집 기록상 철회 표시 없음':'철회 여부 미확인')+' · 확인 기록: '+limited(source.retraction_status,1000),article);
  node('p','문헌 분류: '+limited(source.article_class,1000)+' · 모델·시험 범위: '+limited(source.model_scope,3000),article);
  node('small','텍스트 SHA256: '+(typeof source.text_sha256==='string'&&/^[a-f0-9]{64}$/.test(source.text_sha256)?source.text_sha256:'유효한 본문 해시 미기록'),article);
 }
 if(!sources.length)node('p','문헌 출처 기록 없음 · 수집 성공이나 문헌 부재를 뜻하지 않습니다.',panel);
 if(sources.length>100)node('p','최초 100개 출처만 표시합니다. 전체는 근거 JSON을 확인하세요.',panel);
 node('h3','자동 검토 제안 · 사람 검토 전',panel);
 for(const claim of claims.slice(0,100)){
  const article=node('article','',panel);node('p','LLM 제안 · 사람 검토 필요',article);node('p',limited(claim.claim),article);
  const sourceId=text(claim.source_id,''),count=ids.get(sourceId)||0;
  node('p','출처: '+(sourceId||'미기록')+' · '+(count===1?'출처 연결됨':count>1?'출처 ID 중복 · 연결 미확인':'출처 연결 미확인'),article);
  node('pre',limited(claim.quote,20000),node('blockquote','',article));
  if(Number.isSafeInteger(claim.quote_start)&&claim.quote_start>=0&&Number.isSafeInteger(claim.quote_end)&&claim.quote_end>=claim.quote_start)node('p','기록된 문자 위치 ['+claim.quote_start+', '+claim.quote_end+')',article);
  node('p','원문 포함 검사 기록: '+(claim.grounding==='EXACT_QUOTE_PRESENT'?'EXACT_QUOTE_PRESENT · 내용의 타당성 판정과 별도':'미확인'),article);
  node('p','종·모델: '+limited(claim.organism,1000)+' · 근거 유형: '+limited(claim.evidence_type,1000)+' · 관계: '+limited(claim.relation,1000),article);
  node('p','해석 이유: '+limited(claim.reasoning),article);node('p','제한점: '+limited(claim.limitations),article);
 }
 if(!claims.length)node('p','표시할 검토 제안이 없습니다. LLM을 실행하지 않았거나 근거가 부족할 수 있습니다.',panel);
 if(claims.length>100)node('p','최초 100개 제안만 표시합니다. 전체는 근거 JSON을 확인하세요.',panel);
 node('h3','프로그램으로 실행한 결정적 검사',panel);node('p','수집·식별자·발췌문 기록 검사이며 과학적 효능·실험 검증과 구분합니다.',panel);
 node('pre',record(d.checks)?jsonText(d.checks):'결정적 검사 기록 미기록 · 통과로 간주하지 않습니다.',panel);
 if(d.run_verification_links!==null&&d.run_verification_links!==undefined){node('h3','별도 실행 검증 연결 기록 · 문헌 검토와 구분',panel);node('p','입력·참조 등 실제 실행 검증의 범위와 영수증 연결 정보입니다. LLM 문헌 검토가 포즈·순위·실험 검증을 승인하지 않습니다.',panel);node('pre',jsonText(d.run_verification_links),panel)}
 node('h3','사람 검토·실제 실험이 필요한 부분',panel);
 const checklist=Array.isArray(d.human_review_checklist)?d.human_review_checklist.filter(item=>typeof item==='string'&&item.trim()):typeof d.human_review_checklist==='string'&&d.human_review_checklist.trim()?[d.human_review_checklist]:[];
 if(checklist.length){const ul=node('ul','',panel);for(const item of checklist.slice(0,100))node('li',limited(item,6000),ul)}else node('p','구체적 검토 항목 미기록 · 사람 검토는 여전히 필요합니다.',panel);
 const details=node('details','',panel);node('summary','전체 기록과 판단 근거',details);node('pre',jsonText(d),details);
}
async function list(){const out=await request('/api/reports');if(!record(out)||!Array.isArray(out.reports))throw Error('저장 목록 응답 형식 미확인');$('list').replaceChildren();let count=0;for(const row of rows(out.reports)){if(typeof row.id!=='string'||!/^[0-9a-f]{32}$/.test(row.id))continue;const id=row.id,b=node('button',id,$('list'));b.onclick=()=>operation(async()=>showReport(await request('/api/report?id='+encodeURIComponent(id)),id));count++}if(!count)node('p','표시할 저장 문헌 근거가 없습니다.',$('list'))}
$('reports').onclick=()=>operation(list);
$('collect').onsubmit=e=>{e.preventDefault();return operation(async()=>{try{if($('review').checked&&!$('model').value)throw Error('ChatGPT 연결 후 모델을 조회하고 사용할 모델을 선택해 주세요.');if($('review').checked&&!reasoningEfforts.includes($('reasoning').value))throw Error('목록에서 추론 강도를 선택해 주세요.');const out=await request('/api/collect',{query:$('query').value.trim(),taxon:Number($('taxon').value),kind:$('kind').value,action:$('action').value,max_papers:Number($('count').value),review:$('review').checked,model:$('model').value||null,reasoning_effort:$('reasoning').value});if(!record(out)||typeof out.id!=='string'||!/^[0-9a-f]{32}$/.test(out.id))throw Error('생성된 근거 기록 식별자를 확인할 수 없습니다.');showReport(await request('/api/report?id='+encodeURIComponent(out.id)),out.id);await list()}finally{await refreshUsage()}})};
Promise.all([status(),list()]).catch(e=>$('message').textContent=text(e?.message,'연결 상태와 저장 목록을 확인할 수 없습니다.'));
</script></html>'''
