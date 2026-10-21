import {reserveDownload,downloadResponse} from './download-stream.js';
import {RESEARCH_PROGRESS_STALE_MS,researchProgressView,applyResearchProgressState} from './research-progress.js';
import {requestErrorText} from './api-request.js';
import {literatureMarkup} from './literature-evidence.js';
import {personalLiteratureLink} from './personal-literature.js';
import {showResearchComplex} from './structure-viewer.js';
import {researchInputRequest,skinContextRequest,generationRequest,screeningRequest,researchIdentityMarkup,boundChoiceMarkup,boundSelectionRequest,researchSourceUnavailable,identityLabel,PRODUCTS,GOALS} from './research-inputs.js';
import {researchSummaryMarkup,researchRunStatus,researchProgressMarkup,researchRunCounts} from './research-results.js';
export {researchInputRequest,skinContextRequest,researchIdentityMarkup} from './research-inputs.js';
export {seedAssessmentMarkup} from './research-results.js';
const esc=x=>String(x??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const post=(api,path,body)=>api(path,{method:'POST',body:JSON.stringify(body)});
export function syncResearchRunRoute(runId,location=globalThis.location,history=globalThis.history){
 if(location?.hash?.split('?')[0]!=='#discovery'||typeof history?.replaceState!=='function')return false;
 if(runId!==null&&(typeof runId!=='string'||!runId||runId.length>200))return false;
 const hash=runId===null?'#discovery':'#discovery?run='+encodeURIComponent(runId);
 if(location.hash===hash)return false;
 history.replaceState(history.state,'',hash);return true;
}
// Keep request acknowledgements only until current server state proves the action.
// A missing response or unresolved cancellation never proves a successful resume.
export function researchControlNotice(request,run){
 const notice=typeof request?.notice==='string'?request.notice:'';
 if(!notice||request?.pending||!run||typeof run!=='object')return notice;
 const cancellation=run.cancellation||{},cancelState=cancellation.state;
 const cancellationObserved=cancellation.requested===true||['stopping','stopped','cleanup_unverified'].includes(cancelState)||['cancelling','cancelled'].includes(run.status)||['cancelling','cancelled'].includes(run.pipeline_status)||run.lifecycle_state==='cancelled';
 if(request.action==='cancel')return cancellationObserved?'':notice;
 if(request.action!=='resume'||cancellationObserved||cancellation.requested!==false||cancelState!==null&&cancelState!==undefined)return notice;
 const verifiedState=['queued','running','preparing','complete','partial','paused','failed'];
 return run.lifecycle_state==='active'||verifiedState.includes(run.status)||verifiedState.includes(run.pipeline_status)?'':notice;
}
export function researchToxicityPolicyMarkup(capabilities){
 const model=capabilities?.additional_toxicity,available=model?.state==='available'&&model?.provider==='admet_ai'&&model?.model_version==='2.0.1'&&model?.policy==='missing_exact_db_predictions_only'&&Number.isInteger(model.max_candidates)&&model.max_candidates>0&&model.max_candidates<=200&&Number.isInteger(model.budget_seconds)&&model.budget_seconds>0&&model.budget_seconds<=600;
 const limits=available?[Number.isInteger(model.max_candidates)&&model.max_candidates>0?'추가 모델 최대 '+model.max_candidates+'개':null,Number.isInteger(model.budget_seconds)&&model.budget_seconds>0?'추가 모델 한도 '+model.budget_seconds+'초':null].filter(Boolean).join(' · '):'';
 const status=available?'추가 독성 모델: 등록됨 · DB 예측이 없는 정확 구조에 자동 적용':model?.state==='unconfigured'?'추가 독성 모델: 미설정 · 해당 예측은 미평가로 남습니다.':model?.state==='held'?'추가 독성 모델: 설정 확인 필요 · 해당 예측을 보류합니다.':'추가 독성 모델: 등록 상태 미확인 · 실행 결과에서 실제 처리 범위를 확인하세요.';
 return`<p><strong>${esc(status)}</strong></p><p>기존 DB의 정확 구조 예측을 먼저 사용합니다.${available?' 없는 구조는 등록된 CPU 모델로 추가 계산합니다.':''} 지원하지 않는 구조·계산 실패·한도 초과는 미평가로 보존합니다.</p>${available?`<p>${esc(model.provider==='admet_ai'?'ADMET-AI':model.provider||'모델 이름 미기록')} · 버전 ${esc(model.model_version||'미기록')}${limits?' · '+esc(limits):''}</p><p class="research-muted">등록된 로컬 CPU 모델을 사용합니다. 화합물별 계산 결과는 연구 결과에서 확인하세요.</p>`:''}<p class="research-muted">이 설정은 새로 제출하는 연구에 적용합니다. 과거 연구의 미평가 기록을 자동으로 다시 계산하지 않습니다. 입력 물질만 평가하거나 새 구조를 생성하는 경우에도 같은 정책을 적용합니다. 전체 연구의 최대 처리 시간은 아래에서 조정할 수 있습니다.</p>`;
}
const ASSISTANCE_TITLES={input_guidance:'입력 확인 안내',execution_explanation:'실행 전 계산 설명',result_explanation:'결과 읽는 순서',followup:'다음 실험 검토'};
export function researchAssistanceMarkup(response={}){
 const status=response.status||'DETERMINISTIC',ai=status==='AI_FACTS_ORGANIZED';
 const note=ai?'ChatGPT가 저장한 사실을 정리했습니다. 입력·수치·순위·완료 상태는 변경하지 않습니다.':status==='DETERMINISTIC_AI_UNAVAILABLE'?'ChatGPT 응답을 확인하지 못해 기본 안내를 보여줍니다. 계산은 계속 사용할 수 있습니다.':'저장된 근거와 정해진 규칙에 따른 기본 안내입니다. ChatGPT 연결 없이 사용할 수 있습니다.';
 const summary=typeof response.summary==='string'?response.summary:'입력과 조건을 확인한 뒤, 후보의 선정 이유·경보·미평가 항목을 함께 검토해 주세요.';
 const transmission=response.research_data_transmission==='COMPLETED'?'허용한 연구 요약을 ChatGPT에 전달했습니다.':response.research_data_transmission==='ATTEMPTED_OUTCOME_UNKNOWN'?'허용한 연구 요약의 전달을 시도했습니다. 응답 완료 여부는 확인되지 않았습니다.':response.research_data_transmission==='REQUEST_OUTCOME_UNKNOWN'?'설명 요청의 응답을 받지 못해 ChatGPT에 자료가 전달됐는지 확인할 수 없습니다. 자동으로 다시 보내지 않습니다.':'이 안내에서 비공개 연구 자료를 ChatGPT에 전달하지 않았습니다.';
 const sections=(Array.isArray(response.explanation?.sections)?response.explanation.sections:[]).slice(0,8),actions=(Array.isArray(response.explanation?.next_actions)?response.explanation.next_actions:[]).filter(action=>typeof action?.text==='string').slice(0,10);
 const sectionsMarkup=sections.map(section=>`<div class="research-assistance-section"><h4>${esc(typeof section?.title==='string'?section.title:'확인할 근거')}</h4>${Array.isArray(section?.facts)?'<ul>'+section.facts.filter(fact=>typeof fact?.text==='string').slice(0,20).map(fact=>'<li>'+esc(fact.text)+'</li>').join('')+'</ul>':typeof section?.text==='string'?'<p>'+esc(section.text)+'</p>':''}</div>`).join('');
 return`<section class="research-assistance"><h3>${esc(response.task_title||ASSISTANCE_TITLES[response.task]||'연구 안내')}</h3><p class="research-assistance-status">${esc(note)}</p><p class="research-assistance-summary">${esc(summary)}</p>${sectionsMarkup}${actions.length?'<h4>다음 행동</h4><ul>'+actions.map(action=>'<li>'+esc(action.text)+'</li>').join('')+'</ul>':''}<p class="research-muted">${esc(transmission)}</p><details><summary>설명에 사용한 근거와 처리 기록</summary><pre>${esc(JSON.stringify(response,null,2))}</pre></details></section>`;
}
export function discoveryLiteratureLinkMarkup(run){
 const descriptor=run.input_descriptor||run.target_descriptor||{},input=descriptor.input||run.research_input||run.request?.input||{},kind=descriptor.kind||input.kind||descriptor.request?.kind;
 if(run.status!=='complete'||['compound','bound_complex'].includes(kind))return'';
 const query=input.value||descriptor.request?.query||descriptor.gene||descriptor.components?.[0]?.gene||input.query;
 if(typeof query!=='string'||!query.trim()||query.length>160)return'';
 return`<p><a class="inline-link" href="${esc(personalLiteratureLink(query))}">새 문헌 검토 · 기존 분석 유지 →</a></p><p class="research-muted">본인 ChatGPT 플랜을 사용하는 선택적 검토입니다. 문헌 화면에서 표적·종·작용 방향과 전달 동의를 확인한 뒤 직접 실행합니다.</p>`;
}
function runHeading(run,descriptor){
 const input=descriptor.input||descriptor.original_input||descriptor;
 if((descriptor.kind||input.kind)==='bound_complex'&&(input.format||descriptor.format)==='upload'){
  const filename=descriptor.status==='verified'&&descriptor.upload_receipt?.filename;
  return typeof filename==='string'&&filename?'업로드한 결합 구조 · '+filename:'업로드한 결합 구조';
 }
 const name=identityLabel(descriptor);return name==='입력 확인'?run.run_id:name;
}
function runLiteratureMarkup(evidence){
 if(evidence?.status==='not_requested')return'<p class="research-muted">이번 화합물 탐색에서는 단백질 문헌 검토를 요청하지 않았습니다.</p>';
 const labels={NOT_RUN:'문헌 검토: 미실행',HOLD:'문헌 검토: 추가 확인 필요',HOLD_UNSUPPORTED_TARGET:'문헌 검토: 지원 범위 확인 필요',HOLD_NO_SOURCES:'문헌 검토: 제한된 검색에서 문헌 미확보',HOLD_UNAVAILABLE:'문헌 검토: 연결 확인 필요',HOLD_REVIEW_UNAVAILABLE:'문헌 검토: AI 연결 보류',HOLD_PARTIAL_RETRIEVAL:'문헌 검토: 일부 자료 수집 보류',HOLD_INVALID_REVIEW:'문헌 검토: 제안 형식·인용 확인 필요',HOLD_USAGE_LIMIT:'문헌 검토: ChatGPT 사용량으로 보류',HOLD_RATE_LIMIT:'문헌 검토: 일시적 요청 제한',COLLECTED_HUMAN_REVIEW_REQUIRED:'문헌 수집됨 · 사람 검토 필요',REVIEW_PROPOSED_HUMAN_REVIEW_REQUIRED:'문헌 해석 준비됨 · 원문 검토 필요'};
 const status=typeof evidence?.status==='string'&&Object.hasOwn(labels,evidence.status)?labels[evidence.status]:evidence?'문헌 검토 상태 확인 필요':'문헌 근거 기록 없음';
 return`<p class="research-muted">${esc(status)}</p><details class="research-run-literature"><summary>문헌 근거 보기</summary>${literatureMarkup(evidence)}</details>`;
}
let generation=0;
export function receptorPreparationRequest({target,pdb='',chain='',ligand='',ph=7.4,hyalPocket=false}){
 if(typeof ph!=='number'||!Number.isFinite(ph)||ph<0||ph>14)throw Error('준비 pH는 0–14 범위의 숫자로 입력해 주세요.');
 if(hyalPocket){
  if(target.kind!=='protein'||target.taxon_id!==9606||!['HYAL1','Q12794'].includes(target.query.toUpperCase()))throw Error('이 문헌 포켓은 사람 HYAL1(Q12794)에만 적용합니다.');
  if((pdb&&pdb.toUpperCase()!=='2PE4')||(chain&&chain!=='A')||ligand)throw Error('문헌 HYAL1 포켓은 2PE4·사슬 A와 사용합니다. 다른 구조·기준 리간드 선택을 확인해 주세요.');
  return{pdb_id:'2PE4',chain:'A',assembly_id:'1',...(ph!==7.4?{ph}:{}),explicit_pocket:{
   schema:'source-residue-pocket-v1',chain:'A',citation:'PMCID:PMC2666596',
   residues:[[129,'D'],[131,'E'],[202,'Y'],[216,'N'],[245,'S'],[247,'Y'],[265,'R']].map(([number,amino_acid])=>({number,amino_acid}))}};
 }
 if(!pdb){if(ph!==7.4)throw Error('pH 조건을 변경하려면 PDB 구조를 지정해 주세요.');return undefined}
 return{pdb_id:pdb,...(chain?{chain}:{}),...(ligand?{ligand_code:ligand}:{}),...(ph!==7.4?{ph}:{})};
}
export async function renderDiscovery(main,api,role='researcher',{isCurrent:isActive}={}){
 const epoch=++generation,startRoute=globalThis.location?.hash?.split('?')[0],routeIsCurrent=()=>epoch===generation&&(isActive?isActive():startRoute!=='#discovery'||globalThis.location?.hash?.split('?')[0]==='#discovery');
 const snapshots=(await api('/v2/snapshots?view=summary')).snapshots||[];if(!routeIsCurrent())return;
 let kind='compound',format='name',selectedCid=null,uploadId=null,selection=null,ligandSmiles='',resolution=null,resolutionKey=null,imageUrl=null,activeRun=null,pollTimer=null,progressTimer=null,lastLiveRun=null,lastCheckedAt=null,resultRequest=0,inputGeneration=0,resolveRequest=0,resolveButtonOwner=0,submissionInFlight=false,submissionSequence=0,navigationIntent=0,reviewingControl=null;
 let exportView=0,exportNotice=0;
 const selectedRuns=new Set(),controlRequests=new Map(),exportLinks=new Set();
 main.innerHTML=`<div class="page-heading"><div><h1>새 피부 화합물 연구</h1><p>무엇을 알고 있는지부터 시작해, 입력을 확인하고 연구 후보를 살펴봅니다.</p></div></div>
 <ol class="research-steps" aria-label="연구 진행 단계"><li data-research-step="1" aria-current="step">1 입력</li><li data-research-step="2">2 입력 확인</li><li data-research-step="3">3 연구 실행</li><li data-research-step="4">4 진행·결과</li></ol>
 <div class="platform-errors research-errors" role="alert" aria-live="polite"></div>
 <p id="research-export-status" class="research-export-status" role="status" aria-live="polite" aria-atomic="true"></p>
 <p id="research-submit-status" class="research-progress-panel" role="status" aria-live="polite" aria-atomic="true"></p>
 <section id="research-step-1" class="platform-panel"><h2>무엇에서 시작할까요?</h2><form id="research-input-form">
 <fieldset class="research-kind-cards"><legend class="sr-only">입력 종류</legend><label class="research-kind-card"><input type="radio" name="research-kind" value="compound" checked><strong>화합물을 알고 있어요</strong><span>이름으로 찾거나, SMILES 분자 구조를 입력합니다.</span></label><label class="research-kind-card"><input type="radio" name="research-kind" value="protein"><strong>단백질을 알고 있어요</strong><span>피부 연구에서 살펴볼 단백질 이름을 입력합니다.</span></label><label class="research-kind-card"><input type="radio" name="research-kind" value="bound_complex"><strong>결합 구조가 있어요</strong><span>단백질과 화합물이 함께 있는 PDB ID 또는 구조 파일을 사용합니다.</span></label></fieldset>
 <div class="platform-toolbar"><label id="research-format-label">입력 방식<select id="research-format"><option value="name">화합물 이름</option><option value="smiles">SMILES 분자 구조 (선택)</option></select></label><label id="research-value-label"><span id="research-value-title">화합물 이름</span><input id="research-value" required placeholder="예: niacinamide, caffeine" autocomplete="off"></label><label id="research-upload-label" hidden>구조 파일 (PDB/mmCIF)<input id="research-upload" type="file" accept=".pdb,.cif,.mmcif" disabled></label></div><p id="research-input-help">이름을 PubChem에서 확인합니다. 여러 분자를 포함하는 물질군 이름은 정확한 종류를 선택해야 합니다.</p>
 <details id="research-protein-options" hidden><summary>다른 종·작용 방향을 지정하려면</summary><p>기본은 사람 단백질이며, 작용 방향은 아직 정하지 않습니다.</p><div class="platform-toolbar"><label>종 (NCBI taxon ID)<input id="research-taxon" type="number" min="1" value="9606"></label><label>살펴볼 작용<select id="research-action"><option value="any">아직 정하지 않음</option><option value="inhibit">억제 후보</option><option value="activate">활성화 후보</option></select></label></div></details>
 <p class="research-muted">접수 전 입력은 이 화면에서만 유지됩니다. 새로고침하거나 다른 메뉴로 이동하면 다시 입력해야 합니다. 접수한 연구는 실행 기록에 보존되며, 화면을 이동해도 접수된 연구는 유지됩니다.</p><button id="research-resolve" class="button primary" type="submit">입력 확인하기 →</button><button id="research-input-help-button" type="button" class="button">입력 방법 도움말</button><div id="research-input-guidance" aria-live="polite"></div></form></section>
 <section id="research-step-2" class="platform-panel" hidden><div id="research-identity" aria-live="polite"></div><details class="research-assistance-controls"><summary>확인 결과 설명 듣기 (선택)</summary><h3>확인 카드가 어려우면</h3><p>기본 안내로 다음 행동을 확인할 수 있습니다. ChatGPT는 입력을 대신 확정하지 않습니다.</p><label><input type="checkbox" id="research-identity-ai">확인한 이름·식별 결과 요약을 ChatGPT 설명에 사용</label><p id="research-identity-scope" class="research-muted">ChatGPT로 전달할 자료: 없음. 원본 구조 파일·SMILES·좌표는 설명에 보내지 않습니다.</p><button type="button" id="research-explain-input" class="button">확인 결과와 다음 행동 안내</button><div id="research-identity-assistance" aria-live="polite"></div></details><div class="form-actions"><button id="research-back-input" class="button">← 입력 수정</button><button id="research-confirm" class="button primary" disabled>이 입력으로 연구하기 →</button></div></section>
 <section id="research-step-3" class="platform-panel" hidden><h2>어떤 후보를 살펴볼까요?</h2><p id="research-run-description"></p><form id="discovery-form">
 <fieldset class="research-mode-cards"><legend id="research-mode-heading">탐색 방법</legend><label id="research-similarity-option"><input type="radio" name="research-mode" value="similarity" checked><strong id="research-similarity-title">비슷한 화합물 찾기</strong><span id="research-similarity-description">출발 화합물과 비슷한 카탈로그 후보를 살펴봅니다.</span></label><label id="research-hop-option"><input type="radio" name="research-mode" value="scaffold_hopping"><strong id="research-hop-title">다른 골격 후보 찾기</strong><span id="research-hop-description">구조가 달라도 일부 특징을 공유하는 후보를 비교합니다.</span></label><label id="research-assessment-option"><input type="radio" name="research-mode" value="assessment"><strong>입력 화합물만 평가</strong><span>후보를 새로 찾지 않고 피부 관련 지표부터 확인합니다.</span></label></fieldset>
 <h3>피부 평가 조건 <small>모르면 비워 두세요</small></h3><div class="platform-toolbar"><label>연구 목적<input id="research-goal" maxlength="500" placeholder="예: 장벽 개선 후보 탐색"></label><label>적용 부위<input id="research-site" maxlength="200" placeholder="예: 얼굴, 두피 · 미정이면 생략"></label><label>주로 살펴볼 위치<select id="research-skin-goal">${Object.entries(GOALS).map(([k,v])=>`<option value="${k}">${v}</option>`).join('')}</select></label><label>제품 형태<select id="research-product"><option value="">아직 정하지 않음</option>${Object.entries(PRODUCTS).map(([k,v])=>`<option value="${k}">${v}</option>`).join('')}</select></label><label>배합 농도 (%)<input id="research-percent" type="number" min="0" max="100" step="any" placeholder="아직 정하지 않음"></label></div><p>침투가 높을수록 좋은 것은 아닙니다. 제품·농도가 미정이면 전신 노출량을 확정하지 않고, 계산 가능한 피부 지표와 경보를 살펴봅니다.</p>
 <section class="research-options"><h3>이번 연구의 계산 범위</h3><p id="research-screening-routine">기본 루틴: 넓게 탐색하고, 필터·평가를 거쳐 대표 후보 최대 200개를 고릅니다. 저장된 예측·도킹은 계산 조건을 대조해 재사용하고, 부족한 계산만 나누어 진행합니다.</p><p id="research-budget-preview" class="research-scope" role="status"></p><div id="research-model-policy" role="status"><p>추가 독성 모델의 등록 상태를 확인하고 있습니다. 기본 입력 확인과 연구는 계속 사용할 수 있습니다.</p></div><details><summary>후보 수·필터·최대 처리 시간 조정</summary><label>후보를 어디서 찾을까요?<select id="discovery-source"><option value="">공통 카탈로그</option>${snapshots.filter(s=>s.status==='verified').flatMap(s=>Object.keys(s.target_counts||{}).map(t=>`<option value="${esc(s.snapshot_id)}|${esc(t)}">${esc(s.provenance?.marker||t)} · ${Number(s.target_counts[t]).toLocaleString('ko-KR')}개 · ${esc(s.snapshot_id)}</option>`)).join('')}</select></label><p>기존 모집단에서는 화학 구조만 재사용합니다. 과거 표적의 활성·도킹 값을 새로운 입력의 근거로 옮기지 않습니다.</p><div class="platform-toolbar"><label><input id="discovery-purchase" type="checkbox" checked>구매 기록 있는 후보</label><label><input id="discovery-critical" type="checkbox" checked>강한 구조 경보 후보 숨김</label><label><input id="discovery-admet" type="checkbox" checked>강한 독성 예측 후보 숨김</label><label id="discovery-opposing-label"><input id="discovery-opposing" type="checkbox" checked>원하는 작용과 반대인 후보 숨김</label><label>대표 후보 목표<input id="research-representatives" type="number" min="1" max="1000" step="1" value="200"></label><label>DB에서 평가할 최대 후보 수<input id="research-catalog-limit" type="number" min="1" max="10000" step="1" value="1000"></label><label>최대 처리 시간 (시간)<input id="discovery-hours" type="number" min="0.1" max="8" step="0.1" value="2"></label></div><p>미평가 후보는 평가값이 비어 있습니다. 필터를 바꾸면 제외 이유와 함께 다시 볼 수 있습니다.</p>
 <div id="research-receptor-options"><h4>단백질 도킹용 구조 (선택)</h4><div class="platform-toolbar"><label>PDB ID<input id="discovery-pdb" placeholder="모르면 비워 두세요"></label><label>단백질 사슬<input id="discovery-chain" placeholder="예: A"></label><label>기준 리간드 코드<input id="discovery-ligand" placeholder="선택"></label><label>준비 pH<input id="discovery-ph" type="number" min="0" max="14" step="0.1" value="7.4"></label></div><label><input id="discovery-hyal-pocket" type="checkbox">사람 HYAL1 · 2PE4 문헌 잔기 포켓</label><p>구조·포켓을 확인할 수 없으면 도킹을 보류하고 가능한 탐색 결과를 보존합니다.</p></div></details></section>
 <details id="research-generation-options" class="research-options"><summary>새 분자 구조도 만들어 비교하기 (선택)</summary><label><input id="research-generation" type="checkbox">카탈로그 후보 외에 새 구조 생성도 요청</label><p>생성 구조는 구매 가능성·합성 가능성·효능이 확인되지 않은 연구 가설입니다. 카탈로그 후보와 구분해 저장합니다.</p><p>작용기 치환 구조는 별도로 추가되어 총 생성 수가 BRICS 재조합 상한보다 많을 수 있습니다.</p><div class="platform-toolbar"><label>BRICS 재조합 상한<input id="research-products" type="number" min="1" max="2000" value="200"></label><label>3차원 비교 상한<input id="research-max3d" type="number" min="1" max="200" value="50"></label><label>분자당 배치 수<input id="research-conformers" type="number" min="1" max="30" value="10"></label><label>합성 난이도 상한 (1–10)<input id="research-max-sa" type="number" min="1" max="10" step="0.1" value="6"></label></div></details>
 <details class="research-options"><summary>ChatGPT 설명과 자료 전달 설정 (선택)</summary><h3>선택적 ChatGPT 설명</h3><label><input id="research-public-literature" type="checkbox" checked>공개 문헌을 이용한 설명 허용</label><label><input id="research-private-ai" type="checkbox">이번 실행의 비공개 입력·결과를 ChatGPT 도우미에게 전달하는 데 동의</label><p>기본 탐색·계산은 ChatGPT 연결 없이도 진행합니다. 비공개 연구 데이터 전달은 기본으로 꺼져 있습니다. 허용하더라도 원본 업로드·분자 구조·좌표는 이 설명에 보내지 않습니다.</p><p id="research-execution-scope" class="research-scope" role="status">ChatGPT로 전달할 비공개 연구 자료: 없음.</p><button type="button" id="research-explain-execution" class="button">실행할 계산을 쉬운 말로 확인</button><div id="research-execution-assistance" aria-live="polite"></div></details>
 <div class="form-actions"><button type="button" id="research-back-confirm" class="button">← 입력 확인</button><button id="discovery-submit" class="button primary" type="submit">새 연구 실행</button></div></form></section>
 <section id="research-step-4" class="platform-panel" hidden><div id="research-live-status" class="research-progress-panel" role="status" aria-live="polite" aria-atomic="false"></div><div id="research-run-controls" aria-live="polite"></div><div id="research-result-actions" class="form-actions"></div><p id="research-result-export-status" class="research-export-status" role="status" aria-live="polite"></p><div id="research-result-summary"></div><div id="research-result-assistance-controls"></div><button id="research-new" class="button">다른 입력으로 새 연구</button></section>
 <section id="discovery-runs" class="platform-panel"><h2>내 연구 실행 기록</h2><p>‘진행·결과 보기’에서 진행 상태와 평가를 열고, 완료된 연구는 Excel로 받을 수 있습니다.</p><div id="research-history-actions" class="form-actions research-history-actions"><button id="discovery-refresh" class="button">상태 갱신</button><button id="research-export-selected" class="button" disabled>선택한 연구를 Excel로 받기</button></div><p id="research-selection-status" role="status">Excel로 묶을 연구를 아래에서 선택해 주세요.</p><div id="discovery-list"></div></section>`;
 const $=id=>main.querySelector('#'+id),error=e=>{if(!isCurrent())return;const target=main.querySelector('.platform-errors');if(!target)return;target.textContent=requestErrorText(e);reveal(target)};
 const isCurrent=()=>routeIsCurrent()&&Boolean(main.querySelector('#discovery-runs'));
 function reveal(target,anchor=target){if(!target)return;target.setAttribute('tabindex','-1');target.focus({preventScroll:true});anchor.scrollIntoView({block:'start',behavior:'instant'})}
 function clearProgressTimer(){if(progressTimer!==null){clearTimeout(progressTimer);progressTimer=null}}
 function stopPoll(){if(pollTimer!==null){clearTimeout(pollTimer);pollTimer=null}clearProgressTimer()}
 function step(n,navigate=true){if(n!==4){stopPoll();reviewingControl=null}for(let i=1;i<=4;i++)$('research-step-'+i).hidden=i!==n;for(const el of main.querySelectorAll('[data-research-step]')){if(Number(el.dataset.researchStep)===n)el.setAttribute('aria-current','step');else el.removeAttribute('aria-current')}main.querySelector('.platform-errors').textContent='';if(navigate){const section=$('research-step-'+n);reveal(section.querySelector('h2')||section,section)}}
 const input=()=>researchInputRequest({kind,format,value:format==='upload'?uploadId:$('research-value').value,taxon:$('research-taxon').value,action:$('research-action').value,pubchemCid:selectedCid,selection,ligandSmiles});
 function clearExportStatus(output){if(output){output.textContent='';delete output.dataset.state;output.setAttribute('aria-busy','false')}}
 function clearExportNotices(){exportView++;exportNotice++;for(const id of['research-export-status','research-result-export-status'])clearExportStatus($(id));for(const link of exportLinks)link.remove?.();exportLinks.clear()}
 function clearIdentity(){inputGeneration++;resolveRequest++;stopPoll();clearExportNotices();activeRun=null;resultRequest++;$('research-generation').checked=false;$('research-private-ai').checked=false;selectedCid=null;selection=null;ligandSmiles='';resolution=null;resolutionKey=null;$('research-confirm').disabled=true;$('research-input-guidance').innerHTML='';$('research-identity-assistance').innerHTML='';$('research-execution-assistance').innerHTML=''}
 function setKind(next){
  kind=next;format=kind==='bound_complex'?'pdb_id':'name';uploadId=null;clearIdentity();$('research-value').value='';$('research-upload').value='';
  $('research-format').innerHTML=kind==='bound_complex'?'<option value="pdb_id">PDB ID</option><option value="upload">PDB/mmCIF 파일</option>':kind==='protein'?'<option value="name">단백질 이름·유전자명·식별자</option>':'<option value="name">화합물 이름</option><option value="smiles">SMILES 분자 구조 (선택)</option>';
  $('research-format').value=format;$('research-format').disabled=kind==='protein';$('research-format-label').hidden=kind==='protein';$('research-protein-options').hidden=kind!=='protein';$('research-upload-label').hidden=true;$('research-upload').disabled=true;$('research-value-label').hidden=false;$('research-value').disabled=false;$('research-value').required=true;
  $('research-value-title').textContent=kind==='protein'?'단백질 이름·유전자명·식별자':kind==='bound_complex'?'PDB ID':'화합물 이름';
  $('research-value').placeholder=kind==='protein'?'예: TYR, tyrosinase, P14679':kind==='bound_complex'?'예: 1ATP':'예: niacinamide, caffeine';
  $('research-input-help').textContent=kind==='protein'?'사람 단백질이 기본입니다. 다른 종이면 아래 조건을 변경해 주세요. 입력을 임의의 대리 단백질로 바꾸지 않습니다.':kind==='bound_complex'?'결합 화합물을 포함한 구조를 사용해 주세요. 다음 단계에서 단백질 사슬과 결합 화합물을 직접 선택합니다.':'이름 검색 시 입력한 이름이 PubChem에 전달됩니다. 여러 분자를 포함하는 물질군 이름은 정확한 종류 또는 SMILES가 필요합니다.';
 }
 for(const radio of main.querySelectorAll('input[name="research-kind"]'))radio.onchange=()=>setKind(radio.value);
 $('research-format').onchange=()=>{format=$('research-format').value;uploadId=null;clearIdentity();const upload=format==='upload';$('research-upload-label').hidden=!upload;$('research-upload').disabled=!upload;$('research-value-label').hidden=upload;$('research-value').disabled=upload;$('research-value').required=!upload;$('research-value-title').textContent=format==='smiles'?'SMILES 분자 구조':format==='pdb_id'?'PDB ID':'화합물 이름';$('research-value').placeholder=format==='smiles'?'예: NC(=O)c1cccnc1':format==='pdb_id'?'예: 1ATP':'예: niacinamide';$('research-input-help').textContent=format==='smiles'?'입체화학·염·전하를 포함해 정확한 분자 구조를 입력해 주세요. 이름 검색을 위해 구조를 PubChem에 보내지 않습니다.':kind==='bound_complex'?'결합 화합물을 포함한 구조가 필요합니다. 다음 단계에서 사용자가 사슬과 화합물을 선택합니다.':'이름 검색 시 입력한 이름이 PubChem에 전달됩니다.'};
 for(const id of['research-value','research-taxon','research-action'])$(id).oninput=clearIdentity;
 $('research-upload').onchange=()=>{uploadId=null;clearIdentity()};
 const assistanceRequests=new Map();
 async function assist(button,holderId,body,valid=()=>true){
  const target=$(holderId),sequence=(assistanceRequests.get(holderId)||0)+1,label=button.textContent;
  assistanceRequests.set(holderId,sequence);button.disabled=true;button.textContent='안내 확인 중…';target.innerHTML='<p role="status">저장한 근거로 안내를 준비합니다. 기본 입력·계산은 계속 사용할 수 있습니다.</p>';
  try{const response=await post(api,'/v2/research-assistance',body);if(isCurrent()&&assistanceRequests.get(holderId)===sequence&&valid())target.innerHTML=researchAssistanceMarkup(response)}
  catch(err){if(isCurrent()&&assistanceRequests.get(holderId)===sequence&&valid())target.innerHTML=researchAssistanceMarkup({task:body.task,status:'DETERMINISTIC_AI_UNAVAILABLE',research_data_transmission:body.permissions?.research_data?'REQUEST_OUTCOME_UNKNOWN':'NOT_ATTEMPTED',summary:body.task==='input_guidance'?'화합물은 정확한 이름 또는 SMILES, 단백질은 이름과 종, 결합 구조는 PDB ID 또는 파일을 입력하세요. 여러 종류가 있으면 직접 선택해야 합니다.':body.task==='execution_explanation'?'입력과 피부 조건을 확인하세요. 단백질이 없는 화합물 연구는 도킹 없이 진행하며, 모든 후보의 피부·독성 자료를 검토합니다.':body.task==='followup'?'선정 이유·경보·미평가 항목을 검토하고, 구매·순도·실험 농도·대조군을 확인한 뒤 실험 계획을 정하세요.':'후보의 구조, 선정 이유, 피부 평가, 독성, 미평가 항목 순서로 확인하세요. 원자료와 Excel도 함께 검토하세요.',message:requestErrorText(err)})}
  finally{button.disabled=false;button.textContent=label}
 }
 $('research-input-help-button').onclick=()=>assist($('research-input-help-button'),'research-input-guidance',{task:'input_guidance',permissions:{research_data:false,public_literature:false}});
 function identityScope(){ $('research-identity-scope').textContent=$('research-identity-ai').checked?'ChatGPT로 전달할 자료: 확인한 이름·식별 상태 요약. 원본 업로드·SMILES·좌표는 전달하지 않습니다.':'ChatGPT로 전달할 자료: 없음. 원본 구조 파일·SMILES·좌표는 설명에 보내지 않습니다.'; }
 $('research-identity-ai').onchange=identityScope;
 $('research-explain-input').onclick=()=>{const id=resolution?.resolution_id;if(!id){error(Error('먼저 입력 확인을 완료해 주세요.'));return}return assist($('research-explain-input'),'research-identity-assistance',{task:'input_guidance',resolution_id:id,permissions:{research_data:$('research-identity-ai').checked,public_literature:false}},()=>resolution?.resolution_id===id)};
 function executionScope(){ $('research-execution-scope').textContent=$('research-private-ai').checked?'ChatGPT로 전달할 자료: 확인한 입력 요약·아래 계산 조건. 이 허용 범위는 실행에 저장되며, 결과 설명은 따로 요청합니다.':'ChatGPT로 전달할 비공개 연구 자료: 없음. 설치 서버의 기본 안내만 사용합니다.'; }
 $('research-private-ai').onchange=executionScope;
 $('research-explain-execution').onclick=()=>{try{const body=runRequest(),id=body.resolution_id;return assist($('research-explain-execution'),'research-execution-assistance',{task:'execution_explanation',resolution_id:id,permissions:body.ai_permissions,execution:Object.fromEntries(['search_mode','skin_context','generation','compute_budget_seconds','catalog_limit','ai_permissions','screening'].map(key=>[key,body[key]]))},()=>resolution?.resolution_id===id)}catch(err){error(err)}};
 async function showIdentity(result){
  if(!isCurrent())return;
  resolution=result.descriptor?{...result.descriptor,resolution_id:result.resolution_id}:result;
  resolutionKey=JSON.stringify(input());
  $('research-identity-ai').checked=false;identityScope();$('research-identity-assistance').innerHTML='';
  $('research-identity').innerHTML=researchIdentityMarkup(resolution);step(2);
  $('research-confirm').disabled=resolution.status!=='verified'||!resolution.resolution_id;
  if($('research-bound-view'))$('research-bound-view').onclick=()=>showResearchComplex({title:identityLabel(resolution),resolutionId:resolution.resolution_id,api});
  for(const button of main.querySelectorAll('[data-compound-cid]'))button.onclick=async()=>{selectedCid=Number(button.dataset.compoundCid);await resolveInput()};
  if(kind==='bound_complex'){
   $('research-bound-choices').innerHTML=researchSourceUnavailable(resolution)?'':resolution.status==='verified'?'<details><summary>구성이나 정확한 분자 구조를 다시 선택 (선택)</summary>'+boundChoiceMarkup(resolution)+'</details>':boundChoiceMarkup(resolution);
   const form=$('research-bound-form');if(form)form.onsubmit=async event=>{
    event.preventDefault();try{
     const ligands=Array.isArray(resolution.choices)?resolution.choices:resolution.inventory?.ligands||[];
     const index=$('research-bound-ligand').value;if(index==='')throw Error('결합 화합물을 선택해 주세요.');
     const ligand=ligands[Number(index)],chains=resolution.inventory?.protein_chains||[];
     const chosen=[...main.querySelectorAll('input[name="bound-chain"]:checked')].map(el=>chains[Number(el.value)]);
     if(chosen.some(c=>String(c.model_id)!==String(ligand.model_id)))throw Error('결합 화합물과 같은 모델의 단백질 사슬을 선택해 주세요.');
     selection=boundSelectionRequest({modelId:ligand.model_id,assemblyId:$('research-bound-assembly').value,proteinChains:chosen.map(c=>c.chain_id),ligandChain:ligand.chain_id,ligandResidue:ligand.residue_id,ligandCode:ligand.code,altloc:$('research-bound-altloc').value.trim()});
     ligandSmiles=$('research-bound-smiles').value.trim();
     await resolveInput();
    }catch(err){$('research-bound-error').textContent=requestErrorText(err)}
   };
  }
  if(imageUrl){URL.revokeObjectURL(imageUrl);imageUrl=null}
  if(resolution.status==='verified'&&resolution.resolution_id&&['compound','bound_complex'].includes(kind)){
   const id=resolution.resolution_id;
   try{const response=resolution.structure_image?{image_data_url:resolution.structure_image}:await api('/v2/research-inputs/'+encodeURIComponent(id)+'/image');if(!isCurrent()||resolution?.resolution_id!==id)return;const data=response.image_data_url;if(typeof data!=='string'||!/^data:image\/png;base64,[A-Za-z0-9+/=]+$/.test(data))throw Error('구조 이미지 미확인');const holder=$('research-identity-image');if(holder)holder.innerHTML=`<img src="${esc(data)}" alt="확인한 화합물의 2차원 구조">`}
   catch{const holder=$('research-identity-image');if(holder)holder.textContent='구조 그림은 현재 표시하지 못했습니다. 식별 근거에서 확인한 구조를 살펴보세요.'}
  }
 }
 async function resolveInput(){
  const button=$('research-resolve'),sequence=++resolveRequest,inputEpoch=inputGeneration;resolveButtonOwner=sequence;button.disabled=true;main.querySelector('.platform-errors').textContent='';
  const current=()=>isCurrent()&&sequence===resolveRequest&&inputEpoch===inputGeneration;
  try{
   if(kind==='bound_complex'&&format==='upload'&&!uploadId){
    const file=$('research-upload').files?.[0];if(!file)throw Error('PDB/mmCIF 구조 파일을 선택해 주세요.');
    if(file.size>50*1024*1024)throw Error('구조 파일은 50 MiB 이하로 선택해 주세요.');
    const extension=file.name.toLowerCase().split('.').pop();if(!['pdb','cif','mmcif'].includes(extension))throw Error('PDB 또는 mmCIF 형식의 파일을 선택해 주세요.');
    button.textContent='구조 파일 확인 중…';
    const text=await file.text();if(!current()||$('research-upload').files?.[0]!==file)return;
    const uploaded=await post(api,'/v2/research-inputs/uploads',{filename:file.name,format:extension==='pdb'?'pdb':'mmcif',text});if(!current()||$('research-upload').files?.[0]!==file)return;
    uploadId=uploaded.upload_id;
    if(!uploadId)throw Error('업로드 식별자를 받지 못했습니다. 파일을 다시 확인해 주세요.');
   }
   const request=input(),key=JSON.stringify(request);button.textContent='입력 정체성 확인 중…';
   const out=await post(api,'/v2/research-inputs/resolve',{input:request});
   if(!current()||JSON.stringify(input())!==key)return;
   await showIdentity(out);
  }catch(err){if(current())error(err)}finally{if(resolveButtonOwner===sequence){button.disabled=false;button.textContent='입력 확인하기 →'}}
 }
 $('research-input-form').onsubmit=event=>{event.preventDefault();return resolveInput()};
 function budgetPreview(){const hours=Number($('discovery-hours').value),limit=Number($('research-catalog-limit').value),assessment=main.querySelector('input[name="research-mode"]:checked')?.value==='assessment';$('research-screening-routine').hidden=assessment;$('research-budget-preview').textContent=assessment?'이번에는 입력 화합물 1개만 평가합니다. DB 후보 탐색·새 구조 생성·도킹은 요청하지 않습니다.':'대표 후보 목표: '+Number($('research-representatives')?.value??200).toLocaleString('ko-KR')+'개 · 전체 연구 최대 처리 시간: '+(Number.isFinite(hours)&&hours>=0.1&&hours<=8?hours+'시간':'입력 확인 필요')+' · DB에서 평가할 최대 후보 수: '+(Number.isInteger(limit)&&limit>=1&&limit<=10000?limit.toLocaleString('ko-KR')+'개':'입력 확인 필요')+'. 실행 시간과 처리 범위는 한도 안에서 달라질 수 있습니다.';}
 $('discovery-hours').oninput=budgetPreview;$('research-catalog-limit').oninput=budgetPreview;if($('research-representatives'))$('research-representatives').oninput=budgetPreview;
 function syncMode(){const mode=main.querySelector('input[name="research-mode"]:checked')?.value||'similarity',disabled=kind==='protein'||mode==='assessment';$('research-generation-options').hidden=disabled;$('research-generation').disabled=disabled;$('research-representatives').disabled=mode==='assessment';budgetPreview()}
 for(const radio of main.querySelectorAll('input[name="research-mode"]'))radio.onchange=syncMode;
 $('research-back-input').onclick=()=>{navigationIntent++;step(1)};
 $('research-back-confirm').onclick=()=>{navigationIntent++;step(2)};
 $('research-new').onclick=()=>{setKind(kind);syncResearchRunRoute(null);step(1)};
 $('research-confirm').onclick=()=>{
  if(!resolution||resolution.status!=='verified'||!resolution.resolution_id||resolutionKey!==JSON.stringify(input())){error(Error('입력이 변경되었습니다. 입력을 다시 확인해 주세요.'));return}
  const protein=kind==='protein';$('research-run-description').textContent=protein?'확인한 단백질에 대한 후보 탐색을 진행합니다. 구조·포켓이 확인되면 도킹을 추가합니다.':kind==='bound_complex'?'선택한 결합 화합물에서 출발해 카탈로그의 유사 후보와 골격 변화 후보를 탐색합니다.':'확인한 화합물을 출발구조로 사용합니다. 활성 측정 대조군으로 간주하지 않으며 단백질을 자동 지정하지 않습니다.';
  $('research-receptor-options').hidden=!protein;$('discovery-opposing-label').hidden=!protein;$('research-generation-options').hidden=protein;
  const bound=kind==='bound_complex';$('research-similarity-option').hidden=bound;$('research-hop-option').hidden=protein;$('research-assessment-option').hidden=protein||bound;
  $('research-mode-heading').textContent=protein?'단백질 표적 기반 탐색':bound?'관측된 결합 구조의 고정 포켓 호핑':'탐색 방법';
  $('research-similarity-title').textContent=protein?'이 단백질의 후보 찾기':'비슷한 화합물 찾기';$('research-similarity-description').textContent=protein?'확인한 단백질의 후보를 찾고, 구조·포켓이 확정되면 도킹합니다.':'출발 화합물과 비슷한 카탈로그 후보를 살펴봅니다.';
  $('research-hop-title').textContent=bound?'같은 단백질 포켓에서 다른 골격 후보 비교':'다른 골격 후보 찾기';$('research-hop-description').textContent=bound?'관측된 결합 화합물로 호핑을 비교하며, 선택한 수용체·좌표계·포켓을 유지해 도킹합니다.':'구조가 달라도 일부 특징을 공유하는 후보를 비교합니다.';
  for(const radio of main.querySelectorAll('input[name="research-mode"]')){radio.checked=radio.value===(bound?'scaffold_hopping':'similarity');radio.disabled=protein?radio.value!=='similarity':bound?radio.value!=='scaffold_hopping':false}
  syncMode();executionScope();budgetPreview();
  step(3);
 };
 if(role==='viewer'){$('discovery-submit').disabled=true;$('discovery-submit').title='연구자 또는 관리자 권한이 필요합니다.'}
 function runRequest(){
   if(!resolution?.resolution_id||resolution.status!=='verified'||resolutionKey!==JSON.stringify(input()))throw Error('먼저 입력을 다시 확인해 주세요.');
   const mode=main.querySelector('input[name="research-mode"]:checked')?.value||'similarity',source=$('discovery-source').value.split('|');
   if(kind==='bound_complex'&&mode!=='scaffold_hopping')throw Error('결합 구조 입력은 관측된 포켓을 유지하는 골격 탐색으로 실행합니다.');
   if(kind==='protein'&&mode!=='similarity')throw Error('단백질 입력은 확인한 표적의 후보 탐색으로 실행합니다.');
   const body={resolution_id:resolution.resolution_id,search_mode:mode,skin_context:skinContextRequest({product:$('research-product').value,appliedPercent:$('research-percent').value,goal:$('research-skin-goal').value,researchGoal:$('research-goal').value,applicationSite:$('research-site').value}),generation:generationRequest({enabled:kind!=='protein'&&mode!=='assessment'&&$('research-generation').checked,products:$('research-products').value,max3d:$('research-max3d').value,conformers:$('research-conformers').value,maxSa:$('research-max-sa').value}),ai_permissions:{public_literature:$('research-public-literature').checked,research_data:$('research-private-ai').checked},filters:{purchasable_only:$('discovery-purchase').checked,hide_critical:$('discovery-critical').checked,hide_opposing:kind==='protein'&&$('discovery-opposing').checked,hide_severe_admet:$('discovery-admet').checked},...(mode!=='assessment'?{screening:screeningRequest($('research-representatives')?.value??200)}:{}),compute_budget_seconds:Number($('discovery-hours').value)*3600,catalog_limit:Number($('research-catalog-limit').value)};
   if(!Number.isInteger(body.catalog_limit)||body.catalog_limit<1||body.catalog_limit>10000)throw Error('DB에서 평가할 후보 수는 1–10,000개의 정수로 입력해 주세요.');
   if(!Number.isFinite(body.compute_budget_seconds)||body.compute_budget_seconds<360||body.compute_budget_seconds>28800)throw Error('최대 처리 시간은 0.1–8시간으로 입력해 주세요.');
   if(source[0]){body.snapshot_id=source[0];body.source_target_id=source[1]}
   if(kind==='protein'){
    const receptor=receptorPreparationRequest({target:{query:input().value,kind:'protein',taxon_id:input().taxon_id},pdb:$('discovery-pdb').value.trim(),chain:$('discovery-chain').value.trim(),ligand:$('discovery-ligand').value.trim(),ph:$('discovery-ph').value.trim()?Number($('discovery-ph').value):NaN,hyalPocket:$('discovery-hyal-pocket').checked});
    if(receptor)body.receptor_request=receptor;
   }
   return body;
 }
 $('discovery-form').onsubmit=async event=>{
  event.preventDefault();const button=$('discovery-submit');if(role==='viewer'||!isCurrent())return;
  const status=$('research-submit-status');
  if(submissionInFlight){status.textContent='연구 접수 여부를 확인하고 있습니다. 같은 요청을 다시 보내지 않습니다. 내 연구 실행 기록에서 접수 상태를 확인할 수 있습니다.';return}
  submissionInFlight=true;const sequence=++submissionSequence,inputEpoch=inputGeneration,intent=navigationIntent,label=button.textContent,form=$('discovery-form');let attempted=false,accepted=false;
  const current=()=>isCurrent()&&sequence===submissionSequence&&inputEpoch===inputGeneration&&intent===navigationIntent;
  status.hidden=false;button.disabled=true;button.textContent='연구 접수 중…';form.setAttribute('aria-busy','true');applyResearchProgressState(status,researchProgressView('submitting'));status.textContent='연구를 서버에 접수하고 있습니다. 화면을 이동해도 접수된 연구는 유지됩니다. 접수한 연구는 내 연구 실행 기록에서 다시 열 수 있습니다.';main.querySelector('.platform-errors').textContent='';reveal(status);
  try{
   const body=runRequest();
   attempted=true;const result=await post(api,'/v2/discovery-runs',body);accepted=true;
   if(!isCurrent())return;await refresh();if(!isCurrent())return;
   status.textContent=current()?'연구가 접수되었습니다. 화면을 나가도 계산은 계속되며 내 연구 실행 기록에서 진행·결과 보기를 누르면 다시 열 수 있습니다.':'이전에 요청한 연구가 접수되었습니다. 수정 중인 입력은 그대로 유지됩니다. 내 연구 실행 기록에서 이전 연구를 확인해 주세요.';
   if(!current())return;
   if(result.run_id){status.hidden=true;await showRun(result.run_id,{navigate:true})}else{step(4);$('research-result-summary').innerHTML='<p>연구가 접수되었습니다. 실행 기록에서 진행 상태와 결과를 확인해 주세요.</p>'}
  }catch(err){if(isCurrent()){status.textContent=attempted&&!accepted?'접수 응답을 확인하지 못했습니다. 같은 실행을 자동으로 다시 보내지 않습니다. 내 연구 실행 기록에서 상태 갱신을 눌러 접수 여부를 먼저 확인해 주세요.':'입력과 조건을 확인해 주세요. 연구는 아직 접수하지 않았습니다.';if(current())error(err)}}finally{submissionInFlight=false;button.disabled=false;button.textContent=label;form.setAttribute('aria-busy','false');status.setAttribute('aria-busy','false');status.setAttribute('data-active','false');status.setAttribute('data-state',accepted?'running':'unknown')}
 };
 async function exportRuns(runIds,format='xlsx',button,ordinaryDownload=false){
  let writer,exportAttempted=false;const label=button?.textContent,view=exportView,notice=++exportNotice;
  const historyRow=button?.closest?.('.research-run'),localOutput=historyRow?.querySelector('[data-run-export-status]')||(button?.closest?.('#research-result-actions')?$('research-result-export-status'):null);
  const attached=node=>!node||main.contains?.(node)!==false;
  // A history export belongs to its row; result/selection exports belong to
  // the current view. A detached button must never receive a ready link.
  const valid=()=>isCurrent()&&attached(button)&&(historyRow?attached(historyRow):view===exportView);
  const ownsNotice=()=>valid()&&view===exportView&&notice===exportNotice;
  const status=(message,state='pending')=>{if(!valid())return;for(const output of new Set([ownsNotice()?$('research-export-status'):null,attached(localOutput)?localOutput:null]))if(output){output.textContent=message;output.dataset.state=state;output.setAttribute('aria-busy',state==='pending'?'true':'false')}};
  try{
   main.querySelector('.platform-errors').innerHTML='';status('내보내기를 준비하고 있습니다.');
   if(button)button.disabled=true;if(!ordinaryDownload&&!api.downloadExport)writer=await reserveDownload(format==='xlsx'?'SkinPipe-results-xlsx.zip':'SkinPipe-results.zip');if(!valid())return;
   exportAttempted=true;status('저장된 연구의 내보내기를 요청하고 있습니다.');
   let out=await post(api,'/v2/exports',{...(runIds.length===1?{run_id:runIds[0]}:{run_ids:runIds}),format,...(format==='xlsx'?{candidate_details:true}:{})});if(!valid())return;const id=out.export_id;
   while(['queued','running'].includes(out.status)){const queued=out.status==='queued';if(button)button.textContent=queued?'내보내기 실행 대기…':'내보내기 준비 중…';status(queued?'내보내기 실행을 기다리고 있습니다. 서버의 앞선 작업이 끝나면 준비를 시작합니다.':'서버에서 내보내기 파일을 준비하고 있습니다.');await new Promise(resolve=>setTimeout(resolve,1500));if(!valid())return;out=await api('/v2/exports/'+id);if(!valid())return}
   if(out.status==='complete'&&!out.result)out=await api('/v2/exports/'+id);
   if(!valid())return;
   if(out.status!=='complete')throw Error('내보내기가 완료되지 않았습니다: '+out.status);
   const filename=out.result.filename;
   if(api.downloadExport){const link=api.downloadExport({exportId:id,filename,anchor:button});if(link&&!historyRow)exportLinks.add(link);status('내보내기 파일이 준비되었습니다. 옆의 저장 링크를 눌러 다운로드하세요. 파일명: '+filename,'ready')}
   else{status('내보내기 파일 데이터를 전달하고 있습니다. 파일명: '+filename);const raw=await api('/v2/exports/'+id+'/download',{rawResponse:true});if(!valid()){await raw.body?.cancel();return}const delivery=await downloadResponse(raw,filename,{writer});writer=null;status(delivery?.mode==='stream'?'파일 데이터 전송이 완료되었습니다. 선택한 저장 위치를 확인하세요. 파일명: '+filename:'브라우저에 다운로드 요청을 전달했습니다. 브라우저 다운로드 목록에서 확인하세요. 파일명: '+filename,'complete');}
   if(ownsNotice())main.querySelector('.platform-errors').innerHTML='';
  }catch(err){
   if(!valid())return;
   if(!ownsNotice()){status('내보내기를 완료하지 못했습니다. '+requestErrorText(err),'failed');return}
   if(err.name==='AbortError'&&!exportAttempted){
    status('직접 저장을 취소했습니다. 내보내기 요청은 아직 접수하지 않았습니다.','cancelled');
    main.querySelector('.platform-errors').innerHTML='<p>직접 저장을 취소했습니다. 내보내기 작업은 아직 시작하지 않았습니다. 원하면 일반 브라우저 다운로드로 받을 수 있습니다 (최대 256 MiB).</p><button type="button" id="research-browser-download" class="button">일반 다운로드로 받기</button>';
    $('research-browser-download').onclick=()=>{if(ownsNotice())return exportRuns(runIds,format,$('research-browser-download'),true)};
    reveal(main.querySelector('.platform-errors'));
   }else if(err.name==='AbortError'){status('다운로드가 중단되었습니다. 접수 여부와 브라우저 다운로드 목록을 확인해 주세요.','failed');error(Error('다운로드가 중단되었습니다. 내보내기 접수 여부를 확인한 뒤 다시 다운로드해 주세요.'))}else{status('내보내기를 완료하지 못했습니다. 오류 안내를 확인해 주세요.','failed');error(err)}
  }finally{if(writer)await writer.abort().catch(()=>{});if(button){button.disabled=false;button.textContent=label}if(isCurrent()&&view===exportView&&notice===exportNotice&&!attached(button))clearExportStatus($('research-export-status'))}
 }
 function openPopulation(summary){sessionStorage.setItem('skinpipe-discovery-result',JSON.stringify({snapshot:summary.snapshot_id,target:summary.target_id,source:'generalized'}));location.hash='research?source=generalized'}
 function runState(run){return run.cancellation?.requested&&run.cancellation.state==='cleanup_unverified'?'cancellation_unverified':run.cancellation?.requested&&run.cancellation.state==='stopping'||run.status==='cancelling'?'cancelling':run.lifecycle_state==='cancelled'?'cancelled':run.status==='needs_data'&&run.progress?.resource_quarantine?.active===true?'cleanup_unverified':run.pipeline_status||run.status}
 function pendingRun(run){return ['queued','running','preparing','cancelling'].includes(runState(run))}
 function runStatus(run){const state=runState(run);return state==='cancellation_unverified'?'취소 요청됨 · 작업 자원 정리 확인 필요':state==='cancelling'?'취소 요청됨 · 계산 종료·자원 정리 중':state==='cancelled'?(run.cancellation?.state==='stopped'&&run.cancellation?.resources_released===true?'취소 완료 · 작업 자원 정리 확인됨':'취소됨 · 작업 자원 정리 확인 필요'):state==='paused'?'일시정지됨 · 작업은 보존됩니다':researchRunStatus(state)}
 const pollInterval=5000,pollDuration=8*60*60*1000+30*60*1000,pollLimit=Math.ceil(pollDuration/pollInterval);
 function resultIsCurrent(runId,sequence){return isCurrent()&&activeRun===runId&&sequence===resultRequest&&!$('research-step-4').hidden}
 function liveUnknown(message){
  clearProgressTimer();const target=$('research-live-status');applyResearchProgressState(target,researchProgressView(null,{connection:'unavailable'}));
  target.textContent=message+(lastLiveRun?' 마지막으로 확인한 상태: '+runStatus(lastLiveRun)+(lastCheckedAt?' ('+lastCheckedAt.toLocaleTimeString('ko-KR',{hour12:false})+')':'')+'.':'');
 }
 function armProgressTimer(runId,sequence){
  clearProgressTimer();progressTimer=setTimeout(()=>{progressTimer=null;if(resultIsCurrent(runId,sequence))liveUnknown('최근 20초 동안 새 상태 응답을 확인하지 못했습니다. 현재 계산 상태는 미확인입니다. 결과 갱신으로 확인해 주세요.');},RESEARCH_PROGRESS_STALE_MS);progressTimer?.unref?.();
 }
 function liveChecking(runId,sequence){
  const target=$('research-live-status'),view=researchProgressView('checking');applyResearchProgressState(target,view);target.innerHTML='<div class="research-progress-heading"><span class="research-progress-icon" aria-hidden="true"></span><div><h2>'+esc(view.title)+'</h2><p>'+esc(view.detail)+'</p></div></div>';armProgressTimer(runId,sequence);
 }
 function liveStatus(run,message,{connection='current'}={}){
  clearProgressTimer();lastLiveRun=run;
  if(connection!=='current'){liveUnknown('최신 실행 상태를 조회하지 못했습니다. 현재 계산 상태는 미확인입니다. '+message);return}
  const checked=new Date();lastCheckedAt=checked;const target=$('research-live-status'),view=researchProgressView(runState(run),{checkedAt:checked.getTime(),now:checked.getTime()});applyResearchProgressState(target,view);
  target.innerHTML='<div class="research-progress-heading"><span class="research-progress-icon" aria-hidden="true">'+esc(view.icon)+'</span><div><h2>'+esc(view.title)+'</h2><p>'+esc(view.detail)+'</p></div></div><p><strong>'+esc(runStatus(run))+'</strong> · '+esc(message)+'</p>'+(view.showProgress?researchProgressMarkup(run.progress,{state:runState(run)}):'')+'<p class="research-progress-checked" aria-live="off">마지막 상태 확인: <time datetime="'+checked.toISOString()+'">'+checked.toLocaleTimeString('ko-KR',{hour12:false})+'</time> · 상태를 마지막으로 조회한 시각입니다. 화면을 이동해도 접수된 연구는 유지됩니다.</p>';
  if(pendingRun(run))armProgressTimer(activeRun,resultRequest);
 }
 function renderRunControls(runId,sequence,run){
  if(!resultIsCurrent(runId,sequence)||reviewingControl?.runId===runId&&reviewingControl.sequence===sequence)return;
  const holder=$('research-run-controls'),request=controlRequests.get(runId),cancellation=run.cancellation||{},released=cancellation.state==='stopped'&&cancellation.resources_released===true,requested=cancellation.requested===true||['cancelling','cancelled'].includes(runState(run));
  holder.setAttribute('aria-busy',request?.pending?'true':'false');
  if(request?.pending){holder.innerHTML='<p role="status">'+(request.action==='cancel'?'취소 요청을 접수하고 있습니다. 아직 계산 종료·자원 반환이 확인된 것은 아닙니다.':'저장한 조건과 남은 계산 한도를 확인해 연구를 이어가고 있습니다.')+'</p><button class="button" disabled>'+(request.action==='cancel'?'취소 요청 중…':'이어가기 요청 중…')+'</button>';return}
  if(request)request.notice=researchControlNotice(request,run);
  const notice=request?.notice?'<p role="status">'+esc(request.notice)+'</p>':'';
  if(requested){
   holder.innerHTML=notice+(released?'<p><strong>작업 프로세스 종료와 자원 정리가 확인되었습니다.</strong> 이미 저장된 결과와 처리 기록은 보존됩니다. 아래에서 결과를 확인하거나 다른 입력으로 새 연구를 시작할 수 있습니다.</p>':cancellation.state==='stopping'||runState(run)==='cancelling'?'<p role="status">취소 요청은 접수되었습니다. 실행 중인 계산을 종료하고 이 연구의 작업 자원을 정리하는 중입니다. 종료 확인 전에는 취소 완료로 표시하지 않습니다. 저장된 결과는 보존됩니다.</p>':'<p role="status">취소 상태이지만 작업 프로세스 종료·자원 정리의 확인 기록이 없습니다. 자원 반환 완료로 해석하지 마세요. 결과 갱신 후에도 같으면 설치 담당자에게 처리 기록과 함께 확인을 요청하세요. 저장된 결과는 보존됩니다.</p>')+(released&&cancellation.resume_allowed===true?'<button id="research-resume-run" class="button">남은 계산 이어가기…</button>':'');
   if(released&&cancellation.resume_allowed===true)$('research-resume-run').onclick=()=>reviewRunControl(runId,sequence,'resume');return;
  }
  if(runState(run)==='budget_deferred'){
   const published=Boolean(run.snapshot_id||run.result_snapshot_id);
   holder.innerHTML=notice+'<section class="notice"><div><h3>처리 시간 한도로 보류되었습니다</h3><p>이미 저장된 계산·체크포인트와 입력·조건은 보존됩니다. 이 연구는 완료된 상태가 아닙니다. '+(published?'이미 등록된 결과는 아래 버튼에서 그 저장 범위만 확인할 수 있습니다.':'최종 결과 스냅샷이 아직 등록되지 않아 전체 후보 표·Excel·오프라인 ZIP을 제공할 수 없습니다. 체크포인트는 최종 결과가 아닙니다.')+'</p><p>더 진행하려면 기존 기록을 유지한 채 별도 연구를 준비하세요. 기존 연구의 시간 한도를 늘리거나 이어서 실행하지 않습니다.</p><ol><li>아래 버튼으로 입력 단계에 돌아가 같은 입력을 다시 확인합니다. 이전 연구의 입력·조건을 자동 복원하지 않습니다. 파일·구조 선택, 후보 출처, 필터, 대표 목표, 생성·수용체 조건을 다시 확인하고 QC 기준은 유지하세요.</li><li>후보 수·필터·최대 처리 시간 조정에서 최대 처리 시간을 기존보다 길게 직접 정하세요(예: 2 → 4시간, 최대 8시간). 이미 8시간이면 이 화면에서 더 늘릴 수 없습니다.</li><li>확인 후 새 연구 실행을 한 번 누르세요. 저장된 계산은 서버가 동일한 계산 조건을 검증한 경우에만 재사용하며, 재사용이나 새 시간 한도 내 완료를 보장하지 않습니다.</li></ol><button type="button" id="research-budget-new" class="button">입력·조건을 다시 확인해 별도 연구 준비</button></div></section>';
   $('research-budget-new').onclick=()=>{if(resultIsCurrent(runId,sequence))$('research-new').onclick()};return;
  }
  if(runState(run)==='cleanup_unverified'){
   holder.innerHTML=notice+'<section class="notice"><div><h3>작업 자원 정리를 확인하지 못해 연구를 보류했습니다</h3><p>입력·처리 기록과 저장된 계산 파일은 보존됩니다. 도킹·QC 단계가 목표에 도달했더라도, 최종 결과로 등록되지 않은 체크포인트는 연구 완료나 Excel 결과로 표시하지 않습니다.</p><p>결과 갱신 후에도 이 상태가 유지되면 내 연구 실행 기록의 처리 기록과 함께 설치 담당자에게 확인을 요청하세요. 정리 확인 전에는 자동으로 재실행하지 않습니다.</p></div></section>';return;
  }
  holder.innerHTML=notice+(['queued','running','preparing','paused'].includes(runState(run))?'<p class="research-muted">화면 이동은 계산을 멈추지 않습니다. 계산을 중단하려면 아래에서 이 연구의 취소를 요청하세요. 저장된 결과는 보존됩니다.</p><button id="research-cancel-run" class="button">이 연구 분석 취소…</button>':'');
  if(['queued','running','preparing','paused'].includes(runState(run)))$('research-cancel-run').onclick=()=>reviewRunControl(runId,sequence,'cancel');
 }
 function reviewRunControl(runId,sequence,action){
  if(!resultIsCurrent(runId,sequence)||controlRequests.get(runId)?.pending)return;
  stopPoll();liveUnknown('확인창이 열려 자동 상태 조회를 잠시 멈췄습니다. 아래에서 계속 진행하거나 요청을 확인해 주세요.');reviewingControl={runId,sequence,action};const holder=$('research-run-controls');
  holder.innerHTML=action==='cancel'?'<div class="notice"><div><strong>이 연구의 분석을 취소할까요?</strong><p>이 연구의 대기·실행 중 계산을 중단하고 작업 자원을 정리합니다. 이미 저장된 결과·처리 기록과 다른 연구는 유지합니다. 종료·정리 확인에는 시간이 걸릴 수 있습니다.</p><div class="form-actions"><button id="research-control-dismiss" class="button">계속 분석하기</button><button id="research-control-confirm" class="button">이 연구 분석 취소</button></div></div></div>':'<div class="notice"><div><strong>저장한 입력·조건으로 남은 계산을 이어갈까요?</strong><p>서버가 종료 상태·입력·남은 계산 한도를 다시 확인합니다. 이미 완료한 결과는 보존합니다.</p><div class="form-actions"><button id="research-control-dismiss" class="button">취소 상태 유지</button><button id="research-control-confirm" class="button">남은 계산 이어가기</button></div></div></div>';
  $('research-control-dismiss').onclick=()=>showRun(runId);$('research-control-confirm').onclick=()=>controlRun(runId,sequence,action);reveal(holder);
 }
 async function controlRun(runId,sequence,action){
  if(!resultIsCurrent(runId,sequence)||controlRequests.get(runId)?.pending)return;
  reviewingControl=null;const request={action,pending:true,notice:''};controlRequests.set(runId,request);stopPoll();renderRunControls(runId,sequence,{status:'running'});
  try{await post(api,'/v2/discovery-runs/'+encodeURIComponent(runId)+'/'+action,{});request.notice=action==='cancel'?'취소 요청을 접수했습니다. 서버가 계산 종료·자원 정리 상태를 확인합니다.':'이어가기 요청을 접수했습니다. 현재 계산 상태를 확인합니다.'}
  catch{request.notice=action==='cancel'?'취소 요청의 응답을 받지 못해 접수 여부·작업 자원 정리 상태를 확인할 수 없습니다. 취소 완료로 해석하지 마세요. 자동으로 다시 요청하지 않습니다. 현재 상태를 조회하며, 연결이 끊겼으면 결과 갱신으로 먼저 확인하세요.':'이어가기 요청의 응답을 받지 못해 접수 여부를 확인할 수 없습니다. 자동으로 다시 요청하지 않습니다. 결과 갱신으로 현재 상태를 먼저 확인하세요.'}
  finally{request.pending=false}
  if(!isCurrent()||activeRun!==runId||$('research-step-4').hidden)return;
  await showRun(runId);if(isCurrent()&&activeRun===runId)await refresh();
 }
 function monitorRun(runId,sequence,{attempt=0,startedAt=Date.now(),failures=0}={}){
  if(!resultIsCurrent(runId,sequence)||reviewingControl?.runId===runId&&reviewingControl.sequence===sequence)return;
  if(attempt>=pollLimit||Date.now()-startedAt>=pollDuration){liveUnknown('이 화면의 자동 상태 확인 한도 (최대 8시간 30분)에 도달했습니다. 서버의 계산을 취소한 것이 아닙니다. 결과 갱신 또는 실행 기록의 상태 갱신으로 현재 상태를 확인해 주세요.');return}
  const delay=Math.min(30000,pollInterval*2**failures);
  pollTimer=setTimeout(async()=>{
   pollTimer=null;if(!resultIsCurrent(runId,sequence))return;
   try{const run=await api('/v2/discovery-runs/'+encodeURIComponent(runId)+'?view=status');if(!resultIsCurrent(runId,sequence))return;
    if(!pendingRun(run)){await showRun(runId,{pollState:{attempt:attempt+1,startedAt}});if(isCurrent()&&activeRun===runId)await refresh();return}
    liveStatus(run,runState(run)==='cancelling'?'5초마다 계산 종료·작업 자원 정리 상태를 확인합니다.':'5초마다 상태를 확인합니다. 완료되면 결과와 Excel 버튼이 자동으로 표시됩니다.');renderRunControls(runId,sequence,run);monitorRun(runId,sequence,{attempt:attempt+1,startedAt,failures:0});
   }catch{if(resultIsCurrent(runId,sequence)){const count=failures+1;if(count>=3){liveUnknown('상태 연결을 연속 3회 확인하지 못해 자동 조회를 멈췄습니다. 현재 계산 상태는 미확인입니다. 저장한 연구는 보존됩니다. 결과 갱신 또는 실행 기록의 상태 갱신을 눌러 다시 확인해 주세요.')}else{liveUnknown('상태 연결을 확인하지 못했습니다 ('+count+'/3회). '+Math.min(30000,pollInterval*2**count)/1000+'초 뒤 다시 조회합니다. 계산 요청을 다시 보내지는 않습니다. 화면을 이동해도 접수된 연구는 유지됩니다.');monitorRun(runId,sequence,{attempt:attempt+1,startedAt,failures:count})}}}
  },delay);pollTimer?.unref?.();
 }
 async function showRun(runId,{navigate=false,pollState={}}={}){
  if(!isCurrent())return;if(navigate)navigationIntent++;
  stopPoll();clearExportNotices();reviewingControl=null;const sequence=++resultRequest;
  try{
   activeRun=runId;lastLiveRun=null;lastCheckedAt=null;step(4,false);$('research-run-controls').innerHTML='';liveChecking(runId,sequence);$('research-result-actions').innerHTML='<button id="research-result-refresh" class="button">결과 갱신</button>';$('research-result-refresh').onclick=()=>refreshRun(runId);$('research-result-assistance-controls').innerHTML='';$('research-result-summary').innerHTML='<p role="status">저장한 결과와 평가 근거를 확인하고 있습니다…</p>';
   const summary=await api('/v2/discovery-runs/'+encodeURIComponent(runId)+'/summary');if(!isCurrent()||activeRun!==runId||sequence!==resultRequest)return;
   if(!submissionInFlight){$('research-submit-status').textContent='';$('research-submit-status').hidden=true}
   $('research-result-summary').innerHTML=researchSummaryMarkup(summary);
   $('research-result-actions').innerHTML=`<button id="research-result-refresh" class="button">결과 갱신</button>${summary.snapshot_id?'<button id="research-open-population" class="button">전체 후보 표 보기</button><button id="research-result-excel" class="button primary">Excel · 전체 후보·개별 물질 탭</button><button id="research-result-zip" class="button">오프라인 ZIP</button>':''}`;
   $('research-result-assistance-controls').innerHTML=`<details class="research-assistance-controls"><summary>결과 설명·다음 실험 안내 (선택)</summary><h3>결과를 읽고 다음 연구 준비하기</h3><p><a class="inline-link" href="#literature?source=generalized">화합물·단백질 공개 문헌 검토 →</a> · 기존 입력은 가져오지 않습니다. 공개 이름을 직접 입력하고 전달에 동의한 뒤 실행합니다.</p><p>계산 결과의 기본 안내는 언제든 사용할 수 있습니다. ChatGPT 설명은 아래 허용 범위를 확인하고 직접 요청합니다.</p><label><input type="checkbox" id="research-result-ai">저장한 결과 요약을 이번 ChatGPT 설명에 사용</label><p id="research-result-scope" class="research-scope" role="status"></p><div class="form-actions"><button id="research-explain-result" class="button">결과 해석 도움말</button><button id="research-followup" class="button">다음 실험 검토 안내</button></div><div id="research-result-assistance" aria-live="polite"></div></details>`;
   $('research-result-refresh').onclick=()=>refreshRun(runId);
   if(summary.snapshot_id){$('research-open-population').onclick=()=>openPopulation(summary);$('research-result-excel').onclick=()=>exportRuns([runId],'xlsx',$('research-result-excel'));$('research-result-zip').onclick=()=>exportRuns([runId],'zip',$('research-result-zip'))}
   const allowed=summary.ai_permissions?.research_data===true,privateBox=$('research-result-ai');privateBox.checked=false;privateBox.disabled=!allowed;
   const resultScope=()=>{$('research-result-scope').textContent=!allowed?'이 실행은 비공개 자료 전달을 허용하지 않았습니다. 결과·후속 연구의 기본 안내만 사용합니다.':privateBox.checked?'ChatGPT로 전달할 자료: 저장한 입력 요약·처리 상태·집계·후보 선정 이유·경보·미평가 항목. 원본 구조·SMILES·좌표는 전달하지 않습니다.':'ChatGPT로 전달할 자료: 없음. 기본 안내만 사용합니다.'};
   privateBox.onchange=resultScope;resultScope();
   for(const [id,task]of[['research-explain-result','result_explanation'],['research-followup','followup']])$(id).onclick=()=>assist($(id),'research-result-assistance',{run_id:runId,task,permissions:{research_data:allowed&&privateBox.checked,public_literature:summary.ai_permissions?.public_literature===true}},()=>activeRun===runId);
   let live=summary,pending=pendingRun(summary),summaryNotice='',liveConnection='current';
   if(pending)try{
    const run=await api('/v2/discovery-runs/'+encodeURIComponent(runId)+'?view=status');if(!resultIsCurrent(runId,sequence))return;live=run;pending=pendingRun(run);
    if(!pending)try{
     const fresh=await api('/v2/discovery-runs/'+encodeURIComponent(runId)+'/summary');if(!resultIsCurrent(runId,sequence))return;
     if(!pendingRun(fresh)){await showRun(runId,{navigate,pollState});return}
     summaryNotice='실행 상태는 '+researchRunStatus(runState(run))+'이지만 결과 요약에 아직 반영되지 않았습니다. 아래 요약은 이전 기록입니다. 결과 갱신을 눌러 최신 결과와 Excel 준비 여부를 확인해 주세요.';
    }catch{if(!resultIsCurrent(runId,sequence))return;summaryNotice='실행 상태는 '+researchRunStatus(runState(run))+'이지만 결과 요약의 최신 조회에 실패했습니다. 아래 요약은 이전 기록입니다. 결과 갱신을 눌러 최신 결과와 Excel 준비 여부를 다시 확인해 주세요.'}
   }catch{if(!resultIsCurrent(runId,sequence))return;liveConnection='unavailable'}
   if(!resultIsCurrent(runId,sequence))return;
   liveStatus(live,summaryNotice||(runState(live)==='cancelling'?'5초마다 계산 종료·작업 자원 정리 상태를 확인합니다.':pending?'5초마다 상태를 확인합니다. 완료되면 결과와 Excel 버튼이 자동으로 표시됩니다.':'아래에서 평가 범위와 미평가 항목을 확인하세요.'),{connection:liveConnection});
   renderRunControls(runId,sequence,live);
   if(summaryNotice)$('research-result-summary').innerHTML='<div class="notice" role="status"><p>'+esc(summaryNotice)+'</p></div>'+$('research-result-summary').innerHTML;
   if(navigate){syncResearchRunRoute(runId);reveal($('research-result-summary').querySelector('h2')||$('research-step-4'),$('research-step-4'))}
   if(pending)monitorRun(runId,sequence,pollState);
  }catch(err){if(isCurrent()&&activeRun===runId&&sequence===resultRequest){$('research-run-controls').setAttribute('aria-busy','false');$('research-run-controls').innerHTML=controlRequests.get(runId)?.notice?'<p role="status">'+esc(controlRequests.get(runId).notice)+'</p>':'';liveUnknown('결과 요약을 불러오지 못했습니다. 연구가 취소되거나 실패했다는 뜻은 아닙니다. 저장한 연구는 보존됩니다. 결과 갱신을 눌러 이 연구의 요약을 다시 확인해 주세요.');$('research-result-summary').innerHTML='<p role="status">이 연구의 결과 요약을 아직 확인하지 못했습니다. 이전 연구의 결과를 대신 표시하지 않습니다.</p>';error(err);try{const run=await api('/v2/discovery-runs/'+encodeURIComponent(runId)+'?view=status');if(!resultIsCurrent(runId,sequence))return;renderRunControls(runId,sequence,run);if(pendingRun(run))monitorRun(runId,sequence,pollState)}catch{if(resultIsCurrent(runId,sequence))$('research-run-controls').innerHTML+='<p role="status">실행 상태도 조회하지 못했습니다. 결과 갱신으로 이 연구의 취소·자원 정리 상태를 먼저 확인하세요.</p>'}}}
 }
 async function refreshRun(runId){await showRun(runId);if(isCurrent()&&activeRun===runId&&!$('research-step-4').hidden)await refresh()}
 async function refresh(){
  try{
   const result=await api('/v2/discovery-runs?view=status');if(!isCurrent())return;
   $('discovery-list').innerHTML=(result.runs||[]).map((run,index)=>{
    const descriptor=run.input_descriptor||run.target_descriptor||run.research_input||run.request?.input||{},available=Boolean(run.snapshot_id);
    return`<article class="research-run"><label><input type="checkbox" data-select-run="${index}" ${selectedRuns.has(run.run_id)?'checked':''} ${available?'':'disabled'}>Excel에 포함</label><div><h3>${esc(runHeading(run,descriptor))}</h3><p>${esc(runStatus(run))} · ${esc(researchRunCounts(run.counts))}</p>${researchProgressMarkup(run.progress,{state:runState(run)})}${runLiteratureMarkup(run.literature_evidence)}${discoveryLiteratureLinkMarkup(run)}<div class="form-actions"><button data-view-run="${index}" class="button">진행·결과 보기</button>${available?`<button data-export-run="${index}" class="button">Excel 받기</button>`:''}<button data-events="${index}" class="button">처리 기록</button></div><p data-run-export-status class="research-export-status" role="status" aria-live="polite"></p><pre id="run-events-${index}"></pre></div></article>`;
   }).join('')||'<p>아직 등록한 연구가 없습니다. 위에서 첫 입력을 확인해 보세요.</p>';
   for(const checkbox of main.querySelectorAll('[data-select-run]'))checkbox.onchange=()=>{const id=result.runs[Number(checkbox.dataset.selectRun)].run_id;if(checkbox.checked)selectedRuns.add(id);else selectedRuns.delete(id);selectionStatus()};
   for(const button of main.querySelectorAll('[data-view-run]'))button.onclick=async()=>{const label=button.textContent;button.disabled=true;button.textContent='불러오는 중…';try{await showRun(result.runs[Number(button.dataset.viewRun)].run_id,{navigate:true})}finally{button.disabled=false;button.textContent=label}};
   for(const button of main.querySelectorAll('[data-export-run]'))button.onclick=()=>exportRuns([result.runs[Number(button.dataset.exportRun)].run_id],'xlsx',button);
   for(const button of main.querySelectorAll('[data-events]'))button.onclick=async()=>{try{const i=Number(button.dataset.events),output=$('run-events-'+i),out=await api('/v2/discovery-runs/'+result.runs[i].run_id+'/events');if(isCurrent()&&output&&$('run-events-'+i)===output)output.textContent=JSON.stringify(out,null,2)}catch(err){error(err)}};
   selectionStatus();
  }catch(err){if(isCurrent())error(err)}
 }
 function selectionStatus(){$('research-export-selected').disabled=!selectedRuns.size;$('research-selection-status').textContent=selectedRuns.size?'Excel에 포함할 연구 '+selectedRuns.size+'개 선택됨 · 선택한 연구를 Excel로 받기 버튼을 누르세요.':'Excel로 묶을 연구를 아래에서 선택해 주세요.';$('research-history-actions').classList?.toggle('is-selected',selectedRuns.size>0)}
 $('discovery-refresh').onclick=refresh;$('research-export-selected').onclick=()=>exportRuns([...selectedRuns],'xlsx',$('research-export-selected'));
 await refresh();
 try{const capabilities=await api('/v2/research-capabilities');if(isCurrent())$('research-model-policy').innerHTML=researchToxicityPolicyMarkup(capabilities)}catch{if(isCurrent())$('research-model-policy').innerHTML=researchToxicityPolicyMarkup(null)}
 const requestedRun=new URLSearchParams(globalThis.location?.hash?.split('?')[1]||'').get('run');if(requestedRun&&isCurrent())await showRun(requestedRun,{navigate:true});
}
