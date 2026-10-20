// This view reports collection receipts and proposed interpretation separately.
// It never promotes LLM prose, quote matching or a source hash to scientific validation.
const SCHEMA='skinpipe-literature-evidence-v1';
const record=value=>value!==null&&typeof value==='object'&&!Array.isArray(value);
const esc=value=>String(value??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const text=(value,fallback='미기록')=>typeof value==='string'&&value.trim()?value:fallback;
const bounded=(value,maximum=12000)=>{const s=text(value);return s.length>maximum?s.slice(0,maximum)+'\n[표시 한도 초과 · 이후 내용은 원본 문헌 근거 기록을 확인하세요.]':s};
const list=value=>Array.isArray(value)?value.filter(record):[];
const statuses={
 HOLD_USAGE_LIMIT:'이 앱의 ChatGPT 사용량 한도 도달 · AI 검토 보류',
 HOLD_RATE_LIMIT:'일시적 요청 속도 제한 · AI 검토 보류',
 NOT_RUN:'문헌 수집 미실행 · 고정된 수집 정책 또는 기록 없음',
 HOLD:'문헌 수집·검토 보류',
 HOLD_UNSUPPORTED_TARGET:'지원 범위 밖 표적 · 문헌 수집 보류',
 HOLD_NO_SOURCES:'제한된 검색에서 문헌을 확보하지 못해 보류',
 HOLD_UNAVAILABLE:'수집 또는 검토 연결 보류',HOLD_REVIEW_UNAVAILABLE:'문헌 검토 연결 실패 · 제안 보류',HOLD_PARTIAL_RETRIEVAL:'문헌 수집 일부 실패 · 검토 보류',
 COLLECTED_HUMAN_REVIEW_REQUIRED:'문헌 수집됨 · 사람 검토 필요',
 HOLD_INVALID_REVIEW:'제안 형식·인용 검사가 충족되지 않아 보류',
 REVIEW_PROPOSED_HUMAN_REVIEW_REQUIRED:'문헌 해석 준비됨 · 원문 검토 필요'
};
const known=(map,key)=>typeof key==='string'&&Object.hasOwn(map,key)?map[key]:null;
const coverage={ABSTRACT_ONLY:'초록만 확보',ABSTRACT_EXCERPT:'초록 일부만 확보',
 OA_FULL_TEXT:'공개 원문 확보',OA_FULL_TEXT_EXCERPT:'공개 원문 일부만 확보',METADATA_ONLY:'메타데이터만 확보 · 본문 없음'};
function literatureURL(value){
 if(typeof value!=='string'||value.length>4096||/[\u0000-\u0020\u007f\\]/.test(value))return null;
 try{
  const url=new URL(value);
  if(url.protocol!=='https:'||url.username||url.password||url.port||url.hash)return null;
  const host=url.hostname.toLowerCase(),path=url.pathname;
  if(host==='doi.org'&&/^\/10\.\d{4,9}\/[^\s]+$/i.test(path))return url.href;
  if(host==='pubmed.ncbi.nlm.nih.gov'&&/^\/\d+\/$/.test(path)&&!url.search)return url.href;
  if(host==='pmc.ncbi.nlm.nih.gov'&&/^\/articles\/PMC\d+\/?$/.test(path)&&!url.search)return url.href;
  if(host==='www.ncbi.nlm.nih.gov'&&/^\/pmc\/articles\/PMC\d+\/?$/.test(path)&&!url.search)return url.href;
  if(['europepmc.org','www.europepmc.org'].includes(host)&&/^\/articles?\/(?:MED\/\d+|PMC\/PMC\d+)\/?$/i.test(path)&&!url.search)return url.href;
  if(host==='www.ebi.ac.uk'&&/^\/europepmc\/webservices\/rest\/(?:search|PMC\d+\/fullTextXML)$/.test(path))return url.href;
 }catch{}
 return null;
}
const citationLink=(label,value)=>{const url=literatureURL(value);return url?`<a class="inline-link" href="${esc(url)}" target="_blank" rel="noopener noreferrer">${esc(label)} ↗</a>`:`<span>${esc(label)} · 허용된 문헌 주소 미기록</span>`};
function checkValue(value){
 if(value===null||value===undefined)return '미기록';
 if(typeof value==='boolean')return value?'true · 해당 출처·원문 대조 기록':'false · 해당 출처·원문 대조 기록';
 if(typeof value==='number')return Number.isFinite(value)?String(value):'유효한 수치 미기록';
 if(typeof value==='string')return bounded(value,3000);
 if(record(value)||Array.isArray(value)){try{return bounded(JSON.stringify(value,null,2),6000)}catch{return '검사 기록 형식 미확인'}}
 return '검사 기록 형식 미확인';
}
const requiredReview='<p><strong>검토 필요</strong> · 논문의 방법·대상·결과와 인용문을 읽고 해석을 검토하세요.</p><p>실험 검증: <code>NOT_EVALUATED</code></p>';
function sourceMarkup(source){
 const id=text(source.source_id),url=literatureURL(source.url);
 const pmid=typeof source.pmid==='number'&&Number.isSafeInteger(source.pmid)?String(source.pmid):text(source.pmid,'');
 const doi=text(source.doi,'');
 const publicationTypes=Array.isArray(source.publication_types)?source.publication_types.filter(x=>typeof x==='string').slice(0,30).join(' · '):'';
 const modelTerms=Array.isArray(source.indexed_model_terms)?source.indexed_model_terms.filter(x=>typeof x==='string').slice(0,30).join(' · '):'';
 const identifiers=[/^\d+$/.test(pmid)?citationLink('PMID '+pmid,'https://pubmed.ncbi.nlm.nih.gov/'+pmid+'/'):pmid?'PMID '+esc(pmid)+' · 식별자 형식 미확인':'PMID 미기록',
  /^10\.\d{4,9}\/\S+$/i.test(doi)?citationLink('DOI '+bounded(doi,500),'https://doi.org/'+doi):doi?'DOI '+esc(bounded(doi,500))+' · 식별자 형식 미확인':'DOI 미기록'];
 const retraction=source.retracted===true?'철회 표시 있음':source.retracted===false?'수집 기록상 철회 표시 없음':'철회 여부 미확인';
 const sha=typeof source.text_sha256==='string'&&/^[a-f0-9]{64}$/.test(source.text_sha256)?source.text_sha256:'유효한 본문 해시 미기록';
 return `<li><p><strong>${esc(bounded(source.title,1500))}</strong><br>출처 ID: <code>${esc(id)}</code> · ${url?citationLink('문헌 원문·초록',source.url):'<span>문헌 링크 미확인 · 외부 이동 비활성</span>'}</p><p>${identifiers.join(' · ')}</p><p>확보 범위: ${esc(known(coverage,source.coverage)||'범위 미확인')} · ${esc(text(source.coverage))}${source.review_text_truncated===true?' · 검토용 본문이 잘려 있음':''}</p><p>${esc(retraction)} · 철회 확인 기록: ${esc(bounded(source.retraction_status,1000))}<br>문헌 분류: ${esc(bounded(source.article_class,1000))} · 모델·시험 범위: ${esc(bounded(source.model_scope,3000))}</p><p>출판 유형: ${esc(bounded(publicationTypes,1500))}<br>색인된 모델 용어: ${esc(bounded(modelTerms,1500))}</p><p>본문 SHA256: <code>${esc(sha)}</code></p></li>`;
}
function claimMarkup(claim,sourceIndex){
 const id=text(claim.source_id,''),matches=sourceIndex.get(id)||[];
 const sourceStatus=matches.length===1?'출처 연결됨':matches.length>1?'출처 ID 중복 · 연결 미확인':'출처 연결 미확인';
 const location=Number.isSafeInteger(claim.quote_start)&&claim.quote_start>=0&&Number.isSafeInteger(claim.quote_end)&&claim.quote_end>=claim.quote_start?` · 기록된 문자 위치 [${claim.quote_start}, ${claim.quote_end})`:'';
 return `<li><p><strong>문헌 해석 · 원문 검토 필요</strong><br>${esc(bounded(claim.claim))}</p><p>출처 ID: <code>${esc(id||'미기록')}</code> · ${sourceStatus}${location}</p><p>첨부 발췌문 (제안에 기록된 원문 문자열)</p><blockquote><pre>${esc(bounded(claim.quote,20000))}</pre></blockquote><p>원문 포함 검사 기록: ${claim.grounding==='EXACT_QUOTE_PRESENT'?'EXACT_QUOTE_PRESENT · 인용문 원문 일치':'미확인'}</p><p>종·모델: ${esc(bounded(claim.organism,1000))} · 근거 유형: ${esc(bounded(claim.evidence_type,1000))} · 관계: ${esc(bounded(claim.relation,1000))}</p><p><strong>해석 이유</strong><br>${esc(bounded(claim.reasoning))}</p><p><strong>제한점</strong><br>${esc(bounded(claim.limitations))}</p></li>`;
}
function quotaTiming(error){
 const messages=[];
 if(record(error)&&typeof error.retry_after_seconds==='number'&&Number.isFinite(error.retry_after_seconds)&&error.retry_after_seconds>=0&&error.retry_after_seconds<=31536000)messages.push('제공된 재시도 대기: '+error.retry_after_seconds+'초');
 if(record(error)&&typeof error.reset_at==='string'&&error.reset_at.length<=120&&/^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(?:\.\d+)?Z$/.test(error.reset_at)&&Number.isFinite(Date.parse(error.reset_at)))messages.push('제공된 한도 갱신 시각: '+error.reset_at);
 return messages;
}
function quotaMarkup(dossier){
 const error=record(dossier.review)&&record(dossier.review.error)?dossier.review.error:null;
 const usage=dossier.status==='HOLD_USAGE_LIMIT'||['usage_exhausted','usage_limit_exceeded','subscription_sharing_usage_limit_exceeded'].includes(error?.code)||error?.usage_exhausted===true;
 const rate=dossier.status==='HOLD_RATE_LIMIT'||['rate_limit_exceeded','rate_limited'].includes(error?.code);
 if(!usage&&!rate)return '';
 const code=typeof error?.code==='string'&&/^[A-Za-z][A-Za-z0-9_]{0,63}$/.test(error.code)?error.code:null;
 return `<aside role="status"><p><strong>${usage?'이 앱의 ChatGPT 사용량 한도 도달 · AI 검토 보류':'일시적 요청 속도 제한 · AI 검토 보류'}</strong>${code?' · '+esc(code):''}</p><p>${usage?'이 앱에서 사용할 수 있는 ChatGPT 사용량 한도에 도달해 AI 검토를 보류했습니다. ':''}수집한 문헌·출처 기록은 보존하며 AI 검토만 보류합니다. 자동 유료 대체 경로로 전환하거나 반복 호출하지 않습니다. 사용 가능해진 뒤 직접 다시 시도하세요.</p>${quotaTiming(error).map(message=>`<p>${esc(message)}</p>`).join('')}</aside>`;
}

function humanChecklist(dossier){
 const value=dossier.human_review_checklist;
 const items=Array.isArray(value)?value.filter(item=>typeof item==='string'&&item.trim()):typeof value==='string'&&value.trim()?[value]:[];
 return items.length?`<details open><summary>사람이 확인할 항목 · ${items.length}개</summary><ul>${items.slice(0,100).map(item=>`<li>${esc(bounded(item,6000))}</li>`).join('')}</ul>${items.length>100?'<p>최초 100개 항목만 표시합니다. 전체는 원본 근거 기록을 확인하세요.</p>':''}</details>`:'';
}
function verificationLinks(value){
 if(value===undefined||value===null)return '';
 return `<details><summary>별도 실행 검증 연결 기록 · 문헌 검토와 구분</summary><p>실행에 연결된 검증 기록의 참조입니다. 이 화면은 연결 정보만 표시합니다. 도킹·작용 방향과 실험 상태는 연결된 원본 기록에서 확인하세요.</p><pre>${esc(checkValue(value))}</pre></details>`;
}

export function literatureMarkup(dossier){
 if(dossier===null||dossier===undefined)return `<section class="platform-panel literature-evidence"><h3>문헌 근거</h3><p>문헌 수집 설정 또는 기록 없음 · 아직 근거를 수집·검토하지 않았습니다.</p>${requiredReview}</section>`;
 if(!record(dossier)||dossier.schema!==SCHEMA)return `<section class="platform-panel literature-evidence"><h3>문헌 근거 · 형식 미확인</h3><p>지원하는 근거 계약을 확인할 수 없어 검증 미완료로 표시합니다. 기록 형식: <code>${esc(record(dossier)?bounded(dossier.schema,200):'미확인')}</code></p>${requiredReview}</section>`;
 const sources=list(dossier.sources),claims=list(dossier.claims),sourceIndex=new Map();
 for(const source of sources){const id=text(source.source_id,'');if(id)sourceIndex.set(id,[...(sourceIndex.get(id)||[]),source])}
 const checks=record(dossier.checks)?Object.entries(dossier.checks):[];
 const status=known(statuses,dossier.status)||'수집·검토 상태 미확인 · 검증 미완료';
 return `<section class="platform-panel literature-evidence"><h3>문헌 근거와 해석 제안</h3><p>${esc(status)}</p>${typeof dossier.reason==='string'&&dossier.reason.trim()?`<p>사유: ${esc(bounded(dossier.reason,6000))}</p>`:''}${requiredReview}${quotaMarkup(dossier)}${humanChecklist(dossier)}${verificationLinks(dossier.run_verification_links)}${dossier.human_review_required!==true?'<p>검토자 확인 상태가 누락됐습니다. 원문을 확인하고 검토 상태를 기록하세요.</p>':''}<details open><summary>수집한 문헌과 확보 범위 · ${sources.length}건</summary>${sources.length?`<ol>${sources.slice(0,100).map(sourceMarkup).join('')}</ol>${sources.length>100?'<p>최초 100개 출처만 표시합니다. 전체는 원본 근거 기록을 확인하세요.</p>':''}`:'<p>문헌 출처 기록이 없습니다. 수집 상태와 오류 내용을 확인하세요.</p>'}</details><details><summary>출처·원문 대조 기록 · ${checks.length}개</summary><p>수집 상태, 출처 식별자와 발췌문을 원문에 대조한 결과입니다.</p>${checks.length?`<dl>${checks.slice(0,100).map(([key,value])=>`<dt>${esc(bounded(key,200))}</dt><dd><pre>${esc(checkValue(value))}</pre></dd>`).join('')}</dl>`:'<p>출처·원문 대조 기록이 없습니다. 처리 기록을 확인하세요.</p>'}</details><details open><summary>문헌 해석 · 검토 전 · ${claims.length}건</summary>${claims.length?`<ol>${claims.slice(0,100).map(claim=>claimMarkup(claim,sourceIndex)).join('')}</ol>${claims.length>100?'<p>최초 100개 제안만 표시합니다. 전체는 원본 근거 기록을 확인하세요.</p>':''}`:'<p>문헌 해석이 저장되지 않았습니다. 처리 상태를 확인하세요.</p>'}</details></section>`;
}
