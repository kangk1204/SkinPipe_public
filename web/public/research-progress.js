// Animation describes an acknowledged state or an in-flight request, never a
// percentage, ETA, or proof of new computation between server observations.
export const RESEARCH_PROGRESS_STALE_MS=20000;
const ACTIVE=new Set(['submitting','checking','queued','preparing','running','cancelling']);
const VIEWS={
 submitting:['running','연구 접수 중','서버의 접수 응답을 기다리고 있습니다. 아직 계산 시작을 확인한 것은 아닙니다.',''],
 checking:['unknown','연구 상태 확인 중','현재 처리 단계와 저장된 결과를 조회하고 있습니다.',''],
 queued:['running','연구 실행 대기','연구를 접수했습니다. 서버에서 계산 시작을 기다리고 있습니다.',''],
 preparing:['running','연구 준비 중','서버가 입력과 계산 자료를 준비하고 있습니다.',''],
 running:['running','연구 진행 중','서버가 연구를 처리하고 있습니다. 확인된 단계와 처리 수를 아래에 표시합니다.',''],
 complete:['complete','연구 완료','아래 평가 결과를 확인하고 Excel로 저장하세요. 미평가·보류 항목도 함께 확인하세요.','✓'],
 completed:['complete','연구 완료','아래 평가 결과를 확인하고 Excel로 저장하세요. 미평가·보류 항목도 함께 확인하세요.','✓'],
 partial:['attention','일부 결과 완료','저장된 결과와 처리하지 못한 단계·이유를 함께 확인하세요.','!'],
 failed:['failed','연구 처리 실패','처리 기록과 오류를 확인하세요. 이미 저장된 결과는 아래에서 검토할 수 있습니다.','!'],
 cancelling:['attention','취소 요청됨 · 자원 정리 중','계산 종료와 작업 자원 정리를 확인하고 있습니다. 종료 확인이 끝날 때까지 기다려 주세요.',''],
 cancelled:['cancelled','연구 취소됨','저장된 결과는 보존됩니다. 작업 자원 정리 상태를 아래에서 확인하세요.','■'],
 cancellation_unverified:['attention','자원 정리 확인 필요','취소 요청은 기록됐지만 작업 자원이 반환됐는지 확인하지 못했습니다.','!'],
 cleanup_unverified:['attention','자원 정리 확인 필요','작업 자원 반환을 확인하지 못했습니다. 처리 기록을 확인하세요.','!'],
 paused:['attention','연구 일시정지','저장된 작업은 보존됩니다. 현재 조건과 이어가기 가능 여부를 확인하세요.','Ⅱ'],
 held:['attention','연구 보류','추가로 확인할 조건과 저장된 결과를 아래에서 확인하세요.','!'],
 needs_data:['attention','추가 확인으로 연구 보류','입력 조건과 처리 기록을 확인하세요. 이미 저장된 자료는 보존됩니다.','!'],
 budget_deferred:['attention','처리 시간 한도로 보류','완료한 결과와 보류된 계산 범위를 확인하세요.','!'],
};
export function researchProgressView(state,{connection='current',checkedAt=null,now=Date.now()}={}){
 const stale=ACTIVE.has(state)&&checkedAt!==null&&(!Number.isFinite(checkedAt)||!Number.isFinite(now)||now<checkedAt||now-checkedAt>=RESEARCH_PROGRESS_STALE_MS);
 if(connection!=='current'||stale)return{state:'unknown',tone:'unknown',title:'현재 상태 미확인',detail:'마지막으로 확인한 기록만 남아 있습니다. 결과 갱신으로 현재 상태를 확인하세요.',icon:'?',active:false,busy:false,showProgress:false};
 const [tone,title,detail,icon]=VIEWS[state]||['unknown','연구 상태 확인 필요','결과 갱신으로 현재 처리 상태를 확인하세요.','?'];
 return{state,tone,title,detail,icon,active:ACTIVE.has(state),busy:ACTIVE.has(state),showProgress:['running','preparing'].includes(state)};
}
export function applyResearchProgressState(element,view){
 element.setAttribute('data-state',view.tone);element.setAttribute('data-phase',view.state);
 element.setAttribute('data-active',String(view.active));element.setAttribute('aria-busy',String(view.busy));
}
