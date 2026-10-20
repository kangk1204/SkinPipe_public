const record=value=>value!==null&&typeof value==='object'&&!Array.isArray(value);
const insist=(ok,message)=>{if(!ok)throw Error(message)};
export function validateResearchSelection(snapshots,selection){
 insist(record(selection)&&Array.isArray(snapshots),'명시한 후보 조회 조건이 없습니다.');
 const matches=snapshots.filter(s=>s.snapshot_id===selection.snapshot);
 insist(matches.length===1&&matches[0].status==='verified',
  '명시한 snapshot이 현재 프로젝트에서 검증되지 않았습니다. 다른 snapshot으로 바꾸지 않습니다.');
 const snapshot=matches[0];
 insist(record(snapshot.target_counts)&&Object.hasOwn(snapshot.target_counts,selection.target),
  '명시한 target이 이 snapshot에 없습니다. 다른 target으로 바꾸지 않습니다.');
 if(selection.protocol!==null&&selection.protocol!==undefined){
  const protocols=(snapshot.provenance?.primary_docking_protocols||[]).filter(p=>p.target_id===selection.target&&p.protocol_id===selection.protocol);
  insist(protocols.length===1,'명시한 protocol이 이 snapshot과 target에 없습니다. 조건을 임의로 만들지 않습니다.');
 }
 return snapshot;
}
