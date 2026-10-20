// Display/selection metadata only. Numeric sorting does not calibrate affinity.
const finite=v=>typeof v==='number'&&Number.isFinite(v);
const record=v=>v!==null&&typeof v==='object'&&!Array.isArray(v);
const engines=new Set(['vina','unidock']);
const hash=v=>typeof v==='string'&&/^[a-f0-9]{64}$/.test(v);
export function dockingInfo(candidate) {
 const summary=candidate.docking_summary;
 if(record(summary))return {actualEngine:engines.has(summary.actual_engine)?summary.actual_engine:null,
  rawCount:Number.isInteger(summary.raw_output_count)&&summary.raw_output_count>0?summary.raw_output_count:0,
  rawEngines:Array.isArray(summary.raw_engines)?summary.raw_engines.filter(e=>engines.has(e)):[],
  protocol:summary.profile_status==='declared_inputs_bound_to_preparation'&&hash(summary.protocol_digest)?summary.protocol_digest:null};
 const trace=candidate.docking_execution,attempts=record(trace)&&Array.isArray(trace.attempts)?trace.attempts:[];
 return {actualEngine:record(trace)&&trace.state==='completed'&&engines.has(trace.actual_engine)&&trace.engine_identity?.engine===trace.actual_engine?trace.actual_engine:null,
  rawCount:attempts.filter(a=>record(a)&&record(a.raw_output_artifact)).length,
  rawEngines:[...new Set(attempts.filter(a=>record(a)&&record(a.raw_output_artifact)&&engines.has(a.engine)).map(a=>a.engine))],
  protocol:null};
}
export function hasDockEvidence(candidate){return finite(candidate.dock_affinity)||dockingInfo(candidate).rawCount>0}
export const engineName=value=>({vina:'Vina',unidock:'UniDock'}[value]||'미확인');
export function dockingCellValue(candidate,key){const info=dockingInfo(candidate);
 if(key==='docking_engine')return info.actualEngine?engineName(info.actualEngine):null;
 if(key==='docking_protocol')return info.protocol;
 return candidate[key];
}
export function dockingOptions(rows){
 const profiles=new Map();
 for(const row of rows){const d=dockingInfo(row);if(!d.protocol)continue;
  const group=profiles.get(d.protocol)||{id:d.protocol,engine:d.actualEngine,count:0};group.count++;profiles.set(d.protocol,group)}
 return [...profiles.values()].sort((a,b)=>a.id.localeCompare(b.id));
}
export function dockingComparison(rows){const scored=rows.filter(r=>finite(r.dock_affinity));
 const profiles=new Set(),actual=new Set();let unknown=0;
 for(const r of scored){const info=dockingInfo(r);if(info.protocol)profiles.add(info.protocol);else unknown++;if(info.actualEngine)actual.add(info.actualEngine)}
 return {count:scored.length,profiles:profiles.size,engines:actual.size,unknown,
  caution:scored.length>0&&(unknown>0||profiles.size>1||actual.size>1)};
}
export function dockingEntries(manifest){
 const entries=(manifest.poses||[]).map(p=>({...p,id:`pose:${p.pose_index}`,kind:'accepted',renderable:true,
  label:`채택 포즈 ${p.pose_index} · ${finite(p.affinity)?p.affinity.toFixed(3):'미확인'} kcal/mol${p.ranked?' · 결과에 반영':''}`}));
 for(const raw of manifest.raw_outputs||[]){
  const base={original_path:raw.original_path,original_sha256:raw.source_sha256,engine:raw.engine,state:raw.state,
   raw_attempt:{engine:raw.engine,state:raw.state,engine_identity:raw.engine_identity,batch_context:raw.batch_context,
    docking_conditions:raw.docking_conditions,conditions_scope:raw.conditions_scope,score_review_policy:raw.score_review_policy}};
  if(!raw.models?.length)entries.push({...base,id:`raw:${raw.attempt_index}`,kind:'raw',renderable:false,
   label:`보관 원본 · ${engineName(raw.engine)} · 3D 표시 불가`,reason:raw.reason});
  for(const m of raw.models||[]){const score=finite(m.affinity)?`${m.affinity.toFixed(3)} kcal/mol`:'점수 미확인';
   const flag=m.score_flag==='finite_outside_operational_window'?' · 운영 구간 밖':m.score_flag==='nonfinite_affinity'?' · 비유한 점수':m.score_binding_verified!==true?' · 점수 연결 미확인':'';
   entries.push({...base,...m,id:`raw:${raw.attempt_index}:${m.pose_index}`,kind:'raw',source_sha256:raw.source_sha256,
    label:`보관 원본 · ${engineName(raw.engine)} · 모델 ${m.model_serial??m.pose_index}${m.model_serial!==null&&m.model_serial!==m.pose_index?` (순서 ${m.pose_index})`:''} · ${score}${flag}${m.renderable?'':' · 3D 표시 불가'}`})
  }
 }
 return entries;
}
