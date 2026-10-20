// Mol* is loaded from the site's own pinned assets when a structure is opened.
import {dockingEntries} from './docking-data.js';
const $ = s => document.querySelector(s);
const esc = x => String(x ?? '').replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
let library, session, generation = 0;
export function structureWarningText(payload){return (payload.warnings||[]).map(w=>`${w.component||'구조'} · ${w.status||'경고'} · ${w.reason||'원본 기록 확인 필요'}`).join('\n');}

function loadLibrary() {
  if (window.molstar?.Viewer) return Promise.resolve(window.molstar);
  if (library) return library;
  library = new Promise((resolve, reject) => {
    if (!$('#molstar-css')) {
      const css = document.createElement('link');
      css.id = 'molstar-css'; css.rel = 'stylesheet'; css.href = 'vendor/molstar/molstar.css';
      document.head.append(css);
    }
    const script = document.createElement('script');
    script.src = 'vendor/molstar/molstar.js?sha256=7d2fe9bc69577f3f7615f2c460ddefd5a7a47bbb425feaa811d2c4e5a611be00';
    script.onload = () => resolve(window.molstar);
    script.onerror = () => { script.remove(); library = null; reject(Error('Mol* 뷰어를 불러오지 못했습니다. 다시 열어 주세요.')); };
    document.head.append(script);
  });
  return library;
}

function saveAsset(asset) {
  if (!asset) return;
  const name = asset.download_file || (asset.format==='pdb' ? asset.source_file.replace(/\.pdbqt$/, '.pdb') : asset.source_file);
  const url = URL.createObjectURL(new Blob([asset.data], {type:'text/plain;charset=utf-8'}));
  const a = document.createElement('a'); a.href = url; a.download = name; a.click();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

function cleanup() {
  generation++;
  if (!session) return;
  session.controller.abort();
  session.viewer?.dispose();
  const focus = session.focus;
  session = null;
  if (focus?.isConnected) focus.focus();
}

async function addStructure(plugin, asset, role, complex) {
  const raw = await plugin.builders.data.rawData({data:asset.data,label:asset.label});
  const trajectory = await plugin.builders.structure.parseTrajectory(raw, asset.format);
  const model = await plugin.builders.structure.createModel(trajectory);
  // A model, rather than a biological assembly, preserves the docking frame.
  const structure = await plugin.builders.structure.createStructure(model, {name:'model',params:{}});
  if (role === 'ligand' || role === 'reference' || role === 'observed') {
    const all = await plugin.builders.structure.tryCreateComponentStatic(structure, 'all', {label:role==='observed'?'관측된 기준 리간드':role==='reference'?'기준 리간드':'도킹 화합물'});
    await plugin.builders.structure.representation.addRepresentation(all, {type:'ball-and-stick',
      color:'element-symbol',colorParams:{carbonColor:{name:'uniform',params:{value:role==='reference'?0xc58b28:0x23766b}}},
      typeParams:{sizeFactor:0.32,ignoreHydrogens:true}});
  } else {
    const polymer = await plugin.builders.structure.tryCreateComponentStatic(structure, 'polymer', {label:'수용체'});
    if (polymer) await plugin.builders.structure.representation.addRepresentation(polymer, {type:'cartoon',
      color:complex?'uniform':'chain-id',colorParams:complex?{value:0x819b91}:{},
      typeParams:{alpha:complex?0.45:1}});
    const other = await plugin.builders.structure.tryCreateComponentStatic(structure, 'non-standard', {label:'저장된 비단백질 성분'});
    if (other) await plugin.builders.structure.representation.addRepresentation(other, {type:'ball-and-stick',color:'element-symbol'});
    if (!polymer && !other) {
      const all = await plugin.builders.structure.tryCreateComponentStatic(structure, 'all');
      await plugin.builders.structure.representation.addRepresentation(all, {type:'ball-and-stick',color:'element-symbol'});
    }
  }
  return structure.obj?.data;
}

async function display(index) {
  const s = session;
  if (!s || s.busy) return;
  s.busy = true;
  s.receptor = s.ligand = s.receptorStructure = s.ligandStructure = null;
  $('#molstar-view').removeAttribute('data-loaded');
  $('#molstar-view').removeAttribute('data-structure-count');
  $('#structure-source').textContent = '';
  $('#structure-warnings').hidden=true;$('#structure-warnings').textContent='';
  const epoch = generation;
  const current = () => session === s && generation === epoch && $('#structure-dialog').open;
  const status = $('#structure-status');
  status.textContent = '구조를 불러오는 중입니다…';
  status.classList.remove('error');
  $('#structure-loading').hidden = false;
  for (const c of document.querySelectorAll('.structure-controls select,.structure-controls button')) c.disabled = true;
  try {
    const entry = s.entries[index];
    s.currentEntry = entry;
    if(s.methods&&$('#structure-method-record'))$('#structure-method-record').textContent=JSON.stringify({
      ...s.methods.record,selected_entry:entry.kind==='raw'?{kind:'raw',attempt:entry.raw_attempt,
        pose_order:entry.pose_index,model_serial:entry.model_serial,score_binding_verified:entry.score_binding_verified}
        :{kind:'accepted',pose_index:entry.pose_index,protocol:entry.protocol}},null,2);
    if(s.viewer) await s.viewer.plugin.clear();
    if(!current()) return;
    $('#structure-title').textContent = !s.complex&&entry.node&&entry.node!==s.title ? s.title+' · '+entry.node : s.title;
    const assets = await s.load(entry, s.controller.signal);
    if (!current()) return;
    $('#structure-warnings').textContent=structureWarningText(assets);
    $('#structure-warnings').hidden=!(assets.warnings||[]).length;
    if (!current()) return;
    const molstar = await loadLibrary();
    if (!current()) return;
    if (!s.viewer) {
      const viewer = await molstar.Viewer.create('molstar-view', {layoutIsExpanded:false,layoutShowControls:false,
        layoutShowRemoteState:false,layoutShowSequence:false,layoutShowLog:false,layoutShowLeftPanel:false,
        viewportShowExpand:false,viewportShowSelectionMode:false,viewportShowAnimation:false,
        volumeStreamingDisabled:true});
      if (!current()) { viewer.dispose(); return; }
      s.viewer = viewer;
      if (!viewer.plugin.canvas3d) throw Error('이 브라우저에서 3D 표시를 시작하지 못했습니다. WebGL 지원을 확인해 주세요.');
      viewer.plugin.canvas3d.setProps({renderer:{backgroundColor:0xf7faf8}});
    }
    await s.viewer.plugin.clear();
    if (!current()) return;
    s.receptor = assets.receptor;
    s.ligand = assets.ligand;
    s.receptorStructure = await addStructure(s.viewer.plugin, assets.receptor, 'receptor', s.complex);
    if (!current()) return;
    s.ligandStructure = assets.ligand ? await addStructure(s.viewer.plugin, assets.ligand, s.observed?'observed':'ligand', true) : null;
    if (!current()) return;
    if(assets.reference)await addStructure(s.viewer.plugin,assets.reference,'reference',true);
    if (!current()) return;
    s.viewer.handleResize();
    if (s.ligandStructure) s.viewer.plugin.managers.camera.focusSphere(s.ligandStructure.boundary.sphere,
      {minRadius:12,extraRadius:8,durationMs:0});
    else s.viewer.plugin.managers.camera.reset(undefined,0);
    status.textContent = s.complex ? `${entry.label}${entry.reason?' · '+entry.reason:''} · 저장 기록의 식별 범위는 계산 조건에서 확인하세요.${assets.ligand?.kind==='raw_docking_projection'?` 가상 점 ${assets.ligand.omitted_glue_points}개 제외 · 표시 결합은 뷰어 추정입니다.`:''}` : entry.label;
    $('#structure-source').textContent = [assets.receptor,assets.ligand,assets.reference].filter(Boolean)
      .map(a => `${a.source_file} · SHA-256 ${a.source_sha256.slice(0,16)}…`).join(' / ');
    $('#molstar-view').dataset.loaded = 'true';
    $('#molstar-view').dataset.selection = String(index);
    $('#molstar-view').dataset.structureCount = String([assets.receptor,assets.ligand,assets.reference].filter(Boolean).length);
  } catch (e) {
    if (current()) { status.textContent = e.message; status.classList.add('error'); }
  } finally {
    if (current()) {
      s.busy = false;
      $('#structure-loading').hidden = true;
      for (const c of document.querySelectorAll('.structure-controls select,.structure-controls button')) c.disabled = false;
      if ($('#focus-binding')) $('#focus-binding').disabled = !s.ligandStructure;
      $('#save-receptor').disabled = !s.receptor;
      if ($('#save-ligand')) $('#save-ligand').disabled = !s.ligand;
      if ($('#save-ligand')) $('#save-ligand').textContent = s.observed?'관측 리간드 SDF 저장':s.ligand?.format==='sdf'?'채택 포즈 SDF 저장':'표시용 PDB 저장';
      if ($('#save-raw')) {$('#save-raw').disabled=!s.currentEntry?.original_path;$('#save-raw').textContent=s.currentEntry?.original_path?'원본 PDBQT 저장':'원본 출력 미저장'}
    }
  }
}

function open({title,subtitle,entries,selected=0,complex=false,observed=false,load,original,methods}) {
  cleanup();
  let dialog = $('#structure-dialog');
  if (!dialog) { dialog=document.createElement('dialog'); dialog.id='structure-dialog';
    dialog.setAttribute('aria-labelledby','structure-title'); document.body.append(dialog); }
  session = {title,entries,complex,observed,load,original,methods,controller:new AbortController(),focus:document.activeElement,viewer:null,busy:false};
  dialog.innerHTML = `<div class="structure-layout"><div class="structure-header"><div><small>${esc(subtitle)}</small>
    <h2 id="structure-title">${esc(title)}</h2></div><button class="icon-button" id="close-structure" aria-label="구조 닫기">×</button></div>
    <div class="structure-controls"><label for="structure-choice">${observed?'관측 결합 구조':complex?'도킹 포즈':'표적 구조'}</label>
    <select id="structure-choice">${entries.map((e,i)=>`<option value="${i}" ${i===selected?'selected':''}>${esc(e.label)}</option>`).join('')}</select>
    ${complex?'<button class="button" id="focus-binding">결합 부위</button>':''}<button class="button" id="reset-structure">전체 구조</button>
    <button class="button" id="save-receptor">수용체 저장</button>${complex?'<button class="button" id="save-ligand">포즈 파일 저장</button><button class="button" id="save-raw">원본 PDBQT 저장</button>':''}</div>
    <p class="structure-note">${observed?'입력 구조에 관측된 리간드 좌표입니다. 결합 부위 버튼으로 주변 포켓을 확대합니다. 선택한 단백질과 기준 리간드의 원래 좌표계를 유지하며, 도킹 결과와 구분합니다.':complex?'채택 SDF와 보관 원본 출력을 구분합니다. 원본의 운영 구간 밖·비유한 점수는 효능이나 비결합 판정이 아닙니다. PDB 표시는 저장 XYZ를 유지하며 결합선은 뷰어가 추정합니다. 결합의 실험적 확인은 분석 근거를 확인하세요.':'분석에 저장된 구조입니다. 준비된 수용체와 PDB 원본을 구분해 확인하세요.'}</p>
    <p id="structure-warnings" class="platform-error" role="alert" hidden></p>
    ${methods?`<details class="structure-methods"><summary>계산 조건과 확인 범위</summary><p>${esc(methods.scope)}</p><pre id="structure-method-record"></pre></details>`:''}
    <div class="structure-canvas"><div id="molstar-view" role="region" aria-label="Mol* 3D 분자 구조"></div><div id="structure-loading"><span class="spinner"></span>구조를 불러오는 중입니다.</div></div>
    <div class="structure-footer"><p id="structure-status" role="status" aria-live="polite"></p><p>드래그: 회전 · 휠/핀치: 확대 · 원자 선택: 정보 보기${complex?' · 수용체: 연한 녹색 / 화합물: 막대 모형':''}</p>
    <details><summary>구조 파일의 출처</summary><p id="structure-source"></p></details><a href="https://molstar.org/" target="_blank" rel="noopener noreferrer">Mol* 5.12.0</a></div></div>`;
  $('#close-structure').onclick = () => dialog.close();
  dialog.onclose = cleanup;
  $('#structure-choice').onchange = e => display(Number(e.target.value));
  $('#reset-structure').onclick = () => session?.viewer?.plugin.managers.camera.reset(undefined,250);
  if ($('#focus-binding')) $('#focus-binding').onclick = () => session?.viewer?.plugin.managers.camera.focusSphere(session.ligandStructure.boundary.sphere,
    {minRadius:12,extraRadius:8,durationMs:250});
  $('#save-receptor').onclick = () => saveAsset(session?.receptor);
  if ($('#save-ligand')) $('#save-ligand').onclick = () => saveAsset(session?.ligand);
  if ($('#save-raw')) $('#save-raw').onclick = async () => {
    const s=session,entry=s?.currentEntry;if(!s?.original||!entry?.original_path)return;
    $('#save-raw').disabled=true;
    try{const asset=await s.original(entry,s.controller.signal);if(session===s)saveAsset(asset)}
    catch(e){if(session===s)$('#structure-status').textContent=e.message}
    finally{if(session===s)$('#save-raw').disabled=!s.currentEntry?.original_path}
  };
  dialog.showModal();
  display(selected);
}

export function showTargetStructure({title,structures,api}) {
  if (!structures.length) return;
  open({title,subtitle:'표적 구조 · Mol*',entries:structures,load:async (entry,signal)=>({receptor:await api(entry.path,{signal})})});
}

export function showResearchComplex({title,resolutionId,api}) {
  open({title,subtitle:'입력 확인 · 원래 관측된 결합 자세',complex:true,observed:true,
    entries:[{label:'관측된 기준 리간드와 선택한 단백질',reason:'입력 구조의 좌표'}],
    load:async(_entry,signal)=>api('/v2/research-inputs/'+encodeURIComponent(resolutionId)+'/structure',{signal})});
}

export function showDockingStructure({title,node,manifest,api}) {
  const entries = dockingEntries(manifest);
  if(!entries.length)return;
  const shaText=async text=>{
    if(!globalThis.crypto?.subtle)throw Error('이 연결에서 구조 파일 해시 검증을 시작할 수 없습니다.');
    return [...new Uint8Array(await crypto.subtle.digest('SHA-256',new TextEncoder().encode(text)))].map(x=>x.toString(16).padStart(2,'0')).join('');
  };
  const selected=Math.max(0,entries.findIndex(e=>e.kind==='accepted'&&e.pose_index===manifest.selected_pose));
  const scope=manifest.identity_scope==='legacy_chemistry_digest_only'
    ?'과거 기록: 수용체 화학적 식별과 SDF 속성만 확인했습니다. 수소 좌표 전체 파일·실제 엔진·검색 조건의 완전한 검증은 미확인입니다.'
    :manifest.identity_scope==='full_receptor_preparation_and_execution'
      ?'수용체 전체 파일·준비·실행 기록의 연결과 원 파일 해시를 확인했습니다. 모델별 표시·점수 연결 여부는 각 항목을 확인하세요. 과거 조건 버전의 상자 정밀도와 미기록 환경은 별도이며, 동일 기록 조건도 친화도 순위 검증은 아닙니다.'
      :manifest.identity_scope==='full_receptor_bytes'
        ?'수용체 전체 파일과 포즈 파일 식별을 확인했습니다. 준비·실제 엔진의 완전한 실행 기록은 미확인입니다.'
        :'구조 식별 검증 범위가 기록되지 않았습니다. 파일·준비·실행 조건을 상세 근거에서 확인해 주세요.';
  open({title,subtitle:`${node} · ${manifest.pdb_id} / ${manifest.chain} · ${manifest.engine} · Mol*`,entries,complex:true,
    selected,methods:{scope,record:{actual_engine:manifest.actual_engine,requested_engine:manifest.requested_engine,
      identity_scope:manifest.identity_scope,declared_protocol_digest:manifest.protocol_digest,
      accepted_or_final_attempt_conditions:manifest.docking_conditions,preparation:manifest.preparation_summary,
      pose_identity_scope:manifest.pose_identity_scope,
      comparison_scope:'Selected raw attempts use their own recorded conditions. A prior GPU attempt may lack full search-input metadata; final CPU conditions do not describe that prior attempt.'}},
    original:async(entry,signal)=>{
      const asset=await api(entry.original_path,{signal});
      if(asset.format!=='pdbqt'||asset.source_sha256!==entry.original_sha256||await shaText(asset.data)!==entry.original_sha256)
        throw Error('원본 도킹 출력이 저장된 해시와 다릅니다. 다시 열어 주세요.');
      return asset;
    },load:async (entry,signal)=>{
      const [receptor,ligand] = await Promise.all([api(manifest.receptor_path,{signal}),entry.renderable&&entry.path?api(entry.path,{signal}):Promise.resolve(null)]);
      if (receptor.source_sha256!==manifest.receptor_sha256 || ligand&&ligand.source_sha256!==entry.source_sha256)
        throw Error('서버의 구조 파일이 변경되었습니다. 이 창을 닫고 결합 구조를 다시 열어 주세요.');
      if(!receptor.projection_sha256||await shaText(receptor.data)!==receptor.projection_sha256)
        throw Error('표시할 수용체 파일의 해시가 다릅니다. 다시 열어 주세요.');
      if(ligand&&await shaText(ligand.data)!==(entry.kind==='raw'?entry.projection_sha256:entry.source_sha256))
        throw Error('표시할 포즈 파일이 저장된 해시와 다릅니다. 다시 열어 주세요.');
      return {receptor,ligand};
    }});
}

// MolScript expression format is provided by the pinned Mol* 5.12.0 builder.
export function correspondenceExpression(indices) {
  if (!Array.isArray(indices) || indices.some(i=>!Number.isInteger(i)||i<0)) throw Error('잘못된 대응 원자 인덱스');
  return {head:{name:'structure-query.generator.atom-groups'},args:{'atom-test':{
    head:{name:'core.set.has'},args:[{head:{name:'core.type.set'},args:[...new Set(indices)]},
      {head:{name:'structure-query.atom-property.core.source-index'}}]}}};
}

export async function showCorrespondence({title,result}) {
  const molstar=await loadLibrary();
  let dialog=document.querySelector('#correspondence-dialog');
  if(dialog?.open)dialog.close();
  if(!dialog){dialog=document.createElement('dialog');dialog.id='correspondence-dialog';document.body.append(dialog)}
  dialog.setAttribute('aria-labelledby','correspondence-title');
  dialog.innerHTML=`<div class="structure-layout"><div class="structure-header"><div><small>스카폴드 대응 · MolStar</small><h2 id="correspondence-title">${esc(title)}</h2></div><button class="icon-button" id="close-correspondence" aria-label="대응 창 닫기">×</button></div><p class="structure-note">기준 구조와 후보의 설명용 3D 정렬입니다. 저장 도킹 포즈나 결합 증거가 아닙니다. 원자 번호는 입력 구조 기준 0부터 시작합니다.</p><div class="structure-controls"><label><input type="checkbox" id="color-mcs" checked> 공통 부분구조 · 파랑</label><label><input type="checkbox" id="color-features"> 3D 약리 특징 대응 · 특징별 색</label><span>미대응: 기준 회색 / 후보 연한 갈색</span><p class="feature-legend">3D 특징 색: Donor 파랑 · Acceptor 자홍 · Aromatic 황토 · Hydrophobe 보라 · PosIonizable 녹색 · NegIonizable 갈색. 특징명·원자 대응 표를 함께 확인하세요.</p></div><div class="structure-canvas"><div id="correspondence-view" role="region" aria-label="기준과 후보의 원자별 대응"></div></div><div class="structure-footer"><p id="correspondence-status" role="status"></p><details open><summary>대응 원자·거리·미대응 특징</summary><div class="table-wrap"><table class="platform-table"><thead><tr><th>특징</th><th>기준 원자</th><th>후보 원자</th><th>거리 Å</th></tr></thead><tbody>${(result.feature_pairs||[]).map(p=>`<tr><td>${esc(p.reference.family)}</td><td>${esc(p.reference.atom_indices.join(', '))}</td><td>${esc(p.candidate.atom_indices.join(', '))}</td><td>${Number(p.distance_angstrom).toFixed(3)}</td></tr>`).join('')}</tbody></table></div><p>공통 원자: ${esc((result.mcs?.atom_pairs||[]).map(p=>`${p.reference_atom}↔${p.candidate_atom}`).join(', '))}</p><p>미대응 특징: 기준 ${(result.unmatched_reference_features||[]).length} / 후보 ${(result.unmatched_candidate_features||[]).length}</p><pre>${esc(JSON.stringify({reference:result.unmatched_reference_features,candidate:result.unmatched_candidate_features,limits:result.limitations},null,2))}</pre></details></div></div>`;
  let viewer=null,closed=false,busy=false;
  const focus=document.activeElement;
  dialog.onclose=()=>{closed=true;viewer?.dispose();if(focus?.isConnected)focus.focus()};
  dialog.querySelector('#close-correspondence').onclick=()=>dialog.close();dialog.showModal();
  if(!result.reference_sdf||!result.candidate_sdf){dialog.querySelector('#correspondence-status').textContent=`3D 대응 미평가: ${result.feature_status}`;return}
  viewer=await molstar.Viewer.create('correspondence-view',{layoutIsExpanded:false,layoutShowControls:false,layoutShowRemoteState:false,layoutShowSequence:false,layoutShowLog:false,layoutShowLeftPanel:false,viewportShowExpand:false,volumeStreamingDisabled:true});
  if(closed){viewer.dispose();return}
  viewer.plugin.canvas3d?.setProps({renderer:{backgroundColor:0xf7faf8}});
  const colors={Donor:0x397bc0,Acceptor:0xc04b70,Aromatic:0xc28a24,Hydrophobe:0x8064bc,LumpedHydrophobe:0x8064bc,PosIonizable:0x00866b,NegIonizable:0x8f4920};
  const draw=async()=>{
    if(busy||closed)return;busy=true;
    const inputs=[...dialog.querySelectorAll('.structure-controls input')];inputs.forEach(x=>x.disabled=true);
    try{
      await viewer.plugin.clear();
      for(const [role,sdf,base] of [['reference',result.reference_sdf,0x8b9599],['candidate',result.candidate_sdf,0xc4a486]]){
        if(closed)return;
        const raw=await viewer.plugin.builders.data.rawData({data:sdf,label:role==='reference'?'기준 화합물':'후보 화합물'});
        const trajectory=await viewer.plugin.builders.structure.parseTrajectory(raw,'sdf');
        const model=await viewer.plugin.builders.structure.createModel(trajectory);
        const structure=await viewer.plugin.builders.structure.createStructure(model,{name:'model',params:{}});
        const add=async(indices,key,color,label)=>{if(!indices.length)return;const component=await viewer.plugin.builders.structure.tryCreateComponentFromExpression(structure,correspondenceExpression(indices),key,{label});if(component)await viewer.plugin.builders.structure.representation.addRepresentation(component,{type:'ball-and-stick',color:'uniform',colorParams:{value:color},typeParams:{sizeFactor:.34,ignoreHydrogens:true}})};
        const all=await viewer.plugin.builders.structure.tryCreateComponentStatic(structure,'all',{label:role});
        await viewer.plugin.builders.structure.representation.addRepresentation(all,{type:'ball-and-stick',color:'uniform',colorParams:{value:base},typeParams:{sizeFactor:.25,ignoreHydrogens:true}});
        if(dialog.querySelector('#color-mcs').checked)await add(result.mcs.atom_pairs.map(p=>p[role+'_atom']),role+'-mcs',0x1f6ec4,'공통 부분구조');
        if(dialog.querySelector('#color-features').checked){const grouped=new Map();for(const pair of result.feature_pairs){const f=pair[role];grouped.set(f.family,[...(grouped.get(f.family)||[]),...f.atom_indices])}for(const [family,indices] of grouped)await add(indices,role+'-'+family,colors[family]||0xcd7436,family)}
      }
      viewer.handleResize();viewer.plugin.managers.camera.reset(undefined,0);
      dialog.querySelector('#correspondence-view').dataset.loaded='true';
      dialog.querySelector('#correspondence-status').textContent=`공통 원자 ${result.mcs.atom_pairs.length}쌍 · 3D 특징 ${result.feature_pairs.length}쌍 · ${result.mcs.complete_search?'공통 구조 검색 완료':'시간 제한 내 부분 대응'}`;
    }catch(e){if(!closed)dialog.querySelector('#correspondence-status').textContent=e.message}
    finally{busy=false;inputs.forEach(x=>x.disabled=false)}
  };
  dialog.querySelector('#color-mcs').onchange=draw;dialog.querySelector('#color-features').onchange=draw;
  await draw();
}


export function showPlatformDocking({title,candidate,api,context}) {
 const entries=(candidate.docking||[]).map((row,index)=>({...row,observation_index:index,pose_index:index+1,kind:'accepted',renderable:true,
  label:`${row.engine||'엔진'} · ${{broad_docking:'넓은 도킹',confirmation:'재확인',receptor_sensitivity:'수용체 민감성'}[row.source_stage||row.stage]||row.source_stage||row.stage||'도킹'} · 시드 ${row.source_seed??row.seed??row.protocol?.seed??'미기록'} · ${Number.isFinite(row.score)?row.score.toFixed(3):'미기록'} kcal/mol`})).filter(r=>r.artifact_valid&&r.pose_sha256);
 if(!entries.length)throw Error('검증된 도킹 포즈가 없습니다.');
 open({title,subtitle:`${candidate.target_id} · 내용 해시 검증 포즈 · MolStar`,entries,complex:true,
  methods:{scope:'수용체·포즈·프로토콜은 선택한 불변 스냅샷에 연결됩니다. 구조 계산 결과이며 실험적 결합 검증은 아닙니다.',record:{snapshot_id:context.snapshot,target_id:context.target}},
  load:async(entry,signal)=>{
   const payload=await api(`/v2/candidates/${encodeURIComponent(candidate.candidate_id)}/docking?snapshot_id=${encodeURIComponent(context.snapshot)}&target_id=${encodeURIComponent(context.target)}&observation_index=${entry.observation_index}`,{signal});
   for(const asset of [payload.receptor,payload.ligand,payload.reference,payload.original_pose].filter(Boolean)){
    const hash=[...new Uint8Array(await crypto.subtle.digest('SHA-256',new TextEncoder().encode(asset.data)))].map(v=>v.toString(16).padStart(2,'0')).join('');
    if(hash!==(asset.projection_sha256||asset.source_sha256))throw Error('표시 구조가 스냅샷의 파일 해시와 다릅니다.');
    asset.source_file=asset.source_file||asset.filename||'structure.'+asset.format;
   }
   entry.original_asset=payload.original_pose||null;entry.original_path=payload.original_pose?'cas:'+entry.pose_sha256:null;
   return payload;
  },original:async entry=>{if(!entry.original_asset)throw Error('원본 포즈가 제공되지 않았습니다.');return entry.original_asset}});
}
