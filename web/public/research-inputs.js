// Input contracts shared by the novice flow and its request-level checks.
export const escapeHtml=value=>String(value??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
export const PRODUCTS={body_lotion:'바디 로션',face_cream:'페이스 크림',hand_cream:'핸드 크림',lip_product:'입술 제품',body_wash:'바디 워시',shampoo:'샴푸'};
export const GOALS={exploratory:'아직 정하지 않음',surface:'피부 표면·장벽',epidermal:'표피',dermal:'진피'};
export function researchInputRequest({kind='compound',format='name',value='',taxon=9606,action='any',pubchemCid=null,selection=null,ligandSmiles=''}={}){
 const text=String(value).trim();
 if(!text)throw Error(kind==='bound_complex'?'PDB ID를 입력하거나 구조 파일을 올려 주세요.':'먼저 이름 또는 분자 구조를 입력해 주세요.');
 if(kind==='protein'){
  const taxonId=Number(taxon);
  if(!Number.isInteger(taxonId)||taxonId<1)throw Error('종의 NCBI taxon ID를 확인해 주세요.');
  if(!['any','activate','inhibit'].includes(action))throw Error('원하는 작용을 확인해 주세요.');
  return{kind,format:'name',value:text,taxon_id:taxonId,desired_action:action};
 }
 if(kind==='compound'){
  if(!['name','smiles'].includes(format))throw Error('화합물 이름 또는 SMILES 입력을 선택해 주세요.');
  const input={kind,format,value:text};
  if(pubchemCid!==null){
   const cid=Number(pubchemCid);
   if(format!=='name'||!Number.isSafeInteger(cid)||cid<1)throw Error('화합물 선택 정보를 다시 확인해 주세요.');
   input.pubchem_cid=cid;
  }
  return input;
 }
 if(kind==='bound_complex'){
  if(!['pdb_id','upload'].includes(format))throw Error('PDB ID 또는 구조 파일을 선택해 주세요.');
  return{kind,format,value:text,...(selection?{selection}:{}),...(ligandSmiles.trim()?{ligand_smiles:ligandSmiles.trim()}:{})};
 }
 throw Error('입력 종류를 다시 선택해 주세요.');
}
export function skinContextRequest({product=null,appliedPercent=null,goal='exploratory',researchGoal='',applicationSite=''}={}){
 const selected=product||null;
 if(selected!==null&&!Object.hasOwn(PRODUCTS,selected))throw Error('제품 형태를 다시 선택해 주세요.');
 if(!Object.hasOwn(GOALS,goal))throw Error('주로 살펴볼 위치를 다시 선택해 주세요.');
 const pct=appliedPercent===null||appliedPercent===''?null:Number(appliedPercent);
 if(pct!==null&&(typeof appliedPercent==='boolean'||!Number.isFinite(pct)||pct<0||pct>100))throw Error('배합 농도는 비워 두거나 0–100% 숫자로 입력해 주세요.');
 return{product:selected,applied_percent:pct,goal,...(researchGoal.trim()?{research_goal:researchGoal.trim()}:{}),...(applicationSite.trim()?{application_site:applicationSite.trim()}:{})};
}
export function generationRequest({enabled=false,products=200,max3d=50,conformers=10,maxSa=6}={}){
 if(!enabled)return{enabled:false};
 const fields={brics_products:Number(products),max_3d:Number(max3d),conformers:Number(conformers),max_sa:Number(maxSa)};
 for(const [key,max]of[['brics_products',2000],['max_3d',200],['conformers',30]])if(!Number.isInteger(fields[key])||fields[key]<1||fields[key]>max)throw Error('새 구조 생성의 수량 조건을 확인해 주세요.');
 if(!Number.isFinite(fields.max_sa)||fields.max_sa<1||fields.max_sa>10)throw Error('합성 난이도 조건은 1–10 사이로 입력해 주세요.');
 return{enabled:true,...fields};
}
export function screeningRequest(target=200){
 const value=Number(target);
 if(typeof target==='boolean'||!Number.isInteger(value)||value<1||value>1000)throw Error('대표 후보 목표는 1–1,000개의 정수로 입력해 주세요.');
 return{target_candidates:value};
}
export function identityLabel(result){
 const input=result?.input||result?.original_input||result||{},kind=result?.kind||input.kind;
 if(kind==='bound_complex'&&input.format==='upload'){
  const filename=result?.status==='verified'&&result?.upload_receipt?.filename;
  return typeof filename==='string'&&filename?'업로드한 결합 구조 · '+filename:'업로드한 결합 구조';
 }
 const name=result?.name||result?.preferred_name||result?.title||result?.gene||input.value||result?.request?.query||result?.target_id;
 if(kind==='compound'&&input.format==='smiles'&&(!name||name===input.value))return'입력한 화합물 (SMILES)';
 return typeof name==='string'&&/^[a-f0-9]{64}$/i.test(name)?'이름이 기록되지 않은 화합물':name||'입력 확인';
}
const e=escapeHtml;
const INPUT_LIMITATIONS={
 'Structural similarity and skin-property predictions do not establish activity or skin benefit.':'구조 유사도와 피부 관련 물성 예측을 기준으로 후보를 비교합니다. 활성 자료는 별도 항목에서 확인하세요.',
 'A bound pose provides a structural hypothesis; it does not establish activity or skin benefit.':'관측된 결합 자세를 기준으로 후보의 형태와 포켓 내 배치를 비교합니다.'
};
export function researchSourceUnavailable(result={}){
 if(result.status==='verified')return false;
 const sourceReasons=['structure_source_unavailable','compound_identity_source_unavailable','protein_name_source_unavailable'];
 const failedReceipt=Array.isArray(result.source_receipts)&&result.source_receipts.some(receipt=>receipt?.status==='source_unavailable');
 return result.status==='source_unavailable'||sourceReasons.includes(result.reason)||result.reason==='structure_identity_unresolved'&&(failedReceipt||['ReadTimeout','ConnectTimeout','Timeout','ConnectionError','SSLError','ProxyError'].includes(result.error_type));
}
export function researchResolutionGuidance(result={}){
 if(researchSourceUnavailable(result)){
  const kind=result.kind||result.input?.kind;
  return kind==='bound_complex'?'공개 구조 제공처의 응답을 받지 못했습니다. 잠시 뒤 다시 확인하거나 PDB/mmCIF 파일을 올려 주세요.':kind==='compound'?'공식 화합물 조회 서비스의 응답을 받지 못했습니다. 잠시 뒤 이름을 다시 확인하거나, 원자료에서 확인한 정확한 SMILES를 입력해 주세요.':kind==='protein'?'공식 단백질 조회 서비스의 응답을 받지 못했습니다. 잠시 뒤 같은 이름을 다시 확인해 주세요.':'공식 조회 서비스의 응답을 받지 못했습니다. 잠시 뒤 같은 입력을 다시 확인해 주세요.';
 }
 const reasons={ligand_exact_chemistry_required:'구조 파일에는 결합 화합물의 정확한 결합 차수·입체화학 정보가 부족합니다. 같은 화합물의 검증된 SMILES를 아래에 입력하고 다시 확인하세요. 모르면 원자료 제공자에게 요청하거나, 같은 공개 구조의 PDB ID로 시작하세요. 좌표나 이름으로 구조를 추측하지 않습니다.',bound_ligand_selection_required:'연구 기준으로 삼을 결합 화합물이 여러 개 있습니다. 화합물 코드·사슬·위치를 확인하고 하나를 직접 선택하세요.',assembly_selection_required:'구조의 조립체가 여러 개입니다. 원래 실험에 해당하는 구성을 고급 선택에서 지정하세요.',ligand_altloc_selection_required:'결합 화합물의 관측 위치가 여러 개입니다. 고급 선택에서 실험에 해당하는 대체 위치를 지정하세요.',bound_ligand_multiple_organic_components:'선택한 화합물에 유기 성분이 여러 개 있습니다. 단일한 연구 기준 구조를 확인해 주세요.'};
 const inputReasons={invalid_smiles:'분자 구조를 읽지 못했습니다. 원자료의 정확한 SMILES를 확인하고 다시 입력하세요.',concrete_structure_required:'한 가지 분자를 확정하지 못했습니다. 혼합물 대신 정확한 분자 구조를 지정해 주세요.',multiple_organic_components_require_selection:'유기 성분이 여러 개입니다. 실제 연구 기준 분자를 확인하고 정확한 구조를 입력하세요.',compound_class_requires_exact_structure:'이 이름은 여러 화합물을 포함하는 분자군입니다. 정확한 종류 이름 또는 검증된 SMILES를 입력하세요. 단백질로 바꾸어 분석하지 않습니다.',identity_not_found:'일치하는 대상을 찾지 못했습니다. 철자·입력 종류를 확인하거나 정확한 구조·식별자를 사용해 주세요.',name_match_limit_requires_exact_structure:'이 이름에 해당하는 물질이 너무 많습니다. 구체적인 이름이나 정확한 SMILES로 범위를 좁혀 주세요.',selected_cid_not_in_name_results:'선택한 화합물이 현재 이름 검색 결과와 맞지 않습니다. 이름을 다시 검색하고 구조를 확인하세요.',ambiguous_compound_name:'같은 이름에 해당하는 화합물이 여러 개입니다. 아래에서 구조와 이름을 확인하고 연구할 물질을 선택하세요.',compound_identity_source_unavailable:'공식 화합물 조회 서비스에 연결하지 못했습니다. 잠시 뒤 다시 확인하거나 정확한 SMILES를 입력하세요.',protein_family_requires_exact_member:'이 이름은 단백질 가족이나 분자군입니다. 연구할 단백질 한 종류의 이름을 입력하세요.',ambiguous_protein_name:'같은 이름에 해당하는 단백질이 여러 개입니다. 원자료의 정확한 구성원·생물종을 확인한 뒤 입력을 구체화하세요.',protein_name_source_unavailable:'공식 단백질 조회 서비스에 연결하지 못했습니다. 잠시 뒤 다시 확인해 주세요.',protein_name_page_incomplete:'검색 결과가 너무 많아 하나로 확정하지 못했습니다. 정확한 단백질 이름으로 범위를 좁혀 주세요.'};
 return reasons[result.reason]||inputReasons[result.reason]||result.message||'입력 이름과 종류를 확인하거나, 정확한 SMILES·구조 파일을 사용해 주세요.';
}
function boundIdentityMarkup(result){
 const selected=result.selection||{},chains=Array.isArray(selected.protein_chains)?selected.protein_chains.join(', '):'미확정';
 const assembly=selected.assembly_id==='as_uploaded'?'업로드한 좌표 그대로':selected.assembly_id==='asymmetric_unit'?'구조의 비대칭 단위':selected.assembly_id?'구조에 기록된 조립체 '+selected.assembly_id:'미확정';
 return`<div class="research-options"><strong>확인할 결합 구조의 구성</strong><p>모델 ${e(selected.model_id||'미확정')} · 단백질 사슬 ${e(chains)} · 결합 화합물 ${e(selected.ligand_code||'미확정')} (사슬 ${e(selected.ligand_chain??'미확정')}, 잔기 ${e(selected.ligand_residue??'미확정')})</p><p>조립체: ${e(assembly)}</p><p>${result.actual_protein_id?'확인한 단백질 식별자: '+e(result.actual_protein_id):'단백질 이름·종은 아직 확정하지 않았습니다. 선택한 사슬의 원래 좌표를 계산에 사용합니다.'}</p><p>화합물의 정확한 화학 구조는 ${result.chemistry_receipts?.basis?'저장한 구조 원자료와의 대조 기록이 있습니다.':'식별 근거에서 확인해 주세요.'} 관측된 결합 좌표만으로 효능이 확인되지는 않습니다.</p>${(result.docking_hold_reasons||[]).length?'<p>구조 탐색은 가능하지만 도킹에는 추가 준비가 필요합니다. 실행 결과에서 구체적인 보류 사유를 확인해 주세요.</p>':''}</div>`;
}
export function researchIdentityMarkup(result){
 const kind=result.kind||result.input?.kind||'protein',compound=kind==='compound',bound=kind==='bound_complex';
 const verified=result.status==='verified',unavailable=researchSourceUnavailable(result),choices=Array.isArray(result.candidates)?result.candidates:Array.isArray(result.choices)?result.choices:[];
 const smiles=result.canonical_smiles||result.smiles||result.ligand?.canonical_smiles||result.ligand?.smiles;
 const status=verified?'확인한 입력':unavailable?'공식 자료 조회를 완료하지 못했습니다':result.status==='selection_required'?'구체적인 항목 선택이 필요합니다':'아직 입력을 확정하지 못했습니다';
 return`<div class="research-identity"><div><span class="eyebrow">${e(status)}</span><h2>${e(identityLabel(result))}</h2><p>${unavailable?'공식 자료 조회를 마친 뒤 이 입력의 계산 대상을 확인합니다.':compound?'이 화합물을 출발구조로 사용합니다. 단백질을 자동으로 선택하지 않습니다.':bound?'선택한 결합 화합물과 단백질 사슬을 바탕으로 연구를 준비합니다.':'이 단백질을 계산 표적으로 사용합니다.'}</p>${result.organism?`<p>종: ${e(result.organism)}</p>`:''}${result.target_id?`<p>식별자: <code>${e(result.target_id)}</code></p>`:''}${result.pubchem_cid?`<p>PubChem CID: ${e(result.pubchem_cid)}</p>`:''}${smiles?`<details><summary>확인한 분자 구조 (SMILES)</summary><code class="research-smiles">${e(smiles)}</code></details>`:''}${bound&&verified?boundIdentityMarkup(result):''}${bound&&verified?'<button type="button" class="button" id="research-bound-view">결합 화합물과 주변 포켓을 3D로 확인</button>':''}</div>${verified&&(compound||bound)?'<div id="research-identity-image" class="research-molecule" aria-label="확인한 화합물의 2차원 구조"></div>':''}</div>
 ${!verified?`<div class="notice"><p>${e(researchResolutionGuidance(result))}</p>${unavailable?'':`<p>${compound?'여러 분자를 포함하는 물질군 이름은 정확한 종류 또는 SMILES가 필요합니다. 표시된 구조를 확인하고 직접 선택해 주세요.':'아래 선택 항목을 확인하거나 입력을 수정해 주세요.'}</p>`}</div>`:''}
 ${!verified&&compound&&choices.length?`<div class="research-choice-list">${choices.filter(c=>Number.isSafeInteger(Number(c.pubchem_cid??c.cid))&&Number(c.pubchem_cid??c.cid)>0).map(c=>`<button type="button" class="research-choice" data-compound-cid="${Number(c.pubchem_cid??c.cid)}"><strong>${e(c.name||c.title||'화합물')}</strong><span>PubChem CID ${e(c.pubchem_cid??c.cid)}</span><code>${e(c.smiles||c.canonical_smiles||'구조를 확인한 후 선택')}</code><span>이 화합물 확인 →</span></button>`).join('')}</div>`:''}
 <div id="research-bound-choices"></div>${(Array.isArray(result.limitations)?result.limitations:[]).map(x=>`<p class="research-muted">${e(Object.hasOwn(INPUT_LIMITATIONS,x)?INPUT_LIMITATIONS[x]:x)}</p>`).join('')}<details class="research-evidence"><summary>식별 근거와 원문 기록</summary><pre>${e(JSON.stringify(result,null,2))}</pre></details>`;
}
export function boundSelectionRequest({modelId,assemblyId='',proteinChains=[],ligandChain,ligandResidue,ligandCode,altloc=''}={}){
 if(modelId===undefined||!Array.isArray(proteinChains)||!proteinChains.length||ligandChain===undefined||ligandResidue===undefined||!ligandCode)throw Error('모델·단백질 사슬·결합 화합물을 선택해 주세요.');
 return{model_id:String(modelId),...(assemblyId?{assembly_id:String(assemblyId)}:{}),protein_chains:proteinChains.map(String),ligand_chain:String(ligandChain),ligand_residue:String(ligandResidue),ligand_code:String(ligandCode),...(altloc?{altloc:String(altloc)}:{})};
}
export function boundChoiceMarkup(result){
 const inventory=result.inventory||{},ligands=Array.isArray(result.choices)?result.choices:(inventory.ligands||[]),chains=inventory.protein_chains||[],selected=result.selection||{},chemistryPending=result.reason==='ligand_exact_chemistry_required';
 const isSelected=l=>selected.ligand_code===l.code&&String(selected.model_id)===String(l.model_id)&&String(selected.ligand_chain)===String(l.chain_id)&&String(selected.ligand_residue)===String(l.residue_id);
 if(!ligands.length)return'<p>선택할 결합 화합물이 없습니다. 결합 화합물을 포함한 PDB/mmCIF 파일을 사용해 주세요.</p>';
 return`<form id="research-bound-form" class="research-options"><h3>어느 결합 화합물을 살펴볼까요?</h3><p>첫 번째 항목은 자동으로 선택하지 않습니다. 실험에 사용한 단백질 사슬과 결합 화합물을 직접 확인해 주세요.</p><label>결합 화합물<select id="research-bound-ligand" required><option value="">선택해 주세요</option>${ligands.map((l,i)=>`<option value="${i}" ${isSelected(l)?'selected':''}>${e(l.code)} · 사슬 ${e(l.chain_id||'(없음)')} · 잔기 ${e(l.residue_id)} · 모델 ${e(l.model_id)}</option>`).join('')}</select></label><fieldset><legend>계산에 포함할 단백질 사슬</legend>${chains.map((c,i)=>`<label><input type="checkbox" name="bound-chain" value="${i}" ${selected.protein_chains?.includes(c.chain_id)&&String(selected.model_id)===String(c.model_id)?'checked':''}>사슬 ${e(c.chain_id||'(없음)')} · ${e(c.residue_count)}개 잔기 · 모델 ${e(c.model_id)}</label>`).join('')||'<p>단백질 사슬을 찾지 못했습니다.</p>'}</fieldset><div class="research-chemistry-help${chemistryPending?' research-chemistry-required':''}"><label>결합 화합물의 정확한 SMILES ${chemistryPending?'(이번 파일에 필요)':'(화학 정보가 부족한 업로드에 사용)'}<input id="research-bound-smiles" placeholder="원자료에서 확인한 정확한 SMILES" ${chemistryPending?'required':''}></label><p>원자 좌표만으로 결합 차수·입체화학을 확정할 수 없는 파일에 필요합니다. 직접 만들거나 ChatGPT가 추측한 구조를 넣지 말고, 같은 화합물의 확인된 구조를 사용하세요.</p></div><details><summary>고급 선택</summary><label>생물학적 조립체<select id="research-bound-assembly"><option value="">아직 지정하지 않음 · 확인 결과에서 기본 구성 확인</option>${(inventory.assembly_ids||[]).map(id=>`<option value="${e(id)}" ${String(selected.assembly_id)===String(id)?'selected':''}>${e(id)}</option>`).join('')}</select></label><p>적용한 조립체는 다음 확인 카드에 표시합니다. 공개 자료의 조립체를 선택하면 구조 확장이 필요하며, 지원되지 않는 경우 도킹을 보류합니다.</p><label>대체 위치 (있는 경우)<input id="research-bound-altloc" placeholder="예: A" value="${e(selected.altloc||'')}"></label></details><button class="button primary" type="submit">선택한 결합 화합물 확인</button><div id="research-bound-error" role="alert"></div></form>`;
}
