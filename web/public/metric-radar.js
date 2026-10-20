// Display-only profiles of saved values. Ranges control drawing, never filtering.
import {helpButton} from './column-help.js';
const esc=value=>String(value??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
export const PROPERTY_AXES=[
 {key:'logP',label:'지용성 · logP',unit:'log 비율',min:-3,max:8,help:'review_logp'},
 {key:'molecular_weight',label:'분자량',unit:'g/mol',min:0,max:800,positive:true,help:'review_mw'},
 {key:'tpsa',label:'극성 · TPSA',unit:'Å²',min:0,max:200,nonnegative:true,help:'review_tpsa'},
 {key:'hydrogen_bond_donors',label:'수소결합 공여체',unit:'개',min:0,max:8,count:true,help:'review_hbd'},
 {key:'hydrogen_bond_acceptors',label:'수소결합 수용체',unit:'개',min:0,max:15,count:true,help:'review_hba'},
 {key:'QED',label:'약물 유사성 · QED',unit:'0–1',min:0,max:1,bounded:true,help:'review_qed'},
];
// Highest component-score bands in filters/topical.py; HBA/QED have no band.
export const SKIN_REFERENCE = {
 "logP":{"min":1,"max":3.5},
 "molecular_weight":{"min":0,"max":350,"minExclusive":true},
 "tpsa":{"min":0,"max":90},
 "hydrogen_bond_donors":{"min":0,"max":3}
};
export const SKIN_REFERENCE_SOURCE='https://github.com/kangk1204/SkinPipe_public/blob/main/skinpipe/skinpipe/filters/topical.py';
export const TOXICITY_AXES=[
 ['Skin_Reaction','피부 감작성'],['AMES','변이원성 · AMES'],['hERG','hERG 채널'],
 ['DILI','간 손상 · DILI'],['ClinTox','임상 독성'],['Carcinogens_Lagunin','발암성'],
].map(([key,label])=>({key,label,unit:'0–1 모델 출력',min:0,max:1,bounded:true,help:'review_tox_'+key}));
const num=value=>typeof value==='number'&&Number.isFinite(value);
const display=value=>value.toLocaleString('ko-KR',{maximumSignificantDigits:6});
const shortReason=reason=>({'지표 정의 미확인':'정의 미확인','품질 경고 · 원자료 확인':'품질 경고','값 오류 · 원자료 확인':'값 오류'})[reason]||reason;
const axisValue=axis=>axis.value===null?shortReason(axis.reason):display(axis.value)+(axis.unit.startsWith('0–1')?' / 1':axis.unit==='log 비율'?'':' '+axis.unit)+(axis.clipped?' †':'');
export function radarProfile(kind,values={}, {allowedKeys=null,blockedKeys=[]}={}){
 const axes=kind==='properties'?PROPERTY_AXES:kind==='toxicity'?TOXICITY_AXES:null;
 if(!axes)throw Error('Unknown radar profile');
 const allowed=allowedKeys===null?null:new Set(allowedKeys),blocked=new Set(blockedKeys);
 return axes.map(axis=>{
  const value=values?.[axis.key];let reason=null;
  if(allowed&&!allowed.has(axis.key))reason='지표 정의 미확인';
  else if(blocked.has(axis.key))reason='품질 경고 · 원자료 확인';
  else if(value===null||value===undefined)reason='미평가';
  else if(!num(value)||(axis.positive&&value<=0)||(axis.nonnegative&&value<0)||(axis.count&&(!Number.isInteger(value)||value<0))||(axis.bounded&&(value<axis.min||value>axis.max)))reason='값 오류 · 원자료 확인';
  return {...axis,value:reason?null:value,reason,normalized:reason?null:Math.max(0,Math.min(1,(value-axis.min)/(axis.max-axis.min))),clipped:!reason&&(value<axis.min||value>axis.max)};
 });
}
let figureId=0;
const point=(index,radius)=>{const angle=-Math.PI/2+index*Math.PI/3;return[230+Math.cos(angle)*radius,186+Math.sin(angle)*radius]};
const coordinate=xy=>xy.map(value=>value.toFixed(2)).join(',');
export function metricRadarMarkup({kind,values={},source='출처 미기록',allowedKeys=null,blockedKeys=[]}){
 const axes=radarProfile(kind,values,{allowedKeys,blockedKeys}),valid=axes.filter(a=>a.value!==null),complete=valid.length===6;
 const id='metric-radar-'+(++figureId),title=kind==='properties'?'분자 물성 프로필':'독성 예측 6항목';
 const references=axes.map(a=>kind==='properties'?SKIN_REFERENCE[a.key]||null:null);
 const referenceNote=kind==='properties'?'주황색 영역·점선은 SkinPipe 피부 물성 점수에서 각 항목이 최고점인 범위입니다. 기준이 있는 4축에 표시하며 HBA·QED의 참고 범위는 미설정입니다.':'';
 const direction=kind==='properties'?'바깥쪽일수록 각 물성값이 큽니다. 적정 범위는 연구 목적·제형과 함께 확인하세요.':'바깥쪽일수록 해당 모델의 독성 양성 출력이 큽니다. 각 항목의 판정 기준은 ?에서 확인하세요.';
 const description=axes.map(a=>a.label+': '+(a.value===null?a.reason:display(a.value)+' '+a.unit)).join('; ');
 const grid=[.25,.5,.75,1].map(r=>`<polygon points="${axes.map((_,i)=>coordinate(point(i,108*r))).join(' ')}"/>`).join('');
 const spokes=axes.map((_,i)=>`<line x1="230" y1="186" x2="${point(i,108)[0]}" y2="${point(i,108)[1]}"/>`).join('');
 const referencePoint=(index,bound)=>point(index,108*(references[index][bound]-axes[index].min)/(axes[index].max-axes[index].min));
 const referenceBand=references.map((ref,i)=>{const next=(i+1)%6;if(!ref||!references[next])return'';return`<polygon class="radar-reference-band" data-reference-axes="${axes[i].key},${axes[next].key}" points="${[referencePoint(i,'max'),referencePoint(next,'max'),referencePoint(next,'min'),referencePoint(i,'min')].map(coordinate).join(' ')}"/>`;}).join('');
 const referenceAxes=references.map((ref,i)=>!ref?'':`<line class="radar-reference-axis" data-reference-key="${axes[i].key}" x1="${referencePoint(i,'min')[0]}" y1="${referencePoint(i,'min')[1]}" x2="${referencePoint(i,'max')[0]}" y2="${referencePoint(i,'max')[1]}"/>`).join('');
 const positions=axes.map((a,i)=>a.value===null?null:point(i,a.normalized*108));
 const segments=positions.flatMap((p,i)=>p&&positions[(i+1)%6]?[`<line class="radar-value-line" x1="${p[0]}" y1="${p[1]}" x2="${positions[(i+1)%6][0]}" y2="${positions[(i+1)%6][1]}"/>`]:[]).join('');
 const polygon=complete?`<polygon class="radar-value-area" points="${positions.map(coordinate).join(' ')}"/>`:'';
 const dots=axes.map((a,i)=>a.value===null?'':`<circle class="radar-value-dot" data-endpoint="${a.key}" data-raw-value="${a.value}" cx="${positions[i][0]}" cy="${positions[i][1]}" r="3.5"><title>${esc(a.label+': '+a.value+' '+a.unit+(a.clipped?' · 표시 범위 밖':''))}</title></circle>`).join('');
 const labels=axes.map((a,i)=>{const [x,y]=point(i,143),anchor=i===0||i===3?'middle':i<3?'start':'end';return`<text x="${x}" y="${y-5}" text-anchor="${anchor}"><tspan class="radar-axis-name">${esc(a.label)}</tspan><tspan x="${x}" dy="18">${esc(axisValue(a))}</tspan>${references[i]?`<tspan class="radar-axis-reference" x="${x}" dy="17">참고 ${esc(references[i].min===0?'≤ '+references[i].max:references[i].min+'–'+references[i].max)}</tspan>`:''}</text>`}).join('');
 const metadata={schema:'skinpipe-display-radar-v1',kind,source:String(source).slice(0,500),reference:kind==='properties'?{policy:'skinpipe-topical-max-component-v1',source:SKIN_REFERENCE_SOURCE,bands:SKIN_REFERENCE}:null,axes:axes.map(({key,value,min,max,reason,clipped})=>({key,value,min,max,reason,clipped}))};
 const svg=`<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 460 370" class="metric-radar-svg" role="img" aria-labelledby="${id}-title" aria-describedby="${id}-description"><title id="${id}-title">${esc(title+' · '+valid.length+'/6개 지표')}</title><desc id="${id}-description">${esc(direction+' '+referenceNote+' '+description)}</desc><metadata>${esc(JSON.stringify(metadata))}</metadata><style>.radar-grid{fill:none;stroke:#d5dedb;stroke-width:1}.radar-reference-band{fill:#f4b84b;fill-opacity:.4;stroke:#a35c00;stroke-width:2.5;stroke-dasharray:6 4;stroke-linejoin:round}.radar-reference-axis{stroke:#a35c00;stroke-width:3;stroke-dasharray:4 3}.metric-radar-svg .radar-axis-reference{fill:#8b4c00;font-size:12px;font-weight:650}.radar-value-area{fill:#397668;fill-opacity:.06;stroke:none}.radar-value-line{stroke:#146a59;stroke-width:2.7;fill:none}.radar-value-dot{fill:#286354;stroke:white;stroke-width:1}.metric-radar-svg text{fill:#263c35;font:13px system-ui,sans-serif}.radar-axis-name{font-weight:600}.radar-empty{fill:#596c64;font-size:14px}</style>${referenceBand}<g class="radar-grid">${grid}${spokes}</g>${referenceAxes}${polygon}${segments}${dots}${labels}${valid.length?'':'<text class="radar-empty" x="230" y="191" text-anchor="middle">미평가</text>'}</svg>`;
 const clipped=axes.some(a=>a.clipped);
 return`<figure class="metric-radar-card" data-radar-kind="${kind}"><figcaption><strong>${title}</strong><span>${valid.length}/6개 지표</span></figcaption><p class="metric-radar-source">${esc(source)}</p><div class="metric-radar-legend">${kind==='properties'?`<span class="metric-radar-legend-item"><span class="metric-radar-swatch" aria-hidden="true"></span><span class="column-heading">참고 범위 · 주황색 영역·점선${helpButton('skin_radar_reference','참고 범위')}</span></span>`:''}<span class="metric-radar-legend-item"><span class="metric-radar-candidate-swatch" aria-hidden="true"></span><span>후보값 · 초록색 실선·점</span></span></div>${svg}${kind==='properties'?`<p class="metric-radar-reading">${esc(referenceNote)}</p>`:''}<p class="metric-radar-reading">${esc(direction)}</p><details class="metric-radar-values"><summary>원값·단위·표시 범위</summary><p>육각형의 축별 표시 범위를 고정해 후보끼리 같은 척도로 볼 수 있습니다. 표시 범위는 그림의 축 설정입니다.${clipped?' †는 표시 범위 밖 값이며, 원값은 아래 표에 표시합니다.':''}</p><table><thead><tr><th>항목</th><th>원값</th><th>단위</th><th>축 표시 범위</th>${kind==='properties'?'<th>피부 참고 범위</th>':''}</tr></thead><tbody>${axes.map(a=>`<tr><th scope="row"><span class="column-heading">${esc(a.label)}${helpButton(a.help,a.label)}</span></th><td title="${esc(a.value??a.reason)}">${esc(a.value===null?a.reason:display(a.value))}</td><td>${esc(a.unit)}</td><td>${a.min} ~ ${a.max}</td>${kind==='properties'?'<td>'+esc(SKIN_REFERENCE[a.key]?(SKIN_REFERENCE[a.key].minExclusive?'0 초과':SKIN_REFERENCE[a.key].min)+' ~ '+SKIN_REFERENCE[a.key].max:'미설정')+'</td>':''}</tr>`).join('')}</tbody></table></details><a class="metric-radar-save" download="SkinPipe-${kind}.svg" href="data:image/svg+xml;charset=utf-8,${encodeURIComponent(svg)}">그래프 SVG 저장</a></figure>`;
}

// The retained review accepts only its recorded ADMET-AI channel; future runs
// pass the verified endpoint definitions from their own result records.
export function recordedReviewRadars(row){
 const records=Array.isArray(row?.admet)?row.admet:[];
 const cards=records.filter(p=>/^admet-ai-2\.0\.1-[a-f0-9]{16}$/.test(p?.model_id||'')&&['success','partial_nonfinite'].includes(p.status)).map(p=>{
  const blocked=Object.values(p.quality_flags||{}).flatMap(x=>Array.isArray(x)?x:[]).filter(x=>typeof x==='string');
  const source='저장 모델 · '+p.model_id;
  return metricRadarMarkup({kind:'properties',values:p.values,source,blockedKeys:blocked})+metricRadarMarkup({kind:'toxicity',values:p.values,source,blockedKeys:blocked});
 });
 return cards.length?`<section class="metric-radar-section" aria-label="물성·독성 레이더"><h3>한눈에 보는 물성·독성</h3><div class="metric-radar-grid">${cards.join('')}</div></section>`:'<p class="metric-radar-reading">레이더에 표시할 저장 모델·유효값이 없습니다. 상세 원자료에서 계산 상태를 확인하세요.</p>';
}
