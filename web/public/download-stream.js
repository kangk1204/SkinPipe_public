// Reserve a user-selected file before waiting for the export job. Without a
// File System Access API, bound the Blob fallback so large exports cannot fill RAM.
function pickerContextAvailable(context) {
 if(context.isSecureContext===false)return false;
 // Sites embeds the app in a cross-origin frame. The picker is forbidden
 // there even when the function exists; avoid opening an unusable dialog.
 // https://wicg.github.io/file-system-access/#api-showsavefilepicker
 try {
  if(context.top&&context.top!==context&&context.top.location.origin!==context.location.origin)return false;
 } catch(error) {
  if(error?.name==='SecurityError')return false;
  throw error;
 }
 return true;
}
export async function reserveDownload(name, picker=globalThis.showSaveFilePicker, context=globalThis) {
 if(typeof picker!=='function'||!pickerContextAvailable(context))return null;
 try {
  const handle=await picker.call(context,{suggestedName:name});
  return await handle.createWritable();
 } catch(error) {
  // An unavailable capability can use the bounded browser download instead.
  // AbortError is the user's cancellation and must not start an export job.
  if(['SecurityError','NotAllowedError','NotSupportedError'].includes(error?.name))return null;
  throw error;
 }
}
export async function downloadResponse(response,name,{writer=null,maxBlobBytes=256*1024*1024,saveBlob}={}) {
 try {
  if (!response.ok) throw Error(`다운로드 실패 (HTTP ${response.status})`);
  if (writer) {
   if(!response.body)throw Error('응답 스트림이 없습니다.');
   await response.body.pipeTo(writer);
   return {mode:'stream'};
  }
  const length=Number(response.headers.get('Content-Length'));
  if(!Number.isSafeInteger(length)||length<=0||length>maxBlobBytes) {
   await response.body?.cancel();
   throw Error('큰 내보내기는 파일 직접 저장을 지원하는 브라우저에서 내려받아 주세요.');
  }
  const reader=response.body.getReader(),parts=[];let received=0;
  while(true){const {value,done}=await reader.read();if(done)break;received+=value.byteLength;
   if(received>maxBlobBytes){await reader.cancel();throw Error('내보내기 크기가 메모리 다운로드 한도를 초과했습니다.');}
   parts.push(value);
  }
  if(received!==length)throw Error('다운로드 길이가 원본 기록과 다릅니다. 다시 시도하세요.');
  const blob=new Blob(parts,{type:response.headers.get('Content-Type')||'application/octet-stream'});
  if(saveBlob)saveBlob(blob,name);else{const url=URL.createObjectURL(blob),a=document.createElement('a');a.href=url;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(url),30000);}
  return {mode:'blob',bytes:received};
 } catch(error) { if(writer)await writer.abort().catch(()=>{});throw error; }
}
