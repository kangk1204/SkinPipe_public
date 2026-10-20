// An acceptance reply can be lost after the server has reserved a job.
// Keep its identifier across view recreation until that exact request is acknowledged.
export function pendingSubmission(storage,key,body,uuid=()=>crypto.randomUUID()){
 const encoded=JSON.stringify(body),raw=storage.getItem(key);let pending;
 try{pending=JSON.parse(raw)}catch{}
 if(!pending||pending.body!==encoded||typeof pending.id!=='string'||!/^[a-f0-9]{32}$/.test(pending.id))
  pending={body:encoded,id:uuid().replaceAll('-','')};
 storage.setItem(key,JSON.stringify(pending));return pending;
}
export function acknowledgeSubmission(storage,key,id,activeJobKey){
 const raw=storage.getItem(key);let pending;try{pending=JSON.parse(raw)}catch{}
 if(pending?.id!==id)return false;
 // Persist the accepted job before releasing its retry identifier.
 if(activeJobKey)storage.setItem(activeJobKey,id);
 storage.removeItem(key);return true;
}
