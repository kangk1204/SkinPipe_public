import assert from 'node:assert/strict';
import {renderDiscovery} from '../web/public/discovery.js';

// A pending request must show the new screen; its late result must not replace
// a different route selected in the meantime.
let finish, current=true, calls=0;
const main={innerHTML:'previous candidate screen'};
const pending=renderDiscovery(main,()=>{calls++;return new Promise(resolve=>{finish=resolve})},'researcher',{isCurrent:()=>current});
assert.equal(calls,1);
assert.match(main.innerHTML,/role="status"/);
assert.doesNotMatch(main.innerHTML,/previous candidate screen/);
current=false;main.innerHTML='newly selected screen';finish({snapshots:[]});
await pending;
assert.equal(main.innerHTML,'newly selected screen');
await renderDiscovery(main,()=>{throw Error('stale route should not request data')},'researcher',{isCurrent:()=>false});
assert.equal(main.innerHTML,'newly selected screen');
console.log('Discovery loading state and late route response verified');
