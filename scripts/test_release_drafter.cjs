// Offline safety regression for the actual inline workflow; no GitHub writes.
const fs=require('fs'),crypto=require('crypto');
const path=require('path');
const workflow=fs.readFileSync(path.join(__dirname,'../.github/workflows/release-drafter.yml'),'utf8');
const scripts=[...workflow.matchAll(/^          script: \|\n((?:(?: {12}.*)?\n)*)/gm)];
if(!scripts.length)throw Error('No inline workflow script found');
const code=scripts.at(-1)[1].split('\n').map(line=>line.replace(/^ {12}/,'')).join('\n');
const run=new (Object.getPrototypeOf(async function(){}).constructor)('github','context','core','require','process',code);
const sha='e1c13fc82a9077720a6b27350a5729993e772a83';
const original='Curated overview\n<!-- release-drafter:start -->\nold\n<!-- release-drafter:end -->\nAssets and checksums unchanged';
async function test(kind){let writes=[];let body=original;
const env={PROPOSED_NOTES:'<!-- release-drafter:start -->\n- Included (#1)\n- After tag (#2)\n<!-- release-drafter:end -->',RELEASE_ID:'123',SNAPSHOT_UPDATED_AT:'T',SNAPSHOT_BODY_HASH:crypto.createHash('sha256').update(body).digest('hex')};
if(kind==='malformed')env.PROPOSED_NOTES='<!-- release-drafter:start -->\n- malformed entry\n<!-- release-drafter:end -->';
if(kind==='markers')env.PROPOSED_NOTES+='<!-- release-drafter:start -->';
const github={rest:{repos:{getCommit:async()=>({data:{sha:kind==='moved-tag'?'changed':sha}}),compareCommitsWithBasehead:async({basehead})=>({data:{status:basehead.startsWith('merge1')?'ahead':'behind'}}),getRelease:async()=>({data:{draft:kind!=='published',tag_name:'v0.1.0a0',updated_at:kind==='changed'?'NEW':'T',body}}),updateRelease:async x=>{writes.push(x);}},pulls:{get:async({pull_number:n})=>({data:{merged:true,merge_commit_sha:'merge'+n,title:kind==='reserved-title'?'<!-- release-drafter:start -->':n===1?'Included':'After tag'}})}}};
const core={summary:{addRaw(){return this},async write(){}}};
let error=null;try{await run(github,{repo:{owner:'test',repo:'test'}},core,require,{env});}catch(e){error=e;}
if(kind==='valid'){if(error)throw error;if(writes.length!==1||writes[0].body.includes('#2')||!writes[0].body.startsWith('Curated overview\n<!--')||!writes[0].body.endsWith('Assets and checksums unchanged'))throw Error('preservation failed');if(Object.keys(writes[0]).sort().join(',')!=='body,owner,release_id,repo')throw Error('unexpected update fields');}
else if(!error||writes.length)throw Error('Fail-closed test failed '+kind);
console.log(kind,'PASS');}
(async()=>{for(const kind of ['valid','changed','published','malformed','markers','moved-tag','reserved-title'])await test(kind)})().catch(e=>{console.error(e);process.exitCode=1});
