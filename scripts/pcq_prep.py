"""usage: pcq_prep.py <file-index> <shards> [MAXORB=150] [LIMIT=0]
Stream one PubChemQC chon300nosalt JSON file from Hugging Face, keep neutral singlets with
<= MAXORB orbitals (basis-count - heavy-atom-count), and split them round-robin into shard files
cand/shard<k>.jsonl. LIMIT>0 keeps only the first LIMIT candidates per shard (for test runs)."""
import sys, os, json, urllib.request, time
idx=int(sys.argv[1]); S=int(sys.argv[2]); MAXORB=int(sys.argv[3]) if len(sys.argv)>3 else 150; LIMIT=int(sys.argv[4]) if len(sys.argv)>4 else 0
API='https://huggingface.co/api/datasets/molssiai-hub/pubchemqc-b3lyp/tree/main/data/b3lyp_pm6_chon300nosalt/train'
RES='https://huggingface.co/datasets/molssiai-hub/pubchemqc-b3lyp/resolve/main/'
files=sorted(x['path'] for x in json.load(urllib.request.urlopen(API)) if x['path'].endswith('.json'))
f=files[idx]; print('file',idx,f,flush=True)
def log(*a): print(time.strftime('%H:%M:%S'),*a,flush=True)
def lines(url):
    off=0; tries=0
    while True:
        try:
            req=urllib.request.Request(url,headers={'Range':f'bytes={off}-'} if off else {})
            resp=urllib.request.urlopen(req,timeout=600)
            if off and resp.status!=206: raise IOError(f'server ignored the Range request (status {resp.status})')
            cr=resp.headers.get('Content-Range'); cl=resp.headers.get('Content-Length')
            total=int(cr.rsplit('/',1)[1]) if cr else (off+int(cl) if cl else None)
            pend=None
            for raw in resp:
                if pend is not None: off+=len(pend); tries=0; yield pend; pend=None
                if not raw.endswith(b'\n'): pend=raw; continue
                off+=len(raw); tries=0; yield raw
            if pend is not None and total is not None and off+len(pend)==total: yield pend+b'\n'; return   # file's last line has no newline
            if pend is None and (total is None or off==total): return
            raise IOError(f'short read: {off} of {total} bytes')   # connection closed early; reconnect from off
        except Exception as e:
            tries+=1; log('reconnect',off,repr(e)[:120])
            if tries>20: raise
            time.sleep(min(60,5*tries))
os.makedirs('cand',exist_ok=True); outs=[open(f'cand/shard{k}.jsonl','w') for k in range(S)]; n=[0]*S
scanned=small=0; buf=None
for raw in lines(RES+f):
    line=raw.decode().rstrip('\n')
    if line=='    {': buf=['{']; continue
    if buf is None: continue
    if line in ('    }','    },'):
        buf.append('}'); r=json.loads(''.join(buf)); buf=None; scanned+=1
        if r.get('multiplicity')!=1 or abs(r.get('charge',1))>1e-9: continue
        if r['basis-count']-r['heavy-atom-count']>MAXORB: continue
        k=small%S; small+=1
        if LIMIT and n[k]>=LIMIT:
            if all(x>=LIMIT for x in n): break
            continue
        outs[k].write(json.dumps(r)+'\n'); n[k]+=1
    else: buf.append(line)
for o in outs: o.close()
json.dump(dict(file=f,index=idx,scanned=scanned,small=small,shards=n),open('cand/summary.json','w'))
log('done',json.dumps(dict(scanned=scanned,small=small,shards=n)))
