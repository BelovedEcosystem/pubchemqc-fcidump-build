"""usage: pcq_prep.py <file-index> <shards|auto> [MAXORB=150] [LIMIT=0] [SUBSET=chon300nosalt] [SKIP=] [PER_SHARD=500] [MAX_SHARDS=64]
Stream one PubChemQC JSON file of the given subset (b3lyp_pm6_<SUBSET>) from Hugging Face, keep neutral singlets
with <= MAXORB orbitals (basis-count - heavy-atom-count), drop CIDs listed in SKIP (a gzipped or plain text file,
one CID per line: molecules an earlier build already processed), and split the rest round-robin into shard files
cand/shard<k>.jsonl. shards=auto picks ceil(candidates/PER_SHARD) shards (1..MAX_SHARDS) so a dense file is not cut
off by the per-job time budget. LIMIT>0 keeps only the first LIMIT candidates per shard (for test runs).
cand/summary.json records the counts, including candidates with P, S, F or Cl ("new_element") vs CHON only.
Optional environment (defaults keep the original behaviour):
  PCQ_EXCLUDE_Z   comma-separated atomic numbers; small neutral singlets containing any of them are dropped (counted)
  PCQ_NO_MULTIFRAG=1  drop small neutral singlets whose PubChem SMILES has several fragments ('.'), counted and listed
                  in <out>/multifrag.jsonl (pcq_convert.py would reject them anyway: "multiple fragments")
  PCQ_FILES       JSON file with the Hugging Face tree listing (avoids one API call per file)
  PCQ_OUTDIR      output directory (default cand)"""
import sys, os, json, urllib.request, time, gzip, math
A=sys.argv+['']*8
idx=int(A[1]); S_ARG=A[2]; MAXORB=int(A[3] or 150); LIMIT=int(A[4] or 0); SUBSET=A[5] or 'chon300nosalt'; SKIP=A[6]
PER=int(A[7] or 500); MAXS=int(A[8] or 64)
skip=set()
if SKIP:
    with (gzip.open(SKIP,'rt') if SKIP.endswith('.gz') else open(SKIP)) as fh: skip={int(l) for l in fh if l.strip()}
NEWZ={9,15,16,17}
EXZ={int(z) for z in os.environ.get('PCQ_EXCLUDE_Z','').split(',') if z.strip()}
NOMF=os.environ.get('PCQ_NO_MULTIFRAG')=='1'; OUT=os.environ.get('PCQ_OUTDIR','cand'); FL=os.environ.get('PCQ_FILES')
API=f'https://huggingface.co/api/datasets/molssiai-hub/pubchemqc-b3lyp/tree/main/data/b3lyp_pm6_{SUBSET}/train'
RES='https://huggingface.co/datasets/molssiai-hub/pubchemqc-b3lyp/resolve/main/'
files=sorted(x['path'] for x in json.load(open(FL) if FL else urllib.request.urlopen(API)) if x['path'].endswith('.json'))
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
os.makedirs(OUT,exist_ok=True); keep=[]; mf=[]
scanned=small=skipped=0; buf=None; small_all=metal=multifrag=0
for raw in lines(RES+f):
    line=raw.decode().rstrip('\n')
    if line=='    {': buf=['{']; continue
    if buf is None: continue
    if line in ('    }','    },'):
        buf.append('}'); r=json.loads(''.join(buf)); buf=None; scanned+=1
        if r.get('multiplicity')!=1 or abs(r.get('charge',1))>1e-9: continue
        if r['basis-count']-r['heavy-atom-count']>MAXORB: continue
        small_all+=1
        if EXZ and set(r['atomic-numbers'])&EXZ: metal+=1; continue
        if NOMF and '.' in r.get('pubchem-isomeric-smiles',''):
            multifrag+=1; mf.append(dict(cid=r['cid'],formula=r['formula'],smiles=r.get('pubchem-isomeric-smiles',''),already_done=r['cid'] in skip)); continue
        small+=1
        if r['cid'] in skip: skipped+=1; continue
        keep.append(json.dumps(r))
    else: buf.append(line)
S=max(1,min(MAXS,math.ceil(len(keep)/PER))) if S_ARG=='auto' else int(S_ARG)
outs=[open(f'{OUT}/shard{k}.jsonl','w') for k in range(S)]; n=[0]*S; newel=0
for i,j in enumerate(keep):
    k=i%S
    if LIMIT and n[k]>=LIMIT: continue
    outs[k].write(j+'\n'); n[k]+=1
    newel+=bool(set(json.loads(j)['atomic-numbers'])&NEWZ)
for o in outs: o.close()
summ=dict(file=f,index=idx,subset=SUBSET,scanned=scanned,small=small,skipped_already_done=skipped,candidates=sum(n),
          new_element=newel,chon_only=sum(n)-newel,nshards=S,shards=n)
if EXZ or NOMF:
    summ.update(small_neutral_singlets_all=small_all,excluded_z=sorted(EXZ),excluded_element=metal,excluded_multifragment=multifrag)
    with open(f'{OUT}/multifrag.jsonl','w') as fh: fh.writelines(json.dumps(x)+'\n' for x in mf)
json.dump(summ,open(f'{OUT}/summary.json','w'))
log('done',json.dumps(summ))
