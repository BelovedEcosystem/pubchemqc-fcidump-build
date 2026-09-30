"""usage: pcq_run.py <records.jsonl> [workers]  one PubChemQC record per line. Resumable via catalog.jsonl."""
import sys, os, json, subprocess, time, hashlib, tempfile, shutil, concurrent.futures as cf, threading
HERE=os.path.dirname(os.path.abspath(__file__)); D=os.environ.get('PCQ_OUT','out'); PY=sys.executable; CHK=f'{HERE}/fcidump_check.py'; lock=threading.Lock()
DEADLINE=float(os.environ.get('DEADLINE','0')); os.makedirs(f'{D}/fcidump',exist_ok=True)
recs=[json.loads(l) for l in open(sys.argv[1]) if l.strip()]; W=int(sys.argv[2]) if len(sys.argv)>2 else 2
cat=f'{D}/catalog.jsonl'; done={json.loads(l)['cid'] for l in open(cat)} if os.path.exists(cat) else set()
def sha(p): return hashlib.sha256(open(p,'rb').read()).hexdigest()
def one(rec):
    with lock:   # claim the CID so a repeated record in the input cannot run twice in the same directory
        if rec['cid'] in done: return
        done.add(rec['cid'])
    if DEADLINE and time.time()>DEADLINE: return
    t0=time.time(); out=f'{D}/fcidump/{rec["cid"]:09d}'; row=dict(cid=rec['cid'],formula=rec['formula'],natoms=rec['atom-count'])
    tf=tempfile.NamedTemporaryFile('w',suffix='.json',delete=False); json.dump(rec,tf); tf.close()
    try:
        r=subprocess.run([PY,f'{HERE}/pcq_convert.py',tf.name,out],capture_output=True,text=True,timeout=7200,
                         env={**os.environ,'OMP_NUM_THREADS':str(max(1,os.cpu_count()//W))})
        last=(r.stdout.strip().splitlines() or [''])[-1] if r.returncode==0 else ((r.stderr.strip().splitlines() or ['?'])[-1])
        if r.returncode:
            row.update(status='rejected' if last.startswith('REJECT') else 'error',reason=last[:300]); raise StopIteration
        row.update(json.load(open(f'{out}/meta.json'))); row.pop('geometry',None); checks={}
        for f in sorted(os.listdir(out)):
            if f.startswith('FCIDUMP_'):
                j=f'{out}/check_{f[:-3]}.json'
                subprocess.run([PY,CHK,f'{out}/{f}','--json',j]+(['--no-mf'] if 'full' in f else []),capture_output=True,timeout=3600)
                c=json.load(open(j)) if os.path.exists(j) else {}
                checks[f]=dict(verdict=c.get('verdict') or c.get('status'),sha256=c.get('info',{}).get('sha256'),sha256_gz=sha(f'{out}/{f}'),bytes=os.path.getsize(f'{out}/{f}'))  # sha256 = uncompressed text
        bad=[f for f,c in checks.items() if c['verdict']!='ACCEPTED']
        if bad: row.update(status='rejected',reason=f'checker did not accept {bad}',checks=checks); raise StopIteration
        row.update(checks=checks,status='ok')
    except StopIteration: pass
    except Exception as e: row.update(status='error',reason=str(e)[:300])
    finally: os.unlink(tf.name)
    if row['status']!='ok': shutil.rmtree(out,ignore_errors=True)
    row['seconds']=round(time.time()-t0,1)
    with lock:
        open(cat,'a').write(json.dumps(row)+'\n')
        print(time.strftime('%H:%M:%S'),row['status'],row['cid'],row['formula'],row.get('norb'),row['seconds'],row.get('reason',''),'' if row.get('stable_external',True) else 'UHF-unstable',flush=True)
with cf.ThreadPoolExecutor(W) as ex: list(ex.map(one,recs))
print('finished',flush=True)
