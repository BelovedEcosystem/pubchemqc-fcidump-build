"""usage: pcq_prep_group.py <first-file> <last-file> <MAXORB> <SUBSET> <SKIP> [PARALLEL=4]
Pooled-mode prep: run pcq_prep.py (one shard each) for every source file in [first, last], PARALLEL at a time, into
cand/f<idx>/. The Hugging Face tree listing is fetched once and shared (PCQ_FILES). A file that fails is retried up to
3 times; files that still fail are listed in cand/group.json (the job itself succeeds so the other files' candidates
are kept; re-dispatch the failed indices). Honours PCQ_EXCLUDE_Z and PCQ_NO_MULTIFRAG (see pcq_prep.py)."""
import sys, os, json, subprocess, time, urllib.request, concurrent.futures as cf
a, b, MAXORB, SUBSET, SKIP = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4], sys.argv[5]
P = int(sys.argv[6]) if len(sys.argv) > 6 else 4
HERE = os.path.dirname(os.path.abspath(__file__)); os.makedirs('cand', exist_ok=True)
API = f'https://huggingface.co/api/datasets/molssiai-hub/pubchemqc-b3lyp/tree/main/data/{SUBSET if SUBSET.startswith('b3lyp_pm6') else 'b3lyp_pm6_'+SUBSET}/train'
for t in range(10):
    try: lst = json.load(urllib.request.urlopen(API, timeout=120)); break
    except Exception as e: print('tree listing retry', t, repr(e)[:120], flush=True); time.sleep(10 * (t + 1))
else: raise SystemExit('could not list the source files')
json.dump(lst, open('cand/files.json', 'w'))
env = {**os.environ, 'PCQ_FILES': os.path.abspath('cand/files.json')}
def one(i):
    for t in range(3):
        out = f'cand/f{i}'
        r = subprocess.run([sys.executable, f'{HERE}/pcq_prep.py', str(i), '1', MAXORB, '0', SUBSET, SKIP, '500', '1'],
                           env={**env, 'PCQ_OUTDIR': out}, capture_output=True, text=True)
        print(f'--- file {i} attempt {t} rc {r.returncode}\n{r.stdout[-1500:]}{r.stderr[-1500:]}', flush=True)
        if r.returncode == 0 and os.path.exists(f'{out}/summary.json'): return i, True
        time.sleep(30)
    return i, False
with cf.ThreadPoolExecutor(P) as ex: res = list(ex.map(one, range(a, b + 1)))
json.dump(dict(first=a, last=b, ok=[i for i, k in res if k], failed=[i for i, k in res if not k]), open('cand/group.json', 'w'))
print(open('cand/group.json').read())
