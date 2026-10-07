"""usage: pcq_pool.py <prep-artifacts-dir> <nshards> <outdir>
Pooled mode: gather every per-file candidate list written by pcq_prep_group.py (<dir>/cand-*/f<idx>/shard0.jsonl),
drop repeated CIDs (first file wins), and split the pool round-robin into <nshards> shard files <outdir>/shard<k>.jsonl.
<outdir>/summary.json sums the per-file counts and lists every file's summary plus any file whose prep failed or
whose summary is missing, so the scan's coverage can be checked."""
import sys, os, json, glob
src, S, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]; os.makedirs(out, exist_ok=True)
files = {}; failed = []; listed = set()
for g in sorted(glob.glob(f'{src}/cand-*/group.json')):
    d = json.load(open(g)); failed += d['failed']; listed |= set(range(d['first'], d['last'] + 1))
for p in glob.glob(f'{src}/cand-*/f*/summary.json'):
    d = json.load(open(p)); files[d['index']] = (d, os.path.dirname(p))
missing = sorted(listed - set(files))
seen = set(); rows = []; dup = 0; mf = []
for i in sorted(files):
    d, dr = files[i]
    for l in open(f'{dr}/shard0.jsonl'):
        if not l.strip(): continue
        c = json.loads(l)['cid']
        if c in seen: dup += 1; continue
        seen.add(c); rows.append(l if l.endswith('\n') else l + '\n')
    if os.path.exists(f'{dr}/multifrag.jsonl'): mf += [dict(json.loads(x), file=i) for x in open(f'{dr}/multifrag.jsonl') if x.strip()]
S = max(1, min(S, len(rows))) if rows else 1
outs = [open(f'{out}/shard{k}.jsonl', 'w') for k in range(S)]; n = [0] * S
for j, l in enumerate(rows): outs[j % S].write(l); n[j % S] += 1
for o in outs: o.close()
with open(f'{out}/multifrag.jsonl', 'w') as fh: fh.writelines(json.dumps(x) + '\n' for x in mf)
key = lambda k: sum(files[i][0].get(k, 0) for i in files)
summ = dict(index='pool', subset=next(iter(files.values()))[0]['subset'] if files else None, files_scanned=len(files),
            files_listed=len(listed), files_failed=sorted(set(failed)), files_missing=missing,
            scanned=key('scanned'), small_neutral_singlets_all=key('small_neutral_singlets_all'), excluded_element=key('excluded_element'),
            excluded_multifragment=key('excluded_multifragment'), small=key('small'), skipped_already_done=key('skipped_already_done'),
            candidates_before_dedup=key('candidates'), duplicate_cids=dup, candidates=len(rows), new_element=None, nshards=S, shards=n,
            per_file={i: {k: files[i][0].get(k) for k in ('file', 'scanned', 'small_neutral_singlets_all', 'excluded_element',
                      'excluded_multifragment', 'small', 'skipped_already_done', 'candidates')} for i in sorted(files)})
summ['new_element'] = sum(bool(set(json.loads(l)['atomic-numbers']) & {9, 15, 16, 17}) for l in rows)
json.dump(summ, open(f'{out}/summary.json', 'w'), indent=1)
print(json.dumps({k: v for k, v in summ.items() if k != 'per_file'}))
