"""usage: pcq_merge.py <dest> <shard-out-dir>...
Merge the out/ directories of several build shards (each with catalog.jsonl and fcidump/<cid>/) into one
directory laid out the way pcq_finalize.py expects: <dest>/catalog.jsonl and <dest>/fcidump/<cid>/.
Idempotent: a CID already in <dest>/catalog.jsonl is skipped. Where two shards both report the same CID, the
first one seen wins and the duplicate is listed in <dest>/merge_duplicates.jsonl."""
import sys, os, json, shutil
dest = sys.argv[1]; srcs = sys.argv[2:]
if not srcs: raise SystemExit(__doc__)
os.makedirs(f'{dest}/fcidump', exist_ok=True)
cat = f'{dest}/catalog.jsonl'
seen = {json.loads(l)['cid'] for l in open(cat)} if os.path.exists(cat) else set()
before = set(seen)   # merged by an earlier call: skipped quietly, not counted as duplicates
added = dup = missing = already = 0
with open(cat, 'a') as out, open(f'{dest}/merge_duplicates.jsonl', 'a') as dups:
    for src in srcs:
        p = f'{src}/catalog.jsonl'
        if not os.path.exists(p): print('no catalog in', src, flush=True); continue
        for l in open(p):
            if not l.strip(): continue
            row = json.loads(l); cid = row['cid']
            if cid in before: already += 1; continue
            if cid in seen: dup += 1; dups.write(json.dumps(dict(cid=cid, src=src, status=row.get('status')))+'\n'); continue
            if row.get('status') == 'ok':
                d = f'{src}/fcidump/{cid:09d}'
                if not os.path.isdir(d): missing += 1; print('ok row without files, skipped:', cid, src, flush=True); continue
                shutil.copytree(d, f'{dest}/fcidump/{cid:09d}', dirs_exist_ok=True)   # a copy interrupted before its catalog line was written is redone
            out.write(json.dumps(row)+'\n'); seen.add(cid); added += 1
print(json.dumps(dict(dest=dest, shards=len(srcs), added=added, already_merged=already, duplicates=dup, ok_rows_missing_files=missing, total=len(seen))), flush=True)
