"""usage: pcq_package.py <out-dir> <library> [<library> ...] [--entry-numbers CSV] [--cids FILE] [--tar]

Package one or more build libraries (the output of pcq_finalize.py: <library>/catalog.jsonl and
<library>/published/<CID, 9 digits>/) into the deposited layout:

  <out>/become-fcidump-pcq/fcidump/pcq<CID>_6-31g_cas8e8o.fcidump   active-space FCIDUMP, uncompressed
  <out>/become-fcidump-pcq/xyz/pcq<CID>.xyz                          PM6 geometry
  <out>/catalog_index.csv                                            one row per entry (columns below)
  <out>/SHA256SUMS.txt                                               SHA-256 of every packaged file
  <out>/become-fcidump-pcq.tar.gz                                    with --tar: reproducible archive of the folder

Only strict-tier entries whose active-space file has the checker verdict ACCEPTED are packaged. Libraries are read in
the order given; a CID that appears in more than one library is taken from the first. Every file is checked against
the hashes recorded by pcq_finalize.py before it is written. --cids restricts the package to the CIDs listed in FILE
(one per line). Column i is the catalog entry number: taken from --entry-numbers (a CSV with columns cid and i) when
given, otherwise numbered 1, 2, ... in the order the entries are packaged. No correlated energy is read or written.

catalog_index.csv columns: i, cid, source, formula, natoms, charge, multiplicity, basis, nbasis, norb, nelec, ms2,
e_rhf_Ha, pubchemqc_b3lyp_eV, rhf_homo_lumo_gap_Ha, cas_edge_degenerate, fcidump_path, xyz_path, bytes, sha256."""
import argparse, csv, gzip, hashlib, io, json, os, tarfile

EV = 27.211386245988          # eV per Eh, used to report the PubChemQC B3LYP energy in eV
TOP = 'become-fcidump-pcq'
MTIME = 1791504000            # fixed member mtime so the archive is reproducible
COLS = ['i', 'cid', 'source', 'formula', 'natoms', 'charge', 'multiplicity', 'basis', 'nbasis', 'norb', 'nelec', 'ms2',
        'e_rhf_Ha', 'pubchemqc_b3lyp_eV', 'rhf_homo_lumo_gap_Ha', 'cas_edge_degenerate', 'fcidump_path', 'xyz_path',
        'bytes', 'sha256']

ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
ap.add_argument('out'); ap.add_argument('libraries', nargs='+')
ap.add_argument('--entry-numbers'); ap.add_argument('--cids'); ap.add_argument('--tar', action='store_true')
a = ap.parse_args()
sha = lambda b: hashlib.sha256(b).hexdigest()
want = {int(l) for l in open(a.cids) if l.strip()} if a.cids else None
num = {int(r['cid']): int(r['i']) for r in csv.DictReader(open(a.entry_numbers))} if a.entry_numbers else None

rows, files, seen = [], {}, set()
for lib in a.libraries:
    for line in open(f'{lib}/catalog.jsonl'):
        r = json.loads(line); cid = r['cid']
        if (want is not None and cid not in want) or cid in seen: continue
        gz, xyz = f'FCIDUMP_pcq{cid}_cas8e8o.gz', f'pcq{cid}.xyz'
        f = r['files'][gz]
        if r.get('tier') != 'strict' or f['checker'] != 'ACCEPTED': continue
        src = f'{lib}/published/{cid:09d}'
        gzb = open(f'{src}/{gz}', 'rb').read(); raw = gzip.decompress(gzb); xb = open(f'{src}/{xyz}', 'rb').read()
        assert sha(gzb) == f['sha256_gz'] and sha(raw) == f['sha256'] and len(raw) == f['bytes_uncompressed'], cid
        assert sha(xb) == r['files'][xyz]['sha256'], cid
        fp, xp = f'fcidump/pcq{cid}_6-31g_cas8e8o.fcidump', f'xyz/{xyz}'
        files[fp], files[xp] = raw, xb; seen.add(cid)
        rows.append(dict(cid=cid, source='pubchemqc', formula=r['formula'], natoms=r['natoms'], charge=0, multiplicity=1,
                         basis=r['basis'], nbasis=r['norb'], norb=f['norb'], nelec=f['nelec'], ms2=f['ms2'],
                         e_rhf_Ha=round(r['e_rhf'], 10), pubchemqc_b3lyp_eV=round(r['e_b3lyp_pubchemqc'] * EV, 6),
                         rhf_homo_lumo_gap_Ha=round(r['homo_lumo_gap'], 10), cas_edge_degenerate=bool(r['cas_edge_degenerate']),
                         fcidump_path=fp, xyz_path=xp, bytes=len(raw), sha256=f['sha256']))
for k, r in enumerate(rows, 1): r['i'] = num[r['cid']] if num is not None else k
rows.sort(key=lambda r: r['i'])
assert len({r['i'] for r in rows}) == len(rows), 'repeated entry numbers'

os.makedirs(a.out, exist_ok=True)
for k, b in sorted(files.items()):
    p = f'{a.out}/{TOP}/{k}'; os.makedirs(os.path.dirname(p), exist_ok=True); open(p, 'wb').write(b)
with open(f'{a.out}/SHA256SUMS.txt', 'w') as fh:
    for k, b in sorted(files.items()): fh.write(f'{sha(b)}  {TOP}/{k}\n')
with open(f'{a.out}/catalog_index.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=COLS, lineterminator='\n'); w.writeheader(); w.writerows(rows)
if a.tar:   # sorted names, USTAR, fixed mtime, uid/gid 0, mode 0644, gzip mtime 0
    with open(f'{a.out}/{TOP}.tar.gz', 'wb') as raw:
        with gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0, compresslevel=9) as g:
            with tarfile.open(fileobj=g, mode='w', format=tarfile.USTAR_FORMAT) as tf:
                for k, b in sorted(files.items()):
                    ti = tarfile.TarInfo(f'{TOP}/{k}'); ti.size = len(b); ti.mtime = MTIME; ti.mode = 0o644
                    ti.uid = ti.gid = 0; ti.uname = ti.gname = ''; tf.addfile(ti, io.BytesIO(b))
print(json.dumps(dict(entries=len(rows), files=len(files))))
