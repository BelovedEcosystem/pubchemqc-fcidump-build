"""Build the strict-tier published library from catalog.jsonl. Idempotent; rerun any time.
Strict tier = status ok AND RHF stable toward UHF (stable_external True).
Every published FCIDUMP is re-read from disk and must reproduce e_rhf to 1e-8 Eh and match both recorded hashes
(sha256 of the uncompressed text, sha256_gz of the .gz archive). Recipe version 2.
usage: PCQ_DIR=<dir with catalog.jsonl and fcidump/> [PCQ_SUBSET=chon300nosalt] pcq_finalize.py"""
import json, os, sys, gzip, shutil, hashlib, csv, platform, tempfile
import numpy as np, scipy, pyscf
from pyscf import ao2mo
from pyscf.tools import fcidump
D = os.environ.get('PCQ_DIR', '/workspace/pubchemqc_v2'); SUBSET = os.environ.get('PCQ_SUBSET', 'chon300nosalt'); L = f'{D}/library'; PUB = f'{L}/published'
SYM = {1:'H',5:'B',6:'C',7:'N',8:'O',9:'F',14:'Si',15:'P',16:'S',17:'Cl',35:'Br'}
TOL = 1e-8
sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()
versions = dict(python=platform.python_version(), pyscf=pyscf.__version__, numpy=np.__version__, scipy=scipy.__version__)
RECIPE = dict(
    recipe_version=2,
    source='PubChemQC B3LYP/6-31G*//PM6 (Hugging Face molssiai-hub/pubchemqc-b3lyp, config ' + SUBSET + '), CC-BY-4.0',
    geometry='PubChemQC PM6-optimized geometry (the B3LYP/6-31G*//PM6 set: B3LYP single points on PM6 structures), used unchanged (Angstrom)',
    structure_check='PySCF B3LYP/6-31G* with xc="b3lypg" (B3LYP with VWN-RPA correlation; the variant that reproduces the PubChemQC energies), Cartesian 6D, grid level 3, conv_tol 1e-9, must match PubChemQC total energy within 1e-4 Eh',
    mean_field='RHF, 6-31G* (spherical), conv_tol 1e-10, must converge; internally stable; stable toward UHF (strict tier)',
    orbitals='canonical RHF orbitals; orthonormality error < 1e-8; HOMO-LUMO gap > 0; degenerate orbitals (same occupation, energies within 1e-5 Eh) fixed by diagonalizing a fixed generic operator inside each block (the Coulomb potential of unit point charges at the atom centroid plus (1.37, 0.71, 0.43) and plus (-0.52, 1.19, -0.88) Bohr, weights 1 and 0.61; eigenvalue order) (cas_edge_degenerate marks molecules where the 8-orbital window cuts through such a block); phase rule: in each orbital the largest-magnitude AO coefficient (the first one within 1e-6 of the maximum) is made positive',
    active_space='8 electrons in 8 orbitals (4 highest occupied + 4 lowest virtual canonical orbitals); inactive core folded into ECORE and 1e integrals',
    full_file='written when the basis has <= 60 orbitals (all orbitals, no frozen core)',
    verification='each file re-read from disk reproduces its one- and two-electron integrals and core energy within 1e-10 and the RHF energy within 1e-8 Eh; CASCI(8e,8o) from the file matches a direct PySCF CASCI within 1e-8 Eh; independent fcidump_check.py verdict ACCEPTED',
    hashes='sha256 = SHA-256 of the uncompressed FCIDUMP text (the identifier; listed in SHA256SUMS against the uncompressed filename, so gunzip then sha256sum matches); sha256_gz = SHA-256 of the .gz archive (listed in SHA256SUMS-gz); archives are written with gzip mtime 0',
    reproducibility='files are not guaranteed byte-identical between runs, even on the same machine and software (floating-point summation order varies at the last digit); repeat runs, including different thread counts, agree within about 1e-9 in every integral, so equivalence is by integral tolerance, not by hash',
    software=versions)

def det_energy(path):
    tmp = None
    if path.endswith('.gz'):
        tmp = tempfile.NamedTemporaryFile(delete=False); tmp.write(gzip.open(path).read()); tmp.close(); p = tmp.name
    else: p = path
    try:
        d = fcidump.read(p, verbose=False)
    finally:
        if tmp: os.unlink(tmp.name)
    n, ne = d['NORB'], d['NELEC']; no = ne // 2
    h1 = np.asarray(d['H1']); g = ao2mo.restore(1, d['H2'], n)
    o = np.arange(no)
    e = d['ECORE'] + 2 * h1[o, o].sum() + 2 * np.einsum('iijj->', g[:no, :no, :no, :no]) - np.einsum('ijji->', g[:no, :no, :no, :no])
    return float(e), n, ne, d.get('MS2', 0)

rows = {}
for l in open(f'{D}/catalog.jsonl'):
    r = json.loads(l); rows[r['cid']] = r
os.makedirs(PUB, exist_ok=True)
published, excluded, problems = [], [], []
for cid, r in sorted(rows.items()):
    if r['status'] != 'ok':
        excluded.append(dict(cid=cid, formula=r['formula'], natoms=r['natoms'], stage='pipeline', reason=r.get('reason', '')[:200])); continue
    if r.get('stable_external') is not True:
        excluded.append(dict(cid=cid, formula=r['formula'], natoms=r['natoms'], stage='strict tier', reason='RHF unstable toward UHF (open-shell instability); file exact but excluded from strictest tier'))
        shutil.rmtree(f'{PUB}/{cid:09d}', ignore_errors=True); continue
    src = f'{D}/fcidump/{cid:09d}'; dst = f'{PUB}/{cid:09d}'
    meta = json.load(open(f'{src}/meta.json'))
    os.makedirs(dst, exist_ok=True)
    files = {}
    ok = True
    for f in sorted(os.listdir(src)):
        if not f.startswith('FCIDUMP_'): continue
        rec = r['checks'][f]; raw = gzip.open(f'{src}/{f}').read() if f.endswith('.gz') else open(f'{src}/{f}', 'rb').read()
        s = hashlib.sha256(raw).hexdigest(); sgz = sha(f'{src}/{f}')
        if s != rec.get('sha256'): problems.append(f'{cid} {f} text sha mismatch'); ok = False; continue
        if sgz != rec.get('sha256_gz'): problems.append(f'{cid} {f} archive sha mismatch'); ok = False; continue
        e, n, ne, ms2 = det_energy(f'{src}/{f}')
        diff = e - r['e_rhf']
        if abs(diff) > TOL or ms2 != 0: problems.append(f'{cid} {f} reread diff {diff:.2e} ms2 {ms2}'); ok = False; continue
        shutil.copy2(f'{src}/{f}', f'{dst}/{f}')
        shutil.copy2(f'{src}/check_{f[:-3]}.json', f'{dst}/check_{f[:-3]}.json')
        files[f] = dict(sha256=s, sha256_gz=sgz, bytes_uncompressed=len(raw), bytes=os.path.getsize(f'{dst}/{f}'), norb=n, nelec=ne, ms2=ms2,
                        checker=r['checks'][f]['verdict'], e_rhf_reread=e, reread_diff=diff)
    if not ok:
        shutil.rmtree(dst, ignore_errors=True)
        excluded.append(dict(cid=cid, formula=r['formula'], natoms=r['natoms'], stage='final verification', reason='; '.join(p for p in problems if p.startswith(f'{cid} ')))); continue
    geo = meta['geometry']
    xyz = f'{dst}/pcq{cid}.xyz'
    with open(xyz, 'w') as fh:
        fh.write(f'{len(geo)}\nPubChem CID {cid} {meta["formula"]} | PubChemQC PM6-optimized geometry (CC-BY-4.0) | Angstrom\n')
        for z, x, y, zz in geo: fh.write(f'{SYM[z]:<2} {x:16.10f} {y:16.10f} {zz:16.10f}\n')
    files[os.path.basename(xyz)] = dict(sha256=sha(xyz), bytes=os.path.getsize(xyz))
    entry = {k: meta.get(k) if k == 'dev_ecore_full' else meta[k] for k in ['cid', 'formula', 'smiles', 'inchi', 'natoms', 'nelec', 'basis', 'norb', 'e_b3lyp_ours', 'e_b3lyp_pubchemqc', 'b3lyp_diff',
                                  'e_rhf', 'e_casci_8e8o', 'e_casci_direct', 'homo_lumo_gap', 'orth_err', 'stable_internal', 'stable_external',
                                  'recipe_version', 'xc', 'degenerate_blocks', 'cas_edge_degenerate', 'dev_h1_cas', 'dev_eri_cas', 'dev_ecore_cas', 'dev_h1_full', 'dev_eri_full', 'dev_ecore_full']}
    entry.update(pubchem_url=f'https://pubchem.ncbi.nlm.nih.gov/compound/{cid}', tier='strict', files=files, recipe=RECIPE)
    json.dump(entry, open(f'{dst}/entry.json', 'w'), indent=1)
    with open(f'{dst}/SHA256SUMS', 'w') as fh:
        for f, v in sorted(files.items()): fh.write(f'{v["sha256"]}  {f[:-3] if f.endswith(".gz") else f}\n')
    with open(f'{dst}/SHA256SUMS-gz', 'w') as fh:
        for f, v in sorted(files.items()):
            if 'sha256_gz' in v: fh.write(f'{v["sha256_gz"]}  {f}\n')
    published.append(entry)

# a directory published by an earlier run whose CID is no longer published (status changed, or gone from the catalog) is removed
for d in os.listdir(PUB):
    if d.isdigit() and int(d) not in {e['cid'] for e in published}: shutil.rmtree(f'{PUB}/{d}', ignore_errors=True)
with open(f'{L}/catalog.jsonl', 'w') as fh:
    for e in published: fh.write(json.dumps({k: v for k, v in e.items() if k != 'recipe'}) + '\n')
cols = ['cid', 'formula', 'natoms', 'nelec', 'norb', 'e_rhf', 'e_casci_8e8o', 'homo_lumo_gap', 'b3lyp_diff', 'smiles', 'inchi', 'pubchem_url']
with open(f'{L}/catalog.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, cols, extrasaction='ignore'); w.writeheader(); [w.writerow(e) for e in published]
with open(f'{L}/excluded.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, ['cid', 'formula', 'natoms', 'stage', 'reason']); w.writeheader(); [w.writerow(e) for e in excluded]
json.dump(RECIPE, open(f'{L}/RECIPE.json', 'w'), indent=1)
summ = dict(candidates=len(rows), published=len(published), excluded=len(excluded),
            excluded_by_stage={s: sum(1 for e in excluded if e['stage'] == s) for s in sorted({e['stage'] for e in excluded})},
            final_verification_problems=problems, max_reread_diff=max((abs(v['reread_diff']) for e in published for v in e['files'].values() if 'reread_diff' in v), default=None),
            full_files=sum(1 for e in published for f in e['files'] if 'full' in f))
json.dump(summ, open(f'{L}/SUMMARY.json', 'w'), indent=1)
print(json.dumps(summ, indent=1))
