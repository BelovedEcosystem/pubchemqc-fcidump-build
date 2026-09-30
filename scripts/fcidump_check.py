#!/usr/bin/env python3
"""Upload checker for customer FCIDUMP files.

usage: fcidump_check.py FILE [--json OUT] [--no-mf]
Exit code 0 = accepted (PASS or WARN only), 1 = rejected (any FAIL).
Accepts plain, .gz or .zip (first member) files.
"""
import sys, os, re, io, json, gzip, zipfile, hashlib, math, time, argparse
import numpy as np, warnings
warnings.filterwarnings("ignore")

REGISTRY = {
    # sha256 of the uncompressed FCIDUMP text -> known published Hamiltonian
}
REG_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "known_fcidumps.json")
if os.path.exists(REG_PATH):
    REGISTRY.update(json.load(open(REG_PATH)))

MAX_BYTES = 30 * 1024**3          # hard cap on uncompressed size
MF_MAX_NORB = 160                 # skip mean-field check above this

class Report:
    def __init__(s): s.checks = []; s.info = {}
    def add(s, name, status, msg):
        s.checks.append({"check": name, "status": status, "detail": msg})
    @property
    def verdict(s):
        st = [c["status"] for c in s.checks]
        return "REJECTED" if "FAIL" in st else ("ACCEPTED_WITH_WARNINGS" if "WARN" in st else "ACCEPTED")

def open_text(path):
    raw = open(path, "rb").read(4)
    if raw[:2] == b"\x1f\x8b":
        return gzip.open(path, "rb").read(), "gzip"
    if raw[:4] == b"PK\x03\x04":
        z = zipfile.ZipFile(path)
        names = [n for n in z.namelist() if not n.endswith("/")]
        return z.read(names[0]), "zip:" + names[0]
    return open(path, "rb").read(), "plain"

def parse_header(text):
    m = re.search(r"&FCI(.*?)(&END|/\s*\n)", text, re.S | re.I)
    if not m:
        return None, None
    body = m.group(1)
    hdr = {}
    for mm in re.finditer(r"([A-Za-z_][A-Za-z_0-9]*)\s*=\s*([-+0-9.,\sEeDd]*?)(?=[A-Za-z_][A-Za-z_0-9]*\s*=|$)", body, re.S):
        hdr[mm.group(1).upper()] = [v for v in re.split(r"[,\s]+", mm.group(2).strip()) if v]
    for mm in re.finditer(r"\b(UHF|IUHF|TUHF)\b(?!\s*=)", body, re.I):
        hdr[mm.group(1).upper()] = ["1"]
    return hdr, m.end()

def canon2(i, j, k, l):
    if i < j: i, j = j, i
    if k < l: k, l = l, k
    a, b = (i, j), (k, l)
    if a < b: a, b = b, a
    return a + b

def check(path, run_mf=True):
    R = Report(); t0 = time.time()
    size = os.path.getsize(path)
    try:
        data, kind = open_text(path)
    except Exception as e:
        R.add("file", "FAIL", f"Could not open file: {e}"); return R, None
    R.info.update(file=os.path.basename(path), container=kind, bytes_on_disk=size, bytes=len(data),
                  sha256=hashlib.sha256(data).hexdigest())
    if len(data) > MAX_BYTES:
        R.add("file", "FAIL", f"File is {len(data)/1e9:.1f} GB uncompressed; limit is {MAX_BYTES/1e9:.0f} GB."); return R, None
    try:
        text = data.decode("ascii")
    except UnicodeDecodeError:
        R.add("file", "FAIL", "File is not plain ASCII text."); return R, None
    R.add("file", "PASS", f"{len(data)/1e6:.1f} MB of text ({kind}).")

    # ---- header
    hdr, end = parse_header(text)
    if hdr is None:
        R.add("header", "FAIL", "No &FCI ... &END header found."); return R, None
    def geti(k, default=None):
        if k not in hdr: return default
        return int(hdr[k][0])
    try:
        norb, nelec, ms2 = geti("NORB"), geti("NELEC"), geti("MS2", None)
    except ValueError:
        R.add("header", "FAIL", "NORB, NELEC or MS2 is not an integer."); return R, None
    if norb is None or nelec is None:
        R.add("header", "FAIL", "Header must give NORB and NELEC."); return R, None
    if ms2 is None:
        ms2 = nelec % 2
        R.add("header", "WARN", f"No MS2 in header; assuming lowest spin (MS2={ms2}).")
    if any(geti(k, 0) for k in ("UHF", "IUHF")) or "TUHF" in hdr:
        R.add("header", "FAIL", "Header marks separate spin-up and spin-down integrals (UHF format). "
              "Only spin-restricted FCIDUMPs are accepted for now."); return R, None
    orbsym = [int(x) for x in hdr.get("ORBSYM", [])]
    if orbsym and len(orbsym) != norb:
        R.add("header", "WARN", f"ORBSYM has {len(orbsym)} entries but NORB={norb}.")
    na, nb = (nelec + ms2) // 2, (nelec - ms2) // 2
    bad = []
    if norb <= 0: bad.append("NORB must be positive")
    if nelec <= 0: bad.append("NELEC must be positive")
    if (nelec + ms2) % 2: bad.append(f"NELEC={nelec} and MS2={ms2} cannot both be right (one must be even and one odd together)")
    if ms2 < 0 or nb < 0: bad.append(f"MS2={ms2} is larger than NELEC={nelec}")
    if na > norb: bad.append(f"{na} spin-up electrons do not fit in {norb} orbitals")
    if bad:
        R.add("header", "FAIL", "; ".join(bad) + "."); return R, None
    S = ms2 / 2
    R.info.update(norb=norb, nelec=nelec, ms2=ms2, n_alpha=na, n_beta=nb, spin_S=S,
                  multiplicity=ms2 + 1, open_shell=bool(ms2))
    R.add("header", "PASS", f"{nelec} electrons in {norb} orbitals, MS2={ms2} "
          f"({'closed shell' if ms2 == 0 else f'{ms2} unpaired electrons, spin {S:g}'}).")

    # ---- body
    body = text[end:]
    body = body.replace("D", "E").replace("d", "e")
    try:
        arr = np.array(body.split(), dtype=np.float64)
    except ValueError:
        R.add("integrals", "FAIL", "Integral section contains non-numeric text."); return R, None
    if arr.size % 5:
        R.add("integrals", "FAIL", f"Integral section is not whole lines of 5 columns ({arr.size} numbers); file may be truncated."); return R, None
    arr = arr.reshape(-1, 5)
    val, idx = arr[:, 0], arr[:, 1:]
    if not np.all(np.isfinite(val)):
        R.add("integrals", "FAIL", f"{int((~np.isfinite(val)).sum())} integral values are NaN or infinite."); return R, None
    if not np.all(idx == np.round(idx)):
        R.add("integrals", "FAIL", "Some orbital indices are not whole numbers."); return R, None
    idx = idx.astype(np.int64)
    if idx.min() < 0 or idx.max() > norb:
        R.add("integrals", "FAIL", f"Orbital index {int(idx.max()) if idx.max()>norb else int(idx.min())} is outside 1..{norb}."); return R, None
    i, j, k, l = idx.T
    two = (i > 0) & (j > 0) & (k > 0) & (l > 0)
    one = (i > 0) & (j > 0) & (k == 0) & (l == 0)
    core = (i == 0) & (j == 0) & (k == 0) & (l == 0)
    eps = (i > 0) & (j == 0) & (k == 0) & (l == 0)
    other = ~(two | one | core | eps)
    if other.any():
        R.add("integrals", "FAIL", f"{int(other.sum())} lines have an index pattern that is not a valid FCIDUMP entry."); return R, None
    if core.sum() > 1:
        R.add("integrals", "WARN", f"{int(core.sum())} core-energy lines; using the last one.")
    ecore = float(val[core][-1]) if core.any() else 0.0
    if not core.any():
        R.add("integrals", "WARN", "No core-energy line (0 0 0 0); assuming 0.")
    if idx.max() < norb:
        R.add("integrals", "WARN", f"Highest orbital index used is {int(idx.max())}, but NORB={norb}.")

    # one-body
    h1 = np.zeros((norb, norb)); seen1 = {}
    conf1 = 0
    for v, a, b in zip(val[one], i[one] - 1, j[one] - 1):
        key = (max(a, b), min(a, b))
        if key in seen1 and abs(seen1[key] - v) > 1e-8: conf1 += 1
        seen1[key] = v; h1[a, b] = h1[b, a] = v
    # two-body: canonical 8-fold keys, vectorised
    a, b, c, d = i[two] - 1, j[two] - 1, k[two] - 1, l[two] - 1
    ij = np.maximum(a, b) * (np.maximum(a, b) + 1) // 2 + np.minimum(a, b)
    kl = np.maximum(c, d) * (np.maximum(c, d) + 1) // 2 + np.minimum(c, d)
    P, Q = np.maximum(ij, kl), np.minimum(ij, kl)
    key = P * (P + 1) // 2 + Q
    v2 = val[two]
    order = np.argsort(key, kind="stable")
    ks, vs = key[order], v2[order]
    dup = ks[1:] == ks[:-1]
    conf2 = int((dup & (np.abs(vs[1:] - vs[:-1]) > 1e-8)).sum())
    ndup = int(dup.sum())
    npair = norb * (norb + 1) // 2
    eri = np.zeros(npair * (npair + 1) // 2)
    eri[key] = v2
    R.info.update(n_lines=int(arr.shape[0]), n_two_body=int(two.sum()), n_one_body=int(one.sum()),
                  core_energy=ecore, max_abs_two_body=float(np.abs(v2).max()) if v2.size else 0.0)
    if conf1 or conf2:
        R.add("symmetry", "FAIL", f"{conf1 + conf2} integrals appear twice with different values under the "
              "standard 8-fold symmetry of real orbitals; the file may use complex orbitals or be corrupt."); return R, None
    msg = f"{int(two.sum()):,} two-electron and {int(one.sum()):,} one-electron integrals; no conflicting duplicates."
    if ndup: msg += f" ({ndup:,} harmless repeated entries.)"
    R.add("symmetry", "PASS", msg)

    # ---- physics checks on the two-electron integrals
    def g(p, q, r, s):
        pq = max(p, q) * (max(p, q) + 1) // 2 + min(p, q)
        rs = max(r, s) * (max(r, s) + 1) // 2 + min(r, s)
        P, Q = max(pq, rs), min(pq, rs)
        return eri[P * (P + 1) // 2 + Q]
    diag_pairs = np.arange(npair)
    Dpq = eri[diag_pairs * (diag_pairs + 1) // 2 + diag_pairs]      # (pq|pq)
    self_rep = np.array([g(p, p, p, p) for p in range(norb)])
    probs = []
    if (self_rep <= 0).any():
        probs.append(f"{int((self_rep <= 0).sum())} orbitals have zero or negative self-repulsion (ii|ii)")
    if (Dpq < -1e-8).any():
        probs.append(f"{int((Dpq < -1e-8).sum())} exchange-type integrals (pq|pq) are negative")
    # Cauchy-Schwarz: |(pq|rs)| <= sqrt((pq|pq)(rs|rs))
    bound = np.sqrt(np.clip(Dpq[P], 0, None) * np.clip(Dpq[Q], 0, None))
    viol = np.abs(v2) - bound
    worst = float(viol.max()) if viol.size else 0.0
    nviol = int((viol > 1e-6).sum())
    if probs:
        R.add("physics", "FAIL", "; ".join(probs) + ". These cannot come from real molecular orbitals."); return R, None
    if worst > 1e-3:
        R.add("physics", "FAIL", f"{nviol} integrals break the Cauchy-Schwarz bound (worst by {worst:.2e}); "
              "the integrals are not self-consistent."); return R, None
    if nviol:
        R.add("physics", "WARN", f"{nviol} integrals exceed the Cauchy-Schwarz bound by up to {worst:.1e} "
              "(likely rounding in the file).")
    else:
        R.add("physics", "PASS", "All self-repulsion and exchange integrals are positive, and every integral "
              "satisfies the Cauchy-Schwarz bound.")

    # ---- size tier
    ndet = math.comb(norb, na) * math.comb(norb, nb)
    if norb <= 20: tier = "small: exact diagonalisation (FCI) is feasible"
    elif norb <= 40: tier = "medium: selected CI or DMRG"
    elif norb <= 100: tier = "large: DMRG-class solve, custom quote"
    else: tier = "very large: beyond routine DMRG, custom quote"
    R.info.update(n_determinants=f"{ndet:.3e}", size_tier=tier)
    R.add("size", "PASS" if norb <= 100 else "WARN", f"{tier} (about {ndet:.1e} determinants).")

    # ---- known published Hamiltonian?
    known = REGISTRY.get(R.info["sha256"])
    if known:
        R.info["known_hamiltonian"] = known
        R.add("provenance", "PASS", f"Byte-identical to a published Hamiltonian: {known['name']}.")
    else:
        R.add("provenance", "PASS", "Not a known published file; treated as a customer Hamiltonian.")
    R.info["parse_seconds"] = round(time.time() - t0, 1)
    ham = dict(norb=norb, nelec=nelec, ms2=ms2, h1=h1, eri=eri, ecore=ecore)
    if run_mf:
        mean_field(R, ham)
    return R, ham

def mean_field(R, H):
    norb = H["norb"]
    if norb > MF_MAX_NORB:
        R.add("sanity_energy", "WARN", f"Skipped mean-field energy for NORB>{MF_MAX_NORB}."); return
    from pyscf import gto, scf, ao2mo, lib
    lib.num_threads(os.cpu_count())
    t0 = time.time()
    mol = gto.M(verbose=0); mol.nelectron = H["nelec"]; mol.spin = H["ms2"]
    mol.incore_anyway = True; mol.nao_nr = lambda *a, **k: norb
    mol.energy_nuc = lambda *a: H["ecore"]
    mf = scf.RHF(mol) if H["ms2"] == 0 else scf.UHF(mol)
    mf.get_hcore = lambda *a: H["h1"]
    mf.get_ovlp = lambda *a: np.eye(norb)
    mf._eri = ao2mo.restore(8, H["eri"], norb)
    mf.max_cycle = 200; mf.conv_tol = 1e-8
    # deterministic guess: occupy orbitals with the lowest one-electron energies
    order = np.argsort(np.diag(H["h1"]))
    na, nb = (H["nelec"] + H["ms2"]) // 2, (H["nelec"] - H["ms2"]) // 2
    da = np.zeros((norb, norb)); db = np.zeros((norb, norb))
    da[order[:na], order[:na]] = 1; db[order[:nb], order[:nb]] = 1
    dm0 = da + db if H["ms2"] == 0 else np.array([da, db])
    e_det = mf.energy_tot(dm0)
    known = R.info.get("known_hamiltonian") or {}
    gs = known.get("reference_guess_occupation")
    if gs and len(gs) == norb:     # known file: start from the published guess to reproduce the published energy
        occ = np.array([{'2': [1, 1], 'a': [1, 0], 'b': [0, 1], '0': [0, 0]}[x] for x in gs], float)
        da = np.diag(occ[:, 0]); db = np.diag(occ[:, 1])
        dm0 = da + db if H["ms2"] == 0 else np.array([da, db])
        mf.conv_tol = 1e-12
    try:
        nmf = mf.newton(); e = nmf.kernel(dm0); conv = bool(nmf.converged)
    except Exception as ex:
        e, conv = None, False
    R.info.update(guess_determinant_energy=float(e_det), mean_field_energy=None if e is None else float(e),
                  mean_field_type="RHF" if H["ms2"] == 0 else "UHF", mean_field_seconds=round(time.time() - t0, 1), mean_field_converged=conv)
    msg = (f"Energy of a simple starting determinant: {e_det:.6f} Eh. "
           + (f"Lowest mean-field ({R.info['mean_field_type']}) energy found: {e:.6f} Eh." if e is not None else "Mean-field did not run."))
    status = "PASS" if conv else "WARN"
    if not conv: msg += " (did not fully converge; this is only a sanity check, not a rejection)"
    if known and e is not None and "reference_mean_field_energy" in known:
        ref = known["reference_mean_field_energy"]; tol = known.get("reference_tolerance", 0.1)
        ok = abs(e - ref) <= tol
        msg += f" Published reference {ref:.6f} Eh; difference {e - ref:+.4f} Eh ({'within' if ok else 'outside'} {tol} Eh)."
        status = "PASS" if ok else "WARN"
    R.add("sanity_energy", status, msg)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("file"); ap.add_argument("--json"); ap.add_argument("--no-mf", action="store_true")
    a = ap.parse_args()
    R, _ = check(a.file, run_mf=not a.no_mf)
    out = {"verdict": R.verdict, "info": R.info, "checks": R.checks}
    print(f"{R.verdict}: {a.file}")
    for c in R.checks: print(f"  [{c['status']}] {c['check']}: {c['detail']}")
    if a.json: json.dump(out, open(a.json, "w"), indent=2, default=str)
    sys.exit(1 if R.verdict == "REJECTED" else 0)

if __name__ == "__main__":
    main()
