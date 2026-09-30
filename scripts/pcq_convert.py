"""Strict FCIDUMP builder for one PubChemQC B3LYP/6-31G*//PM6 record.
usage: pcq_convert.py <record.json> <outdir>
We take only the structure from PubChemQC and redo everything ourselves in PySCF:
  gates      neutral closed-shell, light elements, no formal charges, PubChem identity kept after PM6
  structure  our B3LYP/6-31G*(6D) energy must reproduce PubChemQC's to 1e-4 Eh
  RHF        6-31G* (spherical), tightly converged, internally and externally stable
  orbitals   canonical RHF orbitals with a fixed phase: in each orbital the largest-magnitude AO
             coefficient (first one, within 1e-6) is made positive
  files      8e8o active-space FCIDUMP (+ full FCIDUMP when <= FULL_MAX_NAO orbitals);
             each file is re-read and must reproduce its integrals to 1e-10 and the RHF energy to 1e-8 Eh;
             CASCI from the file must match a direct PySCF CASCI to 1e-8 Eh; gzip written with mtime 0
  degenerate orbitals (same occupation, energies within 1e-5 Eh) are fixed by diagonalizing a fixed
             generic point-charge Coulomb operator inside each block, before the phase rule
Recipe version 2."""
import sys, os, re, json, gzip, shutil, numpy as np
from pyscf import gto, dft, scf, mcscf, fci, ao2mo, lib
from pyscf.tools import fcidump
EV = 27.211386245988
rec = json.load(open(sys.argv[1])); outdir = sys.argv[2]
FULL_MAX_NAO = int(os.environ.get("FULL_MAX_NAO", "60")); BASIS = os.environ.get("BASIS", "6-31g*")
SIMPLE = {1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 35}
Z = [int(z) for z in rec["atomic-numbers"]]; xyz = np.array(rec["coordinates"]).reshape(-1, 3)
def gate(ok, msg):
    if not ok: raise SystemExit(f"REJECT {msg}")
gate(rec["multiplicity"] == 1 and rec.get("pubchem-multiplicity", 1) == 1, "not a singlet")
gate(abs(rec["charge"]) < 1e-9 and rec.get("pubchem-charge", 0) == 0, "charged")
gate(set(Z) <= SIMPLE, f"element outside simple set {sorted(set(Z) - SIMPLE)}")
smi = rec.get("pubchem-isomeric-smiles", "")
gate(not re.search(r"\[[^\]]*[+-][^\]]*\]", smi), "formal charges (zwitterion)")
strip = lambda s: re.split(r"/[tbms]", s)[0]
gate(strip(rec["pubchem-inchi"]) == strip(rec["obabel-inchi"]), "PM6 structure no longer matches the PubChem compound")
gate("." not in smi, "multiple fragments")
atoms = [(z, tuple(c)) for z, c in zip(Z, xyz)]
# structure check: reproduce PubChemQC's own B3LYP/6-31G*(Cartesian d) energy
mol6d = gto.M(atom=atoms, basis="6-31g*", cart=True, unit="Angstrom", verbose=0)
gate(mol6d.nao == rec["basis-count"], f"basis size {mol6d.nao} != {rec['basis-count']}")
# xc named explicitly: "b3lypg" is B3LYP with VWN-RPA correlation (what PySCF >= 2.3 calls "b3lyp"; "b3lyp5" differs by
# ~0.07 Eh for a small organic molecule). Chosen because it reproduces the PubChemQC energies within the 1e-4 Eh gate.
ks = dft.RKS(mol6d, xc="b3lypg"); ks.grids.level = 3; ks.conv_tol = 1e-9
e_b3 = ks.kernel(); e_ref = rec["total-energy"] / EV
gate(ks.converged and abs(e_b3 - e_ref) < 1e-4, f"structure check failed (B3LYP diff {e_b3 - e_ref:.2e} Eh)")
# RHF in the FCIDUMP basis
mol = gto.M(atom=atoms, basis=BASIS, cart=False, unit="Angstrom", verbose=0)
nelecas, ncas = 8, 8
gate(mol.nelectron >= nelecas and mol.nao >= (mol.nelectron - nelecas) // 2 + ncas, f"too small for {nelecas}e{ncas}o active space")
mf = scf.RHF(mol); mf.conv_tol = 1e-10; mf.conv_tol_grad = 1e-6; mf.max_cycle = 200
mf.kernel()
if not mf.converged:
    mf = scf.newton(mf); mf.kernel()
gate(mf.converged, "RHF did not converge")
_, _, st_i, st_e = mf.stability(internal=True, external=True, return_status=True)
gate(st_i, "RHF internally unstable")
if os.environ.get("STRICT_EXTERNAL", "1") == "1":
    gate(st_e, "RHF unstable toward UHF (open-shell character)")
DEG_TOL = 1e-5
def fix_degenerate(C, eps, occ):
    """Degenerate orbitals (same occupation, energies within DEG_TOL) are only defined up to a rotation.
    Fix it: inside each block, diagonalize a fixed generic one-electron operator (Coulomb potential of two
    point charges at fixed offsets from the atom centroid) and order by its eigenvalue."""
    # 1/r potentials of point charges at two fixed off-centre points: unlike dipole/quadrupole terms they carry
    # every multipole order, so they also split pairs like delta orbitals of linear molecules.
    cen = mol.atom_coords().mean(axis=0); W = 0
    for w, off in ((1.0, (1.37, 0.71, 0.43)), (0.61, (-0.52, 1.19, -0.88))):
        with mol.with_rinv_origin(cen + np.array(off)): W = W + w * mol.intor("int1e_rinv")
    C = C.copy(); blocks = 0; worst = np.inf; k = 0; n = C.shape[1]
    while k < n:
        j = k + 1
        while j < n and occ[j] == occ[k] and eps[j] - eps[j - 1] < DEG_TOL: j += 1
        if j - k > 1:
            Cb = C[:, k:j]; lam, U = np.linalg.eigh(Cb.T @ W @ Cb); C[:, k:j] = Cb @ U
            blocks += 1; worst = min(worst, float(np.diff(lam).min()))
        k = j
    return C, blocks, (None if blocks == 0 else worst)
def canon(C):
    C = C.copy()
    for k in range(C.shape[1]):
        a = np.abs(C[:, k]); i = int(np.flatnonzero(a >= a.max() - 1e-6)[0])
        C[:, k] *= (np.sign(C[i, k]) or 1.0)
    return C
Cfix, deg_blocks, deg_split = fix_degenerate(mf.mo_coeff, mf.mo_energy, mf.mo_occ)
gate(deg_blocks == 0 or deg_split > 1e-6, f"degenerate orbitals could not be fixed (split {deg_split})")
mf.mo_coeff = canon(Cfix)  # CASCI below reads orbitals from mf
C, occ, eps = mf.mo_coeff, mf.mo_occ, mf.mo_energy
nao = C.shape[1]; nocc = mol.nelectron // 2
orth = float(np.abs(C.T @ mol.intor("int1e_ovlp") @ C - np.eye(nao)).max())
gate(orth < 1e-8, f"orbitals not orthonormal ({orth:.1e})")
gap = float(eps[nocc] - eps[nocc - 1])
gate(gap > 0, "no HOMO-LUMO gap")
os.makedirs(outdir, exist_ok=True); os.chdir(outdir)
def e_rhf_from(h1, eri, ecore, no):
    return float(ecore + 2 * np.trace(h1[:no, :no]) + 2 * np.einsum("iijj", eri[:no, :no, :no, :no]) - np.einsum("ijji", eri[:no, :no, :no, :no]))
# active space 8e8o around the Fermi level
cas = mcscf.CASCI(mf, ncas, nelecas)
h1, ecore = cas.get_h1eff(C)
Ca = C[:, cas.ncore:cas.ncore + ncas]
h2 = ao2mo.restore(1, ao2mo.full(mol, Ca), ncas)
lo, hi = cas.ncore, cas.ncore + ncas  # active-space edges fall inside a degenerate block?
cas_edge_degenerate = bool((lo > 0 and eps[lo] - eps[lo - 1] < DEG_TOL) or (hi < nao and eps[hi] - eps[hi - 1] < DEG_TOL))
name = f"FCIDUMP_pcq{rec['cid']}_cas8e8o"
fcidump.from_integrals(name, h1, h2, ncas, nelecas, nuc=ecore, ms=0, tol=1e-12)
r = fcidump.read(name, verbose=False)
dev_h1_cas = float(np.abs(np.asarray(r["H1"]) - h1).max())
dev_eri_cas = float(np.abs(ao2mo.restore(1, r["H2"], ncas) - h2).max())
dev_ecore_cas = float(abs(r["ECORE"] - ecore))
gate(max(dev_h1_cas, dev_eri_cas, dev_ecore_cas) < 1e-10, f"active-space file does not reproduce its integrals ({max(dev_h1_cas, dev_eri_cas, dev_ecore_cas):.1e})")
e_cas_file = e_rhf_from(r["H1"], ao2mo.restore(1, r["H2"], r["NORB"]), r["ECORE"], nelecas // 2)
gate(abs(e_cas_file - mf.e_tot) < 1e-8, f"active-space file does not reproduce RHF energy ({e_cas_file - mf.e_tot:.1e})")
e_casci = float(fci.direct_spin1.kernel(r["H1"], ao2mo.restore(1, r["H2"], ncas), ncas, (4, 4), tol=1e-12)[0] + r["ECORE"])
cas.fcisolver.conv_tol = 1e-12
e_casci_direct = float(cas.kernel()[0])
gate(abs(e_casci - e_casci_direct) < 1e-8, f"CASCI from file disagrees with direct CASCI ({e_casci - e_casci_direct:.1e})")
gate(e_casci <= mf.e_tot + 1e-9, "CASCI energy above RHF")
full = nao <= FULL_MAX_NAO; e_full_file = None; dev_h1_full = dev_eri_full = dev_ecore_full = None
if full:
    fname = f"FCIDUMP_pcq{rec['cid']}_full"
    h1f = C.T @ mf.get_hcore() @ C; eri8 = ao2mo.full(mol, C)
    fcidump.from_integrals(fname, h1f, ao2mo.restore(8, eri8, nao), nao, mol.nelectron, nuc=mol.energy_nuc(), ms=0, tol=1e-12)
    r = fcidump.read(fname, verbose=False)
    dev_h1_full = float(np.abs(np.asarray(r["H1"]) - h1f).max())
    dev_eri_full = float(np.abs(ao2mo.restore(1, r["H2"], nao) - ao2mo.restore(1, eri8, nao)).max())
    dev_ecore_full = float(abs(r["ECORE"] - mol.energy_nuc()))
    gate(max(dev_h1_full, dev_eri_full, dev_ecore_full) < 1e-10, f"full file does not reproduce its integrals ({max(dev_h1_full, dev_eri_full, dev_ecore_full):.1e})")
    e_full_file = e_rhf_from(r["H1"], ao2mo.restore(1, r["H2"], nao), r["ECORE"], nocc)
    gate(abs(e_full_file - mf.e_tot) < 1e-8, f"full file does not reproduce RHF energy ({e_full_file - mf.e_tot:.1e})")
meta = dict(cid=rec["cid"], formula=rec["formula"], smiles=smi, inchi=rec["pubchem-inchi"], natoms=mol.natm, nelec=mol.nelectron,
            basis=BASIS, norb=nao, e_b3lyp_ours=float(e_b3), e_b3lyp_pubchemqc=e_ref, b3lyp_diff=float(e_b3 - e_ref),
            e_rhf=float(mf.e_tot), e_rhf_from_cas_file=e_cas_file, e_rhf_from_full_file=e_full_file, e_casci_8e8o=e_casci,
            homo_lumo_gap=gap, orth_err=orth, recipe_version=2, xc="b3lypg", degenerate_blocks=deg_blocks, degenerate_min_split=deg_split, cas_edge_degenerate=cas_edge_degenerate, e_casci_direct=e_casci_direct,
            dev_h1_cas=dev_h1_cas, dev_eri_cas=dev_eri_cas, dev_ecore_cas=dev_ecore_cas, dev_h1_full=dev_h1_full, dev_eri_full=dev_eri_full, dev_ecore_full=dev_ecore_full, stable_internal=bool(st_i), stable_external=bool(st_e), full_written=full,
            geometry=[[int(z)] + [float(v) for v in c] for z, c in zip(Z, xyz)])
for f in list(os.listdir(".")):
    if f.startswith("FCIDUMP_") and not f.endswith(".gz"):
        with open(f, "rb") as a, gzip.GzipFile(f + ".gz", "wb", compresslevel=6, mtime=0) as b: shutil.copyfileobj(a, b)
        os.remove(f)
json.dump(meta, open("meta.json", "w"), indent=1)
print(json.dumps({k: v for k, v in meta.items() if k != "geometry"}))
