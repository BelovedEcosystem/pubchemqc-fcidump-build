# PubChemQC FCIDUMP build

Scripts that turn PubChemQC B3LYP/6-31G* records (Nakata and colleagues, CC-BY-4.0, from MolSSI's public
mirror on Hugging Face, `molssiai-hub/pubchemqc-b3lyp`) into CAS(8e,8o) FCIDUMP files with PySCF, for the
BECOME energy catalog at https://become.belovedecosystem.com/energy/.

- `scripts/pcq_prep.py` streams one source file and keeps neutral singlets with at most 150 orbitals.
- `scripts/pcq_convert.py` rebuilds RHF in 6-31G*, checks it against PubChemQC, and writes the FCIDUMP.
- `scripts/fcidump_check.py` independently checks each FCIDUMP.
- `scripts/pcq_run.py` runs a batch in parallel and writes `catalog.jsonl`.

The `Build PubChemQC FCIDUMPs` workflow runs this across many GitHub runners at once.
- `scripts/pcq_finalize.py` re-checks every file and writes the published library (`SHA256SUMS` for the
  uncompressed text, `SHA256SUMS-gz` for the `.gz` archives).

Recipe version 2 adds: the functional named explicitly (`b3lypg`), degenerate orbitals fixed by a fixed
generic operator, a phase rule on every orbital, gzip written with mtime 0, integral re-read checks, a direct
CASCI cross-check, and a gate rejecting molecules too small for an 8e8o active space. Repeat runs agree within
about 1e-9 in every integral but are not byte-identical, so compare files by tolerance, not by hash.
