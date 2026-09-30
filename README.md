# PubChemQC FCIDUMP build

Scripts that turn PubChemQC B3LYP/6-31G* records (Nakata and colleagues, CC-BY-4.0, from MolSSI's public
mirror on Hugging Face, `molssiai-hub/pubchemqc-b3lyp`) into CAS(8e,8o) FCIDUMP files with PySCF, for the
BECOME energy catalog at https://become.belovedecosystem.com/energy/.

- `scripts/pcq_prep.py` streams one source file and keeps neutral singlets with at most 150 orbitals.
- `scripts/pcq_convert.py` rebuilds RHF in 6-31G*, checks it against PubChemQC, and writes the FCIDUMP.
- `scripts/fcidump_check.py` independently checks each FCIDUMP.
- `scripts/pcq_run.py` runs a batch in parallel and writes `catalog.jsonl`.

The `Build PubChemQC FCIDUMPs` workflow runs this across many GitHub runners at once.
