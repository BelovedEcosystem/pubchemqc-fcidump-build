# PubChemQC FCIDUMP build

Scripts that turn PubChemQC B3LYP/6-31G* records (Nakata and colleagues, CC BY 4.0, from MolSSI's public
mirror on Hugging Face, `molssiai-hub/pubchemqc-b3lyp`) into CAS(8e,8o) FCIDUMP files with PySCF, for the
BECOME energy catalog at https://become.belovedecosystem.com/energy/.

- `scripts/pcq_prep.py` streams one source file of a subset, keeps neutral singlets with at most `maxorb` orbitals
  (basis-function count minus heavy-atom count), applies the optional selection settings below, drops CIDs in an
  optional skip list, and with `shards=auto` sizes the shard count to the candidate count so dense files finish
  inside the per-job time budget.
- `scripts/pcq_prep_group.py` and `scripts/pcq_pool.py` are the pooled mode: many source files are scanned in grouped
  prep jobs, their candidates are pooled (repeated CIDs dropped) and re-sharded.
- `scripts/pcq_convert.py` rebuilds RHF in 6-31G*, checks it against PubChemQC, and writes the FCIDUMP.
- `scripts/fcidump_check.py` independently checks each FCIDUMP.
- `scripts/pcq_run.py` runs a batch in parallel and writes `catalog.jsonl`.
- `scripts/pcq_merge.py` merges the `out/` directories of several shards into one tree (first shard wins on
  a duplicate CID; duplicates are listed in `merge_duplicates.jsonl`).
- `scripts/pcq_finalize.py` re-checks every file and writes the build library: `library/catalog.jsonl` and, for
  each molecule, `library/published/<CID>/` with `FCIDUMP_pcq<CID>_cas8e8o.gz` (plus a full-space file for molecules
  with at most 60 orbitals), `pcq<CID>.xyz`, `entry.json`, `SHA256SUMS` for the uncompressed text and
  `SHA256SUMS-gz` for the `.gz` archives.
- `scripts/pcq_package.py` packages one or more build libraries into the deposited layout: the active-space files
  uncompressed as `become-fcidump-pcq/fcidump/pcq<CID>_6-31g_cas8e8o.fcidump`, the geometries as
  `become-fcidump-pcq/xyz/pcq<CID>.xyz`, the index `catalog_index.csv` (one row per entry), `SHA256SUMS.txt` and,
  with `--tar`, a reproducible `become-fcidump-pcq.tar.gz`.

The `Build PubChemQC FCIDUMPs` workflow runs prep and build across many GitHub runners at once, then a
`finalize` job merges every shard and uploads the build library as the `library` artifact. The deposited catalog is
made from the four sources' libraries with

    python scripts/pcq_package.py <out> <chon300nosalt library> <chnopsfcl300nosalt library> \
        <chnopsfclnakmgca500 library> <b3lyp_pm6 library> --tar [--entry-numbers <CSV of cid,i>]

Column `i` of `catalog_index.csv` is the catalog entry number; without `--entry-numbers` the entries are numbered
in packaging order.

## Selection settings for the published catalog

The published catalog (125,999 molecules) comes from four PubChemQC sources. Each was built with the workflow
inputs below; inputs not listed keep their defaults (`limit=0`, `pool=0`, `minorb=0`, `no_multifrag=0`, empty
`exclude_z`, `allow_z` and `skip`). `shards`, `per_shard`, the pooled shard count and `budget_min` only spread the
work over runners and do not change which molecules are selected.

| Source | `subset` | `file_from`-`file_to` | `maxorb` | Other selection inputs |
|---|---|---|---|---|
| `b3lyp_pm6_chon300nosalt` | `chon300nosalt` | 0-86 | 150 | none |
| `b3lyp_pm6_chnopsfcl300nosalt` | `chnopsfcl300nosalt` | 0-148 | 150 | `skip=data/skip_chnopsfcl300nosalt.txt.gz` |
| `b3lyp_pm6_chnopsfclnakmgca500` | `chnopsfclnakmgca500` | 0-346 | 150 | pooled (`pool=8`), `exclude_z=11,12,19,20` (Na, Mg, K, Ca), `no_multifrag=1`, `skip=data/skip_chnopsfclnakmgca500.txt.gz` |
| `b3lyp_pm6`, source file 1 (CIDs 253,697-507,790) | `b3lyp_pm6` | 1-1 | 200 | `minorb=150`, `allow_z=1,5,6,7,8,9,14,15,16,17,35`, `no_multifrag=1`, `skip=data/skip_b3lyp_pm6.txt.gz`, pooled (`pool=16`) |

Each skip list holds the CIDs that earlier builds had already processed, so a molecule is built once.
`pcq_convert.py` then applies the same gates to every source (neutral closed-shell singlet, elements H, B, C, N, O,
F, Si, P, S, Cl, Br, no formal charges, PubChem identity kept after PM6, single fragment, B3LYP energy reproduced,
stable RHF).

## Method notes

The functional is named explicitly (`b3lypg`), degenerate orbitals are fixed by a fixed generic operator, a phase
rule is applied to every orbital, gzip is written with mtime 0, every file is re-read and its integrals checked (no
correlated energy is computed or stored), and molecules too small for an 8e8o active space are rejected. Repeat runs
agree within about 1e-9 in every integral but are not byte-identical, so compare files by tolerance, not by hash.

## License

The code in this repository is licensed under the Apache License 2.0 (see `LICENSE`). Copyright 2026 Beloved
Ecosystem, Inc. The PubChemQC input data are not part of this repository; they are licensed separately under CC BY 4.0 by
their authors.
