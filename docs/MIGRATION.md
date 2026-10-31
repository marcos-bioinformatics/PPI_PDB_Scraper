# From BeastAffinitySeq8 to v1

This repository is the cleaned-up, reproducible version of the working pipeline
used to build the BeastAffinitySeq8 dataset (developed August 2025 – March
2026). The scientific design and the order of the steps are unchanged; what
changed is how the steps are wired together, plus a few bug fixes.

## Script mapping

| v1 | Seq8 |
|---|---|
| `01_download_pisa.py` | new (the PISA JSON came from a manual bulk download of the PDBe FTP) |
| `02_parse_pisa.py` | `2_Json_scrapping.py` |
| `03_select_assemblies.R` | `master_pisa_out/master_pisa.R` |
| `04_download_structures.py` | `1_ASU_CIF_download.py` + `download_missing_assemblies.sh` (assembly 1 had no script) |
| `05_audit_downloads.py` | `0_audit_db_folders.py` |
| `06_enrich_mmcif.py` | `3_CIF_scrapping_ASU_Assembly.py` |
| `07_quality_filter.R` | `mmcif_enrichment_out/masterr_csv_164k.R` |
| `08_filter_special_residues.py` | `4_filter_special_residues_SR.py` |
| `09_cluster_chains.py` | `5_clustering.py` |
| `10_prune_redundancy.R` | `mmseqs_chain_clustering_out/2_pruning_98.R` |
| `11_build_test_set.R` | `mmseqs_chain_clustering_out/3_test_set_datatable.R` |
| `12_trainval_splits.R` | `mmseqs_chain_clustering_out/4_trainval_splits.R` (+ the validation export done by hand) |
| `analysis/clustering_stats.R` | `mmseqs_chain_clustering_out/1_clustering_stats.R` |
| — | `5_last_filtering.R` (the length cut is now a parameter of step 07) |

## Workflow changes

- **One config file.** All paths and thresholds live in `config/config.yaml`;
  both the Python and the R scripts read it.
- **Snakemake** encodes the order of the steps; no more manual renaming of
  intermediate files between scripts.
- **The ID list is applied first.** PISA data is downloaded only for the listed
  entries, instead of parsing the whole PDB and filtering afterwards.
- **Structures are downloaded on demand.** Only the ASU and assembly files of
  the pairs that survive the PISA parsing are fetched (in Seq8: 164,000 ASU +
  190,214 assemblies downloaded, 69,577 entries / 93,022 assemblies needed).
- **Compressed storage.** mmCIF files are kept as `.cif.gz` and PISA JSON as
  `.json.gz` (about 10x smaller).

## Bug fixes

1. **Chain matching (step 06).** Chains were indexed by `label_asym_id` but PISA
   reports `auth_asym_id`, so whenever both differ (antibodies, entries with
   ligands listed before proteins, ...) a pair could get the sequence and
   coverage of another chain, or none. Chains are now matched by
   `auth_asym_id` of the assembly file (including symmetry copies such as
   `A-2`). On a 1,500-entry sample, 5.7 % of the pairs had a wrong sequence;
   all 160 chains checked against the structure now match their sequence.
2. **Coverage of selenomethionine (step 06).** Only `ATOM` records were counted
   as resolved, so MSE (written as `HETATM`) lowered the coverage of SeMet
   structures. All polymer residues are now counted.
3. **ASU download (step 04).** RCSB now serves files with gzip transfer
   encoding and the compressed `Content-Length`; the old downloader compared it
   with the decompressed size, deleted every file and, on a rerun, deleted the
   files already on disk. Downloads now request the raw bytes and are checked
   by gzip integrity.
4. **PISA pre-filter (step 02).** Interfaces with a missing area or p-value
   passed the filter; they are now dropped.
5. **Missing chains (step 06).** A missing sequence (pandas NaN) was counted as
   present in `has_chain*_asm`.
6. **Test set (step 11).** `setorder()` referred to `.easy_score` instead of
   `easy_score`, so the saved script stopped with an error.
7. Row order is now deterministic in every Python step.

## R scripts

The R scripts keep their original logic. The only edits are mechanical:
reading the config (input/output paths and the thresholds that were literals),
writing outputs with `row.names = FALSE`, sending the histograms of step 07 to
`reports/07_quality_filter_plots.pdf`, the `easy_score` fix above, and the
validation export at the end of step 12.

## Regression against the Seq8 files

The R steps were run on the original Seq8 intermediate tables:

| Step | Result |
|---|---|
| 07 quality filter | identical (130,053 rows) with `max_chain_length: 900`, the value actually used for Seq8 (the saved script said 800) |
| 10 pruning | identical (38,135 rows, same `row_id`s and columns) |
| 11 test set | same bucket sizes; C2 99/100 and 96/100 identical, C1/C3 differ: the saved script is not the exact version that produced the Seq8 test set |
| 12 splits | identical: split_c90/c80/c60, family weights and validation set |
