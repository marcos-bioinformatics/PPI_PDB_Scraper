# From PPI_PDB_Scrapper Seq8 to v1

This repository is the cleaned-up, reproducible version of the working pipeline
used to build the PPI_PDB_Scrapper Seq8 dataset (developed August 2025 – March
2026). The scientific design and the order of the steps are unchanged; what
changed is how the steps are wired together, plus a few bug fixes.

## Script mapping

| v1 | Seq8 |
|---|---|
| `01_download_pisa.py` | new (the PISA JSON came from a manual bulk download of the PDBe FTP) |
| `02_parse_pisa.py` | `2_Json_scrapping.py` |
| `03_select_assemblies.R` | `master_pisa_out/master_pisa.R` |
| `04_download_structures.py` | `download_missing_assemblies.sh` (assembly 1 had no script); `1_ASU_CIF_download.py` replaced by the RCSB Data API |
| `05_audit_downloads.py` | `0_audit_db_folders.py` |
| `06_enrich_mmcif.py` | `3_CIF_scrapping_ASU_Assembly.py` |
| `07_quality_filter.R` | `mmcif_enrichment_out/masterr_csv_164k.R` |
| `08_check_structures.py` | `4_filter_special_residues_SR.py` + contact check (KD-tree adapted from `NeighborSearch.py` / `Test_distance.py` of the ML code) |
| `09_cluster_chains.py` | `5_clustering.py` |
| `10_prune_redundancy.R` | `mmseqs_chain_clustering_out/2_pruning_98.R` |
| `11_build_test_set.R` | `mmseqs_chain_clustering_out/3_test_set_datatable.R` |
| `12_trainval_splits.R` | `mmseqs_chain_clustering_out/4_trainval_splits.R` (+ the validation export done by hand) |
| `13_hard_negatives.py` | new (Seq5 had in-assembly negatives with label conflicts; see below) |
| `14_leakage_audit.py` | new |
| `analysis/clustering_stats.R` | `mmseqs_chain_clustering_out/1_clustering_stats.R` |
| — | `5_last_filtering.R` (the length cut is now a parameter of step 07) |

## Workflow changes

- **One config file.** All paths and thresholds live in `config/config.yaml`;
  both the Python and the R scripts read it.
- **Snakemake** encodes the order of the steps; no more manual renaming of
  intermediate files between scripts.
- **The ID list is applied first.** PISA data is downloaded only for the listed
  entries, instead of parsing the whole PDB and filtering afterwards.
- **Structures are downloaded on demand.** Only the assembly files of the pairs
  that survive the PISA parsing are fetched (in Seq8: 190,214 assemblies
  downloaded, 93,022 needed). The asymmetric-unit files were only used for five
  entry fields (method, resolution, R-free, title, keywords); these now come
  from the RCSB Data API, which saves ~70,000 downloads (~20 GB).
- **Compressed storage.** mmCIF files are kept as `.cif.gz` and PISA JSON as
  `.json.gz` (about 10x smaller).

## New in v1

- **Contact check (step 08).** Every PISA pair must have heavy-atom contacts in
  the assembly structure; this guards against chain-mapping errors.
- **Stale PISA assemblies (step 01).** The PDBe dump keeps files of assemblies
  removed in later entry revisions (e.g. 5u1j); only current assemblies are used.
- **Hard negatives (step 13).** Seq5 labelled every non-interacting chain pair
  as negative, which produced pairs that were positive and negative at the same
  time (same sequences). v1 keeps hetero pairs only, removes cluster-level
  conflicts, requires a minimum distance and assigns each negative to the split
  of its proteins.

## Defaults that differ from Seq8

To reproduce Seq8 exactly, set these values in `config/config.yaml`:

| Parameter | v1 default | Seq8 |
|---|---|---|
| `quality_filter.max_chain_length` | 600 | 900 (the saved script said 800) |
| `quality_filter.min_chain_length` | 39 (at least 40 residues) | 12 (at least 13) |
| `quality_filter.max_solvation_density` | null (off) | -0.0068 |
| `enrichment.homo_similarity` | on an alignment | positional identity, no alignment |
| `clustering.method` | `connected_components` | `set_cover` |
| `test_set.c2_pick` | `random` | `worst` (lowest-quality pair of each cluster) |
| `hard_negatives.conflict_level` | `c40` | — (new in v1) |

`quality_filter.coverage_metric: internal_coverage` and `quality_filter.max_resolution`
are new options, off by default.

## Homology leakage

The Seq8 splits were built on MMseqs2 set-cover clusters. Auditing them
(step 14) shows leaks at the design threshold of every set: 115 of 3,664
validation pairs at 90 %, 8 and 6 of 100 C2 pairs at 70 and 80 %, and 4 of 144
and 4 of 147 C3 chains at 40 and 50 %. v1 clusters by connected components,
which removes them by construction.

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
3. **ASU download.** RCSB now serves files with gzip transfer encoding and the
   compressed `Content-Length`; the old downloader compared it with the
   decompressed size, deleted every file and, on a rerun, deleted the files
   already on disk. Downloads now request the raw bytes and are checked by gzip
   integrity (and ASU files are no longer needed).
4. **ID splitting (step 03).** `strsplit(x, "[,\\s]+")` splits on the letter
   "s" (inside brackets `\s` is not whitespace for R's default regex engine).
   It was harmless with the upper-case RCSB exports used for Seq8, but dropped
   every ID containing an "s" (~9 %) once IDs are lower-case; now
   `[,[:space:]]+`, and step 05 reports any entry lost in step 03.
5. **PISA pre-filter (step 02).** Interfaces with a missing area or p-value
   passed the filter; they are now dropped.
6. **Missing chains (step 06).** A missing sequence (pandas NaN) was counted as
   present in `has_chain*_asm`.
7. **Test set (step 11).** `setorder()` referred to `.easy_score` instead of
   `easy_score`, so the saved script stopped with an error.
8. Row order is now deterministic in every Python step.

## Deep audit (1 November 2026)

Four independent reviews of v1 (acquisition, structure and science,
clustering and orchestration, R scripts) were followed by these changes; the
pipeline was then rerun from a clean clone and checked with
`tests/check_outputs.py`.

- **Selenocysteine / pyrrolysine (step 06).** `U` and `O` are not in BLOSUM62
  and crashed the alignment of any hetero pair containing them (e.g. 1cc1).
- **Special residues (steps 08, 13; inherited from Seq8).** Residues were
  classified by name with Biopython's `is_aa`, which mixes polymer and ligands:
  free ligands such as SAH, SEP or ORN excluded chains, while GFP chromophores,
  UNK and D-amino acids in the polymer were not detected and selenocysteine was
  never counted. Only polymer residues (`label_seq_id`) are now used, for the
  classification and for the contacts; ACE/NH2 caps are ignored.
- **C2 chain cap (step 11; inherited from Seq8).** `cc <- tmp_cc` inside
  `pick_group()` updated a local copy, so the one-use-per-protein cap was not
  enforced between homo and hetero picks: 25 proteins appeared twice in the Seq8
  C2 buckets and 21 of 100 C2-c70 pairs had a protein that no longer appeared in
  training. Fixed with `<<-`; on the Seq8 data both counts drop to 0.
- **Short chains.** Peptides cannot reach the MMseqs2 E-value, so near-identical
  short chains escaped the homology clustering; the default minimum length is
  now 40 residues.
- **Downloads.** Transfers cut mid-way are retried, a corrupt gzip is
  downloaded again, GraphQL queries are retried, steps 01/04 stop when too many
  files fail, PISA files updated by PDBe are refreshed, the redundant
  `{id}_assembly{n}.json` is no longer downloaded (half of the PISA requests),
  and ID lists with BOM, UTF-16, quotes, semicolons or entity IDs (`4HHB_1`) are
  read correctly. Obsolete entries are reported as such.
- **Hard negatives (step 13).** Missing PISA files stop the step instead of
  hiding contacts; structure checks run before the per-entry de-duplication;
  negatives are capped per pair cluster across entries; splits are assigned per
  pair-cluster family at 40 % (`hard_negatives.family_level`, the conflict level).
- **Leakage audit (step 14).** Searches in both directions, compares the test
  set with training and validation, audits every split scheme and the hard
  negatives.
- **Snakemake.** Scripts are declared as inputs (code changes trigger reruns),
  step 01 reruns when its config sections change, step 13 declares its PISA and
  structure inputs, and unsupported `--config`/`--configfile` options stop the
  run instead of being ignored. YAML values such as `1e-3` are read as numbers in
  both Python and R.
- **Step 10 performance.** The exact-pair cap and the cluster pruning were
  rewritten with one global ranking instead of one `group_modify` call per
  group: 39 s instead of 9-17 min on the Seq8 data, with an output identical to
  the Seq8 `pruned98` table (same rows, order and values).
- **Other R fixes.** `null` pruning parameters fall back to the function
  defaults, `test_set.seed: null` and an empty training pool stop with an
  error, under-filled test buckets print a warning, helper columns are no longer
  written, empty strings are read as NA in step 12, and `exptl_method: null`
  disables the method filter.

## R scripts

The R scripts keep their original logic. The only edits are mechanical:
reading the config (input/output paths and the thresholds that were literals),
writing outputs with `row.names = FALSE`, sending the histograms of step 07 to
`reports/07_quality_filter_plots.pdf`, the `easy_score` fix above, and the
validation export at the end of step 12, dropping the helper column `.area`
from the output of step 10 and making its comparison run in `pair_cluster_95`
optional (`pruning.compare_pair_cluster_95`), which halves the run time.

## Regression against the Seq8 files

The R steps were run on the original Seq8 intermediate tables:

| Step | Result |
|---|---|
| 07 quality filter | identical (130,053 rows) with `max_chain_length: 900`, the value actually used for Seq8 (the saved script said 800) |
| 10 pruning | identical (38,135 rows, same `row_id`s and columns) |
| 11 test set | same bucket sizes; C2 99/100 and 96/100 identical, C1/C3 differ: the saved script is not the exact version that produced the Seq8 test set |
| 12 splits | identical: split_c90/c80/c60, family weights and validation set |
