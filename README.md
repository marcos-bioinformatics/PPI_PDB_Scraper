# BeastAffinity

A reproducible pipeline that turns a list of PDB entries into a curated,
non-redundant dataset of **protein–protein interfaces** with PISA energetics,
structure-quality metrics and homology-aware train / validation / test splits.

Starting from any set of IDs exported from the RCSB Advanced Search, it:

1. downloads the PISA interface analysis of every entry from PDBe,
2. keeps protein–protein interfaces as candidate chain pairs,
3. downloads only the biological-assembly mmCIF files those pairs need, plus their
   entry metadata (method, resolution, R-free) from the RCSB Data API,
4. adds sequence, coverage and experimental-quality metadata to every pair,
5. filters by quality, removes chains with non-standard residues and checks
   that every PISA interface has real atomic contacts,
6. clusters all chains with MMseqs2 at eight identity thresholds,
7. prunes redundancy and builds a test set (C1/C2/C3) and cluster-aware train/val splits,
8. optionally adds **hard negatives**: chains of the same complex that are not in contact.

## Pipeline

Python steps do the heavy extraction (downloads, JSON and mmCIF parsing,
clustering); R steps take the data-analysis decisions (filters, pruning, splits).

| Step | Script | What it does | Main output |
|---|---|---|---|
| 01 | `01_download_pisa.py` | Reads the ID lists; downloads PISA JSON from PDBe | `raw/pisa/` |
| 02 | `02_parse_pisa.py` | One row per protein pair in each PISA interface | `02_master_pisa.csv` |
| 03 | `03_select_assemblies.R` | Restricts to the ID list; lists the assemblies needed | `03_needed_pairs.csv` |
| 04 | `04_download_structures.py` | Downloads the needed assembly mmCIF files; entry metadata from the RCSB Data API | `raw/assemblies/`, `04_entry_metadata.csv` |
| 05 | `05_audit_downloads.py` | Checks that every needed file is on disk | `05_audit_downloads.json` |
| 06 | `06_enrich_mmcif.py` | Sequences, coverage, homo/hetero, method, resolution, R-free | `06_enriched.csv` |
| 07 | `07_quality_filter.R` | Exploratory analysis + quality filter + gold/silver tier | `07_quality_filtered.csv` |
| 08 | `08_check_structures.py` | Drops chains with non-standard amino acids and pairs without atomic contacts | `08_checked.csv` |
| 09 | `09_cluster_chains.py` | MMseqs2 chain clustering at 98–40 % identity | `09_clustered.csv` |
| 10 | `10_prune_redundancy.R` | Exact-duplicate soft cap + pruning in `pair_cluster_98` | `10_pruned.csv` |
| 11 | `11_build_test_set.R` | Test set: C1 (train-like), C2 (new pairs), C3 (unseen proteins) | `processed/test.csv` |
| 12 | `12_trainval_splits.R` | Cluster-aware train/val splits and family weights | `processed/trainval_splits.csv`, `processed/validation.csv` |
| 13 | `13_hard_negatives.py` | Non-contacting chain pairs of the same assembly, assigned to splits | `processed/hard_negatives.csv` |

`analysis/clustering_stats.R` is an optional exploratory report on the clusters
(`snakemake --cores 1 clustering_stats`).

```mermaid
flowchart LR
  ids[RCSB ID lists] --> s01[01 PISA download] --> s02[02 parse PISA] --> s03[03 select assemblies]
  s03 --> s04[04 mmCIF download] --> s05[05 audit] --> s06[06 enrich]
  s06 --> s07[07 quality filter] --> s08[08 structure checks] --> s09[09 MMseqs2]
  s09 --> s10[10 prune] --> s11[11 test set] --> s12[12 train/val splits] --> s13[13 hard negatives]
```

## Installation

```bash
conda env create -f environment.yml
conda activate beastaffinity
```

The environment provides Python 3.11, R 4.x (dplyr, tidyr, data.table, yaml),
MMseqs2 and Snakemake.

## Usage

1. Run a query in the [RCSB Advanced Search](https://www.rcsb.org/search/advanced),
   click **Download IDs** and put the `.txt` file(s) in `pdb_ids/`.
2. Review the parameters in `config/config.yaml` (every threshold lives there).
3. Run the pipeline from the repository root:

```bash
snakemake -n              # dry run: list the steps that would run
snakemake --cores 8       # run everything
```

Useful Snakemake commands:

```bash
snakemake --cores 8 --until quality_filter           # stop after step 07
snakemake --cores 8 --forcerun download_structures   # rerun a step (e.g. retry failed downloads)
snakemake --cores 8 --config cfg=config/test.yaml    # use another config file
snakemake --dag | dot -Tsvg > dag.svg                # draw the workflow
```

Snakemake only reruns what is missing or out of date: if you change a section
of the config (for instance `quality_filter`), the next run repeats that step
and everything downstream of it. Every script can also be run by hand from the
repository root, e.g. `python scripts/06_enrich_mmcif.py --config config/config.yaml`
or `Rscript scripts/07_quality_filter.R config/config.yaml`.

### Quick test

`config/test.yaml` runs the whole pipeline on the 1,500 IDs in `tests/pdb_ids/`
(about 0.5 GB of downloads, ~6 minutes) into `data_test/`:

```bash
snakemake --cores 8 --config cfg=config/test.yaml
```

## Choosing the parameters

Every threshold is in `config/config.yaml`. The defaults reproduce the
BeastAffinitySeq8 dataset, but most of them are a matter of what the dataset
is for; the notes below explain what each group means.

**PISA pre-filter (`pisa`, step 02).** PISA reports `interface_area` as the
solvent-accessible area buried *per side* (half of the total buried surface),
so 800 Å² means about 1,600 Å² buried in total, above the typical range of
crystal-packing contacts. `p_value` is PISA's ΔiG P-value: the probability of
finding a more hydrophobic patch of the same size by chance; values below 0.5
suggest an interaction-specific interface. PISA-lite does not provide the
complexation significance score (CSS) of classic PISA, so specificity relies
on area, P-value and solvation energy. These two values are only a pre-filter
(they decide which structures are downloaded); the final thresholds are set in
`quality_filter` and can only be stricter. Only the assemblies an entry has
today are used (`only_current_assemblies`): the PDBe dump can keep assemblies
removed in later revisions.

**Quality filter (`quality_filter`, step 07).** Experimental method, R-free,
minimum coverage of both chains (resolved / SEQRES residues), chain length
range, interface area range, P-value and solvation density
(`solvation_energy / interface_area`, nearly independent of interface size).
`max_chain_length` is often limited by the model you plan to train.

**Structure checks (`special_residues`, `contacts`, step 08).** Chains with
non-standard amino acids other than the allowed ones are dropped. Every pair
must also have at least `min_residue_contacts` residue pairs with heavy atoms
within `cutoff` Å; this catches PISA interfaces that cannot be mapped onto the
structure.

**Redundancy (`pruning`, step 10).** Many entries contain several copies of
the same complex (non-crystallographic symmetry, several assemblies), and close
homologues appear across entries. Identical sequence pairs are soft-capped to
`exact_K_homo` / `exact_K_heter` rows, and pair clusters (`pair_key`, 98 %
identity by default) with more than `activate_n` rows are pruned to
`K_redund_*` rows when their interface areas are almost identical, or to
`K_var_*` rows chosen across interface-area bins otherwise. Keeping several
copies (2 homo / 4 hetero by default) preserves close variants such as point
mutants; set the caps to 1 for a strictly non-redundant dataset.

**Test set and splits (`test_set`, `splits`, steps 11–12).** C1 rows are
train-like, C2 rows are new pairings of known proteins (pair clusters at 70 and
80 %), C3 rows contain proteins unseen in training (chain clusters at 40 and
50 %). Train/validation splits are made by pair cluster, stratified by
homo/hetero, with a family weight of 1/√(cluster size).

**Hard negatives (`hard_negatives`, step 13).** Two protein chains of the
same biological assembly with no PISA interface between them. Homo pairs are
excluded by default (a copy of a protein not touching another copy usually
contradicts a homo interface elsewhere in the complex). A negative is dropped
if its pair of chain clusters (`conflict_level`) is in contact in any
downloaded assembly, if it fails the quality thresholds of step 07, if a chain
has forbidden residues or if the chains are closer than `min_distance` Å. One
negative is kept per entry and sequence pair, and at most `max_per_entry` per
entry, so that ribosomes and capsids do not dominate. Each negative inherits
the split of its proteins (test if any chain is in the test set or in a
held-out C3 cluster).

## Outputs

All files are written under `data/` (configurable with `data_dir`):

```
data/
├── raw/         PISA JSON and assembly mmCIF files (gzipped), with download manifests
├── interim/     tables of steps 01–11
├── processed/   test.csv, trainval_splits.csv, validation.csv, hard_negatives.csv
├── reports/     audits, parse errors, quality-filter plots
└── logs/        one log per step
```

Disk space: about 25 KB of PISA JSON and 0.27 MB of mmCIF per assembly, i.e.
~0.5 GB for the 1,500-ID test and roughly 30 GB for a list of ~160,000 IDs.

## Data sources

- **RCSB PDB** — assembly mmCIF files (`files.rcsb.org`), the Data API
  (`data.rcsb.org`, current assemblies and entry metadata) and the Advanced Search.
  Burley et al., *Nucleic Acids Res.* (2023).
- **PDBe PISA** — interface and assembly analysis
  (`ftp.ebi.ac.uk/pub/databases/msd/pdb-assemblies-analysis`).
  Krissinel & Henrick, *J. Mol. Biol.* 372, 774–797 (2007).
- **MMseqs2** — Steinegger & Söding, *Nat. Biotechnol.* 35, 1026–1028 (2017).
- **Biopython** — Cock et al., *Bioinformatics* 25, 1422–1423 (2009).

## Roadmap

On 21 July 2027 the wwPDB switches to extended PDB IDs (`pdb_` + 8 characters)
and a reorganised archive. Version 2 will accept both ID formats and use the new
archive layout. See [docs/MIGRATION.md](docs/MIGRATION.md) for the history of the pipeline.
