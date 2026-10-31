# BeastAffinity

A reproducible pipeline that turns a list of PDB entries into a curated,
non-redundant dataset of **protein–protein interfaces** with PISA energetics,
structure-quality metrics and homology-aware train / validation / test splits.

Starting from any set of IDs exported from the RCSB Advanced Search, it:

1. downloads the PISA interface analysis of every entry from PDBe,
2. keeps protein–protein interfaces as candidate chain pairs,
3. downloads only the mmCIF files those pairs need (asymmetric unit + biological assemblies),
4. adds sequence, coverage and experimental-quality metadata to every pair,
5. filters by quality and removes chains with non-standard residues,
6. clusters all chains with MMseqs2 at eight identity thresholds,
7. prunes redundancy and builds a test set (C1/C2/C3) and cluster-aware train/val splits.

## Pipeline

Python steps do the heavy extraction (downloads, JSON and mmCIF parsing,
clustering); R steps take the data-analysis decisions (filters, pruning, splits).

| Step | Script | What it does | Main output |
|---|---|---|---|
| 01 | `01_download_pisa.py` | Reads the ID lists; downloads PISA JSON from PDBe | `raw/pisa/` |
| 02 | `02_parse_pisa.py` | One row per protein pair in each PISA interface | `02_master_pisa.csv` |
| 03 | `03_select_assemblies.R` | Restricts to the ID list; lists the assemblies needed | `03_needed_pairs.csv` |
| 04 | `04_download_structures.py` | Downloads the needed ASU and assembly mmCIF files | `raw/asu/`, `raw/assemblies/` |
| 05 | `05_audit_downloads.py` | Checks that every needed file is on disk | `05_audit_downloads.json` |
| 06 | `06_enrich_mmcif.py` | Sequences, coverage, homo/hetero, method, resolution, R-free | `06_enriched.csv` |
| 07 | `07_quality_filter.R` | Exploratory analysis + quality filter + gold/silver tier | `07_quality_filtered.csv` |
| 08 | `08_filter_special_residues.py` | Drops chains with non-standard amino acids | `08_special_residues.csv` |
| 09 | `09_cluster_chains.py` | MMseqs2 chain clustering at 98–40 % identity | `09_clustered.csv` |
| 10 | `10_prune_redundancy.R` | Exact-duplicate soft cap + pruning in `pair_cluster_98` | `10_pruned.csv` |
| 11 | `11_build_test_set.R` | Test set: C1 (train-like), C2 (new pairs), C3 (unseen proteins) | `processed/test.csv` |
| 12 | `12_trainval_splits.R` | Cluster-aware train/val splits and family weights | `processed/trainval_splits.csv`, `processed/validation.csv` |

`analysis/clustering_stats.R` is an optional exploratory report on the clusters
(`snakemake --cores 1 clustering_stats`).

```mermaid
flowchart LR
  ids[RCSB ID lists] --> s01[01 PISA download] --> s02[02 parse PISA] --> s03[03 select assemblies]
  s03 --> s04[04 mmCIF download] --> s05[05 audit] --> s06[06 enrich]
  s06 --> s07[07 quality filter] --> s08[08 special residues] --> s09[09 MMseqs2]
  s09 --> s10[10 prune] --> s11[11 test set] --> s12[12 train/val splits]
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
(about 0.7 GB of downloads, ~6 minutes) into `data_test/`:

```bash
snakemake --cores 8 --config cfg=config/test.yaml
```

## Outputs

All files are written under `data/` (configurable with `data_dir`):

```
data/
├── raw/         PISA JSON (gzipped) and mmCIF files, with download manifests
├── interim/     tables of steps 01–11
├── processed/   test.csv, trainval_splits.csv, validation.csv
├── reports/     audits, parse errors, quality-filter plots
└── logs/        one log per step
```

Disk space: about 25 KB of PISA JSON per assembly and 0.3 MB per mmCIF file,
i.e. ~0.7 GB for the 1,500-ID test and roughly 50 GB for a list of ~160,000 IDs.

## Data sources

- **RCSB PDB** — mmCIF files (`files.rcsb.org`) and the Advanced Search.
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
