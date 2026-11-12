# PPI_PDB_Scraper

> **About this project.** PPI_PDB_Scraper is an independent, personal project.
> I conceived it and developed it by hand over 2025 and early 2026, to have a
> reliable source of protein–protein interface structures for training models
> related to protein function. Once the method was mature, I audited it in
> depth with Claude and rebuilt it as a reproducible Snakemake workflow, so that
> the whole dataset can be regenerated from a list of PDB IDs and a single
> configuration file. I am now sharing it with the community.

A reproducible methodology to scrape **protein–protein interfaces (PPIs) from the
Protein Data Bank** and turn them into a curated, non-redundant dataset for
machine learning. It combines PDBe-PISA interface energetics with structural
quality control, MMseqs2 homology clustering, leakage-free train / validation /
test splits and hard negatives. The workflow is orchestrated with Snakemake
(Python + R) and driven by a single configuration file.

## What it does

From a list of PDB IDs (for example, exported from the RCSB Advanced Search):

1. downloads the **PDBe-PISA** interface analysis of every entry and keeps the
   protein–protein interfaces as candidate chain pairs;
2. downloads the biological-assembly mmCIF files those pairs need, with entry
   metadata, membrane and antibody annotations and UniProt accessions from the
   **RCSB Data API**;
3. filters by experimental and interface quality (method, R-free, coverage,
   chain length, interface area, PISA P-value) and checks every pair in the
   coordinates (non-standard residues, atomic contacts, clashes);
4. removes small interfaces that **PRODIGY-cryst**, a classifier independent of
   PISA, calls crystal contacts;
5. clusters all chains with **MMseqs2** (95–30 % identity, connected components)
   and caps redundant copies;
6. builds a test set with the three classes of Park and Marcotte (C1 new pairs
   of seen proteins, C2 one unseen protein, C3 unseen proteins) and
   cluster-aware train/validation splits with family weights;
7. adds **hard negatives** (chains of one assembly that are not in contact),
   checked against **IntAct**;
8. audits the homology leakage between training and every evaluation set.

Every curation decision, with its reason, precedents, effect and biases, is in
[docs/DECISIONS.md](docs/DECISIONS.md).

## Reference run

164,000 X-ray entries (resolution ≤ 3.0 Å), run on 10–11 November 2026 on a
laptop under WSL2 (20 cores, 14 GB RAM):

| Stage | Result |
|---|---|
| PISA (steps 01–02) | 245,757 assemblies, 564,347 interfaces → 204,075 candidate pairs from 69,580 entries |
| Structures (03–05) | first assembly of 66,807 entries: 155,899 pairs |
| Quality and structure checks (06–08b) | 105,440 pairs, 31,722 unique chain sequences |
| Redundancy cap (10) | 43,710 pairs |
| Final dataset (11–13) | **39,747 train/val pairs**, **750 test pairs** (C1 150, C2 300, C3 300), **492 hard negatives**, from 32,989 entries |
| Leakage audit (14) | 0 leaks by sequence in the 11 audited sets; the 31 output checks pass |
| Time | about 2 h 20 min with the files on disk; downloading from scratch added about 6 h for PISA and 1 h 20 min for the structures (run of 2 November 2026, every assembly) |
| Disk | about 22 GB of raw files, 1.2 GB of intermediate tables, 66 MB of final tables |

## Pipeline

![The steps of PPI_PDB_Scraper (01–14 and 08b): four phases, the scripts they run (Python in blue, R in purple) and which step reads the files of which](docs/pipeline.png)

*Each step is a Snakemake rule with the script it runs. Python steps do the
extraction (downloads, parsing, geometry, clustering); R steps make the
data-analysis decisions (filters, redundancy, splits).*

| Step | Script | What it does | Main output |
|---|---|---|---|
| 01 | `01_download_pisa.py` | Downloads PISA JSON from PDBe for the ID lists | `raw/pisa/` |
| 02 | `02_parse_pisa.py` | One row per protein pair in each PISA interface | `02_master_pisa.csv` |
| 03 | `03_select_assemblies.R` | Keeps the first assembly of each entry (by default) and lists the files needed | `03_needed_pairs.csv` |
| 04 | `04_download_structures.py` | Assembly mmCIF files; entry and entity metadata from the RCSB Data API | `raw/assemblies/`, `04_entry_metadata.csv` |
| 05 | `05_audit_downloads.py` | Checks that every needed file is on disk | `05_audit_downloads.json` |
| 06 | `06_enrich_mmcif.py` | Sequences, coverage, homomeric/heteromeric, method, resolution, R-free | `06_enriched.csv` |
| 07 | `07_quality_filter.R` | Descriptive analysis, scope and quality filter, gold/silver tier | `07_quality_filtered.csv` |
| 08 | `08_check_structures.py` | Non-standard residues, atomic contacts and clashes | `08_checked.csv` |
| 08b | `08b_prodigy_filter.py` | Drops pairs below 1,000 Å² that PRODIGY-cryst calls crystal contacts | `08b_checked.csv` |
| 09 | `09_cluster_chains.py` | MMseqs2 chain clustering at 95–30 % identity | `09_clustered.csv` |
| 10 | `10_prune_redundancy.R` | Cap on identical sequence pairs | `10_pruned.csv` |
| 11 | `11_build_test_set.R` | Test set: C1, C2 and C3 | `processed/test.csv` |
| 12 | `12_trainval_splits.R` | Train/validation splits by pair cluster, family weights | `processed/trainval_splits.csv` |
| 13 | `13_hard_negatives.py` | Non-contacting chain pairs, checked against IntAct | `processed/hard_negatives.csv` |
| 14 | `14_leakage_audit.py` | Homology leakage of every evaluation set | `reports/14_leakage_audit.json` |

Optional rules: `prodigy_cryst` (PRODIGY-cryst call for every final positive),
`leakage_control` (leakage audit with an independent, global-alignment
measure) and `clustering_stats`.

## Installation

Linux x86-64 or WSL2, [conda](https://conda-forge.org/download/) or miniforge,
about 22 GB of disk for the reference dataset (0.5 GB for the quick test) and
network access to PDBe, RCSB and IntAct.

```bash
git clone https://github.com/marcos-bioinformatics/PPI_PDB_Scraper.git
cd PPI_PDB_Scraper
conda env create -f environment.yml            # Python 3.11, R, MMseqs2, Snakemake
conda env create -f envs/prodigy_cryst.yml     # PRODIGY-cryst for step 08b (Python 3.8)
conda run -n prodigy_cryst pip install --no-deps prodigy-cryst==1.0.1
conda activate ppi_pdb_scraper
```

`environment.lock.yml` pins the exact versions of the reference run. Without the
PRODIGY-cryst environment, set `prodigy.filter: false` to skip step 08b. The
optional leakage control has its own environment (`envs/leakage_control.yml`).

## Configuration

The pipeline has many parameters, all set in `config/config.yaml`. Some of the
default values follow precedents from deep-learning models or the literature;
these are collected, with the reason, effect and biases of every curation
decision, in [docs/DECISIONS.md](docs/DECISIONS.md). Even so, each user should
choose these values for their own problem.

[docs/METHODS.md](docs/METHODS.md) explains what each parameter does, and a
configuration can inherit from the default one with `base: config/config.yaml`
and change only what it needs (as `config/test.yaml` does).

## Usage

1. Put the ID list(s) in `pdb_ids/` (the `.txt` files of *Download IDs* in the
   RCSB Advanced Search) and, optionally, the query in `scope.input_query`. The
   repository includes the 17 lists of the reference run.
2. Review the parameters in `config/config.yaml` (see *Configuration*).
3. Run from the repository root:

```bash
snakemake -n                       # dry run: what would run and why
snakemake --cores 8                # the whole pipeline
python tests/check_outputs.py      # 31 consistency and leakage checks
snakemake --cores 8 prodigy_cryst  # optional analyses: prodigy_cryst, leakage_control
```

Snakemake reruns only what is missing or out of date: changing a section of the
config or a script repeats that step and everything after it.
`snakemake --cores 8 --until quality_filter` stops after a step, and
`--config cfg=<file>` uses another configuration. [docs/NOTES.md](docs/NOTES.md) has more tips and the
troubleshooting.

**Quick test.** `config/test.yaml` runs the whole pipeline on 1,500 IDs
(`tests/pdb_ids/`, about 0.5 GB of downloads) into `data_test/`, in about
10 minutes:

```bash
snakemake --cores 8 --config cfg=config/test.yaml
python tests/check_outputs.py --config config/test.yaml
```

Step 11 warns that C1 has fewer pairs than requested: 1,500 IDs are too few to
fill it.

## Outputs

All files are written under `data/` (configurable with `data_dir`):

```
data/
├── raw/         PISA JSON and assembly mmCIF files (gzipped), with download manifests
├── interim/     tables of steps 01–11
├── processed/   test.csv, trainval_splits.csv, validation.csv, hard_negatives.csv
├── reports/     audits, parse errors, quality-filter plots, provenance.json
└── logs/        one log per step
```

Each row of the final tables carries the PISA interface metrics, what was
measured in the structure, the entry metadata, the annotations of each chain,
the chain and pair cluster IDs at every threshold, the test class or the split
columns, the gold/silver `tier` and `family_weight_c40`.
`reports/provenance.json` records the ID lists, the configuration, the code
commit and the software versions of each run.

## Documentation

- [docs/METHODS.md](docs/METHODS.md): what each step does, what its parameters
  mean, and how strict the clustering and leakage criteria are.
- [docs/DECISIONS.md](docs/DECISIONS.md): every curation decision, with its
  reason, precedents, effect and biases.
- [docs/COLUMNS.md](docs/COLUMNS.md): the columns of the final tables.
- [docs/NOTES.md](docs/NOTES.md): known issues, Snakemake tips and troubleshooting.
- [docs/MIGRATION.md](docs/MIGRATION.md): how this version relates to the
  working versions (2025–2026) and the bugs fixed on the way.

## Data sources and terms of use

The datasets the pipeline builds keep the terms of their sources:

- **PDB** (RCSB PDB: assembly mmCIF files, Data API) — [CC0 1.0](https://www.wwpdb.org/about/usage-policies).
  Burley et al., *Nucleic Acids Res.* (2023).
- **PDBe PISA** — PDBe data under [CC0](https://www.ebi.ac.uk/pdbe/about/public-data-access-statement),
  with the [EMBL-EBI terms of use](https://www.ebi.ac.uk/about/terms-of-use).
  Krissinel & Henrick, *J. Mol. Biol.* 372, 774–797 (2007).
- **IntAct** (PSICQUIC, step 13) — CC BY 4.0; cite IntAct when you use the
  `intact_*` columns. del Toro et al., *Nucleic Acids Res.* 50, D648–D653 (2022).
- **SIFTS**, **PDBTM** and **SAbDab** annotations, through RCSB.

Software: MMseqs2 (Steinegger & Söding, 2017), PRODIGY-cryst (Elez et al.,
2018; Jiménez-García et al., 2019), Biopython (Cock et al., 2009), Snakemake
(Mölder et al., 2021); for the optional control, EMBOSS (Rice et al., 2000) and
parasail (Daily, 2016). The test classes follow Park and Marcotte (2012).
[docs/DECISIONS.md](docs/DECISIONS.md) lists every reference.

## Roadmap

On 21 July 2027 the wwPDB switches to extended PDB IDs (`pdb_` + 8 characters)
and a reorganised archive. Version 2 will accept both ID formats and the new
archive layout.

## Citation, license and contact

If you use PPI_PDB_Scraper, please cite it as described in
[CITATION.cff](CITATION.cff) (GitHub shows it under *Cite this repository*).
The code is released under the [MIT License](LICENSE). Questions and bug
reports: [GitHub issues](https://github.com/marcos-bioinformatics/PPI_PDB_Scraper/issues).
