# Notes: known issues, Snakemake tips and troubleshooting

## Known issues

- The homo filter of step 13 only removes identical sequences, not the
  `is_homo` definition of step 06 (the conflict filter usually catches the rest).
- The hard negatives come from few large complexes (35 % from proteasome and
  24 % from ribosome entries in the reference run; the 289 test negatives come
  from 77 families at 40 %, 54 % of them from these two kinds of complex) and
  none is homomeric: evaluate them per family and against heteromeric positives
  only. The IntAct check finds no contradiction rather than confirming them:
  almost all IntAct records for these pairs are expanded from n-ary complexes.
- C1 has 96 homomeric and 54 heteromeric pairs in the reference run, not 75 and
  75: only 54 heteromeric pair clusters leave both proteins with other partners
  in the train/validation pool, and homomeric pairs fill the class.
- The splits are by pair cluster, so a protein can appear in training and in
  validation with different partners. A *seen* protein (test classes) is one
  with a row in the train/validation pool, so in a given split it may be only
  in validation. Validation is mostly not "known proteins with new partners":
  in `split_c40`, 77 % of the validation pairs have neither protein (chain
  cluster at 40 %) in training.
- In `split_c40`, one influenza hemagglutinin family of 852 rows makes up two
  thirds of the heteromeric validation rows: a cluster that does not fit the
  10 % target is still taken if that leaves the total closer to it.
- `is_homo` means the same entity or an aligned identity of at least 0.80, so
  it includes pairs of a protein and its truncated form (e.g. the p66–p51
  heterodimer of HIV-1 reverse transcriptase), and heteromeric pairs can share
  a UniProt accession (e.g. hemagglutinin HA1–HA2, cleaved from one precursor).
- The cap on identical sequence pairs (step 10) does not look at the
  interface: of the 1,839 sequence pairs with interfaces of different area in
  one entry (95 % homomeric), 1,279 keep rows of only one of them.
- The optional PRODIGY-cryst annotation does not score pairs with a
  selenocysteine (16 in the reference run), a position modelled both as
  selenomethionine and methionine (1) or a duplicated atom in the file (1).
- Steps 01–04 query live services: a rerun of those steps can download updated
  PISA files and metadata while keeping older mmCIF files. To freeze a dataset,
  keep `data/raw/` and `data/interim/` and do not rerun steps 01–04.

## Working with Snakemake

Snakemake only reruns what is missing or out of date. Each rule receives the
config sections it uses as `params` and declares its script (and the shared
`common.py` / `common.R`) as inputs, so changing a section of the config or
editing a script repeats that step and everything downstream of it. Editing
`scripts/common.py` therefore reruns the whole pipeline, downloads included;
after a change that does not affect the early steps, `snakemake --touch --until
<rule>` marks them as up to date. A config can inherit from another one with
`base: config/config.yaml`; `--config` only accepts `cfg=<file>`, and a target
goes before it. The path of the config file is part of every rule's command,
so running with another config file on a `data_dir` that already holds outputs
makes Snakemake rerun every step from 01, downloads included ("Code has changed
since last execution"). Either edit `config/config.yaml` in place, give the new
config its own `data_dir`, or add `--rerun-triggers mtime params input
software-env`, with which only the steps whose config sections changed rerun. Every script can also be run by hand from the repository root,
e.g. `python scripts/06_enrich_mmcif.py --config config/config.yaml` or
`Rscript scripts/07_quality_filter.R config/config.yaml`.

```bash
snakemake --cores 8 --forcerun download_structures    # rerun a step (e.g. retry failed downloads)
snakemake --cores 8 --touch --until check_structures  # mark steps 01-08 as up to date without running them
snakemake --dag | dot -Tsvg > dag.svg                 # draw the workflow
```

## Troubleshooting

- **"Output … has older modification time than input … clock skew" (WSL2).**
  The WSL2 clock can jump back after the computer sleeps. Run the same
  `snakemake` command again; it resumes from the step that stopped.
- **Failed downloads.** Each file is retried `download_retries` times. Steps 01
  and 04 stop if more than `max_failed_fraction` of the files still fail (run
  the same command again to retry them; files already on disk are kept) and
  report the remaining failures in their manifests.
- **Stopping a run.** On SIGTERM, Snakemake lets the running job finish (which
  can take hours during a download) and exits without losing anything; Ctrl-C
  stops at once and removes the incomplete outputs of that step.
