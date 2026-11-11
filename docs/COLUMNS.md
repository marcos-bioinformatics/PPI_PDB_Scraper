# Output columns

Data dictionary of the final tables. Terms (pair, homomeric/heteromeric, chain
cluster, pair cluster, family, class, scheme, tier) are those of the glossary at
the top of [METHODS.md](METHODS.md). Paths are relative to `data_dir` (`data/`).

| Table | Rows | Columns |
|---|---|---|
| `processed/test.csv` | test positives (classes C1, C2, C3) | 120 |
| `processed/trainval_splits.csv` | every positive not in the test set, with its split in each scheme | 121 |
| `processed/validation.csv` | the rows of `trainval_splits.csv` marked `val` in `splits.validation_split` (`split_c40`) | 123 |
| `processed/hard_negatives.csv` | hard negatives | 59 |
| `reports/prodigy_cryst.csv` | optional (`snakemake prodigy_cryst`): PRODIGY-cryst call for every positive | 10 |

- Every positive table has the same 117 columns, one row per pair (one PISA
  interface between two protein chains). `test.csv` adds the class columns
  (`test_bucket`, `test_subbucket`, `test_unseen_chains`); `trainval_splits.csv`
  adds the split and weight columns (`split_c60`, `split_c40`,
  `cluster_size_train_c40`, `family_weight_c40`); `validation.csv` has those and
  `test_bucket`, `test_subbucket` (both `VAL`).
- `molA_` and `molB_` are the two chains of the pair: in the positives, the
  first and second molecule of the PISA interface (PISA's order, not sorted); in
  the hard negatives, the two chains in alphabetical order of chain ID. Such
  columns are described once for both chains.
- Names with a level follow the configuration: `cXX` are the levels of
  `clustering.thresholds` (95, 90, 60, 40, 30), `split_cXX` the levels of
  `splits.schemes` (60, 40), and the weight columns the scheme of
  `splits.validation_split` (40).
- *Step* is the script that creates the column (`scripts/NN_*`); *opt.* is
  `scripts/annotate_prodigy_cryst.py`. *Empty* means an empty CSV field (NA when
  read). Booleans are written `TRUE`/`FALSE` in the positives (R) and
  `True`/`False` in the hard negatives (pandas).

## Identifiers

| Column | Step | Description |
|---|---|---|
| `row_id` | 10 | Unique key of a positive: `pdb_id\|assembly_id\|interface_id\|pair_index_in_interface\|molA_chain_id\|molB_chain_id` (e.g. `6a34\|1\|1\|1\|A\|B`). |
| `pdb_id` | 02 | PDB ID, lower case. |
| `assembly_id` | 02 | Biological assembly number (with `assemblies.selection: first`, the lowest current assembly ID of the entry). |
| `interface_id` | 02 | PISA interface ID within the assembly. |
| `pair_index_in_interface` | 02 | Always 1 (one pair per PISA interface); part of `row_id`. |
| `molA_chain_id`, `molB_chain_id` | 02 | Author chain ID in the assembly file (`auth_asym_id`); symmetry copies carry a suffix, e.g. `A-2`. |
| `chain_instance_id_A`, `chain_instance_id_B` | 10 | Chain instance: `pdb_id\|assembly_id\|chain_id`. |
| `chain_pair_id` | 10 | The two chain instance IDs in sorted order, joined with `\|\|` (order-invariant). |
| `pair_seq_exact` | 10 | The two chain sequences (`molA_seq`, `molB_seq`) in sorted order, joined with `\|\|`; equal values mean identical sequence pairs (the key of the exact-pair cap of step 10). |
| `source_file` | 02 | Path of the PISA-lite JSON file the row was parsed from (`raw/pisa/<xy>/<id>_assembly<n>_interfaces.json.gz`). |

## Test classes and splits

| Column | Step | Description |
|---|---|---|
| `test_bucket` | 11, 12 | Class of a test pair: `C3` (both proteins unseen), `C2` (one unseen), `C1` (both seen, pair unseen), with "seen" at `test_set.level` (40 %). `VAL` in `validation.csv`. |
| `test_subbucket` | 11, 12 | Class and kind (from `is_homo`): `C3_homo`, `C3_hetero`, `C2_hetero`, `C1_homo`, `C1_hetero`. `VAL` in `validation.csv`. |
| `test_unseen_chains` | 11 | Chains whose chain cluster at `test_set.level` is held out of training: `AB` (C3), `A` or `B` (C2); empty for C1. |
| `split_c60`, `split_c40` | 12 | Split of the pair in each scheme (pair clusters at 60 % and 40 %): `train` or `val`. Independent splits. Always `val` for `split_c40` in `validation.csv`. |
| `cluster_size_train_c40` | 12 | Rows of the pair's family (`pair_cluster_c40`) that are `train` in `split_c40`; 1 if the family has no training rows (every `val` row). |
| `family_weight_c40` | 12 | `1 / cluster_size_train_c40 ^ splits.weight_exponent` (0.5 by default: a family of n rows weighs √n in total); 1 for `val` rows. |

## PISA interface

Copied from the interface block of the PISA-lite file (step 02); one row per
protein–protein interface with exactly two protein molecules.

| Column | Step | Description |
|---|---|---|
| `interface_area` | 02 | Solvent-accessible area buried per side (Å², half of the total buried area); 800 < value ≤ 5,400 after step 07. |
| `p_value` | 02 | PISA ΔiG P-value (0–1): probability of a more hydrophobic surface patch of the same size by chance; ≤ 0.4 after step 07. |
| `solvation_energy` | 02 | PISA ΔiG, solvation free energy gain of interface formation (kcal/mol; negative = hydrophobic interface). |
| `stabilization_energy` | 02 | PISA interface energy that adds hydrogen-bond and salt-bridge terms to ΔiG (kcal/mol); in this dataset ≈ `solvation_energy` − 0.44 × H-bonds − 0.37 × salt bridges (fit, not documented in the file). |
| `molA_interface_residues`, `molB_interface_residues` | 02 | Residues of the chain with buried surface area > 0 in this interface (count of positive values in PISA's per-residue `buried_surface_areas`). |
| `number_interface_residues_raw` | 02 | PISA-lite's `number_interface_residues`, as is. Not an interface count: it equals the length of the second molecule's (molB) residue list in the PISA file, i.e. about the modelled residues of molB (= `molB_resolved_count` in 86 % of rows). Do not use as interface size. |
| `number_hydrogen_bonds` | 02 | Interchain hydrogen bonds (PISA). |
| `number_salt_bridges` | 02 | Interchain salt bridges (PISA). |
| `number_disulfide_bonds` | 02 | Interchain disulfide bonds (PISA). |
| `number_covalent_bonds` | 02 | Interchain covalent bonds (PISA). |
| `number_other_bonds` | 02 | PISA "other bonds": interchain atom pairs listed in the file's `other_bonds` block (2.0–4.0 Å in a sample), i.e. van der Waals contacts. |
| `protein_molecule_count` | 02 | Protein molecules in the interface (always 2). |
| `other_molecule_count` | 02 | Non-protein molecules in the interface (always 0). |
| `molA_molecule_id`, `molB_molecule_id` | 02 | PISA molecule ID within the interface (1 and 2). |
| `molA_molecule_class`, `molB_molecule_class` | 02 | PISA molecule class (always `Protein`). |

## PISA assembly

Copied from the assembly block of the same PISA-lite file (PISA's assembly
analysis of the whole biological assembly, including ligands). All but
`assembly_interface_count` are empty when PISA-lite gives no assembly values (82
of 40,616 positives in the current run).

| Column | Step | Description |
|---|---|---|
| `assembly_dissociation_energy` | 02 | ΔG_diss, free energy of dissociation of the assembly (kcal/mol; > 0 = stable in PISA's model). For one-interface assemblies ≈ −`stabilization_energy` − `assembly_entropy`. |
| `assembly_solvation_energy_gain` | 02 | ΔG_int, solvation free energy gain of assembly formation (kcal/mol). |
| `assembly_entropy` | 02 | TΔS_diss, rigid-body entropy change on dissociation (kcal/mol). |
| `assembly_dissociation_area` | 02 | Area of the interfaces broken on dissociation (Å²). |
| `assembly_accessible_surface_area` | 02 | Solvent-accessible surface area of the assembly (Å²). |
| `assembly_buried_surface_area` | 02 | Total buried surface area of the assembly, both sides of every interface (Å²). |
| `assembly_formula` | 02 | PISA formula: upper case = macromolecules, lower case = ligands, `(n)` = copies, e.g. `A(2)a(4)`. |
| `assembly_composition` | 02 | Chains (author IDs) and ligands (`[code](copies)`) of the assembly, e.g. `AB[HEZ](4)`. |
| `assembly_interface_count` | 02 | Number of interfaces PISA lists for the assembly (any molecule type). |

## Chains (per chain, `molA_` / `molB_`)

From the assembly mmCIF file (first model; step 06) and the RCSB entity
annotations (step 04). "Resolved" = a SEQRES position (`_pdbx_poly_seq_scheme`)
with at least one non-hydrogen atom of occupancy > 0.

| Column | Step | Description |
|---|---|---|
| `molA_entity_id_mmcif`, `molB_entity_id_mmcif` | 06 | Entity ID of the chain in the assembly file (`molA_entity_id`, `molB_entity_id` in the hard negatives). |
| `molA_seq`, `molB_seq` | 06 | Canonical one-letter sequence of the entity (`_entity_poly.pdbx_seq_one_letter_code_can`), upper case, terminal `X` (caps such as ACE/NH2) removed. |
| `molA_seq_len`, `molB_seq_len` | 06 | Length of `seq` (residues); 40–1,000 after step 07. |
| `molA_seqres_count_chain`, `molB_seqres_count_chain` | 06 | SEQRES positions of the chain (`_pdbx_poly_seq_scheme`); can exceed `seq_len` by the terminal caps or `X` removed from `seq`. Positives only. |
| `molA_resolved_count`, `molB_resolved_count` | 06 | Resolved SEQRES positions. Positives only. |
| `molA_missing_residues`, `molB_missing_residues` | 06 | `seqres_count_chain` − `resolved_count`. Positives only. |
| `molA_coverage`, `molB_coverage` | 06 | `resolved_count / seqres_count_chain` (fraction 0–1); > 0.86 after step 07 with `coverage_metric: coverage`. |
| `molA_internal_coverage`, `molB_internal_coverage` | 06 | Resolved residues / SEQRES positions from the first to the last resolved one (fraction 0–1): disordered tags and tails do not count, internal gaps do. |
| `molA_missing_nterm`, `molB_missing_nterm` | 06 | Unresolved SEQRES positions before the first resolved residue. |
| `molA_missing_cterm`, `molB_missing_cterm` | 06 | Unresolved SEQRES positions after the last resolved residue. |
| `molA_uniprot_ids`, `molB_uniprot_ids` | 04 | UniProt accessions of the entity (RCSB, SIFTS mapping), `;`-separated; empty if none. |
| `molA_membrane_sources`, `molB_membrane_sources` | 04 | Membrane-protein annotation sources of the entity among PDBTM, OPM, MemProtMD, mpstruc (`;`-separated); empty if none. Pairs with a source in `scope.exclude_membrane` (PDBTM) are excluded, so only the others appear. |
| `molA_is_antibody`, `molB_is_antibody` | 04 | The entity has a SAbDab antibody annotation in RCSB. Flagged only, unless `scope.exclude_antibodies`. |

## Entry metadata

RCSB Data API, one value per entry (step 04, joined in step 06).

| Column | Step | Description |
|---|---|---|
| `exptl_method` | 04 | Experimental method of the entry: `quality_filter.exptl_method` if the entry has it, else its first method. Always `X-RAY DIFFRACTION` after step 07. |
| `exptl_methods` | 04 | Every method of the entry, `;`-separated (e.g. `X-RAY DIFFRACTION;NEUTRON DIFFRACTION`). |
| `resolution` | 04 | High-resolution limit of the refinement of `exptl_method` (`_refine.ls_d_res_high`, Å); empty if RCSB has none (not filtered unless `quality_filter.max_resolution` is set). |
| `r_free` | 04 | R-free of that refinement (`_refine.ls_R_factor_R_free`, fraction 0–1); < 0.28 after step 07 (entries without R-free are excluded). |
| `title` | 04 | Entry title (`_struct.title`). |
| `keywords` | 04 | `_struct_keywords.pdbx_keywords`; can be empty. |
| `deposit_date` | 04 | Deposition date (YYYY-MM-DD). |
| `initial_release_date` | 04 | First release date (YYYY-MM-DD); compared with `scope.max_release_date`. |
| `revision_date` | 04 | Date of the latest revision (YYYY-MM-DD). |
| `revision` | 04 | Current revision, `major.minor`. Read it as text with care: the tables were read as numbers on the way, so `1.0` appears as `1` (positives) and `1.10` as `1.1`. |

## Homomeric / heteromeric

| Column | Step | Description |
|---|---|---|
| `is_homo_entity` | 06 | Both chains are the same entity (`molA_entity_id_mmcif` = `molB_entity_id_mmcif`). |
| `seq_identity` | 06 | 1 for the same entity; otherwise identical residues of a global alignment with free end gaps (BLOSUM62, gap open −10, extend −0.5) divided by the length of the shorter sequence (fraction 0–1), so a truncated construct of the same protein scores 1. |
| `is_homo` | 06 | Homomeric pair: `is_homo_entity` or `seq_identity` ≥ `enrichment.homo_similarity` (0.80). |

## Quality tier

| Column | Step | Description |
|---|---|---|
| `has_entry_meta` | 06 | RCSB entry metadata found. Always TRUE (other rows are dropped in step 06). |
| `has_chainA_asm`, `has_chainB_asm` | 06 | The chain and its sequence were found in the assembly file. Always TRUE (other rows are dropped in step 06). |
| `tier` | 07 | `gold` if `p_value` < `quality_filter.gold_max_p_value` (0.2), else `silver` (0.2 ≤ P ≤ 0.4). |

## Structure checks

Step 08, on the assembly file: first model, polymer residues only, heavy atoms
(not H or D) of occupancy > 0, one conformer per residue (the alternate location
with the highest mean occupancy). The 20 standard amino acids are always
allowed; `special_residues.ignore` (ACE, NH2) are neither counted nor forbidden;
any other residue not in `special_residues.allowed` (MSE, SEC) excludes the
chain.

| Column | Step | Description |
|---|---|---|
| `molA_has_special`, `molB_has_special` | 08 | The chain has at least one residue of `special_residues.allowed`. |
| `molA_number_special`, `molB_number_special` | 08 | Number of modelled residues of `special_residues.allowed` in the chain. |
| `molA_special_residues`, `molB_special_residues` | 08 | Allowed special residues and counts, `NAME:count` joined with `;` (e.g. `MSE:4`); empty if none. |
| `molA_special_status`, `molB_special_status` | 08 | `ok`, or `forbidden_or_error` (forbidden residue or structure error; such rows are excluded). Always `ok`. |
| `molA_special_error`, `molB_special_error` | 08 | Forbidden residues (`NAME:count;...`) or the structure error (`chain_not_found`, `parse_error:<exception>`). Always empty. |
| `has_special`, `number_special`, `special_residues` | 08 | The same three values for both chains together (counts summed). |
| `n_atom_contacts` | 08 | Interchain heavy-atom pairs within `contacts.cutoff` (≤ 5 Å). |
| `n_residue_contacts` | 08 | Residue pairs (one per chain) with at least one heavy-atom pair within `contacts.cutoff`; ≥ `contacts.min_residue_contacts` (1) in the positives. |
| `min_interchain_distance` | 08 | Shortest heavy-atom distance between the two chains (Å, 3 decimals). |
| `molA_interface_gaps`, `molB_interface_gaps` | 08 | Unresolved SEQRES positions of the chain within ±`contacts.interface_gap_window` (2) positions of an interface residue (a residue with a heavy atom within `contacts.cutoff` of the other chain, not the PISA interface). Reported only, unless `contacts.max_interface_gaps` is set. |
| `n_clash_atoms` | 08 | Heavy atoms of either chain closer than `contacts.clash_distance` (< 1.7 Å) to the other chain (sum over both chains). |
| `max_clash_fraction` | 08 | Larger of the two chains' fractions of heavy atoms in clashes (0–1, 4 decimals); pairs above `contacts.max_clash_fraction` (0.30) are excluded. |

## Clusters

Step 09: MMseqs2 all-vs-all search of the unique chain sequences; two
sequences are linked at XX % when their identity is ≥ XX % with
`clustering.coverage` (0.9) of either sequence aligned and E-value ≤ 1e-3;
clusters are the connected components (`clustering.method:
connected_components`).

| Column | Step | Description |
|---|---|---|
| `molA_cluster_c95`, `_c90`, `_c60`, `_c40`, `_c30`; `molB_cluster_c95` … `_c30` | 09 | Chain cluster at that identity: ID of its representative sequence (first 16 hex characters of the SHA-1 of the sequence; for connected components, the smallest ID in the cluster). |
| `pair_cluster_c95`, `_c90`, `_c60`, `_c40`, `_c30` | 09 | Pair cluster: the two chain cluster IDs in sorted order joined with `\|\|` (order-invariant). `pair_cluster_c40` is the family. |

## Hard negatives only

`hard_negatives.csv` (step 13): two protein chains of the same assembly with no
PISA interface between them. Columns shared with the positives, with the same
meaning: `pdb_id`, `assembly_id`, `molA_chain_id`, `molB_chain_id`, the chain
columns except `seqres_count_chain`, `resolved_count` and `missing_residues`
(the entity ID is `molA_entity_id`, `molB_entity_id`), the entry metadata,
`n_atom_contacts`, `n_residue_contacts`, `min_interchain_distance` and the
cluster columns (taken from the clustering of the positives in step 09; every
negative chain has a sequence clustered there). There are no PISA, homo/hetero,
tier, special-residue, clash or gap columns and no `row_id`.

| Column | Step | Description |
|---|---|---|
| `label` | 13 | Always 0 (negative). |
| `n_protein_chains_assembly` | 13 | Protein (polypeptide) chains in the assembly file, symmetry copies included. |
| `n_atom_contacts`, `n_residue_contacts` | 13 | As in step 08; always 0, since `hard_negatives.min_distance` (6 Å) exceeds `contacts.cutoff` (5 Å). |
| `min_interchain_distance` | 13 | As in step 08 (Å); ≥ `hard_negatives.min_distance` (6 Å) and ≤ `max_distance` if set. |
| `intact_evidences` | 13 | IntAct evidences (MITAB 2.7 rows, PSICQUIC) summed over every pair of UniProt accessions of the two chains (different accessions only). Empty if there is no such pair (a chain without UniProt accession, or both chains with only the same one) or a query failed. Absent if `hard_negatives.intact.enabled` is false. |
| `intact_direct_evidences` | 13 | The subset of those evidences that are binary (no spoke or matrix expansion) and of type direct interaction (MI:0407 or a PSI-MI descendant). Empty in the same cases. With `intact.drop_direct`, negatives with a value > 0 are dropped, so it is 0 or empty. |
| `split_c60`, `split_c40` | 13 | Split in each scheme of step 12: `test` if a chain occurs only in test positives or belongs to a chain cluster held out by step 11; else `val` if a chain occurs only in validation positives of the scheme; else `train`. All negatives of a family (`hard_negatives.family_level`, c40) take the most restrictive split of the group (test > val > train). |
| `split` | 13 | Copy of the scheme of `splits.validation_split` (`split_c40`). |

## PRODIGY-cryst report

`reports/prodigy_cryst.csv` (optional): one row per positive of `test.csv` and
`trainval_splits.csv` (so validation rows too), sorted by `row_id`. An
annotation; no step filters on it. Structures are prepared as PRODIGY-cryst
does (first model, one conformer per atom, HETATM residues such as MSE and
water removed), keeping only the two chains of the pair.

| Column | Step | Description |
|---|---|---|
| `row_id`, `pdb_id`, `assembly_id`, `molA_chain_id`, `molB_chain_id` | opt. | Copied from the positive; join on `row_id`. |
| `prodigy_class` | opt. | PRODIGY-cryst call: `BIO` (biological) or `XTAL` (crystal contact); empty if not scored. |
| `prodigy_p_bio` | opt. | Classifier probability of `BIO` (0–1, 4 decimals); empty if not scored. |
| `prodigy_ics` | opt. | Interchain residue contacts: residue pairs with any atoms within 5.0 Å (PRODIGY-cryst definition); empty if not scored. |
| `prodigy_link_density` | opt. | `prodigy_ics` / (residues of A in contact × residues of B in contact); empty if not scored. |
| `prodigy_error` | opt. | Why a pair was not scored: `non-standard residue <NAME>` (in ATOM records), `parse <exception>: …`, or the PRODIGY-cryst exception (truncated to 120 characters); empty if scored. |

## Notes

**Meanings not verified in the code.** These columns are copied from PISA-lite
without a definition in the file; their descriptions follow PISA's documentation
and checks on the data, not this pipeline's code:
`stabilization_energy` (relation to H-bonds and salt bridges fitted on the
data), `number_other_bonds` (inferred from the bond list in the file),
`number_covalent_bonds`, and the `assembly_*` values (`dissociation_energy`,
`solvation_energy_gain`, `entropy`, `dissociation_area`,
`accessible_surface_area`, `buried_surface_area`, `formula`, `composition`).

**Constant columns** in the final tables (kept for compatibility):

- Forced by the code or the filters, in every positive table:
  `pair_index_in_interface` (1), `protein_molecule_count` (2),
  `other_molecule_count` (0), `molA_molecule_class`/`molB_molecule_class`
  (`Protein`), `has_entry_meta`, `has_chainA_asm`, `has_chainB_asm` (TRUE),
  `molA_special_status`/`molB_special_status` (`ok`),
  `molA_special_error`/`molB_special_error` (empty), `exptl_method`
  (`X-RAY DIFFRACTION`, also in the hard negatives).
- Constant in the current run but not forced: `molA_molecule_id` (1),
  `molB_molecule_id` (2) (PISA's numbering), `number_disulfide_bonds` (0),
  `number_covalent_bonds` (0), `assembly_id` (1 in every table).
- `validation.csv`: `test_bucket`, `test_subbucket` (`VAL`), `split_c40`
  (`val`), `cluster_size_train_c40`, `family_weight_c40` (1).
- `hard_negatives.csv`: `label` (0), `n_atom_contacts`, `n_residue_contacts` (0),
  `intact_direct_evidences` (0 or empty).
