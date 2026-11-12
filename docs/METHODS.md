# Methods and parameters

What each step does, what its parameters mean and why the defaults are what
they are. Most thresholds depend on the intended use of the dataset; all of
them are in `config/config.yaml`. Every curation decision, with its reason,
precedents, effect and biases, is recorded in [DECISIONS.md](DECISIONS.md).

**Terms used in this document and in the output columns.**

- *Pair*: two protein chains of one biological assembly, with a PISA interface
  (a positive) or without one (a hard negative).
- *Homomeric / heteromeric pair* (`is_homo`): the two chains are the same entity
  or align with at least `enrichment.homo_similarity` identity (0.80), or not.
  Homomeric therefore includes a protein paired with a truncated form of itself
  (e.g. the p66–p51 heterodimer of HIV-1 reverse transcriptase), and a
  heteromeric pair can join two chains of one UniProt accession (e.g. the HA1
  and HA2 chains of hemagglutinin, cleaved from one precursor). A homomeric pair
  is often two chains of a larger homo-oligomer, not a homodimer.
- *Chain cluster* (`molA_cluster_cXX`, `molB_cluster_cXX`): the connected
  component of chains linked at XX % identity (step 09).
- *Pair cluster* (`pair_cluster_cXX`): the unordered pair of the chain clusters
  of the two chains. A *family* is a pair cluster at 40 %, the level used for the
  family weights and the hard-negative families.
- *Class* (C1, C2, C3; column `test_bucket`, and `test_subbucket` with the
  homo/hetero kind): the novelty class of a test pair (Park and Marcotte, 2012).
- *Scheme* (`split_c60`, `split_c40`): a train/validation split by pair cluster
  at that identity.
- *Tier* (`tier`): *gold* for a PISA ΔiG P-value below 0.2, *silver* from 0.2 to
  the final maximum (0.4).

**PISA pre-filter (`pisa`, step 02).** PISA reports `interface_area` as the
solvent-accessible area buried *per side* (half of the total buried surface),
so 800 Å² means about 1,600 Å² buried in total, above the typical range of
crystal-packing contacts. `p_value` is PISA's ΔiG P-value: the probability of
finding a more hydrophobic patch of the same size by chance; values below 0.5
suggest an interaction-specific interface. PISA-lite does not provide the
complexation significance score (CSS) of classic PISA, so specificity relies
on area and P-value, with thresholds chosen for this pipeline.
These two values are only a pre-filter (they decide which structures are
downloaded); the final thresholds are set in `quality_filter` and can only be
stricter. Only the assemblies an entry has today are used
(`only_current_assemblies`): the PDBe dump can keep assemblies removed in later
revisions, and entries that RCSB no longer serves are listed as obsolete. PDBe
recomputes PISA when an entry is revised, so a rerun of step 01 downloads again
any file that PDBe has updated since the previous run.

**Scope and assemblies (`scope`, `assemblies`, steps 03, 07 and 13).** By
default only the first biological assembly of each entry is used
(`assemblies.selection: first`: its lowest current assembly ID, as in AlphaFold 3
and PINDER); `all` keeps every assembly that PISA analysed. `max_release_date`
keeps only the entries released up to a date, so that a later rerun gives the
same entries. Pairs with a chain that PDBTM annotates as a membrane protein
(RCSB entity annotations; `exclude_membrane`) are excluded. Antibody chains are
flagged from their SAbDab annotation in RCSB (`molA_is_antibody`,
`molB_is_antibody`) and excluded only with `exclude_antibodies: true`. The same
scope applies to the hard negatives, and step 07 logs how many pairs each
criterion removes.

**Quality filter (`quality_filter`, step 07).** Experimental method, R-free,
optional maximum resolution, minimum coverage of both chains, chain length
range, interface area range and P-value. The maximum interface area (5,400 Å²
per side, about the 99th percentile of the candidates) trims the extreme tail
so that a deep-learning model is not trained on extreme values; it is not an
artefact filter, since such large interfaces are mostly intertwined homodimers.
The final P-value (≤ 0.4) and the *gold*/*silver* tiers (P < 0.2) are
specificity classes, so that a model can be trained or evaluated on data of
different specificity; P is not a calibrated probability that an interface is
biological. Method and
resolution are often already restricted by the RCSB query (the reference dataset uses X-ray
structures up to 3.0 Å); cryo-EM entries are excluded by design because they
have no R-free. `coverage_metric` chooses between `coverage` (resolved / SEQRES
residues) and `internal_coverage` (resolved residues between the first and the
last resolved one): with the latter, disordered expression tags and tails are
not penalised, only internal gaps (in the full run it would rescue 91 % of the
chains that fail on `coverage`).
`min_chain_length` (at least 40 residues by default) keeps peptides out: 40
residues is the usual boundary between peptides and proteins in structural
datasets (PISCES lists, PINDER) and peptide binding is a different problem.
`max_chain_length` (1,000) is usually set by the model you plan to train. Homo/hetero pairs are defined by
the same entity or a sequence identity of at least `enrichment.homo_similarity`,
measured on a global alignment over the shorter chain, so a truncated construct
of the same protein counts as homo. Step 07 also writes descriptive statistics
of every variable to its log and histograms to `reports/07_quality_filter_plots.pdf`.

**Structure checks (`special_residues`, `contacts`, step 08).** Only the
polymer residues of each chain are considered (ligands, waters and free amino
acids that share the chain ID are ignored), in the first model, with heavy
atoms of occupancy > 0 and, where a residue has alternate locations, the
conformer with the highest mean occupancy. A chain is dropped if it has a polymer residue that is not one of
the 20 standard amino acids, not in `allowed` and not in `ignore` (terminal caps
ACE and NH2 by default). The default is conservative: only the selenium
residues used for phasing or naturally present (MSE, selenomethionine; SEC,
selenocysteine) are allowed. Modified residues could be mapped to their parent
amino acid, as AlphaFold and the PDB canonical sequence do, but the mapping is
not always clean, so this version excludes the chain instead. In
the full run the excluded chains (5.8 %) carry cysteine modifications that are
usually crystallisation artefacts (CSO, CSD, CSS, CME), catalytic or cofactor
lysines (KCX, LLP), GFP-like chromophores, N-terminal modifications (PCA),
post-translational modifications, D-amino acids or unknown residues (UNK); add
any of them to `allowed` if they are acceptable for your use. Every pair must
also have at least `min_residue_contacts` residue pairs with polymer heavy atoms
within `cutoff` Å (5 Å, the interface definition of AlphaFold 3; a KD-tree
search); this catches PISA interfaces that cannot
be mapped onto the structure. In the full run no pair failed it. Step 08 also
reports, for each chain, the unresolved residues within `interface_gap_window`
positions of an interface residue (`molA_interface_gaps`, `molB_interface_gaps`);
`max_interface_gaps` excludes pairs above a maximum (null by default: report
only). Pairs where a chain has more than `max_clash_fraction` (30 %) of its
atoms within `clash_distance` (1.7 Å) of the other chain are excluded, the
chain-removal rule of AlphaFold 3; every pair reports `n_clash_atoms` and
`max_clash_fraction`.

**Small interfaces (`prodigy`, step 08b).** Biological interfaces and crystal
contacts overlap most at small areas. Among the pairs with `interface_area`
below `max_area` (1,000 Å² per side), step 08b removes those that PRODIGY-cryst,
a contact-based classifier independent of PISA (the one PINDER uses), calls
crystal (XTAL); pairs it cannot score are kept. Its calls are in
`interim/08b_prodigy_band.csv` and its counts in `reports/08b_prodigy_filter.json`.
Above `max_area`, PISA is the only judge of whether an interface is biological.
`snakemake --cores 8 prodigy_cryst` (optional) annotates every final positive
in `reports/prodigy_cryst.csv` (`row_id`, `prodigy_class` BIO/XTAL,
`prodigy_p_bio`). Selenomethionine, written as HETATM records, is scored as
methionine, since PRODIGY-cryst would otherwise remove it with the ligands;
selenocysteine is not converted, so its pairs are not scored. In the reference
run, step 08b removed 4,385 of the 15,798 pairs below 1,000 Å² (28 %; 38 % of
the *silver* and 19 % of the *gold* ones), and PRODIGY-cryst calls 98.7 % of the
final positives biological (96.8 % of the test pairs); without step 08b it
called 95.4 %, and only 75.4 % of the interfaces of 800–1,000 Å² per side.
Above `max_area` it still calls 533 final positives crystal, 479 of them below
1,250 Å². PRODIGY-cryst is itself a
predictor: its authors report 74 % accuracy on the independent DC benchmark,
against 79 % for PISA (Elez *et al.*, 2018).

**Clustering (`clustering`, step 09).** Chains are clustered at the identity
levels of `thresholds` (95, 90, 60, 40 and 30 % by default). The clusters
define the test classes (40 %), the validation splits (60 and 40 %), the family
weights (40 %) and the hard-negative caps and families (95 and 40 %); step 09
stops if a level that the configuration reads is missing.
The default method, `connected_components`, runs one all-vs-all MMseqs2 search
and links two sequences when their identity is above the threshold with the
required alignment coverage in either search direction (the same rule as the
leakage audit); clusters are the connected components (single linkage). Two
sequences of different clusters then never have an MMseqs2 hit that meets the
linking criterion, so cluster-based splits do not leak homologues *as long as
the search sees every hit*: `max_seqs` must exceed the number of prefilter
candidates of the most populated family (antibody chains need thousands). The
price is that large families merge into big clusters (at 40 %, the largest
component of the full run holds 16 % of the sequences), and that two members
of a cluster can be less similar than the threshold through a chain of links.
The alternative `set_cover` (MMseqs2 `cluster` or `linclust`, used by the
working versions; despite the name, with `--cov-mode 1` MMseqs2 switches to
greedy incremental clustering, CD-HIT-like) only guarantees that every member
is similar to its representative: members of different clusters can still be
similar above the threshold: in the working versions it left 3 % of the
validation pairs, 6–8 % of the C2 pairs and about 3 % of the unseen-protein
test chains with a homologue in training at their design threshold
([MIGRATION.md](MIGRATION.md)). Other approaches can be preferable depending on the
question: greedy incremental clustering (CD-HIT), graph-partitioning methods
that balance split sizes while removing cross-split homology (GraphPart),
structure-based clustering (Foldseek, TM-score), which also catches remote
homologues below ~30 % sequence identity, or splits by domain classification
(CATH, SCOP, ECOD).

**Redundancy (`pruning`, step 10).** Many entries contain several copies of
the same complex (non-crystallographic symmetry, several assemblies), and close
homologues appear across entries.

1. *Identical sequence pairs* are capped, across the whole dataset, to the best
   `exact_K_homo` / `exact_K_heter` rows (ranked by coverage, P-value,
   interface area, R-free, resolution and solvation energy), taken from
   different entries where possible, so that the kept copies are not NCS
   copies of one crystal. In the reference run this cap removed 59 % of the
   rows (105,440 → 43,710). It looks at sequences, not at interfaces: when one
   entry has several interfaces of the same sequence pair (the different
   interfaces of a homotetramer, for instance), the kept rows can all be of one
   of them.
2. *Pair-cluster pruning* is optional and off by default (`cluster_pruning`).
   Every distinct sequence pair is kept, point mutants included, and large
   families are down-weighted instead by the family weight of step 12. When it
   is enabled, pair clusters (`pair_key`, `pair_cluster_95` by default: step 10
   names the clusters of step 09 without the `c`, e.g. `pair_cluster_c95`) with
   more than `activate_n` rows are reduced to `K_redund_*` rows when their
   interface areas are almost identical (range within `area_range_pct_thr` % of
   the mean) or to `K_var_*` rows spread over interface-area bins otherwise
   (52,085 → 43,164 with pair clusters at 98 %, in the run of 2 November 2026). This also removes
   mutant series; for a strictly non-redundant set, enable it and set every `K`
   to 1.

**Test set (`test_set`, step 11).** The external test set follows the three
classes of Park and Marcotte (2012): a protein is *seen* if its chain cluster
at `level` (40 % identity, the low-homology limit of AlphaFold 3) has a row in
the train/validation pool (every row not in the test set). Validation is drawn
from that pool later, so in a given split a seen protein may be only in
validation.

- **C3, unseen proteins** (`n_c3`): neither protein is seen. Chain clusters are
  held out, and every row that contains one leaves the training pool; a homo
  pair needs one held-out cluster, a hetero pair two.
- **C2, a new protein with a known partner** (`n_c2`): one protein is held out
  and the other keeps other partners in the pool (hetero by definition).
- **C1, new pairs of seen proteins** (`n_c1`): both proteins are seen, but the
  whole pair cluster at 40 % is removed, so no pair of the pool joins
  homologues of both chains.

Clusters are drawn at random (seeded), only among those with at most
`max_cluster_rows` rows, so that no large family is emptied from training;
there is no preference for rare families and no stratification by length or
area. Each pair cluster contributes one random representative. C1 and C3 aim
at `homo_frac` homomeric pairs (0.5, while about 76 % of the pairs after step 10
are homomeric; C2 is heteromeric by definition); when one kind runs out the
other fills the class, and step 11 warns if a class is short. Report metrics per class and separately for homo and hetero pairs. The
held-out clusters are written to `11_reserved_clusters.csv`, so that step 13
sends every negative with an unseen protein to the test set, and
`reports/11_test_set.json` gives the cost in training rows.

**Train/validation splits (`splits`, step 12).** The remaining pool is split
once per level of `schemes`, by pair cluster (by default at 60 and 40 %
identity: `split_c60`, `split_c40`), sending `val_frac` of the rows to
validation with the same homo/hetero proportion. Splitting by pair cluster
means that a protein A can appear in training and validation, but no training
pair joins a homologue of A with a homologue of its partner at the scheme's
identity: the validation sets measure generalisation to new pairs. At 40 % this
matches what AlphaFold 3 calls a low-homology interface, and most validation
pairs are of new proteins too: in the reference run, 77 % of the `split_c40`
validation pairs have neither chain cluster at 40 % in training, against 40 %
in `split_c60`. Clusters stay whole: a cluster that would overshoot the target
is still taken if that leaves the total closer to it, so one large family can
dominate a stratum (an influenza hemagglutinin family of 852 rows makes up two
thirds of the heteromeric validation rows of `split_c40` in the reference run). The schemes are
independent splits, so train one model per scheme. `validation_split`
(`split_c40`) selects the scheme exported as `validation.csv` and the one the
family weight is computed on: `family_weight_c40` = 1 / n^`weight_exponent`,
where n is the number of rows of the pair cluster at 40 % in training (0.5 by
default, so a family of n rows weighs √n; 1 gives every family the same
weight).

**Hard negatives (`hard_negatives`, step 13).** Two protein chains of the
same biological assembly with no PISA interface between them: negatives for
*direct contact*, not proof of non-interaction (in PSI-MI terms, chains of the
same complex without contact are still a physical association). Pairs of
identical sequences are excluded by default (a copy of a protein not touching
another copy usually contradicts a homo interface elsewhere in the complex). A
negative is dropped if its pair of chain clusters at `conflict_level` (40 % by
default) has a PISA interface in any downloaded assembly. This is a precaution: an
interface between homologues at 40 % does not prove that this pair interacts,
but it makes the negative doubtful (homologues above ~30–40 % identity that do
interact tend to do so in the same way, Aloy et al., *J. Mol. Biol.* 2003;
transferring an interaction itself needs much higher identity). It is also dropped if it fails the quality thresholds of
step 07 (except the interface ones), if a chain has forbidden residues (same
rules as step 08) or if the chains are closer than `min_distance` Å
(`max_distance`, if set, keeps only near-miss pairs). One negative is kept per
entry and sequence pair, at most `max_per_pair_cluster` per pair cluster
(`pair_cluster_key`, 95 % by default) across entries, so that complexes
deposited many times do not dominate, and at most `max_per_entry` per entry, so
that ribosomes and capsids do not either. Each negative gets a split for every
scheme of step 12 (`split_c60`, `split_c40`): test if a protein
appears only in the test set or in a held-out cluster of step 11, validation if a
protein appears only in validation, training otherwise; negatives of the same
pair cluster at `family_level` (40 %, the conflict level) all go to the most
restrictive of those sets, so a family is never split between training and
evaluation. Finally, every UniProt pair of a negative is looked up in IntAct
(PSICQUIC; `intact`): `intact_evidences` counts its evidences and
`intact_direct_evidences` the binary ones of type direct interaction (MI:0407 or
a PSI-MI descendant). Chains of one complex without contact are expected to
have association or physical association evidence; direct evidence contradicts
the negative, which is dropped (`drop_direct`). Answers are cached with their
date in `raw/intact_psicquic.json`. The step stops if a PISA file is missing,
since that would hide contacts.

**Leakage audit (step 14).** Chains are searched in both directions with the
identity, coverage and E-value criteria of step 09: every test class against
the training and validation positives and, separately, against the training
and validation hard negatives; the validation set of each scheme against its
training positives; and the hard negatives of the test set against the
training and validation positives and the non-test negatives. `reports/14_leakage_audit.json` gives the number
of leaks of each set at its design threshold (chain level for C3 and for the
unseen protein of C2; pair level for C2, C1, validation and hard negatives)
and `reports/14_leakage_by_threshold.tsv` the full picture at every threshold.
With connected components all of them should be zero (zero leaks *by
sequence*, with these criteria); `tests/check_outputs.py` fails otherwise. The audit uses the same
tool and criteria as the clustering, so it is a consistency check: it cannot
see remote homologues below the lowest threshold. It is therefore repeated
without the coverage rule, so that a single shared domain in multidomain
proteins also counts (`design_threshold_leaks_without_coverage_rule`); that
second audit is reported, not enforced. The optional independent control
(`analysis/leakage_control.py`) audits the same sets with global-alignment
identities instead; see *Clustering and leakage* in the README.
