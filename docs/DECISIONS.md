# Curation decisions

This document records the curation decisions of PPI_PDB_Scraper: what the
pipeline does, why, which published resources do the same or otherwise, what
each decision costs, and how to change it. Most thresholds are user parameters
in `config/config.yaml`; the values given here are the defaults. The README
describes the pipeline step by step.

Each entry has the same fields:

- **Decision**: what the pipeline does.
- **Reason**: why. When the value is a matter of what the dataset is for, the
  entry says it is the user's choice and which precedent the default follows.
- **Precedent**: published datasets and methods that do the same or otherwise
  (keys in brackets, full references at the end).
- **Effect**: numbers from the full run of 10 November 2026 (164,000 PDB IDs,
  X-ray structures at 3.0 Å or better). *Previous behaviour* gives the numbers
  of the run of 2 November 2026, before the changes of this version.
- **Bias / limitation**: what the decision leaves out or distorts.
- **Configuration**: the parameter in `config/config.yaml` and its default.
- **Status**: `implemented` (in the code and in the full run); `planned`;
  `under discussion` (the direction is set, the design is not).

Blocks: A scope · F sources · I interface · Q structural quality ·
R redundancy · H homology and clustering · T test set · S train/validation
splits · N hard negatives · U leakage audit · P provenance and reproducibility.

## A · Scope

### A-01 · Scope set by an RCSB search, stored with its date
- **Decision:** The scope of the dataset (method, resolution, dates…) is set by an RCSB Advanced Search whose exported ID lists are read from `input_dir`. The query (JSON file or URL) and the date it was run are stored in the configuration and copied, with the SHA-256 of every ID list, to `reports/provenance.json`.
- **Reason:** The ID list cannot be rebuilt without the query and its date. Each result needs "the exact parameters and inputs that were used" [Sandve2013] and data need "detailed provenance" (FAIR R1.2) [Wilkinson2016]. An RCSB query can be stored as a URL-encoded link [Rose2026].
- **Precedent:** Reference datasets publish their scope and archive date: AlphaFold 2 (PDB downloaded 28 August 2019) [Jumper2021]; ProteinMPNN (PDB "as of Aug 02, 2021", X-ray or cryo-EM better than 3.5 Å) [Dauparas2022]; PPIRef (PDB "as of June 20, 2023") [Bushuiev2024]; PINDER (RCSB snapshot of 29 January 2024, release `2024-02`) [Kovtun2024].
- **Effect:** 164,000 IDs, from a search for X-ray entries at resolution ≤ 3.0 Å; `reports/provenance.json` holds the SHA-256 of the 17 ID lists. The query of this run was not stored (`scope.input_query` is empty).
- **Configuration:** `input_dir` (`pdb_ids`); `scope.input_query` and `scope.input_query_date` (new).
- **Status:** `implemented`

### A-02 · Four-character PDB IDs only
- **Decision:** Only 4-character PDB IDs are accepted. Entity or assembly forms such as `4HHB_1` are reduced to the entry ID; any other token is logged in `01_rejected_tokens.tsv` and dropped.
- **Reason:** It is the ID format in use today and the one the PDBe PISA API takes ("The 4-character Protein Data Bank (PDB) accession code") [Leines2026].
- **Effect:** 0 tokens rejected.
- **Bias / limitation:** On 21 July 2027 the wwPDB moves to 12-character IDs (`pdb_` + 8 characters; existing IDs gain the prefix `pdb_0000`) [wwPDBextID]. Version 2 of the pipeline will accept both formats (README, Roadmap).
- **Status:** `implemented`

### A-03 · Optional maximum release date
- **Decision:** A maximum release date can be set; entries released after it are left out. By default there is no cut-off. The initial release date of each entry is stored in the entry metadata (step 04), and the cut-off applies to positives (step 07) and hard negatives (step 13).
- **Reason:** The value is the user's choice; the default keeps every entry available on the download day.
- **Precedent:** Reference models fix a maximum release date: 30 April 2018 for AlphaFold 2 and AlphaFold-Multimer [Jumper2021; Evans2021], 30 September 2021 for AlphaFold 3 [Abramson2024]. ProteinMPNN, PPIRef and PINDER date the PDB snapshot they use (2 August 2021, 20 June 2023, 29 January 2024) [Dauparas2022; Bushuiev2024; Kovtun2024].
- **Effect:** No cut-off; the final positives were released between 31 March 1995 and 11 February 2026.
- **Bias / limitation:** Without a cut-off the content depends on the download day (P-02), and a comparison with a model trained up to a given date needs the cut-off to be set.
- **Configuration:** `scope.max_release_date` (null = no cut-off).
- **Status:** `implemented`

### A-04 · X-ray diffraction only
- **Decision:** Only entries whose experimental method is X-ray diffraction are kept. An entry solved by several methods (e.g. joint X-ray and neutron refinement) counts as X-ray if one of its methods is X-ray diffraction (F-05); the filter is applied to positives (step 07) and hard negatives (step 13).
- **Reason:** The method is the user's choice. The default keeps X-ray diffraction because the quality criteria of the pipeline are crystallographic (R-free, Q-01).
- **Precedent:** X-ray only: the PINDER test set (≤ 3.5 Å) [Kovtun2024] and PPI4DOCK (high-resolution X-ray heterodimers) [YuGuerois2016]. Most recent datasets also admit cryo-EM at ≤ 3.5 Å: ProteinMPNN, DIPS-Plus, PPIRef [Dauparas2022; Morehead2023; Bushuiev2024].
- **Effect:** No pair fails the method filter, because the RCSB search already selected X-ray entries. The 74 pairs of 42 joint X-ray and neutron entries count as X-ray (F-05), including 21 from entries that list the neutron method first, which the previous behaviour excluded.
- **Bias / limitation:** Leaving out cryo-EM biases the dataset towards homomers: since 2010, 69 % of cryo-EM structures have been heteromeric, against 13–15 % of X-ray structures [XuDunbrack2019]. 79.6 % of the train/validation pairs are homomeric (31,636 of 39,747). EM accounted for 41 % of the entries released in 2025 (7,233 of 17,599) [RCSBstats]. Admitting cryo-EM needs a quality criterion other than R-free, such as map–model FSC or Q-score, which the PDB already reports [Lawson2021].
- **Configuration:** `quality_filter.exptl_method` (`"X-RAY DIFFRACTION"`). Changing it alone does not admit methods without R-free: step 07 always applies `max_r_free` and drops entries without R-free (Q-02).
- **Status:** `implemented`

### A-05 · No resolution cut-off of its own
- **Decision:** By default the pipeline applies no resolution cut-off; the cut-off is set in the RCSB search (≤ 3.0 Å in the full run).
- **Reason:** The value is the user's choice. The 3.0 Å of the full run is stricter than usual: 3.5 Å is "a standard threshold in the field" [Morehead2023].
- **Precedent:** 3.5 Å: ProteinMPNN, DIPS-Plus, PPIRef, the PINDER test set and the Dockground templates [Dauparas2022; Morehead2023; Bushuiev2024; Kovtun2024; Collins2022]. 6.0 Å: the Dockground bound set [Collins2022]. 9 Å: AlphaFold-Multimer and AlphaFold 3 [Evans2021; Abramson2024].
- **Effect:** Largest resolution value among the 155,899 enriched pairs: 3.02 Å.
- **Bias / limitation:** Until the query is stored (A-01), the cut-off of a run is recorded only outside the repository.
- **Configuration:** `quality_filter.max_resolution` (null = no cut-off).
- **Status:** `implemented`

### A-06 · Membrane proteins excluded by default
- **Decision:** Pairs with a chain annotated as a membrane protein by any source of a configurable list are excluded, from positives and hard negatives. The annotations are those RCSB attaches to each polymer entity (PDBTM, OPM, MemProtMD, mpstruc); the default list is `[PDBTM]`, and the user can change or empty it.
- **Reason:** Membrane proteins should be left out by default; the value is the user's choice.
- **Precedent:** Dockground annotates membrane proteins in its bound set (15.4 % of 667,331 interfaces) and keeps a separate membrane set, because generic docking methods "are not optimal for the membrane proteins" [Collins2022; Kotthoff2022].
- **Effect:** 3,032 of the 155,899 enriched pairs (1.9 %) have a chain annotated by PDBTM and are excluded; 1,935 of them fail no other criterion of step 07. Previous behaviour: none excluded. Chains annotated only by the other sources stay: in the final data, 190 train/validation pairs, 6 test pairs and 10 hard negatives have a chain annotated by OPM, mpstruc or MemProtMD but not by PDBTM.
- **Bias / limitation:** The general features of transmembrane interfaces "do not differ much from those of soluble proteins" [Duarte2013]; the exclusion removes a whole class of interfaces, and its reach depends on the coverage of the annotation source.
- **Configuration:** `scope.exclude_membrane` (`[PDBTM]`; `[]` keeps them).
- **Status:** `implemented`

### A-07 · No limit on complex size
- **Decision:** Large complexes are not excluded: there is no limit on the number of chains per assembly.
- **Reason:** The value is the user's choice; it can also be set in the RCSB search that defines the scope (A-01).
- **Precedent:** Size limits appear where the model reads the whole complex: AlphaFold 3 (≤ 300 chains in training, cropped to the 20 closest chains) [Abramson2024], ProteinMPNN (< 10,000 residues) [Dauparas2022], DIPS-Plus (≤ 17,500 atoms) [Morehead2023]. PINDER and Dockground set no limit [Kovtun2024; Collins2022].
- **Bias / limitation:** Large complexes can be over-represented: 171 of the 492 hard negatives (35 %) come from entries with "proteasome" in the title and 116 (24 %) from entries whose title names a ribosome or a ribosomal subunit; in the test set, 156 of the 289 (54 %) come from one of the two. The exact-pair cap (R-01) and the caps on negatives (N-11) limit this, but do not remove it.
- **Status:** `implemented`

### A-08 · Antibodies flagged; exclusion is optional
- **Decision:** Antibody chains are flagged in the columns `molA_is_antibody` and `molB_is_antibody`, from the SAbDab annotation that RCSB attaches to each polymer entity. The user decides whether to exclude them; by default they are kept.
- **Reason:** Antibodies must be explicit in the dataset; excluding them is the user's choice.
- **Precedent:** AlphaFold 3 reports protein–protein metrics with and without antibodies [Abramson2024]; QSalign-HET excludes immunoglobulins [DeyLevy2021]; PINDER marks complexes with antibodies (`contains_antibody`, from SAbDab) [Kovtun2024; Schneider2022]; Docking Benchmark 5 keeps antigen–antibody complexes as a separate category [Vreven2015].
- **Effect:** 3,670 of the 155,899 enriched pairs have an antibody chain; so do 3,358 of the 118,044 pairs that pass step 07 (2.8 %) and 2,679 of the 40,497 final positives (6.6 %). None is excluded.
- **Bias / limitation:** The wwPDB notes that quaternary-structure programs such as PISA "do not work with every case, such as antibodies" [wwPDBproc]. Antibody chains are the most populated sequence family (see `clustering.max_seqs`). The flag is only as complete as its source.
- **Configuration:** `scope.exclude_antibodies` (false).
- **Status:** `implemented`

### A-09 · Protein–protein interfaces only
- **Decision:** Only PISA interfaces between two protein molecules are kept. A protein–protein interface is kept even if the assembly also contains DNA, RNA or ligands.
- **Reason:** The dataset covers protein–protein interfaces; other molecules in the assembly are not a reason to drop one.
- **Precedent:** Dockground keeps protein–protein interfaces with DNA/RNA (18.7 %) or a ligand (67.8 %) at the interface, and annotates them [Collins2022].
- **Effect:** 46,262 of the 564,347 PISA interfaces dropped as not protein–protein.
- **Bias / limitation:** Interfaces with a nucleic acid or a ligand at the interface are not annotated.
- **Status:** `implemented`

### A-10 · Homo- and hetero-oligomeric pairs treated alike
- **Decision:** Homo- and hetero-oligomeric pairs go through the same steps; homomeric pairs are neither excluded nor subsampled apart. The test classes aim at half homo and half hetero pairs (C2 is hetero by definition, T-09) and splits are stratified by homo/hetero (S-04). Metrics are reported overall and separately for homo and hetero pairs.
- **Reason:** Both kinds stay in one dataset; metrics are reported overall and per kind because homo- and hetero-oligomeric interfaces differ in amino acid composition and residue–residue preferences [OfranRost2003; JonesThornton1996].
- **Precedent:** PINDER keeps both and prefers heterodimers for its test set [Kovtun2024]; AlphaFold-Multimer and AlphaFold 3 train on both [Evans2021; Abramson2024]; PPI4DOCK keeps heterodimers only [YuGuerois2016].
- **Effect:** 79.6 % of the train/validation pairs are homomeric (31,636 of 39,747); 246 of the 750 test pairs (C2 pairs are heteromeric by definition, and C1 has 96 homomeric and 54 heteromeric pairs, T-05, T-09). Homomeric follows `enrichment.homo_similarity`, so it includes 334 train/validation pairs and 1 test pair of two different entities with aligned identity ≥ 0.80, such as the p66–p51 heterodimer of HIV-1 reverse transcriptase in 6bsh (558 and 441 residues, identical where they align). Conversely, 897 of the 8,111 heteromeric train/validation pairs and 45 of the 504 heteromeric test pairs join two chains with a common UniProt accession; 418 of the former come from entries whose title names hemagglutinin, whose HA1 and HA2 chains are cleaved from one precursor.
- **Bias / limitation:** The PDB is dominated by homomers (14.2 % of entries are heteromers) [XuDunbrack2019], and A-04 adds to this. PINDER notes "the inherent bias towards homodimers" of such data [Kovtun2024].
- **Configuration:** `enrichment.homo_similarity` (0.80: same entity or aligned identity ≥ 0.80 counts as homo); `test_set.homo_frac` (0.5).
- **Status:** `implemented`

### A-11 · Interfaces with a symmetry copy are positives
- **Decision:** An interface between a chain and its symmetry copy in the assembly (e.g. `A`–`A-2`) is a positive like any other.
- **Reason:** The interface between a chain and its symmetry copy is often the biological one. About 14 % of PDB entries need chains from a neighbouring unit cell to build an annotated assembly, and isologous interfaces are exactly C2 symmetric when one monomer is a symmetry copy of the other [XuDunbrack2019]. The common dimer of the PALP family is in the asymmetric unit of only 42 of 109 entries (39 %) [XuDunbrack2011].
- **Precedent:** ProtCID builds symmetry copies [XuDunbrack2011]; AlphaFold 3 and PINDER expand the assembly and treat each copy as a chain [Abramson2024; Kovtun2024].
- **Effect:** 77,440 of the 204,075 candidate pairs (37.9 %) involve a symmetry copy.
- **Bias / limitation:** Over the whole crystal lattice, contacts through crystallographic operators are biological less often than non-crystallographic ones (37.3 % of the latter are biological) [Baskaran2014]. A symmetry-copy positive is only as good as the annotated assembly (I-01).
- **Status:** `implemented`

### A-12 · First biological assembly by default
- **Decision:** By default only the first biological assembly of each entry is used; the alternative is every assembly that PISA analysed. The first assembly is the lowest current assembly ID of the entry (RCSB Data API); the selection is made in step 03.
- **Reason:** The value is the user's choice; the default follows AlphaFold 3 and PINDER, which use the first biological assembly.
- **Precedent:** AlphaFold 3: "The first bioassembly is then expanded" [Abramson2024]; PINDER: "the first biological assembly was generated" [Kovtun2024]. For `all`, AlphaFold-Multimer samples one assembly uniformly from all available ones [Evans2021].
- **Effect:** 155,899 of the 204,075 candidate pairs (76.4 %) are in the first assembly, from 66,807 of the 69,580 entries with candidate pairs; the other 2,773 entries (4.0 %) have candidates only in other assemblies and are left out. Previous behaviour (all assemblies): 15,920 entries (22.9 %) had pairs in more than one assembly.
- **Bias / limitation:** 30 % of PDB entries have more than one assembly: in 83 % of them they are essentially the same, in 17 % they are hypothetical assemblies of different size [XuDunbrack2019]. Published estimates of the accuracy of the first assembly are around 80 % [XuDunbrack2019]. Interfaces present only in other assemblies are lost.
- **Configuration:** `assemblies.selection` (`first`; alternative `all`).
- **Status:** `implemented`

### A-13 · Obsolete entries excluded
- **Decision:** Entries that RCSB no longer serves are marked `obsolete` in the PISA manifest and left out, with all their assemblies.
- **Reason:** An obsolete entry is no longer part of the archive: the wwPDB makes an entry obsolete at the request of those responsible for it (better data or an improved interpretation) or when its publication is retracted [wwPDBobs].
- **Effect:** 32 entries (59 assemblies).
- **Configuration:** `pisa.only_current_assemblies` (true), which also detects obsolete entries.
- **Status:** `implemented`

### A-14 · Entries without PISA data excluded
- **Decision:** IDs with no PISA-lite data are logged (`no_pisa_data` in the PISA manifest) and left out. Step 01 stops if they exceed 20 % of the IDs.
- **Reason:** These entries are better left out. PDBe computes PISA-lite for every annotated assembly and updates it weekly [Leines2026; PDBe2023], so the gaps are technical, not a scientific selection. A share above 20 % points to a wrong ID list or a changed listing format.
- **Effect:** 276 of 164,000 IDs (0.17 %).
- **Configuration:** `pisa.max_missing_fraction` (0.2).
- **Status:** `implemented`

## F · Sources and conflicts between them

### F-01 · Interfaces taken from PDBe PISA-lite, not recomputed
- **Decision:** Interfaces and their energetics come from the PDBe PISA-lite bulk files and are not recomputed.
- **Reason:** PISA [KrissinelHenrick2007] is the standard tool to analyse interfaces in crystal structures. PISA-lite is the PDBe run of PISA 2.x on the annotated assemblies of the PDB, updated weekly with the archive; it does not predict assemblies [Leines2026; PDBe2023]. PISA and EPPIC give the same call for 88 % of PDB interfaces [Baskaran2014], and reach almost the same precision on assemblies (85 % for EPPIC) [Bliven2018].
- **Effect:** 564,347 interfaces parsed from 245,757 assemblies of 163,692 entries; 204,075 candidate pairs after the pre-filter (I-03, I-04).
- **Bias / limitation:** "PISA results are computational predictions, and any conclusions must be validated" [Leines2026]; the only independent cross-check is PRODIGY-cryst, applied to the interfaces below 1,000 Å² per side (I-02). PISA-lite does not provide the complexation significance score (CSS) of classic PISA.
- **Configuration:** `pisa.base_url`.
- **Status:** `implemented`

### F-02 · Current RCSB assemblies prevail over the PDBe dump
- **Decision:** Only the assemblies an entry has today according to the RCSB Data API are used. PISA-lite files of assemblies removed in later revisions (`stale`) are not downloaded.
- **Reason:** The PDBe dump can keep files of assemblies removed in later revisions of an entry (e.g. 5u1j). RCSB and PDBe distribute the same wwPDB archive, which versions every entry [wwPDBproc], so the current version prevails.
- **Effect:** 27 stale assemblies (20 entries) left out. In the run of 2 November 2026, PISA's number of macromolecules equalled the polymer instances of the RCSB assembly with the same ID in 92,625 of 92,755 comparable assemblies (99.86 %); the 130 differences are short or mostly unmodelled polymers that PISA counts as ligands (e.g. a 4-residue peptide in 1ztz, a trinucleotide in 1rrc), and 273 PISA files gave no macromolecule count.
- **Bias / limitation:** The pipeline assumes that PISA-lite assembly IDs are those of the PDB. PISA-lite works only on the annotated assemblies [PDBe2023], but the PDBe protocol describes the API `assemblyid` as "PISA-assigned" [Leines2026]. A check showed that both describe the same assemblies (Effect).
- **Configuration:** `pisa.only_current_assemblies` (true).
- **Status:** `implemented`

### F-03 · Coordinates from the RCSB assembly files
- **Decision:** Coordinates are read from the RCSB biological-assembly files (`{id}-assembly{n}.cif.gz`), not from the asymmetric unit. Only the assemblies that PISA analysed and that hold a candidate pair are downloaded.
- **Reason:** PISA-lite analyses the assemblies annotated in the PDB [Leines2026; PDBe2023], so the RCSB file of the same assembly holds the chains and symmetry copies that its interfaces refer to; in the run of 2 November 2026 the number of macromolecules agreed in 99.86 % of the comparable assemblies (F-02).
- **Effect:** 66,807 assembly files, the first assembly of each entry (A-12).
- **Configuration:** `rcsb.download_url`.
- **Status:** `implemented`

### F-04 · PISA chain IDs matched to `auth_asym_id`
- **Decision:** A PISA chain ID is matched to the `auth_asym_id` of the assembly file, symmetry copies included (e.g. `A-2`).
- **Reason:** PISA reports author chain IDs. The previous version (Seq8) matched by `label_asym_id`, and 5.7 % of the pairs got the sequence of another chain. It is an identifier mapping, not a scientific choice.
- **Effect:** Step 08 finds atomic contacts for every PISA pair (0 pairs without contacts), which supports the mapping; a wrong copy that also touches the partner would pass too. In the final positives the median is 69.5 residue contacts per 1,000 Å² of PISA area, and 27 pairs have fewer than 40 (the lowest, 6cug A–D, has 11 contacts for 999 Å²); they were not reviewed.
- **Status:** `implemented`

### F-05 · Method, resolution and R-free from the matching refinement record
- **Decision:** Method, resolution and R-free come from the RCSB Data API. If one of the entry's methods is `quality_filter.exptl_method`, that is its method, and resolution and R-free come from the `_refine` record whose `pdbx_refine_id` is that method; otherwise the first method and the first `_refine` record are used. All the methods of the entry are kept in `exptl_methods`.
- **Reason:** Minor correction. The PDBx/mmCIF dictionary allows several methods per entry, and `_refine.pdbx_refine_id` "can be used to distinguish the results of joint refinements" [mmCIFdic].
- **Precedent:** RCSB computes `resolution_combined` from `refine.ls_d_res_high` of the record whose `pdbx_refine_id` is the method [RCSBschema].
- **Effect:** 48 entries (94 pairs) have more than one method. The 74 pairs of joint X-ray and neutron entries count as X-ray, including the 21 of entries that list the neutron method first, which the previous behaviour excluded (A-04); 18 of those pass step 07. Previous behaviour: mixed-method entries with X-ray listed first took whichever `_refine` record came first (pairs affected not counted).
- **Configuration:** `quality_filter.exptl_method` (`"X-RAY DIFFRACTION"`).
- **Status:** `implemented`

### F-06 · Canonical entity sequence
- **Decision:** The sequence of each chain is the canonical sequence of its entity (`_entity_poly.pdbx_seq_one_letter_code_can`), with terminal `X` removed. For the homo/hetero alignment of step 06 only, U and O are scored as C and K, and any other letter outside BLOSUM62 as `X`.
- **Reason:** In the canonical sequence, non-standard residues are represented "by the codes of their parents if parent is specified [...] or by letter 'X' if parent is not specified" [mmCIFdic]. Terminal `X` are caps such as ACE or NH2; removing them gives capped and uncapped constructs of a protein the same sequence.
- **Precedent:** AlphaFold 3 maps modified residues to their parent for clustering and treats residues without a parent as unknown [Abramson2024].
- **Bias / limitation:** Any terminal `X` is removed, not only caps. Chains with non-standard residues are handled in step 08 (Q-08).
- **Status:** `implemented`

### F-07 · PISA defines the interface; atomic contacts are checked
- **Decision:** PISA defines the interface. Step 08 checks that each PISA pair has at least one residue pair with heavy atoms within `contacts.cutoff` in the assembly file (5.0 Å, I-13); a pair without contacts is dropped.
- **Reason:** The PISA definition is respected; the check guards against chain-mapping errors (F-04). When a model is trained, the data are preprocessed with atomic tests that check the number and nature of the contacts to build the supervision masks.
- **Precedent:** Datasets define an interaction by atomic contact in the coordinates: AlphaFold 3, heavy atoms < 5 Å [Abramson2024]; PINDER, backbone atoms ≤ 10 Å [Kovtun2024]; PPIRef, heavy atoms ≤ 10 Å [Bushuiev2024].
- **Effect:** 0 pairs dropped. Median of 104 residue contacts per pair; Spearman correlation of 0.96 between PISA interface area and number of residue contacts.
- **Configuration:** `contacts.cutoff` (5.0 Å, see I-13), `contacts.min_residue_contacts` (1).
- **Status:** `implemented`

### F-08 · PISA metrics copied; interface residues counted
- **Decision:** PISA metrics are copied as they are, except the number of interface residues, which is computed per side as the number of residues with buried surface area > 0 in PISA (`molA_interface_residues`, `molB_interface_residues`). The raw PISA field stays in the intermediate tables (`number_interface_residues_raw`) and is not exported.
- **Reason:** Minimal change to a misleading field. In PISA-lite, `number_interface_residues` is the length of the residue list of one molecule (the whole chain), not the number of interface residues. Data must carry "accurate and relevant attributes" (FAIR R1) [Wilkinson2016], and the PDBe protocol reads the field literally ("201 interface residues" for an interface of 766.45 Å²) [Leines2026].
- **Effect:** Median of 42 interface residues per side (interquartile range 33–58) in the 109,825 pairs that pass step 08, about 36 Å² of buried area per interface residue. Previous behaviour: the field held the length of a whole chain.
- **Status:** `implemented`

### F-09 · UniProt accessions of each chain, from SIFTS through RCSB
- **Decision:** Each chain carries the UniProt accessions of its entity (`molA_uniprot_ids`, `molB_uniprot_ids`), from the SIFTS mapping that the RCSB Data API serves per polymer entity. They are annotations; no filter uses them.
- **Reason:** UniProt accessions link each pair to interaction databases such as IntAct, which the cross-check of the hard negatives needs (N-03).
- **Precedent:** Negatome maps its pairs to UniProt and removes those that IntAct reports as interacting [Smialowski2010].
- **Effect:** 91,076 of the 103,683 polymer entities of the entry metadata (87.8 %) have a UniProt accession, and 91.1 % of the chains of the final positives.
- **Bias / limitation:** Entities without a UniProt mapping (for instance synthetic or designed proteins) have no accession, and fusion proteins have several.
- **Status:** `implemented`

## I · Interface: biological versus crystal contact, and PISA thresholds

### I-01 · Accept the PDB-annotated assemblies as biological
- **Decision:** The biological assemblies annotated in the PDB (selection in A-12) are taken as biological units. The type of definition of each assembly (author, software or both) is neither recorded nor used.
- **Reason:** The PDB annotation is used as given, and the limitation is declared.
- **Precedent:** AlphaFold 3 and PINDER also build their complexes from the annotated first biological assembly [Abramson2024; Kovtun2024].
- **Effect:** Not measured.
- **Bias / limitation:** 10–15 % of PDB entries may have incorrect or ambiguous assembly annotations [XuDunbrack2023]; the lower bound for author annotations is about 7 % [Baskaran2014]; the accuracy of the first assembly is commonly around 80 % [XuDunbrack2019]. Software-defined assemblies inherit the error of the predictor: PISA mispredicts 25 % of heterodimers [DeyLevy2021]. Because the definition type is not recorded, software-only assemblies cannot be told apart.
- **Status:** `implemented`

### I-02 · Independent biological/crystal check: PRODIGY-cryst on small interfaces
- **Decision:** Step 08b runs PRODIGY-cryst, a contact-based classifier independent of PISA [JimenezGarcia2019; Elez2018], on every pair below 1,000 Å² per side and removes those it calls crystal (XTAL); pairs it cannot score are kept. Above that area, PISA-lite is the only judge; it does not provide the complexation significance score (CSS) of classic PISA. The optional rule `prodigy_cryst` annotates every final positive (`prodigy_class`, `prodigy_p_bio`, in `reports/prodigy_cryst.csv`).
- **Reason:** The biological character of an interface rests on the PDB annotation (I-01) and on PISA-lite metrics, which are computational predictions that PDBe asks users to validate [Leines2026]. Biological interfaces and crystal contacts overlap most at small areas (I-03): without the filter, PRODIGY-cryst called crystal 24.6 % of the final positives of 800–1,000 Å² per side, 7.1 % of those of 1,000–1,250 Å² and at most 0.8 % above. Below 1,000 Å², a pair is kept only if PRODIGY-cryst also calls it biological.
- **Precedent:** PINDER annotates every interface with PRODIGY-cryst and requires a physiological call in its test set [Kovtun2024]. QSbio gives an error estimate for each assembly so that users can choose the highest-confidence ones [DeyLevy2021]. EPPIC classifies interfaces by evolutionary analysis [Duarte2012].
- **Effect:** Step 08b scored all 15,798 pairs below 1,000 Å² of the 109,825 that pass step 08 and removed the 4,385 it calls crystal (27.8 %): 3,515 homomeric and 870 heteromeric, 1,712 *gold* and 2,673 *silver*; 105,440 rows go on to step 09. It calls crystal 19.4 % of the *gold* and 38.4 % of the *silver* pairs below 1,000 Å². Selenomethionine, written as HETATM records, is scored as methionine; before this conversion, the preparation of PRODIGY-cryst removed it with the other HETATM records, selenomethionine chains lost their methionines, and the filter removed 102 more pairs (4,487). In the final positives, the optional annotation scored 40,479 of 40,497 (18 failures: 16 pairs with a selenocysteine, 1 with a position modelled both as selenomethionine and methionine, 5by8, and 1 with a duplicated atom in the file, 9ebl) and calls 98.7 % of them biological: 99.5 % of the *gold* and 96.2 % of the *silver* pairs; train/validation 98.7 %, test 96.8 % (C1 97.3 %, C2 94.3 %, C3 99.0 %). Previous behaviour (no filter, same run up to step 08, same scoring): 95.4 % of 41,968 scored final positives called biological, 75.4 % of those of 800–1,000 Å² (1,431 pairs called crystal, 78 in the test set), 87.2 % of the test pairs.
- **Bias / limitation:** PRODIGY-cryst is itself a predictor: its authors report an accuracy of 92 % on the MANY set it was trained on (cross-validation), but 74 % on the independent DC benchmark of 80 biological and 80 crystal interfaces, against 79 % for PISA and 81 % for EPPIC [Elez2018]. The filter therefore also removes biological interfaces that it misclassifies and keeps crystal contacts that it calls biological: it adds an independent opinion where PISA-lite is least reliable, not a ground truth. It removes more *silver* than *gold* pairs. Pairs at or above 1,000 Å² are not filtered: the annotation calls 533 final positives crystal, 479 of them at 1,000–1,250 Å². Selenocysteine, written as ATOM records, is not converted, and pairs with it are kept unscored. Crystal dimers cannot be distinguished from real dimers by interface size alone [Bahadur2004]; PISA and EPPIC give the same call for 88 % of PDB interfaces [Baskaran2014].
- **Configuration:** `prodigy.filter` (true), `prodigy.max_area` (1000). PRODIGY-cryst runs in its own environment (`envs/prodigy_cryst.yml`: Python 3.8 and scikit-learn 0.22); with `prodigy.filter: false`, step 09 reads the output of step 08.
- **Status:** `implemented`

### I-03 · Pre-filter: interface area of at least 800 Å² per side
- **Decision:** In step 02, a protein–protein interface becomes a candidate pair only if PISA's `interface_area`, the area buried per side [PISA-iarea], is at least 800 Å².
- **Reason:** The value is the user's choice. The default, 800 Å² per side (about 1,600 Å² in total), is the mean size of the standard-size interfaces of Lo Conte et al. [LoConte1999] and lies well above the average crystal contact (570 Å² in total, about 285 Å² per side) [Janin2007].
- **Precedent:** NOXclass: crystal contacts bury 436 Å² per side on average, and a single cut at 650 Å² per side misclassifies about 7 % of interfaces [Zhu2006]. DIPS-Plus requires at least 500 Å² of buried surface [Morehead2023]. The DCbio/DCxtal benchmark covers 800–2,000 Å² per side [Capitani2016].
- **Effect:** 278,374 of the 564,347 PISA interfaces fall below the cut (counted after removing the interfaces that are not protein–protein).
- **Bias / limitation:** Area is necessary but not sufficient: biological interfaces and crystal contacts overlap between about 800 and 2,000 Å² per side [Duarte2012]. Below 1,000 Å², step 08b removes the pairs that PRODIGY-cryst calls crystal (I-02). Small interfaces are lost, such as those of electron-transfer complexes (450–600 Å² per side) and about half of Lo Conte's standard-size interfaces [Janin2007; LoConte1999]. The dataset represents stable interfaces of standard size or larger, not weak transient interactions.
- **Configuration:** `pisa.min_interface_area` (800)
- **Status:** `implemented`

### I-04 · Pre-filter: PISA ΔiG P-value of at most 0.5
- **Decision:** In step 02, an interface becomes a candidate pair only if its ΔiG P-value is at most 0.5.
- **Reason:** PISA's own definition. P is the probability of a lower ΔiG when interface atoms are picked at random from the protein surface up to the observed area. P < 0.5 means a more hydrophobic interface than expected, possibly interaction-specific; P > 0.5 suggests a packing artefact; 0.5 is the neutral point [PISA-ipv].
- **Effect:** 35,636 interfaces that pass the area cut are removed; 204,075 candidate pairs remain.
- **Bias / limitation:** P reflects solvation only: ΔiG does not include hydrogen bonds or salt bridges across the interface [PISA-idg]. Polar interfaces, such as those of weak homodimers [Nooren2003b], are disfavoured.
- **Configuration:** `pisa.max_p_value` (0.5)
- **Status:** `implemented`

### I-05 · The pre-filter decides which structures are downloaded
- **Decision:** Only assemblies with at least one candidate pair after the pre-filter (I-03, I-04) are downloaded, so the final thresholds of step 07 can only be stricter than the pre-filter.
- **Reason:** The pre-filter uses permissive values (PISA's neutral P = 0.5 and 800 Å²). It limits the download without fixing the final selection, and the final thresholds can be tightened without downloading again.
- **Effect:** 66,807 assembly files downloaded, the first assembly of each of the 66,807 entries with candidate pairs in it (A-12), out of 69,580 candidate entries.
- **Bias / limitation:** Relaxing a final threshold below the pre-filter requires re-running steps 02–04. Assemblies without any pre-filtered interface are not downloaded, which also restricts the hard negatives (N-03). The pre-filter keeps area ≥ 800 Å² and step 07 keeps area > 800 Å²; the 2 interfaces with exactly 800.00 Å² pass the first and fail the second.
- **Configuration:** `pisa.min_interface_area` (800), `pisa.max_p_value` (0.5)
- **Status:** `implemented`

### I-06 · Final interface area: above 800 and up to 5,400 Å² per side
- **Decision:** Step 07 keeps pairs with 800 < `interface_area` ≤ 5,400 Å² per side.
- **Reason:** The lower limit is that of I-03. The upper limit is pragmatic: it keeps extreme percentiles out of the data used to train a deep-learning model. 5,400 Å² is about the 99th percentile of the 204,075 candidates (5,352 Å²). It is not an artefact filter.
- **Precedent:** No published upper limit was found; DIPS-Plus, PINDER and AlphaFold 3 set none [Morehead2023; Kovtun2024; Abramson2024].
- **Effect:** 1,364 rows fail the area range (1,134 fail only this criterion): 1,362 above 5,400 Å² and the 2 at exactly 800.00 Å² (I-05).
- **Bias / limitation:** The cut removes interfaces that the literature treats as biological: EPPIC calls every interface above 2,200 Å² per side biological [Duarte2012], and homodimer interfaces reach about 5,000 Å² per side [Janin2007]. Of the 1,362 enriched pairs above 5,400 Å², 82.8 % are homo pairs, and the 1,194 that meet the length and P limits have a median P of 0.029; intertwined associations are common among homo-oligomers [MacKinnon2013].
- **Configuration:** `quality_filter.min_interface_area` (800), `quality_filter.max_interface_area` (5400)
- **Status:** `implemented`

### I-07 / I-08 · Final P-value of at most 0.4, and gold/silver tiers at P < 0.2
- **Decision:** Step 07 keeps pairs with P ≤ 0.4 and labels them *gold* if P < 0.2 and *silver* otherwise.
- **Reason:** Arbitrary specificity classes, defined to show a deep-learning model data of different specificity and to see whether that changed the results; the differences were small. P is not a calibrated probability that an interface is biological.
- **Precedent:** None for 0.4 or 0.2. The direction (lower P, more specific) follows PISA's definition [PISA-ipv].
- **Effect:** 16,776 rows with 0.4 < P ≤ 0.5 fail the cut (12,656 fail only this criterion). In train/validation, 29,989 of 39,747 pairs (75.4 %) are *gold*; in the test set, 449 of 750.
- **Bias / limitation:** P rewards hydrophobicity. Interfaces of complexes are about as apolar as the average protein surface, homodimer interfaces are more apolar [Janin2007; LoConte1999], and weak homodimers are more polar [Nooren2003b]. *Gold* means "more hydrophobic than random", not "biological with more certainty". After step 08, 67.0 % of homo and 61.5 % of hetero pairs are *gold*. Before the filter of step 08b, PRODIGY-cryst (I-02) called 97.6 % of the *gold* and 88.7 % of the *silver* final positives biological, so the tiers do order specificity.
- **Configuration:** `quality_filter.max_p_value` (0.4), `quality_filter.gold_max_p_value` (0.2)
- **Status:** `implemented`

### I-09 · No ΔiG cut-off; solvation-density filter deprecated
- **Decision:** No cut-off on PISA's ΔiG is applied. The optional solvation-density filter (ΔiG / interface area) has been removed from step 07 and from the configuration.
- **Reason:** The solvation-density filter was an idea that was tested and is now deprecated.
- **Precedent:** No published cut-off on ΔiG or on ΔiG per unit area for separating biological interfaces from crystal contacts was found. ΔiG covers solvation only [PISA-idg], PISA energies are not accurate enough to be the only classification criterion [Capitani2016], and the P-value already compares ΔiG with surface patches of the same area [PISA-ipv].
- **Effect:** None (filter removed). The Seq8 cut (density < −0.0068) would have removed 9,356 of the 109,825 pairs that pass step 08 (8.5 %).
- **Status:** `implemented`

### I-10 · Interfaces with missing or invalid area or P-value are dropped
- **Decision:** Step 02 drops interfaces without an area or a P-value, with a non-finite area, or with P outside [0, 1].
- **Reason:** Without these values the criteria of I-03 and I-04 cannot be applied; P is a probability in [0, 1] by definition [PISA-ipv].
- **Effect:** 0 interfaces.
- **Status:** `implemented`

### I-11 · Interfaces with a number of molecules other than two are dropped
- **Decision:** Step 02 drops any PISA interface that does not have exactly two molecules.
- **Reason:** In PISA an interface is always formed by two structures: its area is computed from the two isolated structures and the pair, and divided by two [PISA-iarea]; PDBe describes it as the interface between two chains [PISA-tut]. The filter is defensive.
- **Precedent:** AlphaFold 3 defines interfaces as pairs of chains [Abramson2024]; DIPS-Plus takes every pair of chains as a binary complex [Morehead2023]; PINDER is built from dimers [Kovtun2024].
- **Effect:** 0 interfaces.
- **Status:** `implemented`

### I-12 · The interface is the unit of the dataset
- **Decision:** Each PISA interface between two protein chains is one row. Copies related by non-crystallographic symmetry and equivalent assemblies of the same entry stay as separate rows until pruning in step 10.
- **Reason:** Kept as is and declared. Redundancy among copies is controlled by sequence in step 10 (R-01, R-05, R-07).
- **Precedent:** AlphaFold 3 groups interfaces by the pair of chain clusters and samples them with a weight that depends on cluster size [Abramson2024]; PINDER clusters interfaces and keeps one representative per cluster in its test set [Kovtun2024]. Xu and Dunbrack recommend grouping identical interfaces, such as each copy of the same homodimer interface in a tetrameric asymmetric unit [XuDunbrack2019]. In PPIRef, 322,454 interfaces reduce to 45,553 after near-duplicates are removed [Bushuiev2024].
- **Effect:** With the first assembly only (A-12), no entry contributes pairs from more than one assembly; previously 15,920 entries (22.9 % of candidate entries) did. The cap of R-01 reduces the 105,440 rows left after step 08b (I-02) to 43,710.
- **Bias / limitation:** Before pruning, the number of rows is not the number of distinct interactions, and copies of one interface within an entry are not independent observations.
- **Status:** `implemented`

### I-13 · Atomic contact: heavy atoms within 5 Å (AlphaFold 3 criterion)
- **Decision:** Two chains are in contact when a pair of their polymer heavy atoms is within 5 Å, the interface definition of AlphaFold 3. Every pair also reports its clashes (`n_clash_atoms`, atoms of either chain closer than 1.7 Å to the other; `max_clash_fraction`, the larger fraction of the two chains), and pairs where a chain has more than 30 % of its atoms in clashes are excluded, the chain-removal rule of AlphaFold 3.
- **Reason:** Alignment with the interface and clash criteria of AlphaFold 3 [Abramson2024]. The clash rule is a safeguard: the review of the full run found only local clashes (Effect).
- **Precedent:** AlphaFold 3: minimum heavy-atom separation below 5 Å [Abramson2024]. CAPRI/DockQ: any pair of heavy atoms within 5 Å [Basu2016]. 5 Å is the most common cut-off for interface residues [Esmaielbeiki2016] and the best value of a 3.5–12.5 Å sweep for separating biological from crystal interfaces [Elez2018].
- **Effect:** 0 pairs fail the contact check (I-14) or the clash rule. 584 of the 109,825 pairs that pass step 08 (327 entries) have atoms closer than 1.7 Å, 428 of them with a symmetry copy; the largest clashing fraction is 5.2 % (1jkm) and the median 0.1 %, so they are local clashes, not overlapping chains; 214 of these pairs are in the final positives. In the run of 2 November 2026, 118 of the 294 pairs closer than 1 Å came from two icosahedral capsids (1b35, 2iz9) and 81 involved chains with partial occupancy.
- **Bias / limitation:** Local clashes below the 30 % rule are kept and only reported; the EPPIC benchmark curation was stricter and discarded entries with more than 5 clashes below 1.5 Å [Duarte2012].
- **Configuration:** `contacts.cutoff` (5.0), `contacts.clash_distance` (1.7), `contacts.max_clash_fraction` (0.30; null reports only)
- **Status:** `implemented`

### I-14 · At least one residue pair in contact
- **Decision:** Step 08 keeps a pair if at least one residue pair is in contact (I-13).
- **Reason:** An existence test. Interface size is already guaranteed by PISA's area (I-03, I-06); a PISA interface without atomic contacts points to a chain-mapping problem.
- **Precedent:** AlphaFold 3 accepts a single heavy-atom pair below 5 Å [Abramson2024], and EPPIC one atom within 6 Å [Duarte2012]. Where contacts are the size criterion, more is required: at least 10 residues in contact in 3D Complex [Levy2006]; PINDER leaves systems with fewer than 7 interface residues out of clustering [Kovtun2024].
- **Effect:** 0 pairs excluded. The minimum observed is 11 residue pairs; the median, 104.
- **Configuration:** `contacts.min_residue_contacts` (1)
- **Status:** `implemented`

## Q · Structural quality

### Q-01 · R-free below 0.28
- **Decision:** Step 07 keeps pairs from entries with R-free < 0.28.
- **Reason:** The value is the user's choice. R-free is the cross-validated measure of agreement between model and data [Brunger1992; Read2011]; the default matches CrysFormer and is 0.01 below the CSAR cut.
- **Precedent:** CrysFormer: R-free ≤ 0.28 [Dun2023]. CSAR 2010: entries with R-free > 0.29 removed [Dunbar2011]. PLINDER test set: R-free ≤ 0.45 [Durairaj2024]. Large training sets (AlphaFold 2 and 3, PINDER) filter by resolution only [Jumper2021; Abramson2024; Kovtun2024].
- **Effect:** 1 row fails the criterion: its entry has no R-free (Q-02). All 66,806 X-ray entries with R-free have R-free below 0.28, which indicates that the RCSB query of the run (A-01) already applied the cut.
- **Bias / limitation:** An absolute cut ignores the dependence of R-free on resolution; the wwPDB ranks R-free against entries of similar resolution [Gore2017; wwPDB-VR].
- **Configuration:** `quality_filter.max_r_free` (0.28)
- **Status:** `implemented`

### Q-02 · Entries without R-free are excluded explicitly
- **Decision:** Entries without R-free are excluded, and step 07 records how many. Until now they were excluded silently, because a missing value does not pass the R-free filter.
- **Reason:** The R-free criterion (Q-01) cannot be evaluated without a value; these entries stay out.
- **Effect:** 1 pair of 1 entry excluded for lacking R-free.
- **Bias / limitation:** Entries deposited without R-free are lost. PISCES keeps them [Wang2005], and wwPDB validation reports give an R-free recalculated from the structure factors [Gore2017; wwPDB-VR].
- **Configuration:** `quality_filter.max_r_free` (0.28)
- **Status:** `implemented`

### Q-03 · Coverage above 0.86 in both chains
- **Decision:** Step 07 keeps a pair only if both chains have a coverage above 0.86 (metric in Q-04).
- **Reason:** The value is the user's choice.
- **Precedent:** Filtering by chain completeness is common, with varied values: Top2018 discards chains less than 60 % complete after residue filtering [Williams2022]; ProteinFlow tolerates up to 30 % missing at the ends and 10 % inside the chain [ProteinFlow]; the PINDER test set requires no atom gaps within 4 Å of the interface [PINDER-code]. None uses 0.86.
- **Effect:** 18,856 rows fail the criterion (15,424 fail only this one).
- **Bias / limitation:** With the default metric, disordered tags and tails count as much as internal gaps (Q-04).
- **Configuration:** `quality_filter.min_coverage` (0.86)
- **Status:** `implemented`

### Q-04 · Coverage metric, and a check of gaps at the interface
- **Decision:** By default, coverage is resolved residues divided by SEQRES residues; `internal_coverage` (resolved residues between the first and the last resolved one) remains available. Step 08 also reports, for each chain of a pair, the unresolved SEQRES residues within 2 positions of one of its interface residues (`molA_interface_gaps`, `molB_interface_gaps`); a maximum can be set to exclude pairs.
- **Reason:** The metric stays a user choice; completeness at the interface is reported directly and, by default, excludes nothing.
- **Precedent:** ProteinFlow measures internal gaps after disregarding missing ends [ProteinFlow]; PISCES offers lists without chain breaks [PISCES-web]; the PINDER test set requires no atom gaps within 4 Å of the interface [PINDER-code]; Hahn et al. recommend a global and a local quality filter [Hahn2021].
- **Effect:** 22,076 of the 224,091 chains of the enriched pairs (9.9 %) have coverage ≤ 0.86; 20,149 of them (91 %) have an internal coverage above 0.86 and would pass with `internal_coverage`. 44,663 of the 109,825 pairs that pass step 08 (40.7 %) have an unresolved residue within 2 positions of an interface residue, and so do 17,272 of the 40,497 final positives (42.7 %).
- **Bias / limitation:** The global metric penalises disordered expression tags and tails, and does not tell whether the missing residues lie at the interface.
- **Configuration:** `quality_filter.coverage_metric` (`coverage`); `contacts.interface_gap_window` (2); `contacts.max_interface_gaps` (null: report only)
- **Status:** `implemented`

### Q-05 · A resolved residue needs a heavy atom with occupancy above zero
- **Decision:** In step 06, a polymer residue counts as resolved if it has at least one heavy atom (ATOM or HETATM) with occupancy > 0 in the first model; steps 08 and 13 also ignore atoms with zero occupancy, so they neither count as residues nor create contacts. Until now occupancy was not checked.
- **Reason:** Correction to the standard: regions with zero occupancy should never be considered part of the experimental model [Wlodawer2008], and the PDBx dictionary groups zero-occupancy residues with unobserved ones (`pdbx_unobs_or_zero_occ_residues`) [PDBx].
- **Effect:** Compared with the previous behaviour on the same chains, 1,082 chains of 511 entries lose 5,539 resolved residues (at most 110 in one chain), and 60 of them fall to a coverage of 0.86 or less.
- **Status:** `implemented`

### Q-06 · Minimum chain length of 40 residues
- **Decision:** Step 07 keeps a pair only if both chains have more than 39 residues (at least 40), counted on the canonical sequence of the entity.
- **Reason:** 40 residues is the boundary between peptides and proteins used by PISCES and PINDER [PISCES-web; Kovtun2024]. The earlier README argument, that shorter chains cannot reach the MMseqs2 E-value of step 09, does not hold as written: with the step 09 parameters, identical fragments reach E ≤ 10⁻³ from about 18–20 residues [Test-MMseqs2].
- **Precedent:** PISCES lists: 40–10,000 residues [PISCES-web]. PINDER: chains under 40 residues are peptides, and the test set requires at least 40 [Kovtun2024; PINDER-code]. CrysFormer: at least 40 [Dun2023]. Top2018: at least 38 [Williams2022]. DIPS-Plus: more than 50 [Morehead2023].
- **Effect:** 2,509 rows have a chain of fewer than 40 residues (1,110 fail only this criterion, and 8 more fail only the two length limits, Q-07).
- **Bias / limitation:** Peptide–protein interfaces are excluded.
- **Configuration:** `quality_filter.min_chain_length` (39, i.e. length > 39)
- **Status:** `implemented`

### Q-07 · Maximum chain length: 1,000 residues by default
- **Decision:** Step 07 keeps a pair only if both chains have at most `max_chain_length` residues; the default changes from 600 to 1,000.
- **Reason:** The limit depends on the model to be trained and is the user's choice.
- **Precedent:** Bernett et al. restrict proteins to 50–1,000 residues because of GPU limits [Bernett2024]. PISCES and ProteinFlow allow up to 10,000 [PISCES-web; ProteinFlow].
- **Effect:** 910 rows have a chain longer than 1,000 residues (412 fail only this criterion, and 8 fail only the two length limits, Q-06). 3,849 of the 118,044 pairs that pass step 07 have a chain of 601 to 1,000 residues, which the previous default (600) excluded.
- **Bias / limitation:** Long chains are excluded rather than cropped; AlphaFold 2 crops them and samples chains according to their length [Jumper2021].
- **Configuration:** `quality_filter.max_chain_length` (1000)
- **Status:** `implemented`

### Q-08 · Non-standard residues: only MSE and SEC are accepted
- **Decision:** A chain is accepted only if all its polymer residues are standard amino acids, MSE (selenomethionine) or SEC (selenocysteine); any other polymer residue excludes the chain. SEM is removed from the allowed list: in the Chemical Component Dictionary it is O-benzyl-L-serine, not a selenium residue [CCD].
- **Reason:** Modified residues could be mapped to their parent amino acid, but the mapping sometimes causes problems; in this version they are excluded. MSE and SEC are kept because they are selenium residues, used for phasing or naturally present.
- **Effect:** 10,249 of the 176,931 chains scanned in step 08 carry a forbidden residue (319 residue types; the most frequent are LLP, KCX and CSO), and 8,219 pairs are excluded. MSE occurs in 8,950 chains and SEC in 28. With the previous list (MSE, SEM, SEC), SEM occurred in no chain, so removing it changes no pair.
- **Bias / limitation:** Stricter than common practice, which maps modified residues to the parent amino acid: AlphaFold 2 and 3 convert MSE to MET, and AlphaFold 2 maps other modified residues to their parent or to X [Jumper2021; AF2-code; Abramson2024]; the canonical PDB sequence uses the parent or X [PDBx]; CSAR changed modified residues back to their natural amino acids [Dunbar2011]. No dataset was found that excludes whole chains.
- **Configuration:** `special_residues.allowed` ([MSE, SEC])
- **Status:** `implemented`

### Q-09 · Terminal caps ACE and NH2 are ignored
- **Decision:** In steps 08 and 13, ACE and NH2 residues are neither counted nor treated as forbidden. In step 06, the terminal X that such caps leave in the canonical sequence is stripped, so capped and uncapped constructs of a protein share one sequence (F-06).
- **Reason:** The Chemical Component Dictionary classifies ACE (acetyl group) and NH2 (amino group) as non-polymer components without a parent amino acid: they are terminal caps, not amino acids [CCD].
- **Effect:** Not measured.
- **Configuration:** `special_residues.ignore` ([ACE, NH2])
- **Status:** `implemented`

### Q-10 · UNK and D-amino acids exclude the chain
- **Decision:** A chain with any UNK residue or any D-amino acid is excluded, under the rule of Q-08.
- **Reason:** Neither has a standard parent: UNK is an unknown residue, written as X in the canonical sequence [PDBx], and D-amino acids are D-peptide-linking components without a parent in the Chemical Component Dictionary [CCD].
- **Effect:** Included in the 10,249 chains of Q-08.
- **Bias / limitation:** Stricter than AlphaFold 3, which removes only chains made entirely of unknown residues [Abramson2024].
- **Configuration:** `special_residues.allowed` (UNK and D-amino acids not listed)
- **Status:** `implemented`

### Q-11 · A forbidden residue excludes the whole chain
- **Decision:** A forbidden residue (Q-08, Q-10) excludes the whole chain, wherever it lies, even far from the interface.
- **Reason:** A chain-level rule, kept in this version and declared as a limitation.
- **Effect:** 8,219 pairs excluded in step 08.
- **Bias / limitation:** Comparable resources judge quality locally: Hahn et al. recommend a global and a local filter [Hahn2021]; PLINDER applies its criteria to the ligand and the pocket [Durairaj2024]; PINDER checks gaps near the interface [PINDER-code]; CSAR checked that modified residues lay far from the binding site [Dunbar2011].
- **Status:** `implemented`

### Q-12 · Only the first model
- **Decision:** Steps 06, 08 and 13 read only the first model of each file.
- **Reason:** Standard practice; crystal structures generally have a single model [Biopython].
- **Precedent:** AlphaFold 2 processes only the first model [AF2-code]; DIPS-Plus picks the first model [Morehead2023].
- **Effect:** Not measured.
- **Status:** `implemented`

### Q-13 · Alternative conformations: highest-occupancy conformer in every step
- **Decision:** Where a residue has alternative locations, steps 08 and 13 keep the conformer with the highest mean occupancy (ties: the first in the file); step 06 counts a residue as resolved if any of its conformers has a heavy atom with occupancy > 0 (Q-05), so the three steps see the same residues. Until now steps 08 and 13 read only conformer A.
- **Reason:** Correction to the standard: AlphaFold 2 and 3 keep the alternative location with the largest occupancy [Jumper2021; Abramson2024], as Biopython does by default [Biopython]; conformer A is not always the most occupied.
- **Effect:** No chain fails to be read in step 08. With the previous behaviour, 88 chains had no atom in conformer A: they passed step 06 and led to the exclusion of 51 pairs in step 08.
- **Bias / limitation:** Only one conformer is kept; stricter datasets drop residues with alternative conformations (Top2018) or require none in the binding pocket (PLINDER) [Williams2022; Durairaj2024].
- **Status:** `implemented`

### Q-14 · Hydrogens are ignored
- **Decision:** Hydrogen and deuterium atoms are ignored when counting resolved residues (step 06) and in the residue and contact checks (steps 08 and 13).
- **Reason:** Standard practice: AlphaFold 3 removes hydrogens [Abramson2024]; CSAR ignores them because they are "usually a vestigial result of refinement" [Dunbar2011]; CrysFormer removes them [Dun2023].
- **Effect:** Not measured.
- **Configuration:** `enrichment.ignore_hydrogens` (true; step 06 only)
- **Status:** `implemented`

### Q-15 · Rows without entry metadata or chain sequence are dropped
- **Decision:** Step 06 drops rows without entry metadata (method, resolution, R-free) or without the sequence of either chain.
- **Reason:** The step 07 filters need the metadata, and the clustering and the leakage audit need the sequences.
- **Effect:** 0 rows: 155,899 pairs in and out.
- **Status:** `implemented`

## R · Redundancy

### R-01 · Cap the rows of each identical sequence pair
- **Decision:** Rows that share the same unordered pair of chain sequences, anywhere in the dataset, are capped at 2 rows for homo pairs and 4 for hetero pairs, keeping the best-ranked rows (R-02). The kept rows come from distinct PDB entries where possible: the best-ranked row of each entry goes first, so the kept copies are not NCS copies within one entry.
- **Reason:** Identical copies add no new information, and removing them is standard practice [Hobohm1992; Wang2005]. A few copies are kept because the same interface observed in more than one crystal form is stronger evidence of a biological contact [XuDunbrack2011]. The caps of 2 and 4 are pragmatic values and user parameters.
- **Precedent:** PISCES keeps one chain per group, or one entry in its by-entry mode [Wang2005]. PLINDER keeps *n* representatives per component, in its test set [Durairaj2024]. No published precedent was found for the values 2 and 4, or for halving the cap for homo pairs.
- **Effect:** 105,440 → 43,710 rows (59 % removed), counted after the filter of step 08b (I-02). The 105,440 rows contain 28,185 distinct sequence pairs; the most repeated pair has 1,170 rows. All 13,549 pairs that keep more than one row keep as many distinct entries as possible; in 6,194 of them the pair occurs in a single entry. Previous behaviour: 141,523 → 52,085 rows.
- **Bias / limitation:** When a sequence pair occurs in a single entry, its kept rows are still copies from that entry, which are not independent observations. The cap looks at sequences, not at interfaces: 1,839 of the 28,185 sequence pairs (6.5 %, 95 % of them homomeric) have, within one entry, interfaces whose areas differ by more than 20 % (for instance the different interfaces of a homotetramer), and for 1,279 of them (70 %) all the kept rows are of one of those interfaces, judged by area.
- **Configuration:** `pruning.exact_dedup_mode` (`softcap`; also `keep1`, `none`), `pruning.exact_K_homo` (2), `pruning.exact_K_heter` (4).
- **Status:** `implemented`

### R-02 · Ranking used to choose the kept rows
- **Decision:** Inside each redundancy group (R-01, R-05), rows are ranked by highest minimum coverage of the two chains, then lowest PISA P-value, largest interface area, lowest R-free, lowest resolution and most negative ΔiG. The top rows are kept.
- **Reason:** A design choice of the pipeline. Completeness comes first, the same principle PDBe uses to rank structures [PDBe-API]; the order then favours the most reliable (lowest P-value) and largest interface of the group, and crystallographic quality only breaks the remaining ties.
- **Precedent:** PDBe `best_structures`: UniProt coverage, then resolution [PDBe-API]. PISCES: highest resolution [Wang2005]. Dockground: best resolution, then fewest missing residues [Kotthoff2022].
- **Bias / limitation:** No published set uses the PISA P-value or the interface area to choose a representative; resolution and R-free act only as late tie-breakers. The ranking does not tell apart the different interfaces of one sequence pair, so the kept rows can all be of one of them (R-01).
- **Configuration:** None; the order is fixed in `rank_quality()` of step 10.
- **Status:** `implemented`

### R-03 · Pair-cluster pruning: optional, at 95 %
- **Decision:** The second redundancy stage, which prunes large pair clusters (R-04, R-05), is optional and off by default. When it is enabled, rows are grouped by pair cluster at 95 % identity (`pair_cluster_95`, H-07), so point mutants and variants with a few substitutions of the same complex fall in the same group.
- **Reason:** Pruning whole groups discards point mutants and conformational variants that some tasks need (R-06); large families are down-weighted instead (R-04, S-06). When pruning is wanted, the 95 % level follows the RCSB clusters.
- **Precedent:** RCSB publishes sequence clusters at 100, 95, 90, 70, 50 and 30 % [RCSB-clust]; at ≥ 95 % identity two proteins have practically the same structure, Δ = 0.40·e^(1.87H) Å, about 0.44 Å at 95 % identity [Chothia1986].
- **Effect:** Pruning off: the 43,710 rows left by R-01 are kept; they form 17,216 pair clusters at 95 %. Previous behaviour (pruning on, at 98 %): stage 2 reduced 52,085 rows to 43,164.
- **Configuration:** `pruning.cluster_pruning` (false), `pruning.pair_key` (`pair_cluster_95`).
- **Status:** `implemented`

### R-04 · Large families are down-weighted, not pruned
- **Decision:** By default no pair cluster is pruned: large families keep all their distinct pairs and are down-weighted by the family weight of S-06. If pair-cluster pruning is enabled (R-03), clusters with more than `activate_n` rows are pruned as in R-05.
- **Reason:** A weight keeps the data and corrects the over-representation, as comparable training sets do.
- **Precedent:** AF2 accepts chains with probability inverse to the size of their cluster [Jumper2021], and AF3 weights by 1/N_clust [Abramson2024]. No published activation threshold for pruning was found.
- **Effect:** No pair cluster is pruned. Of the 17,216 pair clusters at 95 %, 1,966 have more than 8 rows after step 08b (71,045 rows, 67.4 %); after the exact-pair cap of R-01, which runs first, 560 do (11,051 of the 43,710 rows, 25.3 %), and those are the ones pruning would act on. Previous behaviour (pruning on): stage 2 reduced 52,085 rows to 43,164.
- **Bias / limitation:** The weight is a column of the output: a model that ignores it sees large families at full weight.
- **Configuration:** `pruning.cluster_pruning` (false), `pruning.activate_n` (8, only with pruning on); `splits.weight_exponent` (0.5).
- **Status:** `implemented`

### R-05 · Keep interface-area diversity in large pair clusters
- **Decision:** When pair-cluster pruning is enabled (R-03), in each pair cluster selected by R-04, if the interface areas are almost equal (range ≤ 4 % of the mean), the 2 (homo) or 4 (hetero) best-ranked rows are kept; otherwise 4 (homo) or 8 (hetero) rows are kept, spread over up to four interface-area quantile bins (the best row of each bin, then the next best). The cluster type follows the majority of its rows; mixed clusters use the hetero values.
- **Reason:** The values are the user's choice. The rule keeps geometric diversity, approximated by the interface area, inside heavily repeated complexes.
- **Precedent:** Keeping the distinct binding modes of a pair of homologues is established practice: 3did groups interfaces into "interaction topologies" by interface similarity [Mosca2014], and ProtCID clusters the interfaces common to several crystal forms [XuDunbrack2011]. Homologues usually, but not always, interact in the same way [Aloy2003].
- **Effect:** Not applied (pruning off). Previous behaviour (pruning on, at 98 %): 52,085 → 43,164 rows; all 19,215 pair clusters at 98 % remained.
- **Bias / limitation:** Area is a proxy for the binding mode: two different binding modes can have the same area. A measure of interface similarity (e.g. iDist [Bushuiev2024]) would group binding modes directly. The 4 % range and the caps 2/4/4/8 have no published precedent.
- **Configuration:** `pruning.area_range_pct_thr` (4.0), `pruning.K_redund_homo` (2), `pruning.K_redund_heter` (4), `pruning.K_var_homo` (4), `pruning.K_var_heter` (8).
- **Status:** `implemented` (applies only when `pruning.cluster_pruning` is enabled)

### R-06 · Point mutants are kept
- **Decision:** Point mutants and variants with a few substitutions are distinct sequence pairs, so the exact-pair cap (R-01) keeps them, and with pair-cluster pruning off (R-03) they all stay. Enabling the pruning keeps only a few per complex.
- **Reason:** For tasks such as the effect of mutations on binding (ΔΔG), mutants are the data [Jankauskaite2019]. They must not cross a split: variants of one complex on both sides of a split cause leakage, as in SKEMPI splits by PDB entry [Bushuiev2024]; the splits by pair cluster keep them on one side (S-01).
- **Precedent:** SKEMPI 2.0 keeps all 7,085 mutations [Jankauskaite2019]; PPIRef replaces the SKEMPI splits by PDB entry with leak-free splits [Bushuiev2024].
- **Bias / limitation:** Mutant series enlarge the families of well-studied complexes; the family weight (S-06) limits their influence. At ≥ 98 % identity the structure is practically the same [Chothia1986], so for a purely structural task they add little.
- **Configuration:** `pruning.cluster_pruning` (false).
- **Status:** `implemented`

### R-07 · NCS copies and equivalent assemblies are not collapsed per entry
- **Decision:** NCS copies and equivalent assemblies of one entry are not collapsed per entry; they are reduced only through the sequence-based cap of R-01 (and R-05 when pair-cluster pruning is enabled). The unit of the dataset is the interface (I-12).
- **Reason:** NCS copies share the pair of sequences, so the cap of R-01 already limits them.
- **Precedent:** The literature treats the entry or the crystal form as the unit of independent observation: PISCES offers by-entry culling and culling within each entry [Wang2005]; ProtCID counts crystal forms, not copies [XuDunbrack2011]; AF3 expands only the first bioassembly [Abramson2024].
- **Effect:** The 43,710 rows after step 10 come from 34,942 entries, and the 40,497 final positives from 32,989. Previous behaviour: the 43,164 rows after pruning came from 29,391 entries.
- **Bias / limitation:** The number of rows is not the number of independent observations. The changes of R-01 (distinct entries) and A-12 (first assembly by default) reduce this dependence but do not remove it: when a sequence pair occurs in a single entry, its kept rows are still copies from that entry.
- **Configuration:** None.
- **Status:** `implemented`

### R-08 · Cap identical copies, keep the rest and weight families
- **Decision:** By default the pipeline caps identical sequence pairs (R-01) and keeps every distinct pair, with a family weight per row (S-06) for weighting or sampling at training time. Pair-cluster pruning (R-03 to R-05) can be enabled for a smaller, less redundant set, and `exact_dedup_mode: none` keeps even the identical copies.
- **Reason:** The mode is the user's choice. The default keeps the data and leaves the rebalancing to training, as comparable training sets do; the output tables also carry the chain and pair cluster IDs at every threshold for per-cluster sampling.
- **Precedent:** Comparable training sets keep all members and sample per cluster: AF2 (acceptance ∝ 1/cluster size at 40 %) [Jumper2021], AF3 (weight ∝ 1/N_clust) [Abramson2024], ProteinMPNN (one random member of each 30 % cluster per draw) [Dauparas2022; ProteinMPNN-code], and PINDER, which credits interface cluster-based sampling of its training split with strong generalization improvements [Kovtun2024]. Non-redundant reference sets remove redundancy [Hobohm1992; Wang2005; Bushuiev2024], and a strictly non-redundant training set slightly improved the prediction of difficult cases [Hamp2015].
- **Effect:** 105,440 rows after step 08b → 43,710 (−59 %); 48,311 → 34,942 entries; homomeric share 78.4 % → 75.8 %. Previous behaviour (pruning on): 141,523 → 43,164 rows (−70 %).
- **Bias / limitation:** Without the weight, large families dominate training and metrics; pruning instead discards data [Teufel2023], including the conformational variability of repeated complexes.
- **Configuration:** `pruning.exact_dedup_mode` (`softcap`), `pruning.cluster_pruning` (false), `splits.weight_exponent` (0.5).
- **Status:** `implemented`

## H · Homology and clustering

### H-01 · Homology by sequence, not by structure
- **Decision:** Homology between chains is measured by sequence identity with MMseqs2 [Steinegger2017], for the clustering (step 09) and for the leakage audit (step 14), with the same linking rule. No structure- or domain-based clustering is used.
- **Reason:** Sequence identity is the standard homology criterion of comparable training sets [Jumper2021; Abramson2024; Dauparas2022] and of the PDB's own clusters [RCSB-clust].
- **Precedent:** Structure-aware alternatives: PINDER combines Foldseek and MMseqs2 clusters [Kovtun2024]; PPIRef removes near-duplicate interfaces with iDist [Bushuiev2024]; Dockground removes redundancy with MM-align [Kotthoff2022].
- **Effect:** 0 leaks by sequence at the design threshold in the 11 audited sets. Structural leakage is not measured.
- **Bias / limitation:** The guarantee holds at sequence level only, and interfaces can be similar when sequences are not: about 90 % of native interfaces have a close structural neighbour even without similarity between the monomers [Gao2010], and a sequence-based split of DIPS leaves a near-duplicate in training for 53 % of the test [Bushuiev2024]. "0 leaks" therefore means "0 leaks by sequence at the design threshold". Structure-based clustering (Foldseek [vanKempen2024]), interface similarity (iDist [Bushuiev2024]) or domain classification (ECOD [Cheng2014]) can be used instead or in addition.
- **Configuration:** None.
- **Status:** `implemented`

### H-02 · Identity thresholds
- **Decision:** Chains are clustered at the identity levels of `clustering.thresholds`, by default 95, 90, 60, 40 and 30 %. Step 09 checks that the list holds every level the configuration reads: the novelty level of the test set, the split schemes, the hard-negative levels and the pruning key if pruning is enabled.
- **Reason:** The levels are the user's choice; the default follows the levels of the precedents, plus 60 % for one of the validation schemes (S-02).
- **Precedent:** RCSB: 100, 95, 90, 70, 50 and 30 % [RCSB-clust]; AF2 and AF3: 40 % [Jumper2021; Abramson2024]; ProteinMPNN and PINDER: 30 % [Dauparas2022; Kovtun2024]. The twilight zone of sequence alignments lies at 20–35 % identity [Rost1999].
- **Effect:** From 31,722 unique sequences: 18,734 clusters at 95 %, 16,760 at 90 %, 12,847 at 60 %, 9,354 at 40 % and 6,821 at 30 %. Previous behaviour: eight levels, from 98 to 40 %.
- **Configuration:** `clustering.thresholds` ([0.95, 0.90, 0.60, 0.40, 0.30]).
- **Status:** `implemented`

### H-03 · Linking criterion: alignment identity
- **Decision:** Two sequences are linked at a threshold when the MMseqs2 identity `fident` (identical positions / aligned columns) is at least the threshold and the coverage rule (H-04) holds.
- **Reason:** `easy-search` computes the real alignment identity, not an estimate from the score [MMseqs-src; MMseqs-wiki]. Identity over the aligned columns of a local alignment can mislead without a coverage control [Teufel2023]; the pipeline applies one (H-04).
- **Precedent:** GraphPart divides by the length of the shorter sequence [Teufel2023]; AF3 uses the percentage of residues of the evaluated chain that are identical to the training chain [Abramson2024].
- **Bias / limitation:** Alignments of well-conserved fragments overestimate identity [Wang2003]; the coverage rule limits this effect but does not remove it.
- **Configuration:** `clustering.thresholds`.
- **Status:** `implemented`

### H-04 · Coverage rule: 90 % of either chain
- **Decision:** A link also requires the alignment to cover ≥ 90 % of either chain (query or target); every pair is searched in both directions. The clustering (step 09) and the leakage audit (step 14) share this rule. Step 14 repeats the audit without the coverage rule and reports it (it does not enforce it).
- **Reason:** Coverage of either chain links a fragment or an isolated domain to its full-length protein, the use MMseqs2 documents for its coverage mode 1 [MMseqs-wiki]. Either direction is accepted because MMseqs2 may keep only a shorter local alignment in one direction. The value is the user's choice: how strictly training and evaluation data must be separated depends on what a model has to generalise to, and a lower coverage trades construct leaks for larger, chained clusters (Effect).
- **Precedent:** RCSB requires 80 % coverage of the member only [RCSB-clust]; the Linclust benchmark uses a minimum coverage of 90 % [Steinegger2018]. For the audit, AF3 measures identity over the evaluated chain [Abramson2024], and GraphPart, by default, over the length of the Needleman–Wunsch alignment, gaps included (over the shorter sequence in its MMseqs2 mode), without a coverage threshold [Teufel2023].
- **Effect:** Without the coverage rule, the audit would count as leaks at the design threshold 61 of the 435 C3 chains, 35 of the 292 unseen C2 chains, 7 of the 300 C2 pairs, 16 of the 150 C1 pairs, 455 of the 3,976 validation pairs of `split_c60`, 717 of the 4,393 of `split_c40` and 4 of the 289 test negatives; against the negatives, 1 C3 chain and 2 unseen C2 chains. These are links between chains that share a region but cover each other by less than 90 % (Bias). Simulations of steps 09–14 with a lower coverage, run before selenomethionine was converted in step 08b (I-02) and measured with the control of U-02 (by GraphPart's identity / by AlphaFold 3's with a significant alignment): at 0.90, 6 / 22 of the 431 C3 chains, 3 / 21 of the 286 unseen C2 chains and 9 / 84 of the 3,989 `split_c40` validation pairs; at 0.80, 2 / 9 of 429, 2 / 9 of 289 and 0 / 225 of 4,007, and step 14 flags two validation rows of one homomeric sequence pair at exactly 40 % identity that the search of step 09 had not linked; at 0.50, 1 / 4 of 430, 0 / 9 of 289 and 0 / 31 of 4,009, but the largest chain cluster at 40 % grows from 5,065 to 8,592 sequences, the largest training family from 3,670 to 9,300 rows, and the test negatives fall from 342 to 206. In the reference run, with the conversion, the control finds at 0.90 3 / 20 of 435, 3 / 13 of 292 and 43 / 95 of 4,393 (U-02); its `split_c40` validation set is a different draw (S-04).
- **Bias / limitation:** For leakage-free splits the rule is permissive: two multidomain chains that share an identical domain but cover each other by less than 90 % (e.g. two ~300-residue chains sharing a ~150-residue domain) are not linked at any threshold. Such a chain can count as unseen in the test set (C2, C3) while sharing a domain with training; by the AF3 definition it would count as high homology. The rule also leaves unlinked two constructs of the same protein when each carries residues the other lacks (expression tags, linkers, termini) that exceed 10 % of its length; a construct that only adds residues to the other is linked (U-02).
- **Configuration:** `clustering.coverage` (0.90), `clustering.cov_mode` (1; 0 requires both chains).
- **Status:** `implemented`

### H-05 · MMseqs2 search parameters
- **Decision:** One all-vs-all MMseqs2 `easy-search` with E-value ≤ 10⁻³, sensitivity `-s 7.5`, `--max-seqs 20000`, minimum identity equal to the lowest threshold (30 %) and no coverage filter at search time (`-c 0`; coverage is applied afterwards, H-04).
- **Reason:** 10⁻³ is the MMseqs2 default E-value and 7.5 the "sensitive" end of its sensitivity scale [MMseqs-src]. `--max-seqs` caps the prefilter candidates per query and "affects sensitivity" [MMseqs-src; MMseqs-wiki]; it must exceed the largest family, or links are lost and clusters leak.
- **Precedent:** MMseqs2 defaults: `-e 0.001`, `--max-seqs 300` [MMseqs-src]. PLINDER: E < 0.01, max seqs 5,000 [Durairaj2024]. GraphPart aligns all against all without prefilter [Teufel2023].
- **Effect:** With `max_seqs` 20,000, the 11 audited sets have 0 leaks by sequence; the largest cluster has 5,891 sequences (at 30 %; 5,066 at 40 %). In the run of 2 November 2026, steps 09–14 with `max_seqs: 1000` left 15 of 4,195 validation pairs (90 % scheme) and 13 of 44 test hard negatives with a homologue in their reference set at the design threshold.
- **Bias / limitation:** The search still depends on the k-mer prefilter; an alignment without prefilter (`--prefilter-mode 2`) would remove the dependence on `max_seqs` [MMseqs-src; Teufel2023]. Homologues below the lowest threshold (30 % identity) are not searched.
- **Configuration:** `clustering.evalue` (1.0e-3), `clustering.search_sensitivity` (7.5), `clustering.max_seqs` (20000); the search floor is the lowest value of `clustering.thresholds`.
- **Status:** `implemented`

### H-06 · Connected components (single linkage)
- **Decision:** At each threshold, the clusters are the connected components (single linkage) of the graph of linked sequences.
- **Reason:** Two sequences in different components are never linked above the threshold, so cluster-based splits do not leak by sequence, provided the search reports every hit (H-05). Greedy methods (set cover, CD-HIT, Linclust) only guarantee that each member is similar to its representative and "should not be used for homology partitioning" [Teufel2023; Steinegger2018].
- **Precedent:** PLINDER splits train and validation by connected components [Durairaj2024]; PINDER masks all cluster members and transitive hits [Kovtun2024]; GraphPart and the Blue and Cobalt algorithms also separate by graph [Teufel2023; Petti2022].
- **Effect:** 0 leaks by sequence in the 11 audited sets. Largest component: 1,258 sequences (4.0 %) at 90 %, 4,552 (14.3 %) at 60 %, 5,066 (16.0 %) at 40 % and 5,891 (18.6 %) at 30 %. The previous version (set-cover clusters) left 115 of 3,664 validation pairs with a homologue in training at 90 %.
- **Bias / limitation:** Single linkage can merge families into giant components through a few bridging sequences [Petti2022]; in DIPS one component of near-duplicates holds 36 % of the interfaces [Bushuiev2024]. Two members of a cluster can be less similar than the threshold.
- **Configuration:** `clustering.method` (`connected_components`). The alternative `set_cover` is in fact a greedy MMseqs2 clustering: with the default engine (`mmseqs cluster`) and `--cov-mode 1`, MMseqs2 clusters greedily by representative (CD-HIT-like), not by set cover [MMseqs-wiki]; the `linclust` engine is also representative-based [Steinegger2018]. The name is inherited from the earlier version.
- **Status:** `implemented`

### H-07 · Pair cluster: unordered pair of chain clusters
- **Decision:** At each threshold, the pair cluster of an interface is the unordered pair of the clusters of its two chains (`min||max`). Pair clusters are the groups for the optional pruning (R-03), the C1 test class, the train/validation splits and the caps on hard negatives.
- **Reason:** It is the definition AF3 uses for interface clustering in training and evaluation ("a join on the cluster IDs of the constituent chains") [Abramson2024]. The evaluation of pair-input predictions depends on whether each of the two proteins was seen in training [Park2012].
- **Precedent:** PINDER labels each system with paired cluster IDs [Kovtun2024]; Bernett et al. assign pairs to the blocks of a proteome partition by their two proteins [Bernett2024]; DataSAIL formalises such two-dimensional splits [Joeres2025].
- **Effect:** The 105,440 rows after step 08b form 17,216 pair clusters at 95 %; the 39,747 train/validation rows form 11,067 pair clusters at 60 % and 8,066 at 40 %.
- **Bias / limitation:** A split by pair cluster does not separate single proteins: a chain can appear in training and validation with different partners (S-05).
- **Configuration:** None (columns `pair_cluster_cXX`).
- **Status:** `implemented`

### H-08 · Cluster only the chains that pass steps 07–08
- **Decision:** Only the chains of the pairs that pass the quality and structure filters (steps 07–08) are clustered; clusters are not computed over the whole PDB.
- **Reason:** Partitioning the dataset at hand is standard practice [Teufel2023; Joeres2025]; AF3 also measures homology against its own filtered training set [Abramson2024].
- **Precedent:** AF2 and ProteinMPNN instead use clusters of the whole PDB computed before filtering [Jumper2021; ProteinMPNN-code].
- **Effect:** 31,722 unique sequences (210,880 chain instances in 105,440 rows) are clustered.
- **Bias / limitation:** Cluster IDs change when the filters change, so they are not stable between runs or comparable with other datasets, and they say nothing about homology with external data, such as the pre-training set of a model.
- **Configuration:** None.
- **Status:** `implemented`

## T · Test set

### T-01 · Three test classes by what the model has seen
- **Decision:** The external test set has three classes. A protein is *seen* if its chain cluster at 40 % identity has a row in the train/validation pool (every row not in the test set). C3: neither protein seen; C2: exactly one seen; C1: both seen, but no pair of the pool joins homologues of both chains. Metrics are reported per class.
- **Reason:** Performance on pairs depends on whether their proteins were seen in training, so it has to be reported separately for each class.
- **Precedent:** Park and Marcotte: "predictive performances should be reported separately for each distinct test class", with C1 (both proteins in training), C2 (one) and C3 (none) [Park2012]. Hamp and Rost apply the classes by sequence similarity [Hamp2015], and with similarity minimised, sequence-based predictors fall to random performance [Bernett2024].
- **Effect:** 150 C1, 300 C2 and 300 C3 pairs, with 0 leaks by sequence in each class at 40 % (U-01). Previous behaviour: a train-like bucket, new pairs at 70 and 80 % and unseen proteins at 40 and 50 %; 37 % (70 %) and 48 % (80 %) of the new-pair test had a training pair homologous at 40 % in both chains.
- **Bias / limitation:** Validation is drawn from the same pool (S-03), so in a given scheme a *seen* protein may be only in validation: in `split_c40`, 139 of the 150 C1 pairs have both proteins in training, 5 one and 6 neither, and 26 of the 300 C2 pairs have their seen protein only in validation (`split_c60`: 124, 6 and 20; 30).
- **Status:** `implemented`

### T-02 · Order: held-out proteins, then new pairs, then validation
- **Decision:** Step 11 first holds out chain clusters (C3 and C2), then removes whole pair clusters for C1; step 12 splits the rest into training and validation, and step 13 sends every negative with a held-out protein to the test set.
- **Reason:** Holding out a protein removes every row that contains it, so it has to come first; otherwise C1 or the validation set could have used those proteins.
- **Precedent:** PINDER and PLINDER build the test set first and then clean the training set [Kovtun2024; Durairaj2024].
- **Status:** `implemented`

### T-03 · Novelty at 40 % identity
- **Decision:** A protein or a pair is new when no training chain or pair is linked to it at 40 % identity (connected components of step 09, H-06).
- **Reason:** 40 % is the novelty level of the reference models; the level is a parameter.
- **Precedent:** AlphaFold 3 calls chains and interfaces *low homology* below 40 % identity to training [Abramson2024]; AlphaFold 2 removes test chains with a training template above 40 % [Jumper2021]; the gold standard of Bernett et al. separates proteins at 40 % [Bernett2024]. PINDER and ProteinMPNN use 30 % [Kovtun2024; Dauparas2022].
- **Bias / limitation:** Novelty is by sequence only (H-01), with the coverage rule of H-04.
- **Configuration:** `test_set.level` (0.40); must be one of `clustering.thresholds`.
- **Status:** `implemented`

### T-04 · Clusters drawn at random, with a size cap
- **Decision:** Chain clusters (C2, C3) and pair clusters (C1) are drawn at random with a fixed seed, among those with at most 50 rows. C3 homo pairs need one held-out cluster, C3 hetero pairs both, and C2 pairs one, chosen so that the partner stays seen.
- **Reason:** Random draws keep the test representative of the families of the dataset (T-10), without the bias towards rare families of the previous version; the cap keeps large families in training, the aim of partitioning methods that lose as little data as possible.
- **Precedent:** GraphPart and the algorithms of Petti and Eddy separate homologues while "keeping as many sequences as possible" [Teufel2023; Petti2022].
- **Effect:** 659 chain clusters at 40 % are held out (C3 and C2) and 150 pair clusters are taken for C1; the test set and the rows it removes from training total 3,963 of the 43,710 rows (9.1 %), and the train/validation pool keeps 39,747. A simulation on the run of 2 November 2026 had estimated about 7.5 %.
- **Bias / limitation:** Families larger than the cap never enter the test set.
- **Configuration:** `test_set.max_cluster_rows` (50), `test_set.seed` (7).
- **Status:** `implemented`

### T-05 · C1: new pairs of seen proteins at 40 %
- **Decision:** A C1 pair is taken with its whole pair cluster at 40 %, only if both of its chain clusters keep other partners in training.
- **Reason:** It is the pair-level novelty that AlphaFold 3 requires for a low-homology interface, the same level as the conflict filter of the hard negatives (N-07).
- **Precedent:** AlphaFold 3 leaves an interface out of its low-homology set if both chains have more than 40 % identity to two chains of the same training complex [Abramson2024]. Close homologues that interact "almost invariably interact the same way" [Aloy2003].
- **Effect:** 150 C1 pairs, 96 homomeric and 54 heteromeric: heteromeric pair clusters whose two proteins keep other partners in the train/validation pool run out as the class is filled, and homo pairs fill the class (T-09). Before the filter of step 08b (I-02), 75 homo and 75 hetero.
- **Status:** `implemented`

### T-06 · One representative per pair cluster
- **Decision:** Each test pair represents a different pair cluster at 40 %.
- **Reason:** Copies and close homologues of one complex would weigh the metrics towards it.
- **Precedent:** PINDER limits its test set to one system per cluster, and AlphaFold 3 samples one interface per cluster [Kovtun2024; Abramson2024].
- **Status:** `implemented`

### T-07 · Representative chosen at random
- **Decision:** The representative of each pair cluster is a random row of the required kind (homo or hetero).
- **Reason:** Choosing by quality would mix structure quality with novelty; a random row has the quality of the training data.
- **Precedent:** AlphaFold 3 samples one interface per cluster without a quality criterion [Abramson2024]; a test set that does not represent the population of interest does not generalise [Park2012].
- **Status:** `implemented`

### T-08 · No train-like bucket
- **Decision:** The test set has no bucket drawn from the most populated families.
- **Reason:** A train-like test measures what a random split measures, which does not generalise to new proteins; validation (S-02) already covers pairs close to training.
- **Precedent:** Park and Marcotte show that random cross-validation estimates "should not be expected to generalize" [Park2012]; with random splits, models learn sequence similarity and node degree [Bernett2024].
- **Effect:** Previous behaviour: 259 of the 300 train-like pairs had a homologous pair at ≥ 98 % in training.
- **Status:** `implemented`

### T-09 · Half homo, half hetero when possible
- **Decision:** Each class aims at half homo and half hetero pairs; when one kind runs out, the other fills the class. C2 is hetero by definition. Metrics are reported overall and separately for homo and hetero pairs (A-10).
- **Reason:** About four fifths of the train/validation pairs are homo (A-10); balancing the classes keeps hetero pairs visible in the test metrics.
- **Precedent:** AlphaFold-Multimer reports homomers and heteromers separately [Evans2021]; homo- and hetero-oligomeric interfaces differ [JonesThornton1996]; PINDER biases its test set towards heterodimers [Kovtun2024].
- **Effect:** C3: 150 homomeric and 150 heteromeric pairs; C1: 96 homomeric and 54 heteromeric (T-05); C2: 300 heteromeric. 246 of the 750 test pairs are homomeric.
- **Configuration:** `test_set.homo_frac` (0.5).
- **Status:** `implemented`

### T-10 · No stratification by length or interface area
- **Decision:** Test pairs are not stratified by chain length or interface area.
- **Reason:** No PPI benchmark was found that stratifies its test set this way; random draws follow the distribution of the families of the dataset, one pair per pair cluster (T-06).
- **Effect:** With one pair per family, 67 % of them heteromeric (T-09), the test has smaller interfaces and more *silver* pairs than the train/validation rows: a median of 1,345 against 1,619 Å² per side, 18.7 % against 10.8 % below 1,000 Å², and 40.1 % against 24.6 % *silver*. One random train/validation pair per family at 40 % narrows the gap without closing it (mean of 50 draws): its heteromeric pairs have a median of 1,433 Å² (test 1,313 Å²; all train/validation heteromeric rows 1,707 Å²) and 41.1 % *silver* (test 43.8 %; all heteromeric rows 22.9 %). Previous behaviour: strata proportional to the candidates, which made the unseen-protein interfaces smaller on average.
- **Bias / limitation:** The test reflects the families of the dataset, not its rows.
- **Status:** `implemented`

### T-11 · Class sizes: 300, 300 and 150
- **Decision:** 300 C3 pairs (150 homo, 150 hetero), 300 C2 pairs and 150 C1 pairs.
- **Reason:** The values are the user's choice; C1 is the scarcest class (T-05).
- **Precedent:** Comparable structural test sets range from 180 to about 2,000 pairs (PINDER-AF2 180, PINDER-S 250, PINDER-XL 1,955), and AlphaFold 3 validates on 600 protein–protein interfaces [Kovtun2024; Abramson2024].
- **Bias / limitation:** With 150 pairs, the 95 % interval of a proportion near 0.5 is about ±8 percentage points (binomial).
- **Configuration:** `test_set.n_c3` (300), `test_set.n_c2` (300), `test_set.n_c1` (150).
- **Status:** `implemented`

### T-12 · Class names follow Park and Marcotte
- **Decision:** C1, C2 and C3 mean what they mean in Park and Marcotte (both, one or no protein shared with training), defined by homology as in Hamp and Rost.
- **Reason:** The previous version used the same names with other meanings.
- **Precedent:** [Park2012; Hamp2015; Bernett2024]. DataSAIL formalises such two-dimensional splits [Joeres2025].
- **Status:** `implemented`

## S · Train/validation splits and weights

### S-01 · Split by pair cluster, not by row or chain
- **Decision:** Step 12 assigns whole pair clusters to training or validation. A pair cluster is the unordered pair of the clusters of its two chains (H-07), so all rows of the same pair cluster go to the same side.
- **Reason:** Homologous pairs stay on the same side: no training pair joins homologues of both chains of a validation pair at the scheme's identity.
- **Precedent:** AlphaFold 3 defines an interface cluster as the pair of cluster IDs of its chains (at 40 %) and uses it to sample and to evaluate [Abramson2024]; PINDER gives each PPI the pair of cluster IDs of its chains [Kovtun2024]; GraphPart shows that homologues must be kept in the same partition, and that methods that only control the cluster representative should not be used for partitioning [Teufel2023].
- **Effect:** 11,067 and 8,066 pair clusters at 60 and 40 % among the 39,747 train/validation rows; 0 leaks at pair level in both validation sets (U-01).
- **Bias / limitation:** Splitting by pair does not keep proteins apart: a protein can be in training and in validation with different partners (S-05).
- **Status:** `implemented`

### S-02 · Validation schemes at 60 and 40 % pair-cluster identity
- **Decision:** Step 12 builds one independent train/validation split by pair cluster per level of `splits.schemes`, by default 60 and 40 % (`split_c60`, `split_c40`); the user can choose other levels. Each scheme is a separate split, so one model is trained per scheme. `split_c40` is exported as `validation.csv`.
- **Reason:** A validation set split at 90 % stays close to training: most of its pairs had a homologous pair in training at 40 % (see Effect). Lower defaults bring validation closer to the low-homology validation used for model selection in AlphaFold 3. The thresholds are the user's choice.
- **Precedent:** Which method appears best "crucially depends on the sequence similarity between the test and the training set" [Hamp2015]; PLINDER builds its similarity graphs at 50, 70, 95 and 100 % [Durairaj2024]; AlphaFold 3 selects models on low-homology interfaces (no training complex with two chains above 40 % identity to the interface chains) [Abramson2024]; PINDER credits part of its gains to "the diverse and less leaky validation split" [Kovtun2024].
- **Effect:** In `split_c60`, 1,933 of the 3,976 validation pairs (48.6 %) have a training pair homologous at 40 % in both chains; in `split_c40`, none by construction. Previous behaviour (schemes at 90, 80 and 60 %): 2,669 of 4,185 (64 %) in `split_c90`, 2,092 of 4,184 (50 %) in `split_c80` and 1,662 of 4,188 (40 %) in `split_c60`.
- **Bias / limitation:** The separation is by sequence only, with the coverage rule of step 09: two multidomain proteins that share a single domain are not linked (U-02).
- **Configuration:** `splits.schemes` ([0.60, 0.40]); each level must be in `clustering.thresholds` (step 09 checks it). `splits.validation_split` (`split_c40`) selects the scheme exported as `validation.csv` and the one the family weight is computed on (S-06).
- **Status:** `implemented`

### S-03 · Send 10 % of the rows to validation
- **Decision:** Each scheme sends about 10 % of the train/validation rows to validation. The target is counted in rows, not in clusters.
- **Reason:** A 90/10 split of clusters is common practice (see Precedent). The value is the user's choice.
- **Precedent:** PLINDER uses a "random 90/10 split of component clusters" [Durairaj2024]; GraphPart's example keeps 10 % for testing and 5 % for validation [Teufel2023]; ProteinMPNN keeps 1,464 of 25,361 clusters (5.8 %) for validation [Dauparas2022]; the gold standard of Bernett et al. puts 59,260 of 274,500 PPIs (21.6 %) in validation [Bernett2024].
- **Effect:** 3,976 validation rows in `split_c60` (10.0 % of 39,747) and 4,393 in `split_c40` (11.1 %). Clusters stay whole, so the target is approached cluster by cluster (S-04): in `split_c40`, one family of 852 rows (influenza hemagglutinin, whose HA1 and HA2 chains share one chain cluster at 40 %) brings the heteromeric stratum to 1,267 validation rows against a target of 851.
- **Configuration:** `splits.val_frac` (0.10).
- **Status:** `implemented`

### S-04 · Stratify homo/hetero per cluster and select clusters greedily
- **Decision:** Each pair cluster gets a stratum by majority (homo if at least half of its rows are homomeric). Within each stratum, clusters are shuffled with a fixed seed and taken greedily until the stratum's validation target, in rows, is reached.
- **Reason:** Validation keeps the homo/hetero proportion of the pool while clusters stay whole.
- **Precedent:** GraphPart assigns clusters to partitions "such that the number of sequences of any label is approximately the same across partitions" [Teufel2023]; DataSAIL offers class-stratified splits without leakage [Joeres2025]. The majority rule and the greedy selection are implementation details.
- **Effect:** 79.6 % of the 39,747 train/validation rows are homomeric, and 78.6 % (`split_c60`) and 77.7 % (`split_c40`) of the validation sets. The largest pair cluster at 40 % (3,670 rows) is in training; the largest validation cluster has 852 rows in `split_c40` (19.4 % of its rows, S-03) and 198 in `split_c60` (5.0 %).
- **Bias / limitation:** A cluster that mixes homomeric and heteromeric pairs counts in its majority stratum. A cluster that does not fit is still taken if that leaves the stratum closer to its target, so one large family can push a stratum well past it (S-03). The draw depends on the whole list of clusters: the 102 pairs that the conversion of selenomethionine kept in step 08b (I-02) left only 282 of the previous 3,989 `split_c40` validation rows in validation.
- **Configuration:** `splits.seed` (48; the first scheme uses 48, the next 49, and so on).
- **Status:** `implemented`

### S-05 · Allow a protein in training and validation with different partners
- **Decision:** Because the splits are by pair cluster (S-01), a protein can appear in training and in validation with different partners. `docs/NOTES.md` declares it among the known issues.
- **Reason:** This is acceptable: validation can then include cases close to training, known proteins with new partners.
- **Precedent:** Park and Marcotte show that estimates from test pairs whose proteins are both in training do not generalise to new proteins, and accept that setting provided "this limitation should then be noted" [Park2012]. Bernett et al. share no protein between training and test in their gold standard, and show that when proteins are shared, models learn sequence similarity and node degree [Bernett2024].
- **Effect:** In `split_c40`, 899 of the 4,393 validation pairs have both chain clusters at 40 % in training, 103 one and 3,391 (77 %) neither; in `split_c60`, 2,095, 295 and 1,586 of 3,976.
- **Bias / limitation:** Validation mixes the three situations of the test classes (T-01) in proportions that depend on the scheme, so its metrics average over them; most `split_c40` validation pairs are of proteins absent from training.
- **Status:** `implemented`

### S-06 · Weight rows by 1/√ of their family size at 40 %
- **Decision:** Each train/validation row gets `family_weight_c40` = 1/n^0.5, where n is the number of rows of its pair cluster at 40 % in the training side of the exported scheme (`split_c40`; n = 1 for clusters absent from training). The weight is an output column; using it is the user's choice.
- **Reason:** Large families are down-weighted without being removed (R-04). With 1/√n a family of n rows weighs √n, between no weighting (n) and 1/n (1). The family level, 40 %, is the one AlphaFold 2 and 3 use for their clusters.
- **Precedent:** AlphaFold 2 accepts each chain "with the probability inverse to the size of the cluster" (clusters at 40 %) [Jumper2021]; AlphaFold 3 samples chains and interfaces with weights ∝ 1/N_clust [Abramson2024]; ProteinMPNN picks one random member of each 30 % cluster per epoch, which is equivalent to 1/N [Dauparas2022]. The exponent 0.5 has a precedent outside proteins: XLM samples languages with α = 0.5 [Lample2019]. No protein work using 1/√N was found.
- **Effect:** The largest training family at 40 % (3,670 rows) gets 0.017 per row and weighs about 61 times a singleton; the median weight is 0.32. Previous behaviour (weight on pair clusters at 90 % in `split_c90`): the largest training cluster (198 rows) got 0.071 per row, so it weighed about 14 times a singleton; with 1/n both would weigh the same.
- **Bias / limitation:** Single-linkage clusters at 40 % can be large (the largest component at 40 % holds 16.0 % of the sequences), so a few families get very small weights per row.
- **Configuration:** `splits.validation_split` (`split_c40`), `splits.weight_exponent` (0.5; 1 gives 1/n).
- **Status:** `implemented`

## N · Hard negatives

### N-01 · Define a hard negative as two chains of one assembly without a PISA interface
- **Decision:** A hard negative is a pair of protein chains of the same biological assembly with no PISA interface between them, of any size. These are negatives for direct contact, not for interaction.
- **Reason:** A negative means "no evidence of contact in the downloaded structures", not a proven non-interaction. In PSI-MI, two chains of the same complex without contact are still a *physical association* (MI:0915: "within the same physical complex … not necessarily in direct contact") [PSI-MI]; IntAct curates an n-ary crystal structure as a physical association [IntActManual].
- **Precedent:** Drew et al. use the same criterion: in the biological assembly analysed with PISA, "Protein pairs with no interface area were considered not contacting, comprising benchmark negative examples" [Drew2017]. Negatome's PDB subset takes chain pairs of biological units with more than two chains and all inter-chain distances above 8 Å [Smialowski2010; Blohm2014]; DirectContacts2 uses co-complexed pairs that do not share an interface as test negatives [Claussen2026].
- **Effect:** 1,611,626 candidate pairs before the filters N-04 to N-11; 492 hard negatives at the end.
- **Bias / limitation:** Missing contact can come from non-observed density, truncation or modification, and 6.5–8 % of Negatome's structural pairs were later contradicted by PPI data [Smialowski2010]. Fang et al. treat such "mediated" pairs as positives and find coevolution signal in them [Fang2024].
- **Status:** `implemented`

### N-02 · Leave pairs with a sub-threshold interface unlabelled
- **Decision:** A chain pair with a PISA interface that fails the interface thresholds (area, p-value) is neither a positive nor a negative: any PISA interface, of any size, excludes the pair from the negatives.
- **Reason:** Interface size alone does not separate biological from crystal contacts: crystal contacts above 800 Å² form "crystal dimers that cannot be distinguished from real dimers on the basis of the interface size or symmetry" [Bahadur2004]. A small or weak interface is therefore ambiguous.
- **Precedent:** Negatome also drops chain pairs that are nearest neighbours, even without contact [Smialowski2010]. By contrast, Drew et al. call any interface area above 0 Å² a positive, with no grey zone [Drew2017].
- **Status:** `implemented`

### N-03 · Draw negatives from the downloaded assemblies and cross-check them with IntAct
- **Decision:** Negatives come only from the assemblies downloaded for the positives (those with at least one positive candidate). Each negative is then looked up in IntAct through PSICQUIC, for every pair of UniProt accessions of its chains (F-09): `intact_evidences` counts all evidences and `intact_direct_evidences` the binary ones (not spoke- or matrix-expanded) of type direct interaction (MI:0407) or any of its PSI-MI descendants, read from OLS. Negatives with direct evidence are dropped. The answers are cached with their date.
- **Reason:** The step reuses the structures and PISA files already downloaded, and the conflict filter (N-07) only sees these assemblies, so an external source is needed to catch negatives contradicted by known interactions. In PSI-MI, two chains of the same complex without contact are a physical association (MI:0915) [PSI-MI], so association evidence is expected for these negatives and only direct-interaction evidence contradicts them.
- **Precedent:** Negatome removes the pairs that IntAct reports as interacting: 8 % in Negatome 1.0 and 5 % of the structural pairs in Negatome 2.0 [Smialowski2010; Blohm2014]. These resources draw negatives from all multi-chain complexes, not only from those with a positive [Smialowski2010; Blohm2014; Claussen2026].
- **Effect:** 474 of the 492 negatives have a UniProt accession on both chains, 472 of them a pair of different accessions; their 229 distinct pairs were all queried on 10 November 2026, without failures. 50 pairs (141 negatives) have IntAct evidence: association (MI:0914), physical association (MI:0915) and proximity (MI:2364), plus one direct interaction (MI:0407) for O75865–Q9Y296 that comes from a spoke expansion of an n-ary interaction and therefore does not count. No negative is dropped. 173 of the 174 IntAct records for these pairs are expanded from n-ary complexes (spoke or matrix), so the check finds no contradiction rather than confirming that the chains do not interact. Candidates come from the 66,807 downloaded assemblies out of 245,757 PISA assemblies.
- **Bias / limitation:** Complexes without a positive candidate contribute no negatives. IntAct covers curated literature, mostly model organisms: a negative without evidence is not proven non-interacting, and the result depends on the IntAct release at the query date.
- **Configuration:** `hard_negatives.intact.enabled` (true), `hard_negatives.intact.drop_direct` (true); cache in `paths.intact_cache`.
- **Status:** `implemented`

### N-04 · Exclude pairs of identical sequences
- **Decision:** Pairs of chains with identical sequences are not negatives. The filter uses sequence identity, not the homo definition of step 06.
- **Reason:** A copy of a protein that does not touch another copy usually contradicts a homo interface elsewhere in the complex.
- **Precedent:** Negatome drops pairs that map to the same UniProt accession [Smialowski2010].
- **Effect:** 1,611,626 → 771,147 candidates.
- **Bias / limitation:** 79.6 % of the train/validation positives are homomeric and no negative joins identical sequences, so "identical sequences → positive" is a shortcut a model can learn [Park2012; Bernett2024]. Near-identical chains pass this filter; the conflict filter (N-07) usually removes them.
- **Configuration:** `hard_negatives.include_homo` (false).
- **Status:** `implemented`

### N-05 · Use only proteins present in the positives
- **Decision:** Both chains of a negative must have a sequence present in the clustered positives (step 09).
- **Reason:** The negative shares the protein universe, the cluster IDs and the splits of the positives.
- **Precedent:** If a protein appears only in positives, a model learns the protein and not the pair; negatives should use the same proteins as the positives, ideally preserving each protein's degree [Yu2010; Park2012; Bernett2024].
- **Effect:** 771,147 → 315,257 candidates.
- **Bias / limitation:** The same proteins do not mean the same degree: with 492 negatives, almost all proteins still appear only in positives (N-13).
- **Status:** `implemented`

### N-06 · Apply the quality filters of the positives
- **Decision:** Negatives pass the entry and chain thresholds of step 07 (method, R-free, resolution, coverage, chain length) and its scope filters (release date, membrane proteins, antibodies on request), but not the interface thresholds.
- **Reason:** Both classes are filtered alike, so that quality alone does not separate them. Coverage also guards against contacts hidden in unresolved regions.
- **Precedent:** Negatives filtered differently from positives make the task artificially easier [BenHur2006]. DirectContacts2 requires ≥ 80 % coverage so that no pair is "interacting through unresolved protein domains" [Claussen2026].
- **Effect:** 315,257 → 304,165 candidates.
- **Bias / limitation:** Coverage is resolved residues over SEQRES, that is, of the crystallised construct and not of the full protein. A single-domain construct can hide the real contact; Negatome names truncation as one cause of missing contacts [Smialowski2010].
- **Configuration:** The `quality_filter` thresholds of step 07.
- **Status:** `implemented`

### N-07 · Drop a negative if its pair of clusters at 40 % is in contact
- **Decision:** A negative is dropped if its pair of chain clusters at 40 % identity has a PISA interface in any downloaded assembly.
- **Reason:** A precaution. An interface between homologues at 40 % does not prove that this pair interacts, but it makes the negative doubtful; losing negatives is preferred to mislabelling them.
- **Precedent:** Dropping negatives contradicted by a known interaction is what Negatome does with IntAct [Smialowski2010]. No dataset was found that filters negatives by homology, so the 40 % value has no precedent. Close homologues (30–40 % identity or more) that both interact "almost invariably interact the same way" [Aloy2003]: this concerns the mode of interaction, not whether a homologue interacts. Transferring an interaction (interologs) needs a joint identity above 80 % [Yu2004].
- **Effect:** 304,165 → 23,671 candidates (92 % removed).
- **Bias / limitation:** Very conservative: "few pairs interact at JI < 40%" [Yu2004]. Only the downloaded assemblies are checked (N-03), and chains whose sequence is not among the positives clustered in step 09 keep their own sequence as key, so a contact between their homologues counts only for identical sequences.
- **Configuration:** `hard_negatives.conflict_level` (c40; any clustering level, or `exact`).
- **Status:** `implemented`

### N-08 · Apply the special-residue rules of step 08
- **Decision:** A negative is dropped if either chain has a polymer residue not allowed in step 08 (the same rules as for the positives, Q-08).
- **Reason:** Both classes get the same filters, so that chain chemistry alone does not separate them.
- **Precedent:** Filtering negatives differently from positives changes their distribution and makes the task artificially easier [BenHur2006].
- **Effect:** 23,671 → 23,544 candidates.
- **Configuration:** `special_residues.allowed` and `special_residues.ignore` (shared with step 08).
- **Status:** `implemented`

### N-09 · Require at least 6 Å between the two chains
- **Decision:** The minimum heavy-atom distance between the two chains must be at least 6 Å. There is no maximum by default; `max_distance` can keep only near-miss pairs.
- **Reason:** The value is the user's choice, and the default is restrictive: it leaves a margin above heavy-atom contact cut-offs, the pipeline's own (I-13) and the 5 Å of CAPRI and DockQ [Basu2016].
- **Precedent:** No published 6 Å minimum for negatives was found. Negatome requires Cβ–Cβ distances above 8 Å and also drops nearest neighbours [Smialowski2010]; Drew et al. use no distance, only a zero interface area [Drew2017].
- **Effect:** 23,544 → 23,511 candidates. In the final negatives the median minimum distance is 28.4 Å, and 55 are closer than 12 Å.
- **Bias / limitation:** No literature supports near-miss negatives (`max_distance`); Negatome avoids them on purpose by dropping nearest neighbours [Smialowski2010].
- **Configuration:** `hard_negatives.min_distance` (6.0 Å), `hard_negatives.max_distance` (null).
- **Status:** `implemented`

### N-10 · Keep one negative per entry and sequence pair
- **Decision:** Among the copies of the same sequence pair in one entry, only the best-covered one is kept.
- **Reason:** Copies of the same pair are redundant, and redundancy inflates estimates; the usual practice counts protein pairs, not copies.
- **Precedent:** Negatome and DirectContacts2 keep one record per protein pair [Smialowski2010; Claussen2026]; a strictly non-redundant training slightly improved difficult cases [Hamp2015].
- **Effect:** 23,511 → 6,816 candidates.
- **Status:** `implemented`

### N-11 · Cap negatives per pair cluster (4) and per entry (5)
- **Decision:** At most 4 negatives per pair cluster at 95 % across entries, then at most 5 per entry, sampled at random with a fixed seed.
- **Reason:** Otherwise the same complex deposited many times, and large complexes, would dominate. The values are pragmatic.
- **Precedent:** Pairs without contact concentrate in large complexes: complexes of up to 5 chains are usually fully connected and larger ones are not [Hakes2007]; Fang et al. count 3,273 mediated against 414 direct pairs in their PDB benchmark (*E. coli*) [Fang2024]. No published cap was found; DirectContacts2 handles the problem in the split, grouping complexes by the similarity of their subunits [Claussen2026].
- **Effect:** 6,816 → 807 after the cap per pair cluster → 492 after the cap per entry. The median assembly of a negative has 28 protein chains.
- **Bias / limitation:** The caps are not enough: 171 negatives (35 %) still come from entries with "proteasome" in the title and 116 (24 %) from entries whose title names a ribosome or a ribosomal subunit; in the test set, 156 of the 289 (54 %) come from one of the two.
- **Configuration:** `hard_negatives.max_per_pair_cluster` (4), `hard_negatives.pair_cluster_key` (c95), `hard_negatives.max_per_entry` (5), `hard_negatives.seed` (13).
- **Status:** `implemented`

### N-12 · Assign negatives to splits by their chains, one family at a time
- **Decision:** For every split scheme, a negative goes to the test set if a chain appears only in the test set or belongs to a chain cluster held out by step 11 (the unseen proteins of C2 and C3); to validation if a chain appears only in validation; otherwise to training. All negatives of the same pair cluster at 40 % then take the most restrictive split of the group (test > validation > training). Proteins shared by training and the test set do not pull a negative into the test set.
- **Reason:** For a test pair it matters whether its proteins resemble proteins of the training negatives, not only of the positives [Hamp2015]; negatives therefore follow the separation of the positives, and a family of negatives is never split between training and evaluation.
- **Precedent:** Bernett et al. apply the same partition of the proteome to positives and negatives [Bernett2024]; DirectContacts2 keeps complexes, pairs and proteins disjoint between training and test [Claussen2026].
- **Effect:** In `split_c40`, 155 training, 48 validation and 289 test negatives (`split_c60`: 137, 66 and 289). The 289 test negatives come from 77 families at 40 %, the largest with 37; 133 come from proteasome entries. 0 leaks at pair level at 40 % (U-01).
- **Bias / limitation:** The test negatives are concentrated in few families and none joins identical sequences (2 join chains of one cluster at 95 %): evaluate them per family and against heteromeric positives only. Most negatives (289 of 492) are in the test set, because a negative with a protein held out by step 11, or of the same family as one, goes there.
- **Configuration:** `hard_negatives.family_level` (c40, the same as `conflict_level`).
- **Status:** `implemented`

### N-13 · Keep the negatives unbalanced (about 1:81)
- **Decision:** All hard negatives that pass the filters are kept, and their number is not balanced against the positives.
- **Reason:** For now, this is the number of hard negatives available with these filters.
- **Precedent:** Interacting pairs are under 1 % of all protein pairs [Smialowski2010] and about 0.25 % of yeast pairs [BenHur2006]. Benchmarks use 1:1 [BenHur2006; Bernett2024] or more negatives than positives: 1:1, 1:4, 1:9 and 1:332 in Dunham and Ganapathiraju [Dunham2021]; 1,472 positives against 47,246 negatives in the training set of DirectContacts2 [Claussen2026]. Which method looks best depends on the expected ratio of negatives to positives [Hamp2015].
- **Effect:** 492 negatives against 39,747 train/validation positives, about 1:81; within training and validation (`split_c40`), 203 negatives, about 1:196.
- **Bias / limitation:** A declared limitation. The ratio is inverted with respect to every precedent, and almost all proteins appear only in positives, which is the degree shortcut [Bernett2024]. The hard negatives are a difficult evaluation set, not the full negative class; sampling negatives for evaluation and for training are different problems [Park2011]. Training a binary classifier would need random negatives that preserve each protein's degree [Yu2010; Bernett2024].
- **Status:** `implemented`

## U · Leakage audit

### U-01 · Define a leak as a homologue at the set's design threshold
- **Decision:** Step 14 counts as a leak an evaluation chain or pair with a homologue in its reference data at the set's design threshold. At chain level, the chain has a homologous reference chain; at pair level, a reference pair joins homologues of both chains. C3 is audited at chain level (both chains), C2 at chain level (its unseen chain) and at pair level, C1 at pair level, all at `test_set.level`; the validation sets at pair level, at their scheme's identity; the test negatives at pair level, at `conflict_level`. Each test class is also audited against the training and validation hard negatives (U-03).
- **Reason:** The level follows what each set is meant to measure: a new protein (chain level) or a new pair (pair level).
- **Precedent:** AlphaFold 3 evaluates chains with less than 40 % identity to training, and leaves out an interface if both of its chains have more than 40 % identity to two chains of the same training complex [Abramson2024].
- **Effect:** 0 leaks in the 11 audited sets at their design threshold. `reports/14_leakage_by_threshold.tsv` gives every set at every clustering threshold.
- **Status:** `implemented`

### U-02 · Audit with the clustering criteria, plus an optional independent control
- **Decision:** Step 14 searches with MMseqs2 in both directions and links chains with the identity, coverage and E-value criteria of step 09; `tests/check_outputs.py` enforces zero leaks by this criterion. The optional rule `leakage_control` (`analysis/leakage_control.py`) audits the same sets against the same reference data with global alignments (EMBOSS needleall with GraphPart's defaults) and two identity definitions: identical positions over the alignment length, gaps included (GraphPart's default), and over the length of the evaluation chain (AlphaFold 3), the latter counted only when the two chains also have a significant local alignment (MMseqs2, E ≤ 10⁻³, any coverage). The control reports; it does not change the data.
- **Reason:** The clustering uses connected components and the audit looks for direct hits between sets, so step 14 is a second search that catches links the clustering search missed or a split that crosses a link (a hit that both searches miss stays unseen); only a measure that the splits were not built with shows what the linking rule leaves out. How strict the separation must be depends on what a model has to generalise to, so the control measures it and the threshold stays a choice of each project (H-04).
- **Precedent:** GraphPart evaluates separation with global Needleman–Wunsch identities and finds that MMseqs2 and CD-HIT reductions "can leave a considerable amount of cross-partition NW identities above the threshold" [Teufel2023]; PLINDER also measures leakage with metrics not used to split [Durairaj2024]; PINDER adds structural similarity (Foldseek, iAlign) [Kovtun2024].
- **Effect:** Step 14: 0 leaks at the design thresholds in the 11 audited sets. The control, by GraphPart's identity: 3 of the 435 C3 chains, 3 of the 292 unseen C2 chains, 1 of the 150 C1 pairs, 2 of the 300 C2 pairs, 3 of the 3,976 `split_c60` and 43 of the 4,393 `split_c40` validation pairs, and none against the negatives. Each has an MMseqs2 alignment with its reference (E ≤ 4 × 10⁻⁹, identity ≥ 0.46) that step 09 does not link; all but one of these alignments cover less than 90 % of both chains, mostly a shared region, in some cases the same protein expressed as another construct (expression tags, linkers, termini). By AlphaFold 3's identity with a significant alignment: 20 C3 chains, 13 unseen C2 chains, 1 C1 and 4 C2 pairs, 39 and 95 validation pairs and 4 of the 289 test negatives, plus 2 C3 and 2 unseen C2 chains against the negatives (H-04); without the significance requirement, short chains reach 40 % by chance inside much longer ones (82 C3 chains, 64 unseen C2 chains and 228 `split_c40` validation pairs). 11 C3 chains and 18 unseen C2 chains share a UniProt accession with a training or validation chain; 2 of each are among the leaks by AlphaFold 3's identity. In the run of 2 November 2026, with `max_seqs` = 1,000, step 14 itself found real leaks: 2 in the test set, 15, 5 and 5 in the three validation sets, and 13 of 44 test negatives.
- **Bias / limitation:** The control measures sequence only: remote homologues below about 30 % identity and structural similarity are not measured (H-01). Its parasail prefilter (margin 0.35) was checked against needleall on every pair of the test configuration only. The identity of AlphaFold 3 is applied with a semi-global alignment, since the paper does not state the aligner.
- **Configuration:** Optional rule `leakage_control`, in its own environment (`envs/leakage_control.yml`: parasail, EMBOSS 6.6.0, MMseqs2); about 1.5 h on 20 cores for the full run.
- **Status:** `implemented`

### U-03 · Compare each evaluation set with everything the model trains on
- **Decision:** The test set is compared with all training and validation positives, each validation set with the training positives of its scheme, and the test negatives with the training and validation positives and the non-test negatives. Each test bucket is also compared with the training and validation hard negatives, at its design level (sets `*_vs_negatives`). The validation negatives are not audited: their separation rests on the conflict filter (N-07) and the split by families (N-12).
- **Reason:** For a test pair it matters whether its proteins resemble proteins of the training negatives, not only of the positives [Hamp2015].
- **Precedent:** GraphPart requires that no pair of sequences from different partitions exceeds the threshold [Teufel2023]; Bernett et al. use the same partition for positives and negatives [Bernett2024].
- **Effect:** 6,444 evaluation chains searched against 29,019 reference chains, in both directions; the test classes have 0 leaks against the training and validation negatives.
- **Status:** `implemented`

### U-04 · Fail the output check on any leak at a design threshold
- **Decision:** `tests/check_outputs.py` fails if any audited set has a leak at its design threshold.
- **Reason:** Separation takes precedence over keeping data; the criterion is zero violations.
- **Precedent:** GraphPart: "As we strictly require c to be 0 … achieving separation takes precedence" [Teufel2023]; Petti and Eddy build splits in which each test sequence is less than p % identical to any training sequence [Petti2022].
- **Effect:** The script runs 31 checks, 12 of them on leakage (one checks that every audited set is present) and 1 on the PRODIGY-cryst filter (I-02); all pass.
- **Status:** `implemented`

## P · Provenance and reproducibility

### P-01 · Download manifests with dates and revisions
- **Decision:** The PISA manifest (step 01) stores the PDBe file date of each assembly; the structures manifest (step 04) stores the date each mmCIF file was downloaded (`file_date`); the entry metadata stores the deposit and release dates and the current revision of each entry (`revision`, `revision_date`), from the RCSB Data API. `reports/provenance.json` records the ID lists, the query, the resolved configuration, the code commit and the software versions of each run.
- **Reason:** Data need "detailed provenance" (FAIR R1.2) [Wilkinson2016], and each result should keep track of "the exact parameters and inputs that were used" [Sandve2013].
- **Precedent:** AlphaFold 2 dates its PDB download (28 August 2019) and gives the version of each sequence database [Jumper2021]; QSbio-HET dates its PISA ("as of April 2017") and EPPIC ("downloaded [...] on October 17th, 2019") data [DeyLevy2021].
- **Effect:** `reports/provenance.json` of the reference run records the configuration, the SHA-256 of the 17 ID lists, the software versions (with the R packages and the PRODIGY-cryst environment) and the code commit (6f1c494, without uncommitted changes), from which steps 11–14 ran on 12 November 2026 after the redundant and constant columns were dropped from the exported tables; the remaining columns are identical row by row. Steps 06–10 and the optional rules `prodigy_cryst` and `leakage_control` ran from commit 0192cfa on 11 November; the later commits change only comments in their code, and the optional rules read none of the dropped columns, so they were marked up to date (`snakemake --touch`). Steps 01–05, run on 10 November, were marked up to date instead of querying the live services again: their code differs only in comments, console messages and the HTTP User-Agent, and in how step 03 reads its input, which gives a byte-identical output.
- **Bias / limitation:** The assembly files carry no revision history, so the revision recorded is the one current when step 04 ran; a file downloaded earlier may predate it (compare `file_date` with `revision_date`).
- **Status:** `implemented`

### P-02 · Reruns of steps 01–04 follow the live archive
- **Decision:** Steps 01–04 query live services. A rerun downloads again any PISA file that PDBe has updated and queries the entry metadata again, but keeps the mmCIF files already on disk. To freeze a dataset, `data/raw/` and `data/interim/` are kept and steps 01–04 are not rerun.
- **Reason:** PDBe recomputes PISA-lite weekly with the archive [Leines2026] and the wwPDB versions every entry [wwPDBproc]. Archiving intermediate results ("Record All Intermediate Results", rule 5) is what makes a run reproducible [Sandve2013].
- **Precedent:** PINDER publishes versioned releases (current version `2024-02`) [Kovtun2024]; Dockground rebuilds its bound set weekly with the PDB [Collins2022].
- **Effect:** `data/raw/` takes 27 GB and `data/interim/` 1.2 GB.
- **Bias / limitation:** A rerun can mix new PISA files with older mmCIF files, whose coordinates may no longer be those PISA analysed (`docs/NOTES.md`, Known issues).
- **Status:** `implemented`

### P-03 · Fixed random seeds
- **Decision:** The random steps use fixed seeds from the configuration: 7 for the test set (step 11), 48 for the train/validation splits (step 12; one derived seed per scheme) and 13 for the hard negatives (step 13). Step 11 stops if its seed is null.
- **Reason:** "For Analyses That Include Randomness, Note Underlying Random Seeds" (rule 6) [Sandve2013].
- **Precedent:** AlphaFold 2 trains five models "using different random seeds" [Jumper2021].
- **Configuration:** `test_set.seed` (7), `splits.seed` (48), `hard_negatives.seed` (13).
- **Status:** `implemented`

### P-04 · Exact software environment exported
- **Decision:** `environment.yml` gives version ranges (e.g. `mmseqs2>=15`, `snakemake-minimal>=8,<10`) so that the environment installs on other systems. `environment.lock.yml` (`conda env export`) gives the exact versions the pipeline was run with, and `reports/provenance.json` records the versions of Python, R, MMseqs2 and the main packages of each run.
- **Reason:** "Archive the Exact Versions of All External Programs Used" (rule 3) [Sandve2013]. Providing the output of `conda env export` lets anyone obtain "the exact version of the software used, as well as all its dependencies" [Gruning2018].
- **Precedent:** AlphaFold 2 states the versions of hh-suite (v.3.0-beta.3) and of each database (UniRef90 v.2020_01, Uniclust30 v.2018_08) [Jumper2021]; AlphaFold 3 states HMMER v3.3 [Abramson2024].
- **Effect:** The full run used MMseqs2 18.8cc5c, recorded in `09_mmseqs/summary.json` and `reports/provenance.json`.
- **Bias / limitation:** Installing from `environment.yml` can bring a newer MMseqs2 within the range, which can change the clusters and, with them, the splits; the lock file avoids it.
- **Status:** `implemented`

## References

- **Abramson2024** · Abramson J, et al. (2024). Accurate structure prediction of biomolecular interactions with AlphaFold 3. *Nature* 630:493–500. doi:10.1038/s41586-024-07487-w (Supplementary Information §2.1, 2.5.1, 2.5.3, 2.5.4, 5.8, 6.1, 6.2).
- **AF2-code** · AlphaFold 2 source code, `alphafold/data/mmcif_parsing.py` and `alphafold/common/residue_constants.py`. https://github.com/google-deepmind/alphafold (main branch).
- **Aloy2003** · Aloy P, Ceulemans H, Stark A, Russell RB (2003). The relationship between sequence and interaction divergence in proteins. *J Mol Biol* 332:989–998. doi:10.1016/j.jmb.2003.07.006
- **Bahadur2004** · Bahadur RP, Chakrabarti P, Rodier F, Janin J (2004). A dissection of specific and non-specific protein–protein interfaces. *J Mol Biol* 336:943–955. doi:10.1016/j.jmb.2003.12.073.
- **Baskaran2014** · Baskaran K, Duarte JM, Biyani N, Bliven S, Capitani G (2014). A PDB-wide, evolution-based assessment of protein-protein interfaces. *BMC Struct Biol* 14:22. doi:10.1186/s12900-014-0022-0.
- **Basu2016** · Basu S, Wallner B (2016). DockQ: a quality measure for protein–protein docking models. *PLoS One* 11:e0161879. doi:10.1371/journal.pone.0161879.
- **BenHur2006** · Ben-Hur A, Noble WS. 2006. Choosing negative examples for the prediction of protein-protein interactions. *BMC Bioinformatics* 7(Suppl 1):S2. doi:10.1186/1471-2105-7-S1-S2
- **Bernett2024** · Bernett J, Blumenthal DB, List M (2024). Cracking the black box of deep sequence-based protein–protein interaction prediction. *Brief Bioinform* 25:bbae076. doi:10.1093/bib/bbae076 (PMC10939362)
- **Biopython** · Biopython Tutorial and Cookbook, chapter Bio.PDB. https://biopython.org/docs/latest/Tutorial/chapter_pdb.html.
- **Bliven2018** · Bliven S, Lafita A, Parker A, Capitani G, Duarte JM (2018). Automated evaluation of quaternary structures from protein crystals. *PLoS Comput Biol* 14:e1006104. doi:10.1371/journal.pcbi.1006104
- **Blohm2014** · Blohm P, Frishman G, Smialowski P, et al. 2014. Negatome 2.0: a database of non-interacting proteins derived by literature mining, manual annotation and protein structure analysis. *Nucleic Acids Res* 42:D396–D400. doi:10.1093/nar/gkt1079
- **Brunger1992** · Brünger AT (1992). Free R value: a novel statistical quantity for assessing the accuracy of crystal structures. *Nature* 355:472–475. doi:10.1038/355472a0.
- **Bushuiev2024** · Bushuiev A, Bushuiev R, Kouba P, et al. (2024). Learning to design protein–protein interactions with enhanced generalization. *ICLR 2024*. arXiv:2310.18515
- **Capitani2016** · Capitani G, Duarte JM, Baskaran K, Bliven S, Somody JC (2016). Understanding the fabric of protein crystals: computational classification of biological interfaces and crystal contacts. *Bioinformatics* 32:481–489. doi:10.1093/bioinformatics/btv622.
- **CCD** · Chemical Component Dictionary, via the RCSB Data API (https://data.rcsb.org/rest/v1/core/chemcomp/{ID}), consulted on 10 November 2026.
- **Cheng2014** · Cheng H, et al. (2014). ECOD: an evolutionary classification of protein domains. *PLoS Comput Biol* 10:e1003926. doi:10.1371/journal.pcbi.1003926
- **Chothia1986** · Chothia C, Lesk AM (1986). The relation between the divergence of sequence and structure in proteins. *EMBO J* 5:823–826. doi:10.1002/j.1460-2075.1986.tb04288.x
- **Claussen2026** · Claussen ER, Woodcock-Girard MD, Fischer SN, Drew K. 2026. DirectContacts2: a wiring diagram of human physical protein interactions. *Nat Commun* 17:9037. doi:10.1038/s41467-026-75863-3
- **Collins2022** · Collins KW, Copeland MM, Kotthoff I, Singh A, Kundrotas PJ, Vakser IA (2022). Dockground resource for protein recognition studies. *Protein Sci* 31:e4481. doi:10.1002/pro.4481
- **Dauparas2022** · Dauparas J, Anishchenko I, Bennett N, et al. 2022. Robust deep learning–based protein sequence design using ProteinMPNN. *Science* 378(6615):49–56. doi:10.1126/science.add2187
- **DeyLevy2021** · Dey S, Levy ED (2021). PDB-wide identification of physiological hetero-oligomeric assemblies based on conserved quaternary structure geometry. *Structure* 29:1303–1311.e3. doi:10.1016/j.str.2021.07.012.
- **Drew2017** · Drew K, Müller CL, Bonneau R, Marcotte EM. 2017. Identifying direct contacts between protein complex subunits from their conditional dependence in proteomics datasets. *PLoS Comput Biol* 13:e1005625. doi:10.1371/journal.pcbi.1005625
- **Duarte2012** · Duarte JM, Srebniak A, Schärer MA, Capitani G (2012). Protein interface classification by evolutionary analysis. *BMC Bioinformatics* 13:334. doi:10.1186/1471-2105-13-334.
- **Duarte2013** · Duarte JM, Biyani N, Baskaran K, Capitani G (2013). An analysis of oligomerization interfaces in transmembrane proteins. *BMC Struct Biol* 13:21. doi:10.1186/1472-6807-13-21
- **Dun2023** · Dun C, Pan Q, Jin S, et al. (2023). CrysFormer: Protein structure prediction via 3d Patterson maps and partial structure attention. arXiv:2310.03899 (preprint).
- **Dunbar2011** · Dunbar JB, Smith RD, Yang CY, et al. (2011). CSAR benchmark exercise of 2010: selection of the protein-ligand complexes. *J Chem Inf Model* 51:2036–2046. doi:10.1021/ci200082t.
- **Dunham2021** · Dunham B, Ganapathiraju MK. 2021. Benchmark evaluation of protein-protein interaction prediction algorithms. *Molecules* 27:41. doi:10.3390/molecules27010041
- **Durairaj2024** · Durairaj J, Adeshina Y, Cao Z, et al. (2024). PLINDER: The protein-ligand interactions dataset and evaluation resource. *bioRxiv*. doi:10.1101/2024.07.17.603955.
- **Elez2018** · Elez K, Bonvin AMJJ, Vangone A (2018). Distinguishing crystallographic from biological interfaces in protein complexes: role of intermolecular contacts and energetics for classification. *BMC Bioinformatics* 19(Suppl 15):438. doi:10.1186/s12859-018-2414-9.
- **Esmaielbeiki2016** · Esmaielbeiki R, Krawczyk K, Knapp B, Nebel JC, Deane CM (2016). Progress and challenges in predicting protein interfaces. *Brief Bioinform* 17:117–131. doi:10.1093/bib/bbv027.
- **Evans2021** · Evans R *et al.* (2021). Protein complex prediction with AlphaFold-Multimer. *bioRxiv*. doi:10.1101/2021.10.04.463034 (v2)
- **Fang2024** · Fang T, Szklarczyk D, Hachilif R, von Mering C. 2024. Enhancing coevolutionary signals in protein-protein interaction prediction through clade-wise alignment integration. *Sci Rep* 14:6009. doi:10.1038/s41598-024-55655-9
- **Gao2010** · Gao M, Skolnick J (2010). Structural space of protein–protein interfaces is degenerate, close to complete, and highly connected. *PNAS* 107:22517–22522. doi:10.1073/pnas.1012820107
- **Gore2017** · Gore S, Sanz García E, Hendrickx PMS, et al. (2017). Validation of structures in the Protein Data Bank. *Structure* 25:1916–1927. doi:10.1016/j.str.2017.10.009.
- **Gruning2018** · Grüning B *et al.* (2018). Practical computational reproducibility in the life sciences. *Cell Syst* 6:631–635. doi:10.1016/j.cels.2018.03.014
- **Hahn2021** · Hahn DF, Bayly CI, Bruce Macdonald HE, et al. (2021). Best practices for constructing, preparing, and evaluating protein-ligand binding affinity benchmarks. arXiv:2105.06222v2.
- **Hakes2007** · Hakes L, Robertson DL, Oliver SG, Lovell SC. 2007. Protein interactions from complexes: a structural perspective. *Comp Funct Genomics* 2007:49356. doi:10.1155/2007/49356
- **Hamp2015** · Hamp T, Rost B (2015). More challenges for machine-learning protein interactions. *Bioinformatics* 31:1521–1525. doi:10.1093/bioinformatics/btu857
- **Hobohm1992** · Hobohm U, Scharf M, Schneider R, Sander C (1992). Selection of representative protein data sets. *Protein Sci* 1:409–417. doi:10.1002/pro.5560010313
- **IntActManual** · IntAct Curation Manual, sections 4.6.9 *negative* and 5.4.2–5.4.3. https://intact-portal.gitbook.io/intact-portal-doc/docs/curation-manual
- **Janin2007** · Janin J, Rodier F, Chakrabarti P, Bahadur RP (2007). Macromolecular recognition in the Protein Data Bank. *Acta Cryst D* 63:1–8. doi:10.1107/S090744490603575X.
- **Jankauskaite2019** · Jankauskaitė J, Jiménez-García B, Dapkūnas J, Fernández-Recio J, Moal IH (2019). SKEMPI 2.0: an updated benchmark of changes in protein–protein binding energy, kinetics and thermodynamics upon mutation. *Bioinformatics* 35:462–469. doi:10.1093/bioinformatics/bty635
- **JimenezGarcia2019** · Jiménez-García B, Elez K, Koukos PI, Bonvin AMJJ, Vangone A (2019). PRODIGY-crystal: a web-tool for classification of biological interfaces in protein complexes. *Bioinformatics* 35(22):4821–4823. doi:10.1093/bioinformatics/btz437
- **Joeres2025** · Joeres R, Blumenthal DB, Kalinina OV (2025). Data splitting to avoid information leakage with DataSAIL. *Nat Commun* 16:3337. doi:10.1038/s41467-025-58606-8 (PMC11978981)
- **JonesThornton1996** · Jones S, Thornton JM (1996). Principles of protein-protein interactions. *Proc Natl Acad Sci USA* 93:13–20. doi:10.1073/pnas.93.1.13
- **Jumper2021** · Jumper J, et al. (2021). Highly accurate protein structure prediction with AlphaFold. *Nature* 596:583–589. doi:10.1038/s41586-021-03819-2 (Supplementary Information §1.2.5).
- **Kotthoff2022** · Kotthoff I, Kundrotas PJ, Vakser IA (2022). DOCKGROUND membrane protein-protein set. *PLoS One* 17:e0267531. doi:10.1371/journal.pone.0267531 (PMC9113569)
- **Kovtun2024** · Kovtun D *et al.* (2024). PINDER: the protein interaction dataset and evaluation resource. *bioRxiv*. doi:10.1101/2024.07.17.603980 (v1); documentation: https://pinder-org.github.io/pinder/readme.html
- **KrissinelHenrick2007** · Krissinel E, Henrick K (2007). Inference of macromolecular assemblies from crystalline state. *J Mol Biol* 372:774–797. doi:10.1016/j.jmb.2007.05.022
- **Lample2019** · Lample G, Conneau A. 2019. Cross-lingual language model pretraining. arXiv:1901.07291
- **Lawson2021** · Lawson CL *et al.* (2021). Cryo-EM model validation recommendations based on outcomes of the 2019 EMDataResource challenge. *Nat Methods* 18:156–164. doi:10.1038/s41592-020-01051-w
- **Leines2026** · Leines GD, Magaña P, Nair S, Fleming J, Varadi M, Krissinel E, Velankar S (2026). PDBe PISA: enhanced analysis of macromolecular interactions. *Curr Protoc* 6:e70446. doi:10.1002/cpz1.70446.
- **Levy2006** · Levy ED, Pereira-Leal JB, Chothia C, Teichmann SA (2006). 3D complex: a structural classification of protein complexes. *PLoS Comput Biol* 2:e155. doi:10.1371/journal.pcbi.0020155.
- **LoConte1999** · Lo Conte L, Chothia C, Janin J (1999). The atomic structure of protein–protein recognition sites. *J Mol Biol* 285:2177–2198. doi:10.1006/jmbi.1998.2439.
- **MacKinnon2013** · Mackinnon SS, Malevanets A, Wodak SJ (2013). Intertwined associations in structures of homooligomeric proteins. *Structure* 21:638–649. doi:10.1016/j.str.2013.01.019.
- **mmCIFdic** · wwPDB. PDBx/mmCIF dictionary v5 (items `_exptl.method`, `_refine.pdbx_refine_id`, `_entity_poly.pdbx_seq_one_letter_code_can`; category `pdbx_audit_revision_history`). https://mmcif.wwpdb.org/
- **MMseqs-src** · soedinglab/MMseqs2 source code, master branch: `src/commons/Parameters.cpp`, `src/workflow/EasySearch.cpp`, `src/alignment/Alignment.cpp`. https://github.com/soedinglab/MMseqs2
- **MMseqs-wiki** · MMseqs2 User Guide (GitHub wiki, `Home.md`, accessed 10 November 2026). https://github.com/soedinglab/MMseqs2/wiki
- **Morehead2023** · Morehead A, Chen C, Sedova A, Cheng J (2023). DIPS-Plus: The enhanced database of interacting protein structures for interface prediction. *Sci Data* 10:509. doi:10.1038/s41597-023-02409-3.
- **Mosca2014** · Mosca R, Céol A, Stein A, Olivella R, Aloy P (2014). 3did: a catalog of domain-based interactions of known three-dimensional structure. *Nucleic Acids Res* 42:D374–D379. doi:10.1093/nar/gkt887
- **Nooren2003b** · Nooren IMA, Thornton JM (2003). Structural characterisation and functional significance of transient protein–protein interactions. *J Mol Biol* 325:991–1018. doi:10.1016/S0022-2836(02)01281-0.
- **OfranRost2003** · Ofran Y, Rost B (2003). Analysing six types of protein-protein interfaces. *J Mol Biol* 325:377–387. doi:10.1016/s0022-2836(02)01223-8
- **Park2011** · Park Y, Marcotte EM. 2011. Revisiting the negative example sampling problem for predicting protein-protein interactions. *Bioinformatics* 27:3024–3028. doi:10.1093/bioinformatics/btr514
- **Park2012** · Park Y, Marcotte EM (2012). Flaws in evaluation schemes for pair-input computational predictions. *Nat Methods* 9:1134–1136. doi:10.1038/nmeth.2259 (PMC3531800)
- **PDBe-API** · PDBe REST API, endpoint `/mappings/best_structures`; OpenAPI specification https://www.ebi.ac.uk/pdbe/api/v2/openapi.json (accessed 10 November 2026).
- **PDBe2023** · PDBe (25 May 2023). Improved macromolecular interactions data with PISA-lite. https://www.ebi.ac.uk/pdbe/node/38590
- **PDBx** · wwPDB PDBx/mmCIF dictionary v5. https://mmcif.wwpdb.org (items `_entity_poly.pdbx_seq_one_letter_code_can`, `_atom_site.occupancy`, `_atom_site.label_alt_id`, `pdbx_unobs_or_zero_occ_residues`).
- **Petti2022** · Petti S, Eddy SR (2022). Constructing benchmark test sets for biological sequence analysis using independent set algorithms. *PLoS Comput Biol* 18:e1009492. doi:10.1371/journal.pcbi.1009492 (PMC8929697)
- **PINDER-code** · PINDER source code, `pinder/data/config.py`, class `ClusterConfig`. https://github.com/pinder-org/pinder (main branch).
- **PISA-iarea** · PDBe PISA help, "Interface area". https://www.ebi.ac.uk/msd-srv/prot_int/pi_ilist_iarea.html.
- **PISA-idg** · PDBe PISA help, "ΔiG". https://www.ebi.ac.uk/msd-srv/prot_int/pi_ilist_idg.html.
- **PISA-ipv** · PDBe PISA help, "ΔiG P-value". https://www.ebi.ac.uk/msd-srv/prot_int/pi_ilist_ipv.html.
- **PISA-tut** · PDBe, PDBePISA tutorial (PDF). https://www.ebi.ac.uk/pdbe/sites/default/files/documents/service_tutorials/PDBePISA.pdf.
- **PISCES-web** · PISCES server and precompiled lists. https://dunbrack.fccc.edu/pisces/ and https://dunbrack.fccc.edu/pisces/download/ (consulted on 10 November 2026).
- **ProteinFlow** · Kozlova E, Valentin A, Khadhraoui A, Nakhaee-Zadeh Gutierrez D (2023). ProteinFlow: a Python library to pre-process protein structure data for deep learning applications. *bioRxiv*. doi:10.1101/2023.09.25.559346. Default values taken from the code (README and `proteinflow/cli.py`, https://github.com/adaptyvbio/ProteinFlow).
- **ProteinMPNN-code** · dauparas/ProteinMPNN, `training/utils.py` and `training/README.md` (GitHub, main branch). https://github.com/dauparas/ProteinMPNN
- **PSI-MI** · PSI-MI ontology (HUPO-PSI), terms MI:0407, MI:0915, MI:0914 and MI:0403, in OLS4. https://www.ebi.ac.uk/ols4/ontologies/mi
- **RCSB-clust** · RCSB PDB, "Sequence-based clustering" (accessed 10 November 2026). https://www.rcsb.org/docs/grouping-structures/sequence-based-clustering
- **RCSBschema** · RCSB PDB Data API, `entry` schema. https://data.rcsb.org/rest/v1/schema/entry (accessed 10 November 2026)
- **RCSBstats** · RCSB PDB Search API, query of 10 November 2026: entries released in 2025 by experimental method (17,599 entries, 7,233 of them EM). Data, not literature.
- **Read2011** · Read RJ, Adams PD, Arendall WB, et al. (2011). A new generation of crystallographic validation tools for the Protein Data Bank. *Structure* 19:1395–1412. doi:10.1016/j.str.2011.08.006.
- **Rose2026** · Rose Y, Brown R, Voigt M, Bhikadiya C, Bittrich S, Duarte JM, Segura J, Green RK, Burley SK (2026). RCSB Protein Data Bank: next-generation advanced search for exploration of experimental structures and computed structure models. *Protein Sci* 35:e70731. doi:10.1002/pro.70731
- **Rost1999** · Rost B (1999). Twilight zone of protein sequence alignments. *Protein Eng* 12:85–94. doi:10.1093/protein/12.2.85
- **Sandve2013** · Sandve GK, Nekrutenko A, Taylor J, Hovig E (2013). Ten simple rules for reproducible computational research. *PLoS Comput Biol* 9:e1003285. doi:10.1371/journal.pcbi.1003285
- **Schneider2022** · Schneider C, Raybould MIJ, Deane CM (2022). SAbDab in the age of biotherapeutics: updates including SAbDab-nano, the nanobody structure tracker. *Nucleic Acids Res* 50:D1368–D1372. doi:10.1093/nar/gkab1050
- **Smialowski2010** · Smialowski P, Pagel P, Wong P, et al. 2010. The Negatome database: a reference set of non-interacting protein pairs. *Nucleic Acids Res* 38:D540–D544. doi:10.1093/nar/gkp1026
- **Steinegger2017** · Steinegger M, Söding J (2017). MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. *Nat Biotechnol* 35:1026–1028. doi:10.1038/nbt.3988
- **Steinegger2018** · Steinegger M, Söding J (2018). Clustering huge protein sequence sets in linear time. *Nat Commun* 9:2542. doi:10.1038/s41467-018-04964-5 (PMC6026198)
- **Test-MMseqs2** · In-house check, not a publication: MMseqs2 18.8cc5c `easy-search` with the step 09 parameters (`--min-seq-id 0.4 -c 0 -e 1e-3 -s 7.5`) against the 33,641 unique sequences of the run of 2 November 2026; 200 fragments per length, recovered at E ≤ 10⁻³ in 0 % (10 aa), 90 % (15 aa), 99 % (18 aa) and 100 % (≥ 25 aa) of cases.
- **Teufel2023** · Teufel F, Gíslason MH, Almagro Armenteros JJ, Johansen AR, Winther O, Nielsen H (2023). GraphPart: homology partitioning for biological sequence analysis. *NAR Genom Bioinform* 5:lqad088. doi:10.1093/nargab/lqad088 (PMC10578201)
- **vanKempen2024** · van Kempen M, Kim SS, Tumescheit C, et al. (2024). Fast and accurate protein structure search with Foldseek. *Nat Biotechnol* 42:243–246. doi:10.1038/s41587-023-01773-0
- **Vreven2015** · Vreven T *et al.* (2015). Updates to the integrated protein-protein interaction benchmarks: docking benchmark version 5 and affinity benchmark version 2. *J Mol Biol* 427:3031–3041. doi:10.1016/j.jmb.2015.07.016
- **Wang2003** · Wang G, Dunbrack RL (2003). PISCES: a protein sequence culling server. *Bioinformatics* 19:1589–1591. doi:10.1093/bioinformatics/btg224
- **Wang2005** · Wang G, Dunbrack RL (2005). PISCES: recent improvements to a PDB sequence culling server. *Nucleic Acids Res* 33:W94–W98. doi:10.1093/nar/gki402 (PMC1160163)
- **Wilkinson2016** · Wilkinson MD *et al.* (2016). The FAIR Guiding Principles for scientific data management and stewardship. *Sci Data* 3:160018. doi:10.1038/sdata.2016.18
- **Williams2022** · Williams CJ, Richardson DC, Richardson JS (2022). The importance of residue-level filtering and the Top2018 best-parts dataset of high-quality protein residues. *Protein Sci* 31:290–300. doi:10.1002/pro.4239.
- **Wlodawer2008** · Wlodawer A, Minor W, Dauter Z, Jaskolski M (2008). Protein crystallography for non-crystallographers, or how to get the best (but not more) from published macromolecular structures. *FEBS J* 275:1–21. doi:10.1111/j.1742-4658.2007.06178.x.
- **wwPDB-VR** · wwPDB. User guide to the X-ray validation report (2016 version). https://www.wwpdb.org/validation/2016/XrayValidationReportHelp.
- **wwPDBextID** · wwPDB. Transition to extended IDs and the PDB Beta Archive. https://wwpdb.org/documentation/new-format-for-pdb-ids
- **wwPDBobs** · wwPDB. Statement on retraction of PDB entries. https://www.wwpdb.org/documentation/UAB
- **wwPDBproc** · wwPDB. Deposition policies and wwPDB biocuration procedures. https://www.wwpdb.org/documentation/procedure
- **XuDunbrack2011** · Xu Q, Dunbrack RL (2011). The protein common interface database (ProtCID)—a comprehensive database of interactions of homologous proteins in multiple crystal forms. *Nucleic Acids Res* 39:D761–D770. doi:10.1093/nar/gkq1059
- **XuDunbrack2019** · Xu Q, Dunbrack RL (2019). Principles and characteristics of biological assemblies in experimentally determined protein structures. *Curr Opin Struct Biol* 55:34–49. doi:10.1016/j.sbi.2019.03.006.
- **XuDunbrack2023** · Xu Q, Dunbrack RL (2023). The protein common assembly database (ProtCAD)—a comprehensive structural resource of protein complexes. *Nucleic Acids Res* 51:D466–D478. doi:10.1093/nar/gkac937.
- **Yu2004** · Yu H, Luscombe NM, Lu HX, et al. 2004. Annotation transfer between genomes: protein-protein interologs and protein-DNA regulogs. *Genome Res* 14:1107–1118. doi:10.1101/gr.1774904
- **Yu2010** · Yu J, Guo M, Needham CJ, Huang Y, Cai L, Westhead DR. 2010. Simple sequence-based kernels do not predict protein-protein interactions. *Bioinformatics* 26:2610–2614. doi:10.1093/bioinformatics/btq483
- **YuGuerois2016** · Yu J, Guerois R (2016). PPI4DOCK: large scale assessment of the use of homology models in free docking over more than 1000 realistic targets. *Bioinformatics* 32:3760–3767. doi:10.1093/bioinformatics/btw533
- **Zhu2006** · Zhu H, Domingues FS, Sommer I, Lengauer T (2006). NOXclass: prediction of protein–protein interaction types. *BMC Bioinformatics* 7:27. doi:10.1186/1471-2105-7-27.
