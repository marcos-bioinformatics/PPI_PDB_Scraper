library(tidyr)
library(dplyr)
library(data.table)

source(file.path("scripts", "common.R"))
cfg <- load_config()

# Load master_pisa.csv with all the interfaces
# IDs as text (a PDB ID such as 1e10 or a chain named T would otherwise be read as a number or a logical)
df = read.csv(cfg$paths$master_pisa, colClasses = text_col_classes(cfg$paths$master_pisa))

# Filter to the PDB IDs extracted and downloaded from the PDB

# df already loaded:
setDT(df)
df[, pdb_id := tolower(pdb_id)]

# Load the canonical list of PDB IDs (generated in step 01)
txts <- cfg$paths$ids

# Read tokens robustly (separators: spaces, commas, tabs, newlines)
read_ids <- function(f) {
  x <- readLines(f, warn = FALSE)
  x <- unlist(strsplit(x, "[,[:space:]]+"))
  x <- tolower(trimws(x))
  x[nzchar(x)]
}

pdb_ok <- unique(unlist(lapply(txts, read_ids)))
# Keep 4-character IDs in PDB format (first character a digit)
pdb_ok <- pdb_ok[nchar(pdb_ok) == 4 & grepl("^[0-9][a-z0-9]{3}$", pdb_ok)]

# Filter master
df_sel <- df[pdb_id %in% pdb_ok]

# Assembly selection: "first" keeps only the first biological assembly of each
# entry (its lowest current assembly ID, recorded by step 01), as AlphaFold 3 and
# PINDER do; "all" keeps every assembly that PISA analysed.
selection <- cfg$assemblies$selection
if (!selection %in% c("first", "all")) stop("assemblies.selection must be 'first' or 'all'")
if (selection == "first") {
  man <- fread(cfg$paths$pisa_manifest, colClasses = list(character = "pdb_id"))
  man[, first_assembly := suppressWarnings(as.integer(first_assembly))]
  first_asm <- unique(man[!is.na(first_assembly), .(pdb_id = tolower(pdb_id), assembly_id = first_assembly)])
  n_before <- uniqueN(df_sel$pdb_id)
  rows_before <- nrow(df_sel)
  df_sel <- df_sel[first_asm, on = .(pdb_id, assembly_id), nomatch = NULL]
  cat(sprintf("assemblies.selection = first: %d of %d candidate rows kept; %d of %d entries have candidates in their first assembly\n",
              nrow(df_sel), rows_before, uniqueN(df_sel$pdb_id), n_before))
}

# Save
fwrite(df_sel, ensure_parent(cfg$paths$master_selected))

need_pairs <- unique(df_sel[, .(pdb_id, assembly_id)])
fwrite(need_pairs, ensure_parent(cfg$paths$needed_pairs))

# How many PDB entries have >1 assembly in the master?
multi <- need_pairs[, .N, by=pdb_id][N > 1]

# if you already have need_pairs in memory
need_pairs_extra <- need_pairs[assembly_id != 1]

