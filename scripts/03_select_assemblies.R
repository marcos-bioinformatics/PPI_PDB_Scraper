library(tidyr)
library(dplyr)
library(data.table)

source(file.path("scripts", "common.R"))
cfg <- load_config()

# Carga de master_pisa.csv con todos los interfaces
df = read.csv(cfg$paths$master_pisa)

# FIltrado con los 164k PDB_IDs extraidos y descargados del PDB

# df ya cargado:
setDT(df)
df[, pdb_id := tolower(pdb_id)]

# Cargar la lista canónica de PDB IDs (generada en el paso 01)
txts <- cfg$paths$ids

# Lee tokens robustamente (separadores: espacios, comas, tabs, saltos)
read_ids <- function(f) {
  x <- readLines(f, warn = FALSE)
  x <- unlist(strsplit(x, "[,[:space:]]+"))
  x <- tolower(trimws(x))
  x[nzchar(x)]
}

pdb_ok <- unique(unlist(lapply(txts, read_ids)))
# Filtrar a 4 chars y formato PDB (primer char dígito)
pdb_ok <- pdb_ok[nchar(pdb_ok) == 4 & grepl("^[0-9][a-z0-9]{3}$", pdb_ok)]

# Filtrar master
df_164k <- df[pdb_id %in% pdb_ok]
length(unique(df_164k$pdb_id))

# Guardar
fwrite(df_164k, ensure_parent(cfg$paths$master_selected))

need_pairs <- unique(df_164k[, .(pdb_id, assembly_id)])
fwrite(need_pairs, ensure_parent(cfg$paths$needed_pairs))

# ¿Cuántos pdb tienen >1 assembly en el master?
multi <- need_pairs[, .N, by=pdb_id][N > 1]
nrow(multi)

# ¿Qué assembly_ids distintos aparecen?
sort(unique(need_pairs$assembly_id))

# si ya tienes need_pairs en memoria
need_pairs_extra <- need_pairs[assembly_id != 1]
nrow(need_pairs_extra)                 # <- objetivo real de descargas extra
need_pairs_extra[, .N, by=assembly_id][order(assembly_id)]

