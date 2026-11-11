# Shared helpers for the R steps of the PPI_PDB_Scraper pipeline.
# Usage at the top of each R script (run from the repository root):
#   source(file.path("scripts", "common.R"))
#   cfg <- load_config()

# YAML 1.1 reads values such as 1e-3 as strings: turn numeric strings into numbers.
.coerce_numbers <- function(x) {
  if (is.list(x)) return(lapply(x, .coerce_numbers))
  if (is.character(x) && length(x) == 1 &&
      grepl("^[-+]?([0-9]+\\.?[0-9]*|\\.[0-9]+)([eE][-+]?[0-9]+)?$", trimws(x))) {
    return(as.numeric(x))
  }
  x
}

.load_raw_config <- function(path, seen = character(0)) {
  key <- normalizePath(path, mustWork = TRUE)
  if (key %in% seen) stop("Config inheritance cycle: ", paste(c(seen, key), collapse = " -> "))
  cfg <- yaml::read_yaml(path)
  if (!is.null(cfg$base)) {
    base <- cfg$base
    # relative to the working directory, or else to the config file itself
    if (!file.exists(base) && file.exists(file.path(dirname(path), base))) base <- file.path(dirname(path), base)
    parent <- .load_raw_config(base, c(seen, key))
    cfg$base <- NULL
    cfg <- modifyList(parent, cfg)
  }
  cfg
}

# Reads the YAML config given as first command-line argument (Snakemake passes
# it), or config/config.yaml when the script is run by hand from the repo root.
load_config <- function(path = NULL) {
  if (is.null(path)) {
    args <- commandArgs(trailingOnly = TRUE)
    path <- if (length(args) > 0) args[1] else "config/config.yaml"
  }
  cfg <- .coerce_numbers(.load_raw_config(path))
  data_dir <- if (is.null(cfg$data_dir)) "data" else cfg$data_dir
  cfg$paths <- lapply(cfg$paths, function(p) gsub("{data_dir}", data_dir, p, fixed = TRUE))
  cfg
}

# Optional integer parameter: NULL (null in the YAML) stays NULL, so the
# function receiving it can apply its own default.
int_or_null <- function(x) if (is.null(x)) NULL else as.integer(x)

# Identifier-like columns that must be read as text: a chain named "T", a PDB ID such as
# 1e10 or a revision such as 1.10 would otherwise become logical or numeric. (A chain named
# "NA" is still read as missing: the R steps keep the default na.strings, because missing
# values written by write.csv travel through the Python steps as the text "NA".)
TEXT_COLUMNS <- c("pdb_id", "molA_chain_id", "molB_chain_id", "revision")

# colClasses for the TEXT_COLUMNS present in a CSV file (for read.csv and fread).
text_col_classes <- function(path) {
  present <- intersect(TEXT_COLUMNS, names(data.table::fread(path, nrows = 0)))
  setNames(rep("character", length(present)), present)
}

ensure_parent <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  path
}
