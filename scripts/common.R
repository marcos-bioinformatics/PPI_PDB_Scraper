# Shared helpers for the R steps of the BeastAffinity pipeline.
# Usage at the top of each R script (run from the repository root):
#   source(file.path("scripts", "common.R"))
#   cfg <- load_config()

.load_raw_config <- function(path) {
  cfg <- yaml::read_yaml(path)
  if (!is.null(cfg$base)) {
    base <- .load_raw_config(cfg$base)
    cfg$base <- NULL
    cfg <- modifyList(base, cfg)
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
  cfg <- .load_raw_config(path)
  data_dir <- if (is.null(cfg$data_dir)) "data" else cfg$data_dir
  cfg$paths <- lapply(cfg$paths, function(p) gsub("{data_dir}", data_dir, p, fixed = TRUE))
  cfg
}

ensure_parent <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  path
}
