# =============================================================================
# Train/validation split + family weighting protocol
#
# This script builds cluster-aware train/validation splits from the final
# trainval pool obtained after:
#   (i) structural/biophysical filtering,
#   (ii) exact pair deduplication,
#   (iii) light pruning at pair_cluster_98,
#   (iv) removal of the external test set (C1/C2/C3).
#
# The goal is to evaluate model development under validation schemes with
# increasing homology separation, while preserving a realistic homo/hetero
# composition and avoiding cluster leakage between train and validation.
#
# -----------------------------------------------------------------------------
# OVERVIEW
# -----------------------------------------------------------------------------
#
# Three independent train/validation splits are created:
#
#   1) split_c90  -> validation defined at pair_cluster_90
#   2) split_c80  -> validation defined at pair_cluster_80
#   3) split_c60  -> validation defined at pair_cluster_60
#
# These three splits provide progressively harder validation settings:
#
#   - split_c90 is the main / default split and is the closest to the training
#     distribution while still being cluster-aware.
#
#   - split_c80 is an intermediate robustness split, useful to test whether
#     model improvements survive a stricter pair-homology separation.
#
#   - split_c60 is the hardest validation split among the three and is mainly
#     intended as a stronger generalization / robustness check rather than the
#     default split for day-to-day model development.
#
# -----------------------------------------------------------------------------
# SPLIT CONSTRUCTION PRINCIPLE
# -----------------------------------------------------------------------------
#
# Each split is built at the level of complete pair clusters, not individual
# rows. This is critical: all rows belonging to the same pair_cluster_X are
# assigned together to either train or validation. This guarantees:
#
#   - zero pair-cluster overlap between train and validation
#   - cluster-disjoint validation sets
#   - no leakage of homologous pairs across train/val at the chosen threshold
#
# Validation is targeted to approximately 10% of rows (VAL_FRAC = 0.10).
#
# -----------------------------------------------------------------------------
# HOMO / HETERO STRATIFICATION
# -----------------------------------------------------------------------------
#
# To preserve the global composition of the dataset, the split is stratified by
# homomeric vs heteromeric character. However, stratification is performed at
# the cluster level, not the row level:
#
#   - for each pair cluster, the fraction of rows labeled is_homo01 == 1 is
#     computed
#   - the cluster is assigned a majority-based stratum:
#
#         cluster_stratum = 1 if cluster_homo_frac >= 0.5
#         cluster_stratum = 0 otherwise
#
# Clusters are then selected for validation within each stratum until the
# target validation size is approximately reached.
#
# This preserves the homo/hetero balance while keeping cluster assignments
# atomic and leakage-free.
#
# -----------------------------------------------------------------------------
# GREEDY CLUSTER SELECTION
# -----------------------------------------------------------------------------
#
# Within each stratum, clusters are shuffled with a fixed seed and selected
# greedily until the validation target is reached (or slightly approximated).
# The decision to include or skip a cluster is based on closeness to the target
# number of validation rows, rather than simply taking clusters in arbitrary
# order.
#
# This produces:
#
#   - validation fractions very close to the desired target
#   - strong homo/hetero balance preservation
#   - fully cluster-disjoint train and validation subsets
#
# -----------------------------------------------------------------------------
# FAMILY WEIGHTING
# -----------------------------------------------------------------------------
#
# In addition to the train/validation splits, the script computes a family
# weight based on pair_cluster_90:
#
#   family_weight_c90 = 1 / sqrt(cluster_size_train_c90)
#
# where cluster_size_train_c90 is the size of the corresponding pair_cluster_90
# measured using ONLY the TRAIN subset of split_c90.
#
# This is important:
#
#   - validation rows are not used to define training-family sizes
#   - weighting reflects redundancy structure actually seen during training
#   - large families are downweighted without being removed
#
# If a pair_cluster_90 appears only in validation and not in training (rare but
# possible), its training cluster size is set to 1 to avoid missing values or
# infinite weights.
#
# Interpretation:
#
#   - singleton train cluster      -> weight = 1.0
#   - cluster of size 4            -> weight = 0.5
#   - cluster of size 9            -> weight ~ 0.333
#   - cluster of size 100          -> weight = 0.1
#
# This weighting is intended to reduce family dominance while preserving
# examples from recurrent families, and is consistent with the overall strategy
# adopted for dataset curation:
#
#   exact dedupe + light pruning + inverse-family downweighting
#
# -----------------------------------------------------------------------------
# WHY THREE SPLITS?
# -----------------------------------------------------------------------------
#
# The three splits support different experimental goals:
#
#   split_c90:
#     Main model-development split. Best default choice for tuning, ablations,
#     and model selection.
#
#   split_c80:
#     Robustness split. Tests whether improvements remain under stricter
#     pair-homology separation.
#
#   split_c60:
#     Hard split. Useful to probe stronger generalization behavior, but not
#     necessarily the preferred split for routine iteration.
#
# In practice:
#
#   - split_c90 + family_weight_c90 is the primary training configuration
#   - split_c80 and split_c60 are secondary validation settings for sensitivity
#     / robustness analysis
#
# -----------------------------------------------------------------------------
# EXPECTED AUDITS
# -----------------------------------------------------------------------------
#
# The split script should verify at least the following properties:
#
#   1) Row counts:
#        - train / val counts
#        - realized validation fraction
#
#   2) Homo/hetero composition:
#        - global homo fraction
#        - validation homo fraction
#
#   3) Quality comparability:
#        - coverage
#        - p_value
#        - r_free
#        - interface_area
#
#   4) Cluster integrity:
#        - unique clusters in train
#        - unique clusters in val
#        - zero overlap between train and val at the chosen cluster threshold
#
#   5) Validation cluster-size profile:
#        - number of validation clusters
#        - cluster-size quantiles (p50, p90, p95, p99, max)
#        - singleton fraction
#        - fraction of validation rows belonging to clusters of size > 1
#
#   6) Weight distribution:
#        - summary of cluster_size_train_c90
#        - summary of family_weight_c90
#
# -----------------------------------------------------------------------------
# FINAL INTERPRETATION
# -----------------------------------------------------------------------------
#
# The resulting objects have the following roles:
#
#   split_c90 / split_c80 / split_c60:
#     Cluster-aware train/validation labels for three homology thresholds.
#
#   cluster_size_train_c90:
#     Size of each pair_cluster_90 measured only on the TRAIN subset of
#     split_c90.
#
#   family_weight_c90:
#     Per-row family weight derived from the corresponding train-side
#     pair_cluster_90 size.
#
# Recommended usage:
#
#   - default training / model selection:
#         use split_c90 + family_weight_c90
#
#   - robustness checks:
#         repeat with split_c80 and split_c60
#
# This setup allows clean comparison of:
#
#   - model architecture changes
#   - effect of family downweighting
#   - dependence of performance on homology separation strength
#
# while keeping the validation protocol leakage-free and biologically
# interpretable.
# =============================================================================


library(data.table)

source(file.path("scripts", "common.R"))
cfg <- load_config()

IN_CSV  <- cfg$paths$trainval_pool
OUT_CSV <- ensure_parent(cfg$paths$trainval_splits)

VAL_FRAC <- cfg$splits$val_frac
SEED_MASTER <- as.integer(cfg$splits$seed)

to_homo01 <- function(x) {
  x_chr <- tolower(trimws(as.character(x)))
  out <- rep(NA_integer_, length(x_chr))
  out[x_chr %in% c("1","true","t","yes","y")] <- 1L
  out[x_chr %in% c("0","false","f","no","n")] <- 0L
  out
}

make_split_by_cluster_stratified <- function(
    DT,
    cluster_col,
    out_col,
    strat_col = "is_homo01",
    val_frac = 0.10,
    seed = 1L
) {
  stopifnot(is.data.table(DT))
  stopifnot(cluster_col %in% names(DT))
  stopifnot(strat_col %in% names(DT))
  stopifnot(val_frac > 0 && val_frac < 1)
  
  set.seed(seed)
  
  DT[, (out_col) := "train"]
  
  # trabajar a nivel cluster, no a nivel fila
  clust_dt <- DT[!is.na(get(cluster_col)), .(
    cluster_n = .N,
    cluster_homo_frac = mean(get(strat_col), na.rm = TRUE)
  ), by = cluster_col]
  
  setnames(clust_dt, cluster_col, "cluster_id")
  clust_dt[, cluster_stratum := fifelse(cluster_homo_frac >= 0.5, 1L, 0L)]
  
  val_clusters_all <- character(0)
  
  for (s in c(0L, 1L)) {
    sub <- clust_dt[cluster_stratum == s]
    if (nrow(sub) == 0L) next
    
    # target en filas, no en número de clusters
    target_val <- round(val_frac * sum(sub$cluster_n))
    
    sub <- sub[sample(.N)]
    current_n <- 0L
    val_clusters <- character(0)
    
    for (i in seq_len(nrow(sub))) {
      n_i <- sub$cluster_n[i]
      cid <- sub$cluster_id[i]
      
      take_it <- FALSE
      
      if (current_n < target_val) {
        if ((current_n + n_i) <= target_val) {
          take_it <- TRUE
        } else {
          diff_take <- abs((current_n + n_i) - target_val)
          diff_skip <- abs(current_n - target_val)
          if (diff_take <= diff_skip || length(val_clusters) == 0L) {
            take_it <- TRUE
          }
        }
      }
      
      if (take_it) {
        val_clusters <- c(val_clusters, cid)
        current_n <- current_n + n_i
      }
    }
    
    if (length(val_clusters) == 0L && nrow(sub) > 0L) {
      val_clusters <- sub$cluster_id[1]
    }
    
    val_clusters_all <- c(val_clusters_all, val_clusters)
  }
  
  DT[get(cluster_col) %in% val_clusters_all, (out_col) := "val"]
  
  invisible(DT)
}

audit_split <- function(DT, split_col, cluster_col, strat_col = "is_homo01") {
  list(
    counts = DT[, .N, by = split_col],
    homo_frac_global = DT[, mean(get(strat_col), na.rm = TRUE)],
    homo_frac_val = DT[get(split_col) == "val", mean(get(strat_col), na.rm = TRUE)],
    val_frac_real = DT[, mean(get(split_col) == "val")],
    unique_clusters_train = DT[get(split_col) == "train", uniqueN(get(cluster_col))],
    unique_clusters_val = DT[get(split_col) == "val", uniqueN(get(cluster_col))],
    quality = DT[, .(
      n = .N,
      cov_min_pct_median = median(pmin(as.numeric(molA_coverage_pct), as.numeric(molB_coverage_pct), na.rm = TRUE), na.rm = TRUE),
      p_value_median = median(as.numeric(p_value), na.rm = TRUE),
      r_free_median = median(as.numeric(r_free), na.rm = TRUE),
      area_median = median(as.numeric(interface_area), na.rm = TRUE)
    ), by = split_col]
  )
}

df <- fread(IN_CSV)
setDT(df)

stopifnot("row_id" %in% names(df))
stopifnot(uniqueN(df$row_id) == nrow(df))

# usar is_homo01 si ya existe; si no, derivarlo
if (!("is_homo01" %in% names(df))) {
  stop("Falta is_homo01 en el input.")
}
df[, is_homo01 := to_homo01(is_homo01)]

# Splits
make_split_by_cluster_stratified(
  df,
  cluster_col = "pair_cluster_90",
  out_col = "split_c90",
  strat_col = "is_homo01",
  val_frac = VAL_FRAC,
  seed = SEED_MASTER + 0L
)

make_split_by_cluster_stratified(
  df,
  cluster_col = "pair_cluster_80",
  out_col = "split_c80",
  strat_col = "is_homo01",
  val_frac = VAL_FRAC,
  seed = SEED_MASTER + 1L
)

make_split_by_cluster_stratified(
  df,
  cluster_col = "pair_cluster_60",
  out_col = "split_c60",
  strat_col = "is_homo01",
  val_frac = VAL_FRAC,
  seed = SEED_MASTER + 2L
)

# Weighting basado en pair_cluster_90, calculado SOLO con TRAIN de split_c90
SPLIT_FOR_WEIGHT <- "split_c90"
WEIGHT_CLUSTER_COL <- "pair_cluster_90"

stopifnot(SPLIT_FOR_WEIGHT %in% names(df))
stopifnot(WEIGHT_CLUSTER_COL %in% names(df))

train_dt <- df[get(SPLIT_FOR_WEIGHT) == "train", .(pair_cluster_90)]
cl_sizes <- train_dt[, .N, by = pair_cluster_90]
setnames(cl_sizes, "N", "cluster_size_train_c90")

df <- merge(df, cl_sizes, by = "pair_cluster_90", all.x = TRUE, sort = FALSE)
df[is.na(cluster_size_train_c90), cluster_size_train_c90 := 1L]

df[, family_weight_c90 := 1.0 / sqrt(cluster_size_train_c90)]

# Guardar
fwrite(df, OUT_CSV)

cat("OK ->", OUT_CSV, "\n")

# Validation set: rows marked "val" in the chosen split
df_val <- df[get(cfg$splits$validation_split) == "val"]
df_val[, `:=`(test_bucket = "VAL", test_subbucket = "VAL")]
fwrite(df_val, ensure_parent(cfg$paths$validation))
cat("OK ->", cfg$paths$validation, "(", nrow(df_val), "rows )\n")

cat("\n=== AUDIT split_c90 ===\n")
a90 <- audit_split(df, "split_c90", "pair_cluster_90")
print(a90$counts)
cat("homo global:", a90$homo_frac_global, "\n")
cat("homo val:", a90$homo_frac_val, "\n")
cat("val frac real:", a90$val_frac_real, "\n")
cat("unique clusters train:", a90$unique_clusters_train, "\n")
cat("unique clusters val:", a90$unique_clusters_val, "\n")
print(a90$quality)

cat("\n=== AUDIT split_c80 ===\n")
a80 <- audit_split(df, "split_c80", "pair_cluster_80")
print(a80$counts)
cat("homo global:", a80$homo_frac_global, "\n")
cat("homo val:", a80$homo_frac_val, "\n")
cat("val frac real:", a80$val_frac_real, "\n")
cat("unique clusters train:", a80$unique_clusters_train, "\n")
cat("unique clusters val:", a80$unique_clusters_val, "\n")
print(a80$quality)

cat("\n=== AUDIT split_c60 ===\n")
a60 <- audit_split(df, "split_c60", "pair_cluster_60")
print(a60$counts)
cat("homo global:", a60$homo_frac_global, "\n")
cat("homo val:", a60$homo_frac_val, "\n")
cat("val frac real:", a60$val_frac_real, "\n")
cat("unique clusters train:", a60$unique_clusters_train, "\n")
cat("unique clusters val:", a60$unique_clusters_val, "\n")
print(a60$quality)

cat("\nFamily weight (c90) resumen:\n")
print(summary(df$family_weight_c90))
print(summary(df$cluster_size_train_c90))


# =============================================================================
# EXTRA AUDITS
# A) overlap de clusters entre train y val
# B) distribución de tamaños de cluster en val
# =============================================================================

audit_cluster_overlap <- function(DT, split_col, cluster_col) {
  tr <- unique(DT[get(split_col) == "train" & !is.na(get(cluster_col)), get(cluster_col)])
  va <- unique(DT[get(split_col) == "val"   & !is.na(get(cluster_col)), get(cluster_col)])
  
  inter <- intersect(tr, va)
  
  data.table(
    split_col = split_col,
    cluster_col = cluster_col,
    n_train_clusters = length(tr),
    n_val_clusters = length(va),
    n_overlap_clusters = length(inter),
    overlap_frac_vs_val = ifelse(length(va) > 0, length(inter) / length(va), NA_real_)
  )
}

audit_val_cluster_sizes <- function(DT, split_col, cluster_col) {
  val_sizes <- DT[get(split_col) == "val" & !is.na(get(cluster_col)),
                  .N, by = cluster_col]
  
  if (nrow(val_sizes) == 0L) {
    return(data.table(
      split_col = split_col,
      cluster_col = cluster_col,
      n_val_clusters = 0L,
      n_val_rows = 0L,
      p0 = NA_real_,
      p50 = NA_real_,
      p90 = NA_real_,
      p95 = NA_real_,
      p99 = NA_real_,
      p100 = NA_real_,
      mean = NA_real_,
      pct_singletons = NA_real_,
      pct_rows_in_gt1 = NA_real_
    ))
  }
  
  qs <- quantile(val_sizes$N, probs = c(0, 0.5, 0.9, 0.95, 0.99, 1), na.rm = TRUE)
  
  data.table(
    split_col = split_col,
    cluster_col = cluster_col,
    n_val_clusters = nrow(val_sizes),
    n_val_rows = sum(val_sizes$N),
    p0 = unname(qs[1]),
    p50 = unname(qs[2]),
    p90 = unname(qs[3]),
    p95 = unname(qs[4]),
    p99 = unname(qs[5]),
    p100 = unname(qs[6]),
    mean = mean(val_sizes$N),
    pct_singletons = mean(val_sizes$N == 1) * 100,
    pct_rows_in_gt1 = sum(val_sizes$N[val_sizes$N > 1]) / sum(val_sizes$N) * 100
  )
}

# ---- correr para c90 / c80 / c60
overlap_audits <- rbindlist(list(
  audit_cluster_overlap(df, "split_c90", "pair_cluster_90"),
  audit_cluster_overlap(df, "split_c80", "pair_cluster_80"),
  audit_cluster_overlap(df, "split_c60", "pair_cluster_60")
), use.names = TRUE, fill = TRUE)

val_cluster_size_audits <- rbindlist(list(
  audit_val_cluster_sizes(df, "split_c90", "pair_cluster_90"),
  audit_val_cluster_sizes(df, "split_c80", "pair_cluster_80"),
  audit_val_cluster_sizes(df, "split_c60", "pair_cluster_60")
), use.names = TRUE, fill = TRUE)

cat("\n=== AUDIT A: CLUSTER OVERLAP TRAIN/VAL ===\n")
print(overlap_audits)

cat("\n=== AUDIT B: VAL CLUSTER SIZE DISTRIBUTION ===\n")
print(val_cluster_size_audits)

# opcional: guardar
#fwrite(overlap_audits, "audit_cluster_overlap_train_val.csv")
#fwrite(val_cluster_size_audits, "audit_val_cluster_size_distribution.csv")

# Mostrar clusters solapados exactos
show_overlap_clusters <- function(DT, split_col, cluster_col) {
  tr <- unique(DT[get(split_col) == "train" & !is.na(get(cluster_col)), get(cluster_col)])
  va <- unique(DT[get(split_col) == "val"   & !is.na(get(cluster_col)), get(cluster_col)])
  inter <- intersect(tr, va)
  data.table(cluster_id = inter)[1:50]
}

ov90 <- show_overlap_clusters(df, "split_c90", "pair_cluster_90")
ov80 <- show_overlap_clusters(df, "split_c80", "pair_cluster_80")
ov60 <- show_overlap_clusters(df, "split_c60", "pair_cluster_60")

print(ov90)
print(ov80)
print(ov60)
