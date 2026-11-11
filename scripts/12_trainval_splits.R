# =============================================================================
# Step 12 - Train/validation splits by pair cluster, and family weights
#
# Input : paths.trainval_pool (step 11: everything not in the test set)
# Output: paths.trainval_splits  one row per pair with a split column per scheme
#         paths.validation       rows marked "val" in splits.validation_split
#
# One independent train/validation split per level in splits.schemes (e.g.
# 0.60 -> split_c60, built on pair_cluster_60). Whole pair clusters go to one
# side, so no training pair joins homologues of both chains of a validation
# pair at that identity. Within each homo/hetero stratum (majority of the
# cluster's rows), clusters are shuffled with a fixed seed (splits.seed + i for
# the i-th scheme) and taken until about splits.val_frac of the rows are in
# validation.
#
# Family weight, computed on the scheme exported as validation
# (splits.validation_split, e.g. split_c40 on pair_cluster_40):
#   cluster_size_train_cXX = rows of the pair cluster in the training side
#                            (1 for clusters absent from training)
#   family_weight_cXX      = 1 / cluster_size_train_cXX ^ splits.weight_exponent
# With the default exponent 0.5 a family of n rows weighs sqrt(n) in total; 1
# gives every family the same weight. Using the weight is up to the user.
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
  
  # work at cluster level, not at row level
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
    
    # target in rows, not in number of clusters
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

df <- fread(IN_CSV, na.strings = c("", "NA"), colClasses = text_col_classes(IN_CSV))
setDT(df)
if (nrow(df) == 0L) stop("The train/val pool is empty: check step 11")

stopifnot("row_id" %in% names(df))
stopifnot(uniqueN(df$row_id) == nrow(df))

# use is_homo01 if it already exists; otherwise derive it
if (!("is_homo01" %in% names(df))) {
  stop("is_homo01 is missing from the input.")
}
df[, is_homo01 := to_homo01(is_homo01)]

# Splits: one scheme per level of splits.schemes
tag_of <- function(x) sprintf("c%02d", as.integer(round(as.numeric(x) * 100)))
schemes <- vapply(unlist(cfg$splits$schemes), tag_of, character(1))
if (length(schemes) == 0L) stop("splits.schemes is empty")
for (i in seq_along(schemes)) {
  col <- paste0("pair_cluster_", sub("^c", "", schemes[i]))
  if (!col %in% names(df)) stop(sprintf("%s missing: add %s to clustering.thresholds", col, unlist(cfg$splits$schemes)[i]))
  make_split_by_cluster_stratified(df, cluster_col = col, out_col = paste0("split_", schemes[i]),
                                   strat_col = "is_homo01", val_frac = VAL_FRAC, seed = SEED_MASTER + i - 1L)
}

# Family weight on the training side of the exported scheme
SPLIT_FOR_WEIGHT <- cfg$splits$validation_split
if (!SPLIT_FOR_WEIGHT %in% names(df)) stop(sprintf("splits.validation_split (%s) is not one of the schemes", SPLIT_FOR_WEIGHT))
w_tag <- sub("^split_", "", SPLIT_FOR_WEIGHT)
WEIGHT_CLUSTER_COL <- paste0("pair_cluster_", sub("^c", "", w_tag))
size_col <- paste0("cluster_size_train_", w_tag)
weight_col <- paste0("family_weight_", w_tag)
cl_sizes <- df[get(SPLIT_FOR_WEIGHT) == "train", .N, by = c(WEIGHT_CLUSTER_COL)]
setnames(cl_sizes, "N", size_col)
df <- merge(df, cl_sizes, by = WEIGHT_CLUSTER_COL, all.x = TRUE, sort = FALSE)
df[is.na(get(size_col)), (size_col) := 1L]
df[, (weight_col) := 1.0 / get(size_col)^as.numeric(cfg$splits$weight_exponent)]

# Save
# The aliases of step 10 (pair_cluster_XX, A_cluster_XX, B_cluster_XX) and the columns that repeat
# another one (is_homo01 and is_homo_80 = is_homo; coverage_pct = 100 x coverage; has_mse_sec =
# has_special, since step 08 drops every chain with a residue that is not allowed) stay internal (the audits
# below use them); the exported tables keep the cluster columns of step 09 (pair_cluster_cXX).
ALIAS_COLS <- grep("^(pair_cluster_[0-9]+|[AB]_cluster_[0-9]+|is_homo01|is_homo_80|mol[AB]_coverage_pct|(mol[AB]_)?has_mse_sec)$", names(df), value = TRUE)
fwrite(df[, !ALIAS_COLS, with = FALSE], OUT_CSV)
cat("OK ->", OUT_CSV, "\n")

# Validation set: rows marked "val" in the chosen split
df_val <- df[get(SPLIT_FOR_WEIGHT) == "val"]
df_val[, `:=`(test_bucket = "VAL", test_subbucket = "VAL")]
fwrite(df_val[, !ALIAS_COLS, with = FALSE], ensure_parent(cfg$paths$validation))
cat("OK ->", cfg$paths$validation, "(", nrow(df_val), "rows )\n")

for (sc in schemes) {
  cat(sprintf("\n=== AUDIT split_%s ===\n", sc))
  a <- audit_split(df, paste0("split_", sc), paste0("pair_cluster_", sub("^c", "", sc)))
  print(a$counts)
  cat("homo global:", a$homo_frac_global, "\n")
  cat("homo val:", a$homo_frac_val, "\n")
  cat("val frac real:", a$val_frac_real, "\n")
  cat("unique clusters train:", a$unique_clusters_train, "\n")
  cat("unique clusters val:", a$unique_clusters_val, "\n")
  print(a$quality)
}

cat(sprintf("\nFamily weight (%s, exponent %s):\n", w_tag, cfg$splits$weight_exponent))
print(summary(df[[weight_col]]))
print(summary(df[[size_col]]))


# =============================================================================
# EXTRA AUDITS
# A) cluster overlap between train and val
# B) cluster size distribution in val
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

# ---- for each scheme
pc_of <- function(sc) paste0("pair_cluster_", sub("^c", "", sc))
overlap_audits <- rbindlist(lapply(schemes, function(sc) audit_cluster_overlap(df, paste0("split_", sc), pc_of(sc))),
                            use.names = TRUE, fill = TRUE)
val_cluster_size_audits <- rbindlist(lapply(schemes, function(sc) audit_val_cluster_sizes(df, paste0("split_", sc), pc_of(sc))),
                                     use.names = TRUE, fill = TRUE)

cat("\n=== AUDIT A: CLUSTER OVERLAP TRAIN/VAL ===\n")
print(overlap_audits)

cat("\n=== AUDIT B: VAL CLUSTER SIZE DISTRIBUTION ===\n")
print(val_cluster_size_audits)

# optional: save
#fwrite(overlap_audits, "audit_cluster_overlap_train_val.csv")
#fwrite(val_cluster_size_audits, "audit_val_cluster_size_distribution.csv")

# Show the exact overlapping clusters
show_overlap_clusters <- function(DT, split_col, cluster_col) {
  tr <- unique(DT[get(split_col) == "train" & !is.na(get(cluster_col)), get(cluster_col)])
  va <- unique(DT[get(split_col) == "val"   & !is.na(get(cluster_col)), get(cluster_col)])
  inter <- intersect(tr, va)
  if (length(inter) > 0) print(head(data.table(cluster_id = inter), 50))
}

for (sc in schemes) show_overlap_clusters(df, paste0("split_", sc), pc_of(sc))
