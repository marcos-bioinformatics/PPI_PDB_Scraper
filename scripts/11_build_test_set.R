suppressPackageStartupMessages({
  library(data.table)
})

source(file.path("scripts", "common.R"))
cfg <- load_config()
ts <- cfg$test_set

set.seed(ts$seed)

# Test split protocol (C1 / C2 / C3)
#
# C1 (easy / train-like):
#   Rows sampled from the remaining pool after defining C2/C3, stratified by
#   homomer/heteromer status and coarse length/interface-area bins.
#   C1 is disjoint from train/val at the row level, but does not enforce strong
#   novelty at the chain- or pair-homology level. It is intended to measure
#   in-distribution / near-train performance.
#
# C2 (new pairs, known proteins):
#   Defined at the pair-homology level using pair MMseqs clusters
#   (here: pair_cluster_70 and pair_cluster_80), with additional blocking of
#   exact sequence-identical pairs via pair_seq_exact.
#   Thus, C2 evaluates generalization to unseen pairings of proteins that still
#   belong to the known protein universe of the dataset.
#
# C3 (unseen proteins):
#   Defined at the chain-homology level using single-chain MMseqs clusters
#   (here: chain clusters at c40 and c50, i.e. A_cluster_40/B_cluster_40 and
#   A_cluster_50/B_cluster_50).
#   For each C3 sub-bucket, all rows containing any selected held-out chain
#   cluster are removed from the train/val pool.
#   Thus, C3 evaluates generalization to proteins not seen during training at
#   the corresponding chain-homology threshold.
#
# In short:
#   - C1: train-like / in-distribution
#   - C2: novelty at the pair-homology level
#   - C3: novelty at the chain-homology level

# =============================================================================
# CONFIG
# =============================================================================
INPUT_CSV <- cfg$paths$pruned

# ---------------------------------------------------------------------------
# Tamaños del test (config test_set; por defecto 700 filas en total)
#   C3 = 100 (c40) + 100 (c50) = 200
#   C2 = 100 (c70) + 100 (c80) = 200
#   C1 = 300
# ---------------------------------------------------------------------------
N_C3_40 <- as.integer(ts$n_c3_40)
N_C3_50 <- as.integer(ts$n_c3_50)
N_C1    <- as.integer(ts$n_c1)

C2_PLAN <- data.table(
  subbucket = c("c70", "c80"),
  pair_col  = c("pair_cluster_70", "pair_cluster_80"),
  n_pairs   = c(as.integer(ts$n_c2_c70), as.integer(ts$n_c2_c80))
)

NBINS_LEN  <- as.integer(ts$nbins_len)
NBINS_AREA <- as.integer(ts$nbins_area)
HOMO_FRAC  <- ts$homo_frac
C3_MULT    <- as.integer(ts$c3_mult)

# C2 constraints
C2_CHAIN_CAP_98_GLOBAL <- as.integer(ts$c2_chain_cap_98_global)
C2_MIN_REMAIN_CHAIN    <- as.integer(ts$c2_min_remain_chain)

# C2: qué fila representa a cada pair cluster y en qué orden se eligen
#   "random" (por defecto): al azar, calidad similar a la de train
#   "best": mejor calidad primero; "worst": peor calidad primero (criterio de Seq8)
C2_PICK <- if (is.null(ts$c2_pick)) "random" else ts$c2_pick
stopifnot(C2_PICK %in% c("random", "best", "worst"))

# =============================================================================
# INPUT
# =============================================================================
df <- as.data.table(read.csv(INPUT_CSV, stringsAsFactors = FALSE))

# =============================================================================
# REQUIRED COLUMNS
# =============================================================================
required_cols <- c(
  "row_id",
  "pair_seq_exact",
  "is_homo01",
  "molA_coverage", "molB_coverage",
  "molA_coverage_pct", "molB_coverage_pct",
  "p_value", "r_free", "resolution", "interface_area",
  "molA_seq_len", "molB_seq_len",
  "pair_cluster_60", "pair_cluster_70", "pair_cluster_80",
  "pair_cluster_98", "pair_cluster_95", "pair_cluster_90",
  "A_cluster_98", "B_cluster_98",
  "A_cluster_50", "B_cluster_50",
  "A_cluster_40", "B_cluster_40"
)

missing <- setdiff(required_cols, names(df))
if (length(missing) > 0L) {
  stop("Faltan columnas requeridas: ", paste(missing, collapse = ", "))
}

# =============================================================================
# HELPERS
# =============================================================================
as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

coerce_is_homo01 <- function(x) {
  x_chr <- tolower(trimws(as.character(x)))
  out <- rep(NA_integer_, length(x_chr))
  out[x_chr %chin% c("1","true","t","yes","y")] <- 1L
  out[x_chr %chin% c("0","false","f","no","n")] <- 0L
  out
}

ntile_vec <- function(x, n) {
  r <- data.table::frank(x, ties.method = "average", na.last = "keep")
  out <- integer(length(x))
  nn <- sum(!is.na(r))
  if (nn > 0L) out[!is.na(r)] <- pmin(n, pmax(1L, ceiling(r[!is.na(r)] * n / nn)))
  out[is.na(r)] <- n + 1L
  out
}

add_bins_dt <- function(DT, nbins_len = 4L, nbins_area = 4L) {
  DT[, .len_max := pmax(as_num(molA_seq_len), as_num(molB_seq_len), na.rm = TRUE)]
  DT[is.infinite(.len_max), .len_max := NA_real_]
  DT[, .area_num_tmp := as_num(interface_area)]
  DT[, len_bin  := as.integer(ntile_vec(.len_max, nbins_len))]
  DT[, area_bin := as.integer(ntile_vec(.area_num_tmp, nbins_area))]
  DT[, c(".len_max",".area_num_tmp") := NULL]
  invisible(DT)
}

add_quality_cols_dt <- function(DT) {
  cov_min_frac <- pmin(as_num(DT$molA_coverage), as_num(DT$molB_coverage), na.rm = TRUE)
  cov_min_frac[is.infinite(cov_min_frac)] <- NA_real_
  
  cov_min_pct <- pmin(as_num(DT$molA_coverage_pct), as_num(DT$molB_coverage_pct), na.rm = TRUE)
  cov_min_pct[is.infinite(cov_min_pct)] <- NA_real_
  
  DT[, .cov_min_frac := cov_min_frac]
  DT[, .cov_min_pct  := cov_min_pct]
  DT[, .p_value_num  := as_num(p_value)]
  DT[, .r_free_num   := as_num(r_free)]
  DT[, .res_num      := as_num(resolution)]
  DT[, .area_num     := as_num(interface_area)]
  invisible(DT)
}

drop_quality_cols_dt <- function(DT) {
  drop_cols <- intersect(
    c(".cov_min_frac",".cov_min_pct",".p_value_num",".r_free_num",".res_num",".area_num"),
    names(DT)
  )
  if (length(drop_cols) > 0L) DT[, (drop_cols) := NULL]
  invisible(DT)
}

strat_sample_dt <- function(DT, n, strata_cols = c("len_bin","area_bin")) {
  if (n <= 0L || nrow(DT) == 0L) return(DT[0])
  
  df2 <- copy(DT)
  df2[, .stratum := do.call(paste, c(.SD, sep = "|")), .SDcols = strata_cols]
  df2 <- df2[!is.na(.stratum)]
  if (nrow(df2) == 0L) return(DT[0])
  
  counts <- df2[, .(avail = .N), by = .stratum]
  total_avail <- sum(counts$avail)
  
  counts[, target := floor(n * avail / total_avail)]
  rem <- n - sum(counts$target)
  
  if (rem > 0L) {
    counts[, slack := avail - target]
    setorder(counts, -slack)
    idx <- seq_len(min(rem, nrow(counts)))
    counts[idx, target := target + 1L]
    counts[, slack := NULL]
  }
  
  # join corregido
  setnames(counts, ".stratum", "stratum_key")
  setnames(df2, ".stratum", "stratum_key")
  df2 <- counts[df2, on = "stratum_key"]
  setnames(df2, "target", ".target")
  
  out <- df2[, {
    k <- unique(.target)
    k <- k[!is.na(k)]
    if (length(k) == 0L || k[1] <= 0L) .SD[0] else .SD[sample.int(.N, min(.N, k[1]))]
  }, by = stratum_key]
  
  if (".target" %in% names(out)) out[, .target := NULL]
  
  if (nrow(out) < n) {
    need <- n - nrow(out)
    rest <- df2[!out, on = "row_id"]
    if (nrow(rest) > 0L) {
      out <- rbindlist(list(out, rest[sample.int(.N, min(.N, need))]), use.names = TRUE, fill = TRUE)
    }
  }
  
  if ("stratum_key" %in% names(out)) out[, stratum_key := NULL]
  unique(out, by = "row_id")
}

sample_balanced_dt <- function(DT, n, homo_frac = 0.5, strata_cols = c("len_bin","area_bin")) {
  if (n <= 0L || nrow(DT) == 0L) return(DT[0])
  
  n_homo  <- floor(n * homo_frac)
  n_heter <- n - n_homo
  
  df_h <- DT[is_homo01 == 1L]
  df_t <- DT[is_homo01 == 0L]
  
  if (nrow(df_h) < n_homo) {
    n_homo <- nrow(df_h)
    n_heter <- n - n_homo
  }
  if (nrow(df_t) < n_heter) {
    n_heter <- nrow(df_t)
    n_homo <- min(n - n_heter, nrow(df_h))
  }
  
  out <- rbindlist(list(
    strat_sample_dt(df_h, n_homo, strata_cols),
    strat_sample_dt(df_t, n_heter, strata_cols)
  ), use.names = TRUE, fill = TRUE)
  
  if (nrow(out) < n) {
    need <- n - nrow(out)
    rest <- DT[!out, on = "row_id"]
    if (nrow(rest) > 0L) {
      out <- rbindlist(list(out, rest[sample.int(.N, min(.N, need))]), use.names = TRUE, fill = TRUE)
    }
  }
  
  unique(out, by = "row_id")
}

update_chain_counts <- function(chain_counts, rows) {
  ch <- unique(c(rows$A_cluster_98, rows$B_cluster_98))
  ch <- ch[!is.na(ch)]
  if (length(ch) == 0L) return(chain_counts)
  
  tab <- table(ch)
  inc <- data.table(chain98 = names(tab), inc = as.integer(tab))
  
  if (nrow(chain_counts) == 0L) {
    return(inc[, .(chain98, n = inc)])
  }
  
  out <- merge(chain_counts, inc, by = "chain98", all = TRUE)
  out[is.na(n), n := 0L]
  out[is.na(inc), inc := 0L]
  out[, n := n + inc]
  out[, inc := NULL]
  out
}

chains_over_cap <- function(chain_counts, cap) {
  if (nrow(chain_counts) == 0L) return(character(0))
  chain_counts[n >= cap, chain98]
}

rank_for_test_harder <- function(DT) {
  x <- copy(DT)
  add_quality_cols_dt(x)
  # Para C2/C3: sacrificamos menor coverage / peor calidad al test
  setorder(
    x,
    .cov_min_pct,
    .cov_min_frac,
    -.p_value_num,
    -.r_free_num,
    -.res_num,
    .area_num
  )
  drop_quality_cols_dt(x)
  x
}

rank_for_test_best <- function(DT) {
  x <- copy(DT)
  add_quality_cols_dt(x)
  # Mejor calidad primero
  setorder(
    x,
    -.cov_min_pct,
    -.cov_min_frac,
    .p_value_num,
    .r_free_num,
    .res_num,
    -.area_num
  )
  drop_quality_cols_dt(x)
  x
}

rank_for_c2 <- function(DT) {
  switch(C2_PICK,
         worst  = rank_for_test_harder(DT),
         best   = rank_for_test_best(DT),
         random = DT[sample(.N)])
}

rank_for_test_easy <- function(DT) {
  x <- copy(DT)
  add_quality_cols_dt(x)
  # Para C1: más train-like
  setorder(
    x,
    -easy_score,
    -.cov_min_pct,
    -.cov_min_frac,
    .p_value_num,
    .r_free_num,
    .res_num,
    -.area_num
  )
  drop_quality_cols_dt(x)
  x
}

# =============================================================================
# BASE DF0
# =============================================================================
df0 <- copy(df)

df0[, row_id := as.character(row_id)]
stopifnot(!any(is.na(df0$row_id)))
stopifnot(uniqueN(df0$row_id) == nrow(df0))

df0[, is_homo01 := coerce_is_homo01(is_homo01)]
add_bins_dt(df0, NBINS_LEN, NBINS_AREA)

id_cols <- c(
  "pair_seq_exact",
  "pair_cluster_60", "pair_cluster_70", "pair_cluster_80",
  "pair_cluster_90", "pair_cluster_95", "pair_cluster_98",
  "A_cluster_98", "B_cluster_98",
  "A_cluster_50", "B_cluster_50",
  "A_cluster_40", "B_cluster_40"
)

df0[, (id_cols) := lapply(.SD, as.character), .SDcols = id_cols]

setindex(df0, row_id)
setindex(df0, pair_seq_exact)

# =============================================================================
# C3: unseen proteins (c40 y c50)
# =============================================================================
build_C3_bucket <- function(df_pool, N, chain_col_A, chain_col_B, subbucket_label) {
  pool <- df_pool[!is.na(get(chain_col_A)) & !is.na(get(chain_col_B))]
  if (nrow(pool) == 0L) {
    return(list(
      selected = pool[0],
      blocked  = data.table(row_id = character(), blocked_reason = character())
    ))
  }
  
  freq <- table(c(pool[[chain_col_A]], pool[[chain_col_B]]))
  pool[, .freqA := as.integer(freq[as.character(get(chain_col_A))])]
  pool[, .freqB := as.integer(freq[as.character(get(chain_col_B))])]
  pool[is.na(.freqA), .freqA := 0L]
  pool[is.na(.freqB), .freqB := 0L]
  pool[, .freq_sum := .freqA + .freqB]
  
  setorder(pool, .freq_sum)
  pool_small <- pool[1:min(.N, N * C3_MULT)]
  pool_small <- rank_for_test_harder(pool_small)
  
  sel <- sample_balanced_dt(pool_small, N, HOMO_FRAC, strata_cols = c("len_bin","area_bin"))
  sel[, `:=`(test_bucket = "C3", test_subbucket = subbucket_label)]
  
  holdout <- unique(c(sel[[chain_col_A]], sel[[chain_col_B]]))
  holdout <- holdout[!is.na(holdout)]
  
  blocked <- df_pool[
    get(chain_col_A) %chin% holdout | get(chain_col_B) %chin% holdout,
    .(row_id, blocked_reason = paste0("blocked_by_C3_", subbucket_label))
  ]
  
  drop_cols <- intersect(c(".freqA",".freqB",".freq_sum"), names(sel))
  if (length(drop_cols) > 0L) sel[, (drop_cols) := NULL]
  
  list(
    selected = unique(sel, by = "row_id"),
    blocked  = unique(blocked, by = "row_id")
  )
}

c3_40 <- build_C3_bucket(df0, N_C3_40, "A_cluster_40", "B_cluster_40", "unseen_chain_c40")
df_pool_after_C3_40 <- df0[!c3_40$blocked, on = "row_id"]

c3_50 <- build_C3_bucket(df_pool_after_C3_40, N_C3_50, "A_cluster_50", "B_cluster_50", "unseen_chain_c50")
df_pool_after_C3 <- df_pool_after_C3_40[!c3_50$blocked, on = "row_id"]

c3_selected <- rbindlist(list(c3_40$selected, c3_50$selected), use.names = TRUE, fill = TRUE)
blocked_c3  <- unique(rbindlist(list(c3_40$blocked, c3_50$blocked), use.names = TRUE, fill = TRUE), by = "row_id")

# =============================================================================
# C2: reps por threshold + bloqueo por pair_cluster_t y pair_seq_exact
# =============================================================================
build_reps_for_threshold_dt <- function(df_pool_in, pair_col, subbucket_label, chain_counts_global) {
  stopifnot(pair_col %in% names(df_pool_in))
  
  tmpA <- df_pool_in[!is.na(get(pair_col)) & !is.na(A_cluster_98), .(chain98 = A_cluster_98, pc = get(pair_col))]
  tmpB <- df_pool_in[!is.na(get(pair_col)) & !is.na(B_cluster_98), .(chain98 = B_cluster_98, pc = get(pair_col))]
  tmp  <- unique(rbindlist(list(tmpA, tmpB), use.names = TRUE))
  chain_partner_counts <- tmp[, .(n_distinct_pairs = .N), by = chain98]
  
  df_elig <- df_pool_in[!is.na(get(pair_col))]
  
  df_elig <- merge(df_elig, chain_partner_counts, by.x = "A_cluster_98", by.y = "chain98", all.x = TRUE)
  setnames(df_elig, "n_distinct_pairs", ".A_npairs")
  df_elig <- merge(df_elig, chain_partner_counts, by.x = "B_cluster_98", by.y = "chain98", all.x = TRUE)
  setnames(df_elig, "n_distinct_pairs", ".B_npairs")
  df_elig[is.na(.A_npairs), .A_npairs := 0L]
  df_elig[is.na(.B_npairs), .B_npairs := 0L]
  
  df_elig <- df_elig[.A_npairs >= 2L & .B_npairs >= 2L]
  
  over <- chains_over_cap(chain_counts_global, C2_CHAIN_CAP_98_GLOBAL)
  if (length(over) > 0L) {
    df_elig <- df_elig[!(A_cluster_98 %chin% over | B_cluster_98 %chin% over)]
  }
  
  reps <- rank_for_c2(df_elig)
  reps <- reps[, .SD[1], by = pair_col]
  reps[, `:=`(test_bucket = "C2", test_subbucket = subbucket_label)]
  
  reps
}

select_C2_reps_dt <- function(reps_all, n_pairs, chain_counts_global) {
  if (nrow(reps_all) == 0L) return(reps_all[0])
  
  cand <- rank_for_c2(reps_all)
  
  want_homo  <- floor(n_pairs * HOMO_FRAC)
  want_heter <- n_pairs - want_homo
  
  cc <- copy(chain_counts_global)
  
  pick_group <- function(is_homo_val, k) {
    if (k <= 0L) return(cand[0])
    sub <- cand[is_homo01 == is_homo_val]
    if (nrow(sub) == 0L) return(sub[0])
    
    out <- sub[0]
    for (i in seq_len(nrow(sub))) {
      r <- sub[i]
      tmp_cc <- copy(cc)
      tmp_cc <- update_chain_counts(tmp_cc, r)
      if (any(tmp_cc$n > C2_CHAIN_CAP_98_GLOBAL)) next
      out <- rbindlist(list(out, r), use.names = TRUE, fill = TRUE)
      cc  <- tmp_cc
      if (nrow(out) >= k) break
    }
    out
  }
  
  out_h <- pick_group(1L, want_homo)
  out_t <- pick_group(0L, want_heter)
  reps_sel <- rbindlist(list(out_h, out_t), use.names = TRUE, fill = TRUE)
  
  if (nrow(reps_sel) < n_pairs) {
    need <- n_pairs - nrow(reps_sel)
    rest <- cand[!reps_sel, on = "row_id"]
    add <- rest[0]
    for (i in seq_len(nrow(rest))) {
      r <- rest[i]
      tmp_cc <- copy(cc)
      tmp_cc <- update_chain_counts(tmp_cc, r)
      if (any(tmp_cc$n > C2_CHAIN_CAP_98_GLOBAL)) next
      add <- rbindlist(list(add, r), use.names = TRUE, fill = TRUE)
      cc  <- tmp_cc
      if (nrow(add) >= need) break
    }
    reps_sel <- rbindlist(list(reps_sel, add), use.names = TRUE, fill = TRUE)
  }
  
  unique(reps_sel, by = "row_id")
}

make_C2_blocked_rows <- function(df_pool_in, pair_col, subbucket_label, reps_sel) {
  chosen_clusters <- unique(reps_sel[[pair_col]])
  chosen_exact    <- unique(reps_sel$pair_seq_exact)
  
  df_pool_in[
    (get(pair_col) %chin% chosen_clusters) | (pair_seq_exact %chin% chosen_exact),
    .(row_id, blocked_reason = paste0("blocked_by_C2_", subbucket_label, "_", pair_col, "_plus_exact"))
  ]
}

C2_sel_list <- list()
c2_diag <- list()

df_pool_work <- copy(df_pool_after_C3)
chain_counts_global <- data.table(chain98 = character(), n = integer())

for (i in seq_len(nrow(C2_PLAN))) {
  subb    <- C2_PLAN$subbucket[i]
  paircol <- C2_PLAN$pair_col[i]
  npairs  <- C2_PLAN$n_pairs[i]
  
  reps_all <- build_reps_for_threshold_dt(df_pool_work, paircol, subb, chain_counts_global)
  reps_sel <- select_C2_reps_dt(reps_all, npairs, chain_counts_global)
  
  chain_counts_global <- update_chain_counts(chain_counts_global, reps_sel)
  
  C2_sel_list[[subb]] <- reps_sel
  
  c2_diag[[subb]] <- data.table(
    subbucket  = subb,
    pair_col   = paircol,
    target_n   = npairs,
    selected_n = nrow(reps_sel),
    candidates = nrow(reps_all)
  )
  
  blocked_i <- make_C2_blocked_rows(df_pool_work, paircol, subb, reps_sel)
  df_pool_work <- df_pool_work[!blocked_i, on = "row_id"]
}

c2_diag_tbl <- rbindlist(c2_diag, use.names = TRUE, fill = TRUE)

blocked_c2 <- unique(
  rbindlist(lapply(seq_len(nrow(C2_PLAN)), function(i) {
    subb    <- C2_PLAN$subbucket[i]
    paircol <- C2_PLAN$pair_col[i]
    make_C2_blocked_rows(df_pool_after_C3, paircol, subb, C2_sel_list[[subb]])
  }), use.names = TRUE, fill = TRUE),
  by = "row_id"
)

df_pool_after_C2 <- df_pool_after_C3[!blocked_c2, on = "row_id"]
c2_selected <- rbindlist(C2_sel_list, use.names = TRUE, fill = TRUE)

# =============================================================================
# C1: easy / train-like
# =============================================================================
pc80_freq <- df_pool_after_C2[!is.na(pair_cluster_80), .(.pc80_n = .N), by = pair_cluster_80]

tmp98 <- rbindlist(list(
  df_pool_after_C2[!is.na(A_cluster_98), .(chain98 = A_cluster_98)],
  df_pool_after_C2[!is.na(B_cluster_98), .(chain98 = B_cluster_98)]
), use.names = TRUE)

chain98_freq <- tmp98[, .(chain98_n = .N), by = chain98]

df_c1_pool <- copy(df_pool_after_C2)
df_c1_pool[pc80_freq, on = "pair_cluster_80", .pc80_n := i..pc80_n]
df_c1_pool[is.na(.pc80_n), .pc80_n := 1L]
df_c1_pool[chain98_freq, on = .(A_cluster_98 = chain98), .A98_n := i.chain98_n]
df_c1_pool[chain98_freq, on = .(B_cluster_98 = chain98), .B98_n := i.chain98_n]
df_c1_pool[is.na(.A98_n), .A98_n := 1L]
df_c1_pool[is.na(.B98_n), .B98_n := 1L]

df_c1_pool[, easy_score := .pc80_n + .A98_n + .B98_n]
df_c1_pool <- rank_for_test_easy(df_c1_pool)

df_c1_pool_pref <- df_c1_pool[1:min(.N, N_C1 * 20L)]
if (nrow(df_c1_pool_pref) < N_C1) df_c1_pool_pref <- df_c1_pool

c1_selected <- sample_balanced_dt(df_c1_pool_pref, N_C1, HOMO_FRAC, strata_cols = c("len_bin","area_bin"))
c1_selected[, `:=`(test_bucket = "C1", test_subbucket = "easy")]
blocked_c1 <- c1_selected[, .(row_id, blocked_reason = "blocked_by_C1_test")]

# =============================================================================
# FINAL: test / blocked / trainval
# =============================================================================
df_test <- rbindlist(list(c3_selected, c2_selected, c1_selected), use.names = TRUE, fill = TRUE)
df_test <- unique(df_test, by = "row_id")

cols <- names(df_test)
if (all(c("row_id","test_bucket","test_subbucket") %in% cols)) {
  rest <- setdiff(cols, c("row_id","test_bucket","test_subbucket"))
  setcolorder(df_test, c("row_id","test_bucket","test_subbucket", rest))
}

df_blocked <- unique(
  rbindlist(list(blocked_c3, blocked_c2, blocked_c1), use.names = TRUE, fill = TRUE),
  by = "row_id"
)

df_trainval_pool <- df0[!df_blocked, on = "row_id"]

# =============================================================================
# AUDITS
# =============================================================================
summary_tbl <- data.table(
  n_total_rows         = nrow(df0),
  requested_total_rows = N_C3_40 + N_C3_50 + sum(C2_PLAN$n_pairs) + N_C1,
  n_test_rows          = nrow(df_test),
  n_blocked_rows       = nrow(df_blocked),
  n_trainval_rows      = nrow(df_trainval_pool),
  n_blocked_not_test   = nrow(df_blocked) - nrow(df_test),
  blocked_per_test     = round(nrow(df_blocked) / max(1L, nrow(df_test)), 3)
)

bucket_counts <- df_test[, .(n_rows = .N), by = .(test_bucket, test_subbucket)]
setorder(bucket_counts, test_bucket, test_subbucket)

df_test[, group := fifelse(is_homo01 == 1L, "homo", "hetero")]
bucket_balance <- df_test[, .(n = .N), by = .(test_bucket, test_subbucket, group)]
bucket_balance[, frac := n / sum(n), by = .(test_bucket, test_subbucket)]
df_test[, group := NULL]

bucket_block_cost <- rbindlist(list(
  data.table(
    test_bucket = "C3",
    test_subbucket = "unseen_chain_c40",
    test_n = nrow(c3_40$selected),
    blocked_n = nrow(c3_40$blocked),
    blocked_per_test = round(nrow(c3_40$blocked) / max(1L, nrow(c3_40$selected)), 3)
  ),
  data.table(
    test_bucket = "C3",
    test_subbucket = "unseen_chain_c50",
    test_n = nrow(c3_50$selected),
    blocked_n = nrow(c3_50$blocked),
    blocked_per_test = round(nrow(c3_50$blocked) / max(1L, nrow(c3_50$selected)), 3)
  ),
  rbindlist(lapply(seq_len(nrow(C2_PLAN)), function(i) {
    subb <- C2_PLAN$subbucket[i]
    paircol <- C2_PLAN$pair_col[i]
    sel <- C2_sel_list[[subb]]
    blk <- make_C2_blocked_rows(df_pool_after_C3, paircol, subb, sel)
    data.table(
      test_bucket = "C2",
      test_subbucket = subb,
      test_n = nrow(sel),
      blocked_n = nrow(blk),
      blocked_per_test = round(nrow(blk) / max(1L, nrow(sel)), 3)
    )
  }), use.names = TRUE, fill = TRUE),
  data.table(
    test_bucket = "C1",
    test_subbucket = "easy",
    test_n = nrow(c1_selected),
    blocked_n = nrow(blocked_c1),
    blocked_per_test = round(nrow(blocked_c1) / max(1L, nrow(c1_selected)), 3)
  )
), use.names = TRUE, fill = TRUE)

coverage_tbl <- df_test[, .(
  n = .N,
  cov_min_pct_median = median(pmin(as_num(molA_coverage_pct), as_num(molB_coverage_pct), na.rm = TRUE), na.rm = TRUE),
  cov_min_pct_p10    = quantile(pmin(as_num(molA_coverage_pct), as_num(molB_coverage_pct), na.rm = TRUE), 0.10, na.rm = TRUE),
  p_value_median     = median(as_num(p_value), na.rm = TRUE),
  r_free_median      = median(as_num(r_free), na.rm = TRUE),
  resolution_median  = median(as_num(resolution), na.rm = TRUE),
  area_median        = median(as_num(interface_area), na.rm = TRUE)
), by = .(test_bucket, test_subbucket)]

trainval_cluster_loss <- data.table(
  metric = c(
    "pair_cluster_98_lost",
    "pair_cluster_95_lost",
    "pair_cluster_90_lost",
    "pair_cluster_80_lost",
    "pair_cluster_70_lost",
    "pair_cluster_60_lost",
    "chain_cluster_98_lost",
    "chain_cluster_50_lost",
    "chain_cluster_40_lost"
  ),
  n_lost = c(
    uniqueN(df0$pair_cluster_98[!is.na(df0$pair_cluster_98)]) - uniqueN(df_trainval_pool$pair_cluster_98[!is.na(df_trainval_pool$pair_cluster_98)]),
    uniqueN(df0$pair_cluster_95[!is.na(df0$pair_cluster_95)]) - uniqueN(df_trainval_pool$pair_cluster_95[!is.na(df_trainval_pool$pair_cluster_95)]),
    uniqueN(df0$pair_cluster_90[!is.na(df0$pair_cluster_90)]) - uniqueN(df_trainval_pool$pair_cluster_90[!is.na(df_trainval_pool$pair_cluster_90)]),
    uniqueN(df0$pair_cluster_80[!is.na(df0$pair_cluster_80)]) - uniqueN(df_trainval_pool$pair_cluster_80[!is.na(df_trainval_pool$pair_cluster_80)]),
    uniqueN(df0$pair_cluster_70[!is.na(df0$pair_cluster_70)]) - uniqueN(df_trainval_pool$pair_cluster_70[!is.na(df_trainval_pool$pair_cluster_70)]),
    uniqueN(df0$pair_cluster_60[!is.na(df0$pair_cluster_60)]) - uniqueN(df_trainval_pool$pair_cluster_60[!is.na(df_trainval_pool$pair_cluster_60)]),
    uniqueN(c(df0$A_cluster_98[!is.na(df0$A_cluster_98)], df0$B_cluster_98[!is.na(df0$B_cluster_98)])) -
      uniqueN(c(df_trainval_pool$A_cluster_98[!is.na(df_trainval_pool$A_cluster_98)], df_trainval_pool$B_cluster_98[!is.na(df_trainval_pool$B_cluster_98)])),
    uniqueN(c(df0$A_cluster_50[!is.na(df0$A_cluster_50)], df0$B_cluster_50[!is.na(df0$B_cluster_50)])) -
      uniqueN(c(df_trainval_pool$A_cluster_50[!is.na(df_trainval_pool$A_cluster_50)], df_trainval_pool$B_cluster_50[!is.na(df_trainval_pool$B_cluster_50)])),
    uniqueN(c(df0$A_cluster_40[!is.na(df0$A_cluster_40)], df0$B_cluster_40[!is.na(df0$B_cluster_40)])) -
      uniqueN(c(df_trainval_pool$A_cluster_40[!is.na(df_trainval_pool$A_cluster_40)], df_trainval_pool$B_cluster_40[!is.na(df_trainval_pool$B_cluster_40)]))
  )
)

# Leak checks
c3_40_ch <- unique(c(c3_40$selected$A_cluster_40, c3_40$selected$B_cluster_40))
c3_40_ch <- c3_40_ch[!is.na(c3_40_ch)]
train_40 <- unique(c(df_trainval_pool$A_cluster_40, df_trainval_pool$B_cluster_40))
train_40 <- train_40[!is.na(train_40)]
n_leak_C3_40 <- sum(c3_40_ch %chin% train_40)

c3_50_ch <- unique(c(c3_50$selected$A_cluster_50, c3_50$selected$B_cluster_50))
c3_50_ch <- c3_50_ch[!is.na(c3_50_ch)]
train_50 <- unique(c(df_trainval_pool$A_cluster_50, df_trainval_pool$B_cluster_50))
train_50 <- train_50[!is.na(train_50)]
n_leak_C3_50 <- sum(c3_50_ch %chin% train_50)

df_C2 <- df_test[test_bucket == "C2"]
n_leak_exact_C2 <- length(intersect(df_C2$pair_seq_exact, df_trainval_pool$pair_seq_exact))

leak_cluster_C2 <- rbindlist(lapply(seq_len(nrow(C2_PLAN)), function(i) {
  subb <- C2_PLAN$subbucket[i]
  paircol <- C2_PLAN$pair_col[i]
  sel <- C2_sel_list[[subb]]
  if (is.null(sel) || nrow(sel) == 0L) {
    return(data.table(subbucket = subb, pair_col = paircol, n_leak_clusters = NA_integer_))
  }
  chosen <- unique(sel[[paircol]])
  left_in_train <- unique(df_trainval_pool[[paircol]])
  data.table(subbucket = subb, pair_col = paircol, n_leak_clusters = sum(chosen %chin% left_in_train, na.rm = TRUE))
}), use.names = TRUE, fill = TRUE)

cat("\n==================== SUMMARY ====================\n")
print(summary_tbl)

cat("\n==================== BUCKET COUNTS ====================\n")
print(bucket_counts)

cat("\n==================== BUCKET BALANCE ====================\n")
print(bucket_balance)

cat("\n==================== BUCKET BLOCK COST ====================\n")
print(bucket_block_cost)

cat("\n==================== C2 DIAGNOSTICS ====================\n")
print(c2_diag_tbl)

cat("\n==================== COVERAGE / QUALITY BY BUCKET ====================\n")
print(coverage_tbl)

cat("\n==================== TRAINVAL CLUSTER LOSS ====================\n")
print(trainval_cluster_loss)

cat("\n==================== SANITY CHECKS ====================\n")
cat("Leak C3 chain_40:", n_leak_C3_40, "\n")
cat("Leak C3 chain_50:", n_leak_C3_50, "\n")
cat("Leak C2 pair_seq_exact:", n_leak_exact_C2, "\n")
print(leak_cluster_C2)

if (n_leak_C3_40 > 0) cat("WARNING: leak C3_40 -> trainval_pool\n")
if (n_leak_C3_50 > 0) cat("WARNING: leak C3_50 -> trainval_pool\n")
if (n_leak_exact_C2 > 0) cat("WARNING: leak exact C2 -> trainval_pool\n")
if (any(leak_cluster_C2$n_leak_clusters > 0, na.rm = TRUE)) {
  cat("WARNING: leak pair_cluster_t C2 -> trainval_pool\n")
}

fwrite(df_test,          ensure_parent(cfg$paths$test))
fwrite(df_trainval_pool, ensure_parent(cfg$paths$trainval_pool))
