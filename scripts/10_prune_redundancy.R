suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

# ============================================================
# Pair-cluster prune v3
# - Exact dedupe opcional
# - Pruning ligero/modificable en pair_cluster_98 o pair_cluster_95
# - Sin has_special
# - Sin nres_range_thr
# - Diversidad geométrica basada en interface_area
# - Ranking por coverage_pct -> p_value -> area -> r_free -> resolution -> solvation
# - Auditorías profundas before/after
# ============================================================

source(file.path("scripts", "common.R"))
cfg <- load_config()
pr <- cfg$pruning

df <- read.csv(cfg$paths$clustered) %>%
  select(-matches("^Unnamed"))

# -----------------------------
# 0) Helpers robustos
# -----------------------------
as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

coerce_bool01 <- function(x) {
  x_chr <- trimws(tolower(as.character(x)))
  dplyr::case_when(
    x_chr %in% c("1","true","t","yes","y")  ~ 1L,
    x_chr %in% c("0","false","f","no","n")  ~ 0L,
    TRUE                                    ~ NA_integer_
  )
}

derive_is_homo01 <- function(df, candidates = c("is_homo", "is_homo_entity", "is_homo_80")) {
  present <- candidates[candidates %in% names(df)]
  if (length(present) == 0) return(rep(NA_integer_, nrow(df)))
  for (col in present) {
    v <- coerce_bool01(df[[col]])
    if (sum(!is.na(v)) > 0) return(v)
  }
  rep(NA_integer_, nrow(df))
}

mk_pair_key <- function(a, b) paste(pmin(a, b), pmax(a, b), sep = "||")

# -----------------------------
# 1) Ensure IDs + aliases
# -----------------------------
ensure_ids_and_clusters <- function(df, use_label_asym = FALSE, thrs = c(98,95,90,80,70,60,50,40)) {
  chainA <- if (use_label_asym) "molA_label_asym_id" else "molA_chain_id"
  chainB <- if (use_label_asym) "molB_label_asym_id" else "molB_chain_id"
  
  stopifnot(all(c("pdb_id","assembly_id","interface_id","pair_index_in_interface", chainA, chainB) %in% names(df)))
  
  df <- df %>%
    mutate(
      pdb_id = tolower(as.character(pdb_id)),
      assembly_id = suppressWarnings(as.integer(assembly_id)),
      interface_id = as.character(interface_id),
      pair_index_in_interface = suppressWarnings(as.integer(pair_index_in_interface))
    )
  
  if (!all(c("chain_instance_id_A","chain_instance_id_B") %in% names(df))) {
    df <- df %>%
      mutate(
        chain_instance_id_A = paste(pdb_id, assembly_id, as.character(.data[[chainA]]), sep = "|"),
        chain_instance_id_B = paste(pdb_id, assembly_id, as.character(.data[[chainB]]), sep = "|")
      )
  }
  
  if (!("chain_pair_id" %in% names(df))) {
    df <- df %>%
      mutate(
        chain_pair_id = mk_pair_key(chain_instance_id_A, chain_instance_id_B)
      )
  }
  
  if (!("row_id" %in% names(df))) {
    df <- df %>%
      mutate(
        row_id = paste(
          pdb_id, assembly_id, interface_id, pair_index_in_interface,
          as.character(.data[[chainA]]), as.character(.data[[chainB]]),
          sep = "|"
        )
      )
  }
  
  for (t in thrs) {
    src <- paste0("pair_cluster_c", t)
    dst <- paste0("pair_cluster_", t)
    if (!(dst %in% names(df)) && (src %in% names(df))) {
      df[[dst]] <- df[[src]]
    }
  }
  
  for (t in thrs) {
    srcA <- paste0("molA_cluster_c", t)
    srcB <- paste0("molB_cluster_c", t)
    dstA <- paste0("A_cluster_", t)
    dstB <- paste0("B_cluster_", t)
    if (!(dstA %in% names(df)) && (srcA %in% names(df))) df[[dstA]] <- df[[srcA]]
    if (!(dstB %in% names(df)) && (srcB %in% names(df))) df[[dstB]] <- df[[srcB]]
  }
  
  df
}

# -----------------------------
# 2) Quality ranking
#    Prioridad:
#    1) coverage mínima (%) entre A/B
#    2) p_value menor
#    3) interface_area mayor
#    4) r_free menor
#    5) resolution menor
#    6) solvation_energy más negativo
# -----------------------------
rank_quality <- function(g) {
  covA_pct <- if ("molA_coverage_pct" %in% names(g)) as_num(g$molA_coverage_pct) else 100 * as_num(g$molA_coverage)
  covB_pct <- if ("molB_coverage_pct" %in% names(g)) as_num(g$molB_coverage_pct) else 100 * as_num(g$molB_coverage)
  
  g %>%
    mutate(
      .covA_pct = covA_pct,
      .covB_pct = covB_pct,
      .cov_min_pct = pmin(.covA_pct, .covB_pct, na.rm = TRUE),
      .p_value_num = as_num(p_value),
      .area_num = as_num(interface_area),
      .r_free_num = as_num(r_free),
      .res_num = as_num(resolution),
      .solv_num = as_num(solvation_energy)
    ) %>%
    mutate(
      .cov_min_pct = ifelse(is.infinite(.cov_min_pct), NA_real_, .cov_min_pct)
    ) %>%
    arrange(
      desc(coalesce(.cov_min_pct, -Inf)),   # mayor coverage mínima
      coalesce(.p_value_num, Inf),          # p_value menor
      desc(coalesce(.area_num, -Inf)),      # área mayor
      coalesce(.r_free_num, Inf),           # r_free menor
      coalesce(.res_num, Inf),              # resolución menor
      coalesce(.solv_num, Inf)              # más negativo = mejor
    ) %>%
    select(-starts_with("."))
}

# -----------------------------
# 3) Exact dedupe (opcional)
# -----------------------------
dedupe_exact_pairs <- function(df,
                               mode = c("none", "keep1", "softcap"),
                               K_homo = 1L,
                               K_heter = 2L) {
  mode <- match.arg(mode)
  
  df <- df %>%
    mutate(
      is_homo01 = derive_is_homo01(df),
      pair_seq_exact = mk_pair_key(molA_seq, molB_seq)
    )
  
  if (mode == "none") return(df)
  
  if (mode == "keep1") {
    out <- df %>%
      group_by(pair_seq_exact) %>%
      group_modify(~ rank_quality(.x) %>% slice_head(n = 1)) %>%
      ungroup()
    return(out)
  }
  
  out <- df %>%
    group_by(pair_seq_exact) %>%
    group_modify(~{
      g <- rank_quality(.x)
      any_homo  <- any(g$is_homo01 == 1L, na.rm = TRUE)
      any_heter <- any(g$is_homo01 == 0L, na.rm = TRUE)
      K_keep <- if (isTRUE(any_homo)) K_homo else if (isTRUE(any_heter)) K_heter else K_heter
      slice_head(g, n = min(K_keep, nrow(g)))
    }) %>%
    ungroup()
  
  out
}

# -----------------------------
# 4) Diversidad geométrica por interface_area
# -----------------------------
select_cover_bins_area <- function(g, K, nbins = 4L, area_col = "interface_area") {
  if (nrow(g) <= K) return(rank_quality(g))
  
  gq <- rank_quality(g) %>%
    mutate(
      .area = as_num(.data[[area_col]])
    )
  
  if (all(is.na(gq$.area)) || dplyr::n_distinct(gq$.area, na.rm = TRUE) < 2) {
    return(gq %>% slice_head(n = K) %>% select(-.area))
  }
  
  nb <- min(nbins, K)
  brks <- unique(as.numeric(quantile(
    gq$.area,
    probs = seq(0, 1, length.out = nb + 1),
    na.rm = TRUE,
    names = FALSE
  )))
  
  if (length(brks) < 3) {
    return(gq %>% slice_head(n = K) %>% select(-.area))
  }
  
  gq <- gq %>%
    mutate(.bin = cut(.area, breaks = brks, include.lowest = TRUE, labels = FALSE))
  
  per_bin <- gq %>%
    group_by(.bin) %>%
    slice_head(n = 1) %>%
    ungroup()
  
  if (nrow(per_bin) >= K) {
    return(per_bin %>% slice_head(n = K) %>% select(-.area, -.bin))
  }
  
  rest <- gq %>%
    filter(!(row_id %in% per_bin$row_id)) %>%
    slice_head(n = K - nrow(per_bin))
  
  bind_rows(per_bin, rest) %>%
    slice_head(n = K) %>%
    select(-.area, -.bin)
}

# -----------------------------
# 5) Auditorías
# -----------------------------
ppi_profile <- function(df, label = "df") {
  chains_inst <- unique(c(df$chain_instance_id_A, df$chain_instance_id_B))
  tibble(
    dataset = label,
    n_rows = nrow(df),
    n_unique_row_id = n_distinct(df$row_id),
    n_unique_chain_pair_id = n_distinct(df$chain_pair_id),
    n_unique_pair_cluster_98 = if ("pair_cluster_98" %in% names(df)) n_distinct(df$pair_cluster_98) else NA_integer_,
    n_chain_instances = length(chains_inst)
  )
}

audit_homo <- function(df, label = "df") {
  v <- derive_is_homo01(df)
  tibble(
    dataset = label,
    n_rows = nrow(df),
    homo_n = sum(v == 1L, na.rm = TRUE),
    heter_n = sum(v == 0L, na.rm = TRUE),
    na_n = sum(is.na(v)),
    homo_frac = homo_n / pmax(homo_n + heter_n, 1),
    heter_frac = heter_n / pmax(homo_n + heter_n, 1)
  )
}

audit_chain_repeats <- function(df, label = "df") {
  v <- c(df$chain_instance_id_A, df$chain_instance_id_B)
  v <- v[!is.na(v) & v != ""]
  if (length(v) == 0) {
    return(tibble(
      dataset = label,
      n_chain_instances = 0L,
      n_repeated_chains = 0L,
      pct_repeated_chains = 0,
      p0 = NA_real_, p50 = NA_real_, p90 = NA_real_, p95 = NA_real_, p99 = NA_real_, p100 = NA_real_,
      mean_count = NA_real_, max_count = NA_integer_
    ))
  }
  tab <- as.data.frame(table(v), stringsAsFactors = FALSE)
  colnames(tab) <- c("chain_instance_id", "count")
  tab$count <- as.integer(tab$count)
  qs <- as.numeric(quantile(tab$count, probs = c(0, .5, .9, .95, .99, 1), names = FALSE))
  tibble(
    dataset = label,
    n_chain_instances = nrow(tab),
    n_repeated_chains = sum(tab$count > 1),
    pct_repeated_chains = 100 * mean(tab$count > 1),
    p0 = qs[1], p50 = qs[2], p90 = qs[3], p95 = qs[4], p99 = qs[5], p100 = qs[6],
    mean_count = mean(tab$count),
    max_count = max(tab$count)
  )
}

audit_pair_exact_repeats <- function(df, label = "df") {
  df2 <- df %>%
    mutate(pair_seq_exact = mk_pair_key(molA_seq, molB_seq)) %>%
    count(pair_seq_exact, name = "n")
  qs <- as.numeric(quantile(df2$n, probs = c(0, .5, .9, .95, .99, 1), names = FALSE))
  tibble(
    dataset = label,
    n_pair_groups = nrow(df2),
    n_repeated_pairs = sum(df2$n > 1),
    pct_repeated_pairs = 100 * mean(df2$n > 1),
    p0 = qs[1], p50 = qs[2], p90 = qs[3], p95 = qs[4], p99 = qs[5], p100 = qs[6],
    mean_n = mean(df2$n),
    max_n = max(df2$n)
  )
}

q_summary <- function(x) {
  x <- as.numeric(x)
  qs <- as.numeric(quantile(x, probs = c(0, .5, .9, .95, .99, 1), na.rm = TRUE, names = FALSE))
  tibble(p0 = qs[1], p50 = qs[2], p90 = qs[3], p95 = qs[4], p99 = qs[5], p100 = qs[6], mean = mean(x, na.rm = TRUE))
}

audit_sizes <- function(df, key, label) {
  stopifnot(key %in% names(df))
  sizes <- df %>% filter(!is.na(.data[[key]])) %>% count(.data[[key]], name = "n")
  bind_cols(
    tibble(
      dataset = label, key = key,
      n_groups = nrow(sizes),
      n_rows = nrow(df),
      pct_groups_gt1 = 100 * mean(sizes$n > 1),
      pct_rows_in_gt1 = 100 * sum(sizes$n[sizes$n > 1]) / pmax(nrow(df), 1)
    ),
    q_summary(sizes$n)
  )
}

audit_pair_clusters_multi <- function(df, label, keys = c("pair_cluster_98","pair_cluster_95","pair_cluster_90","pair_cluster_60","pair_cluster_50","pair_cluster_40")) {
  out <- lapply(keys, function(k) if (k %in% names(df)) audit_sizes(df, k, label) else NULL)
  bind_rows(out)
}

audit_chain_clusters_multi <- function(df, label, keys = c("98","95","90","60","50","40")) {
  get_counts <- function(a_col, b_col) {
    v <- c(df[[a_col]], df[[b_col]])
    v <- v[!is.na(v)]
    as.data.frame(table(v), stringsAsFactors = FALSE) %>%
      setNames(c("cluster","n")) %>%
      mutate(n = as.integer(n))
  }
  
  out <- lapply(keys, function(thr) {
    a_col <- paste0("A_cluster_", thr)
    b_col <- paste0("B_cluster_", thr)
    if (!all(c(a_col, b_col) %in% names(df))) return(NULL)
    sizes <- get_counts(a_col, b_col)
    bind_cols(
      tibble(
        dataset = label, key = paste0("CHAIN_cluster_", thr),
        n_groups = nrow(sizes),
        n_rows = nrow(df),
        pct_groups_gt1 = 100 * mean(sizes$n > 1),
        pct_rows_in_gt1 = 100 * sum(sizes$n[sizes$n > 1]) / pmax(2 * nrow(df), 1)
      ),
      q_summary(sizes$n)
    )
  })
  bind_rows(out)
}

audit_quality_summary <- function(df, label) {
  covA_pct <- if ("molA_coverage_pct" %in% names(df)) as_num(df$molA_coverage_pct) else 100 * as_num(df$molA_coverage)
  covB_pct <- if ("molB_coverage_pct" %in% names(df)) as_num(df$molB_coverage_pct) else 100 * as_num(df$molB_coverage)
  cov_min_pct <- pmin(covA_pct, covB_pct, na.rm = TRUE)
  cov_min_pct[is.infinite(cov_min_pct)] <- NA_real_
  
  tibble(
    dataset = label,
    n_rows = nrow(df),
    area_median = median(as_num(df$interface_area), na.rm = TRUE),
    area_p90 = as.numeric(quantile(as_num(df$interface_area), 0.90, na.rm = TRUE)),
    p_value_median = median(as_num(df$p_value), na.rm = TRUE),
    p_value_p90 = as.numeric(quantile(as_num(df$p_value), 0.90, na.rm = TRUE)),
    cov_min_pct_median = median(cov_min_pct, na.rm = TRUE),
    cov_min_pct_p10 = as.numeric(quantile(cov_min_pct, 0.10, na.rm = TRUE)),
    r_free_median = median(as_num(df$r_free), na.rm = TRUE)
  )
}

top_groups <- function(df, key, n = 20) {
  df %>%
    count(.data[[key]], name = "n") %>%
    arrange(desc(n)) %>%
    slice_head(n = n)
}

# -----------------------------
# 6) Función principal
# -----------------------------
prune_pair_cluster_single_step_v3 <- function(df,
                                              pair_key = "pair_cluster_98",
                                              
                                              # exact dedupe opcional
                                              exact_dedup_mode = c("none", "softcap", "keep1"),
                                              exact_K_homo = 1L,
                                              exact_K_heter = 2L,
                                              
                                              # pruning params (si NULL, usa defaults según pair_key)
                                              activate_n = NULL,
                                              area_range_pct_thr = 4.0,
                                              
                                              K_redund_homo = 2L,
                                              K_redund_heter = 4L,
                                              K_var_homo = NULL,
                                              K_var_heter = NULL,
                                              
                                              homo_majority_thr = 0.5,
                                              area_col = "interface_area",
                                              use_label_asym = FALSE,
                                              drop_unnamed = TRUE,
                                              audit_pair_keys = c("pair_cluster_98","pair_cluster_95","pair_cluster_90","pair_cluster_60","pair_cluster_50","pair_cluster_40"),
                                              audit_chain_thr = c("98","95","90","60","50","40")) {
  
  exact_dedup_mode <- match.arg(exact_dedup_mode)
  
  if (drop_unnamed) df <- df %>% select(-matches("^Unnamed"))
  df <- ensure_ids_and_clusters(df, use_label_asym = use_label_asym)
  
  stopifnot(pair_key %in% names(df))
  stopifnot(area_col %in% names(df))
  
  if (is.null(activate_n)) {
    activate_n <- if (pair_key == "pair_cluster_95") 12L else 8L
  }
  if (is.null(K_var_homo)) {
    K_var_homo <- if (pair_key == "pair_cluster_95") 6L else 4L
  }
  if (is.null(K_var_heter)) {
    K_var_heter <- if (pair_key == "pair_cluster_95") 10L else 8L
  }
  
  before_raw <- df
  
  df <- dedupe_exact_pairs(
    df,
    mode = exact_dedup_mode,
    K_homo = exact_K_homo,
    K_heter = exact_K_heter
  )
  
  df <- df %>%
    mutate(
      is_homo01 = derive_is_homo01(df),
      .area = as_num(.data[[area_col]])
    )
  
  cl_stats <- df %>%
    group_by(.data[[pair_key]]) %>%
    summarise(
      n = n(),
      homo_frac = mean(is_homo01 == 1L, na.rm = TRUE),
      heter_frac = mean(is_homo01 == 0L, na.rm = TRUE),
      cluster_type = case_when(
        homo_frac > homo_majority_thr ~ "homo",
        heter_frac > homo_majority_thr ~ "hetero",
        TRUE ~ "mixed/unknown"
      ),
      area_min = min(.area, na.rm = TRUE),
      area_max = max(.area, na.rm = TRUE),
      area_mean = mean(.area, na.rm = TRUE),
      area_range = area_max - area_min,
      area_range_pct = 100 * area_range / area_mean,
      .groups = "drop"
    ) %>%
    mutate(
      active = n > activate_n,
      redundant = active &
        !is.na(area_range_pct) &
        (area_range_pct <= area_range_pct_thr)
    )
  
  base_cols_no_key <- setdiff(names(df), pair_key)
  
  kept <- df %>%
    left_join(cl_stats, by = setNames(pair_key, pair_key)) %>%
    group_by(.data[[pair_key]]) %>%
    group_modify(~{
      g <- .x
      act <- isTRUE(unique(g$active)[1])
      red <- isTRUE(unique(g$redundant)[1])
      ctype <- unique(g$cluster_type)[1]
      if (length(ctype) == 0 || is.na(ctype)) ctype <- "mixed/unknown"
      
      if (!act) {
        return(g %>% select(any_of(base_cols_no_key)))
      }
      
      if (red) {
        K <- if (ctype == "homo") K_redund_homo else if (ctype == "hetero") K_redund_heter else K_redund_heter
        return(rank_quality(g) %>% slice_head(n = min(K, nrow(g))) %>% select(any_of(base_cols_no_key)))
      }
      
      K <- if (ctype == "homo") K_var_homo else if (ctype == "hetero") K_var_heter else K_var_heter
      chosen <- select_cover_bins_area(g, K = min(K, nrow(g)), nbins = 4L, area_col = area_col)
      chosen %>% select(any_of(base_cols_no_key))
    }) %>%
    ungroup()
  
  after <- kept
  
  pruning_summary <- tibble(
    pair_key = pair_key,
    exact_dedup_mode = exact_dedup_mode,
    activate_n = activate_n,
    area_range_pct_thr = area_range_pct_thr,
    K_redund_homo = K_redund_homo,
    K_redund_heter = K_redund_heter,
    K_var_homo = K_var_homo,
    K_var_heter = K_var_heter,
    
    rows_before_raw = nrow(before_raw),
    rows_after_exact_dedup = nrow(df),
    rows_after_prune = nrow(after),
    removed_total = nrow(before_raw) - nrow(after),
    removed_by_exact = nrow(before_raw) - nrow(df),
    removed_by_prune = nrow(df) - nrow(after),
    
    n_clusters = nrow(cl_stats),
    n_active = sum(cl_stats$active),
    n_redundant = sum(cl_stats$redundant),
    n_variable_active = sum(cl_stats$active & !cl_stats$redundant),
    
    n_homo_clusters = sum(cl_stats$cluster_type == "homo"),
    n_heter_clusters = sum(cl_stats$cluster_type == "hetero"),
    n_mixed_clusters = sum(cl_stats$cluster_type == "mixed/unknown")
  )
  
  audits <- list(
    profile = bind_rows(
      ppi_profile(before_raw, "before_raw"),
      ppi_profile(df, "after_exact_dedup"),
      ppi_profile(after, "after_prune")
    ),
    homo = bind_rows(
      audit_homo(before_raw, "before_raw"),
      audit_homo(df, "after_exact_dedup"),
      audit_homo(after, "after_prune")
    ),
    chain_repeats = bind_rows(
      audit_chain_repeats(before_raw, "before_raw"),
      audit_chain_repeats(df, "after_exact_dedup"),
      audit_chain_repeats(after, "after_prune")
    ),
    pair_exact_repeats = bind_rows(
      audit_pair_exact_repeats(before_raw, "before_raw"),
      audit_pair_exact_repeats(df, "after_exact_dedup"),
      audit_pair_exact_repeats(after, "after_prune")
    ),
    pair_clusters_multi = bind_rows(
      audit_pair_clusters_multi(before_raw, "before_raw", keys = audit_pair_keys),
      audit_pair_clusters_multi(df, "after_exact_dedup", keys = audit_pair_keys),
      audit_pair_clusters_multi(after, "after_prune", keys = audit_pair_keys)
    ),
    chain_clusters_multi = bind_rows(
      audit_chain_clusters_multi(before_raw, "before_raw", keys = audit_chain_thr),
      audit_chain_clusters_multi(df, "after_exact_dedup", keys = audit_chain_thr),
      audit_chain_clusters_multi(after, "after_prune", keys = audit_chain_thr)
    ),
    quality_summary = bind_rows(
      audit_quality_summary(before_raw, "before_raw"),
      audit_quality_summary(df, "after_exact_dedup"),
      audit_quality_summary(after, "after_prune")
    ),
    top_pair_clusters_before = if (pair_key %in% names(before_raw)) top_groups(before_raw, pair_key, 20) else NULL,
    top_pair_clusters_after = if (pair_key %in% names(after)) top_groups(after, pair_key, 20) else NULL
  )
  
  list(
    df_after = after,
    df_after_exact = df,
    cluster_stats = cl_stats,
    pruning_summary = pruning_summary,
    audits = audits
  )
}

# ============================================================
# 7) USO RECOMENDADO
# ============================================================

# A) Pruning local en pair_cluster_98
res98 <- prune_pair_cluster_single_step_v3(
  df,
  pair_key = pr$pair_key,
  
  # exact dedupe
  exact_dedup_mode = pr$exact_dedup_mode,
  exact_K_homo = as.integer(pr$exact_K_homo),
  exact_K_heter = as.integer(pr$exact_K_heter),
  
  # pruning
  activate_n = as.integer(pr$activate_n),
  area_range_pct_thr = pr$area_range_pct_thr,
  
  K_redund_homo = as.integer(pr$K_redund_homo),
  K_redund_heter = as.integer(pr$K_redund_heter),
  K_var_homo = as.integer(pr$K_var_homo),
  K_var_heter = as.integer(pr$K_var_heter)
)

# B) Pruning más amplio en pair_cluster_95
res95 <- prune_pair_cluster_single_step_v3(
  df,
  pair_key = "pair_cluster_95",
  
  # exact dedupe
  exact_dedup_mode = "softcap",
  exact_K_homo = 2L,
  exact_K_heter = 4L,
  
  # pruning
  activate_n = 12L,
  area_range_pct_thr = 4.0,
  
  K_redund_homo = 2L,
  K_redund_heter = 4L,
  K_var_homo = 6L,
  K_var_heter = 10L
)



# outputs
res98$pruning_summary
res95$pruning_summary

res98$audits$profile
res98$audits$homo
res98$audits$pair_exact_repeats
res98$audits$pair_clusters_multi
res98$audits$chain_clusters_multi
res98$audits$quality_summary
res98$audits$top_pair_clusters_before
res98$audits$top_pair_clusters_after

res95$audits$profile
res95$audits$homo
res95$audits$pair_exact_repeats
res95$audits$pair_clusters_multi
res95$audits$chain_clusters_multi
res95$audits$quality_summary
res95$audits$top_pair_clusters_before
res95$audits$top_pair_clusters_after

df_after98 <- res98$df_after
df_after95 <- res95$df_after

write.csv(df_after98, ensure_parent(cfg$paths$pruned), row.names = FALSE)
# write.csv(df_after95, "BeastAffinitySeq8_pruned95_v3.csv", row.names = FALSE)