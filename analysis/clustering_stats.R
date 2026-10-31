library(dplyr)
library(tidyr)

# ----------------------------
# 0) Load + sanitize
# ----------------------------
source(file.path("scripts", "common.R"))
cfg <- load_config()

df <- read.csv(cfg$paths$clustered) %>%
  select(-matches("^Unnamed"))

# ----------------------------
# 1) Ensure IDs + cluster aliases (NO FAIL)
# ----------------------------
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
  
  # chain_instance_id_A/B
  if (!all(c("chain_instance_id_A","chain_instance_id_B") %in% names(df))) {
    df <- df %>%
      mutate(
        chain_instance_id_A = paste(pdb_id, assembly_id, as.character(.data[[chainA]]), sep="|"),
        chain_instance_id_B = paste(pdb_id, assembly_id, as.character(.data[[chainB]]), sep="|")
      )
  }
  
  # chain_pair_id
  if (!("chain_pair_id" %in% names(df))) {
    df <- df %>%
      mutate(
        chain_pair_id = paste(
          pmin(chain_instance_id_A, chain_instance_id_B),
          pmax(chain_instance_id_A, chain_instance_id_B),
          sep="||"
        )
      )
  }
  
  # row_id
  if (!("row_id" %in% names(df))) {
    df <- df %>%
      mutate(
        row_id = paste(
          pdb_id, assembly_id, interface_id, pair_index_in_interface,
          as.character(.data[[chainA]]), as.character(.data[[chainB]]),
          sep="|"
        )
      )
  }
  
  # Pair-cluster aliases: pair_cluster_cXX -> pair_cluster_XX
  for (t in thrs) {
    src <- paste0("pair_cluster_c", t)
    dst <- paste0("pair_cluster_", t)
    if (!(dst %in% names(df)) && (src %in% names(df))) {
      df[[dst]] <- df[[src]]
    }
  }
  
  # Chain-cluster aliases:
  # molA_cluster_cXX -> A_cluster_XX ; molB_cluster_cXX -> B_cluster_XX
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

df <- ensure_ids_and_clusters(df, use_label_asym = FALSE)

# ----------------------------
# 2) Small profiles / helpers
# ----------------------------
ppi_profile <- function(df) {
  chains_inst <- unique(c(df$chain_instance_id_A, df$chain_instance_id_B))
  tibble(
    n_rows = nrow(df),
    n_unique_row_id = n_distinct(df$row_id),
    n_unique_chain_pair_id = n_distinct(df$chain_pair_id),
    n_unique_pair_cluster_98 = if ("pair_cluster_98" %in% names(df)) n_distinct(df$pair_cluster_98) else NA_integer_,
    n_chain_instances = length(chains_inst)
  )
}

coerce_is_homo01 <- function(x) {
  x_chr <- trimws(tolower(as.character(x)))
  dplyr::case_when(
    x_chr %in% c("1","true","t","yes","y")  ~ 1L,
    x_chr %in% c("0","false","f","no","n")  ~ 0L,
    TRUE                                    ~ NA_integer_
  )
}

# ranking: preferir no especiales, p_value mínimo, y luego desempates
rank_and_keep1 <- function(df) {
  df %>%
    mutate(
      special_rank = case_when(
        isTRUE(has_special) ~ 2L,
        isFALSE(has_special) ~ 0L,
        TRUE ~ 1L
      ),
      cov_min = pmin(molA_coverage, molB_coverage, na.rm = TRUE)
    ) %>%
    arrange(
      special_rank,
      coalesce(p_value, Inf),
      coalesce(r_free, Inf),
      coalesce(resolution, Inf),
      coalesce(solvation_energy, Inf),
      desc(coalesce(interface_area, -Inf)),
      desc(coalesce(cov_min, -Inf))
    ) %>%
    slice_head(n = 1) %>%
    select(-special_rank, -cov_min)
}

rank_key <- function(df) {
  df %>%
    mutate(
      special_rank = case_when(
        isFALSE(has_special) ~ 0L,
        is.na(has_special)   ~ 1L,
        isTRUE(has_special)  ~ 2L,
        TRUE                 ~ 1L
      ),
      cov_min = pmin(molA_coverage, molB_coverage, na.rm = TRUE)
    ) %>%
    arrange(
      special_rank,
      coalesce(p_value, Inf),
      coalesce(r_free, Inf),
      coalesce(resolution, Inf),
      coalesce(solvation_energy, Inf),
      desc(coalesce(interface_area, -Inf)),
      desc(coalesce(cov_min, -Inf))
    ) %>%
    select(-special_rank, -cov_min)
}

# ----------------------------
# 3) Dedup 1: exact sequence pair
# ----------------------------
df2 <- df %>%
  mutate(
    seqA_id = molA_seq,
    seqB_id = molB_seq,
    pair_seq_exact = paste(pmin(seqA_id, seqB_id), pmax(seqA_id, seqB_id), sep="||")
  ) %>%
  group_by(pair_seq_exact) %>%
  group_modify(~ rank_and_keep1(.x)) %>%
  ungroup() %>%
  select(-seqA_id, -seqB_id, -pair_seq_exact)

cat("DEDUP por par exacto de secuencia:\n"); print(ppi_profile(df2))
cat("\nFilas eliminadas:", nrow(df) - nrow(df2), "\n")

# ----------------------------
# 4) Dedup 2: "soft cap" per exact pair (K homo/hetero)
# ----------------------------
K_homo  <- 1L
K_heter <- 2L

df_dedup_soft <- df %>%
  mutate(
    is_homo = coerce_is_homo01(is_homo),
    pair_seq_exact = paste(pmin(molA_seq, molB_seq), pmax(molA_seq, molB_seq), sep="||")
  ) %>%
  group_by(pair_seq_exact) %>%
  group_modify(~{
    g <- rank_key(.x)
    any_homo  <- any(g$is_homo == 1L, na.rm = TRUE)
    any_heter <- any(g$is_homo == 0L, na.rm = TRUE)
    K_keep <- if (isTRUE(any_homo)) K_homo else if (isTRUE(any_heter)) K_heter else K_heter
    slice_head(g, n = min(K_keep, nrow(g)))
  }) %>%
  ungroup()

cat("Filas eliminadas:", nrow(df) - nrow(df_dedup_soft), "\n")

# ----------------------------
# 5) Audits (NO FAIL)
# ----------------------------
audit_homo <- function(df, label = "df") {
  x <- coerce_is_homo01(df$is_homo)
  homo_n <- sum(x == 1L, na.rm = TRUE)
  heter_n <- sum(x == 0L, na.rm = TRUE)
  denom <- max(homo_n + heter_n, 1)
  tibble(
    dataset = label,
    n_rows = nrow(df),
    homo_n = homo_n,
    heter_n = heter_n,
    na_n = sum(is.na(x)),
    homo_frac = homo_n / denom,
    heter_frac = heter_n / denom
  )
}

audit_chain_repeats <- function(df, label = "df") {
  stopifnot(all(c("chain_instance_id_A","chain_instance_id_B") %in% names(df)))
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
    mutate(pair_seq_exact = paste(pmin(molA_seq, molB_seq), pmax(molA_seq, molB_seq), sep="||")) %>%
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

before <- df
after  <- df_dedup_soft

bind_rows(audit_homo(before, "before"), audit_homo(after, "after"))
bind_rows(audit_chain_repeats(before, "before"), audit_chain_repeats(after, "after"))
bind_rows(audit_pair_exact_repeats(before, "before"), audit_pair_exact_repeats(after, "after"))

top_chain_repeats <- function(df, n = 20) {
  stopifnot(all(c("chain_instance_id_A","chain_instance_id_B") %in% names(df)))
  v <- c(df$chain_instance_id_A, df$chain_instance_id_B)
  v <- v[!is.na(v) & v != ""]
  if (length(v) == 0) return(tibble(chain_instance_id=character(), count=integer()))
  as.data.frame(table(v), stringsAsFactors = FALSE) %>%
    rename(chain_instance_id = v, count = Freq) %>%
    mutate(count = as.integer(count)) %>%
    arrange(desc(count)) %>%
    slice_head(n = n)
}

top_chain_repeats(before, 20)
top_chain_repeats(after, 20)

# ----------------------------
# 6) Study clusters (fixed to YOUR column names)
# ----------------------------
q_summary <- function(x) {
  x <- as.numeric(x)
  qs <- quantile(x, probs = c(0, .5, .9, .95, .99, 1), na.rm = TRUE, names = FALSE)
  tibble(p0=qs[1], p50=qs[2], p90=qs[3], p95=qs[4], p99=qs[5], p100=qs[6], mean=mean(x, na.rm=TRUE))
}

audit_sizes <- function(df, key, label) {
  stopifnot(key %in% names(df))
  sizes <- df %>% count(.data[[key]], name="n")
  bind_cols(
    tibble(dataset=label, key=key,
           n_groups=nrow(sizes),
           n_rows=nrow(df),
           pct_groups_gt1=100*mean(sizes$n>1),
           pct_rows_in_gt1=100*sum(sizes$n[sizes$n>1])/nrow(df)),
    q_summary(sizes$n)
  )
}

# Uses pair_cluster_95/90/60/50/40 which we created as aliases already
audit_pair_cluster_dual <- function(df, label) {
  out <- list()
  for (k in c(98,95,90,60,50,40)) {
    key <- paste0("pair_cluster_", k)
    if (key %in% names(df)) out[[length(out)+1]] <- audit_sizes(df, key, label)
  }
  bind_rows(out)
}

audit_chain_clusters <- function(df, label) {
  get_counts <- function(k) {
    a <- df[[paste0("A_cluster_", k)]]
    b <- df[[paste0("B_cluster_", k)]]
    v <- c(a,b)
    v <- v[!is.na(v)]
    if (length(v) == 0) return(tibble(cluster=character(), n=integer()))
    as.data.frame(table(v), stringsAsFactors = FALSE) %>%
      rename(cluster=v, n=Freq) %>%
      mutate(n = as.integer(n))
  }
  
  res <- list()
  for (k in c(98,95,90,60,50,40)) {
    if (!all(c(paste0("A_cluster_",k), paste0("B_cluster_",k)) %in% names(df))) next
    sizes <- get_counts(k)
    if (nrow(sizes) == 0) next
    res[[length(res)+1]] <- bind_cols(
      tibble(dataset=label, key=paste0("CHAIN_cluster_",k),
             n_groups=nrow(sizes),
             n_rows=nrow(df),
             pct_groups_gt1=100*mean(sizes$n>1),
             pct_rows_in_gt1=100*sum(sizes$n[sizes$n>1]) / (2*nrow(df))),
      q_summary(sizes$n)
    )
  }
  bind_rows(res)
}

# Define df_dedup24/df_dedup12 so later blocks don't crash
df_dedup24 <- df_dedup_soft
df_dedup12 <- df_dedup_soft

pair_a_before <- audit_pair_cluster_dual(df, "before")
pair_a_24     <- audit_pair_cluster_dual(df_dedup24, "after_2_4")
pair_a_12     <- audit_pair_cluster_dual(df_dedup12, "after_1_2")
pair_a <- bind_rows(pair_a_before, pair_a_24, pair_a_12)
pair_a

chain_a <- bind_rows(
  audit_chain_clusters(df, "before"),
  audit_chain_clusters(df_dedup24, "after_2_4"),
  audit_chain_clusters(df_dedup12, "after_1_2")
)
chain_a
