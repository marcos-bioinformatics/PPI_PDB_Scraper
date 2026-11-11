library(dplyr)
library(tidyr)
library(data.table)

source(file.path("scripts", "common.R"))
cfg <- load_config()
qf <- cfg$quality_filter

# The histograms of this script are saved to a PDF in the reports folder
pdf(ensure_parent(cfg$paths$quality_plots))

df = read.csv(cfg$paths$enriched)

# Data Analysis

# df must be loaded (data.frame or data.table)
setDT(df)

# ----------------------------
# 1) NAs per column
# ----------------------------
na_report <- data.table(
  column = names(df),
  n_na   = vapply(df, function(x) sum(is.na(x)), integer(1)),
  pct_na = vapply(df, function(x) mean(is.na(x)) * 100, numeric(1))
)[order(-pct_na)]

print(na_report)

# To see only the columns with NAs:
na_report_nonzero <- na_report[n_na > 0]
print(na_report_nonzero)

# ----------------------------
# 2) Summary of quantitative variables
# ----------------------------
# Detect numeric columns (integer/numeric included)
num_cols <- names(df)[vapply(df, is.numeric, logical(1))]

num_summary <- rbindlist(lapply(num_cols, function(cn) {
  x <- df[[cn]]
  data.table(
    column = cn,
    n      = length(x),
    n_na   = sum(is.na(x)),
    min    = suppressWarnings(min(x, na.rm = TRUE)),
    q1     = suppressWarnings(as.numeric(quantile(x, 0.25, na.rm = TRUE))),
    median = suppressWarnings(median(x, na.rm = TRUE)),
    mean   = suppressWarnings(mean(x, na.rm = TRUE)),
    q3     = suppressWarnings(as.numeric(quantile(x, 0.75, na.rm = TRUE))),
    max    = suppressWarnings(max(x, na.rm = TRUE)),
    sd     = suppressWarnings(sd(x, na.rm = TRUE))
  )
}), fill = TRUE)


# NA removal
# How many NAs per column (and % of nrow)
na_por_col <- df %>%
  summarise(across(
    everything(),
    list(n_na = ~sum(is.na(.)), pct_na = ~mean(is.na(.)) * 100)
  )) %>%
  pivot_longer(
    everything(),
    names_to = c("col", ".value"),
    names_pattern = "^(.*)_(n_na|pct_na)$"
  ) %>%
  arrange(desc(n_na))

print(na_por_col, n=100)

cols_clave <- c("molA_coverage","molB_coverage","r_free","exptl_method",
                "molA_seq_len","molB_seq_len","solvation_energy","interface_area")

df2 <- df %>%
  filter(if_all(all_of(cols_clave), ~ !is.na(.)))


# STUDY OF ENERGY QUANTILES

# -----------------------------
# Helper: robust coercion of is_homo
# -----------------------------
coerce_is_homo01 <- function(x) {
  x_chr <- trimws(tolower(as.character(x)))
  dplyr::case_when(
    x_chr %in% c("1","true","t","yes","y") ~ 1L,
    x_chr %in% c("0","false","f","no","n") ~ 0L,
    TRUE                                   ~ NA_integer_
  )
}

df2 <- df %>%
  mutate(is_homo01 = coerce_is_homo01(is_homo))

# -----------------------------
# 1) Global quantiles
# -----------------------------
quantile_probs <- c(0, 0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 1)

global_quants <- df2 %>%
  summarise(
    across(
      c(p_value, interface_area, solvation_energy),
      ~ list(quantile(., probs = quantile_probs, na.rm = TRUE))
    )
  )

global_quants


# -----------------------------
# 2) Quantiles by homo / hetero
# -----------------------------
by_homo_quants <- df2 %>%
  mutate(
    group = case_when(
      is_homo01 == 1L ~ "homo",
      is_homo01 == 0L ~ "hetero",
      TRUE            ~ "unknown"
    )
  ) %>%
  group_by(group) %>%
  summarise(
    p_value_q            = list(quantile(p_value, probs = quantile_probs, na.rm = TRUE)),
    interface_area_q     = list(quantile(interface_area, probs = quantile_probs, na.rm = TRUE)),
    solvation_energy_q   = list(quantile(solvation_energy, probs = quantile_probs, na.rm = TRUE)),
    .groups = "drop"
  )

# -----------------------------
# 3) More readable summary in long format
# -----------------------------
make_quant_table <- function(data, var, group_name = NULL) {
  tibble(
    variable = var,
    quantile = quantile_probs,
    value = as.numeric(quantile(data[[var]], probs = quantile_probs, na.rm = TRUE)),
    group = ifelse(is.null(group_name), "global", group_name)
  )
}

quant_table_global <- bind_rows(
  make_quant_table(df2, "p_value"),
  make_quant_table(df2, "interface_area"),
  make_quant_table(df2, "solvation_energy")
)

quant_table_by_homo <- bind_rows(
  make_quant_table(filter(df2, is_homo01 == 1L), "p_value", "homo"),
  make_quant_table(filter(df2, is_homo01 == 1L), "interface_area", "homo"),
  make_quant_table(filter(df2, is_homo01 == 1L), "solvation_energy", "homo"),
  
  make_quant_table(filter(df2, is_homo01 == 0L), "p_value", "hetero"),
  make_quant_table(filter(df2, is_homo01 == 0L), "interface_area", "hetero"),
  make_quant_table(filter(df2, is_homo01 == 0L), "solvation_energy", "hetero")
)

print(quant_table_global, n=50)
print(quant_table_by_homo, n=100)


# -----------------------------
# 4) Additional compact summary
# -----------------------------
summary_compact <- df2 %>%
  mutate(
    group = case_when(
      is_homo01 == 1L ~ "homo",
      is_homo01 == 0L ~ "hetero",
      TRUE            ~ "unknown"
    )
  ) %>%
  group_by(group) %>%
  summarise(
    n = n(),
    p_value_median = median(p_value, na.rm = TRUE),
    p_value_p90 = quantile(p_value, 0.90, na.rm = TRUE),
    
    interface_area_median = median(interface_area, na.rm = TRUE),
    interface_area_p10 = quantile(interface_area, 0.10, na.rm = TRUE),
    interface_area_p90 = quantile(interface_area, 0.90, na.rm = TRUE),
    
    solvation_energy_median = median(solvation_energy, na.rm = TRUE),
    solvation_energy_p10 = quantile(solvation_energy, 0.10, na.rm = TRUE),
    .groups = "drop"
  )

summary_compact

# -----------------------------
# 5) (Optional) look at the size vs energy correlation
# -----------------------------
df2 %>%
  summarise(
    cor_area_solv = cor(interface_area, solvation_energy, use = "complete.obs")
  )

# -----------------------------
# 6) (Optional) quick histograms
# -----------------------------
hist(df2$p_value, breaks = 100, main = "p_value", xlab = "p_value")
hist(df2$interface_area, breaks = 100, main = "interface_area", xlab = "interface_area")
hist(df2$solvation_energy, breaks = 100, main = "solvation_energy", xlab = "solvation_energy")

# Entries without R-free are excluded (there is no cross-validated statistic of the
# refinement); report how many, so the exclusion is explicit.
no_rfree <- df2 %>%
  filter(is.na(r_free), if (is.null(qf$exptl_method)) TRUE else exptl_method == qf$exptl_method)
cat(sprintf("\nExcluded for missing R-free: %d candidate pairs in %d entries\n",
            nrow(no_rfree), n_distinct(no_rfree$pdb_id)))

# Scope (config `scope`): release-date cut-off, membrane proteins and antibodies.
# Antibody chains are always flagged; their pairs are dropped only on request.
sc <- cfg$scope
as_flag <- function(x) toupper(trimws(as.character(x))) %in% c("TRUE", "T", "1")
has_source <- function(x, sources) {
  vapply(strsplit(ifelse(is.na(x), "", as.character(x)), ";", fixed = TRUE),
         function(v) any(v %in% sources), logical(1))
}
membrane_sources <- unlist(sc$exclude_membrane)
df2 <- df2 %>%
  mutate(
    .in_date = if (is.null(sc$max_release_date)) TRUE else
      !is.na(initial_release_date) &
      as.Date(as.character(initial_release_date)) <= as.Date(as.character(sc$max_release_date)),
    .membrane = has_source(molA_membrane_sources, membrane_sources) |
      has_source(molB_membrane_sources, membrane_sources),
    .antibody = as_flag(molA_is_antibody) | as_flag(molB_is_antibody)
  )
cat(sprintf(paste0("Scope: %d pairs released after the cut-off (%s); %d with a membrane-protein chain ",
                   "(%s, excluded); %d with an antibody chain (%s)\n"),
            sum(!df2$.in_date), if (is.null(sc$max_release_date)) "none" else as.character(sc$max_release_date),
            sum(df2$.membrane), if (length(membrane_sources)) paste(membrane_sources, collapse = ", ") else "none",
            sum(df2$.antibody), if (isTRUE(sc$exclude_antibodies)) "excluded" else "flagged only"))

# FILTERING
df_filt <- df2 %>%
  filter(
    .in_date,
    !.membrane,
    if (isTRUE(sc$exclude_antibodies)) !.antibody else TRUE,
    !is.na(r_free),
    # coverage_metric: "coverage" (resolved / SEQRES) or "internal_coverage"
    .data[[paste0("molA_", qf$coverage_metric)]] > qf$min_coverage,
    .data[[paste0("molB_", qf$coverage_metric)]] > qf$min_coverage,
    r_free < qf$max_r_free,
    # optional filters: null in the config = disabled
    if (is.null(qf$max_resolution)) TRUE else resolution <= qf$max_resolution,
    p_value <= qf$max_p_value,
    if (is.null(qf$exptl_method)) TRUE else exptl_method == qf$exptl_method,
    molA_seq_len > qf$min_chain_length, molA_seq_len <= qf$max_chain_length,
    molB_seq_len > qf$min_chain_length, molB_seq_len <= qf$max_chain_length,
    interface_area > qf$min_interface_area,
    interface_area <= qf$max_interface_area
  )

df_filt <- df_filt %>%
  select(-starts_with(".")) %>%
  mutate(tier = if_else(p_value < qf$gold_max_p_value, "gold", "silver"))

table(df_filt$tier)
table(df_filt$is_homo_80)
summary(df_filt$molA_coverage_pct)
summary(df_filt$molB_coverage_pct)
summary(df_filt$interface_area)
summary(df_filt$solvation_energy)


# HISTOGRAMS
# Quick histograms of quality variables
par(mfrow = c(2, 3))  # 2 rows x 3 columns

hist(df_filt$molA_coverage_pct,
     breaks = 50,
     main = "molA_coverage_pct",
     xlab = "Coverage %",
     col = "lightblue",
     border = "white")

hist(df_filt$molB_coverage_pct,
     breaks = 50,
     main = "molB_coverage_pct",
     xlab = "Coverage %",
     col = "lightgreen",
     border = "white")

hist(df_filt$r_free,
     breaks = 50,
     main = "r_free",
     xlab = "r_free",
     col = "lightcoral",
     border = "white")

hist(df_filt$p_value,
     breaks = 50,
     main = "p_value",
     xlab = "p_value",
     col = "khaki",
     border = "white")

hist(df_filt$interface_area,
     breaks = 50,
     main = "interface_area",
     xlab = "Interface area",
     col = "plum",
     border = "white")

par(mfrow = c(1, 1))  # reset

df = df_filt

# ---------
# Normalise is_homo to logical TRUE/FALSE
# ---------
to_bool <- function(x) {
  x <- as.character(x)
  x <- trimws(tolower(x))
  fifelse(
    x %in% c("true", "t", "1", "yes", "y"),
    TRUE,
    fifelse(
      x %in% c("false", "f", "0", "no", "n"),
      FALSE,
      NA
    )
  )
}

df[, is_homo_bool := to_bool(is_homo)]

# Check that it makes sense now
cat("\n=== original is_homo table ===\n")
print(table(df$is_homo, useNA = "ifany"))

cat("\n=== normalised is_homo_bool table ===\n")
print(table(df$is_homo_bool, useNA = "ifany"))

# ---------
# summary function
# ---------
summ_len <- function(x) {
  x <- x[!is.na(x)]
  c(
    n    = length(x),
    min  = min(x),
    mean = mean(x),
    p90  = as.numeric(quantile(x, 0.90)),
    p95  = as.numeric(quantile(x, 0.95)),
    p96  = as.numeric(quantile(x, 0.96)),
    p97  = as.numeric(quantile(x, 0.97)),
    p98  = as.numeric(quantile(x, 0.98)),
    p99  = as.numeric(quantile(x, 0.99)),
    max  = max(x)
  )
}

# ---------
# individual chains
# ---------
all_lengths <- c(df$molA_seq_len, df$molB_seq_len)

cat("\n=== GENERAL (individual chains) ===\n")
print(summ_len(all_lengths))

homo_lengths <- c(df[is_homo_bool == TRUE,  molA_seq_len],
                  df[is_homo_bool == TRUE,  molB_seq_len])

cat("\n=== HOMO (individual chains) ===\n")
print(summ_len(homo_lengths))

hetero_lengths <- c(df[is_homo_bool == FALSE, molA_seq_len],
                    df[is_homo_bool == FALSE, molB_seq_len])

cat("\n=== HETERO (individual chains) ===\n")
print(summ_len(hetero_lengths))



write.csv(df_filt, ensure_parent(cfg$paths$quality_filtered), row.names = FALSE)
invisible(dev.off())
