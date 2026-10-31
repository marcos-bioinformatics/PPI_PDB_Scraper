library(dplyr)
library(tidyr)
library(data.table)

source(file.path("scripts", "common.R"))
cfg <- load_config()
qf <- cfg$quality_filter

# Los histogramas de este script se guardan en un PDF de informes
pdf(ensure_parent(cfg$paths$quality_plots))

df = read.csv(cfg$paths$enriched)

# Data Analysis

# df debe estar cargado (data.frame o data.table)
setDT(df)

# ----------------------------
# 1) NAs por columna
# ----------------------------
na_report <- data.table(
  column = names(df),
  n_na   = vapply(df, function(x) sum(is.na(x)), integer(1)),
  pct_na = vapply(df, function(x) mean(is.na(x)) * 100, numeric(1))
)[order(-pct_na)]

print(na_report)

# Si quieres ver sólo columnas con NA:
na_report_nonzero <- na_report[n_na > 0]
print(na_report_nonzero)

# ----------------------------
# 2) Summary para cuantitativas
# ----------------------------
# Detecta numéricas (incluye integer/numeric)
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


# Eliminación de NAs
# Cuántos NA hay por columna (y % sobre nrow)
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

# ESTUDIO DE CUANTILES ENERGETICOS

library(dplyr)
library(tidyr)

# -----------------------------
# Helper: coerción robusta de is_homo
# -----------------------------
coerce_is_homo01 <- function(x) {
  x_chr <- trimws(tolower(as.character(x)))
  dplyr::case_when(
    x_chr %in% c("1","true","t","yes","y") ~ 1L,
    x_chr %in% c("0","false","f","no","n") ~ 0L,
    TRUE                                   ~ NA_integer_
  )
}

# -----------------------------
# Añadir variable de densidad
# -----------------------------
df2 <- df %>%
  mutate(
    is_homo01 = coerce_is_homo01(is_homo),
    solvation_density = solvation_energy / interface_area
  )

# -----------------------------
# 1) Cuantiles globales
# -----------------------------
quantile_probs <- c(0, 0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 1)

global_quants <- df2 %>%
  summarise(
    across(
      c(p_value, interface_area, solvation_energy, solvation_density),
      ~ list(quantile(., probs = quantile_probs, na.rm = TRUE))
    )
  )

global_quants


# -----------------------------
# 2) Cuantiles por homo / hetero
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
    solvation_density_q  = list(quantile(solvation_density, probs = quantile_probs, na.rm = TRUE)),
    .groups = "drop"
  )

by_homo_quants

# -----------------------------
# 3) Resumen más legible en formato largo
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
  make_quant_table(df2, "solvation_energy"),
  make_quant_table(df2, "solvation_density")
)

quant_table_by_homo <- bind_rows(
  make_quant_table(filter(df2, is_homo01 == 1L), "p_value", "homo"),
  make_quant_table(filter(df2, is_homo01 == 1L), "interface_area", "homo"),
  make_quant_table(filter(df2, is_homo01 == 1L), "solvation_energy", "homo"),
  make_quant_table(filter(df2, is_homo01 == 1L), "solvation_density", "homo"),
  
  make_quant_table(filter(df2, is_homo01 == 0L), "p_value", "hetero"),
  make_quant_table(filter(df2, is_homo01 == 0L), "interface_area", "hetero"),
  make_quant_table(filter(df2, is_homo01 == 0L), "solvation_energy", "hetero"),
  make_quant_table(filter(df2, is_homo01 == 0L), "solvation_density", "hetero")
)

print(quant_table_global, n=50)
print(quant_table_by_homo, n=100)


# -----------------------------
# 4) Resumen compacto adicional
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
    
    solvation_density_median = median(solvation_density, na.rm = TRUE),
    solvation_density_p10 = quantile(solvation_density, 0.10, na.rm = TRUE),
    solvation_density_p90 = quantile(solvation_density, 0.90, na.rm = TRUE),
    .groups = "drop"
  )

summary_compact

# -----------------------------
# 5) (Opcional) mirar correlación tamaño vs energía
# -----------------------------
df2 %>%
  summarise(
    cor_area_solv = cor(interface_area, solvation_energy, use = "complete.obs"),
    cor_area_density = cor(interface_area, solvation_density, use = "complete.obs")
  )

# -----------------------------
# 6) (Opcional) histogramas rápidos
# -----------------------------
hist(df2$p_value, breaks = 100, main = "p_value", xlab = "p_value")
hist(df2$interface_area, breaks = 100, main = "interface_area", xlab = "interface_area")
hist(df2$solvation_energy, breaks = 100, main = "solvation_energy", xlab = "solvation_energy")
hist(df2$solvation_density, breaks = 100, main = "solvation_density", xlab = "solvation_energy / interface_area")

#La densidad energética (solvation_energy / interface_area) es casi independiente del tamaño
# (ver correlación con el tamaño, -0.024)

# FILTRADO
df_filt <- df2 %>%
  mutate(
    solvation_density = solvation_energy / interface_area
  ) %>%
  filter(
    molA_coverage > qf$min_coverage,
    molB_coverage > qf$min_coverage,
    r_free < qf$max_r_free,
    p_value <= qf$max_p_value,
    exptl_method == qf$exptl_method,
    molA_seq_len > qf$min_chain_length, molA_seq_len <= qf$max_chain_length,
    molB_seq_len > qf$min_chain_length, molB_seq_len <= qf$max_chain_length,
    interface_area > qf$min_interface_area,
    interface_area <= qf$max_interface_area,
    solvation_density < qf$max_solvation_density
  )

df_filt <- df_filt %>%
  mutate(tier = if_else(p_value < qf$gold_max_p_value, "gold", "silver"))

table(df_filt$tier)
table(df_filt$is_homo_80)
summary(df_filt$molA_coverage_pct)
summary(df_filt$molB_coverage_pct)
summary(df_filt$interface_area)
summary(df_filt$solvation_energy)


# HISTOGRAMAS
# Histogramas rápidos de variables de calidad
par(mfrow = c(2, 3))  # 2 filas x 3 columnas

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

hist(df_filt$solvation_density,
     breaks = 50,
     main = "solvation_density",
     xlab = "solvation_energy / interface_area",
     col = "lightsalmon",
     border = "white")

par(mfrow = c(1, 1))  # reset

library(data.table)

# df <- fread("tu_csv.csv")
df = df_filt

# ---------
# Normalizar is_homo a lógico TRUE/FALSE
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

# Comprobar que ahora sí tiene sentido
cat("\n=== tabla is_homo original ===\n")
print(table(df$is_homo, useNA = "ifany"))

cat("\n=== tabla is_homo_bool normalizada ===\n")
print(table(df$is_homo_bool, useNA = "ifany"))

# ---------
# función resumen
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
# cadenas individuales
# ---------
all_lengths <- c(df$molA_seq_len, df$molB_seq_len)

cat("\n=== GENERAL (cadenas individuales) ===\n")
print(summ_len(all_lengths))

homo_lengths <- c(df[is_homo_bool == TRUE,  molA_seq_len],
                  df[is_homo_bool == TRUE,  molB_seq_len])

cat("\n=== HOMO (cadenas individuales) ===\n")
print(summ_len(homo_lengths))

hetero_lengths <- c(df[is_homo_bool == FALSE, molA_seq_len],
                    df[is_homo_bool == FALSE, molB_seq_len])

cat("\n=== HETERO (cadenas individuales) ===\n")
print(summ_len(hetero_lengths))





write.csv(df_filt, ensure_parent(cfg$paths$quality_filtered), row.names = FALSE)
invisible(dev.off())
