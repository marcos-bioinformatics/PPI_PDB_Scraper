# =============================================================================
# Step 11 - External test set: the three classes of Park and Marcotte (2012)
#
# Input : paths.pruned (step 10)
# Output: paths.test                  test pairs (test_bucket C1 / C2 / C3)
#         paths.trainval_pool         every row that may be used for training
#         paths.test_reserved         chain clusters held out from training
#         reports_dir/11_test_set.json counts per class and cost in training rows
#
# A protein is "seen" if its chain cluster at test_set.level (40 % by default,
# connected components of step 09) has a row in the training pool.
#   C3  neither protein seen          (new proteins)
#   C2  exactly one protein seen      (a new protein with a known partner)
#   C1  both proteins seen, the pair not: no training pair joins homologues of
#       both chains at the level      (new pairs of known proteins)
# Construction:
#   1. Chain clusters are held out ("reserved"): every row that contains one
#      leaves the training pool. Pairs with both chains reserved are C3
#      candidates, pairs with one are C2 candidates. C3 homo pairs need one
#      reserved cluster, C3 hetero pairs need both, and more C2 candidates are
#      made by reserving one cluster of hetero pairs, until each class has
#      enough candidates.
#   2. C1: pair clusters whose two chain clusters keep other partners in the
#      pool are removed whole.
#   3. One representative per pair cluster at the level, chosen at random.
# Each class aims at test_set.homo_frac homo pairs; when one kind runs out the
# other fills the class (C2 is hetero by definition). Only chain clusters (C2,
# C3) and pair clusters (C1) with at most test_set.max_cluster_rows rows are
# held out, so that no large family is emptied from training. Clusters are
# drawn at random (seeded), with no preference for rare families and no
# stratification by length or area.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# minimal JSON writer for the report (named lists and scalars)
to_json <- function(x, ind = "") {
  if (is.list(x)) {
    inner <- paste0(ind, "  ")
    items <- vapply(names(x), function(k) paste0(inner, '"', k, '": ', to_json(x[[k]], inner)), character(1))
    return(paste0("{\n", paste(items, collapse = ",\n"), "\n", ind, "}"))
  }
  if (is.character(x)) return(paste0('"', x, '"'))
  if (length(x) == 1L && is.na(x)) return("null")
  format(x, scientific = FALSE, trim = TRUE)
}

source(file.path("scripts", "common.R"))
cfg <- load_config()
ts <- cfg$test_set

if (is.null(ts$seed)) stop("test_set.seed cannot be null: the test set would not be reproducible")
set.seed(as.integer(ts$seed))

LEVEL <- sprintf("%02d", as.integer(round(as.numeric(ts$level) * 100)))
N_C3 <- as.integer(ts$n_c3)
N_C2 <- as.integer(ts$n_c2)
N_C1 <- as.integer(ts$n_c1)
HOMO_FRAC <- as.numeric(ts$homo_frac)
CAP <- as.integer(ts$max_cluster_rows)

cols <- c(paste0(c("A_cluster_", "B_cluster_", "pair_cluster_"), LEVEL))
df0 <- fread(cfg$paths$pruned, na.strings = c("", "NA"), colClasses = setNames(rep("character", 3), cols))
missing <- setdiff(c("row_id", "is_homo01", cols), names(df0))
if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "),
                          " (test_set.level must be one of clustering.thresholds)")
if (anyDuplicated(df0$row_id)) stop("Duplicate row_id in the input")

df <- copy(df0)
df[, `:=`(.A = get(cols[1]), .B = get(cols[2]), .p = get(cols[3]))]
df[, `:=`(.a = pmin(.A, .B), .b = pmax(.A, .B), .homo = is_homo01 %in% 1L)]

# rows that contain each chain cluster (on either side)
size_dt <- unique(rbind(df[, .(cl = .A, row_id)], df[, .(cl = .B, row_id)]))[, .N, by = cl]
size_of <- setNames(size_dt$N, size_dt$cl)
pc <- df[, .(a = .a[1], b = .b[1], n = .N, has_homo = any(.homo), has_het = any(!.homo)), by = .p]
pc[, small_a := size_of[a] <= CAP][, small_b := size_of[b] <= CAP]

reserved <- character(0)
res <- function(x) x %chin% reserved
n_c3 <- function(kind) pc[res(a) & res(b) & (if (kind == "homo") has_homo else has_het), .N]
# C2 candidates: one cluster reserved and the other still seen in the pool
seen_now <- function() {
  p <- df[!(.A %chin% reserved | .B %chin% reserved)]
  unique(c(p$.a, p$.b))
}
n_c2 <- function() {
  seen <- seen_now()
  pc[xor(res(a), res(b)) & has_het & ifelse(res(a), b, a) %chin% seen, .N]
}

# reserve clusters from a shuffled candidate list until `enough()` holds or the list ends
reserve_until <- function(cand, pick, enough) {
  for (i in seq_len(nrow(cand))) {
    if (enough()) break
    new <- pick(cand[i])
    if (length(new)) reserved <<- union(reserved, new)
  }
  enough()
}

# --- 1. C3: homo pairs need one reserved cluster, hetero pairs both ------------
target <- c(homo = round(N_C3 * HOMO_FRAC), het = N_C3 - round(N_C3 * HOMO_FRAC))
cand_homo <- pc[a == b & has_homo & small_a][sample(.N)]
cand_het <- pc[has_het & small_a & small_b][sample(.N)]
pick_c3 <- function(r) unique(c(r$a, r$b))
ok_homo <- reserve_until(cand_homo, pick_c3, function() n_c3("homo") >= target[["homo"]])
if (!ok_homo) target[["het"]] <- N_C3 - n_c3("homo")
ok_het <- reserve_until(cand_het, pick_c3, function() n_c3("het") >= target[["het"]])
if (!ok_het && ok_homo) {
  target[["homo"]] <- N_C3 - n_c3("het")
  invisible(reserve_until(cand_homo, pick_c3, function() n_c3("homo") >= target[["homo"]]))
}

# --- 2. C2: reserve the smaller cluster of hetero pairs with neither reserved -------
cand_c2 <- pc[a != b & has_het & (small_a | small_b)][sample(.N)]
pick_c2 <- function(r) {
  if (res(r$a) || res(r$b)) return(character(0))
  live <- pc[!res(a) & !res(b)]
  d <- unique(rbind(live[, .(cl = a, .p)], live[, .(cl = b, .p)]))[, .N, by = cl]
  deg_live <- setNames(d$N, d$cl)
  # reserve a small cluster whose partner keeps other pairs (so that it stays seen)
  opts <- c(r$a, r$b)[c(r$small_a && deg_live[r$b] >= 2L, r$small_b && deg_live[r$a] >= 2L)]
  if (!length(opts)) return(character(0))
  opts[which.min(size_of[opts])]
}
invisible(reserve_until(cand_c2, pick_c2, function() n_c2() >= N_C2))

# --- 3. pool after the reservation; C2 needs its known protein in the pool ------
blocked_res <- df[.A %chin% reserved | .B %chin% reserved, row_id]
pool <- df[!row_id %chin% blocked_res]
seen <- unique(c(pool$.a, pool$.b))
c2_pcs <- pc[xor(res(a), res(b)) & has_het & ifelse(res(a), b, a) %chin% seen]
c3_pcs <- pc[res(a) & res(b)]

# --- 4. C1: pair clusters whose two chain clusters keep other partners -----------
pcR <- pool[, .(a = .a[1], b = .b[1], n = .N, has_homo = any(.homo), has_het = any(!.homo)), by = .p]
deg_dt <- unique(rbind(pcR[, .(cl = a, .p)], pcR[, .(cl = b, .p)]))[, .N, by = cl]
deg <- setNames(deg_dt$N, deg_dt$cl)
c1_sel <- list(homo = character(0), het = character(0))
take_c1 <- function(cand, kind, k) {
  for (i in seq_len(nrow(cand))) {
    if (length(c1_sel[[kind]]) >= k) break
    r <- cand[i]
    if (r$.p %chin% unlist(c1_sel)) next
    cl <- unique(c(r$a, r$b))
    if (all(deg[cl] >= 2L)) {
      deg[cl] <<- deg[cl] - 1L
      c1_sel[[kind]] <<- c(c1_sel[[kind]], r$.p)
    }
  }
}
t1 <- c(homo = round(N_C1 * HOMO_FRAC), het = N_C1 - round(N_C1 * HOMO_FRAC))
c1_homo <- pcR[a == b & has_homo & n <= CAP][sample(.N)]
c1_het <- pcR[has_het & n <= CAP][sample(.N)]
take_c1(c1_homo, "homo", t1[["homo"]])
take_c1(c1_het, "het", N_C1 - length(c1_sel$homo))
if (length(c1_sel$het) < t1[["het"]]) take_c1(c1_homo, "homo", N_C1 - length(c1_sel$het))

# --- 5. one representative per pair cluster ----------------------------------------
draw <- function(p_ids, kind, k) {
  p_ids <- if (length(p_ids) > k) sample(p_ids, k) else p_ids
  rows <- df[.p %chin% p_ids & (if (kind == "homo") .homo else !.homo)]
  rows[, .SD[sample(.N, 1L)], by = .p]
}
c3_homo_ids <- c3_pcs[has_homo == TRUE, .p]
c3_h <- draw(c3_homo_ids, "homo", min(target[["homo"]], length(c3_homo_ids)))
c3_het_ids <- setdiff(c3_pcs[has_het == TRUE, .p], c3_h$.p)
c3_t <- draw(c3_het_ids, "het", N_C3 - nrow(c3_h))
c2_r <- draw(c2_pcs$.p, "het", N_C2)
c1_h <- draw(c1_sel$homo, "homo", length(c1_sel$homo))
c1_t <- draw(c1_sel$het, "het", length(c1_sel$het))

label <- function(dt, bucket, kind) {
  dt[, `:=`(test_bucket = bucket, test_subbucket = paste0(bucket, "_", kind),
            test_unseen_chains = paste0(fifelse(.A %chin% reserved, "A", ""), fifelse(.B %chin% reserved, "B", "")))]
}
df_test <- rbindlist(list(label(c3_h, "C3", "homo"), label(c3_t, "C3", "hetero"), label(c2_r, "C2", "hetero"),
                          label(c1_h, "C1", "homo"), label(c1_t, "C1", "hetero")), use.names = TRUE)

# --- 6. training pool: no reserved cluster, no C1 pair cluster ---------------------
blocked <- union(blocked_res, df[.p %chin% unlist(c1_sel), row_id])
df_trainval_pool <- df0[!row_id %chin% blocked]
if (nrow(df_trainval_pool) == 0L) stop("The train/val pool is empty: reduce the test set sizes (test_set)")

# --- 7. checks and report -----------------------------------------------------------
tp_cl <- unique(c(df_trainval_pool[[cols[1]]], df_trainval_pool[[cols[2]]]))
stopifnot(!any(reserved %chin% tp_cl))
stopifnot(!any(df_test[test_bucket == "C1", .p] %chin% df_trainval_pool[[cols[3]]]))
stopifnot(all(df_test[test_bucket == "C1", .a %chin% tp_cl & .b %chin% tp_cl]))
stopifnot(all(df_test[test_bucket == "C2", fifelse(.A %chin% reserved, .B, .A) %chin% tp_cl]))
stopifnot(all(df_test[test_bucket == "C3", .A %chin% reserved & .B %chin% reserved]))

counts <- df_test[, .N, by = .(test_bucket, test_subbucket)][order(test_bucket, test_subbucket)]
requested <- c(C1 = N_C1, C2 = N_C2, C3 = N_C3)
for (b in names(requested)) {
  got <- df_test[test_bucket == b, .N]
  if (got < requested[[b]]) cat(sprintf("WARNING: %s has %d of the %d requested pairs\n", b, got, requested[[b]]))
}
report <- list(
  level = as.numeric(ts$level), max_cluster_rows = CAP, homo_frac = HOMO_FRAC,
  requested = as.list(requested),
  selected = setNames(as.list(counts$N), counts$test_subbucket),
  reserved_chain_clusters = length(reserved),
  candidates = list(C3_pair_clusters = nrow(c3_pcs), C2_pair_clusters = nrow(c2_pcs),
                    C1_pair_clusters_taken = length(unlist(c1_sel))),
  rows_in = nrow(df0), rows_removed_from_training = length(blocked),
  rows_removed_pct = round(100 * length(blocked) / nrow(df0), 1),
  trainval_pool_rows = nrow(df_trainval_pool)
)
cat("\n=== TEST SET ===\n"); print(counts)
cat(sprintf("Reserved chain clusters at %s %%: %d; training rows removed: %d of %d (%.1f %%)\n",
            LEVEL, length(reserved), length(blocked), nrow(df0), report$rows_removed_pct))

aux <- grep("^\\.", names(df_test), value = TRUE)
df_test[, (aux) := NULL]
# The aliases of step 10 (pair_cluster_XX, A_cluster_XX, B_cluster_XX) and the columns that repeat
# another one (is_homo01 and is_homo_80 = is_homo; coverage_pct = 100 x coverage; has_mse_sec =
# has_special, since step 08 drops every chain with a residue that is not allowed) stay internal:
# the exported table keeps the cluster columns of step 09 (molA_cluster_cXX, pair_cluster_cXX).
alias <- grep("^(pair_cluster_[0-9]+|[AB]_cluster_[0-9]+|is_homo01|is_homo_80|mol[AB]_coverage_pct|(mol[AB]_)?has_mse_sec)$", names(df_test), value = TRUE)
df_test[, (alias) := NULL]
setcolorder(df_test, c("row_id", "test_bucket", "test_subbucket", "test_unseen_chains"))
fwrite(df_test, ensure_parent(cfg$paths$test))
fwrite(df_trainval_pool, ensure_parent(cfg$paths$trainval_pool))
fwrite(data.table(level = as.numeric(ts$level), cluster = sort(reserved)), ensure_parent(cfg$paths$test_reserved))
writeLines(to_json(report), ensure_parent(file.path(cfg$paths$reports_dir, "11_test_set.json")))
