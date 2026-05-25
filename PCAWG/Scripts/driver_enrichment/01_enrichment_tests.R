library(dplyr)
library(tidyr)
library(purrr)
out_dir      <- "~/Documents/material-tickTack-2026/PCAWG/Plot/"

# ══════════════════════════════════════════════════════════════════════════════
# STATISTICAL ENRICHMENT: mutation events × sample classes (WGD / HM / Classic)
#
# Tests per (cancer_type, gene):
#   A. Global Fisher (RxC, simulated)  : all classes × all events
#   B. Pairwise Fisher (RxC, simulated): HM–Classic, HM–WGD, WGD–Classic
#   C. Per-event 2×2 Fisher            : each event enriched in HM vs rest?
#   D. Per-event pairwise 2×2 Fisher   : each event × all three class pairs
#
# FDR: two columns —
#   p_adj_local  : BH within each (IntoGen_cancer_type, test_type)
#   p_adj_global : BH across all tests of the same test_type (all cancer types)
#
# Pan-cancer: all tumor types pooled as IntoGen_cancer_type == "pan_cancer".
# ══════════════════════════════════════════════════════════════════════════════

MIN_CLASS_N  <- 3   # min samples per class to enter a test
FDR_METHOD   <- "BH"

EVENT_LEVELS <- c("CNA_driver", "CI_M", "M", "WT")
CLASS_LEVELS <- c("WGD", "HM", "Classic")
CLASS_PAIRS  <- list(c("HM", "Classic"), c("HM", "WGD"), c("WGD", "Classic"))

Drivers  <- readRDS("~/Documents/material-tickTack-2026/PCAWG/Data/Drivers.rds")

# ── 1. Collapse to one event per (cancer_type, sample, gene) ─────────────────
# Priority: CNA_driver > CI_M > M > WT  (mirrors build_heatmap_df)
event_matrix <- Drivers %>%
  filter(!is.na(mutation_status), !is.na(class)) %>%
  mutate(
    mutation_status = factor(mutation_status, levels = EVENT_LEVELS),
    class           = factor(class, levels = c(CLASS_LEVELS, "other"))
  ) %>%
  group_by(IntoGen_cancer_type, gene, sample) %>%
  arrange(mutation_status, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

# Subset to the three classes of interest
event_stat <- event_matrix %>%
  filter(class %in% CLASS_LEVELS) %>%
  mutate(class = factor(as.character(class), levels = CLASS_LEVELS))

# ── 2. Helper: build 3-row (WGD/HM/Classic) × 4-col (events) count matrix ────
make_count_mat <- function(df) {
  base <- expand.grid(
    class           = CLASS_LEVELS,
    mutation_status = EVENT_LEVELS,
    stringsAsFactors = FALSE
  )
  counts <- df %>%
    transmute(
      class           = as.character(class),
      mutation_status = as.character(mutation_status)
    ) %>%
    count(class, mutation_status)

  full_tbl <- base %>%
    left_join(counts, by = c("class", "mutation_status")) %>%
    replace_na(list(n = 0L))

  mat <- matrix(
    full_tbl$n,
    nrow     = length(CLASS_LEVELS),
    ncol     = length(EVENT_LEVELS),
    dimnames = list(CLASS_LEVELS, EVENT_LEVELS)
  )
  mat
}

# ── 3. Safe Fisher tests ──────────────────────────────────────────────────────

# RxC Fisher with Monte Carlo simulation (for global and pairwise 4-event tests)
safe_fisher_rc <- function(mat) {
  tryCatch(
    fisher.test(mat, simulate.p.value = TRUE, B = 2000)$p.value,
    error = function(e) NA_real_
  )
}

# 2×2 Fisher: returns p-value, OR, and 95% CI
# Arguments: counts as (class1_event, class1_no_event, class2_event, class2_no_event)
safe_fisher_2x2 <- function(a, b, c_, d) {
  tbl <- matrix(c(a, b, c_, d), nrow = 2,
                dimnames = list(c("cls1", "cls2"), c("event", "no_event")))
  if (any(rowSums(tbl) == 0) || all(tbl == 0)) {
    return(list(p_value = NA_real_, or = NA_real_,
                ci_low  = NA_real_, ci_high = NA_real_))
  }
  tryCatch({
    res <- fisher.test(tbl)
    list(p_value = res$p.value,
         or      = as.numeric(res$estimate),
         ci_low  = res$conf.int[1],
         ci_high = res$conf.int[2])
  }, error = function(e)
    list(p_value = NA_real_, or = NA_real_,
         ci_low  = NA_real_, ci_high = NA_real_))
}

# ── 4. Core test function for one (cancer_type, gene) group ───────────────────
# .x is the slice of event_stat for that group (already one row per sample).
run_gene_tests <- function(.x) {
  mat   <- make_count_mat(.x)
  n_cls <- setNames(rowSums(mat), rownames(mat))  # samples per class
  usable <- names(n_cls)[n_cls >= MIN_CLASS_N]

  if (length(usable) < 2) return(tibble())

  rows <- list()

  # ── A. Global (all usable classes × events present in ≥1 usable class) ──────
  m_g <- mat[usable, , drop = FALSE]
  m_g <- m_g[, colSums(m_g) > 0, drop = FALSE]
  if (nrow(m_g) >= 2 && ncol(m_g) >= 2) {
    rows[["global"]] <- tibble(
      test_type  = "global",
      comparison = paste(sort(usable), collapse = "_vs_"),
      event      = NA_character_,
      p_value    = safe_fisher_rc(m_g),
      odds_ratio = NA_real_, ci_low = NA_real_, ci_high = NA_real_
    )
  }

  # ── B. Pairwise all-events (2 classes × 4 events) ────────────────────────────
  for (pair in CLASS_PAIRS) {
    if (!all(pair %in% usable)) next
    m_p <- mat[pair, , drop = FALSE]
    m_p <- m_p[, colSums(m_p) > 0, drop = FALSE]
    if (ncol(m_p) < 2) next
    key <- paste("pair", pair[1], pair[2], sep = "_")
    rows[[key]] <- tibble(
      test_type  = "pairwise_all_events",
      comparison = paste(pair, collapse = "_vs_"),
      event      = NA_character_,
      p_value    = safe_fisher_rc(m_p),
      odds_ratio = NA_real_, ci_low = NA_real_, ci_high = NA_real_
    )
  }

  # ── C. Per-event 2×2: HM vs rest (Classic + WGD pooled) ─────────────────────
  # OR > 1  →  event enriched IN HM
  # OR < 1  →  event depleted in HM (enriched in Classic/WGD)
  if ("HM" %in% usable) {
    rest      <- intersect(c("Classic", "WGD"), usable)
    hm_total  <- n_cls["HM"]
    rst_total <- sum(n_cls[rest])
    for (ev in EVENT_LEVELS) {
      hm_ev  <- mat["HM", ev]
      rst_ev <- sum(mat[rest, ev])
      res    <- safe_fisher_2x2(hm_ev, hm_total - hm_ev,
                                rst_ev, rst_total - rst_ev)
      rows[[paste("hm_rest", ev, sep = "_")]] <- tibble(
        test_type  = "per_event_HM_vs_rest",
        comparison = "HM_vs_rest",
        event      = ev,
        p_value    = res$p_value,
        odds_ratio = res$or,
        ci_low     = res$ci_low,
        ci_high    = res$ci_high
      )
    }
  }

  # ── D. Per-event pairwise 2×2 ────────────────────────────────────────────────
  # OR > 1  →  event enriched in class1 relative to class2
  for (pair in CLASS_PAIRS) {
    if (!all(pair %in% usable)) next
    c1 <- pair[1]; c2 <- pair[2]
    for (ev in EVENT_LEVELS) {
      res <- safe_fisher_2x2(
        mat[c1, ev], n_cls[c1] - mat[c1, ev],
        mat[c2, ev], n_cls[c2] - mat[c2, ev]
      )
      key <- paste("pair_ev", c1, c2, ev, sep = "_")
      rows[[key]] <- tibble(
        test_type  = "per_event_pairwise",
        comparison = paste(pair, collapse = "_vs_"),
        event      = ev,
        p_value    = res$p_value,
        odds_ratio = res$or,
        ci_low     = res$ci_low,
        ci_high    = res$ci_high
      )
    }
  }

  # ── E. WGD + HM pooled vs Classic ──────────────────────────────────────────
  # OR > 1  →  event enriched in the non-classic pool (WGD+HM) vs Classic
  wgdhm_in_data <- intersect(c("WGD", "HM"), usable)
  n_wgdhm       <- sum(n_cls[wgdhm_in_data])

  if (length(wgdhm_in_data) > 0 && "Classic" %in% usable &&
      n_wgdhm >= MIN_CLASS_N) {

    pooled_wgdhm  <- colSums(mat[wgdhm_in_data, , drop = FALSE])
    classic_row   <- mat["Classic", ]
    m_pool        <- rbind(WGDHM = pooled_wgdhm, Classic = classic_row)
    n_classic_use <- unname(n_cls["Classic"])

    # 2×4 test
    m_pool_nz <- m_pool[, colSums(m_pool) > 0, drop = FALSE]
    if (ncol(m_pool_nz) >= 2) {
      rows[["wgdhm_classic_all"]] <- tibble(
        test_type  = "pairwise_all_events",
        comparison = "WGDHM_vs_Classic",
        event      = NA_character_,
        p_value    = safe_fisher_rc(m_pool_nz),
        odds_ratio = NA_real_, ci_low = NA_real_, ci_high = NA_real_
      )
    }

    # Per-event 2×2
    for (ev in EVENT_LEVELS) {
      res <- safe_fisher_2x2(
        pooled_wgdhm[ev], n_wgdhm       - pooled_wgdhm[ev],
        classic_row[ev],  n_classic_use  - classic_row[ev]
      )
      rows[[paste("wgdhm_classic_ev", ev, sep = "_")]] <- tibble(
        test_type  = "per_event_pairwise",
        comparison = "WGDHM_vs_Classic",
        event      = ev,
        p_value    = res$p_value,
        odds_ratio = res$or,
        ci_low     = res$ci_low,
        ci_high    = res$ci_high
      )
    }
  }

  bind_rows(rows) %>%
    mutate(
      n_WGD     = unname(n_cls["WGD"]),
      n_HM      = unname(n_cls["HM"]),
      n_Classic = unname(n_cls["Classic"]),
      n_WGDHM   = unname(sum(n_cls[c("WGD", "HM")]))
    )
}

# ── 5. Run per cancer type ────────────────────────────────────────────────────
message("Running per-cancer-type enrichment tests ...")

results_per_ct <- event_stat %>%
  group_by(IntoGen_cancer_type, gene) %>%
  group_modify(~ run_gene_tests(.x)) %>%
  ungroup()

# ── 6. Run pan-cancer (all tumor types pooled) ───────────────────────────────
message("Running pan-cancer enrichment tests ...")

results_pancancer <- event_stat %>%
  group_by(gene) %>%
  group_modify(~ run_gene_tests(.x)) %>%
  ungroup() %>%
  mutate(IntoGen_cancer_type = "pan_cancer")

# ── 7. Combine and apply FDR correction ──────────────────────────────────────
results_all <- bind_rows(results_per_ct, results_pancancer)

# p_adj_local : BH within each (cancer_type, test_type)
results_all <- results_all %>%
  group_by(IntoGen_cancer_type, test_type) %>%
  mutate(p_adj_local = p.adjust(p_value, method = FDR_METHOD)) %>%
  ungroup()

# p_adj_global: BH across all cancer types for each test_type
# (pan_cancer rows are included in the global correction pool)
results_all <- results_all %>%
  group_by(test_type) %>%
  mutate(p_adj_global = p.adjust(p_value, method = FDR_METHOD)) %>%
  ungroup()

# ── 8. Count and proportion table ────────────────────────────────────────────
count_table <- bind_rows(
  event_stat %>%
    mutate(scope = "per_cancer_type"),
  event_stat %>%
    mutate(IntoGen_cancer_type = "pan_cancer", scope = "pan_cancer")
) %>%
  group_by(IntoGen_cancer_type, gene, class, mutation_status) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(IntoGen_cancer_type, gene, class) %>%
  mutate(class_total = sum(n), prop = n / class_total) %>%
  ungroup() %>%
  rename(event = mutation_status)

# ── 9. HM-signature summary ───────────────────────────────────────────────────
# For each (gene, event): how many per-cancer-type tests show HM enrichment?
hm_signature_summary <- results_all %>%
  filter(
    IntoGen_cancer_type != "pan_cancer",
    test_type           == "per_event_HM_vs_rest",
    !is.na(odds_ratio),
    odds_ratio          >  1,
    p_adj_local         < 0.05
  ) %>%
  group_by(gene, event) %>%
  summarise(
    n_cancer_types_enriched = n(),
    cancer_types            = paste(sort(IntoGen_cancer_type), collapse = ", "),
    median_OR               = median(odds_ratio, na.rm = TRUE),
    min_p_adj_local         = if (all(is.na(p_adj_local))) NA_real_ else min(p_adj_local, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_cancer_types_enriched), event, desc(median_OR))

# Also flag pan-cancer significant hits
hm_pancancer_sig <- results_all %>%
  filter(
    IntoGen_cancer_type == "pan_cancer",
    test_type           == "per_event_HM_vs_rest",
    !is.na(odds_ratio),
    p_adj_local         < 0.05   # local == global for pan_cancer test_type group
  ) %>%
  select(gene, event, p_value, p_adj_local, odds_ratio, ci_low, ci_high,
         n_WGD, n_HM, n_Classic) %>%
  arrange(p_adj_local, desc(abs(log2(odds_ratio))))

# ── 10. Save outputs ──────────────────────────────────────────────────────────
out_dir_stats <- file.path(out_dir, "statistical_enrichment")
dir.create(out_dir_stats, showWarnings = FALSE, recursive = TRUE)

saveRDS(results_all,           file.path(out_dir_stats, "enrichment_results.rds"))
saveRDS(count_table,           file.path(out_dir_stats, "event_counts.rds"))
saveRDS(hm_signature_summary,  file.path(out_dir_stats, "hm_signature_summary.rds"))
saveRDS(hm_pancancer_sig,      file.path(out_dir_stats, "hm_pancancer_significant.rds"))

write.csv(results_all,          file.path(out_dir_stats, "enrichment_results.csv"),
          row.names = FALSE)
write.csv(count_table,          file.path(out_dir_stats, "event_counts.csv"),
          row.names = FALSE)
write.csv(hm_signature_summary, file.path(out_dir_stats, "hm_signature_summary.csv"),
          row.names = FALSE)
write.csv(hm_pancancer_sig,     file.path(out_dir_stats, "hm_pancancer_significant.csv"),
          row.names = FALSE)

# ── 11. Diagnostic print ──────────────────────────────────────────────────────
message("\n══ Results summary ════════════════════════════════════════════")
message("Total test rows:  ", nrow(results_all))
message("Unique genes tested:  ",
        n_distinct(results_all$gene[results_all$IntoGen_cancer_type != "pan_cancer"]))
message("Cancer types tested:  ",
        n_distinct(results_all$IntoGen_cancer_type) - 1L, "  (+pan_cancer)")

message("\n── Pan-cancer: events significantly enriched / depleted in HM ──")
message("  (per_event_HM_vs_rest, p_adj_local < 0.05)\n")
hm_pancancer_sig %>%
  mutate(direction = ifelse(odds_ratio > 1, "enriched in HM", "depleted in HM")) %>%
  select(gene, event, direction, odds_ratio, ci_low, ci_high,
         p_adj_local, n_HM, n_Classic, n_WGD) %>%
  print(n = 40)

message("\n── Cross-cancer HM signature: top gene×event pairs ──")
message("  (genes enriched in HM in multiple cancer types)\n")
print(hm_signature_summary, n = 30)

message("\nOutput saved to: ", out_dir_stats)
message("Files written:")
for (f in c("enrichment_results.csv", "event_counts.csv",
            "hm_signature_summary.csv", "hm_pancancer_significant.csv",
            "enrichment_results.rds", "event_counts.rds",
            "hm_signature_summary.rds", "hm_pancancer_significant.rds")) {
  message("  ", file.path(out_dir_stats, f))
}
