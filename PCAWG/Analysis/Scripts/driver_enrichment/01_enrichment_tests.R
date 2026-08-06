library(dplyr)
library(tidyr)
library(purrr)

out_dir <- "../../Plot/"

# ══════════════════════════════════════════════════════════════════════════════
# STATISTICAL ENRICHMENT: mutation events × sample classes (WGD / HM / Classic)
#
# KEY CHANGE vs previous version:
# Events are tested as INDEPENDENT binary indicators (is event X present
# in this sample for this gene?) rather than collapsing to one event per
# sample via priority ordering. This prevents CNA_driver from masking CI_M
# in samples that carry both.
#
# WT is kept as a category read directly from mutation_status (amplified
# gene with no somatic mutation), NOT derived as absence of other events.
#
# Zero-count cells are flagged with zero_count = TRUE so they remain
# visible in plots even though OR and p-value are NA.
# ══════════════════════════════════════════════════════════════════════════════

MIN_CLASS_N  <- 3
FDR_METHOD   <- "BH"

EVENT_LEVELS <- c("CNA_driver", "CI_M", "M", "WT")
CLASS_LEVELS <- c("WGD", "HM", "Classic")
CLASS_PAIRS  <- list(c("HM", "Classic"), c("HM", "WGD"), c("WGD", "Classic"))

Drivers <- readRDS("../../Data/Drivers.rds")

# ── 1. Build binary event table ───────────────────────────────────────────────
# One row per (cancer_type, gene, sample, class, event).
# A sample appears in multiple event rows if it carries multiple events.
# WT is read directly from mutation_status, not derived.

sample_events <- Drivers %>%
  filter(!is.na(mutation_status), !is.na(class)) %>%
  filter(class %in% CLASS_LEVELS) %>%
  mutate(class = factor(as.character(class), levels = CLASS_LEVELS)) %>%
  select(IntoGen_cancer_type, gene, sample, class, mutation_status) %>%
  distinct()

# For each (cancer_type, gene, sample, class), determine which events
# are present as binary indicators.
sample_class <- sample_events %>%
  group_by(IntoGen_cancer_type, gene, sample, class) %>%
  summarise(
    has_CNA_driver = any(mutation_status == "CNA_driver"),
    has_CI_M       = any(mutation_status == "CI_M"),
    has_M          = any(mutation_status == "M"),
    has_WT         = any(mutation_status == "WT"),   # read from data, not derived
    .groups = "drop"
  )

# Pivot to long format
event_binary <- sample_class %>%
  pivot_longer(
    cols      = starts_with("has_"),
    names_to  = "event",
    values_to = "present"
  ) %>%
  mutate(
    event = recode(event,
                   has_CNA_driver = "CNA_driver",
                   has_CI_M       = "CI_M",
                   has_M          = "M",
                   has_WT         = "WT"),
    event = factor(event, levels = EVENT_LEVELS)
  )

event_stat <- event_binary %>%
  filter(class %in% CLASS_LEVELS) %>%
  mutate(class = factor(as.character(class), levels = CLASS_LEVELS))

# ── 2. Helper: build count matrix and class sizes ─────────────────────────────
make_count_mat <- function(df) {
  
  # total samples per class (each sample counted once regardless of events)
  class_sizes <- df %>%
    distinct(sample, class) %>%
    count(class, name = "n_samples")
  
  # samples carrying each event per class
  counts <- df %>%
    filter(present) %>%
    distinct(sample, class, event) %>%
    count(class, event, name = "n")
  
  base <- expand.grid(
    class = CLASS_LEVELS,
    event = EVENT_LEVELS,
    stringsAsFactors = FALSE
  )
  
  full_tbl <- base %>%
    left_join(counts,      by = c("class", "event")) %>%
    left_join(class_sizes, by = "class") %>%
    replace_na(list(n = 0L, n_samples = 0L))
  
  mat <- matrix(
    full_tbl$n,
    nrow     = length(CLASS_LEVELS),
    ncol     = length(EVENT_LEVELS),
    dimnames = list(CLASS_LEVELS, EVENT_LEVELS)
  )
  
  n_cls <- setNames(
    class_sizes$n_samples[match(CLASS_LEVELS, class_sizes$class)],
    CLASS_LEVELS
  )
  n_cls[is.na(n_cls)] <- 0L
  
  list(mat = mat, n_cls = n_cls)
}

# ── 3. Safe Fisher tests ──────────────────────────────────────────────────────

safe_fisher_rc <- function(mat) {
  tryCatch(
    fisher.test(mat, simulate.p.value = TRUE, B = 2000)$p.value,
    error = function(e) NA_real_
  )
}

# Modified: returns zero_count flag alongside stats.
# zero_count = TRUE when at least one class in the 2x2 has zero event carriers.
# OR and p-value are NA in that case but the row is retained.
safe_fisher_2x2 <- function(a, b, c_, d) {
  zero_count <- (a == 0L || c_ == 0L)
  tbl <- matrix(c(a, b, c_, d), nrow = 2,
                dimnames = list(c("cls1", "cls2"), c("event", "no_event")))
  if (any(rowSums(tbl) == 0) || all(tbl == 0)) {
    return(list(p_value    = NA_real_,
                or         = NA_real_,
                ci_low     = NA_real_,
                ci_high    = NA_real_,
                zero_count = zero_count))
  }
  tryCatch({
    res <- fisher.test(tbl)
    list(p_value    = res$p.value,
         or         = as.numeric(res$estimate),
         ci_low     = res$conf.int[1],
         ci_high    = res$conf.int[2],
         zero_count = zero_count)
  }, error = function(e)
    list(p_value    = NA_real_,
         or         = NA_real_,
         ci_low     = NA_real_,
         ci_high    = NA_real_,
         zero_count = zero_count))
}

# ── 4. Core test function for one (cancer_type, gene) group ──────────────────

run_gene_tests <- function(.x) {
  res_list <- make_count_mat(.x)
  mat      <- res_list$mat
  n_cls    <- res_list$n_cls
  usable   <- names(n_cls)[n_cls >= MIN_CLASS_N]
  
  if (length(usable) < 2) return(tibble())
  
  rows <- list()
  
  # ── A. Global ────────────────────────────────────────────────────────────
  m_g <- mat[usable, , drop = FALSE]
  m_g <- m_g[, colSums(m_g) > 0, drop = FALSE]
  if (nrow(m_g) >= 2 && ncol(m_g) >= 2) {
    rows[["global"]] <- tibble(
      test_type  = "global",
      comparison = paste(sort(usable), collapse = "_vs_"),
      event      = NA_character_,
      p_value    = safe_fisher_rc(m_g),
      odds_ratio = NA_real_,
      ci_low     = NA_real_,
      ci_high    = NA_real_,
      zero_count = FALSE
    )
  }
  
  # ── B. Pairwise all-events ────────────────────────────────────────────────
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
      odds_ratio = NA_real_,
      ci_low     = NA_real_,
      ci_high    = NA_real_,
      zero_count = FALSE
    )
  }
  
  # ── C. Per-event 2×2: HM vs rest ─────────────────────────────────────────
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
        ci_high    = res$ci_high,
        zero_count = res$zero_count
      )
    }
  }
  
  # ── D. Per-event pairwise 2×2 ─────────────────────────────────────────────
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
        ci_high    = res$ci_high,
        zero_count = res$zero_count
      )
    }
  }
  
  # ── E. WGD + HM pooled vs Classic ────────────────────────────────────────
  wgdhm_in_data <- intersect(c("WGD", "HM"), usable)
  n_wgdhm       <- sum(n_cls[wgdhm_in_data])
  
  if (length(wgdhm_in_data) > 0 && "Classic" %in% usable &&
      n_wgdhm >= MIN_CLASS_N) {
    
    pooled_wgdhm  <- colSums(mat[wgdhm_in_data, , drop = FALSE])
    classic_row   <- mat["Classic", ]
    m_pool        <- rbind(WGDHM = pooled_wgdhm, Classic = classic_row)
    n_classic_use <- unname(n_cls["Classic"])
    
    m_pool_nz <- m_pool[, colSums(m_pool) > 0, drop = FALSE]
    if (ncol(m_pool_nz) >= 2) {
      rows[["wgdhm_classic_all"]] <- tibble(
        test_type  = "pairwise_all_events",
        comparison = "WGDHM_vs_Classic",
        event      = NA_character_,
        p_value    = safe_fisher_rc(m_pool_nz),
        odds_ratio = NA_real_,
        ci_low     = NA_real_,
        ci_high    = NA_real_,
        zero_count = FALSE
      )
    }
    
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
        ci_high    = res$ci_high,
        zero_count = res$zero_count
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

# ── 6. Run pan-cancer ─────────────────────────────────────────────────────────
message("Running pan-cancer enrichment tests ...")

results_pancancer <- event_stat %>%
  group_by(gene) %>%
  group_modify(~ run_gene_tests(.x)) %>%
  ungroup() %>%
  mutate(IntoGen_cancer_type = "pan_cancer")

# ── 7. Combine and apply FDR correction ──────────────────────────────────────
results_all <- bind_rows(results_per_ct, results_pancancer)

# FDR applied only to non-NA p-values; zero_count rows retain NA
results_all <- results_all %>%
  group_by(IntoGen_cancer_type, test_type) %>%
  mutate(p_adj_local = p.adjust(p_value, method = FDR_METHOD)) %>%
  ungroup()

results_all <- results_all %>%
  group_by(test_type) %>%
  mutate(p_adj_global = p.adjust(p_value, method = FDR_METHOD)) %>%
  ungroup()

# ── 8. Count and proportion table ────────────────────────────────────────────
count_table <- bind_rows(
  event_stat %>% mutate(scope = "per_cancer_type"),
  event_stat %>% mutate(IntoGen_cancer_type = "pan_cancer", scope = "pan_cancer")
) %>%
  filter(present) %>%
  distinct(IntoGen_cancer_type, gene, sample, class, event, scope) %>%
  group_by(IntoGen_cancer_type, gene, class, event) %>%
  summarise(n = n(), .groups = "drop") %>%
  left_join(
    bind_rows(
      event_stat %>% mutate(scope = "per_cancer_type"),
      event_stat %>% mutate(IntoGen_cancer_type = "pan_cancer", scope = "pan_cancer")
    ) %>%
      distinct(IntoGen_cancer_type, gene, sample, class) %>%
      group_by(IntoGen_cancer_type, gene, class) %>%
      summarise(class_total = n(), .groups = "drop"),
    by = c("IntoGen_cancer_type", "gene", "class")
  ) %>%
  mutate(
    prop  = n / class_total,
    event = factor(event, levels = EVENT_LEVELS)
  )

# ── 9. HM-signature summary ───────────────────────────────────────────────────
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
    min_p_adj_local         = min(p_adj_local, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_cancer_types_enriched), event, desc(median_OR))

hm_pancancer_sig <- results_all %>%
  filter(
    IntoGen_cancer_type == "pan_cancer",
    test_type           == "per_event_HM_vs_rest",
    !is.na(odds_ratio),
    p_adj_local         < 0.05
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
print(hm_signature_summary, n = 30)

message("\n── CI_M in per_event_pairwise (verification) ──")
results_all %>%
  filter(test_type == "per_event_pairwise", event == "CI_M") %>%
  count(comparison, zero_count) %>%
  print()

message("\n── Zero-count rows per event and comparison (pan-cancer) ──")
results_all %>%
  filter(IntoGen_cancer_type == "pan_cancer",
         test_type == "per_event_pairwise",
         zero_count) %>%
  select(gene, event, comparison, n_WGD, n_HM, n_Classic) %>%
  print(n = 30)

message("\nOutput saved to: ", out_dir_stats)

