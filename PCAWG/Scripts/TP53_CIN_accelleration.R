rm(list = ls())
library(tidyverse)
library(ggplot2)

fmt_p <- function(test) {
  p <- test$p.value
  if (p < 0.001) "p < 0.001" else paste0("p = ", round(p, 3))
}

DRIVERS_PATH <- "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Data/Drivers.rds"
SEGMENTS_PATH <- "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Data/Segments.rds"
GENE <- "TP53"

# drivers_df %>% dplyr::group_by(gene) %>% 
#   dplyr::summarise(n = sum(!is.na(mult_estimate))) %>% 
#   dplyr::arrange(-n)

drivers_df <- readRDS(DRIVERS_PATH) %>%
  dplyr::filter(gene == GENE) %>%
  dplyr::select(segment_id, gene, karyotype, sample_id, mutation_status,
                mult_estimate, timed, clock_mean) %>%
  dplyr::rename(sample = sample_id)

segment_df <- readRDS(SEGMENTS_PATH) %>%
  dplyr::select(segment_name, karyotype, chr, clock_mean, sample, ttype)

segment_df %>% dplyr::sample_n(10)
drivers_df %>% dplyr::sample_n(10)

# ── SANITY CHECK ────────────────────────────────────────────────────────────
# FIX 1: Verify that segment clock_mean is truly per-CNA, not a sample summary.
# If every segment within a sample shares the same clock_mean, the variable
# is a sample-level average and cannot be used to split segments into
# "pre-TP53" vs "post-TP53" temporal bins. The analysis collapses entirely
# in that case — all segments would land on the same side of the split.
clock_variance_check <- segment_df %>%
  group_by(sample) %>%
  summarise(n_unique_clocks = n_distinct(clock_mean), .groups = "drop")

cat("Samples with only one unique clock_mean value:",
    sum(clock_variance_check$n_unique_clocks == 1),
    "out of", nrow(clock_variance_check), "\n")

# If the above shows that most samples have n_unique_clocks == 1,
# the segment clock cannot be used as a per-CNA timestamp.
# In that case you need a different measure of CIN timing 
# (e.g. segment-level posterior timing from tickTack, if available).
# The rest of this script assumes the check passes (per-CNA clocks exist).

## 1. Segment lengths --------------------------------------------------------
seg <- segment_df %>%
  separate(segment_name, c("nm_chr", "start", "end"),
           sep = "_", convert = TRUE, remove = FALSE) %>%
  mutate(length = end - start) %>%
  rename(seg_clock = clock_mean)

## 2. Timed TP53 mutant samples ----------------------------------------------
tp53 <- drivers_df %>%
  filter(mutation_status == "M",
         timed, !is.na(mult_estimate), !is.na(clock_mean),
         clock_mean > 0, clock_mean < 1) %>%
  transmute(
    sample,
    c_tp53 = clock_mean,
    mult   = round(mult_estimate)
  )

## 3. CIN rate (genome gained per unit molecular time) ----------------------
# FIX 2: Replace the fixed 1 Mb pseudocount with a sample-scaled one.
# A fixed pseudocount dominates samples with low total CIN burden,
# artificially shrinking their log_ratio toward 0.
# We use 1% of the sample's total gained length as the pseudocount instead.
cin_rate <- function(df, split) {
  total_len <- sum(df$length, na.rm = TRUE)
  pseudo    <- max(total_len * 0.01, 1e5)   # 1% of total, min 100 kb
  pre  <- sum(df$length[df$seg_clock <  split], na.rm = TRUE) + pseudo
  post <- sum(df$length[df$seg_clock >= split], na.rm = TRUE) + pseudo
  tibble(
    rate_pre      = pre  / split,
    rate_post     = post / (1 - split),
    log_ratio     = log2(rate_post / rate_pre),
    total_len_mb  = total_len / 1e6          # keep for matching / QC
  )
}

mut_rates <- seg %>%
  inner_join(tp53, by = "sample") %>%
  group_by(sample, ttype, mult, c_tp53) %>%
  group_modify(~ cin_rate(.x, .y$c_tp53)) %>%
  ungroup()

## 4. WT background: matched empirical control ─────────────────────────────
# FIX 3: Match WT samples on tumor type AND total CIN burden.
# The original code matched only on tumor type, but TP53-mutant tumors
# tend to have globally higher CIN. Without burden-matching, the
# post/pre ratio in WT samples reflects an intrinsically different CIN
# landscape, not just the absence of TP53 mutation.
set.seed(42)

# Compute total CIN burden for WT samples
wt_seg_burden <- seg %>%
  semi_join(
    drivers_df %>% filter(mutation_status == "WT") %>% distinct(sample),
    by = "sample"
  ) %>%
  group_by(sample) %>%
  summarise(total_len_mb = sum(length, na.rm = TRUE) / 1e6,
            ttype        = first(ttype),
            .groups      = "drop")

# Mutant pool with burden info, for matching
mut_pool <- mut_rates %>%
  select(sample, ttype, c_tp53, total_len_mb)

# For each WT sample: find the mutant in the same ttype whose burden is
# closest, then use that mutant's c_tp53 as the random split.
# This preserves the ttype-specific CIN timing distribution while
# controlling for the amount of total genomic instability.
find_closest_split <- function(wt_burden, pool_ttype) {
  if (nrow(pool_ttype) == 0) return(median(mut_pool$c_tp53, na.rm = TRUE))
  idx <- which.min(abs(pool_ttype$total_len_mb - wt_burden))
  pool_ttype$c_tp53[idx]
}

wt_splits <- wt_seg_burden %>%
  rowwise() %>%
  mutate(
    split = find_closest_split(
      total_len_mb,
      filter(mut_pool, ttype == .data$ttype)
    )
  ) %>%
  ungroup() %>%
  select(sample, ttype, split, total_len_mb)

wt_rates <- seg %>%
  inner_join(wt_splits, by = c("sample", "ttype")) %>%
  group_by(sample, ttype, split) %>%          # drop total_len_mb from grouping
  group_modify(~ cin_rate(.x, .y$split)) %>%
  ungroup() %>%
  # The group_modify recomputes total_len_mb inside cin_rate;
  # drop the duplicated column brought in from wt_splits
  select(-total_len_mb)

## 5. Statistical tests ──────────────────────────────────────────────────────

# FIX 4: Stratify the primary test by tumor type.
# Pooling all types confounds the comparison — OV is near-universally
# TP53-mutant with high CIN; luminal breast is the opposite. A simple
# Wilcoxon on the pooled data partly recovers ttype differences, not TP53.
# We use a linear model with ttype as a covariate instead.
res_primary <- bind_rows(
  mutate(mut_rates, grp = "mut"),
  mutate(wt_rates,  grp = "wt", mult = NA_integer_, c_tp53 = NA_real_)
)

# Tumor-type stratified comparison (ttype as fixed effect)
lm_primary <- lm(log_ratio ~ grp + ttype, data = res_primary)
summary(lm_primary)

# Keep the Wilcoxon as a non-parametric complement (it is robust to
# non-normality), but note it does not control for tumor type.
test_primary <- wilcox.test(log_ratio ~ grp, data = res_primary,
                            alternative = "greater")
cat("Wilcoxon (unadjusted):", fmt_p(test_primary), "\n")

# FIX 5: Correct the mult == 2 "conservative" reasoning.
# mult == 2 means TP53 was ALREADY PRESENT before the CN gain, so c_tp53
# is an UPPER BOUND on TP53 mutation time. The true TP53 event happened
# earlier, meaning the true "post-TP53" window is LONGER than measured.
# This makes our log_ratio UNDERESTIMATE the true acceleration — so finding
# log_ratio > 0 here is genuinely conservative: we'd expect an even larger
# ratio with the true (earlier) split.
# mult == 1 means TP53 happened AFTER the gain, so c_tp53 is a LOWER BOUND.
# In that case the post-TP53 window is shorter than measured, and log_ratio
# could be OVERESTIMATED. mult == 1 is therefore the anti-conservative case.
mut_conservative <- mut_rates %>% filter(mult == 2)   # genuinely conservative
mut_anticonserv  <- mut_rates %>% filter(mult == 1)   # lower bound: be cautious

test_conservative <- wilcox.test(mut_conservative$log_ratio, mu = 0,
                                 alternative = "greater")
test_anticonserv  <- wilcox.test(mut_anticonserv$log_ratio,  mu = 0,
                                 alternative = "greater")

cat("mult == 2 (conservative):", fmt_p(test_conservative), "\n")
cat("mult == 1 (anti-conservative):", fmt_p(test_anticonserv), "\n")

## 6. Plots ──────────────────────────────────────────────────────────────────

# --- Plot 1: Mutant vs WT CIN acceleration (tumor-type aware) -------------
comparison_df <- bind_rows(
  mutate(mut_rates, grp = "TP53 mutant"),
  mutate(wt_rates,  grp = "Wild-type (WT)")
)

p1_label <- data.frame(
  grp   = "TP53 mutant",
  label = fmt_p(test_primary),
  y     = max(comparison_df$log_ratio, na.rm = TRUE) * 0.95
)

p1 <- ggplot(comparison_df, aes(x = grp, y = log_ratio, fill = grp)) +
  geom_violin(alpha = 0.4, color = NA, adjust = 1.2) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8, linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.7) +
  geom_text(data = p1_label, aes(label = label, y = y),
            x = 1.5, vjust = -0.5, size = 3.5, color = "grey30",
            inherit.aes = FALSE) +
  # Facet by tumor type to show the stratified pattern
  facet_wrap(~ ttype, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c("TP53 mutant"    = "#d95f02",
                               "Wild-type (WT)" = "#7570b3")) +
  labs(
    title    = "CIN rate acceleration around TP53 mutation time",
    subtitle = "Stratified by tumor type; WT matched on ttype + CIN burden",
    x        = NULL,
    y        = expression(log[2] * "(rate"[post] * " / rate"[pre] * ")")
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position  = "bottom",
        plot.title       = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        strip.text       = element_text(size = 9))
p1

# --- Plot 2: Multiplicity — with correct labeling of conservatism ----------
mult_labels <- c(
  "1" = "mult = 1\nTP53 after gain\n(c_tp53: lower bound → may overestimate)",
  "2" = "mult = 2\nTP53 before gain\n(c_tp53: upper bound → conservative)"
)

p2_label <- data.frame(
  mult_bin = factor("2", levels = c("1", "2")),
  label    = fmt_p(test_conservative),
  y        = max(mut_rates$log_ratio, na.rm = TRUE) * 0.95
)

p2 <- mut_rates %>%
  mutate(mult_bin = factor(as.character(mult), levels = c("1", "2"))) %>%
  ggplot(aes(x = mult_bin, y = log_ratio, fill = mult_bin)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5, linewidth = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 1.2, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.7) +
  geom_text(data = p2_label, aes(x = mult_bin, y = y, label = label),
            vjust = -0.5, size = 3.5, color = "grey30", inherit.aes = FALSE) +
  scale_fill_brewer(palette = "Set2") +
  scale_x_discrete(labels = mult_labels) +
  labs(
    title    = "CIN acceleration in TP53 mutants by multiplicity",
    subtitle = paste(
      "mult = 2: upper bound on TP53 time → log2 ratio underestimates true acceleration (conservative).",
      "mult = 1: lower bound on TP53 time → log2 ratio may overestimate acceleration.",
      sep = "\n"
    ),
    x = "Mutation multiplicity (timing window)",
    y = expression(log[2] * "(rate"[post] * " / rate"[pre] * ")")
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position  = "none",
        plot.title       = element_text(face = "bold"),
        panel.grid.minor = element_blank())
p2

