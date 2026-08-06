rm(list = ls())
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
})

has_lme4 <- requireNamespace("lme4",  quietly = TRUE)
has_coin <- requireNamespace("coin",  quietly = TRUE)
has_broom <- requireNamespace("broom", quietly = TRUE)

# ── Section 0: Setup ──────────────────────────────────────────────────────────

# DRIVERS_PATH  <- "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Data/Drivers.rds"
# SEGMENTS_PATH <- "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Data/Segments.rds"
# PLOT_DIR      <- "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Plot/TP53_CIN/"
DRIVERS_PATH  <- "~/GitHub/material-tickTack-2026/PCAWG/Data/Drivers.rds"
SEGMENTS_PATH <- "~/GitHub/material-tickTack-2026/PCAWG/Data/Segments.rds"
PLOT_DIR      <- "~/GitHub/material-tickTack-2026/PCAWG/Plot/TP53_CIN/"

GENE      <- "TP53"
BURST_TOL <- 0      # n_unique_clocks <= BURST_TOL+1 → classified as burst
K_PERM    <- 200    # permutation replicates for the gradual-stratum null
HALDANE   <- 0.5   # continuity correction for count-based log-ratio

if (!dir.exists(PLOT_DIR)) dir.create(PLOT_DIR, recursive = TRUE)

fmt_p <- function(p) {
  if (inherits(p, "htest")) p <- p$p.value
  if (p < 0.001) "p < 0.001" else paste0("p = ", round(p, 3))
}

lr_count <- function(clocks, split) {
  pre  <- sum(clocks <  split)
  post <- sum(clocks >= split)
  log2(((post + HALDANE) / (1 - split)) / ((pre + HALDANE) / split))
}

run_mult_test <- function(df, m) {
  d <- filter(df, mult == m)
  if (nrow(d) < 5) {
    cat(sprintf("  mult == %d: only %d samples, skipping.\n", m, nrow(d)))
    return(invisible(NULL))
  }
  tt <- wilcox.test(d$delta, mu = 0, alternative = "greater")
  cat(sprintf("  mult == %d: %s  (n = %d, median delta = %.3f)\n",
              m, fmt_p(tt$p.value), nrow(d), median(d$delta)))
}

theme_nature <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(
      axis.line        = element_line(linewidth = 0.4, colour = "black"),
      axis.ticks       = element_line(linewidth = 0.3, colour = "black"),
      axis.text        = element_text(colour = "black", size = base_size - 1),
      axis.title       = element_text(colour = "black", size = base_size),
      strip.background = element_blank(),
      strip.text       = element_text(size = base_size - 1, colour = "black"),
      legend.text      = element_text(size = base_size - 1),
      legend.title     = element_text(size = base_size),
      legend.key.size  = unit(0.35, "cm"),
      plot.margin      = margin(4, 4, 4, 4)
    )
}

# ── Section 1: Load, define groups, classify ──────────────────────────────────

drivers_all <- readRDS(DRIVERS_PATH)

segment_df <- readRDS(SEGMENTS_PATH) %>%
  dplyr::select(segment_name, karyotype, chr, clock_mean, sample, ttype, class)

# All samples that appear in the timing data
all_samples <- segment_df %>% distinct(sample)

# TP53 entries in Drivers (non-diploid-WT by construction)
tp53_drivers <- drivers_all %>%
  filter(gene == GENE) %>%
  dplyr::select(segment_id, gene, karyotype, sample_id, mutation_status,
                mult_estimate, timed, clock_mean) %>%
  dplyr::rename(sample = sample_id) %>%
  mutate(mutation_status = if_else(mutation_status == "CNA_driver", "M", mutation_status))

# Mutant group: TP53 is M (or CNA_driver recoded to M)
mut_samples <- tp53_drivers %>%
  filter(mutation_status == "M") %>%
  distinct(sample)

# WT group:
#   - Primary:   absent from Drivers entirely → diploid WT for TP53
#   - Secondary: in Drivers with mutation_status == "WT" (on a CNA segment, not mutated)
#
# The Drivers table stores a gene only when it is NOT diploid WT, so absence = diploid WT.
wt_samples <- all_samples %>%
  left_join(tp53_drivers %>% distinct(sample, mutation_status), by = "sample") %>%
  filter(is.na(mutation_status) | mutation_status == "WT") %>%
  mutate(wt_type = if_else(is.na(mutation_status), "diploid_absent", "on_CNA_WT"))

# Per-sample segment summary
seg_summary <- segment_df %>%
  separate(segment_name, into = c("chr_dup", "start", "end"),
           sep = "_", convert = TRUE, remove = FALSE) %>%
  select(-chr_dup) %>%
  mutate(length = end - start) %>%
  group_by(sample) %>%
  summarise(
    n_unique_clocks = n_distinct(clock_mean),
    total_cin_mb    = sum(length, na.rm = TRUE) / 1e6,
    cin_clock       = if_else(n_distinct(clock_mean) == 1L,
                              first(clock_mean), NA_real_),
    ttype           = first(ttype),
    class           = first(class),
    .groups         = "drop"
  ) %>%
  mutate(is_burst = n_unique_clocks <= BURST_TOL + 1L)

# Full segment table with lengths for the gradual analysis
seg <- segment_df %>%
  separate(segment_name, into = c("chr_dup", "start", "end"),
           sep = "_", convert = TRUE, remove = FALSE) %>%
  select(-chr_dup) %>%
  mutate(length    = end - start,
         seg_clock = clock_mean)

cat("════════════════════════════════════════════════════════════════════════\n")
cat("SECTION 1 — Sample classification\n")
cat("════════════════════════════════════════════════════════════════════════\n")
cat(sprintf("Total samples in timing data:  %d\n", nrow(all_samples)))
cat(sprintf("TP53 mutant samples:           %d\n", nrow(mut_samples)))
cat(sprintf("TP53 WT samples (total):       %d\n", nrow(wt_samples)))
cat(sprintf("  diploid-absent:              %d\n", sum(wt_samples$wt_type == "diploid_absent")))
cat(sprintf("  on-CNA WT:                  %d\n", sum(wt_samples$wt_type == "on_CNA_WT")))
cat(sprintf("Burst samples:                 %d  (n_unique_clocks == 1)\n", sum(seg_summary$is_burst)))
cat(sprintf("Gradual samples:               %d  (n_unique_clocks  > 1)\n", sum(!seg_summary$is_burst)))

cat("\nMutant / WT in each stratum:\n")
seg_summary %>%
  mutate(grp = case_when(
    sample %in% mut_samples$sample ~ "mut",
    sample %in% wt_samples$sample  ~ "wt",
    TRUE ~ "other"
  )) %>%
  count(is_burst, grp) %>%
  pivot_wider(names_from = grp, values_from = n, values_fill = 0L) %>%
  print()

# ── Section 2 (PRIMARY): Multiplicity-based relative ordering ─────────────────
#
# mult_estimate directly encodes the temporal order of TP53 SNV vs the CNA gain:
#   mult == 2: SNV on both chromosome copies → predates the duplication
#              → TP53 was mutated BEFORE the CNA gain
#   mult == 1: SNV on one copy → postdates the duplication
#              → TP53 was mutated AFTER the CNA gain
#
# This test applies to both burst and gradual samples and requires only that
# the TP53 mutation is timed (timed == TRUE, mult_estimate not NA).

cat("\n════════════════════════════════════════════════════════════════════════\n")
cat("SECTION 2 (PRIMARY) — Multiplicity: does TP53 precede CIN events?\n")
cat("════════════════════════════════════════════════════════════════════════\n")

mult_df <- tp53_drivers %>%
  filter(mutation_status == "M", timed, !is.na(mult_estimate)) %>%
  group_by(sample) %>%
  summarise(mult = round(median(mult_estimate, na.rm = TRUE)), .groups = "drop") %>%
  left_join(select(seg_summary, sample, total_cin_mb, ttype, class, is_burst), by = "sample")

cat(sprintf("Timed TP53 mutant samples contributing to ordering test: %d\n", nrow(mult_df)))
cat(sprintf("  mult == 2 (TP53 precedes CIN): %d  (%.1f%%)\n",
            sum(mult_df$mult == 2), 100 * mean(mult_df$mult == 2)))
cat(sprintf("  mult == 1 (CIN precedes TP53): %d  (%.1f%%)\n",
            sum(mult_df$mult == 1), 100 * mean(mult_df$mult == 1)))

# 2.1 Global binomial test: is prop(mult == 2) > 0.5?
cat("\n── Test 2.1: Binomial test prop(mult == 2) > 0.5 ──\n")
binom_all <- binom.test(sum(mult_df$mult == 2), nrow(mult_df),
                        p = 0.5, alternative = "greater")
cat(fmt_p(binom_all), "\n")
cat(sprintf("95%% CI for prop(mult==2): [%.3f, %.3f]\n",
            binom_all$conf.int[1], binom_all$conf.int[2]))

# 2.2 Per tumor-type proportions and CIs
cat("\n── Test 2.2: Per-tumor-type ordering ──\n")
mult_by_ttype <- mult_df %>%
  group_by(ttype) %>%
  summarise(
    n           = n(),
    n_mult2     = sum(mult == 2),
    prop_mult2  = mean(mult == 2),
    ci_lo       = binom.test(sum(mult == 2), n())$conf.int[1],
    ci_hi       = binom.test(sum(mult == 2), n())$conf.int[2],
    p_binom     = binom.test(sum(mult == 2), n(), p = 0.5,
                             alternative = "greater")$p.value,
    .groups     = "drop"
  ) %>%
  arrange(desc(prop_mult2))
print(mult_by_ttype, n = Inf)

# 2.3 Does TP53-first (mult==2) associate with larger CIN burden?
cat("\n── Test 2.3: CIN burden by relative ordering ──\n")
lm_mult_burden <- lm(total_cin_mb ~ factor(mult) + ttype, data = mult_df)
coef_mult <- coef(summary(lm_mult_burden))["factor(mult)2", , drop = FALSE]
cat("LM total_cin_mb ~ mult + ttype   (mult==2 vs mult==1):\n")
print(coef_mult)
wt_mult_burden <- wilcox.test(total_cin_mb ~ factor(mult), data = mult_df, alternative = "greater")
cat("Wilcoxon (mult==2 burden > mult==1):", fmt_p(wt_mult_burden), "\n")

# Figure 1: PRIMARY multiplicity figure
clr_mult <- c("1" = "#66c2a5", "2" = "#fc8d62")

fig1a <- mult_by_ttype %>%
  filter(n >= 5) %>%
  mutate(ttype = fct_reorder(ttype, prop_mult2)) %>%
  ggplot(aes(x = prop_mult2, y = ttype)) +
  geom_col(fill = "#fc8d62", alpha = 0.8, width = 0.65) +
  geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi), width = 0.25, linewidth = 0.5) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  scale_x_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(x = "SNV before CNA gain (mult = 2)", y = NULL) +
  theme_nature()

fig1b <- mult_df %>%
  mutate(mult_f = factor(mult, levels = c("1", "2"),
                         labels = c("mult = 1", "mult = 2"))) %>%
  ggplot(aes(x = mult_f, y = total_cin_mb, fill = factor(mult))) +
  geom_violin(alpha = 0.4, color = NA, adjust = 1.2) +
  geom_boxplot(width = 0.22, outlier.shape = NA, alpha = 0.85, linewidth = 0.5) +
  annotate("text", x = 1.5,
           y = max(mult_df$total_cin_mb, na.rm = TRUE) * 0.97,
           label = fmt_p(wt_mult_burden), size = 2.8, colour = "grey30") +
  scale_fill_manual(values = clr_mult) +
  labs(x = NULL, y = "CIN burden (Mb)") +
  theme_nature()
  

fig1 <- fig1a + fig1b

ggsave(paste0(PLOT_DIR, "fig1_mult_ordering.pdf"),
       plot = fig1, width = 14, height = 6)
cat("\nFig 1 saved.\n")

# ── Section 3: Burst-stratum timing analyses ──────────────────────────────────

cat("\n════════════════════════════════════════════════════════════════════════\n")
cat("SECTION 3 — Burst-stratum timing (samples with a single CNA clock)\n")
cat("════════════════════════════════════════════════════════════════════════\n")

burst_summary <- seg_summary %>%
  filter(is_burst) %>%
  mutate(grp = case_when(
    sample %in% mut_samples$sample ~ "mut",
    sample %in% wt_samples$sample  ~ "wt",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(grp), !is.na(cin_clock))

cat(sprintf("Burst samples with defined cin_clock: %d  (mut=%d, wt=%d)\n",
            nrow(burst_summary),
            sum(burst_summary$grp == "mut"),
            sum(burst_summary$grp == "wt")))

# 3A: CIN burst timing by mutation status
cat("\n── Test 3A: CIN burst timing (mutant vs WT) ──\n")
lm3a <- lm(cin_clock ~ grp + ttype, data = burst_summary)
cat("LM cin_clock ~ grp + ttype:\n")
print(coef(summary(lm3a))["grpwt", , drop = FALSE])
wt3a <- wilcox.test(cin_clock ~ grp, data = burst_summary, alternative = "less")
cat("Wilcoxon (mut clock < wt clock):", fmt_p(wt3a), "\n")
cat("Medians:\n")
burst_summary %>% group_by(grp) %>%
  summarise(median_cin_clock = median(cin_clock), n = n(), .groups = "drop") %>% print()

# 3B: Does TP53 multiplicity predict when the CIN burst occurs?
#
# NOTE on validity: mult_estimate is derived from VAF and CN context independently
# of the timing clock, so comparing cin_clock between mult groups is NOT tautological
# (unlike the tp53_clock vs cin_clock scatter, which was circular).
#
# Among timed TP53-mutant burst samples:
#   mult == 2 → TP53 SNV predates the CNA gain → TP53 was lost BEFORE the burst
#   mult == 1 → TP53 SNV postdates the CNA gain → burst happened BEFORE TP53 was lost
#
# If TP53 loss enables or accelerates the burst: samples where TP53 was lost first
# (mult == 2) might have their burst at a different molecular time than samples where
# the burst preceded TP53 loss (mult == 1).

mult_burst <- mult_df %>%
  filter(is_burst, mult %in% c(1, 2)) %>%
  left_join(select(seg_summary, sample, cin_clock), by = "sample") %>%
  filter(!is.na(cin_clock))

cat(sprintf("\n── Test 3B: Does mult predict CIN burst timing? (n = %d timed burst mutants) ──\n",
            nrow(mult_burst)))
cat(sprintf("  mult == 2 (TP53 before burst): n = %d,  median cin_clock = %.3f\n",
            sum(mult_burst$mult == 2),
            median(mult_burst$cin_clock[mult_burst$mult == 2], na.rm = TRUE)))
cat(sprintf("  mult == 1 (burst before TP53): n = %d,  median cin_clock = %.3f\n",
            sum(mult_burst$mult == 1),
            median(mult_burst$cin_clock[mult_burst$mult == 1], na.rm = TRUE)))

lm3b <- lm(cin_clock ~ factor(mult) + ttype, data = mult_burst)
cat("LM cin_clock ~ factor(mult) + ttype  (effect of mult==2 vs mult==1):\n")
print(coef(summary(lm3b))["factor(mult)2", , drop = FALSE])
wt3b <- wilcox.test(cin_clock ~ factor(mult), data = mult_burst)
cat("Wilcoxon (two-sided):", fmt_p(wt3b), "\n")

# Figure 2: Burst-stratum — two panels
clr_grp  <- c("mut" = "#d95f02", "wt" = "#7570b3")
clr_mult <- c("1" = "#66c2a5", "2" = "#fc8d62")

fig2a <- burst_summary %>%
  mutate(grp_label = if_else(grp == "mut", "TP53 mut", "WT")) %>%
  ggplot(aes(x = grp_label, y = cin_clock, fill = grp)) +
  geom_violin(alpha = 0.4, color = NA, adjust = 1.2) +
  geom_boxplot(width = 0.22, outlier.shape = NA, alpha = 0.85, linewidth = 0.5) +
  annotate("text", x = 1.5,
           y = max(burst_summary$cin_clock, na.rm = TRUE) * 0.97,
           label = fmt_p(wt3a), size = 2.8, colour = "grey30") +
  scale_fill_manual(values = clr_grp) +
  labs(x = NULL, y = "CIN clock") +
  theme_nature()
  

fig2b <- mult_burst %>%
  mutate(mult_f = factor(mult, levels = c("1", "2"),
                         labels = c("mult = 1", "mult = 2"))) %>%
  ggplot(aes(x = mult_f, y = cin_clock, fill = mult_f)) +
  geom_violin(alpha = 0.4, color = NA, adjust = 1.2) +
  geom_boxplot(width = 0.22, outlier.shape = NA, alpha = 0.85, linewidth = 0.5) +
  geom_jitter(width = 0.1, alpha = 0.3, size = 1.0) +
  annotate("text", x = 1.5,
           y = max(mult_burst$cin_clock, na.rm = TRUE) * 0.97,
           label = fmt_p(wt3b), size = 2.8, colour = "grey30") +
  scale_fill_manual(values = clr_mult) +
  labs(x = NULL, y = "CIN clock") +
  theme_nature()
  

fig2 <- fig2a | fig2b

ggsave(paste0(PLOT_DIR, "fig2_burst_timing.pdf"),
       plot = fig2, width = 14, height = 6)
cat("\nFig 2 saved.\n")

# ── Section 4: Gradual-stratum timing analyses ────────────────────────────────

cat("\n════════════════════════════════════════════════════════════════════════\n")
cat("SECTION 4 — Gradual-stratum timing (samples with multiple CNA clocks)\n")
cat("════════════════════════════════════════════════════════════════════════\n")

# Timed TP53 mutants in gradual samples only
tp53_timed_grad <- tp53_drivers %>%
  filter(mutation_status == "M", timed, !is.na(mult_estimate),
         clock_mean > 0, clock_mean < 1) %>%
  semi_join(filter(seg_summary, !is_burst), by = "sample") %>%
  group_by(sample) %>%
  summarise(c_tp53 = mean(clock_mean),
            mult    = round(median(mult_estimate)),
            .groups = "drop")

n_grad_mut_total <- mut_samples %>%
  semi_join(filter(seg_summary, !is_burst), by = "sample") %>% nrow()
cat(sprintf("Gradual mutant samples:           %d\n", n_grad_mut_total))
cat(sprintf("Timed gradual mutant samples:     %d  (%.1f%% retained)\n",
            nrow(tp53_timed_grad),
            100 * nrow(tp53_timed_grad) / max(n_grad_mut_total, 1)))

# WT controls for gradual stratum: match on ttype + CIN burden
grad_wt_burden <- seg_summary %>%
  filter(!is_burst) %>%
  semi_join(wt_samples, by = "sample") %>%
  select(sample, total_cin_mb, ttype)

grad_mut_pool <- seg_summary %>%
  filter(!is_burst) %>%
  semi_join(tp53_timed_grad, by = "sample") %>%
  left_join(tp53_timed_grad, by = "sample") %>%
  select(sample, c_tp53, total_cin_mb, ttype)

set.seed(42)
wt_splits_grad <- grad_wt_burden %>%
  group_by(ttype) %>%
  group_modify(function(df, key) {
    pool <- dplyr::filter(grad_mut_pool, ttype == key$ttype)
    df$split <- if (nrow(pool) == 0) {
      median(grad_mut_pool$c_tp53, na.rm = TRUE)
    } else {
      pool$c_tp53[vapply(df$total_cin_mb,
                         \(b) which.min(abs(pool$total_cin_mb - b)),
                         integer(1))]
    }
    df
  }) %>%
  ungroup() %>%
  select(sample, ttype, split)

splits_all_grad <- bind_rows(
  transmute(tp53_timed_grad, sample, split = c_tp53, grp = "mut"),
  transmute(wt_splits_grad,  sample, split,          grp = "wt")
)

# Restrict seg to gradual samples
seg_grad <- seg %>% semi_join(filter(seg_summary, !is_burst), by = "sample")

per_sample_grad <- seg_grad %>%
  inner_join(splits_all_grad, by = "sample") %>%
  group_by(sample, ttype, grp, split) %>%
  summarise(
    pre_n  = sum(seg_clock <  first(split)),
    post_n = sum(seg_clock >= first(split)),
    .groups = "drop"
  ) %>%
  filter(split > 0, split < 1) %>%
  mutate(
    rate_pre  = (pre_n  + HALDANE) / split,
    rate_post = (post_n + HALDANE) / (1 - split),
    log_ratio = log2(rate_post / rate_pre)
  )

# 4A: Poisson GLMM
cat("\n── Test 4A: Poisson GLMM (event count, pre vs post TP53) ──\n")
long_grad <- per_sample_grad %>%
  pivot_longer(c(pre_n, post_n), names_to = "side", values_to = "events") %>%
  mutate(
    side     = factor(if_else(side == "post_n", "post", "pre"), levels = c("pre", "post")),
    grp      = factor(grp, levels = c("wt", "mut")),
    expo     = if_else(side == "post", 1 - split, split),
    log_expo = log(expo)
  )

if (has_lme4) {
  m4a <- lme4::glmer(
    events ~ side * grp + ttype + offset(log_expo) + (1 | sample),
    family  = poisson,
    data    = long_grad,
    control = lme4::glmerControl(optimizer = "bobyqa")
  )
  key_rows <- grep("^side|sidepost:grpmut", rownames(coef(summary(m4a))))
  print(coef(summary(m4a))[key_rows, , drop = FALSE])
  cat("Key effect: 'sidepost:grpmut'  (positive → TP53 mutants accumulate CNA\n",
      "            faster after TP53 loss than WT accumulates after matched split)\n")
  if (lme4::isSingular(m4a))
    message("Note: singular fit — random-effect variance collapsed to zero.")
} else {
  message("lme4 not installed → fixed-effects GLM fallback.")
  m4a <- glm(events ~ side * grp + ttype + offset(log_expo),
             family = poisson, data = long_grad)
  key_rows <- grep("^side|sidepost:grpmut", rownames(coef(summary(m4a))))
  print(coef(summary(m4a))[key_rows, , drop = FALSE])
}

# 4B: Stratified rank test (van Elteren)
cat("\n── Test 4B: Stratified rank test (blocked by ttype) ──\n")
informative_grad <- per_sample_grad %>%
  group_by(ttype) %>%
  filter(n_distinct(grp) == 2, n() >= 2) %>%
  ungroup() %>%
  mutate(grp   = factor(grp, levels = c("wt", "mut")),
         ttype = droplevels(factor(ttype)))

cat(sprintf("Informative tumor types: %d; samples: %d\n",
            n_distinct(informative_grad$ttype), nrow(informative_grad)))

if (has_coin && n_distinct(informative_grad$ttype) >= 1 && nrow(informative_grad) >= 4) {
  st4b <- coin::wilcox_test(log_ratio ~ grp | ttype,
                            data = informative_grad, alternative = "greater")
  cat("van Elteren (mut > wt):", fmt_p(coin::pvalue(st4b)), "\n")
} else {
  st4b <- wilcox.test(log_ratio ~ grp, data = per_sample_grad, alternative = "greater")
  cat("Wilcoxon (unadjusted, mut > wt):", fmt_p(st4b$p.value), "\n")
}

# 4C: Permutation null
cat("\n── Test 4C: Permutation null ──\n")
clocks_by_sample_grad <- seg_grad %>%
  filter(sample %in% tp53_timed_grad$sample) %>%
  group_by(sample) %>%
  summarise(clocks = list(seg_clock), .groups = "drop")

split_pool_grad <- tp53_timed_grad$c_tp53
set.seed(42)
null_tbl_grad <- tp53_timed_grad %>%
  left_join(clocks_by_sample_grad, by = "sample") %>%
  rowwise() %>%
  mutate(
    obs_lr    = lr_count(clocks, c_tp53),
    null_mean = mean(vapply(sample(split_pool_grad, K_PERM, replace = TRUE),
                            \(s) lr_count(clocks, s), numeric(1))),
    delta     = obs_lr - null_mean
  ) %>%
  ungroup() %>%
  select(sample, mult, obs_lr, null_mean, delta)

test_delta <- wilcox.test(null_tbl_grad$delta, mu = 0, alternative = "greater")
cat("All gradual timed mutants, delta > 0:", fmt_p(test_delta), "\n")
cat(sprintf("Median delta = %.3f\n", median(null_tbl_grad$delta, na.rm = TRUE)))

# 4D: Multiplicity substratification
cat("\n── Test 4D: Multiplicity substratification (gradual) ──\n")
cat("mult == 2 → conservative upper bound (underestimates true post window)\n")
cat("mult == 1 → point estimate (cautious)\n")
run_mult_test(null_tbl_grad, 2)
run_mult_test(null_tbl_grad, 1)

# Figure 3: Gradual-stratum
fig3a <- per_sample_grad %>%
  mutate(grp_label = if_else(grp == "mut", "TP53 mut", "WT")) %>%
  ggplot(aes(x = grp_label, y = log_ratio, fill = grp)) +
  geom_violin(alpha = 0.4, color = NA, adjust = 1.2) +
  geom_boxplot(width = 0.22, outlier.shape = NA, alpha = 0.85, linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  scale_fill_manual(values = clr_grp) +
  labs(x = NULL,
       y = expression(log[2] * "(rate"[post] * " / rate"[pre] * ")")) +
  theme_nature()
  

fig3b <- null_tbl_grad %>%
  pivot_longer(c(obs_lr, null_mean), names_to = "kind", values_to = "lr") %>%
  mutate(kind = recode(kind, obs_lr = "Observed", null_mean = "Null")) %>%
  ggplot(aes(x = kind, y = lr, fill = kind)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.7, linewidth = 0.5) +
  geom_jitter(width = 0.12, alpha = 0.25, size = 1.0) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  annotate("text", x = 1.5, y = Inf, vjust = 1.5,
           label = fmt_p(test_delta), size = 2.8, colour = "grey30") +
  scale_fill_brewer(palette = "Set2") +
  labs(x = NULL,
       y = expression(log[2] * "(rate"[post] * " / rate"[pre] * ")")) +
  theme_nature()
  

fig3c <- null_tbl_grad %>%
  filter(!is.na(mult)) %>%
  mutate(mult_f = factor(mult, levels = c("1", "2"),
                         labels = c("mult = 1", "mult = 2"))) %>%
  ggplot(aes(x = mult_f, y = delta, fill = factor(mult))) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.7, linewidth = 0.5) +
  geom_jitter(width = 0.12, alpha = 0.25, size = 1.0) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  scale_fill_manual(values = clr_mult) +
  labs(x = NULL, y = "Delta (obs - null mean)") +
  theme_nature()
  

fig3 <- fig3a | fig3b | fig3c

ggsave(paste0(PLOT_DIR, "fig3_gradual_acceleration.pdf"),
       plot = fig3, width = 18, height = 6)
cat("\nFig 3 saved.\n")

# ── Section 5 (Supplementary): CIN burden comparison, timing-free ─────────────

cat("\n════════════════════════════════════════════════════════════════════════\n")
cat("SECTION 5 (Supplementary) — CIN burden: all mutants vs all WT\n")
cat("════════════════════════════════════════════════════════════════════════\n")

burden_df <- seg_summary %>%
  mutate(grp = case_when(
    sample %in% mut_samples$sample ~ "mut",
    sample %in% wt_samples$sample  ~ "wt",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(grp))

cat(sprintf("Samples: mut = %d, wt = %d\n",
            sum(burden_df$grp == "mut"), sum(burden_df$grp == "wt")))
cat("(WT includes all diploid-absent samples — the expanded WT group)\n")

lm5 <- lm(total_cin_mb ~ grp + ttype, data = burden_df)
cat("\nLM total_cin_mb ~ grp + ttype:\n")
print(coef(summary(lm5))["grpwt", , drop = FALSE])

wt5 <- wilcox.test(total_cin_mb ~ grp, data = burden_df, alternative = "greater")
cat("Wilcoxon (mut burden > wt burden):", fmt_p(wt5), "\n")

burden_df %>% group_by(grp) %>%
  summarise(median_mb = median(total_cin_mb), n = n(), .groups = "drop") %>% print()

# Per-ttype Wilcoxon tests for figure annotation
burden_pvals <- burden_df %>%
  group_by(ttype) %>%
  filter(n_distinct(grp) == 2, n() >= 4) %>%
  summarise(
    p_val = tryCatch(
      wilcox.test(total_cin_mb ~ grp, alternative = "greater")$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(p_label = sapply(p_val, fmt_p))

fig4 <- burden_df %>%
  mutate(grp_label = if_else(grp == "mut", "TP53 mut", "WT")) %>%
  left_join(burden_pvals, by = "ttype") %>%
  ggplot(aes(x = grp_label, y = total_cin_mb, fill = grp)) +
  geom_violin(alpha = 0.4, color = NA, adjust = 1.2) +
  geom_boxplot(width = 0.22, outlier.shape = NA, alpha = 0.85, linewidth = 0.5) +
  geom_text(
    data = burden_pvals,
    aes(x = 1.5, y = Inf, label = p_label),
    vjust = 1.3, size = 2.2, colour = "grey30", inherit.aes = FALSE
  ) +
  facet_wrap(~ ttype, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = clr_grp) +
  labs(x = NULL, y = "CIN burden (Mb)") +
  theme_nature(base_size = 9) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust = 1))

ggsave(paste0(PLOT_DIR, "fig4_burden_supplementary.pdf"),
       plot = fig4, width = 18, height = 12)
cat("\nFig 4 saved.\n")

# ── Section 6: Temporal ordering landscape — CI_M genes and pan-driver mult ───
#
# Three new analyses given that Test 3B (mult predicting burst timing) was not
# significant:
#   6A. Pan-driver multiplicity landscape: extend Section 2 to all M genes.
#   6B. clock_rank of CI_M vs M drivers in gradual samples (CI_M lacks mult_estimate).
#   6C. Within-sample TP53 vs CI_M temporal ordering via clock_mean comparison.

cat("\n════════════════════════════════════════════════════════════════════════\n")
cat("SECTION 6A — Pan-driver multiplicity landscape (all M genes)\n")
cat("════════════════════════════════════════════════════════════════════════\n")

all_M_mult_raw <- drivers_all %>%
  dplyr::select(-sample) %>%
  dplyr::rename(sample = sample_id) %>%
  filter(mutation_status == "M", timed, !is.na(mult_estimate)) %>%
  group_by(gene, sample) %>%
  summarise(mult = round(median(mult_estimate, na.rm = TRUE)), .groups = "drop")

all_M_mult <- all_M_mult_raw %>%
  group_by(gene) %>%
  filter(n() >= 10) %>%
  summarise(
    n          = n(),
    n_early    = sum(mult == 2),
    prop_early = n_early / n,
    lo         = binom.test(n_early, n)$conf.int[1],
    hi         = binom.test(n_early, n)$conf.int[2],
    p_val      = binom.test(n_early, n, p = 0.5, alternative = "greater")$p.value,
    .groups    = "drop"
  ) %>%
  mutate(fdr = p.adjust(p_val, method = "BH")) %>%
  arrange(desc(prop_early))

cat(sprintf("Genes passing n >= 10 timed samples: %d\n", nrow(all_M_mult)))
print(all_M_mult, n = Inf)

fig6a <- all_M_mult %>%
  mutate(
    gene    = fct_reorder(gene, prop_early),
    sig     = case_when(fdr < 0.1 ~ "FDR < 0.1", p_val < 0.05 ~ "p < 0.05", TRUE ~ "NS"),
    is_tp53 = gene == "TP53"
  ) %>%
  ggplot(aes(x = prop_early, y = gene, colour = sig)) +
  geom_segment(aes(x = lo, xend = hi, yend = gene), linewidth = 0.5) +
  geom_point(aes(size = n, shape = is_tp53)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  scale_colour_manual(
    values = c("FDR < 0.1" = "#e41a1c", "p < 0.05" = "#ff7f00", "NS" = "grey60"),
    name = NULL
  ) +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 18), guide = "none") +
  scale_size_continuous(range = c(2, 5), name = "n") +
  scale_x_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(x = "SNV before CNA gain (mult = 2)", y = NULL) +
  theme_nature() +
  theme(axis.text.y = element_text(size = 8))

ggsave(paste0(PLOT_DIR, "fig6a_pan_driver_mult.pdf"),
       plot = fig6a, width = 10, height = max(5, nrow(all_M_mult) * 0.35 + 2))
cat("\nFig 6A saved.\n")

cat("\n════════════════════════════════════════════════════════════════════════\n")
cat("SECTION 6B — clock_rank of CI_M vs driver SNVs in gradual samples\n")
cat("════════════════════════════════════════════════════════════════════════\n")

gradual_samples <- seg_summary %>% filter(!is_burst) %>% pull(sample)

# Compute per-sample clock_rank from segment timing values
sample_clock_ranks <- segment_df %>%
  filter(sample %in% gradual_samples) %>%
  group_by(sample) %>%
  mutate(clock_rank = dense_rank(clock_mean),
         n_clocks   = n_distinct(clock_mean)) %>%
  ungroup() %>%
  distinct(sample, clock_mean, clock_rank, n_clocks)

driver_rank_df <- drivers_all %>%
  dplyr::select(-sample) %>%
  dplyr::rename(sample = sample_id) %>%
  filter(
    mutation_status %in% c("M", "CI_M"),
    timed,
    sample %in% gradual_samples
  ) %>%
  left_join(sample_clock_ranks, by = c("sample", "clock_mean")) %>%
  filter(!is.na(clock_rank)) %>%
  mutate(rank_cat = if_else(clock_rank == 1L, "earliest", "later"))

cat(sprintf("Gradual samples: %d\n", length(gradual_samples)))
cat(sprintf("  M entries with clock_rank:   %d (unique samples: %d)\n",
            sum(driver_rank_df$mutation_status == "M"),
            n_distinct(driver_rank_df$sample[driver_rank_df$mutation_status == "M"])))
cat(sprintf("  CI_M entries with clock_rank: %d (unique samples: %d)\n",
            sum(driver_rank_df$mutation_status == "CI_M"),
            n_distinct(driver_rank_df$sample[driver_rank_df$mutation_status == "CI_M"])))

rank_tab <- driver_rank_df %>%
  count(mutation_status, rank_cat) %>%
  group_by(mutation_status) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()
print(rank_tab)

# Fisher test: is prop(earliest) higher for M than CI_M?
fish_tab <- rank_tab %>%
  select(mutation_status, rank_cat, n) %>%
  pivot_wider(names_from = rank_cat, values_from = n, values_fill = 0L) %>%
  tibble::column_to_rownames("mutation_status") %>%
  as.matrix()

ft6b       <- NULL
ft6b_label <- "NA"
cat("\n── Test 6B.1: Fisher test prop(earliest) M vs CI_M ──\n")
if (all(dim(fish_tab) == c(2, 2))) {
  ft6b       <- fisher.test(fish_tab)
  ft6b_label <- fmt_p(ft6b$p.value)
  cat(ft6b_label, "\n")
  cat(sprintf("OR = %.2f  [%.2f, %.2f]\n",
              ft6b$estimate, ft6b$conf.int[1], ft6b$conf.int[2]))
} else {
  cat("Contingency table not 2x2 — check levels.\n")
}

cat("\n── Test 6B.2: Wilcoxon clock_rank distribution M vs CI_M ──\n")
wt6b <- wilcox.test(clock_rank ~ mutation_status, data = driver_rank_df)
cat(fmt_p(wt6b), "\n")

fig6b <- rank_tab %>%
  mutate(mutation_status = factor(mutation_status, levels = c("M", "CI_M"))) %>%
  ggplot(aes(x = mutation_status, y = prop, fill = rank_cat)) +
  geom_col(width = 0.55, colour = "white", linewidth = 0.3) +
  annotate("text", x = 1.5, y = Inf, vjust = 1.5,
           label = ft6b_label, size = 2.8, colour = "grey30") +
  scale_fill_manual(values = c("earliest" = "#4292c6", "later" = "#d9d9d9"), name = NULL) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "Proportion") +
  theme_nature()

ggsave(paste0(PLOT_DIR, "fig6b_clock_rank_M_vs_CIM.pdf"),
       plot = fig6b, width = 6, height = 5)
cat("\nFig 6B saved.\n")

cat("\n════════════════════════════════════════════════════════════════════════\n")
cat("SECTION 6C — Within-sample ordering: TP53 vs CI_M (clock_mean comparison)\n")
cat("════════════════════════════════════════════════════════════════════════\n")

tp53_clocks6c <- drivers_all %>%
  dplyr::select(-sample) %>%
  dplyr::rename(sample = sample_id) %>%
  filter(gene == "TP53", mutation_status == "M", timed) %>%
  group_by(sample) %>%
  summarise(tp53_clock = mean(clock_mean, na.rm = TRUE), .groups = "drop")

cim_clocks6c <- drivers_all %>%
  dplyr::select(-sample) %>%
  dplyr::rename(sample = sample_id) %>%
  filter(mutation_status == "CI_M", timed) %>%
  group_by(sample) %>%
  summarise(
    cim_clock_min = min(clock_mean, na.rm = TRUE),
    cim_clock_med = median(clock_mean, na.rm = TRUE),
    n_cim         = n(),
    .groups       = "drop"
  )

both6c <- inner_join(tp53_clocks6c, cim_clocks6c, by = "sample") %>%
  left_join(select(seg_summary, sample, ttype, is_burst), by = "sample")

cat(sprintf("Samples with both timed TP53 and ≥1 timed CI_M mutation: %d\n", nrow(both6c)))
cat(sprintf("  burst: %d  |  gradual: %d\n",
            sum(both6c$is_burst, na.rm = TRUE), sum(!both6c$is_burst, na.rm = TRUE)))

prop_tp53_first <- mean(both6c$tp53_clock < both6c$cim_clock_min, na.rm = TRUE)
cat(sprintf("Prop(TP53 precedes earliest CI_M): %.1f%% (%d / %d samples)\n",
            100 * prop_tp53_first,
            sum(both6c$tp53_clock < both6c$cim_clock_min, na.rm = TRUE),
            sum(!is.na(both6c$tp53_clock) & !is.na(both6c$cim_clock_min))))

wt6c       <- NULL
wt6c_label <- "n < 5"
cat("\n── Test 6C.1: Wilcoxon signed-rank (paired) TP53 clock < CI_M earliest clock ──\n")
if (nrow(both6c) >= 5) {
  wt6c       <- wilcox.test(both6c$tp53_clock, both6c$cim_clock_min,
                            paired = TRUE, alternative = "less")
  wt6c_label <- fmt_p(wt6c)
  cat(wt6c_label, "\n")
} else {
  cat("Too few paired samples (< 5) — test not run.\n")
}

cat("\n── Test 6C.2: Gradual-only (removes burst where all clocks collapse) ──\n")
both6c_gradual <- both6c %>% filter(!is_burst)
cat(sprintf("Gradual samples only: %d\n", nrow(both6c_gradual)))
if (nrow(both6c_gradual) >= 5) {
  wt6c_gradual <- wilcox.test(both6c_gradual$tp53_clock, both6c_gradual$cim_clock_min,
                               paired = TRUE, alternative = "less")
  cat(fmt_p(wt6c_gradual), "\n")
  cat(sprintf("Prop(TP53 first) in gradual: %.1f%%\n",
              100 * mean(both6c_gradual$tp53_clock < both6c_gradual$cim_clock_min, na.rm = TRUE)))
} else {
  cat("Too few gradual paired samples — test not run.\n")
}

p6c_label <- wt6c_label
fig6c <- both6c %>%
  ggplot(aes(x = tp53_clock, y = cim_clock_min, colour = ttype, shape = is_burst)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_point(size = 2.0, alpha = 0.75) +
  scale_shape_manual(values = c("TRUE" = 17, "FALSE" = 16),
                     labels = c("TRUE" = "burst", "FALSE" = "gradual"),
                     name = NULL) +
  annotate("text", x = 0.05, y = 0.94,
           label = sprintf("n = %d; %s", nrow(both6c), p6c_label),
           hjust = 0, size = 2.8, colour = "grey30") +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  labs(x = "TP53 clock", y = "Earliest CI_M clock", colour = NULL) +
  theme_nature()

ggsave(paste0(PLOT_DIR, "fig6c_tp53_vs_cim_ordering.pdf"),
       plot = fig6c, width = 7, height = 6)
cat("\nFig 6C saved.\n")

cat("\n════════════════════════════════════════════════════════════════════════\n")
cat("Done. Plots written to:", PLOT_DIR, "\n")
cat("════════════════════════════════════════════════════════════════════════\n")
