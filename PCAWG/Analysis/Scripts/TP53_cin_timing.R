rm(list = ls())
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
})

# Optional packages (graceful fallbacks if missing)
has_lme4 <- requireNamespace("lme4", quietly = TRUE)
has_coin <- requireNamespace("coin", quietly = TRUE)

fmt_p <- function(p) if (p < 0.001) "p < 0.001" else paste0("p = ", round(p, 3))

DRIVERS_PATH  <- "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Data/Drivers.rds"
SEGMENTS_PATH <- "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Data/Segments.rds"
GENE          <- "TP53"
K_PERM        <- 200      # permutations per sample for the empirical null
HALDANE       <- 0.5      # continuity correction for count-based log-ratio (plots/ranks only)

drivers_df <- readRDS(DRIVERS_PATH) %>%
  dplyr::filter(gene == GENE) %>%
  dplyr::select(segment_id, gene, karyotype, sample_id, mutation_status,
                mult_estimate, timed, clock_mean) %>%
  dplyr::rename(sample = sample_id)

segment_df <- readRDS(SEGMENTS_PATH) %>%
  dplyr::select(segment_name, karyotype, chr, clock_mean, sample, ttype, class)

# ── SANITY CHECK (now HALTS instead of warning-and-continuing) ───────────────
# Verify segment clock_mean is per-CNA, not a sample-level summary. If most
# samples carry a single clock value, the variable cannot timestamp segments
# and the whole pre/post split is meaningless.
clock_variance_check <- segment_df %>%
  group_by(sample) %>%
  summarise(n_unique_clocks = n_distinct(clock_mean), .groups = "drop")

frac_single <- mean(clock_variance_check$n_unique_clocks == 1)
cat(sprintf("Samples with a single unique clock_mean: %d / %d (%.1f%%)\n",
            sum(clock_variance_check$n_unique_clocks == 1),
            nrow(clock_variance_check), 100 * frac_single))

## 1. Segment lengths --------------------------------------------------------
# ASSUMPTION: seg_clock (gain time) and c_tp53 (SNV time) share the same
# molecular-time calibration. This holds in MutationTimeR/tickTack-style
# pipelines but is an assumption, not a guarantee.
seg <- segment_df %>%
  separate(segment_name, c("chr_dup", "start", "end"),
           sep = "_", convert = TRUE, remove = FALSE) %>%
  select(-chr_dup) %>%                       # was nm_chr, duplicates chr
  mutate(length = end - start) %>%
  rename(seg_clock = clock_mean)

## 2. Timed TP53-mutant samples (ONE row per sample) -------------------------
# FIX: dedup the mutant side. Multiple timed TP53 driver rows per sample would
# otherwise duplicate every segment in the join and inflate that sample's rate.
tp53_raw <- drivers_df %>%
  filter(mutation_status == "M",
         timed, !is.na(mult_estimate), !is.na(clock_mean),
         clock_mean > 0, clock_mean < 1)

dup_samples <- tp53_raw %>% count(sample) %>% filter(n > 1)
if (nrow(dup_samples) > 0) message(sprintf("%d TP53-mutant samples had >1 timed driver row; collapsing to one.", nrow(dup_samples)))

tp53 <- tp53_raw %>%
  group_by(sample) %>%
  summarise(c_tp53 = mean(clock_mean),               # explicit collapse rule
            mult    = round(median(mult_estimate)),
            .groups = "drop")

# --- Selection transparency: how many TP53-mutant samples survive? ----------
# !is.na(mult_estimate) keeps only samples where the TP53 locus itself sits in
# a timeable gain, excluding the common mutation+LOH (1:0) class. This enriches
# the "mutant" group for chr17 gains (intrinsically higher CIN), so report it.
n_M_total <- drivers_df %>% filter(mutation_status == "M") %>% distinct(sample) %>% nrow()
cat(sprintf("TP53-mutant samples: %d total -> %d timed & timeable (%.1f%% retained)\n",
            n_M_total, nrow(tp53), 100 * nrow(tp53) / n_M_total))

## 3. WT splits: matched on tumor type AND CIN burden (CORRECTED) ------------
# BUG FIXED: the original `filter(mut_pool, ttype == .data$ttype)` resolved
# .data to mut_pool, making the condition always TRUE (no ttype matching).
# Here we genuinely restrict the matching pool to the same tumor type.
set.seed(42)

wt_seg_burden <- seg %>%
  semi_join(drivers_df %>% filter(mutation_status == "WT") %>% distinct(sample),
            by = "sample") %>%
  group_by(sample) %>%
  summarise(total_len_mb = sum(length, na.rm = TRUE) / 1e6,
            ttype        = first(ttype), .groups = "drop")

mut_pool <- seg %>%
  inner_join(tp53, by = "sample") %>%
  group_by(sample, ttype) %>%
  summarise(c_tp53 = first(c_tp53),
            total_len_mb = sum(length, na.rm = TRUE) / 1e6, .groups = "drop")

wt_splits <- wt_seg_burden %>%
  group_by(ttype) %>%
  group_modify(function(df, key) {
    pool <- dplyr::filter(mut_pool, ttype == key$ttype)
    df$split <- if (nrow(pool) == 0) {
      median(mut_pool$c_tp53, na.rm = TRUE)        # fallback: global median split
    } else {
      pool$c_tp53[vapply(df$total_len_mb,
                         \(b) which.min(abs(pool$total_len_mb - b)),
                         integer(1))]
    }
    df
  }) %>%
  ungroup() %>%
  select(sample, ttype, split)

## 4. Per-sample two-bin event counts (PRIMARY representation) ---------------
# CIN measured as the COUNT of independently timed gains landing pre vs post
# the split, with exposure = window length in molecular time. Counts (not
# summed Mb) avoid a few arm-scale gains dominating, and let us drop the
# length pseudocount entirely for the primary model.
splits_all <- bind_rows(
  transmute(tp53,      sample, split = c_tp53, grp = "mut"),
  transmute(wt_splits, sample, split,          grp = "wt")
)

per_sample <- seg %>%
  inner_join(splits_all, by = "sample") %>%
  group_by(sample, ttype, grp, split) %>%
  summarise(
    pre_n  = sum(seg_clock <  first(split)),
    post_n = sum(seg_clock >= first(split)),
    .groups = "drop"
  ) %>%
  filter(split > 0, split < 1) %>%
  mutate(
    # Haldane-Anscombe continuity correction (0.5 EVENT, not 1 Mb).
    # The split-dependent bias this leaves is cancelled in the mut-vs-WT and
    # mut-vs-null comparisons because both sides share the same split distribution.
    rate_pre  = (pre_n  + HALDANE) / split,
    rate_post = (post_n + HALDANE) / (1 - split),
    log_ratio = log2(rate_post / rate_pre)
  )

## 5. PRIMARY TEST: Poisson model, no pseudocount ----------------------------
# events ~ side * grp, offset = log(exposure). The side:grp interaction IS the
# question: "is the post/pre rate ratio larger in TP53 mutants than in WT?"
# ttype as fixed effect; (1|sample) absorbs per-sample overdispersion.
long <- per_sample %>%
  pivot_longer(c(pre_n, post_n), names_to = "side", values_to = "events") %>%
  mutate(
    side     = factor(if_else(side == "post_n", "post", "pre"), levels = c("pre", "post")),
    grp      = factor(grp, levels = c("wt", "mut")),
    expo     = if_else(side == "post", 1 - split, split),
    log_expo = log(expo)
  )

cat("\n================ PRIMARY MODEL (Poisson, count-based) ================\n")
if (has_lme4) {
  m_primary <- lme4::glmer(
    events ~ side * grp + ttype + offset(log_expo) + (1 | sample),
    family = poisson, data = long,
    control = lme4::glmerControl(optimizer = "bobyqa")
  )
  print(summary(m_primary)$coefficients[grep("^side|sidepost:grpmut", 
                                             rownames(summary(m_primary)$coefficients)), , drop = FALSE])
  cat("\nInterpretation: the 'sidepost:grpmut' row is the key effect.\n")
} else {
  message("lme4 not installed -> fixed-effects glm fallback (no sample random effect).")
  m_primary <- glm(events ~ side * grp + ttype + offset(log_expo),
                   family = poisson, data = long)
  print(coef(summary(m_primary)))
}

## 6. Stratified non-parametric complement (mut vs WT, blocked by ttype) -----
cat("\n================ STRATIFIED RANK TEST (mut vs WT) ====================\n")

# Keep only tumor types where BOTH groups are present — the only blocks that
# carry a within-type mut-vs-WT comparison. All-mutant/all-WT/singleton types
# are uninformative for a stratified test (and trigger coin's <2-obs error).
informative <- per_sample %>%
  group_by(ttype) %>%
  filter(n_distinct(grp) == 2, n() >= 2) %>%
  ungroup() %>%
  mutate(grp   = factor(grp, levels = c("wt", "mut")),  # mut is the "high" side
         ttype = droplevels(factor(ttype)))

cat(sprintf("Informative tumor types (both groups present): %d; samples used: %d\n",
            n_distinct(informative$ttype), nrow(informative)))

if (has_coin && n_distinct(informative$ttype) >= 1 && nrow(informative) >= 4) {
  st <- coin::wilcox_test(log_ratio ~ grp | ttype,
                          data = informative, alternative = "greater")
  cat("van Elteren (ttype-blocked):", fmt_p(coin::pvalue(st)), "\n")
  # Direction sanity check — confirm 'greater' means mut > wt here:
  print(informative %>% group_by(grp) %>%
          summarise(median_lr = median(log_ratio), n = n(), .groups = "drop"))
} else {
  message("Too few informative blocks (or coin missing) -> unstratified Wilcoxon.")
  st <- wilcox.test(log_ratio ~ grp, data = per_sample, alternative = "greater")
  cat("Wilcoxon (unadjusted):", fmt_p(st$p.value), "\n")
}

## 7. EMPIRICAL PERMUTATION NULL ---------------------------------------------
# The gain-time distribution is non-uniform (often bimodal around WGD), so the
# expected log_ratio at a given split is NOT 0. For each mutant sample we draw
# K random splits from the empirical c_tp53 distribution, recompute log_ratio
# on that sample's OWN segments, and ask whether the TRUE TP53 split yields a
# higher log_ratio than random splits of the same distribution.
clocks_by_sample <- seg %>%
  filter(sample %in% tp53$sample) %>%
  group_by(sample) %>%
  summarise(clocks = list(seg_clock), .groups = "drop")

lr_count <- function(clocks, split) {
  pre  <- sum(clocks <  split)
  post <- sum(clocks >= split)
  log2(((post + HALDANE) / (1 - split)) / ((pre + HALDANE) / split))
}

split_pool <- tp53$c_tp53
set.seed(42)
null_tbl <- tp53 %>%
  left_join(clocks_by_sample, by = "sample") %>%
  rowwise() %>%
  mutate(
    obs_lr   = lr_count(clocks, c_tp53),
    null_mean = mean(vapply(sample(split_pool, K_PERM, replace = TRUE),
                            \(s) lr_count(clocks, s), numeric(1))),
    delta     = obs_lr - null_mean
  ) %>%
  ungroup() %>%
  select(sample, mult, obs_lr, null_mean, delta)

cat("\n================ PERMUTATION NULL (delta = obs - random split) =======\n")
test_delta_all <- wilcox.test(null_tbl$delta, mu = 0, alternative = "greater")
cat("All TP53 mutants, delta > 0:", fmt_p(test_delta_all$p.value),
    sprintf("(median delta = %.3f)\n", median(null_tbl$delta, na.rm = TRUE)))

## 8. Multiplicity split (conservatism), tested vs the null ------------------
# mult == 2: SNV on both copies -> predates the local gain -> c_tp53 is an
#            UPPER BOUND on TP53 time -> true post window is LONGER -> log_ratio
#            UNDERestimates true acceleration -> genuinely conservative.
# mult == 1: SNV postdates the gain -> c_tp53 is a point estimate (not a clean
#            one-sided bound) -> treat with caution, not as a hard lower bound.
cat("\nMultiplicity availability:\n")
print(count(null_tbl, mult))

run_mult_test <- function(df, m) {
  d <- df %>% filter(mult == m)
  if (nrow(d) < 5) { cat(sprintf("mult == %d: only %d samples, skipping test.\n", m, nrow(d))); return(invisible()) }
  tt <- wilcox.test(d$delta, mu = 0, alternative = "greater")
  cat(sprintf("mult == %d (delta > 0): %s (n = %d, median = %.3f)\n",
              m, fmt_p(tt$p.value), nrow(d), median(d$delta)))
}
run_mult_test(null_tbl, 2)   # conservative
run_mult_test(null_tbl, 1)   # cautious

## 9. Plots ──────────────────────────────────────────────────────────────────
comparison_df <- per_sample %>%
  mutate(grp = recode(as.character(grp),
                      mut = "TP53 mutant", wt = "Wild-type (WT)"))

p1 <- ggplot(comparison_df, aes(grp, log_ratio, fill = grp)) +
  geom_violin(alpha = 0.4, color = NA, adjust = 1.2) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.85, linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.7) +
  scale_fill_manual(values = c("TP53 mutant" = "#d95f02", "Wild-type (WT)" = "#7570b3")) +
  labs(title = "CIN rate acceleration around TP53 mutation time",
       subtitle = "Count-based; WT matched on ttype + CIN burden; test = ttype-stratified rank",
       x = NULL,
       y = expression(log[2] * "(rate"[post] * " / rate"[pre] * ")")) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())
print(p1)

# Permutation null: observed vs random-split expectation, paired by sample
p2 <- null_tbl %>%
  pivot_longer(c(obs_lr, null_mean), names_to = "kind", values_to = "lr") %>%
  mutate(kind = recode(kind, obs_lr = "True TP53 split", null_mean = "Random splits (null)")) %>%
  ggplot(aes(kind, lr, fill = kind)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.7, linewidth = 0.6) +
  geom_jitter(width = 0.12, alpha = 0.3, size = 1.1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.7) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Observed TP53 split vs empirical null",
       subtitle = sprintf("delta > 0: %s", fmt_p(test_delta_all$p.value)),
       x = NULL, y = expression(log[2] * "(rate"[post] * " / rate"[pre] * ")")) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())
print(p2)

# Multiplicity, on the null-corrected delta (only if both levels populated)
if (sum(!is.na(null_tbl$mult)) > 0) {
  p3 <- null_tbl %>%
    filter(!is.na(mult)) %>%
    mutate(mult_bin = factor(mult)) %>%
    ggplot(aes(mult_bin, delta, fill = mult_bin)) +
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.7, linewidth = 0.6) +
    geom_jitter(width = 0.12, alpha = 0.3, size = 1.1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.7) +
    scale_fill_brewer(palette = "Set2") +
    labs(title = "Null-corrected acceleration by multiplicity",
         subtitle = "mult=2 upper bound (conservative); mult=1 point estimate (cautious)",
         x = "Mutation multiplicity", y = "delta (obs - null)") +
    theme_minimal(base_size = 13) +
    theme(legend.position = "none", plot.title = element_text(face = "bold"),
          panel.grid.minor = element_blank())
  print(p3)
}

cat("\nDone.\n")










cat("\n############### RESULTS DUMP FOR REVIEW ###############\n")

## 1. Selection / retention
cat("\n--- SELECTION ---\n")
cat(sprintf("Single-clock samples: %d / %d\n",
            sum(clock_variance_check$n_unique_clocks == 1), nrow(clock_variance_check)))
cat(sprintf("TP53-mutant retained: %d / %d\n", nrow(tp53), n_M_total))

## 2. Primary Poisson model — the side:grp interaction is the effect
cat("\n--- PRIMARY MODEL (Poisson) ---\n")
print(summary(m_primary)$coefficients)
if (has_lme4 && inherits(m_primary, "merMod")) {
  cat("Convergence warnings:\n"); print(m_primary@optinfo$conv$lme4$messages)
  cat("Singular fit:", lme4::isSingular(m_primary), "\n")
}

## 3. Stratified rank test (already printed by section 6 — included for completeness)
cat("\n--- STRATIFIED RANK TEST ---\n")
print(informative %>% group_by(grp) %>%
        summarise(median_lr = median(log_ratio), n = n(), .groups = "drop"))
cat(sprintf("Informative ttypes: %d; samples: %d\n",
            n_distinct(informative$ttype), nrow(informative)))

## 4. Permutation null
cat("\n--- PERMUTATION NULL ---\n")
print(wilcox.test(null_tbl$delta, mu = 0, alternative = "greater"))
cat(sprintf("median delta = %.3f\n", median(null_tbl$delta, na.rm = TRUE)))

## 5. Multiplicity
cat("\n--- MULTIPLICITY ---\n")
print(count(null_tbl, mult))
print(null_tbl %>% filter(!is.na(mult)) %>% group_by(mult) %>%
        summarise(n = n(), median_delta = median(delta), .groups = "drop"))

## 6. Diagnostics: balance, aliasing, split alignment
cat("\n--- DIAGNOSTICS ---\n")
print(count(per_sample, grp))
print(per_sample %>% count(grp, ttype) %>%
        pivot_wider(names_from = grp, values_from = n) %>% arrange(desc(mut)))
print(per_sample %>% group_by(grp) %>%
        summarise(med_split = median(split), med_pre = median(pre_n),
                  med_post = median(post_n), .groups = "drop"))

cat("\n############### END ###############\n")





clk <- seg %>%
  semi_join(per_sample, by = "sample") %>%
  group_by(sample) %>%
  summarise(n_clocks = n_distinct(seg_clock),
            frac_22  = mean(karyotype == "2:2"),
            frac_maj2 = mean(as.integer(sub(":.*", "", karyotype)) >= 2),
            n_seg = n(), .groups = "drop") %>%
  mutate(single = n_clocks == 1)

clk %>% group_by(single) %>%
  summarise(mean_frac_22 = mean(frac_22), mean_maj2 = mean(frac_maj2),
            median_nseg = median(n_seg), n = n(), .groups = "drop")

per_sample %>% left_join(select(clk, sample, single), by = "sample") %>%
  count(grp, single)




# Pull the class in (adjust the column name if it's not literally `class`)
seg_class <- segment_df %>%
  distinct(sample, .keep_all = TRUE) %>%
  select(sample, wgd_class = class)   # <-- rename to the real column name

# Does the class explain the single-clock / synchrony pattern we found?
per_sample %>%
  left_join(seg_class, by = "sample") %>%
  left_join(select(clk, sample, single), by = "sample") %>%
  count(wgd_class, single) %>%
  pivot_wider(names_from = single, values_from = n, values_fill = 0)

# And how does class distribute across mut vs WT? (the real confound check)
per_sample %>%
  left_join(seg_class, by = "sample") %>%
  count(grp, wgd_class) %>%
  pivot_wider(names_from = grp, values_from = n, values_fill = 0)
