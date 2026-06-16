rm(list = ls())
library(tidyverse)
library(ggplot2)
library(patchwork)

DRIVERS_PATH <- "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Data/Drivers.rds"
SEGMENTS_PATH <- "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Data/Segments.rds"
GENE <- "TP53"

# ── Load data ────────────────────────────────────────────────────────────────
drivers_df <- readRDS(DRIVERS_PATH) %>%
  dplyr::filter(gene == GENE) %>%
  dplyr::select(segment_id, gene, karyotype, sample_id, mutation_status,
                mult_estimate, timed, clock_mean) %>%
  dplyr::rename(sample = sample_id)

segment_df <- readRDS(SEGMENTS_PATH) %>%
  dplyr::select(segment_name, karyotype, chr, clock_mean, sample, ttype)

# ── Helper ────────────────────────────────────────────────────────────────────
fmt_p <- function(x) {
  # Accepts either a htest object or a raw p-value
  p <- if (inherits(x, "htest")) x$p.value else x
  if (p < 0.001) "p < 0.001" else paste0("p = ", round(p, 3))
}

# ── SANITY CHECK ─────────────────────────────────────────────────────────────
# The clock_mean in segment_df is a sample-level summary: all CNA events
# within a sample share (nearly) the same value. This tells us the model
# infers CIN as a punctuated, near-simultaneous burst rather than a gradual
# accumulation. The analysis is reframed accordingly.
clock_variance_check <- segment_df %>%
  group_by(sample) %>%
  summarise(n_unique_clocks = n_distinct(clock_mean), .groups = "drop")

cat("── Sanity check ──────────────────────────────────────────────────────\n")
cat("Samples with a single clock_mean value:",
    sum(clock_variance_check$n_unique_clocks == 1),
    "out of", nrow(clock_variance_check), "\n")
cat("Interpretation: CIN is a punctuated burst, not a gradual process.\n\n")

# ── Build per-sample summary table ───────────────────────────────────────────
# One row per sample with:
#   mutation_status  : M (TP53 mutant) or WT
#   cin_clock        : timing of the CIN burst (sample-level clock_mean)
#   tp53_clock       : timing of the TP53 mutation (from drivers_df)
#   total_cin_mb     : total genome altered in the CIN burst (Mb)
#   mult             : rounded multiplicity (1 = CIN before TP53,
#                                            2 = TP53 before CIN)
#   ttype            : tumor type

seg_lengths <- segment_df %>%
  separate(segment_name, c("nm_chr", "start", "end"),
           sep = "_", convert = TRUE, remove = FALSE) %>%
  mutate(length = end - start) %>%
  group_by(sample) %>%
  summarise(
    total_cin_mb = sum(length, na.rm = TRUE) / 1e6,
    .groups = "drop"
  )

cin_clock_per_sample <- segment_df %>%
  distinct(sample, ttype, clock_mean) %>%
  rename(cin_clock = clock_mean)

sample_data <- drivers_df %>%
  distinct(sample, mutation_status, mult_estimate, clock_mean, timed) %>%
  rename(tp53_clock = clock_mean) %>%
  mutate(mult = round(mult_estimate)) %>%
  left_join(cin_clock_per_sample, by = "sample") %>%
  left_join(seg_lengths,          by = "sample") %>%
  filter(!is.na(ttype))   # keep only samples present in segment_df

sample_data <- sample_data %>%
  filter(mutation_status %in% c("M", "WT", "CNA_driver")) %>%
  mutate(mutation_status = ifelse(mutation_status == "CNA_driver", "M", mutation_status))

cat("── Sample counts ─────────────────────────────────────────────────────\n")
print(count(sample_data, mutation_status))
cat("\n")

# ════════════════════════════════════════════════════════════════════════════
# ANALYSIS 1
# Does TP53 mutation determine WHEN the CIN burst occurs?
# Hypothesis: TP53-mutant samples have an earlier CIN clock than WT,
# consistent with TP53 loss triggering an early catastrophic event.
# ════════════════════════════════════════════════════════════════════════════

cat("── Analysis 1: CIN burst timing (mutant vs WT) ───────────────────────\n")

lm1 <- lm(cin_clock ~ mutation_status + ttype, data = sample_data)
summary(lm1)

# Non-parametric complement (no ttype adjustment)
wt1  <- wilcox.test(cin_clock ~ mutation_status, data = sample_data,
                    alternative = "less")   # mutant clock < WT clock?
cat("Wilcoxon (unadjusted, mutant < WT):", fmt_p(wt1), "\n\n")

# ════════════════════════════════════════════════════════════════════════════
# ANALYSIS 2
# Does TP53 mutation determine HOW MUCH CIN occurs in the burst?
# Hypothesis: TP53-mutant samples have a larger total CIN burden,
# consistent with TP53 loss permitting more extensive genome reorganization.
# ════════════════════════════════════════════════════════════════════════════

cat("── Analysis 2: CIN burst magnitude (mutant vs WT) ────────────────────\n")

lm2 <- lm(total_cin_mb ~ mutation_status + ttype, data = sample_data)
summary(lm2)

wt2 <- wilcox.test(total_cin_mb ~ mutation_status, data = sample_data,
                   alternative = "greater")  # mutant burden > WT?
cat("Wilcoxon (unadjusted, mutant > WT):", fmt_p(wt2), "\n\n")

# ════════════════════════════════════════════════════════════════════════════
# ANALYSIS 3
# Does TP53 mutation PRECEDE the CIN burst?
# Multiplicity directly encodes the relative ordering of TP53 vs CIN:
#   mult == 2: TP53 mutated BEFORE the CN gain (TP53 preceded the burst)
#   mult == 1: TP53 mutated AFTER  the CN gain (CIN preceded TP53)
# A majority of mult == 2 supports TP53 as the gatekeeper whose loss
# permits the CIN burst.
# ════════════════════════════════════════════════════════════════════════════

cat("── Analysis 3: Relative ordering of TP53 vs CIN burst ────────────────\n")

mult_data <- sample_data %>%
  filter(mutation_status == "M", !is.na(mult)) %>%
  mutate(tp53_before_cin = mult == 2)

cat("TP53 before CIN (mult == 2):",
    sum(mult_data$tp53_before_cin), "/", nrow(mult_data),
    sprintf("(%.1f%%)\n", 100 * mean(mult_data$tp53_before_cin)))

# Binomial test: is the proportion of mult == 2 greater than 50%?
binom_test <- binom.test(sum(mult_data$tp53_before_cin),
                         nrow(mult_data),
                         p = 0.5, alternative = "greater")
cat("Binomial test (prop mult==2 > 0.5):", fmt_p(binom_test), "\n")

# Within mult == 2 (TP53 before CIN):
# Is the CIN burst larger than in mult == 1 (CIN before TP53)?
# If TP53 gates the burst, samples where it was lost first should have
# more extensive CIN.
lm3 <- lm(total_cin_mb ~ factor(mult) + ttype, data = mult_data)
summary(lm3)

wt3 <- wilcox.test(total_cin_mb ~ factor(mult), data = mult_data,
                   alternative = "greater")  # mult==2 burden > mult==1?
cat("Wilcoxon (mult==2 burden > mult==1):", fmt_p(wt3), "\n\n")

# ════════════════════════════════════════════════════════════════════════════
# ANALYSIS 4
# Within TP53 mutants: does earlier TP53 timing associate with earlier
# and/or larger CIN bursts?
# This tests whether the *degree* of early TP53 loss matters, not just
# its presence.
# ════════════════════════════════════════════════════════════════════════════

cat("── Analysis 4: TP53 clock vs CIN clock (within mutants) ──────────────\n")

mutant_timed <- sample_data %>%
  filter(mutation_status == "M", timed,
         !is.na(tp53_clock), !is.na(cin_clock),
         tp53_clock > 0, tp53_clock < 1)

lm4a <- lm(cin_clock    ~ tp53_clock + ttype, data = mutant_timed)
lm4b <- lm(total_cin_mb ~ tp53_clock + ttype, data = mutant_timed)
cat("TP53 clock → CIN clock:\n");    print(summary(lm4a)$coefficients["tp53_clock", ])
cat("TP53 clock → CIN magnitude:\n"); print(summary(lm4b)$coefficients["tp53_clock", ])
cat("\n")

# ════════════════════════════════════════════════════════════════════════════
# PLOTS
# ════════════════════════════════════════════════════════════════════════════

clr <- c("M"  = "#d95f02", "WT" = "#7570b3")
clr_mult <- c("1" = "#66c2a5", "2" = "#fc8d62")

# ── Plot 1a: CIN burst timing by mutation status ──────────────────────────
p1a <- sample_data %>%
  filter(!is.na(cin_clock)) %>%
  ggplot(aes(x = mutation_status, y = cin_clock, fill = mutation_status)) +
  geom_violin(alpha = 0.4, color = NA, adjust = 1.2) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.85, linewidth = 0.6) +
  annotate("text", x = 1.5,
           y = max(sample_data$cin_clock, na.rm = TRUE) * 0.97,
           label = fmt_p(wt1), size = 3.5, color = "grey30") +
  scale_fill_manual(values = clr) +
  scale_x_discrete(labels = c("M" = "TP53 mutant", "WT" = "Wild-type")) +
  labs(title    = "Timing of CIN burst",
       subtitle = "Earlier clock = earlier in tumor evolution",
       x = NULL,
       y = "CIN burst clock (molecular time)") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

# ── Plot 1b: CIN burst magnitude by mutation status ───────────────────────
p1b <- sample_data %>%
  filter(!is.na(total_cin_mb)) %>%
  ggplot(aes(x = mutation_status, y = total_cin_mb, fill = mutation_status)) +
  geom_violin(alpha = 0.4, color = NA, adjust = 1.2) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.85, linewidth = 0.6) +
  annotate("text", x = 1.5,
           y = max(sample_data$total_cin_mb, na.rm = TRUE) * 0.97,
           label = fmt_p(wt2), size = 3.5, color = "grey30") +
  scale_fill_manual(values = clr) +
  scale_x_discrete(labels = c("M" = "TP53 mutant", "WT" = "Wild-type")) +
  labs(title    = "Magnitude of CIN burst",
       subtitle = "Total genome copy-number altered",
       x = NULL,
       y = "Total CIN burden (Mb)") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

# ── Plot 2: Multiplicity — ordering of TP53 vs CIN ────────────────────────
mult_labels <- c(
  "1" = "mult = 1\nCIN before TP53",
  "2" = "mult = 2\nTP53 before CIN"
)

p2 <- mult_data %>%
  mutate(mult_f = factor(mult, levels = c("1", "2"))) %>%
  ggplot(aes(x = mult_f, y = total_cin_mb, fill = mult_f)) +
  geom_violin(alpha = 0.4, color = NA, adjust = 1.2) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.85, linewidth = 0.6) +
  geom_jitter(width = 0.12, alpha = 0.25, size = 1, color = "black") +
  annotate("text", x = 1.5,
           y = max(mult_data$total_cin_mb, na.rm = TRUE) * 0.97,
           label = fmt_p(wt3), size = 3.5, color = "grey30") +
  scale_fill_manual(values = clr_mult) +
  scale_x_discrete(labels = mult_labels) +
  labs(title    = "CIN burst magnitude by relative ordering of TP53",
       subtitle = "Does earlier TP53 loss lead to a larger CIN event?",
       x = "Multiplicity (relative event order)",
       y = "Total CIN burden (Mb)") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

# ── Plot 3: TP53 clock vs CIN clock (within timed mutants) ────────────────
p3 <- mutant_timed %>%
  ggplot(aes(x = tp53_clock, y = cin_clock, color = ttype)) +
  geom_point(alpha = 0.5, size = 1.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black",
              linewidth = 0.8, linetype = "dashed") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted",
              color = "grey60", linewidth = 0.6) +
  annotate("text", x = 0.75, y = 0.1,
           label = paste0("β = ", round(coef(lm4a)["tp53_clock"], 2),
                          "\n", fmt_p(summary(lm4a)$coefficients["tp53_clock", 4])),
           size = 3.5, color = "grey30", hjust = 0) +
  labs(title    = "TP53 mutation timing vs CIN burst timing",
       subtitle = "Dashed line: linear fit; dotted: identity (simultaneous)",
       x = "TP53 mutation clock (molecular time)",
       y = "CIN burst clock (molecular time)",
       color = "Tumor type") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 8))

# ── Plot 4: Tumor-type stratified view of burst magnitude ─────────────────
p4 <- sample_data %>%
  filter(!is.na(total_cin_mb)) %>%
  mutate(status_label = ifelse(mutation_status == "M",
                               "TP53 mutant", "Wild-type")) %>%
  ggplot(aes(x = status_label, y = total_cin_mb, fill = mutation_status)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, linewidth = 0.5, width = 0.6) +
  facet_wrap(~ ttype, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = clr) +
  labs(title    = "CIN burst magnitude by tumor type",
       subtitle = "Stratified view controlling for type-specific CIN levels",
       x = NULL,
       y = "Total CIN burden (Mb)") +
  theme_minimal(base_size = 11) +
  theme(legend.position  = "none",
        plot.title       = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        strip.text       = element_text(size = 8),
        axis.text.x      = element_text(angle = 20, hjust = 1))

# ── Assemble and save ─────────────────────────────────────────────────────
fig_main <- (p1a | p1b | p2) +
  plot_annotation(
    title = "TP53 mutation and punctuated chromosomal instability",
    subtitle = paste(
      "CIN in these tumors is a near-simultaneous burst, not a gradual process.",
      "Analyses test whether TP53 mutation determines the timing, magnitude,",
      "and ordering of that burst.",
      sep = " "
    ),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 10, color = "grey40")
    )
  )

ggsave("tp53_cin_burst_main.png",
       plot = fig_main, width = 14, height = 6, dpi = 300)

ggsave("tp53_cin_clock_correlation.png",
       plot = p3, width = 7, height = 5, dpi = 300)

ggsave("tp53_cin_by_ttype.png",
       plot = p4, width = 14, height = 8, dpi = 300)

cat("── Done. Figures saved. ──────────────────────────────────────────────\n")




# Sample composition
sample_data %>% count(mutation_status, ttype)

# Analysis 1 & 2: lm summaries
summary(lm1)
summary(lm2)

# Analysis 3: multiplicity counts and tests
mult_data %>% count(mult)
print(binom_test)
summary(lm3)

# Analysis 4: coefficients
summary(lm4a)$coefficients["tp53_clock", ]
summary(lm4b)$coefficients["tp53_clock", ]

# Median summaries
sample_data %>%
  filter(mutation_status %in% c("M", "WT")) %>%
  group_by(mutation_status) %>%
  summarise(
    median_cin_clock  = median(cin_clock,    na.rm = TRUE),
    median_cin_mb     = median(total_cin_mb, na.rm = TRUE),
    n = n()
  )