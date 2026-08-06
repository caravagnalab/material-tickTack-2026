# ==============================================================================
# 03_compare_exposures.R
# ------------------------------------------------------------------------------
# Stage 3 of the time-stratified signature analysis.
#
# This script is NOT a job-array task; it loops over every tumour type that
# has a paired BASCULE fit from stage 02 and produces three things per
# tumour type, plus a cohort-wide summary:
#
#   1. Paired per-signature deltas (all − early) with a Wilcoxon signed-rank
#      test per signature. Signatures with FDR-significant positive delta are
#      "late-acting"; negative delta are "early-acting".
#   2. Per-sample cosine distance between the early and all exposure vectors,
#      summarising how much the spectrum shifts with the addition of late
#      mutations.
#   3. Presence/absence flips: signatures crossing a min-exposure threshold
#      in one set but not the other.
#
# Outputs land in <out_root>/_summary/ and per-type under <out_root>/<tt>/.
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyverse)
  library(ggplot2)
})

# ── Tunables ──────────────────────────────────────────────────────────────────
# CLEAN_EPS: exposures below this are treated as noise (zeroed). After cleaning,
# per-sample exposures are renormalised so each row sums to 1. The same
# threshold is reused for presence/absence flips so "cleaned out" and "absent"
# mean the same thing.
CLEAN_EPS    <- 0.05
PRESENCE_EPS <- CLEAN_EPS   # kept as alias for readability in presence_flips()
FDR_ALPHA    <- 0.05

# ── Paths ─────────────────────────────────────────────────────────────────────
samples_rds <- "/orfeo/cephfs/scratch/cdslab/gsantacatterina/material-tickTack-2026/PCAWG/Data/Samples.rds"
out_root    <- "/orfeo/cephfs/scratch/cdslab/gsantacatterina/material-tickTack-2026/PCAWG/Analysis_results/Signatures_time/"

samples_info <- readRDS(samples_rds)
ttypes       <- unique(samples_info$ttype)

summary_dir <- file.path(out_root, "_summary")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Helpers
# ==============================================================================

# ── extract_exposures ─────────────────────────────────────────────────────────
# Pull the BASCULE per-sample exposure matrix out of a fit object as a tidy
# tibble with columns (sample, signature, exposure). BASCULE's internal layout
# changes between versions; this wrapper centralises the slot access so a
# future schema change only needs one fix.
extract_exposures <- function(fit) {
  # bascule fits expose exposures via get_exposure() in recent versions.
  # Fall back to direct slot access if that helper isn't present.
  exp_mat <- tryCatch(
    bascule::get_exposure(fit, matrix = TRUE),
    error = function(e) fit$exposure$SBS
  )
  
  as.data.frame(exp_mat) %>%
    tibble::rownames_to_column("sample") %>%
    tidyr::pivot_longer(-sample, names_to = "signature", values_to = "exposure")
}

# ── cosine_distance ───────────────────────────────────────────────────────────
# 1 - cos(angle) between two non-negative vectors, with safe handling of
# all-zero rows (returns NA instead of dividing by zero).
cosine_distance <- function(a, b) {
  num <- sum(a * b)
  den <- sqrt(sum(a^2)) * sqrt(sum(b^2))
  if (den == 0) return(NA_real_)
  1 - num / den
}

# ── clean_exposures ───────────────────────────────────────────────────────────
# BASCULE's `min_exposure` argument prunes during fitting, but
# refine_denovo_signatures() and downstream rescaling can leave tiny residual
# exposures (e.g. 0.002) that are numerical noise rather than real signature
# activity. This helper enforces a hard threshold and renormalises so each
# sample's exposures sum back to 1 — keeping them interpretable as proportions.
#
# A sample whose entire row falls below `eps` is left as all-zeros (rather
# than dividing by zero). Such samples should be rare; if you see them, the
# fit for that sample failed and the row should probably be dropped from the
# comparison.
clean_exposures <- function(exp_long, eps = CLEAN_EPS, renormalise = TRUE) {
  out <- exp_long %>%
    dplyr::mutate(exposure = ifelse(exposure < eps, 0, exposure))
  
  if (renormalise) {
    out <- out %>%
      dplyr::group_by(sample) %>%
      dplyr::mutate(
        row_sum  = sum(exposure),
        exposure = ifelse(row_sum > 0, exposure / row_sum, 0)
      ) %>%
      dplyr::select(-row_sum) %>%
      dplyr::ungroup()
  }
  
  out
}

# ── align_signatures ──────────────────────────────────────────────────────────
# The two fits may discover slightly different de novo signature sets, so the
# exposure tables may have different column sets. We align to the UNION,
# filling missing entries with 0 (a signature absent from one fit contributes
# zero exposure to every sample in that fit). This is the right thing to do
# for the delta comparison.
align_signatures <- function(early_long, all_long) {
  all_sigs <- union(unique(early_long$signature), unique(all_long$signature))
  all_samples <- union(unique(early_long$sample), unique(all_long$sample))
  
  scaffold <- tidyr::expand_grid(sample = all_samples, signature = all_sigs)
  
  e <- scaffold %>%
    dplyr::left_join(early_long, by = c("sample", "signature")) %>%
    dplyr::mutate(exposure = tidyr::replace_na(exposure, 0)) %>%
    dplyr::rename(exposure_early = exposure)
  
  a <- scaffold %>%
    dplyr::left_join(all_long, by = c("sample", "signature")) %>%
    dplyr::mutate(exposure = tidyr::replace_na(exposure, 0)) %>%
    dplyr::rename(exposure_all = exposure)
  
  dplyr::inner_join(e, a, by = c("sample", "signature"))
}

# ── per_signature_wilcoxon ────────────────────────────────────────────────────
# Paired Wilcoxon signed-rank per signature, with BH correction across
# signatures within a tumour type. Returns a tidy result table including
# median delta (effect size; Wilcoxon p-values alone are misleading).
per_signature_wilcoxon <- function(paired_long) {
  paired_long %>%
    dplyr::group_by(signature) %>%
    dplyr::summarise(
      n_samples   = dplyr::n(),
      median_early = median(exposure_early),
      median_all   = median(exposure_all),
      median_delta = median(exposure_all - exposure_early),
      # wilcox.test errors out if all deltas are zero — wrap in tryCatch.
      p_value = tryCatch(
        stats::wilcox.test(exposure_all, exposure_early,
                           paired = TRUE, exact = FALSE)$p.value,
        error = function(e) NA_real_
      ),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      p_adj       = p.adjust(p_value, method = "BH"),
      significant = !is.na(p_adj) & p_adj < FDR_ALPHA,
      direction   = dplyr::case_when(
        !significant       ~ "ns",
        median_delta > 0   ~ "late-acting",   # higher in `all`
        median_delta < 0   ~ "early-acting",  # higher in `early`
        TRUE               ~ "ns"
      )
    ) %>%
    dplyr::arrange(p_adj)
}

# ── presence_flips ────────────────────────────────────────────────────────────
# Count samples where a signature is present in one set but not the other,
# using PRESENCE_EPS as the threshold. Useful for sparse signatures where
# medians hide on/off behaviour.
presence_flips <- function(paired_long, eps = PRESENCE_EPS) {
  paired_long %>%
    dplyr::mutate(
      in_early = exposure_early > eps,
      in_all   = exposure_all   > eps
    ) %>%
    dplyr::group_by(signature) %>%
    dplyr::summarise(
      gained_in_all   = sum(!in_early &  in_all),   # appears once late muts added
      lost_in_all     = sum( in_early & !in_all),   # disappears (rare; usually noise)
      present_in_both = sum( in_early &  in_all),
      absent_in_both  = sum(!in_early & !in_all),
      .groups = "drop"
    )
}

# ── per_sample_cosine ─────────────────────────────────────────────────────────
# One cosine distance per sample, computed over the aligned (union) signature
# space. Returns a tibble (sample, cosine_distance).
per_sample_cosine <- function(paired_long) {
  paired_long %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(
      cosine_distance = cosine_distance(exposure_early, exposure_all),
      .groups = "drop"
    )
}

# ==============================================================================
# Main loop
# ==============================================================================

# Accumulators for cohort-wide summaries.
all_wilcox  <- list()
all_cosines <- list()
all_flips   <- list()

for (tt in ttypes) {
  message("==== ", tt, " ====")
  fits_path <- file.path(out_root, tt, "bascule_paired", "fits.rds")
  if (!file.exists(fits_path)) {
    message("  No paired fit found; skipping.")
    next
  }
  
  fits <- readRDS(fits_path)
  
  # Pull exposures, then clean: zero out anything below CLEAN_EPS and
  # renormalise per sample so rows sum to 1. This removes the long tail of
  # numerically-nonzero-but-biologically-meaningless exposures BASCULE leaves
  # behind after refinement.
  early_long <- extract_exposures(fits$fit_early_ref) %>% clean_exposures()
  all_long   <- extract_exposures(fits$fit_all_ref)   %>% clean_exposures()
  
  paired <- align_signatures(early_long, all_long)
  
  # Keep only samples that appear in BOTH original fits (defensive — the union
  # in align_signatures might have introduced one-sided rows if the fits had
  # different sample sets, which shouldn't happen but is cheap to guard).
  paired <- paired %>%
    dplyr::group_by(sample) %>%
    dplyr::filter(any(exposure_early > 0) | any(exposure_all > 0)) %>%
    dplyr::ungroup()
  
  # Per-signature paired tests.
  wilcox_tt <- per_signature_wilcoxon(paired) %>%
    dplyr::mutate(ttype = tt, .before = 1)
  
  # Per-sample cosine distances.
  cosine_tt <- per_sample_cosine(paired) %>%
    dplyr::mutate(ttype = tt, .before = 1)
  
  # Presence/absence flips.
  flips_tt <- presence_flips(paired) %>%
    dplyr::mutate(ttype = tt, .before = 1)
  
  # Save per-tumour-type artefacts.
  tt_out <- file.path(out_root, tt, "comparison")
  dir.create(tt_out, recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(wilcox_tt, file.path(tt_out, "wilcoxon_per_signature.tsv"))
  readr::write_tsv(cosine_tt, file.path(tt_out, "cosine_per_sample.tsv"))
  readr::write_tsv(flips_tt,  file.path(tt_out, "presence_flips.tsv"))
  saveRDS(paired,             file.path(tt_out, "paired_exposures.rds"))
  
  all_wilcox[[tt]]  <- wilcox_tt
  all_cosines[[tt]] <- cosine_tt
  all_flips[[tt]]   <- flips_tt
}

# ── Cohort-wide concatenations ────────────────────────────────────────────────
wilcox_all  <- dplyr::bind_rows(all_wilcox)
cosines_all <- dplyr::bind_rows(all_cosines)
flips_all   <- dplyr::bind_rows(all_flips)

readr::write_tsv(wilcox_all,  file.path(summary_dir, "wilcoxon_per_signature_all_ttypes.tsv"))
readr::write_tsv(cosines_all, file.path(summary_dir, "cosine_per_sample_all_ttypes.tsv"))
readr::write_tsv(flips_all,   file.path(summary_dir, "presence_flips_all_ttypes.tsv"))

# ── Cohort-level summary plots ────────────────────────────────────────────────
# 1. Heatmap of median delta per (ttype, signature), starred where FDR sig.
p_heatmap <- wilcox_all %>%
  ggplot(aes(x = signature, y = ttype, fill = median_delta)) +
  geom_tile() +
  geom_text(aes(label = ifelse(significant, "*", "")),
            colour = "black", size = 3, vjust = 0.78) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                       midpoint = 0,
                       name = "median(all - early)") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title    = "Per-signature exposure shift: all - early",
       subtitle = "* = FDR < 0.05, paired Wilcoxon signed-rank within tumour type",
       x = NULL, y = NULL)

ggsave(file.path(summary_dir, "delta_heatmap.pdf"), p_heatmap, width = 12, height = 7)

# 2. Boxplot of per-sample cosine distances by tumour type. 
p_cos <- cosines_all %>%
  ggplot(aes(x = reorder(ttype, cosine_distance, FUN = median, na.rm = TRUE),
             y = cosine_distance)) +
  geom_boxplot(outlier.size = 0.6, fill = "grey90") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Per-sample spectrum shift (cosine distance between early and all)",
       x = NULL, y = "1 - cos(early, all)")

ggsave(file.path(summary_dir, "cosine_by_ttype.pdf"), p_cos, width = 10, height = 5)

message("Stage 03 done. Summary outputs under: ", summary_dir)