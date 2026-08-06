# ==============================================================================
# 02_bascule_paired.R
# ------------------------------------------------------------------------------
# Stage 2 of the time-stratified signature analysis.
#
# For one tumour type, runs two BASCULE fits with IDENTICAL hyperparameters,
# IDENTICAL reference catalogue, and IDENTICAL k_list / keep_sigs:
#
#   * fit_early : on the SBS96 matrix built from time < TIME_THRESHOLD only
#   * fit_all   : on the SBS96 matrix built from every mutation
#
# Identical settings matter: any difference in inferred signatures between
# the two fits should be attributable to the data, not to the fitting setup.
#
# Same input/output samples on both fits (stage 01 enforced this), so the
# two exposure tables can be paired sample-by-sample in stage 03.
#
# Submit:
#   sbatch --dependency=afterok:<stage1_jobid> run_bascule_paired.sh
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyverse)
})

# ── Job array index ───────────────────────────────────────────────────────────
ttype_idx <- as.integer(
  Sys.getenv("SLURM_ARRAY_TASK_ID",
             unset = Sys.getenv("SGE_TASK_ID", unset = NA_character_))
)
if (is.na(ttype_idx)) stop("SLURM_ARRAY_TASK_ID / SGE_TASK_ID unset.")

# ── Paths (kept in sync with stage 01) ────────────────────────────────────────
samples_rds <- "/orfeo/cephfs/scratch/cdslab/gsantacatterina/material-tickTack-2026/PCAWG/Data/Samples.rds"
out_root    <- "/orfeo/cephfs/scratch/cdslab/gsantacatterina/material-tickTack-2026/PCAWG/Analysis_results/Signatures_time/"
project     <- "PCAWG"

samples_info <- readRDS(samples_rds)
ttypes       <- unique(samples_info$ttype)
tt           <- ttypes[ttype_idx]
if (is.na(tt)) stop("ttype_idx out of range.")
message("Tumour type: ", tt)

out_data_tt <- file.path(out_root, tt)

# ── BASCULE Python env ────────────────────────────────────────────────────────
reticulate::use_condaenv("/orfeo/scratch/cdslab/ggandolfi/miniconda/envs/bascule-env")
py <- reticulate::import("pybascule")
suppressPackageStartupMessages(library(bascule))

# ── Helper: load a SigProfiler SBS96 matrix as samples-as-rows ────────────────
load_sbs96 <- function(subset_dir, project) {
  sbs_path <- file.path(subset_dir, "output", "SBS", paste0(project, ".SBS96.all"))
  if (!file.exists(sbs_path)) {
    stop("Missing SBS96 matrix: ", sbs_path,
         "\nDid stage 01_partition_and_matrix.R run for this tumour type?")
  }
  m <- read.table(sbs_path, header = TRUE)
  rownames(m) <- m$MutationType
  counts <- m %>%
    dplyr::select(-MutationType) %>%
    t() %>%
    as.data.frame()
  counts[rowSums(counts) > 0, , drop = FALSE]
}

early_counts <- load_sbs96(file.path(out_data_tt, "early"), project)
all_counts   <- load_sbs96(file.path(out_data_tt, "all"),   project)

# ── Align row sets so the two fits see the same samples ───────────────────────
# Stage 01 already enforced this, but a row may have been dropped here by the
# rowSums>0 filter; intersect to be safe.
common_samples <- intersect(rownames(early_counts), rownames(all_counts))
message("Samples in both early & all matrices: ", length(common_samples))

if (length(common_samples) < 5) {
  message("Fewer than 5 paired samples for ", tt, " — BASCULE fits would be unstable. Exiting.")
  quit(save = "no", status = 0)
}

early_counts <- early_counts[common_samples, , drop = FALSE]
all_counts   <- all_counts  [common_samples, , drop = FALSE]

# ── Shared BASCULE settings ───────────────────────────────────────────────────
# Defined once and reused so the two fits are guaranteed identical except for
# the input counts. If you change anything here, change it for BOTH fits.
ref_cat <- list(SBS = bascule::COSMIC_sbs %>% as.data.frame())

bascule_args <- list(
  k_list           = 0,
  reference_cat    = ref_cat,
  keep_sigs        = c("SBS1", "SBS5"),
  hyperparameters  = NULL,
  lr               = 0.005,
  optim_gamma      = 0.1,
  n_steps          = 3000,
  py               = NULL,
  enumer           = "parallel",
  nonparametric    = TRUE,
  autoguide        = FALSE,
  filter_dn        = FALSE,
  min_exposure     = 0.1,
  CUDA             = FALSE,
  compile          = FALSE,
  store_parameters = FALSE,
  store_fits       = TRUE,
  seed_list        = 10
)

fit_one <- function(counts_df, label) {
  message("--- BASCULE fit: ", label, " ---")
  do.call(bascule::fit, c(list(counts = list(SBS = counts_df)), bascule_args))
}

fit_early <- fit_one(early_counts, "early")
fit_all   <- fit_one(all_counts,   "all")

# Refine de novo on both — same call, same args.
fit_early_ref <- refine_denovo_signatures(fit_early)
fit_all_ref   <- refine_denovo_signatures(fit_all)

# ── Save outputs ──────────────────────────────────────────────────────────────
out_dir <- file.path(out_data_tt, "bascule_paired")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(
  list(
    fit_early     = fit_early,
    fit_all       = fit_all,
    fit_early_ref = fit_early_ref,
    fit_all_ref   = fit_all_ref,
    samples       = common_samples
  ),
  file.path(out_dir, "fits.rds")
)

# Exposure plots for visual inspection.
pdf(file.path(out_dir, "exposures.pdf"), width = 12, height = 5)
print(plot_exposures(fit_early_ref, sample_name = TRUE))
print(plot_exposures(fit_all_ref,   sample_name = TRUE))
dev.off()

message("Stage 02 (paired BASCULE) done for tumour type: ", tt)
