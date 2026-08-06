# ==============================================================================
# 00_setup.R
# ------------------------------------------------------------------------------
# Shared setup for the PCAWG signature-analysis pipeline.
#
# This script is `source()`d at the top of 01_matrix_generation.R and
# 02_bascule_inference.R. It is responsible for everything those two pipeline
# stages have in common:
#
#   * reading the HPC job-array index (which tumour type to process)
#   * loading the PCAWG sample metadata
#   * filtering it down to the chosen tumour type
#   * defining the output directory for that tumour type
#
# It does NOT load analysis-specific dependencies (SigProfiler, BASCULE,
# reticulate envs); those belong to the scripts that actually use them, so
# that each stage fails fast and cleanly if its own environment is broken.
#
# Expected usage:
#   sbatch --array=1-30 run_matrix.sh    # launches 01 for each tumour type
#   sbatch --array=1-30 run_bascule.sh   # launches 02 for each tumour type
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyverse)
})

# ── Read the job-array index ──────────────────────────────────────────────────
# SLURM exports SLURM_ARRAY_TASK_ID; SGE/PBS export SGE_TASK_ID. We accept
# either so the same script runs on both schedulers without modification.
ttype_idx <- as.integer(
  Sys.getenv("SLURM_ARRAY_TASK_ID",
             unset = Sys.getenv("SGE_TASK_ID",
                                unset = NA_character_))
)

if (is.na(ttype_idx)) {
  stop("Could not read job array index: ",
       "SLURM_ARRAY_TASK_ID and SGE_TASK_ID are both unset. ",
       "Submit this script as a job array.")
}

# ── Paths ─────────────────────────────────────────────────────────────────────
# Centralised here so the two pipeline stages stay in sync. If the project
# is moved or renamed, this is the only block that needs editing.
samples_rds <- "/orfeo/cephfs/scratch/cdslab/gsantacatterina/material-tickTack-2026/PCAWG/Data/Samples.rds"
out_data    <- "/orfeo/cephfs/scratch/cdslab/gsantacatterina/material-tickTack-2026/PCAWG/Analysis_results/Signatures/"
fit_dir     <- "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Fit/inference_results_5ncomponents/"
project     <- "PCAWG"

dir.create(out_data, recursive = TRUE, showWarnings = FALSE)

# ── Load and filter sample metadata ───────────────────────────────────────────
samples_info <- readRDS(samples_rds)

ttypes <- unique(samples_info$ttype)
tt     <- ttypes[ttype_idx]

if (is.na(tt)) {
  stop("ttype_idx = ", ttype_idx,
       " is out of range. There are only ", length(ttypes), " tumour types.")
}

message("ttype_idx: ", ttype_idx, " | tumour type: ", tt)

# Keep only the samples for this tumour type, and assign a stable per-class
# synthetic ID (S_<class><n>) — SigProfiler and BASCULE both rely on this
# `sampleId` downstream, so it has to be generated identically by both
# pipeline stages. That's exactly why it lives in the shared setup.
samples_info <- samples_info %>%
  dplyr::filter(ttype == tt) %>%
  dplyr::select(sample, ttype, class) %>%
  dplyr::group_by(class) %>%
  dplyr::mutate(sampleId = paste0("S_", class, dplyr::row_number())) %>%
  dplyr::ungroup()

# ── Output directory for this tumour type ─────────────────────────────────────
# Defined here (not inside each pipeline stage) so stage 01 can write the
# SigProfiler matrix and stage 02 can read it from the same location without
# either script having to redefine the path.
out_data_tt <- file.path(out_data, tt)
dir.create(out_data_tt, recursive = TRUE, showWarnings = FALSE)
