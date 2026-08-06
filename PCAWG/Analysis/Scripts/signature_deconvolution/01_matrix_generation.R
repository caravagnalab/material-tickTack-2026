# ==============================================================================
# 01_matrix_generation.R
# ------------------------------------------------------------------------------
# Pipeline stage 1: SigProfiler SBS96 matrix generation.
#
# For one tumour type (selected by the job-array index), this script:
#   1. reads every sample's CNAqc fit object
#   2. extracts SNV calls and reshapes them to SigProfiler's expected layout
#   3. writes the combined table to <out_data_tt>/PCAWG.txt
#   4. calls SigProfilerMatrixGeneratorR, which produces the SBS96 matrix
#      that stage 02 (BASCULE) will consume.
#
# Submitted as a job array, one task per tumour type:
#   sbatch --array=1-30 run_matrix.sh
# ==============================================================================

# Shared setup: defines tt, samples_info, project, fit_dir, out_data_tt.
source("00_setup.R")

# ── Python / SigProfiler environment ──────────────────────────────────────────
# SigProfilerMatrixGeneratorR is an R wrapper around a Python tool, so we
# point reticulate at the conda env that has the Python deps installed.
# This is stage-01 specific; stage 02 uses a different env, which is why
# these reticulate calls live in the per-stage script and not in 00_setup.R.
reticulate::use_python("/orfeo/scratch/area/lvaleriani/myconda/bin")
reticulate::py_config()

suppressPackageStartupMessages({
  library(SigProfilerMatrixGeneratorR)
})

# ── Per-sample extraction ─────────────────────────────────────────────────────
# Reads one CNAqc RDS and returns SNVs in SigProfiler's long format.
# Notes:
#   * chr is stored as "chrN" in CNAqc; SigProfiler expects bare "N", hence
#     the separate()/select(-tmp) dance.
#   * ref != alt filter removes any spurious non-variant rows.
#   * distinct() guards against duplicate mutation entries.
get_sig_profiler_table <- function(i) {
  sample_name <- samples_info$sample[i]
  sample_id   <- samples_info$sampleId[i]
  message("  [", i, "/", nrow(samples_info), "] ", sample_name, " -> ", sample_id)

  cnaqc  <- readRDS(file.path(fit_dir, paste0(sample_name, ".rds")))
  genome <- cnaqc$reference_genome

  # Rewrite the per-mutation sample column to use the synthetic sampleId so
  # all downstream files key on the same identifier across the pipeline.
  cnaqc$mutations$sample <- sample_id

  cnaqc$mutations %>%
    dplyr::select(chr, from, to, ref, alt, sample) %>%
    dplyr::distinct() %>%
    tidyr::separate(col = chr, sep = "chr", into = c("tmp", "chr")) %>%
    dplyr::select(-tmp) %>%
    dplyr::mutate(
      Project  = project,
      Genome   = genome,
      mut_type = "SNP",
      Type     = "SOMATIC",
      ID       = sample,
      Sample   = sample
    ) %>%
    dplyr::rename(chrom = chr, pos_start = from, pos_end = to) %>%
    dplyr::select(Project, Sample, ID, Genome, mut_type,
                  chrom, pos_start, pos_end, ref, alt, Type) %>%
    dplyr::filter(ref != alt) %>%
    dplyr::distinct()
}

# ── Build the combined table for this tumour type ─────────────────────────────
sigprofiler_table <- lapply(seq_len(nrow(samples_info)), get_sig_profiler_table) %>%
  dplyr::bind_rows()

# ── Write the input file and run SigProfiler ──────────────────────────────────
# SigProfilerMatrixGeneratorR reads the .txt from `matrix_path` and writes
# its outputs (including PCAWG.SBS96.all) under <matrix_path>/output/SBS/.
input_txt <- file.path(out_data_tt, paste0(project, ".txt"))

write.table(sigprofiler_table,
            file      = input_txt,
            quote     = FALSE,
            sep       = "\t",
            row.names = FALSE)

SigProfilerMatrixGeneratorR(
  project     = project,
  genome      = "GRCh37",
  matrix_path = out_data_tt,
  plot        = FALSE
)

message("Stage 01 done for tumour type: ", tt)
