# ==============================================================================
# 01_partition_and_matrix.R
# ------------------------------------------------------------------------------
# Stage 1 of the time-stratified signature analysis.
#
# Input:  one cohort-wide txt with rows = mutations, columns including
#         sample, chr, from, to, ref, alt, karyotype, time. `time` is the
#         CNAqc-derived molecular-clock estimate for the segment carrying
#         the mutation.
#
# What this script does, per tumour type (job-array task):
#   1. Reads the cohort table and the sample metadata (for ttype + class).
#   2. For each sample of this tumour type, builds two mutation sets:
#         * early : rows with time < TIME_THRESHOLD
#         * all   : every row (any time, including diploid / NA-time rows)
#   3. Drops samples whose `early` set has fewer than MIN_EARLY_MUTS rows
#      from BOTH sets, so the two matrices have identical sample rows.
#   4. Writes two SigProfiler input files (early / all) and runs
#      SigProfilerMatrixGeneratorR on each, producing two SBS96 matrices.
#
# Outputs land under:
#   <out_root>/<tt>/early/  and  <out_root>/<tt>/all/
# Each containing a PCAWG.txt and a SigProfilerMatrixGeneratorR output tree.
#
# Submit:
#   sbatch --array=1-30 run_partition.sh
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyverse)
  library(data.table)
})

# ── Tunables ──────────────────────────────────────────────────────────────────
# Fixed across samples, as requested. Bumped to a const so it's visible at the
# top of the file rather than buried in a filter() call.
TIME_THRESHOLD  <- 0.6
MIN_EARLY_MUTS  <- 100

# ── Job array index ───────────────────────────────────────────────────────────
ttype_idx <- as.integer(
  Sys.getenv("SLURM_ARRAY_TASK_ID",
             unset = Sys.getenv("SGE_TASK_ID", unset = NA_character_))
)
if (is.na(ttype_idx)) {
  stop("SLURM_ARRAY_TASK_ID / SGE_TASK_ID unset. Submit as a job array.")
}

# ── Paths ─────────────────────────────────────────────────────────────────────
cohort_txt   <- "/orfeo/scratch/cdslab/antonelloa/material-tickTack-2026/PCAWG/Fit/SBS_data/alpha_muts.txt"
samples_rds  <- "/orfeo/cephfs/scratch/cdslab/gsantacatterina/material-tickTack-2026/PCAWG/Data/Samples.rds"
out_root     <- "/orfeo/cephfs/scratch/cdslab/gsantacatterina/material-tickTack-2026/PCAWG/Analysis_results/Signatures_time/"
project      <- "PCAWG"

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

# ── Load metadata + pick tumour type for this task ────────────────────────────
samples_info <- readRDS(samples_rds)
ttypes       <- unique(samples_info$ttype)
tt           <- ttypes[ttype_idx]
if (is.na(tt)) stop("ttype_idx ", ttype_idx, " out of range (", length(ttypes), " tumour types).")
message("Tumour type: ", tt)

samples_info_tt <- samples_info %>%
  dplyr::filter(ttype == tt) %>%
  dplyr::select(sample, ttype, class) %>%
  dplyr::group_by(class) %>%
  dplyr::mutate(sampleId = paste0("S_", class, dplyr::row_number())) %>%
  dplyr::ungroup()

samples_this_tt <- samples_info_tt$sample

# ── Load the cohort mutation table ────────────────────────────────────────────
# fread handles the optional leading row-number column the original snippet
# showed; we drop it if present. Adjust if your real file uses a different
# delimiter or has no header.
muts <- data.table::fread(cohort_txt)
if (colnames(muts)[1] %in% c("V1", "")) muts[, 1 := NULL]

# Restrict early to this tumour type's samples.
muts_tt <- muts %>%
  dplyr::filter(sample %in% samples_this_tt) %>%
  dplyr::mutate(
    # Normalise chr to bare "1", "2", ... — SigProfiler expects no "chr" prefix.
    chrom = sub("^chr", "", chr)
  )

message("Cohort rows for ", tt, ": ", nrow(muts_tt),
        " across ", dplyr::n_distinct(muts_tt$sample), " samples.")

# ── Per-sample count check on the EARLY subset ────────────────────────────────
# Build the early set once and use it to decide which samples survive. Samples
# below the floor are dropped from BOTH `early` and `all` matrices so the two
# resulting matrices have the same row set — that's a precondition for the
# paired comparison done in stage 03.
early_counts <- muts_tt %>%
  dplyr::filter(!is.na(time), time < TIME_THRESHOLD) %>%
  dplyr::count(sample, name = "n_early")

passing_samples <- early_counts %>%
  dplyr::filter(n_early >= MIN_EARLY_MUTS) %>%
  dplyr::pull(sample)

dropped <- setdiff(samples_this_tt, passing_samples)
message("Samples passing the ", MIN_EARLY_MUTS,
        "-mutation early floor: ", length(passing_samples),
        " / ", length(samples_this_tt),
        " (dropped: ", length(dropped), ")")

# Persist the attrition record so reviewers can audit it later.
out_data_tt <- file.path(out_root, tt)
dir.create(out_data_tt, recursive = TRUE, showWarnings = FALSE)
saveRDS(
  list(
    threshold      = TIME_THRESHOLD,
    min_early_muts = MIN_EARLY_MUTS,
    passing        = passing_samples,
    dropped        = dropped,
    early_counts   = early_counts
  ),
  file.path(out_data_tt, "sample_attrition.rds")
)

if (length(passing_samples) == 0) {
  message("No samples pass the early-mutation floor for ", tt, ". Exiting cleanly.")
  quit(save = "no", status = 0)
}

# ── Reusable: long mutation table -> SigProfiler-format table ─────────────────
# `subset_df` is a per-sample mutation table; `subset_name` is "early" or "all"
# and is only used for logging. Returns a data.frame ready to write.
to_sigprofiler_format <- function(subset_df, samples_info_tt, project,
                                  genome = "GRCh37") {
  subset_df %>%
    # Use synthetic sampleId so downstream signature outputs key on the same
    # IDs we use in the main pipeline. left_join keeps mutations whose sample
    # is in samples_info_tt.
    dplyr::inner_join(samples_info_tt %>% dplyr::select(sample, sampleId),
                      by = "sample") %>%
    dplyr::select(sampleId, chrom, from, to, ref, alt) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      Project  = project,
      Genome   = genome,
      mut_type = "SNP",
      Type     = "SOMATIC",
      ID       = sampleId,
      Sample   = sampleId
    ) %>%
    dplyr::rename(pos_start = from, pos_end = to) %>%
    dplyr::select(Project, Sample, ID, Genome, mut_type,
                  chrom, pos_start, pos_end, ref, alt, Type) %>%
    dplyr::filter(ref != alt) %>%
    dplyr::distinct()
}

# ── Build the two subsets ─────────────────────────────────────────────────────
muts_passing <- muts_tt %>% dplyr::filter(sample %in% passing_samples)

early_df <- muts_passing %>%
  dplyr::filter(!is.na(time), time < TIME_THRESHOLD)

all_df   <- muts_passing
# Note: `all_df` keeps mutations with NA time (e.g. diploid segments). If you
# want to restrict `all` to "all mutations that had a time estimate", swap to:
#   all_df <- muts_passing %>% dplyr::filter(!is.na(time))

# ── SigProfiler env (matches stage 01 of the main pipeline) ───────────────────
reticulate::use_python("/orfeo/scratch/area/lvaleriani/myconda/bin")
reticulate::py_config()
suppressPackageStartupMessages(library(SigProfilerMatrixGeneratorR))

# ── Run SigProfiler twice: early, then all ────────────────────────────────────
run_sigprofiler_subset <- function(subset_df, subset_name) {
  message("--- SigProfiler: ", subset_name, " (", nrow(subset_df), " rows) ---")
  subset_dir <- file.path(out_data_tt, subset_name)
  dir.create(subset_dir, recursive = TRUE, showWarnings = FALSE)

  table_sp <- to_sigprofiler_format(subset_df, samples_info_tt, project)

  write.table(
    table_sp,
    file      = file.path(subset_dir, paste0(project, ".txt")),
    quote     = FALSE,
    sep       = "\t",
    row.names = FALSE
  )

  SigProfilerMatrixGeneratorR(
    project     = project,
    genome      = "GRCh37",
    matrix_path = subset_dir,
    plot        = FALSE
  )
}

run_sigprofiler_subset(early_df, "early")
run_sigprofiler_subset(all_df,   "all")

message("Stage 01 (partition + matrices) done for tumour type: ", tt)
