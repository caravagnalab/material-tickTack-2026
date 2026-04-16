rm(list=ls())
.libPaths("~/R/orfeo_R_4.4/")
library(dplyr)
library(tickTack)

args <- commandArgs(trailingOnly = TRUE)
sample_idx <- as.integer(args[1])

data_dir = "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Fit/data/"
output_dir = "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Fit/results/"
inference_results_dir = "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Fit/inference_results_5ncomponents/"

info_smooth_samples <- readRDS(paste0(output_dir,"00_info_fit_smoot5Mb_1Mbml_15mm_0.4pi.RDS"))
fit_samples_info = info_smooth_samples %>% filter(is_fittable_after_smoothing == TRUE)
path_sample <- fit_samples_info$cnaqc_path

tolerance = 0.0001 #0.001
min_mutations_number = 15
max_distance_CNAqc = 5e6 #1e7 in GEL
min_segment_length = 1e6 #1e7 in GEL


lines <- readLines(failed_runs_path)
section1_start <- which(lines == "=== No available fit object ===") + 1

no_fit_object_samples <- lines[section1_start:(section2_start - 2)]
no_fit_object_samples <- no_fit_object_samples[no_fit_object_samples != ""]

cat("Samples with no result: many mutations -->  too slow inference:", length(no_fit_object_samples), "\n")



# Re-fit samples with no fit object available
sample_name = no_fit_object_samples[sample_idx]
cat(sprintf("Fitting missing sample: %s\n", sample_name))

f <- info_smooth_samples %>% filter(sample == sample_name) %>% pull(cnaqc_path)

x <- readRDS(f)
min_segment_length = 1e6
x$mutations <- x$snvs

  fit_tickTack <- tickTack::fit_tickTack(x = x,
                                         tolerance = tolerance,
                                         min_mutations_number = min_mutations_number,
                                         max_distance_smooth = max_distance_CNAqc,
                                         min_segment_length = min_segment_length,
                                         n_components = 5)

saveRDS(fit_tickTack, paste0(inference_results_dir, sample_name, ".rds"))
cat(sprintf("  -> Done: %s\n", sample_name))
