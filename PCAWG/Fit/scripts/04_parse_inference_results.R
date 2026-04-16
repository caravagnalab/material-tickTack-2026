rm(list=ls())
.libPaths("~/R/orfeo_R_4.4/")
library(dplyr)
library(tickTack)

data_dir = "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Fit/data/"
inference_results_dir = "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Fit/inference_results_5ncomponents/"
output_dir = "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Fit/results/"

info_smooth_samples <- readRDS(paste0(output_dir,"00_info_fit_smoot5Mb_1Mbml_15mm_0.4pi.RDS"))
fittable_samples = info_smooth_samples %>% filter(is_fittable_after_smoothing == TRUE) %>% pull(sample)
fitted_samples = tools::file_path_sans_ext(list.files(inference_results_dir))
path_sample <- paste0(inference_results_dir, fitted_samples, ".rds")

###### save table with the info for the timed samples in the main Data directory PCAWG/Data ######
table_samples <- info_smooth_samples %>% filter(sample %in% fitted_samples)
saveRDS(table_samples, paste0("/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Data/Samples.rds"))

###### inspect the non fit samples (5 samples) #####
no_fit_available_samplename <- setdiff(unique(fittable_samples), unique(fitted_samples))

failed_runs_path <- paste0(output_dir, "02_failed_runs_after_check_5ncomponents.txt")
writeLines(character(0), failed_runs_path)

# Write samples with no fit object available
if (length(no_fit_available_samplename) > 0) {
  write("=== No available fit object ===", file = failed_runs_path, append = TRUE)
  write(no_fit_available_samplename, file = failed_runs_path, append = TRUE)
  write("", file = failed_runs_path, append = TRUE)
}

write("=== Available fit object but no results (draws_and_summary is NULL) ===", file = failed_runs_path, append = TRUE)

info_segment_summary <- lapply(seq_along(1:length(path_sample)), function(i){
  cat(sprintf("Processing sample %d/%d: %s\n", i, length(path_sample), fitted_samples[i]))
  
  x <- readRDS(path_sample[i])
  info_x <- info_smooth_samples %>% filter(sample == fitted_samples[i])
  results_timing <- x$results_timing
  
  info_x$ncomponents <- length(results_timing$bic_values)
  best_K = results_timing$best_K
  info_x$best_K = best_K
  
  if (is.null(results_timing$draws_and_summary)) {
    cat(sprintf("  -> FAILED: draws_and_summary is NULL\n"))
    write(fitted_samples[i], file = failed_runs_path, append = TRUE)
    return(NULL)
  }
  
  segment_clock_table <- results_timing$draws_and_summary[[best_K]]$summarized_results %>% 
    dplyr::select(segment_name, karyotype, chr, clock_mean, clock_low, clock_high)
  
  info_x_repeated <- info_x[rep(1, nrow(segment_clock_table)), ]
  combined_table <- bind_cols(segment_clock_table, info_x_repeated)
  
  cat(sprintf("  -> Done: %d segments found\n", nrow(segment_clock_table)))
  return(combined_table)
  
}) %>% bind_rows()

saveRDS(info_segment_summary, paste0(output_dir, "04_parsed_inference_results_5ncomponents.rds"))
saveRDS(info_segment_summary, paste0("/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Data/Segments_withoutArmLevelAnn.rds"))

