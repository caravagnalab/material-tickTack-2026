rm(list=ls())
.libPaths("~/R/orfeo_R_4.4/")
library(dplyr)
library(tickTack)

args <- commandArgs(trailingOnly = TRUE)
chunk_id <- as.integer(args[1])
total_chunks <- as.integer(args[2])

cat("subjob", chunk_id, "of", total_chunks, "\n")   

data_dir = "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Fit/data/"
inference_results_dir = "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Fit/inference_results/"

info_smooth_samples <- readRDS(paste0(data_dir,"00_info_fit_smoot5Mb_1Mbml_15mm_0.4pi.RDS"))
fit_samples_info = info_smooth_samples %>% filter(is_fittable_after_smoothing == TRUE)
path_sample <- fit_samples_info$cnaqc_path

num_files <- length(path_sample)
chunk_size <- ceiling(num_files / total_chunks)
start_index <- (chunk_id - 1) * chunk_size + 1
end_index <- min(chunk_id * chunk_size, num_files)
chunk_files <- path_sample[start_index:end_index]            

message(sprintf("Processing files %d to %d out of %d", start_index, end_index, num_files) )

tolerance = 0.0001 #0.001 
min_mutations_number = 15
max_distance_CNAqc = 5e6 #1e7 in GEL 
min_segment_length = 1e6 #1e7 in GEL 


lapply(chunk_files, function(f){
  x <- readRDS(f)
  min_segment_length = 1e6
  
  x$mutations <- x$snvs

  fit_tickTack <- tickTack::fit_tickTack(x = x,
                                         tolerance = tolerance,
                                         min_mutations_number = min_mutations_number,
                                         max_distance_smooth = max_distance_CNAqc, #1e7 in GEL
                                         min_segment_length = min_segment_length)
  
  # res_timing <- fit_tickTack$results_timing
  # best_K <- res_timing$best_K
  # summary <- res_timing$draws_and_summary[[1]]$summarized_results
  
  saveRDS(fit_tickTack, paste0(inference_results_dir,basename(dirname(f)),".rds"))
  print("done")
})
