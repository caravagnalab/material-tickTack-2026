.libPaths("~/R/orfeo_R_4.4/")

library(tickTack)
library(dplyr)
library(data.table)
library(parallel)

base_dir <- "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/simulations/"
source(paste0(base_dir, "/utils.R"))

dir_PCAWG_smooth <- "/orfeo/cephfs/scratch/cdslab/scocomello/material_tickTack/PCAWG/fit_final/"
output_dir <- "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/simulations/data/"

MIN_MUTATIONS <- 10

sample_dir_path <- list.files(dir_PCAWG_smooth, full.names = TRUE)

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", detectCores() - 1))

results_list <- mclapply(seq_along(sample_dir_path), function(i) {
  
  file_path <- sample_dir_path[i]
  sample_id <- basename(file_path)
  
  message(sprintf("Processing sample %d / %d: %s", i, length(sample_dir_path), sample_id))
  
  x_after_inference <- readRDS(file_path)
  
  purity <- x_after_inference$metadata$purity
  median_DP <- median(x_after_inference$mutations$DP)
  
  cna_pre_fit <- x_after_inference$cna %>%
    mutate(len = as.integer(.data$to) - as.integer(.data$from)) %>%
    filter(len >= 1e06)
  
  accepted_data <- tickTack:::prepare_input_data(
    x_after_inference$mutations,
    cna_pre_fit,
    purity,
    possible_k = c("2:1", "2:2", "2:0"),
    alpha = 0.05,
    min_mutations_number = MIN_MUTATIONS
  )
  
  cna <- accepted_data$accepted_cna
  n_seg <- nrow(cna)
  
  coords <- tstrsplit(cna$segment_name, "_", fixed = TRUE)
  cna$len <- as.integer(coords[[3]]) - as.integer(coords[[2]])
  
  seg_vec <- accepted_data$input_data$seg_assignment
  seg_counts <- data.table(segment_id = seg_vec)[, .(n_mut = .N), by = segment_id]
  
  seg_counts_joined <- merge(seg_counts, cna, by = "segment_id", all.x = TRUE)
  seg_counts_joined[, density_mut_per_segment := n_mut / len]
  
  data.frame(
    purity = purity,
    median_DP = median_DP,
    median_density_mut = median(seg_counts_joined$density_mut_per_segment, na.rm = TRUE),
    n_seg = n_seg,
    sample_id = sample_id,
    median_n_mut = median(seg_counts_joined$n_mut, na.rm = TRUE),
    median_seg_length = median(seg_counts_joined$len, na.rm = TRUE)
  )
  
}, mc.cores = n_cores)

# Combine results
info_parameters_PCAWG <- rbindlist(results_list)

# Save output
saveRDS(
  info_parameters_PCAWG,
  file = file.path(output_dir, "00_info_parameters_PCAWG_selectedWithMinMut10.rds")
)

# 
# ########## get info of fitted samples before filtering and smoothing ###########
# sample_dir_path <- list.files(dir_PCAWG_smooth)
# 
# 
# results_list <- vector("list", length(sample_dir_path))
# 
# for (i in seq_along(sample_dir_path)) {
# 
#   r <- sample_dir_path[i]
#   message(sprintf("Processing sample %d / %d: %s", i, length(sample_dir_path), r))
# 
#   x <- readRDS(paste0(dir_PCAWG_smooth, r))
# 
#   info_sample <- tibble(
#     purity = x$metadata$purity,
#     median_DP = median(x$before_smoothing$snvs$DP)
#   )
# 
#   cna <- x$before_smoothing$cna %>%
#     mutate(karyo = paste0(Major,":",minor)) %>%
#     filter(karyo %in% c("2:0", "2:2", "2:1"))
# 
#   n_seg <- nrow(cna)
# 
#   cna <- cna %>% mutate(len = to - from)
# 
#   seg_vec <- x$results_timing$data$input_data$seg_assignment
# 
#   seg_counts <- tibble(segment_id = seg_vec) %>%
#     count(segment_id, name = "n_mut")
# 
#   seg_counts_joined <- seg_counts %>%
#     left_join(cna, by = "segment_id") %>%
#     mutate(density_mut_per_segment = n_mut / len)
# 
#   sample_info_CN_mut <- tibble(
#     purity = info_sample$purity,
#     median_DP = info_sample$median_DP,
#     median_density_mut = median(seg_counts_joined$density_mut_per_segment, na.rm = TRUE),
#     n_seg = n_seg,
#     sample_id = r,
#     median_n_mut = median(seg_counts_joined$n_mut, na.rm = TRUE),
#     median_seg_length = median(seg_counts_joined$len, na.rm = TRUE)
#   )
# 
#   results_list[[i]] <- sample_info_CN_mut
# 
#   rm(
#     x_after_inference,
#     cna,
#     seg_vec,
#     seg_counts,
#     seg_counts_joined,
#     coords,
#     from,
#     to,
#     info_sample
#   )
#   gc(verbose = FALSE)
# }
# 
# info_parameters_PCAWG <- bind_rows(results_list)
# 
# saveRDS(info_parameters_PCAWG, paste0(output_dir, "/00_info_parameters_PCAWG_before_smoothing.rds"))
