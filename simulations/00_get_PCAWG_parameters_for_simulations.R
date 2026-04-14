.libPaths("~/R/orfeo_R_4.4/")
library(tickTack)
library(ggplot2)
library(dplyr)
dir_PCAWG_smooth = "/orfeo/cephfs/scratch/cdslab/scocomello/material_tickTack/PCAWG/fit_final/"
output_dir = "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/simulations/data/"

sample_dir_path <- list.files(dir_PCAWG_smooth)


results_list <- vector("list", length(sample_dir_path))

for (i in seq_along(sample_dir_path)) {
  
  r <- sample_dir_path[i]
  message(sprintf("Processing sample %d / %d: %s", i, length(sample_dir_path), r))
  
  x_after_inference <- readRDS(paste0(dir_PCAWG_smooth, r))
  
  info_sample <- tibble(
    purity = x_after_inference$metadata$purity, 
    median_DP = median(x_after_inference$mutations$DP)
  )
  
  cna <- x_after_inference$results_timing$data$accepted_cna
  n_seg <- nrow(cna)
  
  coords <- strsplit(cna$segment_name, "_")
  from <- as.integer(sapply(coords, `[`, 2))
  to   <- as.integer(sapply(coords, `[`, 3))
  cna$len <- to - from
  
  seg_vec <- x_after_inference$results_timing$data$input_data$seg_assignment
  
  seg_counts <- tibble(segment_id = seg_vec) %>%
    count(segment_id, name = "n_mut")
  
  seg_counts_joined <- seg_counts %>%
    left_join(cna, by = "segment_id") %>%
    mutate(density_mut_per_segment = n_mut / len)
  
  sample_info_CN_mut <- tibble(
    purity = info_sample$purity,
    median_DP = info_sample$median_DP,
    median_density_mut = median(seg_counts_joined$density_mut_per_segment, na.rm = TRUE),
    n_seg = n_seg,
    sample_id = r,
    median_n_mut = median(seg_counts_joined$n_mut, na.rm = TRUE),
    median_seg_length = median(seg_counts_joined$len, na.rm = TRUE)
  )
  
  results_list[[i]] <- sample_info_CN_mut
  
  rm(
    x_after_inference,
    cna,
    seg_vec,
    seg_counts,
    seg_counts_joined,
    coords,
    from,
    to,
    info_sample
  )
  gc(verbose = FALSE)
}

info_parameters_PCAWG <- bind_rows(results_list)

saveRDS(info_parameters_PCAWG, paste0(output_dir, "/00_info_parameters_PCAWG.rds"))




########## get info of fitted samples before filtering and smoothing ###########
sample_dir_path <- list.files(dir_PCAWG_smooth)


results_list <- vector("list", length(sample_dir_path))

for (i in seq_along(sample_dir_path)) {
  
  r <- sample_dir_path[i]
  message(sprintf("Processing sample %d / %d: %s", i, length(sample_dir_path), r))
  
  x <- readRDS(paste0(dir_PCAWG_smooth, r))
  
  info_sample <- tibble(
    purity = x$metadata$purity, 
    median_DP = median(x$before_smoothing$snvs$DP)
  )
  
  cna <- x$before_smoothing$cna %>% 
    mutate(karyo = paste0(Major,":",minor)) %>% 
    filter(karyo %in% c("2:0", "2:2", "2:1"))
  
  n_seg <- nrow(cna)
  
  cna <- cna %>% mutate(len = to - from)
  
  seg_vec <- x$results_timing$data$input_data$seg_assignment
  
  seg_counts <- tibble(segment_id = seg_vec) %>%
    count(segment_id, name = "n_mut")
  
  seg_counts_joined <- seg_counts %>%
    left_join(cna, by = "segment_id") %>%
    mutate(density_mut_per_segment = n_mut / len)
  
  sample_info_CN_mut <- tibble(
    purity = info_sample$purity,
    median_DP = info_sample$median_DP,
    median_density_mut = median(seg_counts_joined$density_mut_per_segment, na.rm = TRUE),
    n_seg = n_seg,
    sample_id = r,
    median_n_mut = median(seg_counts_joined$n_mut, na.rm = TRUE),
    median_seg_length = median(seg_counts_joined$len, na.rm = TRUE)
  )
  
  results_list[[i]] <- sample_info_CN_mut
  
  rm(
    x_after_inference,
    cna,
    seg_vec,
    seg_counts,
    seg_counts_joined,
    coords,
    from,
    to,
    info_sample
  )
  gc(verbose = FALSE)
}

info_parameters_PCAWG <- bind_rows(results_list)

saveRDS(info_parameters_PCAWG, paste0(output_dir, "/00_info_parameters_PCAWG_before_smoothing.rds"))
