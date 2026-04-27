.libPaths("~/R/orfeo_R_4.4/")
library(tickTack)
library(ggplot2)
library(dplyr)

dir_PCAWG = "/orfeo/cephfs/scratch/cdslab/scocomello/data/clonal_analysis_PCAWG/"
output_dir = "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/simulations/data/"

sample_dir_path <- list.files(dir_PCAWG)
results_list <- vector("list", length(sample_dir_path))

for (i in seq_along(sample_dir_path)) {
  
  r <- sample_dir_path[i]
  message(sprintf("Processing sample %d / %d: %s", i, length(sample_dir_path), r))
  
  x <- readRDS(paste0(dir_PCAWG, r, "/fit.rds"))
  
  info_sample <- tibble(
    purity = x$purity,
    median_DP = median(x$snvs$DP)
  )
  
  cna <- x$cna
  n_seg_complex <- nrow(cna)
  cna <- cna %>%
    mutate(karyotype = paste0(Major, ":", minor)) %>%
    filter(karyotype %in% c("2:0", "2:1", "2:2"))
  n_seg <- nrow(cna)
  
  cna <- cna %>% mutate(len = to - from)
  
  snvs_info <- cna %>%
    select(segment_id, len, karyotype) %>%
    left_join(x$snvs, by = "segment_id") %>%
    group_by(segment_id) %>%
    mutate(
      n_mut = sum(!is.na(DP)),  # count only real SNV rows, not NA-filled ones
      mutation_density = n_mut / len
    ) %>%
    ungroup() %>% select(segment_id, len, n_mut, mutation_density) %>% distinct()
  
  # Add sample-level metadata to every row
  snvs_info <- snvs_info %>%
    mutate(
      sample_id     = r,
      purity        = info_sample$purity,
      median_DP     = info_sample$median_DP,
      n_seg         = n_seg,
      n_seg_complex = n_seg_complex
    )
  
  results_list[[i]] <- snvs_info
  
  rm(x, cna, snvs_info, info_sample)
  gc(verbose = FALSE)
}

info_parameters_PCAWG <- bind_rows(results_list)
saveRDS(info_parameters_PCAWG, paste0(output_dir, "/00_A_info_parameters_singleseg_PCAWG_before_smoothing.rds"))