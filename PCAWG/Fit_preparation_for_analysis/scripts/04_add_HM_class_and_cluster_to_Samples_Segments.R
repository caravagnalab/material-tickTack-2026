rm(list=ls())
.libPaths("~/R/orfeo_R_4.4/")
library(dplyr)
library(ggplot2)
library(patchwork)

FIT_DIR <-"/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Fit/inference_results_5ncomponents/"
PLOT_DIR <-"/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Fit_preparation_for_analysis//plot/"
RES_FINAL_DIR <- "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Fit_preparation_for_analysis/results/"


info_fit_HM <- readRDS(paste0(RES_FINAL_DIR, "00_HM_class_info_fit.rds"))
Samples <- readRDS("/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Data/Samples.rds")

df1 <- info_fit_HM %>% dplyr::select(-c(class))
identical(df1, Samples)
saveRDS(info_fit_HM, "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Data/Samples.rds")



info_HM_segments <- readRDS(paste0(RES_FINAL_DIR, "00_HM_max_clusters.rds"))
Segments <- readRDS("/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Data/Segments.rds")


info_HM_segments = info_HM_segments %>% 
  group_by(sample) %>% 
  mutate(HM_cluster = ifelse(n_cna_events_per_timing_group == max(n_cna_events_per_timing_group), T, F)) %>%
  ungroup()

Segments %>% filter(wgd_status != "wgd") %>% 
  group_by(sample) %>% dplyr::mutate(clock_rank = dense_rank(clock_mean)) %>% ungroup() %>% 
  group_by(sample, clock_rank) %>% summarise(mean(clock_mean))


