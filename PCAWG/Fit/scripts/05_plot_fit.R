rm(list=ls())
.libPaths("~/R/orfeo_R_4.4/")
library(tickTack, lib.loc = "~/R/orfeo_R_4.4/")
library(dplyr)
library(patchwork)
library(ggplot2)

data_dir = "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Fit/data/"
inference_results_dir = "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Fit/inference_results_5ncomponents/"
output_dir = "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Fit/results/"
plot_dir = "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Fit/plot/"

info_smooth_samples <- readRDS(paste0(output_dir,"00_info_fit_smoot5Mb_1Mbml_15mm_0.4pi.RDS"))
fittable_samples = info_smooth_samples %>% filter(is_fittable_after_smoothing == TRUE) %>% pull(sample)
fitted_samples = tools::file_path_sans_ext(list.files(inference_results_dir))
path_sample <- paste0(inference_results_dir, fitted_samples, ".rds")

source("/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Fit/scripts/utils.R")
###### save table with the info for the timed samples in the main Data directory PCAWG/Data ######
info_fit <- readRDS(paste0("/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Data/Samples.rds"))
info_fit_best_K  <- readRDS(paste0("/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Data/Segments.rds"))

info_fit <- info_fit %>% left_join(info_fit_best_K %>% dplyr::select(sample, best_K) %>% distinct(),
                                   by = join_by(sample))

pdf(paste0(output_dir, "all_fit_VAFecdf_ncomp5.pdf"), width = 20, height = 8)

lapply(1:nrow(info_fit), function(i){
  
  fit <- info_fit[i, ]
  
  res <- readRDS(paste0(inference_results_dir,fit$sample,".rds"))
  res$ttype <- fit$ttype
  
  p <- tickTack::plot_cnaqc_choose_K(res,fit$best_K, add_VAF_ecdf = TRUE) + plot_annotation(
    caption = paste0("sample: ", fit$sample,"; ", fit$wgd_status),
    theme = theme(
      plot.caption = element_text(size = 12, hjust = 0) )) # & theme(text = element_text(size = 15))
  
  p <- p & theme(
    text = element_text(size = 14),
    # 1. Main Titles for each sub-plot
    plot.title = element_text(size = 14, face = "bold"),
    
    # 2. Axis Titles (the "x" and "y" labels)
    axis.title = element_text(size = 10),
    axis.title.x = element_text(size = 12), # Specific to X
    axis.title.y = element_text(size = 13), # Specific to Y
    
    # 3. Axis Text (the numbers/categories on the ticks)
    axis.text = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 10),
    
    # 4. Legend Text and Title
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    
  ) 
  
  print(p + theme(aspect.ratio = 1))
  ggsave(paste0(plot_dir,"/",fit$sample,"_ecdf.pdf"), plot = p, width = 13, height =10, units="in", dpi=300)
  
})

dev.off()