rm(list=ls())
.libPaths("~/R/orfeo_R_4.4/")
library(dplyr)
library(ggplot2)
library(patchwork)

FIT_DIR <-"/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Fit/inference_results_5ncomponents/"
PLOT_DIR <-"/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Fit_preparation_for_analysis//plot/"
RES_FINAL_DIR <- "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Fit_preparation_for_analysis/results/"

Samples <- readRDS("/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Data/Samples.rds")
####### info fit of WGD samples #####
WGD_samples <- Samples %>% filter (wgd_status == "wgd")
info_fit_WGD <- Samples %>% filter(sample %in% (WGD_samples
                                                                                                      %>%pull(sample)))

###### info fit of classical + HM ########
other_samples <- Samples %>% filter (wgd_status == "no_wgd")
# Samples_name <- Samples %>% filter(sample %in% (other_samples %>% pull(sample)))


##### collect features to be used in finding HM ####
Segments <- readRDS("/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Data/Segments.rds")
Segments <- Segments %>% filter(sample %in% other_samples$sample)
info_fit <- Segments %>% tidyr::separate(segment_name, into = c("chr", "from", "to"), sep = "_", remove = FALSE) %>%
  mutate(from = as.integer(from)) %>% 
  mutate(to = as.integer(to))

total_bps_in_genome <- sum(CNAqc::chr_coordinates_hg19[1:22,]$length)

info_fit = info_fit %>% group_by(sample) %>% dplyr::mutate(clock_rank = dense_rank(clock_mean)) %>%
  ungroup %>% 
  dplyr::group_by(sample, clock_rank) %>% 
  dplyr::mutate(frac_genome_affected_per_timing_group = sum((to - from)/total_bps_in_genome)) %>% 
  dplyr::mutate(length_genome_affected_per_timing_group = sum(to - from)) %>% 
  dplyr::ungroup()


info_HM = lapply(unique (info_fit$sample), function (s) {
  print (s)
  info_fit_single <- info_fit %>% ungroup %>%
    dplyr::filter(sample == s) %>%
    dplyr::select(sample, 
                  ttype, 
                  ploidy, 
                  wgd_status, 
                  karyotype, 
                  clock_rank, clock_mean, 
                  chr, frac_genome_affected_per_timing_group, from, to) %>%
    dplyr::mutate(n_cna_per_sample = n(), n_clusters_per_sample = max(clock_rank)) 
  
  info_fit_single <- info_fit_single %>% dplyr::group_by(clock_rank) %>%
    dplyr::mutate(median_segment_length_per_timing_group = median(to-from)) %>%
    dplyr::mutate(mean_segment_length_per_timing_group = mean(to-from)) %>%
    dplyr::mutate(n_chr_affected_per_timing_group = length(unique (chr)) ) %>%
    dplyr::mutate(n_cna_events_per_timing_group = n()) %>%
    dplyr::select(sample, 
                  ttype, 
                  wgd_status, 
                  ploidy, 
                  n_cna_per_sample, 
                  n_clusters_per_sample, 
                  n_chr_affected_per_timing_group, 
                  clock_rank, clock_mean, 
                  frac_genome_affected_per_timing_group, 
                  n_cna_events_per_timing_group,
                  mean_segment_length_per_timing_group, 
                  mean_segment_length_per_timing_group) %>%
    dplyr::distinct () 
}) %>% do.call("bind_rows", .)

saveRDS(info_HM, paste0(RES_FINAL_DIR, "00_HM_max_clusters.rds"))

################################################################
################################################################



################################################################
########## select the group of segments of the clock which has the highest number of CN ###########
# info_HM <- readRDS(paste0(RES_FINAL_DIR, "06_A_info_HM.rds"))

info_HM_maxgroup = info_HM %>% 
  group_by(sample) %>% 
  filter(n_cna_events_per_timing_group == max(n_cna_events_per_timing_group)) %>%
  ungroup()
################################################################

######################## KMEANS CLUSTERING #####################
km_res = kmeans(info_HM_maxgroup %>%
                  dplyr::select( n_chr_affected_per_timing_group,
                          frac_genome_affected_per_timing_group,
                          n_cna_events_per_timing_group,
                          ploidy), 2)$cluster

info_HM_maxgroup$class = as.factor(km_res)
cluster_sizes <- table(info_HM_maxgroup$class)
classic_cluster <- names(which.max(cluster_sizes))
info_HM_maxgroup <- info_HM_maxgroup %>% 
  mutate(HM_class = factor(ifelse(class == classic_cluster, "Classic", "HM"), 
                           levels = c("Classic", "HM")))

table <- info_HM_maxgroup

png(paste0(PLOT_DIR,"00_genomic_features_byKmeansClass_HM_correlation.png"), width = 1200, height = 1000, res = 150)

if(T){
  grp = table$HM_class
  grp_levels = levels(grp)
  cols = c("steelblue","firebrick")
  names (cols) = grp_levels
  panel.scatter = function(x, y, group, ... ){
    # col = as. numeric(group)
    points(x,y,col=cols[group],pch=19,cex=0.8, ... )
  }
  panel.density = function(x, group, ... ){
    usr = par("usr"); on.exit(par(usr))
    par(usr=c(usr[1:2],0,1.5))
    for (g in levels(group) ){
      d=density(x[group == g])
      lines (d$x, d$y, col=cols[g],#col = which(levels(group) == g),
             lwd=2)
      
    }
  }
  
  op <- par(mar = c(4, 4, 2, 8))
  pairs(
    table[, c("n_chr_affected_per_timing_group",
              "frac_genome_affected_per_timing_group",
              "n_cna_events_per_timing_group",
              "ploidy")],
    diag.panel = function(x, ... ) panel.density(x, table$HM_class),
    lower.panel = function(x, y, ... ) panel.scatter(x, y, table$HM_class),
    upper.panel = function(x, y, ... ) panel.scatter(x, y, table$HM_class)
  )
  legend(
    "bottom",
    inset = c(-0.15,0),
    legend = levels(table$HM_class),
    col = cols,#seq_along(levels(table$is_HM) ),
    pch = 19,
    pt.cex=1,
    xpd=T,
    bty = "n"
  )
}

dev.off()


my_colors <- c("Classic" = "darkseagreen", "HM" = "brown4")


p1 = table %>% ggplot(aes(x= ploidy, color=HM_class, fill=HM_class, alpha = 0.3)) + 
  geom_density() + 
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme_minimal() 
theme_minimal() 
p2 = table %>% ggplot(aes(x = n_chr_affected_per_timing_group, color=HM_class, fill=HM_class, alpha = 0.3) ) +
  geom_density() +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme_minimal() 
p3 = table %>% ggplot(aes(x = frac_genome_affected_per_timing_group, color=HM_class, fill=HM_class, alpha = 0.3)) +
  geom_density() + 
  theme_minimal() + 
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme_minimal() 
p4 = table %>% ggplot(aes(x = n_cna_events_per_timing_group, color=HM_class, fill=HM_class, alpha = 0.3)) +
  geom_density() + 
  xlim(0,100) +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme_minimal() 

unique_plot <- (p1 + p2) / (p3 + p4) + 
  plot_layout(guides = "collect") + 
  plot_annotation(
    title = "Genomic Distribution by HM Class",
    subtitle = "Analysis of Ploidy, Chromosomes, and CNA Events",
    tag_levels = 'A' # Automatically adds A, B, C, D labels to each panel
  )

ggsave(paste0(PLOT_DIR,"00_genomic_features_byKmeansClass_HM.pdf"), plot = unique_plot, width = 9, height =5, units="in", dpi=300)
################################################################
################################################################




################################################################
######################## attach the classification in the two original tables #############


HM_samples <- info_HM_maxgroup %>% filter(HM_class == "HM")
Classic_samples <- info_HM_maxgroup %>% filter(HM_class == "Classic")
WGD_samples

Samples <- readRDS("/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Data/Samples.rds")

info_fit <- Samples %>% 
  mutate( class = case_when(sample %in% (HM_samples %>% 
                                           pull(sample) %>% 
                                           unique()) ~ "HM",
                            sample %in% (Classic_samples %>% 
                                           pull(sample) %>% 
                                           unique()) ~ "Classic",
                            sample %in% (WGD_samples %>% 
                                           pull(sample) %>% 
                                           unique()) ~ "WGD",
  ))

saveRDS(info_fit, paste0(RES_FINAL_DIR, "00_HM_class_info_fit.rds"))


# 
# 
# ################################################################
# ################# copy plots of fit in a pdf file for each of the three groups ##############
# ################################################################
# # if (!require("qpdf")) install.packages("qpdf")
# library(qpdf)
# 
# # HM
# source_folder <- PLOT_FIT_DIR
# dest_folder   <- "/orfeo/cephfs/scratch/cdslab/scocomello/material_tickTack/PCAWG/plot_final_per_class/HM"
# output_name   <- "Combined_HM_Plots.pdf"
# 
# if (!dir.exists(dest_folder)) dir.create(dest_folder, recursive = TRUE)
# 
# all_pdfs <- list.files(source_folder, pattern = "\\.pdf$", full.names = FALSE)
# matches  <- all_pdfs[tools::file_path_sans_ext(all_pdfs) %in% HM_samples$sample]
# 
# file_paths_source <- file.path(source_folder, matches)
# file_paths_dest   <- file.path(dest_folder, matches)
# 
# success <- file.copy(from = file_paths_source, to = dest_folder, overwrite = TRUE)
# 
# copied_files <- file_paths_dest[success]
# 
# qpdf::pdf_combine(
#   input = copied_files,
#   output = file.path(dest_folder, output_name)
# )
# 
# 
# # WGD
# source_folder <- PLOT_FIT_DIR
# dest_folder   <- "/orfeo/cephfs/scratch/cdslab/scocomello/material_tickTack/PCAWG/plot_final_per_class/WGD"
# output_name   <- "Combined_WGD_Plots.pdf"
# 
# if (!dir.exists(dest_folder)) dir.create(dest_folder, recursive = TRUE)
# 
# all_pdfs <- list.files(source_folder, pattern = "\\.pdf$", full.names = FALSE)
# matches  <- all_pdfs[tools::file_path_sans_ext(all_pdfs) %in% WGD_samples$sample]
# 
# file_paths_source <- file.path(source_folder, matches)
# file_paths_dest   <- file.path(dest_folder, matches)
# 
# success <- file.copy(from = file_paths_source, to = dest_folder, overwrite = TRUE)
# 
# copied_files <- file_paths_dest[success]
# 
# qpdf::pdf_combine(
#   input = copied_files,
#   output = file.path(dest_folder, output_name)
# )
# 
# 
# 
# # Classic
# source_folder <- PLOT_FIT_DIR
# dest_folder   <- "/orfeo/cephfs/scratch/cdslab/scocomello/material_tickTack/PCAWG/plot_final_per_class/Classic"
# output_name   <- "Combined_Classic_Plots.pdf"
# 
# if (!dir.exists(dest_folder)) dir.create(dest_folder, recursive = TRUE)
# 
# all_pdfs <- list.files(source_folder, pattern = "\\.pdf$", full.names = FALSE)
# matches  <- all_pdfs[tools::file_path_sans_ext(all_pdfs) %in% Classic_samples$sample]
# 
# file_paths_source <- file.path(source_folder, matches)
# file_paths_dest   <- file.path(dest_folder, matches)
# 
# success <- file.copy(from = file_paths_source, to = dest_folder, overwrite = TRUE)
# 
# copied_files <- file_paths_dest[success]
# 
# qpdf::pdf_combine(
#   input = copied_files,
#   output = file.path(dest_folder, output_name)
# )
# 
# 
# ################################################################
# ################################################################
# 
# 
