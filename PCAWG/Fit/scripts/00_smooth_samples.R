.libPaths("~/R/orfeo_R_4.4/")
library("dplyr")
library("ggplot2")
source("/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Fit/scripts/utils.R")
data_dir = "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Fit/data/"
output_dir = "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Fit/results/"
info_WGD_PCAWG <- readRDS(paste0(data_dir,"info_WGD_PCAWG.rds"))

info_df <- readRDS(paste0(data_dir,"info_cnaqc_samples.rds"))
info_df <- info_df %>% 
  select(-ploidy) %>% left_join(info_WGD_PCAWG %>% 
                                  select(- purity), by=join_by("sample"=="sample_id"))

n_samples = nrow(info_df %>% 
                   dplyr::select(sample) %>% 
                   unique())

n_tumour_types = nrow(info_df %>%
                       dplyr::select(ttype) %>% distinct ())

info_df %>%
  ggplot(aes(x = ttype)) +
  geom_bar(fill = "steelblue") +
  theme_minimal() +
  labs(x = "Cancer type", y = "Count")+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) 

# smoothinp parameter
max_distance_CNAqc = 5e6 #1e7 in GEL 
min_segment_length = 1e6 #1e7 in GEL 
min_purity = 0.4
min_mutations_number_per_segment = 15 


fittable_flags <- lapply(1:(nrow(info_df)), function(idx){
  cat(sprintf ("Completion = %.1f%%\n", idx / nrow(info_df) * 100))
  fit_x <- tryCatch(readRDS (paste0(info_df[idx, ]$cnaqc_path)), error = function(e) NULL)
  
  fit_x$mutations <- fit_x$snvs
  cnaqc_x <- CNAqc::smooth_segments(fit_x, maximum_distance = max_distance_CNAqc) 
  cnaqc_x$snvs <- cnaqc_x$mutations
  fittable <- is_fittable(
    cnaqc_x,
    min_segment_length = min_segment_length,
    min_purity = min_purity,
    min_mutations_number_per_segment = min_mutations_number_per_segment,
    possible_k = c("2:1",
                   "2:2",
                   "2:0")
  )
  return (fittable)
  
}) %>% unlist()

info_df$is_fittable_after_smoothing <- fittable_flags

###################### add IntoGene annotation ######################### 

tab <- read.csv(paste0(data_dir,"table_Intogene_annotation_mapped.tsv"), sep = "\t")
tab <- tab %>% rename(IntoGen_cancer_type = CANCER_TYPE,
                      IntoGen_cancer_name = CANCER_NAME)

info_df <- info_df %>% left_join(tab, by = c("cancer_type", "cancer_type_short"))

## missing annotation for 8 CMDI, 37 LIRI and 69 PBCA samples ##
info_df %>% filter(is.na(IntoGen_cancer_name),
                   is_fittable_after_smoothing) %>% dplyr::select(ttype) %>% table

saveRDS(info_df, paste0(output_dir,"/00_info_fit_smoot5Mb_1Mbml_15mm_0.4pi.RDS"))
