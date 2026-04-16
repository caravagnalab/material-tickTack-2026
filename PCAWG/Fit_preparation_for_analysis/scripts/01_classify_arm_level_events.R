rm(list=ls())
.libPaths("~/R/orfeo_R_4.4/")
library(dplyr)
library(ggplot2)
library(parallel)
library(tibble)
library(stringr)

main_path <-"/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Fit/inference_results_5ncomponents/"
inference_output_dir <-"/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Fit/results/"
summary_samples <- readRDS(paste0(inference_output_dir, "04_parsed_inference_results_5ncomponents.rds"))
output_dir <- "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Fit_preparation_for_analysis/results/"





chromosomes_table = CNAqc:::get_reference('hg19') %>%
  mutate(p_arm_length = centromerStart - from, q_arm_length = to - centromerEnd)

# summary_samples = summary_samples %>% 
#   dplyr::select(sample, IntoGen_cancer_type) %>% 
#   dplyr::distinct()


# Classes: whole chromosome, arm-event (p-q), focal event
classify_segment = function(s, data){
  #print(segment)
  segment = data[s,]
  #print(segment)
  if (segment$chr %in% chromosomes_table$chr){
    targeted_chr =  chromosomes_table %>% filter(chr == segment$chr) 
    return(segment %>% 
      dplyr::mutate(len= to-from) %>% 
      dplyr::mutate(
        class = case_when(
          ## Amplifications 
          # whole chromosome : the CNA involves more than 80% of the total length of the chromosome
          len >= targeted_chr$length * .8 ~ 'whole_chromosome',
          # p-arm : The CNA is on the p-arm and involves more than 80% of its length
          (from + targeted_chr$from < targeted_chr$centromerStart ) & (len >= targeted_chr$p_arm_length * .8) ~ 'p_arm',
          # q-arm : The CNA is on the q-arm and involves more than 80% of its length
          (from + targeted_chr$from > targeted_chr$centromerEnd)  & (len >= targeted_chr$q_arm_length * .8) ~ 'q_arm',
          # focal : The CNA is on one arm and involves less than 80% of its length
          .default = 'focal'
        )
      ))
  }else{
    segment$class = NA
    return(segment)
  }
  
}

##### classify the events at arm level excluding the segments in X and Y chromosomes ####

sample_paths = list.files(main_path)
sample_paths = na.omit(sample_paths) 

# sp = sample_paths[1]
# sp = "09508a0d-ebe0-4fa1-b7b2-1710814181cd.rds"

df_new = parallel::mclapply(sample_paths, function(sp) {
  print(which(sp == sample_paths) / length(sample_paths) * 100)
  print(sp)
  
  sample_name = str_replace(sp, ".rds", "")
  if (!(sample_name %in% summary_samples$sample)) return(NULL)
  summary_sample_info <- summary_samples %>% filter(sample == sample_name)

  res = summary_sample_info %>% rowwise() %>% 
    tidyr::separate(segment_name, sep = "_", into = c("chr", "from", "to"), convert = TRUE) %>% 
    dplyr::mutate(rank = dense_rank(clock_mean)) %>% 
    dplyr::filter(chr != "chrNA")
  
  res_annotated = lapply(1:nrow(res), classify_segment, data = res) %>% do.call("bind_rows", .)
  
  
  dplyr::bind_cols(res_annotated%>% select(class), summary_samples %>% dplyr::filter(sample == sample_name,
                                                                                     chr != "chrNA"))
}, mc.cores = 1)

df_new = bind_rows(df_new)
colnames(df_new)[1] <- "arm_level_class"
# colnames(df_new)[16] <- "class"

saveRDS(df_new, paste0(output_dir,"01_arm_level_events_table_BIC.rds"))

df_new <- readRDS(paste0(output_dir,"01_arm_level_events_table_BIC.rds"))
saveRDS(df_new, paste0("/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Data/Segments.rds"))
