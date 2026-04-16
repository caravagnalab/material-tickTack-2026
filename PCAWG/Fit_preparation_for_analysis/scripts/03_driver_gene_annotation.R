
#################################### original version ###############################

rm(list=ls())
.libPaths("~/R/orfeo_R_4.4/")
library(dplyr)
library(ggplot2)
library(mclust)
library(cluster)
library(biomaRt)

main_path <-"/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Fit/inference_results_5ncomponents/"
data_dir <-"/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Fit_preparation_for_analysis/data/"
inference_output_dir <-"/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Fit/results/"
output_dir <- "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Fit_preparation_for_analysis/results/"

summary_segments <- readRDS(paste0(output_dir,"01_arm_level_events_table_BIC.rds"))
info_fit <- readRDS("/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Data/Samples.rds")
PCAWG_drivers <- read.csv(paste0(data_dir, "/TableS3_panorama_driver_mutations_ICGC_samples.public.tsv"), sep = "\t")
driver_list_IntoGene_CI <- readRDS("/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Fit_preparation_for_analysis/data/driver_list_IntoGene_CI.rds")
gene_coords_IntoGene_CI <- readRDS("/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Fit_preparation_for_analysis/data/gene_coords_IntoGene_CI.rds")

source("/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Fit_preparation_for_analysis/scripts/utils.R")


summary_gene_annotation <- lapply(1:nrow(info_fit), function(idx) {

  cat(sprintf("Processing sample %d / %d (%.1f%%)\n", idx, nrow(info_fit), idx/nrow(info_fit)*100))

  # sample_name <- "03c3c692-8a86-4843-85ae-e045f0fa6f88"
  # info_single <- info_fit%>%filter(sample == sample_name)

  info_single <- info_fit[idx, ]
  summary_segments_single <- summary_segments %>% filter(sample == info_single$sample)

  # sample_info <- info_single %>% dplyr:: select(best_K, CNAs, mean_mut_per_segment, median_mut_per_segment,
  # is_WGD, ttype, reference_genome, ploidy, purity, cancer_type, cancer_type_short, class)

  # select the smoothed version of the fit otherwise the merhe with the timed segments would not be possible
  single_fit <- readRDS(paste0(main_path,info_single$sample,".rds"))

  snvs_prepared <- single_fit$mutations

  # PCAWG DRIVERS
  PCAWG_drivers_SNV <- PCAWG_drivers %>%
    filter(sample_id == info_single$sample) %>%
    filter(pos != "x") %>%
    mutate(mutation_id = paste0("chr",chr,":",pos))

  # PCAWG CNA DRIVERS
  PCAWG_drivers_CNA <- PCAWG_drivers %>%
    filter(sample_id == info_single$sample) %>%
    filter(top_category == "CNA") %>%
    mutate(mutation_id = paste0("chr",chr,":",pos))

  # Geni mutati, da controllare molteplicità della mutazione
  driver_snvs = snvs_prepared %>% mutate(mutation_id = paste0(chr,":",from),
                                         gene = if (!"gene" %in% colnames(.)) NA_character_ else gene) %>%
    filter(mutation_id %in% PCAWG_drivers_SNV$mutation_id) %>%
    dplyr::select(segment_id, chr, from, gene, karyotype, VAF, NV, DP) %>%
    mutate( sample_id = info_single$sample, mutatation_status = "M")
  # uncomment 169- 186 to account for not annotated drivers
  # , mutation_call = "in PCAWG driver annotation")

  # # PCAWG DRIVERS not in snvs table: Geni mutati, da controllare molteplicità della mutazione
  # driver_snvs_mut_id = snvs_prepared %>% mutate(mutation_id = paste0(chr,":",from),
  #                                        gene = if (!"gene" %in% colnames(.)) NA_character_ else gene) %>%
  #   filter(mutation_id %in% PCAWG_drivers_SNV$mutation_id)
  # driver_snvs_missing_mut_id = setdiff(PCAWG_drivers_SNV %>% pull(mutation_id) , driver_snvs_mut_id %>% pull(mutation_id) )
  # if( length(driver_snvs_missing_mut_id) > 0  ){
  #   driver_snvs_missing <- PCAWG_drivers_SNV %>%
  #     filter(mutation_id %in% driver_snvs_missing_mut_id) %>%
  #     dplyr::mutate(segment_id = NA_character_) %>%
  #     mutate(from = as.numeric(pos)) %>%
  #     dplyr::select( chr, from, gene, segment_id) %>%
  #     dplyr:: mutate(karyotype = NA_character_, VAF = NA, NV = NA, DP = NA) %>%
  #     mutate( sample_id = info_single$sample, mutatation_status = "M", mutation_call = "in PCAWG driver annotation but missing")
  #
  #   driver_snvs <- bind_rows(driver_snvs, driver_snvs_missing)
  # }


  driver_cna = gene_coords_IntoGene_CI %>%
    filter(hgnc_symbol %in% PCAWG_drivers_CNA$gene) %>%
    mutate(mutation_status = "WT")
  # uncomment next line to account for not annotated drivers
  # , mutation_call = "in PCAWG")  # trovare le coordinate e dinferire karyo


  #driver_cna_final = data.frame()
  if(nrow(driver_cna) > 0){
    driver_cna_final = lapply(1:nrow(driver_cna), function(d){
      # print(d)
      cna_gene <- single_fit$cna %>% dplyr::select(chr, from, to, segment_id, Major, minor) %>% filter(chr == paste0("chr", driver_cna$chromosome_name[d]),
                                                                                                       from <= driver_cna$start_position[d],
                                                                                                       to >= driver_cna$end_position[d])
      if(nrow(cna_gene) > 0){
        cbind(
          driver_cna[d,],
          cna_gene,
          PCAWG_drivers_CNA$sample_id[d]
        )
      } else {
        cbind(
          driver_cna[d, ][0, ] ,
          cna_gene,
          PCAWG_drivers_CNA$sample_id[d][0]
        )
      }

    }) %>% Reduce(rbind, .)

  } else {
    driver_cna_final <- data.frame(
      ensembl_gene_id = character(),
      hgnc_symbol = character(),
      chromosome_name = character(),
      start_position = numeric(),
      end_position = numeric(),
      mutation_status = character(),
      mutatation_call = character(),
      chr = character(),
      from = numeric(),
      to = numeric(),
      segment_id = character(),
      Major = numeric(),
      minor = numeric(),
      sample_id = character(),
      stringsAsFactors = FALSE
    )
  }


  # CI
  ci_positions = gene_coords_IntoGene_CI %>% filter(hgnc_symbol %in% driver_list_IntoGene_CI$CI_genes)
  # CI mutations (exclude tail)
  ci_mutations = lapply(1:nrow(ci_positions), function(cig){
    annotations = snvs_prepared %>% filter(chr == paste0("chr", ci_positions[cig,]$chr),
                                           as.numeric(from) >= as.numeric(ci_positions[cig,]$start_position),
                                           as.numeric(to) <= as.numeric(ci_positions[cig,]$end_position)) %>%
      dplyr::select(segment_id, chr, from, karyotype, sample, VAF, NV, DP)
    if (nrow(annotations)>0){
      gene = ci_positions$hgnc_symbol[cig]
      return(
        cbind(
          snvs_prepared %>% filter(chr == paste0("chr", ci_positions[cig,]$chr),
                                   as.numeric(from) >= as.numeric(ci_positions[cig,]$start_position),
                                   as.numeric(to) <= as.numeric(ci_positions[cig,]$end_position)) %>%
            dplyr::select(segment_id, chr, from, karyotype, sample, VAF, NV, DP),
          gene
        ))
    } else {
      ci_mutations <- data.frame(
        segment_id = character(),
        chr = character(),
        from = numeric(),
        karyotype = character(),
        sample = character(),
        VAF = numeric(),
        NV = numeric(),
        DP = numeric(),
        gene = character(),
        stringsAsFactors = FALSE
      )
    }
  }) %>% Reduce(rbind, .)

  ## Geni WT amplificati
  geni_amplificati_su_seg_timati = lapply(1:nrow(summary_segments_single), function(s){
    # print(s)
    from_s= strsplit(summary_segments_single$segment_name[s],"_")[[1]][2]
    to_s= strsplit(summary_segments_single$segment_name[s],"_")[[1]][3]

    # seg_gene <- (gene_coords %>% mutate(chr = paste0("chr", chromosome_name)) %>% filter(
    #   chr %in% summary_segments_single$chr[s],
    #   start_position >= from_s,
    #   end_position <= to_s) %>% dplyr::select(hgnc_symbol))
    #
    seg_gene <- (gene_coords_IntoGene_CI %>% mutate(chr = paste0("chr", chromosome_name)) %>% filter(
      chr %in% summary_segments_single$chr[s],
      as.numeric(start_position) >= as.numeric(from_s),
      as.numeric(end_position) <= as.numeric(to_s)) %>% dplyr::select(hgnc_symbol))

    if(nrow(seg_gene) > 0) {
      cbind(
        summary_segments_single[s,],
        seg_gene
      )
    } else {
      cbind(
        summary_segments_single[s,][0, ] ,
        seg_gene
      )
    }

  }) %>% Reduce(rbind, .) #%>% dpl
  geni_amplificati_su_seg_timati = geni_amplificati_su_seg_timati %>% dplyr::select(segment_name, karyotype, hgnc_symbol, sample) %>% mutate(mutation_status = "WT")


  ## TABELLA TOTALE
  driver_snvs # annotati per PCAWG, amplificati o no, da annotare la multiplicity per amplificazioni timate
  driver_cna_final # annotati per PCAWG, da annotare la multiplicity per amplificazioni timate -> unire con info della mutazione
  ci_mutations  # mutazioni di geni CI che possono essere o non essere in CNA
  geni_amplificati_su_seg_timati

  if(nrow(driver_snvs) > 0){
    driver_snvs <- filter_tail_mutations(driver_snvs, info_single$purity, exclude = FALSE)
  }
  if (nrow(ci_mutations) > 0){
    ci_mutations <- filter_tail_mutations(ci_mutations, info_single$purity, exclude = TRUE)
  }

  all_genes = c(driver_snvs$gene, driver_cna_final$hgnc_symbol,
                ci_mutations$gene, geni_amplificati_su_seg_timati$hgnc_symbol) %>% unique()


  drivers_df = lapply(1:length(all_genes), function(g){
    # print(g)
    final_df <- tibble::tibble(
      segment_id = character(),
      gene = character(),
      karyotype = character(),
      sample_id = character(),
      NV = numeric(),
      DP = numeric(),
      mutatation_status = character(),
      mult_estimate = numeric(),
      timed = logical()
    )

    is_mutated = nrow(driver_snvs %>% filter(gene == all_genes[g])) > 0
    is_driver_CNA = nrow(driver_cna_final %>% filter(hgnc_symbol == all_genes[g])) > 0
    is_CI_mutated = nrow(ci_mutations %>% filter(gene == all_genes[g])) > 0
    is_on_timed_segment = nrow(geni_amplificati_su_seg_timati %>% filter(hgnc_symbol == all_genes[g])) > 0 # check that there are some cases with M in P ann but mmissing + in timed segments

    if (is_mutated) {
      final_df = driver_snvs %>% filter(gene == all_genes[g]) %>%
        dplyr::select(segment_id, gene, karyotype, sample_id, NV, DP, mutatation_status, mult_estimate)
      #   # uncomment next line to account for not annotated drivers
      # , mutation_call)
      if (is_on_timed_segment) {
        final_df$timed = T
        # uncomment next lines 350-357 to account for not annotated drivers
        #   if (any(final_df$mutation_call == "in PCAWG driver annotation but missing")){
        #   info_seg_time = geni_amplificati_su_seg_timati %>% filter(hgnc_symbol == all_genes[g]) %>% rowwise %>%
        #     mutate(segment_id = paste0(strsplit(segment_name, "_")[[1]][1], ":",
        #                                strsplit(segment_name, "_")[[1]][2], ":",
        #                                strsplit(segment_name, "_")[[1]][3], ":",karyotype,":1"))
        #   final_df <- final_df %>% mutate(karyotype = ifelse( mutation_call == "in PCAWG driver annotation but missing", info_seg_time$karyotype, karyotype),
        #                                   segment_id = ifelse( mutation_call == "in PCAWG driver annotation but missing", info_seg_time$segment_id, segment_id))
        # }
      }else{ final_df$timed = F}
    }

    if (!is_mutated & is_driver_CNA){
      final_df = driver_cna_final %>% filter(hgnc_symbol == all_genes[g]) %>% mutate(gene = hgnc_symbol, karyotype = paste0(Major, ":", minor), sample_id = `PCAWG_drivers_CNA$sample_id[d]`,
                                                                                     NV = NA, DP = NA, mutatation_status = "CNA_driver") %>%
        dplyr::select(segment_id, gene, karyotype, sample_id, NV, DP, mutatation_status) %>%
        mutate(mult_estimate = NA)
      #   # uncomment next line to account for not annotated drivers
      # , mutation_call = NA)
      if (is_on_timed_segment) final_df$timed = T else final_df$timed = F
    }

    if (is_CI_mutated & !is_mutated & !is_driver_CNA){
      final_df = ci_mutations %>% filter(gene == all_genes[g]) %>% mutate(mutatation_status = "CI_M") %>%
        mutate(sample_id = sample) %>%
        dplyr::select(segment_id, gene, karyotype, sample_id, NV, DP, mutatation_status, mult_estimate)
      #   # uncomment next 2 lines to account for not annotated drivers
      # %>%
      #   mutate(mutation_call = NA)
      if (is_on_timed_segment) final_df$timed = T else final_df$timed = F
    }

    if (is_on_timed_segment & !is_CI_mutated & !is_mutated & !is_driver_CNA){
      final_df = geni_amplificati_su_seg_timati %>% filter(hgnc_symbol == all_genes[g]) %>%
        mutate(segment_id = paste0(strsplit(segment_name, "_")[[1]][1], ":",
                                   strsplit(segment_name, "_")[[1]][2], ":",
                                   strsplit(segment_name, "_")[[1]][3], ":",karyotype,":1"),
               gene = hgnc_symbol, sample_id = sample, NV=NA, DP=NA, timed =T, mutatation_status=mutation_status) %>%
        dplyr::select(segment_id, gene, karyotype, sample_id, NV, DP, mutatation_status, timed) %>%
        mutate(mult_estimate = NA)
      #   # uncomment next line to account for not annotated drivers
      # , mutation_call = NA)
    }

    if (nrow(final_df)==0) {
      final_df <- tibble::tibble(
        segment_id = character(),
        gene = character(),
        karyotype = character(),
        sample_id = character(),
        NV = numeric(),
        DP = numeric(),
        mutatation_status = character(),
        mult_estimate = numeric(),
        timed = logical()
        #   # uncomment next line to account for not annotated drivers
        # ,
        # mutation_call = character()
      )
    }
    #print(colnames(final_df))
    final_df %>% mutate(ttype = info_single$ttype,
                        cancer_type_short = info_single$cancer_type_short,
                        cancer_type = info_single$cancer_type)
    return(final_df)

  }) %>% Reduce(rbind,.)


}) %>% Reduce(rbind,.)


# merge with timed segments info and annotate class of gene
gene_annotation <- summary_gene_annotation %>% rename(mutation_status = mutatation_status )

summary_segments_merge <- summary_segments %>% dplyr::select(-c(karyotype))

colnames(summary_segments_merge)
colnames(gene_annotation)

gene_annotation <- gene_annotation %>% rowwise() %>% mutate(segment_name = paste0(strsplit(segment_id,split=":")[[1]][1],"_",
                                                                                  strsplit(segment_id,split=":")[[1]][2],"_",
                                                                                  strsplit(segment_id,split=":")[[1]][3]),
                                                            sample = sample_id
) %>% left_join(summary_segments_merge, by = c("segment_name", "sample"))


saveRDS(gene_annotation, paste0(output_dir,"CHECK_tmp_03_gene_annotation_table.rds"))

####################################################################################

# filter per ttype
Intogen_drivers <- read.csv(paste0(data_dir, "2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv"), sep = "\t")

gene_annotation <- readRDS(paste0(output_dir,"CHECK_tmp_03_gene_annotation_table.rds"))

drivers_Int = Intogen_drivers %>% dplyr::select(SYMBOL, CANCER_TYPE)
ttypes = drivers_Int$CANCER_TYPE %>% unique()

annotation_new = lapply(ttypes, function(tt){
  genes_tt = drivers_Int %>% filter(CANCER_TYPE == tt) %>% pull(SYMBOL) %>% unique()

  gene_annotation %>%
    filter(IntoGen_cancer_type == tt) %>%
    mutate(to_keep = ifelse( (gene %in% genes_tt) | (mutation_status!="WT"), T, F)) %>%
    filter(to_keep==T) %>% dplyr::select(!to_keep)
}) %>% Reduce(rbind, .)

saveRDS(annotation_new, paste0(output_dir,"03_gene_annotation_table.rds"))


gene_annotation_table <- readRDS("/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Fit_preparation_for_analysis/results/03_gene_annotation_table.rds")
saveRDS(gene_annotation_table, paste0("/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Data/Drivers.rds"))



# 
# ## get table for BTM
# gene_annotation_new_notation <- readRDS(paste0(RES_FINAL_DIR,"06_Cb_gene_annotation_table.rds"))
# 
# BTM_table <- gene_annotation_new_notation %>% filter(timed == TRUE)
# saveRDS(BTM_table, paste0(RES_FINAL_DIR,"06_Cb_BTM_table.rds"))


# 
# 
# #################################### speed up version ###############################
# rm(list=ls())
# .libPaths("~/R/orfeo_R_4.4/")
# library(dplyr)
# library(ggplot2)
# library(mclust)
# library(cluster)
# library(biomaRt)
# library(data.table)
# library(future.apply)
# 
# main_path  <- "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Fit/inference_results_5ncomponents/"
# data_dir   <- "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Fit_preparation_for_analysis/data/"
# output_dir <- "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Fit_preparation_for_analysis/results/"
# 
# summary_segments        <- readRDS(paste0(output_dir, "01_arm_level_events_table_BIC.rds"))
# info_fit                <- readRDS("/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Data/Samples.rds")
# PCAWG_drivers           <- read.csv(paste0(data_dir, "TableS3_panorama_driver_mutations_ICGC_samples.public.tsv"), sep = "\t")
# driver_list_IntoGene_CI <- readRDS(paste0(data_dir, "driver_list_IntoGene_CI.rds"))
# gene_coords_IntoGene_CI <- readRDS(paste0(data_dir, "gene_coords_IntoGene_CI.rds"))
# 
# source("/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Fit_preparation_for_analysis/scripts/utils.R")
# 
# # ── Pre-processing: split large tables by sample before the loop ──────────────
# message("Pre-processing: splitting tables by sample...")
# PCAWG_drivers_by_sample    <- split(PCAWG_drivers,    PCAWG_drivers$sample_id)
# summary_segments_by_sample <- split(summary_segments, summary_segments$sample)
# 
# # Pre-compute CI gene set once
# ci_gene_set <- driver_list_IntoGene_CI$CI_genes
# 
# # Convert gene coords to data.table for fast range joins
# gc_dt <- as.data.table(gene_coords_IntoGene_CI)
# gc_dt[, chr_full := paste0("chr", chromosome_name)]
# 
# # ── Set up parallelism ────────────────────────────────────────────────────────
# n_cores <- 1
# message(sprintf("Running with %d cores...", n_cores))
# plan(multisession, workers = n_cores)
# 
# # ── Per-sample function ───────────────────────────────────────────────────────
# process_sample <- function(idx) {
#   tryCatch({
#     info_single <- info_fit[idx, ]
#     sname       <- info_single$sample
#     message(sprintf("Processing sample %d / %d : %s", idx, nrow(info_fit), sname))
# 
#     # O(1) lookup instead of filtering full tables
#     summary_segments_single <- summary_segments_by_sample[[sname]]
#     drivers_sample          <- PCAWG_drivers_by_sample[[sname]]
# 
#     single_fit    <- readRDS(paste0(main_path, sname, ".rds"))
#     snvs_prepared <- single_fit$mutations
# 
#     # Add mutation_id and gene column once upfront
#     snvs_prepared <- snvs_prepared %>%
#       mutate(
#         mutation_id = paste0(chr, ":", from),
#         gene        = if (!"gene" %in% colnames(.)) NA_character_ else gene
#       )
# 
#     # ── PCAWG drivers for this sample ──────────────────────────────────────────
#     if (is.null(drivers_sample) || nrow(drivers_sample) == 0) {
#       PCAWG_drivers_SNV <- data.frame()
#       PCAWG_drivers_CNA <- data.frame()
#     } else {
#       PCAWG_drivers_SNV <- drivers_sample %>%
#         filter(pos != "x") %>%
#         mutate(mutation_id = paste0("chr", chr, ":", pos))
#       PCAWG_drivers_CNA <- drivers_sample %>%
#         filter(top_category == "CNA") %>%
#         mutate(mutation_id = paste0("chr", chr, ":", pos))
#     }
# 
#     # ── Driver SNVs ────────────────────────────────────────────────────────────
#     if (nrow(PCAWG_drivers_SNV) > 0) {
#       driver_snvs <- snvs_prepared %>%
#         filter(mutation_id %in% PCAWG_drivers_SNV$mutation_id) %>%
#         dplyr::select(segment_id, chr, from, gene, karyotype, VAF, NV, DP) %>%
#         mutate(sample_id = sname, mutatation_status = "M")
#     } else {
#       driver_snvs <- data.frame()
#     }
# 
#     # ── Driver CNAs ────────────────────────────────────────────────────────────
#     driver_cna <- gene_coords_IntoGene_CI %>%
#       filter(hgnc_symbol %in% PCAWG_drivers_CNA$gene) %>%
#       mutate(mutation_status = "WT")
# 
#     if (nrow(driver_cna) > 0 && nrow(single_fit$cna) > 0) {
#       cna_dt  <- as.data.table(single_fit$cna)[, .(chr, from, to, segment_id, Major, minor)]
#       dcna_dt <- as.data.table(driver_cna)
#       dcna_dt[, chr_cna := paste0("chr", chromosome_name)]
#       dcna_dt[, start_position := as.numeric(start_position)]
#       dcna_dt[, end_position   := as.numeric(end_position)]
#       cna_dt[,  from := as.numeric(from)]
#       cna_dt[,  to   := as.numeric(to)]
# 
#       # Non-equi join: gene overlaps CNA segment
#       driver_cna_final <- cna_dt[dcna_dt,
#                                  on = .(chr = chr_cna, from <= start_position, to >= end_position),
#                                  nomatch = 0
#       ] %>% as.data.frame()
# 
#       if (nrow(driver_cna_final) > 0) {
#         driver_cna_final$sample_id <- PCAWG_drivers_CNA$sample_id[
#           match(driver_cna_final$hgnc_symbol, PCAWG_drivers_CNA$gene)
#         ]
#       }
#     } else {
#       driver_cna_final <- data.frame(
#         ensembl_gene_id = character(), hgnc_symbol = character(),
#         chromosome_name = character(), start_position = numeric(),
#         end_position    = numeric(),  mutation_status = character(),
#         chr = character(), from = numeric(), to = numeric(),
#         segment_id = character(), Major = numeric(), minor = numeric(),
#         sample_id  = character(), stringsAsFactors = FALSE
#       )
#     }
# 
#     # ── CI mutations ───────────────────────────────────────────────────────────
#     ci_positions <- gene_coords_IntoGene_CI %>%
#       filter(hgnc_symbol %in% ci_gene_set)
# 
#     if (nrow(ci_positions) > 0 && nrow(snvs_prepared) > 0) {
#       snvs_dt <- as.data.table(snvs_prepared)[,
#                                               .(segment_id, chr, from, to, karyotype, sample, VAF, NV, DP)
#       ]
#       snvs_dt[, from := as.numeric(from)]
#       snvs_dt[, to   := as.numeric(to)]
# 
#       cip_dt <- as.data.table(ci_positions)
#       cip_dt[, chr_ci          := paste0("chr", chromosome_name)]
#       cip_dt[, start_position  := as.numeric(start_position)]
#       cip_dt[, end_position    := as.numeric(end_position)]
# 
#       ci_joined <- snvs_dt[cip_dt,
#                            on = .(chr = chr_ci, from >= start_position, to <= end_position),
#                            nomatch = 0
#       ] %>% as.data.frame()
# 
#       if (nrow(ci_joined) > 0) {
#         ci_mutations <- ci_joined %>% rename(gene = hgnc_symbol)
#       } else {
#         ci_mutations <- data.frame()
#       }
#     } else {
#       ci_mutations <- data.frame()
#     }
# 
#     # ── Genes on timed segments ────────────────────────────────────────────────
#     if (!is.null(summary_segments_single) && nrow(summary_segments_single) > 0) {
#       seg_parts <- do.call(rbind, strsplit(summary_segments_single$segment_name, "_"))
#       segs_dt <- data.table(
#         chr          = summary_segments_single$chr,
#         from_s       = as.numeric(seg_parts[, 2]),
#         to_s         = as.numeric(seg_parts[, 3]),
#         segment_name = summary_segments_single$segment_name,
#         karyotype    = summary_segments_single$karyotype,
#         sample       = summary_segments_single$sample
#       )
# 
#       gc_local <- copy(gc_dt)  # avoid modifying shared data.table
#       joined <- gc_local[segs_dt,
#                          on = .(chr_full = chr, start_position >= from_s, end_position <= to_s),
#                          nomatch = 0
#       ] %>% as.data.frame()
# 
#       if (nrow(joined) > 0) {
#         geni_amplificati_su_seg_timati <- joined %>%
#           dplyr::select(segment_name, karyotype, hgnc_symbol, sample) %>%
#           mutate(mutation_status = "WT")
#       } else {
#         geni_amplificati_su_seg_timati <- data.frame()
#       }
#     } else {
#       geni_amplificati_su_seg_timati <- data.frame()
#     }
# 
#     # ── Filter tail mutations ──────────────────────────────────────────────────
#     if (nrow(driver_snvs) > 0)
#       driver_snvs  <- filter_tail_mutations(driver_snvs,  info_single$purity, exclude = FALSE)
#     if (nrow(ci_mutations) > 0)
#       ci_mutations <- filter_tail_mutations(ci_mutations, info_single$purity, exclude = TRUE)
# 
#     # ── Build per-gene summary ─────────────────────────────────────────────────
#     all_genes <- unique(c(
#       if (nrow(driver_snvs) > 0)                   driver_snvs$gene                        else NULL,
#       if (nrow(driver_cna_final) > 0)               driver_cna_final$hgnc_symbol            else NULL,
#       if (nrow(ci_mutations) > 0)                   ci_mutations$gene                       else NULL,
#       if (nrow(geni_amplificati_su_seg_timati) > 0) geni_amplificati_su_seg_timati$hgnc_symbol else NULL
#     ))
# 
#     if (length(all_genes) == 0) return(NULL)
# 
#     # Pre-split sub-tables once to avoid re-filtering inside gene loop
#     ds_by_gene   <- if (nrow(driver_snvs) > 0)                   split(driver_snvs,                   driver_snvs$gene)                   else list()
#     dcna_by_gene <- if (nrow(driver_cna_final) > 0)               split(driver_cna_final,               driver_cna_final$hgnc_symbol)       else list()
#     ci_by_gene   <- if (nrow(ci_mutations) > 0)                   split(ci_mutations,                   ci_mutations$gene)                  else list()
#     amp_by_gene  <- if (nrow(geni_amplificati_su_seg_timati) > 0) split(geni_amplificati_su_seg_timati, geni_amplificati_su_seg_timati$hgnc_symbol) else list()
# 
#     drivers_df <- lapply(all_genes, function(g) {
#       is_mutated      <- !is.null(ds_by_gene[[g]])   && nrow(ds_by_gene[[g]])   > 0
#       is_driver_CNA   <- !is.null(dcna_by_gene[[g]]) && nrow(dcna_by_gene[[g]]) > 0
#       is_CI_mutated   <- !is.null(ci_by_gene[[g]])   && nrow(ci_by_gene[[g]])   > 0
#       is_on_timed_seg <- !is.null(amp_by_gene[[g]])  && nrow(amp_by_gene[[g]])  > 0
# 
#       if (is_mutated) {
#         final_df <- ds_by_gene[[g]] %>%
#           dplyr::select(segment_id, gene, karyotype, sample_id, NV, DP, mutatation_status, mult_estimate) %>%
#           mutate(timed = is_on_timed_seg)
# 
#       } else if (is_driver_CNA) {
#         final_df <- dcna_by_gene[[g]] %>%
#           mutate(
#             gene             = hgnc_symbol,
#             karyotype        = paste0(Major, ":", minor),
#             NV               = NA_real_,
#             DP               = NA_real_,
#             mutatation_status = "CNA_driver",
#             mult_estimate    = NA_real_,
#             timed            = is_on_timed_seg
#           ) %>%
#           dplyr::select(segment_id, gene, karyotype, sample_id, NV, DP, mutatation_status, mult_estimate, timed)
# 
#       } else if (is_CI_mutated) {
#         final_df <- ci_by_gene[[g]] %>%
#           mutate(
#             mutatation_status = "CI_M",
#             sample_id         = sample,
#             timed             = is_on_timed_seg
#           ) %>%
#           dplyr::select(segment_id, gene, karyotype, sample_id, NV, DP, mutatation_status, mult_estimate, timed)
# 
#       } else if (is_on_timed_seg) {
#         amp    <- amp_by_gene[[g]]
#         parts  <- do.call(rbind, strsplit(amp$segment_name, "_"))
#         final_df <- amp %>%
#           mutate(
#             segment_id        = paste0(parts[,1],":",parts[,2],":",parts[,3],":",karyotype,":1"),
#             gene              = hgnc_symbol,
#             sample_id         = sample,
#             NV                = NA_real_,
#             DP                = NA_real_,
#             timed             = TRUE,
#             mutatation_status = mutation_status,
#             mult_estimate     = NA_real_
#           ) %>%
#           dplyr::select(segment_id, gene, karyotype, sample_id, NV, DP, mutatation_status, mult_estimate, timed)
# 
#       } else {
#         return(NULL)
#       }
# 
#       final_df %>% mutate(
#         ttype             = info_single$ttype,
#         cancer_type_short = info_single$cancer_type_short,
#         cancer_type       = info_single$cancer_type
#       )
#     })
# 
#     data.table::rbindlist(drivers_df, fill = TRUE) %>% as.data.frame()
# 
#   }, error = function(e) {
#     message(sprintf("ERROR in sample %d (%s): %s", idx, info_fit$sample[idx], e$message))
#     message(traceback())
#     return(NULL)
#   })
# }
# 
# # ── Run ───────────────────────────────────────────────────────────────────────
# # STEP 1: test on 3 samples sequentially first to catch any logic errors
# message("=== STEP 1: Testing on 3 samples sequentially ===")
# test_results <- lapply(1:3, process_sample)
# failed_tests <- which(sapply(test_results, is.null))
# if (length(failed_tests) > 0) {
#   stop(sprintf("Sequential test failed on sample indices: %s — fix errors before running in parallel.",
#                paste(failed_tests, collapse = ", ")))
# } else {
#   message("Sequential test passed. Launching parallel run...")
# }
# 
# # STEP 2: full parallel run
# message(sprintf("=== STEP 2: Full parallel run on %d samples ===", nrow(info_fit)))
# summary_gene_annotation_list <- future_lapply(
#   seq_len(nrow(info_fit)),
#   process_sample,
#   future.seed = TRUE
# )
# 
# # Check for failed samples
# failed_idx <- which(sapply(summary_gene_annotation_list, is.null))
# if (length(failed_idx) > 0) {
#   message(sprintf("WARNING: %d samples returned NULL (failed or empty): indices %s",
#                   length(failed_idx), paste(failed_idx, collapse = ", ")))
# }
# 
# summary_gene_annotation <- data.table::rbindlist(
#   Filter(Negate(is.null), summary_gene_annotation_list),
#   fill = TRUE
# ) %>% as.data.frame()
# 
# # ── Merge with timed segment info ─────────────────────────────────────────────
# message("Merging with segment timing info...")
# 
# gene_annotation        <- summary_gene_annotation %>% rename(mutation_status = mutatation_status)
# summary_segments_merge <- summary_segments %>% dplyr::select(-karyotype)
# 
# # Vectorised segment_name extraction (replaces slow rowwise strsplit)
# gene_annotation <- gene_annotation %>%
#   mutate(
#     segment_name = gsub(":", "_", sub("^([^:]+:[^:]+:[^:]+):.*", "\\1", segment_id)),
#     sample       = sample_id
#   ) %>%
#   left_join(summary_segments_merge, by = c("segment_name", "sample"))
# 
# saveRDS(gene_annotation, paste0(output_dir, "tmp_03_gene_annotation_table_parallel.rds"))
# message("Saved tmp_03_gene_annotation_table.rds")
# 
# ── Filter per cancer type ────────────────────────────────────────────────────
# message("Filtering by cancer type...")
# Intogen_drivers <- read.csv(
#   paste0(data_dir, "2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv"), sep = "\t"
# )
# gene_annotation <- readRDS(paste0(output_dir, "tmp_03_gene_annotation_table.rds"))
# 
# genes_by_tt <- split(Intogen_drivers$SYMBOL, Intogen_drivers$CANCER_TYPE)
# 
# annotation_new <- gene_annotation %>%
#   filter(!is.na(IntoGen_cancer_type)) %>%
#   rowwise() %>%
#   filter(
#     gene %in% (genes_by_tt[[IntoGen_cancer_type]] %||% character(0)) |
#       mutation_status != "WT"
#   ) %>%
#   ungroup()
# 
# saveRDS(annotation_new, paste0(output_dir, "03_gene_annotation_table_parallel.rds"))
# message("Done. Saved 03_gene_annotation_table_parallel.rds")
# 
