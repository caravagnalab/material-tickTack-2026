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

###################################### get coordinates of genes of interest ########################################

# create gene list for gene coords

CI_drivers <- readRDS(paste0(data_dir, "CI_drivers.rds"))
Intogen_drivers <- read.csv(paste0(data_dir, "2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv"), sep = "\t") 
list_of_drivers <- c(CI_drivers, (Intogen_drivers$SYMBOL %>% unique()) )


# to match the hg19 for PCAWG use the GRCh37 mart
ensembl <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  GRCh = 37
)

gene_coords <- getBM(
  attributes = c("ensembl_gene_id",
                 "hgnc_symbol",
                 "chromosome_name",
                 "start_position",
                 "end_position"),
  filters = "hgnc_symbol",
  values = list_of_drivers,
  mart = ensembl
)

saveRDS(gene_coords,paste0(data_dir, "gene_coords.rds"))

missing_annotation <- setdiff(list_of_drivers, gene_coords$hgnc_symbol)
missing_annotation_Intogene <- setdiff(Intogen_drivers$SYMBOL, gene_coords$hgnc_symbol)
missing_annotation_CI <- setdiff(CI_drivers, gene_coords$hgnc_symbol)

################## add info for missing annotations genes #####################

gene_coords <- readRDS(paste0(data_dir, "gene_coords.rds"))
missing_annotation_gene_names <- read.csv(paste0(data_dir, "missing_driver_annotation_info.tsv"), sep = "\t")

missing_annotation_ensembl <- missing_annotation_gene_names %>% pull(Ensembl.gene.id)

ensembl <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  GRCh = 37
)

saveRDS(ensembl,paste0(data_dir, "ensembl.rds"))


gene_coords_missing_drivers <- getBM(
  attributes = c("ensembl_gene_id",
                 "hgnc_symbol",
                 "chromosome_name",
                 "start_position",
                 "end_position"),
  filters = "ensembl_gene_id",   
  values = missing_annotation_ensembl,
  mart = ensembl
)

new_gene_symbol_for_missing_genes <- gene_coords_missing_drivers

gene_coords_all <- bind_rows(gene_coords , gene_coords_missing_drivers %>% mutate(chromosome_name = as.character(chromosome_name)))
saveRDS(gene_coords_all,paste0(data_dir, "gene_coords_IntoGene_CI.rds"))

new_names <- gene_coords_missing_drivers %>% left_join(missing_annotation_gene_names, by= join_by(ensembl_gene_id == Ensembl.gene.id))
new_drivers_symbols <- list(IntoGene_drivers = c( setdiff(gene_coords$hgnc_symbol, CI_drivers) , new_names %>% filter(Original.name %in% missing_annotation_Intogene) %>% pull(hgnc_symbol)),
                            CI_genes = c(setdiff(gene_coords$hgnc_symbol, Intogen_drivers$SYMBOL), new_names %>% filter(Original.name %in% missing_annotation_CI) %>% pull(hgnc_symbol)))

saveRDS(new_drivers_symbols,paste0(data_dir, "driver_list_IntoGene_CI.rds"))

