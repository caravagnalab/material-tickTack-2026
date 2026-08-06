# Process data to get input for signatures
library(dplyr)
library(RESOLVE)

### Input data
# pcawg_fits_dir = "~/dati_Orfeo/scocomello/material-tickTack-2026/PCAWG/Fit/inference_results_5ncomponents"
pcawg_fits_dir = "../Fit/inference_results_5ncomponents"

pcawg_fits = list.files(pcawg_fits_dir)

# samples = readRDS("~/dati_Orfeo/antonelloa/material-tickTack-2026/PCAWG/Data/Samples.rds")
# segments = readRDS("~/dati_Orfeo/antonelloa/material-tickTack-2026/PCAWG/Data/Segments.rds")
samples = readRDS("../Data/Samples.rds")
segments = readRDS("../Data/Segments.rds")


ttypes = samples %>% filter(!is.na(IntoGen_cancer_type)) %>% pull(IntoGen_cancer_type) %>% unique()
ttypes <- "OVT"

save_matrix = function(nmf_matrix,path,file_name){
  
  nmf_matrix <- t(nmf_matrix)
  
  in_df <- nmf_matrix %>% as.data.frame() %>% tibble::rownames_to_column("Mutation Types")
  
  data.table::fwrite(data.table::as.data.table(in_df), 
                     file = file.path(path,file_name), 
                     sep = "\t")
  
}

### Create input for signature deconvolution
for (t in ttypes){
  print(t)
  samples_tt = samples %>% filter(!is.na(IntoGen_cancer_type)) %>% filter(IntoGen_cancer_type == t) %>% pull(sample)
  if (length(samples_tt) > 2){
  ### Process data to call signatures on all clusters together
  cna_data = lapply(samples_tt, function(s){
    x = readRDS(paste0(pcawg_fits_dir, "/",s, ".rds"))
    cna_data = x$cna %>% select("chr","from","to","Major","minor")
    cna_data$sample = s
    cna_data
  }) %>% Reduce(rbind, .)
  chrom = strsplit(cna_data$chr, "chr") %>% unlist()
  chrom = chrom[chrom != ""]
  cna_data_new = 
    data.frame(
      "sample" = cna_data$sample,
      "chrom" = chrom,
      "start" = cna_data$from,
      "end" = cna_data$to,
      "major" = cna_data$Major,
      "minor" = cna_data$minor
    )
  cna_counts = RESOLVE::getCNCounts(cna_data_new)
  
  
  ### Process data to call signatures separately for each cluster (clusters with enough CNA events)
  # cna_data_cluster = lapply(samples_tt, function(s){
  #   x = readRDS(paste0(pcawg_fits_dir, "/",s, ".rds"))
  #   cna_data = x$cna %>% select("chr","from","to","Major","minor")
  #   segments %>% filter(sample == s) %>% select(segment_name, clock_mean)
  #   cna_data$sample = s
  #   cna_data
  # }) %>% Reduce(rbind, .)
  # cna_counts_by_cluster = RESOLVE::getCNCounts(cna_data_cluster)
  
  
  tt_dir = paste0("../Analysis_results/CNA_Signatures/raw_data/", t)
  dir.create(tt_dir)
  saveRDS(cna_counts, file = paste0(tt_dir, "/CNV_counts.rds"))
  save_matrix(cna_counts,tt_dir,"CNV_counts.txt")
  }
  
}

## Run signature deconvolution
library(reticulate)
# use_condaenv(
#   "/Users/aliceantonello/Library/r-miniconda/envs/sigprofiler_arm64",
#   required = TRUE
# )


use_condaenv( 
  "/u/cdslab/scocomello/scratch/miniconda/envs/sigprofiler", 
  required = TRUE
)

# py_install("SigProfilerExtractor", envname = "sigprofiler_env", pip = TRUE)
# py_install("SigProfilerExtractor", envname = "sigprofiler", pip = TRUE)

sig_ex <- import("SigProfilerExtractor.sigpro")

create_sigprofiler_dirs <- function(output, context = "CNV48", max_sigs = 10) {
  subdirs <- c("Signatures", "Activities", "Plots", "Solution_Stats",
               "Decomposed_Solution", "Suggested_Solution")
  
  dirs_to_create <- c()
  
  for (n_sig in 1:max_sigs) {
    base <- file.path(output, context, "All_Solutions",
                      paste0(context, "_", n_sig, "_Signatures"))
    for (d in subdirs) dirs_to_create <- c(dirs_to_create, file.path(base, d))
  }
  
  suggested_base <- file.path(output, context, "Suggested_Solution",
                              paste0(context, "_De-Novo_Solution"))
  for (d in subdirs) dirs_to_create <- c(dirs_to_create, file.path(suggested_base, d))
  
  dirs_to_create <- c(dirs_to_create,
                      file.path(output, context, "Plots"),
                      file.path(output, context, "Seeds"),
                      file.path(output, context, "Matrices")
  )
  
  for (d in dirs_to_create) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    Sys.chmod(d, mode = "0755")  # ensure write permission
  }
  
  # Also chmod the full output tree
  system(paste("chmod -R 755", shQuote(output)))
}

counts_folder = "../Analysis_results/CNA_Signatures/raw_data"

for (t in ttypes){
  samples_tt = samples %>% filter(!is.na(IntoGen_cancer_type)) %>% filter(IntoGen_cancer_type == t) %>% pull(sample)
  if (length(samples_tt) > 2){
    x <- paste0(counts_folder,"/",t,"/CNV_counts.txt")
    
    output = path.expand(paste0("../Analysis_results/CNA_Signatures/Sigprofiler/",t))
    create_sigprofiler_dirs(output, context = "CNV48", max_sigs = 10)
    
    sig_ex$sigProfilerExtractor(
      input_type   = "matrix",
      output      = output,   # output folder
      input_data   = x,       # path to your matrix
      reference_genome = "GRCh37",           # or "GRCh37"
      opportunity_genome = "GRCh37",
      context_type = "CNV48",                # CNV48 for copy number
      minimum_signatures = 1L,
      maximum_signatures = 10L,              # adjust to your sample size
      nmf_replicates = 100L,                 # higher = more stable, slower
      cpu = -1L                              # -1 uses all available cores
    )
    
    # # Best solution signatures
    # sigs <- read.table(
    #   "SigProfiler_output/CNV48/Suggested_Solution/CNV48_De-Novo_Solution/Signatures/CNV48_De-Novo_Signatures.txt",
    #   header = TRUE, sep = "\t", row.names = 1
    # )
    # 
    # # Sample exposures (how much each signature contributes per sample)
    # exposures <- read.table(
    #   "SigProfiler_output/CNV48/Suggested_Solution/CNV48_De-Novo_Solution/Activities/CNV48_De-Novo_Activities_refit.txt",
    #   header = TRUE, sep = "\t", row.names = 1
    # )
    # 
    # head(sigs)
    # head(exposures)
  }
}
















