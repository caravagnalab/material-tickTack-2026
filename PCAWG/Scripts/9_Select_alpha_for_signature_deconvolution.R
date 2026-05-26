suppressPackageStartupMessages({
  library(tidyverse)
  # library(ggrepel)
  library(dplyr)
  library(posterior)
})

# ## Parse command line arguments
# args <- commandArgs(trailingOnly = TRUE)
# usage <- "
# Usage: Rscript 9 Select alpha for signature deconvolution.R \\
#   --sample    <sample>
# "
# # Helper: extract value after a flag
# get_arg <- function(args, flag, default = NULL) {
#   idx <- which(args == flag)
#   if (length(idx) == 0) return(default)
#   if (idx == length(args)) stop(paste("Flag", flag, "requires a value."))
#   args[idx + 1]
# }
# if ("--help" %in% args || "-h" %in% args) {
#   cat(usage)
#   quit(status = 0)
# }


## Import data
fits_dir = "/orfeo/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Fit/inference_results/"
Samples = list.files(fits_dir)
# input_path    <- get_arg(args, "--input", default = 1)
all_samples_df = lapply(Samples, function(s){
  # x = readRDS(paste0(fits_dir, Samples[input_path]))
  # fce8d8c6-f2a0-43a8-9a7a-b9c519a6686c
  # x = readRDS(paste0(fits_dir, "fce8d8c6-f2a0-43a8-9a7a-b9c519a6686c.rds"))
  x = readRDS(paste0(fits_dir, s))
  
  mutations = x$mutations 
  timed_segments = x[["results_timing"]][["data"]][["accepted_cna"]] %>% 
    separate(segment_name, into = c("chr", "from", "to"), sep = "_") %>%
    mutate(segment_id = paste0(chr,":",from,":",to,":",karyotype,":1")) %>%
    pull(segment_id)
  mutations = mutations %>% filter(segment_id %in% timed_segments)
  purity = x$purity
  
  ## Filter non tail mutations 
  get_clonal_peaks = function(k, purity) {
    multiplicities <- strsplit(k, ":") %>% unlist() %>% as.numeric()
    major <- multiplicities[1]
    n_tot <- sum(multiplicities)
    
    # get only Major and 1
    multiplicities <- c(1, major)
    
    peaks <- unique(multiplicities) * purity / (n_tot * purity + 2 * (1 - purity))
    return(sort(peaks))
  }
  # segment_id = "chr1:170045273:174114731:2:1:1" 
  filter_mutations = function(mutations, segment_id, purity, alpha = 0.05){
    
    k = paste0(strsplit(segment_id, ":")[[1]][4], ":",strsplit(segment_id, ":")[[1]][5])
    seg_chr = strsplit(segment_id, ":")[[1]][1]
    seg_from = strsplit(segment_id, ":")[[1]][2] %>% as.double()
    seg_to = strsplit(segment_id, ":")[[1]][3]  %>% as.double()
    
    peaks <- get_clonal_peaks(k, purity)
    
    segment_mutations <- mutations %>%
        dplyr::filter(chr == seg_chr, from > seg_from, to < seg_to) %>%
        tidyr::drop_na(DP)
      
      accepted_mutations <- data.frame()
      if (nrow(segment_mutations) > 0) {
        # Check if mutation is inside CI
        probs <- c(alpha/2, 1 - alpha/2)
        
        DP <- segment_mutations$DP
        NV <- segment_mutations$NV
        
        accepted_idx <- lapply(1:length(DP), function(i) {
          #fisso i picchi per il segmento e vedo se tutte le mutazioni ricadono in almeno uno dei due intervalli intorno ai picchiche sono diversi a seconda del valore di DP per la specifica mutazione
          for (p in peaks) {
            quantiles <- stats::qbinom(probs, DP[i], p)
            if ((NV[i] >= quantiles[1]) && (NV[i] <= quantiles[2])) {
              return(i)
            }
          }
        }) %>% unlist()
        
        # Get only good mutations
        accepted_mutations <- data.frame(DP = DP[accepted_idx], NV = NV[accepted_idx])
      }
      accepted_mutations = paste0(accepted_mutations$NV, accepted_mutations$DP)
      segment_mutations %>% mutate(nv_dp = paste0(NV,DP)) %>% filter(nv_dp %in% accepted_mutations)
  }
  filtered_muts = lapply(mutations$segment_id %>% unique(), function(s){
    filter_mutations(mutations, s, purity, alpha = 0.05)
  }) %>% Reduce(rbind, .)
  
  best_K = x$results_timing$best_K
  segment_clock_table <- x$results_timing$draws_and_summary[[best_K]]$summarized_results #%>%
    #dplyr::select(segment_name, karyotype, chr, clock_mean, clock_low, clock_high)
  
  # x$results_timing$draws_and_summary[[best_K]]$draws 
  segment_posterior = lapply(1:nrow(segment_clock_table), function(s){
    segment = segment_clock_table$segment_name[s]
    index = segment_clock_table$segment_id[s]
    segment_tau_median = segment_clock_table %>% filter(segment_name == segment) %>% pull(clock_median)
    segment_tau_index = strsplit(strsplit(names(segment_tau_median), "\\[")[[1]][2], "\\]")[[1]][1] %>% as.double()
    
    draws_df = as_draws_df(x$results_timing$draws_and_summary[[best_K]]$draws) #("tau", format = "matrix")    
    
    theta_single = draws_df[[paste0("theta[", index, ",", segment_tau_index, ",1]")]] %>% median()
    theta_double = draws_df[[paste0("theta[", index, ",", segment_tau_index, ",2]")]] %>% median()
    
    data.frame(
      "segment_id"= segment, 
      "theta_1"= theta_single,
      "theta_2"= theta_double
    )
  }) %>% Reduce(rbind,.)
  
  ### Generate final table
  # Sample id, Segment id, Karyo, Intogen_ttype, ref, alt, from, to, NV, DP, timing
  segments = filtered_muts$segment_id %>% unique()
  alphas_df = lapply(segments, function(s){
    segment_id_tick_tack = paste0(strsplit(s, ":")[[1]][1], "_", strsplit(s, ":")[[1]][2], "_",strsplit(s, ":")[[1]][3])
    k = paste0(strsplit(s, ":")[[1]][4], ":", strsplit(s, ":")[[1]][5])
    muts_on_s = filtered_muts %>% 
      filter(mutation_type == "single base substitution") %>%
      select(sample, segment_id, karyotype, ref, alt, chr, from, to, NV, DP) %>%
      filter(segment_id == s)
    muts_on_s$time = segment_clock_table %>% filter(segment_name == segment_id_tick_tack) %>% pull(clock_mean)
    theta_1 = segment_posterior %>% filter(segment_id == segment_id_tick_tack) %>% pull(theta_1)
    theta_2 = segment_posterior %>% filter(segment_id == segment_id_tick_tack) %>% pull(theta_2)
    peaks = get_clonal_peaks(k, purity) %>% sort()
    low_peak = peaks[1]
    high_peak = peaks[2]
    # Assign the mutation to the alpha or beta peak according to the likelihood
    muts_on_s %>% rowwise() %>%
      mutate(
        ll_single_copy = theta_1*dbinom(x = NV, size = DP, prob = low_peak),
        ll_double_copy = theta_2*dbinom(x = NV, size = DP, prob = high_peak)
      ) %>%
      mutate(is_alpha = ifelse(ll_double_copy > ll_single_copy, T,F)) %>%
      filter(is_alpha) %>% select(!c("ll_single_copy", "ll_double_copy", "is_alpha"))
  }) %>% Reduce(rbind, .)
  
  alphas_df
  
  # write.table(alphas_df, 
  #           file = paste0("/orfeo/scratch/cdslab/antonelloa/material-tickTack-2026/PCAWG/Fit/SBS_data/", 
  #                         muts_on_s$sample %>% unique(), ".txt"),
  #           row.names = F, quote = F, sep = "\t")
}) %>% Reduce(rbind, .)


write.table(all_samples_df, 
            file = paste0("/orfeo/scratch/cdslab/antonelloa/material-tickTack-2026/PCAWG/Fit/SBS_data/alpha_muts.txt"),
            row.names = F, quote = F, sep = "\t")






