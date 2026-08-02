setwd("../material-tickTack-2026/PCAWG/Scripts/")
Segments = readRDS("../Data/Segments.rds")

Segments %>% select(sample, best_K) %>% distinct() %>% select(best_K) %>% table()

fit_list <- list.files("/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Fit/inference_results_5ncomponents")

for (i in 1:50){
  
  row_i = Segments %>% filter(sample == tools::file_path_sans_ext(basename(fit_list[i])))
  fit <- readRDS(paste0("/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Fit/inference_results_5ncomponents/",fit_list[i]))
  
  # fit$results_timing$data$input_data$seg_assignment
  mut_number <- as.data.frame(table(fit$results_timing$data$input_data$seg_assignment))$Freq
  print(mean(mut_number))
  
  row_i$best_K
  n_cna_simple = length(mut_number)
  
}




# use this to perform clustering 
plot_vaf_ecdf = function(x, K) {
  input_data = x$results_timing$data$input_data
  cna_data = x$results_timing$draws_and_summary[[K]]$summarized_results
  
  cna_data_with_vaf = lapply(cna_data$segment_id, function(i) {
    idxs = which(input_data$seg_assignment == i)
    DP = input_data$DP[idxs]
    NV = input_data$NV[idxs]
    VAF = NV / DP
    
    cna_data %>% dplyr::filter(segment_id == i) %>%
      dplyr::bind_cols(dplyr::tibble(VAF = VAF))
  }) %>% do.call(dplyr::bind_rows, .)
  
  cna_data_with_vaf %>%
    dplyr::mutate(clock_rank = dplyr::dense_rank(clock_mean)) %>%
    ggplot2::ggplot(aes(x = VAF, color = factor(clock_rank), group = segment_name)) +
    ggplot2::stat_ecdf() +
    ggplot2::facet_grid(clock_rank~karyotype) +
    ggplot2::theme_light() +
    ggplot2::scale_x_continuous(limits = c(-0.01,1.01)) +
    labs(fill = "Cluster", col = "Cluster")
}
