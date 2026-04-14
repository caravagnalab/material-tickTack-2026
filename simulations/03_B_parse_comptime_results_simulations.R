base_dir = "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/simulations/"
#!/usr/bin/env Rscript
.libPaths("~/R/orfeo_R_4.4/")
library(dplyr)
library(purrr)
library(tibble)
library(tidyr)

methods <- c("variational", "mcmc")

parameteres_combination <- list.files(paste0(base_dir, "/results_simulations/"))

results_all_K <- map_dfr(parameteres_combination, function(sim_name) {
  
  sim_path <- paste0(base_dir, "/results_simulations/", sim_name)
  
  sim_params <- strsplit(sim_name, "_")[[1]]
  n_clusters <- as.numeric(sim_params[2])
  n_events <- sim_params[3]
  purity <- sim_params[4]
  coverage <- sim_params[5]
  mutation_density <- sim_params[6]
  
  iter_folders <- list.files(sim_path)
  if (length(iter_folders) == 0) {
    return(tibble(sim_name = sim_name, status = "NO_SUBFOLDERS"))
  }
  
  map_dfr(iter_folders, function(iter.i) {
    
    iter_path <- paste0(sim_path, "/", iter.i)
    timing_file <- paste0(iter_path, "/timing_results.rds")
    inference_file <- paste0(iter_path, "/inference_results.rds")
    
    # Read timing results if available
    time_info <- if (file.exists(timing_file)) readRDS(timing_file) else NULL
    
    # Read inference results if available
    inference_results <- if (file.exists(inference_file)) readRDS(inference_file) else NULL
    
    map_dfr(methods, function(m) {
      
      # Extract failed Ks for this method
      failed_K <- if (!is.null(inference_results)) {
        if (m == "variational") inference_results$failed_K_variational else inference_results$failed_K_mcmc
      } else NA
      
      # Get the timing info per K
      if (!is.null(time_info) && any(time_info$method == m)) {
        time_per_K <- time_info %>% filter(method == m) %>%
          dplyr::select(K, time_seconds)
      } else {
        time_per_K <- tibble(K = NA, time_seconds = NA)
      }
      
      # Build a tibble per K
      time_per_K %>%
        mutate(
          sim_name = sim_name,
          iter = iter.i,
          inference = m,
          status = if (!is.null(inference_results)) "OK" else "RESULT_NOT_AVAILABLE",
          failed = ifelse(K %in% failed_K, TRUE, FALSE),
          n_clusters = n_clusters,
          n_events = n_events,
          purity = purity,
          coverage = coverage,
          mutation_density = mutation_density
        )
    })
  })
})

# Save the detailed table
saveRDS(results_all_K, paste0(base_dir, "results_summary/03_B_results_comptime_simulations_per_K.rds"))
