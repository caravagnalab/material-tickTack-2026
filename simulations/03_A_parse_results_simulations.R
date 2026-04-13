base_dir = "/orfeo/cephfs/scratch/cdslab/scocomello/material_tickTack_2026/simulations/"
#!/usr/bin/env Rscript
.libPaths("~/R/orfeo_R_4.4/")

library(dplyr)
library(purrr)
library(tibble)
library(mclust)
library(tidyr)

get_clusters <- function(tau) {
  as.numeric(factor(round(tau, 6)))
}

cor_safe <- function(x, y) {
  if (is.null(y) || all(is.na(y)) || sd(y) == 0 || sd(x) == 0) return(NA)
  cor(x, y)
}

methods <- c(
  "var_BIC", "var_AIC", "var_ICL", "var_LOO"
  # , "mcmc_BIC", "mcmc_AIC", "mcmc_ICL", "mcmc_LOO"
)

parameteres_combination <- list.files(paste0(base_dir, "/results_simulations/"))

results_all <- map_dfr(parameteres_combination, function(sim_name) {
  
  sim_path <- paste0(base_dir, "/results_simulations/", sim_name)
  
  sim_params <- strsplit(sim_name, "_")[[1]]
  
  n_clusters <- as.numeric(sim_params[2])
  n_events <- sim_params[3]
  purity <- sim_params[4]
  coverage <- sim_params[5]
  mutation_density <- sim_params[6]
  
  iter_folders <- list.files(sim_path)
  
  if (length(iter_folders) == 0) {
    return(tibble(
      sim_name = sim_name,
      status = "NO_SUBFOLDERS"
    ))
  }
  
  map_dfr(iter_folders, function(iter.i) {
    
    iter_path <- paste0(sim_path, "/", iter.i)
    
    sim_file <- paste0(iter_path, "/sim.rds")
    timing_file <- paste0(iter_path, "/timing_results.rds")
    inference_file <- paste0(iter_path, "/inference_results.rds")
    
    if (!file.exists(inference_file)) {
      return(tibble(
        sim_name = sim_name,
        iter = iter.i,
        inference = NA,
        criterion = NA,
        # inference = strsplit(m, "_")[[1]][1],
        # criterion = strsplit(m, "_")[[1]][2],
        status = "RESULT_NOT_AVAILABLE",
        K = NA,
        correct_K = NA,
        correlation = NA,
        adj_RI = NA,
        failed_K = NA,
        n_clusters = n_clusters,
        n_events = n_events,
        purity = purity,
        coverage = coverage,
        mutation_density = mutation_density
      ))
    }
    
    inference_results <- readRDS(inference_file)
    
    inference_results <- inference_results %>%
      mutate(
        cluster_var_BIC = get_clusters(`tau_var_BIC <- ...`),
        cluster_var_AIC = get_clusters(`tau_var_AIC <- ...`),
        cluster_var_ICL = get_clusters(`tau_var_ICL <- ...`),
        cluster_var_LOO = get_clusters(`tau_var_LOO <- ...`)
        # ,
        
        # cluster_mcmc_BIC = get_clusters(`tau_mcmc_BIC <- ...`),
        # cluster_mcmc_AIC = get_clusters(`tau_mcmc_AIC <- ...`),
        # cluster_mcmc_ICL = get_clusters(`tau_mcmc_ICL <- ...`),
        # cluster_mcmc_LOO = get_clusters(`tau_mcmc_LOO <- ...`)
      )
    tau_cols <- paste0("tau_", methods, " <- ...")
    true_K <- n_clusters
    
    map_dfr(methods, function(m) {
      
      tau_col <- paste0("tau_", m, " <- ...")
      cluster_col <- paste0("cluster_", m)

      tau_values <- inference_results[[tau_col]]
      cluster_values <- inference_results[[cluster_col]]
      
      K_est <- length(unique(cluster_values))
      
      failed_K_value <- if (grepl("^var", m)) {
        paste(unique(inference_results$failed_K_variational), collapse = ",")
      } else {
        paste(unique(inference_results$failed_K_mcmc), collapse = ",")
      }
      
      tibble(
        sim_name = sim_name,
        iter = iter.i,
        inference = strsplit(m, "_")[[1]][1],
        criterion = strsplit(m, "_")[[1]][2],
        status = "OK",
        K = K_est,
        correct_K = (K_est == true_K),
        correlation = cor_safe(inference_results$true_tau, unname(tau_values)),
        adj_RI = adjustedRandIndex(
          inference_results$true_tau_cluster,
          cluster_values
        ),
        failed_K = failed_K_value,
        n_clusters = n_clusters,
        n_events = n_events,
        purity = purity,
        coverage = coverage,
        mutation_density = mutation_density
      )
    })
  })
})


saveRDS(results_all, paste0(base_dir,"results_summary/03_A_results_simulations.rds"))




# ###################### if inference results is not available only take the variational one if it exists #######################
# test_mcmc = FALSE
# if (test_mcmc == TRUE){
#   library(dplyr)
#   library(purrr)
#   library(tibble)
#   library(mclust)
#   library(tidyr)
#   
#   build_from_ms_var <- function(ms_var, true_tau) {
#     
#     df <- tibble(
#       segment_idx = ms_var$seg_assignment_AIC$segment_idx,
#       true_tau = true_tau
#     )
#     
#     criteria <- c("AIC", "BIC", "ICL", "LOO")
#     
#     for (crit in criteria) {
#       tau_vals <- ms_var[[paste0("seg_assignment_", crit)]]$tau_assignment
#       df[[paste0("tau_var_", crit, " <- ...")]] <- tau_vals
#     }
#     
#     df$failed_K_variational <- paste(ms_var$failed_K, collapse = ",")
#     
#     df
#   }
#   
#   get_clusters <- function(tau) {
#     as.numeric(factor(round(tau, 6)))
#   }
#   
#   cor_safe <- function(x, y) {
#     if (is.null(y) || all(is.na(y)) || sd(y) == 0 || sd(x) == 0) return(NA)
#     cor(x, y)
#   }
#   
#   methods <- c(
#     "var_BIC", "var_AIC", "var_ICL", "var_LOO",
#     "mcmc_BIC", "mcmc_AIC", "mcmc_ICL", "mcmc_LOO"
#   )
#   
#   parameteres_combination <- list.files(paste0(base_dir, "/results_simulations/"))
#   
#   results_all_var <- map_dfr(parameteres_combination, function(sim_name) {
#     
#     sim_path <- paste0(base_dir, "/results_simulations/", sim_name)
#     
#     sim_params <- strsplit(sim_name, "_")[[1]]
#     
#     n_clusters <- as.numeric(sim_params[2])
#     n_events <- sim_params[3]
#     purity <- sim_params[4]
#     coverage <- sim_params[5]
#     mutation_density <- sim_params[6]
#     
#     iter_folders <- list.files(sim_path)
#     
#     if (length(iter_folders) == 0) {
#       return(tibble(
#         sim_name = sim_name,
#         status = "NO_SUBFOLDERS"
#       ))
#     }
#     
#     map_dfr(iter_folders, function(iter.i) {
#       
#       iter_path <- paste0(sim_path, "/", iter.i)
#       
#       sim_file <- paste0(iter_path, "/sim.rds")
#       inference_file <- paste0(iter_path, "/inference_results.rds")
#       ms_var_path <- paste0(iter_path, "/ms_results_variational.rds")
#       
#       inference_results <- NULL
#       status_global <- "OK"
#       
#       if (file.exists(inference_file)) {
#         
#         inference_results <- readRDS(inference_file)
#         
#       } else if (file.exists(ms_var_path)) {
#         
#         ms_var <- readRDS(ms_var_path)
#         sim <- readRDS(sim_file)
#         
#         inference_results <- build_from_ms_var(ms_var, sim$true_tau)
#         
#         # add missing columns for compatibility
#         inference_results$true_tau_cluster <- sim$taus_clust
#         inference_results$failed_K_mcmc <- NA
#         
#         status_global <- "ONLY_VARIATIONAL"
#         
#       } else {
#         
#         return(tibble(
#           sim_name = sim_name,
#           iter = iter.i,
#           status = "NO_RESULTS_AVAILABLE"
#         ))
#       }
#       
#       for (m in methods) {
#         tau_col <- paste0("tau_", m, " <- ...")
#         cluster_col <- paste0("cluster_", m)
#         
#         if (tau_col %in% colnames(inference_results)) {
#           inference_results[[cluster_col]] <- get_clusters(inference_results[[tau_col]])
#         }
#       }
#       
#       true_K <- n_clusters
#       
#       map_dfr(methods, function(m) {
#         
#         tau_col <- paste0("tau_", m, " <- ...")
#         cluster_col <- paste0("cluster_", m)
#         
#         tau_exists <- tau_col %in% colnames(inference_results)
#         cluster_exists <- cluster_col %in% colnames(inference_results)
#         
#         if (!tau_exists || !cluster_exists) {
#           return(tibble(
#             sim_name = sim_name,
#             iter = iter.i,
#             inference = strsplit(m, "_")[[1]][1],
#             criterion = strsplit(m, "_")[[1]][2],
#             status = paste0("NOT_AVAILABLE_", status_global),
#             K = NA,
#             correct_K = NA,
#             correlation = NA,
#             adj_RI = NA,
#             failed_K = NA,
#             n_clusters = n_clusters,
#             n_events = n_events,
#             purity = purity,
#             coverage = coverage,
#             mutation_density = mutation_density
#           ))
#         }
#         
#         tau_values <- inference_results[[tau_col]]
#         cluster_values <- inference_results[[cluster_col]]
#         
#         K_est <- length(unique(cluster_values))
#         
#         failed_K_value <- if (grepl("^var", m)) {
#           inference_results$failed_K_variational[1]
#         } else {
#           inference_results$failed_K_mcmc[1]
#         }
#         
#         tibble(
#           sim_name = sim_name,
#           iter = iter.i,
#           inference = strsplit(m, "_")[[1]][1],
#           criterion = strsplit(m, "_")[[1]][2],
#           status = status_global,
#           K = K_est,
#           correct_K = (K_est == true_K),
#           correlation = cor_safe(inference_results$true_tau, unname(tau_values)),
#           adj_RI = if (!all(is.na(inference_results$true_tau_cluster))) {
#             adjustedRandIndex(
#               inference_results$true_tau_cluster,
#               cluster_values
#             )
#           } else NA,
#           failed_K = failed_K_value,
#           n_clusters = n_clusters,
#           n_events = n_events,
#           purity = purity,
#           coverage = coverage,
#           mutation_density = mutation_density
#         )
#       })
#     })
#   })
#   
#   saveRDS(results_all_var, paste0(base_dir,"results_summary/results_all_simulations_variational.rds"))
#   
# }
