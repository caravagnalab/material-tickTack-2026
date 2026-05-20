base_dir <- "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/simulations/"
#!/usr/bin/env Rscript
.libPaths("~/R/orfeo_R_4.4/")

# Model-based clustering based on parameterized finite Gaussian mixture models. 
# Models are estimated by EM algorithm initialized by hierarchical model-based agglomerative clustering. 
# The optimal model is then selected according to BIC.

library(dplyr)
library(purrr)
library(tibble)
library(mclust)

# ── helpers ───────────────────────────────────────────────────────────────────

cor_safe <- function(x, y) {
  if (is.null(y) || all(is.na(y)) || sd(y, na.rm = TRUE) == 0 || sd(x, na.rm = TRUE) == 0) return(NA)
  cor(x, y, use = "complete.obs")
}

get_clusters <- function(tau) {
  as.numeric(factor(round(tau, 6)))
}

cluster_mclust <- function(tau_vals, K_max = 5) {
  tau_vals <- unname(tau_vals)
  valid    <- !is.na(tau_vals)
  
  if (sum(valid) < 2) {
    return(list(clusters  = rep(NA, length(tau_vals)),
                K         = NA,
                tau_clust = rep(NA, length(tau_vals))))
  }
  
  # all values identical — no variance to model, treat as single cluster
  if (sd(tau_vals[valid]) == 0) {
    return(list(clusters  = rep(1L, length(tau_vals)),
                K         = 1L,
                tau_clust = rep(tau_vals[valid][1], length(tau_vals))))
  }
  
  fit <- tryCatch(
    withCallingHandlers(
      Mclust(tau_vals[valid], G = 1:K_max, verbose = FALSE),
      warning = function(w) invokeRestart("muffleWarning")
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(list(clusters  = rep(NA, length(tau_vals)),
                K         = NA,
                tau_clust = rep(NA, length(tau_vals))))
  }
  
  cluster_labels        <- rep(NA_real_, length(tau_vals))
  cluster_labels[valid] <- fit$classification
  cluster_means         <- tapply(tau_vals[valid], fit$classification, mean)
  tau_clust             <- rep(NA_real_, length(tau_vals))
  tau_clust[valid]      <- cluster_means[fit$classification]
  
  list(clusters = cluster_labels, K = fit$G, tau_clust = tau_clust)
}

# ── main loop ─────────────────────────────────────────────────────────────────

parameteres_combination <- list.files(paste0(base_dir, "/results_simulations/"))
all_rows <- vector("list", length(parameteres_combination))

for (i in seq_along(parameteres_combination)) {
  
  sim_name   <- parameteres_combination[i]
  sim_path   <- paste0(base_dir, "/results_simulations/", sim_name)
  sim_params <- strsplit(sim_name, "_")[[1]]
  
  n_clusters       <- as.numeric(sim_params[2])
  n_events         <- sim_params[3]
  purity           <- sim_params[4]
  coverage         <- sim_params[5]
  mutation_density <- sim_params[6]
  true_K           <- n_clusters
  
  iter_folders <- list.files(sim_path)
  
  if (length(iter_folders) == 0) {
    all_rows[[i]] <- tibble(sim_name = sim_name, status = "NO_SUBFOLDERS")
    next
  }
  
  sim_rows <- vector("list", length(iter_folders))
  
  for (j in seq_along(iter_folders)) {
    
    iter.i         <- iter_folders[j]
    inference_file <- paste0(sim_path, "/", iter.i, "/merged_res_competing_methods.rds")
    
    na_row <- function(method) tibble(
      sim_name = sim_name, iter = iter.i, method = method,
      status = "RESULT_NOT_AVAILABLE",
      K = NA, correct_K = NA, MSE = NA, MSE_clust = NA,
      correlation = NA, adj_RI = NA,
      n_clusters = n_clusters, n_events = n_events,
      purity = purity, coverage = coverage,
      mutation_density = mutation_density
    )
    
    if (!file.exists(inference_file)) {
      sim_rows[[j]] <- bind_rows(lapply(
        c("Var-BIC", "ticktack_Single", "AmpTimeR", "MutTimeR"), na_row))
      next
    }
    
    res          <- readRDS(inference_file)
    true_tau     <- res$true_tau
    true_cluster <- res$true_tau_cluster
    
    # ── Var-BIC ───────────────────────────────────────────────
    tau_BIC       <- res[["tau_var_BIC"]]
    cl_BIC        <- get_clusters(tau_BIC)
    cl_means_BIC  <- tapply(tau_BIC, cl_BIC, mean, na.rm = TRUE)
    tau_clust_BIC <- cl_means_BIC[cl_BIC]
    
    # ── competing methods ─────────────────────────────────────
    competing <- list(
      list(label = "ticktack_Single", col = "tau_single"),
      list(label = "AmpTimeR",        col = "tau_AmpTimeR"),
      list(label = "MutTimeR",        col = "tau_MutTimeR")
    )
    
    rows_competing <- lapply(competing, function(m) {
      tau_vals <- res[[m$col]]
      cl       <- cluster_mclust(tau_vals, K_max = 5)
      tibble(
        sim_name = sim_name, iter = iter.i, method = m$label, status = "OK",
        K           = cl$K,
        correct_K   = if (is.na(cl$K)) NA else (cl$K == true_K),
        MSE         = mean((tau_vals      - true_tau)^2, na.rm = TRUE),
        MSE_clust   = mean((cl$tau_clust  - true_tau)^2, na.rm = TRUE),
        correlation = cor_safe(true_tau, tau_vals),
        adj_RI      = if (anyNA(cl$clusters)) NA else
          adjustedRandIndex(true_cluster, cl$clusters),
        n_clusters = n_clusters, n_events = n_events,
        purity = purity, coverage = coverage, mutation_density = mutation_density
      )
    })
    
    sim_rows[[j]] <- bind_rows(
      tibble(
        sim_name = sim_name, iter = iter.i, method = "Var-BIC", status = "OK",
        K           = length(unique(cl_BIC)),
        correct_K   = (length(unique(cl_BIC)) == true_K),
        MSE         = mean((tau_BIC       - true_tau)^2, na.rm = TRUE),
        MSE_clust   = mean((tau_clust_BIC - true_tau)^2, na.rm = TRUE),
        correlation = cor_safe(true_tau, tau_BIC),
        adj_RI      = adjustedRandIndex(true_cluster, cl_BIC),
        n_clusters = n_clusters, n_events = n_events,
        purity = purity, coverage = coverage, mutation_density = mutation_density
      ),
      rows_competing
    )
    
    rm(res, true_tau, true_cluster, tau_BIC, cl_BIC, tau_clust_BIC,
       rows_competing)
    gc(verbose = FALSE)
  }
  
  all_rows[[i]] <- bind_rows(sim_rows)
  rm(sim_rows)
  gc(verbose = FALSE)
  
  message("Done: ", sim_name, " (", i, "/", length(parameteres_combination), ")")
}

results_all <- bind_rows(all_rows)

saveRDS(results_all,
        paste0(base_dir, "results_summary/07_results_competing_methods_BIC.rds"))

message("Saved: results_summary/07_results_competing_methods_BIC.rds")
