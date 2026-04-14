# rm(list=ls())
# n_clocks = 2
# n_events = 10
# purity = 0.8
# coverage = 80
# mutation_density = 2.e-07

base_dir = "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/simulations/"
#!/usr/bin/env Rscript
.libPaths("~/R/orfeo_R_4.4/")
library(tickTack)
library(patchwork)
library(cmdstanr)
library(posterior)
library(dplyr)
library(tidyr)
library(ggplot2)
args <- commandArgs(TRUE)
n_clocks <- as.double(args[1])
n_events <- as.double(args[2])
pi = purity <- as.double(args[3])
coverage <- as.double(args[4])
mutation_density <- as.numeric(args[5])

source(paste0(base_dir,"/utils.R"))

MIN_MUTATIONS = 10

epsilon <- 0.16
tolerance <- 0.0001
max_attempts <- 2
seed <- 1234
set.seed(seed)

print (paste0("n_clocks: ",n_clocks)) 
print (paste0("n_events: ",n_events))
print(paste0("purity: ",purity)) 
print(paste0("coverage: ",coverage))
print(paste0("mutation_density: ",mutation_density)) 
print(paste0("epsilon: ",epsilon))
print(paste0("tolerance: ",tolerance)) 
print(paste0("max_attempts: ",max_attempts)) 
print(paste0("seed: ",seed))

##################################################################################################
main_dir = paste0(base_dir,"/results_simulations/sim_",
                  n_clocks, "_", n_events, "_", purity, "_", coverage, "_", mutation_density, "/")

if (!dir.exists(main_dir)) {
  dir.create(main_dir)  
}

for (i.iter in 1:10) {
  sub_dir = paste0(main_dir, i.iter)
  if (!dir.exists(sub_dir)) {
    dir.create(sub_dir)
    
    sim = simulate_dataset(n_events, n_clocks, mutation_density, purity, coverage, sigma_tau = .01, min_dist = epsilon)
    saveRDS(sim, paste0(sub_dir, "/sim.rds"))
    
    p1 <- plot_sim(sim)
    ggplot2::ggsave(paste0(sub_dir,"/sim.png"),plot = p1, width = 10, height = 8)

    data <- prepare_tickTack_input_data(sim, purity, MIN_MUTATIONS, alpha = 0) 
    input_data = data$input_data
    saveRDS(data, paste0(sub_dir, "/data.rds"))
  
    
    m5 = cmdstanr::cmdstan_model(paste0(base_dir,"model/v5.stan"))

  ##################################### Inference variational H ###################################################
    # with n_clocks = 3 , n_events = 10, purity = 0.8, coverage = 60, mutation_density = 1e-6 the variational fit does not converge for K 0 5 
    # see how many times does this happens and in which conditions 
    
    library(loo)
    
    fits_and_time_var <- model_selection_fit(K = 5, input_data, m5, inference_type = "variational")
    fits_var <- fits_and_time_var$fits
    timing_var <- fits_and_time_var$timing
    
    ms_results_var = model_selection(fits_var, input_data, K = 5)
    best_k_var_BIC = ms_results_var$best_K_BIC
    best_fit_var_BIC = fits_var[best_k_var_BIC]
    seg_ass_var_BIC <- ms_results_var$seg_assignment_BIC
    
    p_fit <- ggplot(data = data.frame(x = seg_ass_var_BIC$tau_assignment, y = sim$true_taus), 
                    aes(x = x, y = y, alpha = 0.3)) + theme_minimal() +
      xlab("tau inferred") +
      ylab(" tau simulated") +
      ggtitle("Inferred tau hierarchical") +
      geom_point() 
    ggsave(paste0(sub_dir, "/fit_h_variational_BIC.png"), plot = p_fit, width = 10, height = 8)
    
    saveRDS(fits_var, paste0(sub_dir, "/fits_variational.rds"))
    saveRDS(ms_results_var, paste0(sub_dir, "/ms_results_variational.rds"))
    saveRDS(timing_var, paste0(sub_dir, "/timing_variational.rds"))
    
  ################################ Inference variational single segment #####################################
    
    m_single = cmdstanr::cmdstan_model(paste0(base_dir,"model/tickTack_noclust.stan"))
    
    start_time <- Sys.time()
    fit_single = m_single$variational(data = input_data, tol_rel_obj = 1e-4)
    end_time <- Sys.time()
    
    elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))

    tau_per_segment = fit_single$draws("tau", format = "matrix") %>% colMeans()
    
    timing_single <- rbind(data.frame(
      method = "mcmc",
      K = "all",
      time_seconds = elapsed
    ))
    
    p_fit_single <- ggplot(data = data.frame(x = tau_per_segment, y = sim$true_taus), 
                    aes(x = x, y = y, alpha = 0.3)) + theme_minimal() +
      xlab("tau inferred") +
      ylab(" tau simulated") +
      ggtitle("Inferred tau single segment") +
      geom_point() 
    
    p_final <- (p_fit + p_fit_single)
    ggsave(paste0(sub_dir, "/fit_HBIC_vs_single.png"), plot = p_final, width = 10, height = 4)
    

    ################################## save the merged results ######################################
    timing_all <- rbind(timing_var, timing_single)
    
    saveRDS(timing_all, paste0(sub_dir, "/timing_results.rds"))
  
    merged_res <- dplyr::tibble(
      segment_idx = seq_len(n_events),
      true_tau = sim$true_taus,
      true_tau_cluster = sim$taus_clust,
      tau_var_BIC <- ms_results_var$seg_assignment_BIC$tau_assignment,
      tau_var_AIC <- ms_results_var$seg_assignment_AIC$tau_assignment,
      tau_var_ICL <- ms_results_var$seg_assignment_ICL$tau_assignment,
      tau_var_LOO <- ms_results_var$seg_assignment_LOO$tau_assignment,
      
      tau_single = unname(tau_per_segment),

      failed_K_variational = paste(ms_results_var$failed_K, collapse = ";"),
      # failed_K_mcmc = paste(ms_results_mcmc$failed_K, collapse = ";"))
    )
    saveRDS(merged_res, paste0(sub_dir, "/inference_results.rds"))
    
  }
  
}
