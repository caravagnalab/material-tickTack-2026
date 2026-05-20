rm(list=ls())
# n_clocks = 2
# n_events = 10
# purity = 0.8
# coverage = 80
# mutation_density = 2.e-06

base_dir = "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/simulations/"

.libPaths("~/R/orfeo_R_4.4/")
library(tickTack)
library(patchwork)
library(cmdstanr)
library(posterior)
library(dplyr)
library(tidyr)
library(ggplot2)
library(loo)

args <- commandArgs(TRUE)
n_clocks        <- as.double(args[1])
n_events        <- as.double(args[2])
pi = purity     <- as.double(args[3])
coverage        <- as.double(args[4])
mutation_density <- as.numeric(args[5])

source(paste0(base_dir, "/utils.R"))

MIN_MUTATIONS <- 10
K_max         <- 3
epsilon       <- 0.16
tolerance     <- 0.0001
max_attempts  <- 2
seed          <- 1234
set.seed(seed)

print(paste0("n_clocks: ",        n_clocks))
print(paste0("n_events: ",        n_events))
print(paste0("purity: ",          purity))
print(paste0("coverage: ",        coverage))
print(paste0("mutation_density: ", mutation_density))
print(paste0("epsilon: ",         epsilon))
print(paste0("tolerance: ",       tolerance))
print(paste0("max_attempts: ",    max_attempts))
print(paste0("seed: ",            seed))

main_dir <- paste0(base_dir, "/results_simulations/sim_",
                   n_clocks, "_", n_events, "_", purity, "_", coverage, "_", mutation_density, "/")

if (!dir.exists(main_dir)) dir.create(main_dir)

m5       <- cmdstanr::cmdstan_model(paste0(base_dir, "model/v5.stan"))
m_single <- cmdstanr::cmdstan_model(paste0(base_dir, "model/tickTack_noclust.stan"))

for (i.iter in 1:10) {
  
  sub_dir <- paste0(main_dir, i.iter)
  if (dir.exists(sub_dir)) next
  dir.create(sub_dir)
  
  # ---- Simulate -------------------------------------------------------
  sim <- simulate_dataset(n_events, n_clocks, mutation_density, purity, coverage,
                          sigma_tau = .01, min_dist = epsilon)
  saveRDS(sim, paste0(sub_dir, "/sim.rds"))
  
  data       <- prepare_tickTack_input_data(sim, purity, MIN_MUTATIONS, alpha = 0)
  input_data <- data$input_data
  saveRDS(data, paste0(sub_dir, "/data.rds"))
  
  # ---- Variational fit with retry on failure --------------------------
  fits_and_time_var <- NULL
  attempt_seeds     <- c(seed, seed + 100, seed + 200)   # original + 2 retries
  
  for (attempt in seq_along(attempt_seeds)) {
    
    attempt_seed <- attempt_seeds[attempt]
    message(sprintf("[iter %d] model_selection_fit attempt %d (seed = %d)",
                    i.iter, attempt, attempt_seed))
    
    tryCatch({
      set.seed(attempt_seed)
      fits_and_time_var <- model_selection_fit(
        K              = K_max,
        input_data     = input_data,
        m          = m5,
        inference_type = "variational"
      )
    }, error = function(e) {
      message(sprintf("[iter %d] Attempt %d failed: %s", i.iter, attempt, conditionMessage(e)))
      fits_and_time_var <<- NULL
    })
    
    if (!is.null(fits_and_time_var)) break   # success — stop retrying
  }
  
  # If all attempts failed, skip to the next iteration
  if (is.null(fits_and_time_var)) {
    message(sprintf("[iter %d] All %d attempts failed — skipping iteration.",
                    i.iter, length(attempt_seeds)))
    next
  }
  
  # ---- Model selection & summarized results ---------------------------
  fits_var    <- fits_and_time_var$fits
  timing_var  <- fits_and_time_var$timing
  
  ms_results_var  <- model_selection(fits_var, input_data, K = K_max)
  best_k_var_BIC  <- ms_results_var$best_K_BIC
  seg_ass_var_BIC <- ms_results_var$seg_assignment_BIC
  
  fits_and_time_var_all_K <- build_summarized_results(
    fits_and_time_var = fits_and_time_var,
    sim               = sim,
    K_max             = K_max
  )
  
  # ---- Build x object -------------------------------------------------
  cna_renamed <- sim$cn %>%
    dplyr::rename(from = startpos, to = endpos, Major = nMaj1_A, minor = nMin1_A) %>%
    mutate(CCF = 1)
  
  x <- list(
    mutations        = sim$muts %>% dplyr::rename(from = start, to = end),
    cna              = cna_renamed,
    metadata         = tibble(purity = purity),
    reference_genome = "hg19",
    results_timing   = list(
      draws_and_summary = fits_and_time_var_all_K$draws_and_summary,
      data = list(
        accepted_cna = cna_renamed %>%
          mutate(segment_name = paste0(chr, "_", from, "_", to)),
        input_data   = data$input_data
      )
    )
  )
  
  # ---- Plots ----------------------------------------------------------
  pdf(paste0(sub_dir, "/sim_plot.pdf"), width = 20, height = 8)
  
  print(plot_sim(sim))
  
  p <- tickTack::plot_cnaqc_choose_K(x, best_k_var_BIC, add_VAF_ecdf = TRUE)
  print(p)
  p <- tickTack::plot_cnaqc_choose_K(x, best_k_var_BIC, add_VAF_hist = TRUE)
  print(p)
  
  p_fit <- ggplot(
    data = data.frame(x = seg_ass_var_BIC$tau_assignment, y = sim$true_taus),
    aes(x = x, y = y, alpha = 0.3)
  ) +
    theme_minimal() +
    xlab("tau inferred") + ylab("tau simulated") +
    ggtitle("Inferred tau hierarchical") +
    geom_point()
  
  # ---- Single-segment model -------------------------------------------
  start_time <- Sys.time()
  fit_single <- m_single$variational(data = input_data, tol_rel_obj = 1e-4)
  end_time   <- Sys.time()
  elapsed    <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  tau_per_segment <- fit_single$draws("tau", format = "matrix") %>% colMeans()
  
  timing_single <- data.frame(method = "mcmc", K = "all", time_seconds = elapsed)
  
  p_fit_single <- ggplot(
    data = data.frame(x = tau_per_segment, y = sim$true_taus),
    aes(x = x, y = y, alpha = 0.3)
  ) +
    theme_minimal() +
    xlab("tau inferred") + ylab("tau simulated") +
    ggtitle("Inferred tau single segment") +
    geom_point()
  
  print(p_fit + p_fit_single)
  dev.off()
  
  # ---- Save results ---------------------------------------------------
  timing_all <- rbind(timing_var, timing_single)
  
  saveRDS(fits_var,          paste0(sub_dir, "/fits_variational.rds"))
  saveRDS(ms_results_var,    paste0(sub_dir, "/ms_results_variational.rds"))
  saveRDS(timing_var,        paste0(sub_dir, "/timing_variational.rds"))
  saveRDS(timing_all,        paste0(sub_dir, "/timing_results.rds"))
  
  merged_res <- dplyr::tibble(
    segment_idx              = seq_len(n_events),
    true_tau                 = sim$true_taus,
    true_tau_cluster         = sim$taus_clust,
    tau_var_BIC              = ms_results_var$seg_assignment_BIC$tau_assignment,
    tau_var_AIC              = ms_results_var$seg_assignment_AIC$tau_assignment,
    tau_var_ICL              = ms_results_var$seg_assignment_ICL$tau_assignment,
    tau_var_LOO              = ms_results_var$seg_assignment_LOO$tau_assignment,
    tau_single               = unname(tau_per_segment),
    failed_K_variational     = paste(ms_results_var$failed_K, collapse = ";")
  )
  saveRDS(merged_res, paste0(sub_dir, "/inference_results.rds"))
  
  message(sprintf("[iter %d] Done.", i.iter))
}