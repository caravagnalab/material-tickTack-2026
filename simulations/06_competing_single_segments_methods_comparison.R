# rm(list=ls())
# n_clocks = 1
# n_events = 5
# purity = 0.4
# coverage = 40
# mutation_density = 2.e-07
base_dir = "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/simulations/"

.libPaths("~/R/orfeo_R_4.4/")
library(tickTack)
library(patchwork)
library(cmdstanr)
library(posterior)
library(dplyr)
library(tidyr)
library(ggplot2)

args <- commandArgs(TRUE)
n_clocks         <- as.double(args[1])
n_events         <- as.double(args[2])
pi = purity      <- as.double(args[3])
coverage         <- as.double(args[4])
mutation_density <- as.numeric(args[5])

source(paste0(base_dir, "/utils.R"))

MIN_MUTATIONS <- 10
epsilon       <- 0.16
tolerance     <- 0.0001
max_attempts  <- 2
seed          <- 1234
set.seed(seed)

print(paste0("n_clocks: ",         n_clocks))
print(paste0("n_events: ",         n_events))
print(paste0("purity: ",           purity))
print(paste0("coverage: ",         coverage))
print(paste0("mutation_density: ", mutation_density))
print(paste0("epsilon: ",          epsilon))
print(paste0("tolerance: ",        tolerance))
print(paste0("max_attempts: ",     max_attempts))
print(paste0("seed: ",             seed))

main_dir <- paste0(base_dir, "/results_simulations/sim_",
                   n_clocks, "_", n_events, "_", purity, "_", coverage, "_", mutation_density, "/")

if (dir.exists(main_dir)) {
  
  for (i.iter in 1:10) {
    
    sub_dir <- paste0(main_dir, i.iter)
    
    if (!dir.exists(sub_dir)) next
    
    sim <- readRDS(paste0(sub_dir, "/sim.rds"))
    
    mult_path <- paste0(sub_dir, "/sample_1_dpclust_info")
    
    get_multiplicities(sim, purity, mult_path, sub_dir)
    res_AmpTimeR <- safe_run(fit_AmpTimeR(sim, mult_path), "fit_AmpTimeR")
    res_MutTime  <- safe_run(fit_MutTimeR(sim, pi),        "fit_MutTimeR")
    
    if (!is.null(res_AmpTimeR) && !is.null(res_MutTime)) {
      
      res <- readRDS(paste0(sub_dir, "/inference_results.rds"))
      
      merged_res <- dplyr::tibble(
        segment_idx      = seq_len(n_events),
        true_tau         = sim$true_taus,
        true_tau_cluster = sim$taus_clust,
        
        tau_AmpTimeR = if (!is.null(res_AmpTimeR)) {
          res_AmpTimeR$tau
        } else {
          rep(NA_real_, n_events)
        },
        
        tau_MutTimeR = if (!is.null(res_MutTime)) {
          res_MutTime$cn_timed$time
        } else {
          rep(NA_real_, n_events)
        }
      )
      
      merged_res <- merged_res %>%
        left_join(res, by = join_by(true_tau, true_tau_cluster, segment_idx))
      
      saveRDS(merged_res, file.path(sub_dir, "merged_res_competing_methods.rds"))
      
      # ---- Comparison plots -> temp pdf, then merge with existing --------
      existing_pdf <- paste0(sub_dir, "/sim_plot.pdf")
      temp_pdf     <- paste0(sub_dir, "/comparison_plots.pdf")
      
      pdf(temp_pdf, width = 20, height = 8)
      
      # AmpTimeR
      if (!is.null(res_AmpTimeR)) {
        p_AmpTimeR <- ggplot(
          data = data.frame(x = res_AmpTimeR$tau, y = sim$true_taus),
          aes(x = x, y = y, alpha = 0.3)
        ) +
          theme_minimal() +
          xlab("tau inferred") + ylab("tau simulated") +
          ggtitle("Inferred tau — AmpTimeR") +
          geom_point()
      } else {
        p_AmpTimeR <- ggplot() +
          theme_minimal() +
          ggtitle("AmpTimeR — failed")
      }
      
      # MutTimeR
      if (!is.null(res_MutTime)) {
        p_MutTimeR <- ggplot(
          data = data.frame(x = res_MutTime$cn_timed$time, y = sim$true_taus),
          aes(x = x, y = y, alpha = 0.3)
        ) +
          theme_minimal() +
          xlab("tau inferred") + ylab("tau simulated") +
          ggtitle("Inferred tau — MutTimeR") +
          geom_point()
      } else {
        p_MutTimeR <- ggplot() +
          theme_minimal() +
          ggtitle("MutTimeR — failed")
      }
      
      print(p_AmpTimeR + p_MutTimeR)
      
      dev.off()
      
      # # Merge with existing sim_plot.pdf
      # if (file.exists(existing_pdf)) {
      #   qpdf::pdf_combine(
      #     input  = c(existing_pdf, temp_pdf),
      #     output = existing_pdf
      #   )
      #   file.remove(temp_pdf)
      #   message(sprintf("[iter %d] Merged comparison plots into %s", i.iter, existing_pdf))
      # } else {
      #   # No existing pdf — just rename temp as the main one
      #   file.rename(temp_pdf, existing_pdf)
      #   message(sprintf("[iter %d] No existing pdf found — saved comparison plots as %s", i.iter, existing_pdf))
      # }
    }
  }
}