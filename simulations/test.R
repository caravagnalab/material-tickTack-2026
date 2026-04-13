
library(cmdstanr)
library(posterior)
library(dplyr)
library(tidyr)
library(ggplot2)

base_path = "/orfeo/cephfs/scratch/cdslab/scocomello/material_tickTack/simulations_v3/results_simulations/sim_3_30_0.9_40_1e-06/1"
base_path = "/orfeo/cephfs/scratch/cdslab/scocomello/material_tickTack/simulations_v3/results_simulations/sim_5_30_0.9_80_7e-06/1/"

sim = readRDS(file.path(base_path, "sim.rds"))
data = readRDS(file.path(base_path, "res_tickTack_h.rds"))
input_data = data$results$data$input_data

s_id = 15
input_data$peaks[[s_id]]
idxs = which(input_data$seg_assignment == s_id)
input_data$karyotype[[s_id]]
input_data$peaks[[s_id]]
sim$true_taus[s_id]

VAF = input_data$NV[idxs] / input_data$DP[idxs]
hist(VAF, breaks = 20)

cut = .7
N2 = sum(VAF > .7)
N1 = sum(VAF <= .7)

N1 / (N1 + N2)
N2 / (N1 + N2)

hist(input_data$NV[idxs] / input_data$DP[idxs])

sim$true_taus %>% hist(breaks = 100)

# m0 = cmdstanr::cmdstan_model("v0.stan")
# m1 = cmdstanr::cmdstan_model("v1.stan")
# m2 = cmdstanr::cmdstan_model("v2.stan")
# m3 = cmdstanr::cmdstan_model("v3.stan")
# m4 = cmdstanr::cmdstan_model("v4.stan")
m5 = cmdstanr::cmdstan_model("v5.stan")


#input_data$gamma_ord = 10

m = m5

input_data$K = 5

fit = m$variational(data = input_data, tol_rel_obj = 0.0001)
# fit = m$sample(data = input_data, chains = 1, iter_warmup = 1000, iter_sampling = 1000)
# fit = m$pathfinder(data = input_data, num_paths = 4, num_threads = 4, tol_rel_obj = 0.0001)

tau_draws    <- fit$draws("tau", format = "matrix") 

sim$taus_clust %>% table()

tau_draws %>% colMeans() %>% sort()
sim$taus_clust %>% unique() %>% sort()
# tau_draws = draws |>
#   select(starts_with("tau[")) |>
#   summarise(across(everything(), mean))
# 
# tau_draws
# sim$taus_clust %>% unique() %>% sort()



# Model selection

library(loo)

# Fit models for K = 1, 2, 3, 4, ...
candidate_Ks = 1:5
fits <- list()

for (k in candidate_Ks) {
  input_data$K = k
  fit = m$variational(data = input_data, tol_rel_obj = 0.0001)
  #fit = m$sample(data = input_data, chains = 1, iter_warmup = 1000, iter_sampling = 1000)
  fits[[k]] <- fit
}

bic_from_cmdstan <- function(fit, k, n_obs = S) {
  log_lik_mat <- fit$draws("log_lik", format = "matrix")  # (draws, S)
  
  # Total log-lik per draw = sum over segments
  total_ll_per_draw <- rowSums(log_lik_mat)
  mean_ll <- mean(total_ll_per_draw)
  
  n_params <- 2 * k - 1  # (K-1) for pi + K for tau
  bic <- -2 * mean_ll + n_params * log(n_obs)
  bic
}

S = input_data$S
bics <- sapply(candidate_Ks, \(k) bic_from_cmdstan(fits[[k]], k, n_obs = S))

K = best_K = which.min(bics)
#K = best_K = 5

fit = fits[[best_K]]

# Extract draws
tau_draws    <- fit$draws("tau", format = "matrix")          # (draws, K)
probs_draws  <- fit$draws("seg_probs", format = "matrix")    # (draws, S*K) — needs reshape
n_draws <- nrow(tau_draws)

# Reshape seg_probs to (draws, S, K)
probs_arr <- array(probs_draws, dim = c(n_draws, S, K))

# Average probabilities over draws → (S, K)
mean_probs  <- apply(probs_arr, c(2, 3), mean)

# Hard assignment per segment
hard_assign <- apply(mean_probs, 1, which.max)   # length S

# Mean tau per clock
mean_tau <- apply(tau_draws, 2, mean)             # length K
sim$taus_clust %>% unique() %>% sort()
sort(mean_tau)

# Tau assigned to each segment
seg_tau <- mean_tau[hard_assign]                  # length S

plot(seg_tau, sim$true_taus)

abs(seg_tau - sim$true_taus) %>% which.max()

s_id = 13
idxs = which(input_data$seg_assignment == s_id)
input_data$karyotype[[s_id]]
input_data$peaks[[s_id]]
sim$true_taus[s_id]

VAF = input_data$NV[idxs] / input_data$DP[idxs]
hist(VAF, breaks = 20)
