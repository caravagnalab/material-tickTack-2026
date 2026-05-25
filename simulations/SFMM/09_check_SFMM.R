base_dir_SFMM = "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/simulations/SFMM/"
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

set.seed(123)

sim_5 <- readRDS("/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/simulations/simulations_5_components_40segments/results_simulations_5components/sim_5_20_0.8_40_2e-06/7/sim.rds")
sim_3 <- readRDS("/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/simulations/results_simulations/sim_3_20_0.4_70_1e-06/1/sim.rds")


source(paste0(base_dir, "/utils.R"))

MIN_MUTATIONS <- 10
K_max         <- 5
epsilon       <- 0.16
tolerance     <- 0.0001
max_attempts  <- 2
seed          <- 1234
set.seed(seed)


sim = sim_3
purity = 0.4

m5       <- cmdstanr::cmdstan_model(paste0(base_dir, "model/v5.stan"))
m_single <- cmdstanr::cmdstan_model(paste0(base_dir, "model/tickTack_noclust.stan"))

data       <- prepare_tickTack_input_data(sim, purity, MIN_MUTATIONS, alpha = 0)
input_data <- data$input_data



fits_and_time_var <- NULL
      fits_and_time_var <- model_selection_fit(
      K              = K_max,
      input_data     = input_data,
      m          = m5,
      inference_type = "variational"
    )

# ---- Model selection & summarized results ---------------------------
fits_var    <- fits_and_time_var$fits
timing_var  <- fits_and_time_var$timing


# check the effective number of components in loo estimation in LOO-CV approximation
ms_results_var  <- model_selection(fits = fits_var, input_data, K = K_max)
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
p_fit







# SFMM
m_SFMM <- cmdstanr::cmdstan_model(paste0(base_dir, "model/v1_SFMM.stan"))

data       <- prepare_tickTack_input_data(sim, purity, MIN_MUTATIONS, alpha = 0)
input_data <- data$input_data

tol_rel_obj    = 0.0001
vi_fits        = NULL
vi_best_k      = NULL       # best K from VI, required for mcmc_refined
n_chains       = 4
iter_warmup    = 1000
iter_sampling  = 1000
iter_warmup_refined = 500   # reduced warmup for VB-initialised MCMC
jitter_sd      = 0.05

stan_data <- c(input_data, list(alpha = 0.1, K = 20))
# fit_var <- m_SFMM$variational(data = stan_data, tol_rel_obj = tol_rel_obj)


fit_mcmc <- m_SFMM$sample(
  data             = stan_data,
  chains           = 2,
  iter_warmup      = 2000,   # generous warmup for sparse Dirichlet
  iter_sampling    = 1000,
  parallel_chains  = 2,
  adapt_delta      = 0.99,   # high adapt_delta to handle curved geometry
  max_treedepth    = 12      # allow deeper trees near boundaries
)

stan_data_orig <- c(input_data, list(K = 5))
fit_orig <- m5$sample(
  data            = stan_data_orig,
  chains          = 1,
  iter_warmup     = 100,
  iter_sampling   = 100,
  parallel_chains = 1
)


# table(input_data$seg_assignment)
# any(input_data$NV > input_data$DP)
# range(input_data$peaks)
# table(input_data$karyotype)
# range(input_data$seg_assignment)
# fit_diag <- m_SFMM$diagnose(data = stan_data)

safe_init <- lapply(1:2, function(i) list(
  tau = rep(0.5, 5),
  pi  = rep(0.2, 5)
))

fit_mcmc <- m_SFMM$sample(
  data            = stan_data,
  chains          = 2,
  iter_warmup     = 200,
  iter_sampling   = 200,
  parallel_chains = 2,
  adapt_delta     = 0.99,
  max_treedepth   = 12,
  init            = safe_init
)


fit_mcmc <- m_SFMM$sample(
  data            = stan_data,
  chains          = 2,
  iter_warmup     = 200,
  iter_sampling   = 200,
  parallel_chains = 2,
  adapt_delta     = 0.99,
  max_treedepth   = 12,
  init            = 0.1        # <-- much smaller init radius
)

# ---- Model selection & summarized results ---------------------------


# posterior mean of pi across draws
pi_draws <- fit_mcmc$draws("pi", format = "matrix")
pi_mean  <- colMeans(pi_draws)
print(round(pi_mean, 3))
# e.g. [0.421, 0.389, 0.187, 0.002, 0.001, 0.000, ...]
#       active  active  active  dead   dead   dead

# posterior mean of K_eff
k_eff_draws <- fit_mcmc$draws("K_eff", format = "matrix")
hist(k_eff_draws, breaks = 0:10)   # posterior distribution over K_eff


stan_data <- c(input_data, list(alpha = 0.1, K = 6))

fit_var <- m_SFMM$variational(data = stan_data, 
                              tol_rel_obj = tol_rel_obj, 
                              algorithm = "fullrank")
# posterior mean of pi across draws
pi_draws <- fit_var$draws("pi", format = "matrix")
pi_mean  <- colMeans(pi_draws)
print(round(pi_mean, 3))
# e.g. [0.421, 0.389, 0.187, 0.002, 0.001, 0.000, ...]
#       active  active  active  dead   dead   dead
tau_draws <- fit_var$draws("tau", format = "matrix")
tau_mean  <- colMeans(tau_draws)
print(round(tau_mean, 3))


# posterior mean of K_eff
k_eff_draws <- fit_var$draws("K_eff", format = "matrix")
hist(k_eff_draws, breaks = 40 )   # posterior distribution over K_eff











# use VI ad initialization for mcmc
# Step 1: VB for all K (as before)
fits_and_time_var <- model_selection_fit_init_mcmc(
  K              = K_max,
  input_data     = input_data,
  m              = m5,
  inference_type = "variational",
  tol_rel_obj    = tolerance
)
fits_var    <- fits_and_time_var$fits
ms_results_var <- model_selection(fits_var, input_data, K = K_max)
best_k_var  <- ms_results_var$best_K_BIC

# Step 2: MCMC only for best K ± 1, warm-started from VB
fits_and_time_mcmc <- model_selection_fit_init_mcmc(
  K                   = K_max,
  input_data          = input_data,
  m                   = m5,
  inference_type      = "mcmc_refined",
  vi_fits             = fits_var,
  vi_best_k           = best_k_var,
  n_chains            = 4,
  iter_warmup_refined = 500,
  iter_sampling       = 1000,
  jitter_sd           = 0.05
)
fits_mcmc <- fits_and_time_mcmc$fits

ms_results_mcmc = model_selection(fits_mcmc, input_data, K = 5)
best_k_mcmc_BIC = ms_results_mcmc$best_K_BIC
best_fit_mcmc_BIC = fits_mcmc[best_k_mcmc_BIC]
seg_ass_mcmc_BIC <- ms_results_mcmc$seg_assignment_BIC

p_fit <- ggplot(data = data.frame(x = seg_ass_mcmc_BIC$tau_assignment, y = sim$true_taus), 
                aes(x = x, y = y, alpha = 0.3)) + theme_minimal() +
  xlab("tau inferred") +
  ylab(" tau simulated") +
  geom_point() 




