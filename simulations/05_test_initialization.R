base_dir = "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/simulations/"
#!/usr/bin/env Rscript
.libPaths("~/R/orfeo_R_4.4/")
source(paste0(base_dir,"/utils.R"))
library(tickTack)

n_clocks = 5
n_events = 5
purity = 0.8
coverage = 80
mutation_density  = 1.8e-06
iter = 1

MIN_MUTATIONS = 10
tol_rel_obj = 0.0001
  
main_dir = paste0(base_dir,"results_simulations/sim_",
                  n_clocks, "_", n_events, "_", purity, "_", coverage, "_", mutation_density, "/")

sim <- readRDS(paste0(main_dir,iter,"/sim.rds"))
data <- prepare_tickTack_input_data(sim, purity, alpha = 0, MIN_MUTATIONS) 
input_data = data$input_data

m5 = cmdstanr::cmdstan_model(paste0(base_dir,"model/v5.stan"))
m = m5

m5_few_segs = cmdstanr::cmdstan_model(paste0(base_dir,"model/v5_few_segs.stan"))


k = 5
# K = 5
# candidate_Ks = 1:K
# 
# for (k in candidate_Ks) {
    input_data$K = k
    
    init_var = get_initialization_v5(input_data, purity)
    # init_var <- init_var %>% 
      
    fit_no_init = m$variational(data = input_data, 
                          tol_rel_obj = tol_rel_obj)
    
    seg_ass_no_init = extract_seg_assignment(fit_no_init, K_sel = 4, S = 4)
    
    stan_data = input_data
    stan_data$alpha_dirichlet <- 0.5   # sparse: favours peaked pi
    stan_data$gamma_local     <- 0.5
    
    fit = m5_few_segs$variational(data = stan_data, 
                        tol_rel_obj = tol_rel_obj)
    
    fit <- m5_few_segs$sample(
      data          = stan_data,
      # init          = list(init_var),   # your good initialization
      chains        = 4,
      iter_warmup   = 1000,
      iter_sampling = 1000,
      adapt_delta   = 0.95,             # higher for mixture models
      max_treedepth = 12
    )
    
    seg_ass_init = extract_seg_assignment(fit, K_sel = 4, S = 4)
    seg_ass_init

    fits[[k]] <- fit
    
    # timing <- rbind(timing, data.frame(
    #   method = "variational",
    #   K = k,
    #   time_seconds = elapsed
    # ))
  # }



ms_results_var = model_selection(fits_var, input_data, K = 5)
