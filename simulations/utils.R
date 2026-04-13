.libPaths("~/R/rstudio_v3/")
library(AmplificationTimeR)
library(BSgenome.Hsapiens.UCSC.hg19)
library(tidyverse)
library(VariantAnnotation)
library(stringr)
require(transport)
library(cluster)
#library(BSgenome.Hsapiens.UCSC.hg38)
library(tickTack)
library(dpclust3p)
library(fossil)



model_selection_fit <- function(K = 5, input_data,  m, inference_type = "variational", tol_rel_obj = 0.0001){
  # Fit models for K = 1, 2, 3, 4, 5
  candidate_Ks = 1:K
  fits <- list()
  
  timing <- data.frame(
    method = character(),
    K = integer(),
    time_seconds = numeric(),
    stringsAsFactors = FALSE
  )
  
  total_time <- 0
  
  if(inference_type == "variational"){
    for (k in candidate_Ks) {
      input_data$K = k
      
      start_time <- Sys.time()
      
      fit = m$variational(data = input_data, tol_rel_obj = tol_rel_obj)

      end_time <- Sys.time()
      elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
      total_time <- total_time + elapsed
      
      fits[[k]] <- fit
      
      timing <- rbind(timing, data.frame(
        method = "variational",
        K = k,
        time_seconds = elapsed
      ))
    }
  } else {
    for (k in candidate_Ks) {
      input_data$K = k
      
      start_time <- Sys.time()
      
      fit = m$sample(data = input_data, chains = 1, iter_warmup = 1000, iter_sampling = 1000)
      
      end_time <- Sys.time()
      elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
      total_time <- total_time + elapsed
      
      fits[[k]] <- fit
      
      timing <- rbind(timing, data.frame(
        method = "mcmc",
        K = k,
        time_seconds = elapsed
      ))
    }
  }
  
  # Add total row
  timing <- rbind(timing, data.frame(
    method = paste0(inference_type, "_total"),
    K = NA,
    time_seconds = total_time
  ))
  
  return(list(
    fits = fits,
    timing = timing
  ))
  
}


model_selection <- function(fits, input_data, K) {
  
  library(loo)
  
  candidate_Ks <- 1:K
  S <- input_data$S
  
  valid_Ks <- c()
  failed_Ks <- c()
  
  results_list <- list()
  
  # ---------------- Loop ----------------
  
  for (k in candidate_Ks) {
    fit <- fits[[k]]
    is_valid <- TRUE
    
    if (is.null(fit)) {
      is_valid <- FALSE
    } else {
      tryCatch({
        log_lik_mat <- fit$draws("log_lik", format = "matrix")
        if (is.null(log_lik_mat) || nrow(log_lik_mat) == 0) {
          is_valid <- FALSE
        }
      }, error = function(e) {
        is_valid <<- FALSE
      })
    }
    
    if (!is_valid) {
      failed_Ks <- c(failed_Ks, k)
      next
    }
    
    valid_Ks <- c(valid_Ks, k)
    
    # ---------------- Log-likelihood ----------------
    
    log_lik_mat <- fit$draws("log_lik", format = "matrix")
    
    # mean log-likelihood (standard Bayesian approximation)
    total_ll <- rowSums(log_lik_mat)
    mean_ll <- mean(total_ll)
    
    n_params <- 2 * k - 1
    
    # ---------------- AIC / BIC ----------------
    
    AIC_val <- -2 * mean_ll + 2 * n_params
    BIC_val <- -2 * mean_ll + n_params * log(S)
    
    # ---------------- LOO (standard) ----------------
    
    loo_res <- loo(log_lik_mat)
    elpd_loo <- loo_res$estimates["elpd_loo", "Estimate"]
    
    # ---------------- ICL ----------------
    
    probs_draws <- fit$draws("seg_probs", format = "matrix")
    n_draws <- nrow(probs_draws)
    
    probs_arr <- array(probs_draws, dim = c(n_draws, S, k))
    mean_probs <- apply(probs_arr, c(2, 3), mean)
    
    mean_probs[mean_probs == 0] <- 1e-12
    entropy <- -sum(mean_probs * log(mean_probs))
    
    ICL_val <- BIC_val + 2 * entropy
    
    # ---------------- Store ----------------
    
    results_list[[length(results_list) + 1]] <- data.frame(
      K = k,
      AIC = AIC_val,
      BIC = BIC_val,
      ICL = ICL_val,
      LOO = elpd_loo
    )
  }
  
  if (length(valid_Ks) == 0) {
    stop("All fits failed.")
  }
  
  criteria <- do.call(rbind, results_list)
  
  # ---------------- Select best ----------------
  
  best_K_AIC <- criteria$K[which.min(criteria$AIC)]
  best_K_BIC <- criteria$K[which.min(criteria$BIC)]
  best_K_ICL <- criteria$K[which.min(criteria$ICL)]
  best_K_LOO <- criteria$K[which.max(criteria$LOO)]
  
  # ---------------- Extract assignments ----------------
  
  extract_seg_assignment <- function(fit, K_sel) {
    
    tau_draws   <- fit$draws("tau", format = "matrix")
    probs_draws <- fit$draws("seg_probs", format = "matrix")
    
    n_draws <- nrow(tau_draws)
    probs_arr <- array(probs_draws, dim = c(n_draws, S, K_sel))
    
    mean_probs <- apply(probs_arr, c(2, 3), mean)
    hard_assign <- apply(mean_probs, 1, which.max)
    
    mean_tau <- apply(tau_draws, 2, mean)
    
    dplyr::tibble(
      tau_assignment = unlist(mean_tau[hard_assign]),
      segment_idx = 1:S
    )
  }
  
  return(list(
    best_K_AIC = best_K_AIC,
    best_K_BIC = best_K_BIC,
    best_K_ICL = best_K_ICL,
    best_K_LOO = best_K_LOO,
    
    seg_assignment_AIC = extract_seg_assignment(fits[[best_K_AIC]], best_K_AIC),
    seg_assignment_BIC = extract_seg_assignment(fits[[best_K_BIC]], best_K_BIC),
    seg_assignment_ICL = extract_seg_assignment(fits[[best_K_ICL]], best_K_ICL),
    seg_assignment_LOO = extract_seg_assignment(fits[[best_K_LOO]], best_K_LOO),
    
    valid_K = valid_Ks,
    failed_K = failed_Ks,
    
    criteria = criteria
  ))
}





prepare_tickTack_input_data = function(sim, pi, min_mutations_number, alpha) {
  N_events = nrow(sim$cn)
  
  cn <- sim$cn %>% 
    dplyr::rename(Major=nMaj1_A, minor=nMin1_A, from=startpos, to=endpos) %>% 
    dplyr::mutate(CCF = 1)
  
  muts <- sim$muts %>% 
    dplyr::rename(from=start, to=end) %>% 
    dplyr::mutate(CCF = 1)
  
  x = list(
    cna = cn, 
    mutations = muts,
    metadata = dplyr::tibble(sample = "sample", purity=pi)
  )
  
  mutations <- tickTack:::Mutations(x)
  
  if(!nrow(mutations)*ncol(mutations)){
    stop("No mutations have been called on this CNAqc object.")
  }

  # cna <- CNAqc::CNA(x)
  segments <- tickTack:::CNA(x)
  
  if(!nrow(segments)*ncol(segments)){
    stop("No CNA events have been called on this CNAqc object.")
  }
  
  # temporarly set the purity here or give it in input before implementing a getter for the purity
  purity = x$metadata$purity
  
  accepted_data <- tickTack:::prepare_input_data(mutations, segments, purity, possible_k = c("2:1", "2:2", "2:0"), alpha = alpha, min_mutations_number = min_mutations_number)

  accepted_data
}

# accuracy score evaluation

parse_vec <- function(x) {
  as.numeric(strsplit(x, "_")[[1]])
}

pairwise_distortion <- function(real, inf) {
  D_real <- as.matrix(dist(real))
  D_inf  <- as.matrix(dist(inf))
  mean(abs(D_real - D_inf))
}

score_clocks <- function(real_str, inferred_str) {
  
  real <- parse_vec(real_str)
  inf  <- parse_vec(inferred_str)
  
  # --- Value accuracy ---
  mae  <- mean(abs(real - inf))
  nmae <- mae / (max(real) - min(real) + 1e-12)
  
  # --- Pairwise temporal distortion ---
  ptd <- pairwise_distortion(real, inf)
  
  # --- Ordering logic ---
  sd_real <- sd(real)
  sd_inf  <- sd(inf)
  
  if (sd_real == 0 && sd_inf == 0) {
    spearman <- 1
  } else if (sd_real == 0 && sd_inf > 0) {
    spearman <- 0
  } else if (sd_real > 0 && sd_inf == 0) {
    spearman <- 0
  } else {
    spearman <- cor(real, inf, method = "spearman")
  }
  
  return(c(MAE = mae,
           NMAE = nmae,
           Spearman = spearman,
           PTD = ptd))
}









# external evaluation metric 
build_truth_matrix <- function(true_times, tol = 0.02) {
  D <- as.matrix(dist(true_times))
  M_true <- D <= tol
  diag(M_true) <- TRUE
  M_true
}

build_inferred_matrix <- function(lower, upper) {
  n <- length(lower)
  
  low_max  <- pmax(outer(lower, lower, pmax))
  high_min <- pmin(outer(upper, upper, pmin))
  
  overlap <- high_min - low_max
  M_inf <- overlap >= 0
  diag(M_inf) <- TRUE
  M_inf
}

structure_score <- function(M_true, M_inf) {
  n <- nrow(M_true)
  idx <- upper.tri(M_true)
  
  agree <- sum(M_true[idx] == M_inf[idx])
  total <- length(M_true[idx])
  
  agree / total
}

time_error <- function(true_times, mean_inf) {
  mean(abs(true_times - mean_inf))
}

false_merge_penalty <- function(true_times, lower, upper, tol = 0.02) {
  
  D_true <- as.matrix(dist(true_times))
  far_true <- D_true > tol
  
  low_max  <- pmax(outer(lower, lower, pmax))
  high_min <- pmin(outer(upper, upper, pmin))
  overlap  <- (high_min - low_max) >= 0
  
  idx <- upper.tri(far_true)
  
  # coppie che nella verità sono lontane
  far_pairs <- far_true[idx]
  
  # CASO K = 1 (o tempi tutti vicini)
  # Nessuna coppia lontana → nessun false merge possibile
  if (!any(far_pairs)) return(0)
  
  # tra le coppie lontane vere, quante vengono fuse nell'inferenza
  false_merges <- overlap[idx][far_pairs]
  
  mean(false_merges)
}

evaluate_clustering <- function(true_times, mean_inf, lower, upper, tol = 0.02) {
  
  M_true <- build_truth_matrix(true_times, tol)
  M_inf  <- build_inferred_matrix(lower, upper)
  
  structure <- structure_score(M_true, M_inf)
  err_time  <- time_error(true_times, mean_inf)
  false_merges <- false_merge_penalty(true_times, lower, upper, tol)
  
  list(
    structure_agreement = structure,
    mean_time_error = err_time,
    false_merge_rate = false_merges
  )
}

total_loss <- function(eval, w1=1, w2=5, w3=5){
  (1 - eval$structure_agreement) +
    w1 * eval$mean_time_error +
    w2 * eval$false_merge_rate
}

build_index <- function(sim, result_inference){
  info_k <- result_inference$results$draws_and_summary[[2]]$summarized_results
  
  clustering_evaluation <- evaluate_clustering(
    true_times = sim$taus_clust,  # tempi veri
    mean_inf   = info_k$clock_mean,
    lower      = info_k$clock_low,
    upper      = info_k$clock_high
  )
  
  loss <- total_loss(clustering_evaluation)
  
  return(list(clustering_evaluation = clustering_evaluation,
              loss = loss))
}


#### max overlap min dist --> internal evaluation metric 

process_file <- function(s, return_res = FALSE, save_res = TRUE) {
  tryCatch({
    # print(s)
    original_dir = dirname(s)
    
    # if(file.exists(paste0(dirname(s), "/results_model_selection.rds"))){
    #   return(NULL)
    # }
    
    new_dir = paste0(dirname(s),"/plot")
    
    # if(file.exists(paste0(new_dir),"/plot_timing_all_K_h.png")){
    #   return(NULL)
    # }
    
    fit <- readRDS(paste0("",s,""))
    sim <- readRDS(paste0(dirname(s),"/sim.rds"))
    
    
    results <- fit$results
    results_model_selection <- fit$results_model_selection
    summarized_results <- results_model_selection$best_fit$summarized_results
    
    best_K <- results_model_selection$best_K
    model_selection_tibble <- results_model_selection$model_selection_tibble
    entropy <- results_model_selection$entropy_list
    
    K = nrow(results_model_selection$model_selection_tibble)
    
    #########################################################
    cna <- results$data$accepted_cna %>%
      rowwise() %>%
      mutate( from = strsplit(segment_name,split="_")[[1]][2], to = strsplit(segment_name,split="_")[[1]][3]) %>%
      mutate( len = as.numeric(to) - as.numeric(from))
    n_mut <- length(results$data$input_data$NV)
    mean_density_mut <- n_mut/sum(cna$len)
    
    results_model_selection$model_selection_tibble$sample_id = paste0(basename(dirname(dirname(s))),":",basename(dirname(s)))
    results_model_selection$model_selection_tibble$mean_density_mut = mean_density_mut
    #########################################################
    
    
    if (!dir.exists(new_dir)) {
      if( save_res == TRUE){
        dir.create(new_dir)
      }
    }
    
    # info on mutations per segment evaluated only on the segments used for the fit,
    # but includes some mut that have been filtered out, very few
    # accepted_cna <- fit$results$data$accepted_cna
    # mut <- as.data.table(x$mutations)
    # cn  <- as.data.table(accepted_cna)
    # cn[, c("chr", "from", "to") := tstrsplit(segment_name, "_", type.convert = TRUE)]
    #
    # mut[, chr := as.character(chr)]
    # cn[, chr := as.character(chr)]
    #
    # setkey(mut, chr, from, to)
    # setkey(cn,  chr, from, to)
    #
    # mut_in_cn <- foverlaps(mut, cn, nomatch = 0)
    # mut_per_segment <- mut_in_cn[, .N, by = segment_id]
    # colnames(mut_per_segment)[2] <- "n_mutations"
    #
    # all_segments <- cn[, .(segment_id)]
    #
    # mut_per_segment <- merge(all_segments,
    #                          mut_per_segment,
    #                          by = "segment_id",
    #                          all.x = TRUE)
    #
    # mut_per_segment[is.na(n_mutations), n_mutations := 0]
    # median_mut_per_segment <- median(mut_per_segment$n_mutations)
    # mean_mut_per_segment <- mean(mut_per_segment$n_mutations)
    total_mut_filtered <- length(results$data$input_data$NV)
    
    results_model_selection$model_selection_tibble$total_mut_filtered = total_mut_filtered
    
    # get the distances needed to evaluate the goodness of the fit wrt the number of mutations
    Kmax <- nrow(results_model_selection$model_selection_tibble)
    
    mean_sep  <- numeric(Kmax)
    mean_ov   <- numeric(Kmax)
    min_sep   <- numeric(Kmax)
    max_ov    <- numeric(Kmax)
    RI    <- numeric(Kmax)
    tau_inferred    <- numeric(Kmax)
    tau_inferred_low    <- numeric(Kmax)
    tau_inferred_high    <- numeric(Kmax)
    entropy_per_segment <- numeric(Kmax)
    entropy_per_segment_norm <- numeric(Kmax)
    
    for (i in seq_len(Kmax)){
      
      # print(paste0(i))
      best_k = i
      best_fit = fit$results$draws_and_summary[[i]]$summarized_results
      row = results_model_selection$model_selection_tibble[i,]
      
      
      unique_intervals <- unique(best_fit %>%
                                   dplyr::select(clock_low, clock_high, clock_mean))
      
      K <- nrow(unique_intervals)
      
      ## RI ##
      inferred_clusters = best_fit$clock_mean
      simulated_clusters <- sim$taus_clust
      RI[i] <- rand.index(inferred_clusters,simulated_clusters)
      tau_inferred[i] <- paste(inferred_clusters, collapse = "_")
      tau_inferred_low[i] <- paste(best_fit$clock_low, collapse = "_")
      tau_inferred_high[i] <- paste(best_fit$clock_high, collapse = "_")
      
      if (is.null(K) || K < 2) {
        mean_sep[i] <- NA
        mean_ov[i]  <- NA
        min_sep[i]  <- NA
        max_ov[i]   <- NA
        next
      }
      
      intervals <- unique_intervals %>%
        dplyr::mutate(
          mean = clock_mean,
          low  = clock_low,
          high = clock_high
        )
      
      dist_mean_mat     <- matrix(NA_real_, K, K)
      dist_overlap_mat  <- matrix(NA_real_, K, K)
      
      mean_vec  <- intervals$mean
      low_vec   <- intervals$low
      high_vec  <- intervals$high
      
      # Pairwise mean distances
      dist_mean_mat <- abs(outer(mean_vec, mean_vec, "-"))
      
      # Overlap
      low_max  <- pmax(outer(low_vec,  low_vec,  pmax))
      high_min <- pmin(outer(high_vec, high_vec, pmin))
      overlap  <- pmax(0, high_min - low_max)
      
      union_len <- pmax(outer(high_vec, high_vec, pmax)) -
        pmin(outer(low_vec,  low_vec,  pmin))
      
      dist_overlap_mat <- ifelse(union_len == 0, 0, overlap / union_len)
      
      diag(dist_mean_mat)    <- 0
      diag(dist_overlap_mat) <- 1
      
      
      closest_mean_distance <- apply(dist_mean_mat, 1, function(x) sort(x)[2])
      closest_overlap       <- apply(dist_overlap_mat, 1, function(x) sort(x, decreasing=TRUE)[2])
      
      mean_sep[i] <- mean(closest_mean_distance)
      mean_ov[i]  <- mean(closest_overlap)
      min_sep[i]  <- min(closest_mean_distance)
      max_ov[i]   <- max(closest_overlap)
      entropy_per_segment[i] <- entropy$entropy_per_segment_matrix[as.character(i)][[1]][[1]][1]
      entropy_per_segment_norm[i] <- entropy$entropy_per_segment_matrix_norm[as.character(i)][[1]][[1]][1]
      
    }
    
    tbl <- results_model_selection$model_selection_tibble
    tbl$mean_cluster_separation <- mean_sep
    tbl$mean_cluster_overlap    <- mean_ov
    tbl$min_cluster_separation  <- min_sep
    tbl$max_cluster_overlap     <- max_ov
    tbl$RI <- RI
    tbl$real_clock <- paste(sim$true_taus, collapse = "_")
    tbl$clust_clock <- paste(sim$taus_clust, collapse = "_")
    tbl$inferred_clock <- tau_inferred
    tbl$inferred_low <- tau_inferred_low
    tbl$inferred_high <- tau_inferred_high
    tbl$entropy_per_segment <- entropy_per_segment
    tbl$entropy_per_segment_norm <- entropy_per_segment_norm
    
    if (save_res == TRUE){
      saveRDS(tbl, paste0(original_dir,"/results_model_selection.rds"))
    }
    # best_K <- results_model_selection$best_K
    # model_selection_tibble <- results_model_selection$model_selection_tibble
    # entropy <- results_model_selection$entropy_list
    #
    # K = nrow(results_model_selection$model_selection_tibble)
    #
    # p_elbo <- list()
    # for (i in 1:K){
    #   p_elbo[[i]] <- tickTack::plot_elbo_h(results$elbo_iterations[[as.character(i)]]) + ggplot2::ggtitle(paste0("K = ", i))
    # }
    # p_elbo <- gridExtra::grid.arrange(grobs = p_elbo, ncol = 2)  #add global title
    # ggplot2::ggsave(paste0(new_dir,"/plot_elbo.png"),plot = p_elbo, width = 30, height = 30)
    #
    #
    #
    # p <- tickTack::plot_timing_h(results, best_K)
    # ggsave(paste0(new_dir,"/plot_timing_h.png"),plot = p, width = 25, height = 5)
    #
    # # print(results$draws_and_summary[[best_K]]$summarized_results, n = nrow(results$draws_and_summary[[best_K]]$summarized_results))
    # # print(results_model_selection$model_selection_tibble)
    #
    #
    
    # plot_model_selection_inference <- list()
    # 
    # K = Kmax
    # for (i in 1:K){
    #   plot_model_selection_inference[[i]] <- tickTack::plot_timing_h(results, i) + ggplot2::ggtitle(paste0("K = ", i))
    # }
    # plot_model_selection_inference <- gridExtra::grid.arrange(grobs = plot_model_selection_inference, nrow = K)
    # ggsave(paste0(new_dir,"/plot_timing_all_K_h.png"),plot = plot_model_selection_inference, width = 25, height = 30)
    if (return_res == TRUE) {
      return(tbl)
    }
    
  }, error = function(e) {
    message(sprintf("Error processing %s: %s", s, e$message))
    NULL
  })
  
}




safe_run <- function(expr, name) {
  tryCatch(
    expr,
    error = function(e) {
      msg <- paste(Sys.time(), "-", name, "failed with error:", e$message, "\n")
      writeLines(msg, error_log)
      return(NULL)
    }
  )
}


generate_multiplicities = function(k, tau, N_mutations, m=1) {
  if (k == "2:0") {
    # 2:0
    n1 = 2 * m * (1 - tau)
    n2 = m * tau
  } else if (k == "2:1") {
    # 2:1
    n1 = m * (3 - 2*tau)
    n2 = m * tau
  } else if (k == "2:2") {
    # 2:2
    n1 = 4 * m * (1 - tau)
    n2 = 2 * m * tau
  } else {
    stop("k not recognized")
  }
  muts <- c(rmultinom(n = 1, size = N_mutations, prob = c(n1, n2)))
  list(n1 = muts[1], n2 = muts[2])
}

k = "2:0"
tau = .25
N_mutations = 20

n1 = 2 * m * (1 - tau)
n2 = m * tau



muts <- rmultinom(n = 1000, size = N_mutations, prob = c(n1, n2))

20 * (n1 / (n1 + n2))
20 * (n2 / (n1 + n2))

rpois(1000, 20 * (n1 / (n1 + n2)))

N2_sampled = c()
N1_sampled = c()
for (i in 1:1000) {
  N1 = 0
  N2 = 0
  while (N1 + N2 != N_mutations) {
    N1 = rpois(1, 20 * (n1 / (n1 + n2)))
    N2 = rpois(1, 20 * (n2 / (n1 + n2)))
  }
  N2_sampled = c(N2_sampled, N2)
  N1_sampled = c(N1_sampled, N1)
}

hist(N2_sampled)
hist(N1_sampled)

muts[1,] %>% hist()





simulate_dataset = function(N_events, N_clocks, mutation_density, pi, coverage, sigma_tau = .01, min_dist = .1) {
  
  repeat {
    taus_sampled <- runif(N_clocks, 0.01, 0.99)
    if (all(diff(sort(taus_sampled)) >= min_dist)) break
  }
  
  clock_indices <- 1:N_clocks
  if (N_events > N_clocks) {
    clock_indices <- c(clock_indices, sample(1:N_clocks, N_events - N_clocks, replace = TRUE))
  }
  clock_indices <- sample(clock_indices)  

  taus_clust <- taus_sampled[clock_indices]
  taus_events <- sapply(taus_clust, function(t) {
    rnorm(1, mean = t, sd = sigma_tau)
  })
  taus_events <- pmin(pmax(taus_events, 0), 1)       
  taus_clust    # tau assignment per event (original tau, no noise)
  taus <- taus_events   # tau per event with noise
  
  karyotypes = sample(c("2:0", "2:1", "2:2"), N_events, replace = TRUE)
  
  sim_data = lapply(1:N_events, function(idx) {
    segment_length <- 5e7   
    N_mutations <- max(10,rpois(1, mutation_density * segment_length))
    
    m = 1
    tau = taus[idx]
    k = karyotypes[idx]
    
    mults = generate_multiplicities(k = k, tau = tau, N_mutations = N_mutations, m = 1)
    
    start_off = 10000
    
    Major = as.numeric(unlist(strsplit(k, ":"))[1])
    minor = as.numeric(unlist(strsplit(k, ":"))[2])
    
    if(idx <=30) {
      cn = dplyr::tibble(
        chr = 1,
        startpos = 1 + (idx - 1) * segment_length,
        endpos   = idx * segment_length,
        nMaj1_A  = Major,
        nMin1_A  = minor
      ) %>%
        dplyr::mutate(startpos = startpos + start_off, endpos=endpos+start_off)
      
      
      mut_positions <- sort(sample(
        cn$startpos:cn$endpos,
        size = N_mutations,
        replace = FALSE
      ))
      
      mult = dplyr::tibble(
        chr = 1,
        end = mut_positions,
        no.chrs.bearing.mut = c(rep(1, mults$n1), rep(2, mults$n2))
      )
      
      muts = dplyr::tibble(
        chr   = 1,
        start = mut_positions,
        end   = mut_positions,
        ref   = "G",
        alt   = "A"
      )
    } else if(idx <= 60 ) {
      cn = dplyr::tibble(
        chr = 2,
        startpos = 1 + (idx - 30 - 1) * segment_length,
        endpos   = (idx - 30) * segment_length,
        nMaj1_A  = Major,
        nMin1_A  = minor
      ) %>%
        dplyr::mutate(startpos = startpos + start_off, endpos=endpos+start_off)
      
      
      mut_positions <- sort(sample(
        cn$startpos:cn$endpos,
        size = N_mutations,
        replace = FALSE
      ))
      
      mult = dplyr::tibble(
        chr = 2,
        end = mut_positions,
        no.chrs.bearing.mut = c(rep(1, mults$n1), rep(2, mults$n2))
      )
      
      muts = dplyr::tibble(
        chr   = 2,
        start = mut_positions,
        end   = mut_positions,
        ref   = "G",
        alt   = "A"
      )
    } else {
      cn = dplyr::tibble(
        chr = 3,
        startpos = 1 + (idx - 60 - 1) * segment_length,
        endpos   = (idx - 60)* segment_length,
        nMaj1_A  = Major,
        nMin1_A  = minor
      ) %>%
        dplyr::mutate(startpos = startpos + start_off, endpos=endpos+start_off)
      
      
      mut_positions <- sort(sample(
        cn$startpos:cn$endpos,
        size = N_mutations,
        replace = FALSE
      ))
      
      mult = dplyr::tibble(
        chr = 3,
        end = mut_positions,
        no.chrs.bearing.mut = c(rep(1, mults$n1), rep(2, mults$n2))
      )
      
      muts = dplyr::tibble(
        chr   = 3,
        start = mut_positions,
        end   = mut_positions,
        ref   = "G",
        alt   = "A"
      )
    }
    
    # cn = dplyr::tibble(
    #   chr = 1,
    #   startpos = 1 + (idx - 1) * segment_length,
    #   endpos   = idx * segment_length,
    #   nMaj1_A  = Major,
    #   nMin1_A  = minor
    # ) %>%
    #   dplyr::mutate(startpos = startpos + start_off, endpos=endpos+start_off)
    # 
    # 
    # mut_positions <- sort(sample(
    #   cn$startpos:cn$endpos,
    #   size = N_mutations,
    #   replace = FALSE
    # ))
    # 
    # mult = dplyr::tibble(
    #   chr = 1,
    #   end = mut_positions,
    #   no.chrs.bearing.mut = c(rep(1, mults$n1), rep(2, mults$n2))
    # )
    # 
    # muts = dplyr::tibble(
    #   chr   = 1,
    #   start = mut_positions,
    #   end   = mut_positions,
    #   ref   = "G",
    #   alt   = "A"
    # )
    
    # Add DP and
    peaks <- get_clonal_peaks(k, pi)
    
    p_vec <- peaks[mult$no.chrs.bearing.mut]

    muts$DP <- rpois(n = N_mutations, lambda = coverage)
    muts$DP[muts$DP == 0] <- 1
    muts$NV <- rbinom(n = N_mutations, size = muts$DP, prob = p_vec)
    muts$VAF <- muts$NV / muts$DP
    
    list(cn = cn, mult=mult, muts=muts)
  })
  
  sim_data
  
  sim_cn <- lapply(1:N_events, function(i) { sim_data[[i]]$cn }) %>% do.call("rbind", .)
  sim_mult <- lapply(1:N_events, function(i) { sim_data[[i]]$mult }) %>% do.call("rbind", .)
  sim_muts <- lapply(1:N_events, function(i) { sim_data[[i]]$muts }) %>% do.call("rbind", .)
  
  list(cn=sim_cn, mult=sim_mult, muts=sim_muts, true_taus = taus, taus_clust=taus_clust)
}




fit_AmpTimeR = function(sim, mult_path) {
  N_events = nrow(sim$cn)
  lapply(1:N_events, function(idx) {
    print(idx)
    cn = sim$cn[idx, ]
    mult = read.delim(mult_path, sep = "\t", header = TRUE)
    muts = sim$muts %>% dplyr::select(chr, start, end, ref, alt)
    
    segment_time <- AmplificationTimeR::time_amplification(
      cn_data = cn %>% as.data.frame(),
      multiplicity_data = mult %>% 
        dplyr::filter(chr == cn$chr, end >= cn$startpos, end <= cn$endpos) %>% 
        as.data.frame() %>% 
        dplyr::select(chr, end, no.chrs.bearing.mut),
      mutation_data = muts %>% dplyr::filter(chr == cn$chr, end >= cn$startpos, end <= cn$endpos) %>% as.data.frame(),
      muts_type = "SBS1 and SBS5",
      sample_id = "test_sample",
      amplification_chrom = cn$chr,
      amplification_start = cn$startpos,
      amplification_stop = cn$endpos,
      is_WGD = TRUE,
      genome = "hg19"
    )
    
    dplyr::tibble(
      segment_idx = idx, 
      tau = segment_time$t_1_mean_bootstrap, 
      tau_low = segment_time$t_1_lower_ci, 
      tau_high = segment_time$t_1_upper_ci, 
      model = "AmpTimeR"
    )
  }) %>% do.call("bind_rows", .)
}




fit_tickTack_single = function(sim, pi, min_mutations) {
  N_events = nrow(sim$cn)
  pcawg_example$cna
  
  cn <- sim$cn %>% 
    dplyr::rename(Major=nMaj1_A, minor=nMin1_A, from=startpos, to=endpos)
  
  muts <- sim$muts %>% 
    dplyr::rename(from=start, to=end)
  
  tickTack_single_res <- tickTack::fit(cn, mutations = muts, purity = pi, alpha = .05, min_mutations_number = min_mutations, beta_binomial = F)
  tickTack_single_res
}


fit_tickTack_h = function(sim, pi, min_mutations, tolerance, INIT, max_attempts) {
  N_events = nrow(sim$cn)
  
  cn <- sim$cn %>% 
    dplyr::rename(Major=nMaj1_A, minor=nMin1_A, from=startpos, to=endpos) %>% 
    dplyr::mutate(CCF = 1)
  
  muts <- sim$muts %>% 
    dplyr::rename(from=start, to=end) %>% 
    dplyr::mutate(CCF = 1)
  
  x = list(
    cna = cn, 
    mutations = muts,
    metadata = dplyr::tibble(sample = "sample", purity=pi)
  )
  
  x <- tickTack::fit_h(x, max_attempts=max_attempts, INIT=INIT, tolerance = tolerance)
  
  results_simulated <- x$results_timing
  results_model_selection <- tickTack::model_selection_h(results_simulated)
  best_K <- results_model_selection$best_K
  results_tickTack <- list( results = results_simulated, results_model_selection = results_model_selection)
  
  results_tickTack
}

fit_MutTimeR <- function(sim, pi) {
  
  cn <- sim$cn %>% 
    dplyr::rename(Major=nMaj1_A, minor=nMin1_A, from=startpos, to=endpos) %>% 
    dplyr::mutate(CCF = 1)
  
  muts <- sim$muts %>% 
    dplyr::rename(from=start, to=end) %>% 
    dplyr::mutate(CCF = 1)
  
  x = list(
    cna = cn, 
    mutations = muts,
    metadata = dplyr::tibble(sample = "sample", purity=pi)
  )
  
  
  vcf <- convert_to_vcf(x$mutations, x$cna, x$metadata)
  bb <- convert_to_granges(x$cna, x$mutations, x$metadata)
  # clusters <- data.frame(cluster=1:1, proportion=c(pi), n_ssms=c(100))
  
  # mt <- MutationTimeR::mutationTime(vcf, bb, clusters=clusters)
  mt <- MutationTimeR::mutationTime(vcf, bb, n.boot = 10)
  mcols(bb) <- cbind(mcols(bb),mt$T)
  
  res_MutTime <- list(vcf = vcf, cn_timed = bb)  
  res_MutTime
}



convert_to_vcf <- function(mutations, cna, metadata) {
  # Validate input
  if(missing(mutations) || missing(cna) || missing(metadata)) {
    stop("mutations, cna, and metadata are required")
  }
  
  # Create VRanges object
  v <- VRanges(
    seqnames = unlist(lapply(mutations$chr, function(c) {str_replace_all(c, "chr", "")})),
    ranges = IRanges(start = mutations$from, end = mutations$to),
    ref = mutations$ref,
    # alt = mutations$alt,
    totalDepth = mutations$DP,
    altDepth = mutations$NV,
    QUAL = rep(NA_real_, nrow(mutations)),
    FILTER = rep("PASS", nrow(mutations))
  )
  
  # Set sample name
  sampleNames(v) <- metadata$sample
  
  # Convert to VCF
  vcf <- as(v, "VCF")
  
  # Create custom geno header
  geno_header <- DataFrame(
    Number = c("2", "1", "1"),
    Type = c("Integer", "Integer", "String"),
    Description = c(
      "Allelic depths (number of reads in each observed allele)",
      "Total read depth",
      "Variant filters"
    ),
    row.names = c("AD", "DP", "FT")
  )
  
  # Add custom info headers
  info(header(vcf)) <- rbind(
    info(header(vcf)),
    DataFrame(
      Number = 1,
      Type = rep("Integer", 2),
      Description = c("Tumour ref count", "Tumour alt count"),
      row.names = c("t_ref_count", "t_alt_count")
    )
  )
  
  # Set geno header
  geno_header <- DataFrame(
    Number = c("2", "1", "1"),
    Type = c("Integer", "Integer", "String"),
    Description = c(
      "Allelic depths of reference and alternate alleles",
      "Total read depth",
      "Variant filters"
    ),
    row.names = c("AD", "DP", "FT")
  )
  geno(vcf)
  geno(header(vcf)) <- geno_header
  
  # Prepare genotype data
  geno(vcf)$AD <- as.matrix(rep(NA, nrow(mutations)))
  geno(vcf)$DP <- as.matrix(mutations$DP)
  geno(vcf)$FT <- as.matrix(rep(NA, nrow(mutations)))
  
  # Add info fields
  info(vcf)$t_alt_count <- altDepth(v)
  info(vcf)$t_ref_count <- totalDepth(v) - altDepth(v)
  
  return(vcf)
}

convert_to_granges <- function(cna, mutations, metadata) {
  # Validate input
  if(missing(cna) || missing(mutations) || missing(metadata)) {
    stop("All three inputs (cna, mutations, metadata) are required")
  }
  
  # Create GRanges
  cna_granges <- GRanges(
    seqnames = gsub("^chr", "", cna$chr),
    ranges = IRanges(start = cna$from, end = cna$to),
    strand = "*",
    major_cn = cna$Major,
    minor_cn = cna$minor,
    clonal_frequency = metadata$purity
  )
  
  return(cna_granges)
}


plot_heatmap <- function(matrix_data, 
                         x_label = "Columns", 
                         y_label = "Rows", 
                         fill_label = "Value", 
                         low_color = "blue", 
                         high_color = "red", 
                         title = "Heatmap") {
  # Check if input is a matrix
  if (!is.matrix(matrix_data)) {
    stop("Input must be a matrix.")
  }
  
  # Convert matrix to a long-format data frame
  df <- melt(matrix_data)
  colnames(df) <- c("Row", "Column", "Value")
  
  # Create the heatmap
  ggplot(df, aes(x = Column, y = Row, fill = Value)) +
    geom_tile() +
    scale_fill_gradient(low = low_color, high = high_color, name = fill_label) +
    labs(x = x_label, y = y_label, title = title) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 10),
      plot.title = element_text(hjust = 0.5)
    )
}



get_multiplicities <- function(sim, purity, mult_path, subdir){
  
  muts <- sim$muts
  df_counts <- muts %>%
    mutate(WT.count = DP - NV,
           mut.count = NV)
  
  info_counts <- GRanges(
    seqnames = paste0("chr", df_counts$chr),
    ranges = IRanges(start = df_counts$start, end = df_counts$end),
    strand = "*",
    WT.count = df_counts$WT.count,
    mut.count = df_counts$mut.count
  )
  
  cn <- sim$cn%>% 
    mutate(CHR = chr, START=startpos, frac_1A=1)%>%
    dplyr::select("CHR","START","nMaj1_A","nMin1_A")
  
  
  info_counts <- GRanges(
    seqnames = paste0(df_counts$chr),
    ranges = IRanges(start = df_counts$start, end = df_counts$end),
    strand = "*",
    WT.count = df_counts$WT.count,
    mut.count = df_counts$mut.count
  )
  
  cn <- sim$cn %>% 
    #dplyr::select("CHR","START","nMaj1_A","nMin1_A") %>% 
    dplyr::mutate(ntot = nMaj1_A + nMin1_A, frac1_A = 1, nMaj2_A=NA, nMin2_A=NA, frac2_A=NA)
  
  write.table(cn, paste0(subdir,"/sample_1_subclones.txt"), col.names = TRUE, quote = FALSE, row.names = FALSE, sep = "\t")
  
  dpclust3p:::GetDirichletProcessInfo(mult_path,
                                      purity,
                                      info_counts, 
                                      paste0(subdir,"/sample_1_subclones.txt"),
                                      T,
                                      SNP.phase.file = "NA",
                                      mut.phase.file = NULL)
  # rm("./my_pseudo_dpclust_files")
}


compute_metrics <- function(posterior_draws) {
  # Input: 
  # posterior_draws: A list where each element contains posterior draws for a sample
  
  # Output:
  # A matrix with Wasserstein distance, KL divergence, and overlap-based similarity
  
  n <- length(posterior_draws)
  results <- matrix(NA, nrow = n, ncol = n, 
                    dimnames = list(paste0("Sample", 1:n), paste0("Sample", 1:n)))
  
  wasserstein <- results
  kl_divergence <- results
  overlap_similarity <- results
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      # Extract draws for sample i and j
      draws_i <- posterior_draws[[i]]
      draws_j <- posterior_draws[[j]]
      
      # Compute Wasserstein distance
      wasserstein[i, j] <- wasserstein1d(draws_i, draws_j)
      wasserstein[j, i] <- wasserstein[i, j]
      
      # Compute KL divergence using density estimates
      density_i <- density(draws_i)
      density_j <- density(draws_j)
      
      # Interpolate densities to match the support
      common_x <- sort(unique(c(density_i$x, density_j$x)))
      p <- approx(density_i$x, density_i$y, xout = common_x, rule = 2)$y
      q <- approx(density_j$x, density_j$y, xout = common_x, rule = 2)$y
      
      # Normalize to avoid NaN issues
      p <- p / sum(p)
      q <- q / sum(q)
      
      kl_divergence[i, j] <- sum(ifelse(p > 0 & q > 0, p * log(p / q), 0))
      kl_divergence[j, i] <- sum(ifelse(q > 0 & p > 0, q * log(q / p), 0))
      
      # Compute overlap-based similarity
      overlap_similarity[i, j] <- sum(pmin(p, q))
      overlap_similarity[j, i] <- overlap_similarity[i, j]
    }
  }
  
  list(
    Wasserstein = wasserstein,
    KL_Divergence = kl_divergence,
    Overlap_Similarity = overlap_similarity
  )
}

compute_WD <- function(posterior_draws) {
  # Input: 
  # posterior_draws: A list where each element contains posterior draws for a sample
  
  # Output:
  # with Wasserstein distance
  
  n <- length(posterior_draws)
  results <- matrix(NA, nrow = n, ncol = n, 
                    dimnames = list(paste0("Sample", 1:n), paste0("Sample", 1:n)))
  
  wasserstein <- results
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      # Extract draws for sample i and j
      draws_i <- posterior_draws[[i]]
      draws_j <- posterior_draws[[j]]
      
      # Compute Wasserstein distance
      wasserstein[i, j] <- wasserstein1d(draws_i, draws_j)
      wasserstein[j, i] <- wasserstein[i, j]
    }
  }
  
  wasserstein
}


compute_overlap_matrix <- function(intervals) {
  N <- nrow(intervals)
  overlap_matrix <- matrix(0, nrow = N, ncol = N)
  
  for (i in 1:N) {
    for (j in 1:N) {
      if (i == j) {
        overlap_matrix[i, j] <- 1  # An interval fully overlaps with itself
      } else {
        a1 <- intervals[i, 1]
        b1 <- intervals[i, 2]
        a2 <- intervals[j, 1]
        b2 <- intervals[j, 2]
        
        intersection_start <- max(a1, a2)
        intersection_end <- min(b1, b2)
        intersection_length <- max(0, intersection_end - intersection_start)
        
        length1 <- b1 - a1
        length2 <- b2 - a2
        
        
        if (length1 > 0 && length2 > 0) {
          overlap_matrix[i, j] <- intersection_length / min(length1, length2)
        } else {
          overlap_matrix[i, j] <- 0  # Avoid division by zero
        }
      }
    }
  }
  
  return(overlap_matrix)
}




#PROBLEM: A matrix is “computationally singular” when: some clusters have zero variance (a cluster has identical points) 
# or data contains duplicate rows, or k-means converged to a solution where a cluster is assigned a single point.

# The singularity of the within-cluster covariance matrix indicates that the clustering solution produces degenerate clusters 
# (e.g., clusters with very low variance or singleton clusters). Because NbClust relies on matrix inversion for several validity 
# indices, these degenerate solutions make some indices mathematically undefined. Consequently, only indices that remain 
# well-defined under such conditions (e.g., Silhouette, Dunn, Davies–Bouldin) were used for evaluation.


# clustering methods : "kmeans", "pam", "hcut" 
# score methods :"silhouette" "gap_stat" , gap statistics to include a score for n clusters = 1 

library(NbClust)

clustering_methods <- function(intervals, data_type = NULL) {
  
  #--------------- SPECIAL HANDLING FOR AMPLIFICATION TIMER ----------------#
  if (!is.null(data_type) && data_type == "amplificationtimer") {
    
    unique_rows <- unique(intervals)
    
    # Case 1: all zeros
    if (all(intervals == 0)) {
      return(list(
        best_K_kmeans = 1,
        clusters_kmeans = rep(1, nrow(intervals)),
        best_K_pam = 1,
        clusters_pam = rep(1, nrow(intervals)),
        best_K_hcut = 1,
        clusters_hcut = rep(1, nrow(intervals)),
        
        best_gap_kmeans = 1,
        clusters_gap_kmenas = rep(1, nrow(intervals)),
        best_gap_pam = 1,
        clusters_gap_pam = rep(1, nrow(intervals)),
        best_gap_hcut = 1,
        clusters_gap_hcut = rep(1, nrow(intervals))
      ))
    }
    
    # Case 2: Some distinct rows, but no scaling
    if (any(apply(intervals, 2, sd) < 1e-10)) {
      # message("Scaling skipped: near-zero variance detected")
      data_prep <- intervals
    } else {
      data_prep <- scale(intervals)
    }
    unique_rows <- unique(data_prep)
    max_clusters <- min(nrow(unique_rows), nrow(intervals) - 1)
    
  } else {
    
    if (any(apply(intervals, 2, sd) < 1e-10)) {
      # message("⚠️ Scaling skipped: near-zero variance detected")
      data_prep <- intervals
    } else {
      data_prep <- scale(intervals)
    }
    unique_rows <- unique(data_prep)
    max_clusters <- min(nrow(unique_rows), nrow(intervals) - 1)
    
  }
  
  if (NCOL(data_prep)==1){
    best_K_kmeans = 1
    clusters_kmeans = as.numeric(as.factor(data_prep[1]))
  } else if (length(data_prep[,1]) < 2 ) {
    best_K_kmeans = 2 
    clusters_kmeans = as.numeric(as.factor(data_prep[,1]))
  } else {
    kmeans_selection <- factoextra::fviz_nbclust(data_prep, kmeans, method = "silhouette", k.max = max_clusters)
    best_K_kmeans <- which( kmeans_selection$data$y == max(kmeans_selection$data$y))
    km = kmeans(x = data_prep, centers = best_K_kmeans)
    clusters_kmeans = km$cluster
  }
  gap_stat_kmeans <- clusGap(data_prep, FUN = kmeans, nstart = 25,
                             K.max = max_clusters, B = 10)
  best_gap_kmeans <- maxSE(f = gap_stat_kmeans$Tab[,"gap"],
                           SE.f = gap_stat_kmeans$Tab[,"SE.sim"],
                           method = "firstSEmax")
  gap_kmeans <- kmeans(data_prep, centers = best_gap_kmeans, nstart = 25)
  clusters_gap_kmenas = gap_kmeans$cluster
  
  
  if (NCOL(data_prep)==1){
    best_K_pam = 1
    clusters_pam = as.numeric(as.factor(data_prep[1]))
  } else if (length(data_prep[,1]) < 2 ) {
    best_K_pam = 2 
    clusters_pam = as.numeric(as.factor(data_prep[,1]))
  } else {
    res_pam <- factoextra::fviz_nbclust(data_prep, cluster::pam, method = "silhouette", k.max = max_clusters)
    best_K_pam <- which( res_pam$data$y == max(res_pam$data$y))
    pam = cluster::pam(x = data_prep, k = best_K_pam)
    clusters_pam = pam$cluster
  }
  gap_stat_pam <- clusGap(data_prep, FUN = cluster::pam, nstart = 25,
                          K.max = max_clusters, B = 10)
  
  best_gap_pam <- maxSE(f = gap_stat_pam$Tab[,"gap"],
                        SE.f = gap_stat_pam$Tab[,"SE.sim"],
                        method = "firstSEmax")
  gap_pam <- cluster::pam(data_prep, k = best_gap_pam)
  clusters_gap_pam = gap_pam$cluster
  best_gap_pam
  
  if (NCOL(data_prep)==1){
    best_K_hcut = 1
    clusters_pam = as.numeric(as.factor(data_prep[1]))
  } else if (length(data_prep[,1]) < 2 ) {
    best_K_hcut = 2 
    clusters_hcut = as.numeric(as.factor(data_prep[,1]))
  } else {
    res_hcut <- factoextra::fviz_nbclust(data_prep, factoextra::hcut, method = "silhouette", k.max = max_clusters)
    best_K_hcut <- which( res_hcut$data$y == max(res_hcut$data$y))
    hcut = factoextra::hcut(x = data_prep, k = best_K_hcut)
    clusters_hcut = hcut$cluster
  }
  
  gap_stat_hcut <- clusGap(data_prep, FUN = factoextra::hcut, nstart = 25,
                           K.max = max_clusters, B = 10)
  best_gap_hcut <- maxSE(f = gap_stat_hcut$Tab[,"gap"],
                         SE.f = gap_stat_hcut$Tab[,"SE.sim"],
                         method = "firstSEmax")
  gap_hcut <- factoextra::hcut(data_prep, k = best_gap_hcut, nstart = 25)
  clusters_gap_hcut = gap_hcut$cluster
  
  
  return(list(best_K_kmeans = best_K_kmeans,
              clusters_kmeans = clusters_kmeans,
              best_K_pam = best_K_pam,
              clusters_pam = clusters_pam,
              best_K_hcut = best_K_hcut,
              clusters_hcut = clusters_hcut,
              
              best_gap_kmeans = best_gap_kmeans,
              clusters_gap_kmenas = clusters_gap_kmenas,
              best_gap_pam = best_gap_pam, 
              clusters_gap_pam = clusters_gap_pam,
              best_gap_hcut = best_gap_hcut,
              clusters_gap_hcut = clusters_gap_hcut)
  )
  
}




# Cluster AmpTimeR
clusterAmpTimeR = function(fp) {
  print("clusterAmp_timer")
  res = readRDS(paste0(fp, '/res_AmpTimeR.rds'))
  intervals = cbind(res$tau_low, res$tau_high)
  na_indices <- which(apply(intervals, 1, function(x) any(is.na(x))))
  if (length(na_indices)!=0){
    intervals <- intervals[-na_indices, , drop = FALSE]
  }
  similarity_matrix = compute_overlap_matrix(intervals)
  
  result_cluster <- clustering_methods(intervals,"amplificationtimer")
  
  return(  list(na_indices=na_indices,
                best_k_kmeans = result_cluster$best_K_kmeans,
                clusters_kmeans = result_cluster$clusters_kmeans,
                best_K_pam = result_cluster$best_K_pam,
                clusters_pam = result_cluster$clusters_pam,
                best_K_hcut = result_cluster$best_K_hcut,
                clusters_hcut = result_cluster$clusters_hcut,
                
                best_gap_kmeans = result_cluster$best_gap_kmeans,
                clusters_gap_kmenas = result_cluster$clusters_gap_kmenas,
                best_gap_pam = result_cluster$best_gap_pam,
                clusters_gap_pam = result_cluster$clusters_gap_pam,
                best_gap_hcut = result_cluster$best_gap_hcut,
                clusters_gap_hcut = result_cluster$clusters_gap_hcut))
}



# Cluster MutTimeR
clusterMutTimeR = function(fp) {
  print("clusterMutTimer")
  res = readRDS(paste0(fp, '/res_MutTime.rds'))
  intervals = cbind(res$cn_timed$time.lo, res$cn_timed$time.up)
  similarity_matrix = compute_overlap_matrix(intervals)
  
  result_cluster <- clustering_methods(intervals)
  
  return(  list(best_k_kmeans = result_cluster$best_K_kmeans,
                clusters_kmeans = result_cluster$clusters_kmeans,
                best_K_pam = result_cluster$best_K_pam,
                clusters_pam = result_cluster$clusters_pam,
                best_K_hcut = result_cluster$best_K_hcut,
                clusters_hcut = result_cluster$clusters_hcut,
                
                best_gap_kmeans = result_cluster$best_gap_kmeans,
                clusters_gap_kmenas = result_cluster$clusters_gap_kmenas,
                best_gap_pam = result_cluster$best_gap_pam,
                clusters_gap_pam = result_cluster$clusters_gap_pam,
                best_gap_hcut = result_cluster$best_gap_hcut,
                clusters_gap_hcut = result_cluster$clusters_gap_hcut))
}




# Cluster tickTack single
cluster_tickTack_single = function(fp) {
  print("clustertickTack_single")
  res = readRDS(paste0(fp, '/res_tickTack_single.rds'))
  draws = lapply(res$inference_results$segment %>% unique(), function(idx) {
    res$inference_result %>%
      dplyr::filter(segment == idx) %>%
      dplyr::filter(tau <= 1, tau >= 0) %>%
      pull(tau)
  })
  similarity_matrix = compute_WD(draws)
  similarity_matrix[is.na(similarity_matrix)] = 0
  
  mean_mat <- sapply(draws, mean)
  sd_mat   <- sapply(draws, sd)
  q05_mat  <- sapply(draws, quantile, 0.05)
  names(q05_mat) <- NULL
  q95_mat  <- sapply(draws, quantile, 0.95)
  names(q95_mat) <- NULL
  
  feature_matrix <- cbind(mean_mat, sd_mat, q05_mat, q95_mat)
  
  result_cluster <- clustering_methods(feature_matrix)
  
  return(  list(best_k_kmeans = result_cluster$best_K_kmeans,
                clusters_kmeans = result_cluster$clusters_kmeans,
                best_K_pam = result_cluster$best_K_pam,
                clusters_pam = result_cluster$clusters_pam,
                best_K_hcut = result_cluster$best_K_hcut,
                clusters_hcut = result_cluster$clusters_hcut,
                
                best_gap_kmeans = result_cluster$best_gap_kmeans,
                clusters_gap_kmenas = result_cluster$clusters_gap_kmenas,
                best_gap_pam = result_cluster$best_gap_pam,
                clusters_gap_pam = result_cluster$clusters_gap_pam,
                best_gap_hcut = result_cluster$best_gap_hcut,
                clusters_gap_hcut = result_cluster$clusters_gap_hcut))
}

cluster_tickTack_h = function(fp) {
  res = readRDS(paste0(fp, '/res_tickTack_h.rds'))
  res$results_model_selection$best_fit$summarized_results$clock_mean
}


plot_sim <- function(sim){
  plot_data <- data.frame(
    taus = sim$taus_clust,    # base cluster assignment
    taus_noise = sim$true_taus  # true tau values with noise
  )
  
  label_data <- plot_data %>%
    group_by(taus) %>%
    summarise(count = n(), .groups = "drop")
  
  p <- ggplot(plot_data, aes(x = taus)) +
    geom_histogram(
      bins = 80,
      fill = "#FDAE61",
      color = "black",
      alpha = 0.3
    ) +
    geom_rug(aes(x = taus_noise), sides = "t", color = "black", alpha = 0.7) +
    geom_text(
      data = label_data,
      aes(x = taus, y = count, label = round(taus, 2)),
      vjust = -0.5,
      size = 4,
      fontface = "bold"
    ) +
    expand_limits(y = max(label_data$count) * 1.2) +
    xlim(0, 1) +
    theme_light(base_size = 14) +
    labs(x = "Tau Clusters", y = "Count")
  
  return(p)
}

