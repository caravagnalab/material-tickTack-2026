data {
  int<lower=1> S;
  int<lower=1> K;
  int<lower=1> N;
  array[S] int<lower=0,upper=1> karyotype;
  array[N] int<lower=1,upper=S> seg_assignment;
  array[S,2] real<lower=0,upper=1> peaks;
  array[N] int<lower=0> NV;
  array[N] int<lower=1> DP;
  
  // Hyperparameters passed as data for flexibility
  real<lower=0> alpha_dirichlet;   // global pi concentration (try 1.0 or 0.1)
  real<lower=0> gamma_local;       // local w concentration  (try 0.5)
}

parameters {
  simplex[K] pi;

  // FIX 1: ordered tau breaks label switching definitively
  // The model now has a canonical ordering: tau[1] < tau[2] < ... < tau[K]
  ordered[K] tau_raw;              
}

transformed parameters {
  // Rescale ordered tau to (0.001, 0.999) after ordering
  vector<lower=0.001, upper=0.999>[K] tau;
  tau = inv_logit(tau_raw) * 0.998 + 0.001;

  array[S,K,2] real<lower=0,upper=1> theta;
  for (s in 1:S) {
    for (k in 1:K) {
      if (karyotype[s] == 1) {
        theta[s,k,1] = (3 - 2*tau[k]) / (3 - tau[k]);
        theta[s,k,2] =      tau[k]     / (3 - tau[k]);
      } else {
        theta[s,k,1] = (2 - 2*tau[k]) / (2 - tau[k]);
        theta[s,k,2] =      tau[k]     / (2 - tau[k]);
      }
    }
  }

  // FIX 2: pre-compute per-segment log-likelihoods for each cluster
  // Shape: seg_lp_data[s,k] = sum of binomial lpmf for all mutations in s under cluster k
  // This separates data likelihood from mixing weights (cleaner + faster)
  array[S,K] real seg_data_lp;
  for (s in 1:S) {
    for (k in 1:K) {
      seg_data_lp[s,k] = 0.0;
      for (i in 1:N) {
        if (seg_assignment[i] == s) {
          seg_data_lp[s,k] += log_sum_exp(
            log(theta[s,k,1]) + binomial_lpmf(NV[i] | DP[i], peaks[s,1]),
            log(theta[s,k,2]) + binomial_lpmf(NV[i] | DP[i], peaks[s,2])
          );
        }
      }
    }
  }
}

model {
  // --- Priors ---
  pi  ~ dirichlet(rep_vector(alpha_dirichlet, K));
  
  // Prior on tau_raw (unconstrained); beta(1,1) on tau ≈ normal(0,2) on tau_raw
  tau_raw ~ normal(0, 2);

  // --- Likelihood ---
  // FIX 3: use pi as a prior weight, but let per-segment data dominate
  for (s in 1:S) {
    vector[K] seg_lp;
    for (k in 1:K) {
      seg_lp[k] = log(pi[k]) + seg_data_lp[s,k];
    }
    target += log_sum_exp(seg_lp);
  }
}

generated quantities {
  array[S] simplex[K] seg_probs;
  array[S] int seg_assignment_hard;
  vector[S] log_lik;

  for (s in 1:S) {
    vector[K] seg_lp;
    for (k in 1:K) {
      seg_lp[k] = log(pi[k]) + seg_data_lp[s,k];
    }
    log_lik[s]              = log_sum_exp(seg_lp);
    seg_probs[s]            = softmax(seg_lp);
    seg_assignment_hard[s]  = sort_indices_desc(seg_probs[s])[1];
  }
}
