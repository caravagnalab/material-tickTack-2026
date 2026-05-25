data {
  int S;
  int K;        // K_max: upper bound, set generously (e.g. 10)
  int N;
  array[S] int karyotype;
  array[N] int seg_assignment;
  array[S,2] real peaks;
  array[N] int NV;
  array[N] int DP;
  
  // Sparse Dirichlet concentration:
  // alpha << 1 (e.g. 0.01) encourages sparsity
  // alpha = 1  recovers the original flat prior
  real<lower=0> alpha;
}

parameters {
  simplex[K] pi;                                    // sparse mixture weights
  vector<lower=0.001, upper=0.999>[K] tau;          // one tau per clock
}

transformed parameters {
  array[S,K,2] real<lower=0,upper=1> theta;
  for (s in 1:S) {
    for (k in 1:K) {
      if (karyotype[s] == 1) {
        theta[s,k,1] = (3 - 2*tau[k]) / (3 - tau[k]);
        theta[s,k,2] = tau[k] / (3 - tau[k]);
      } else {
        theta[s,k,1] = (2 - 2*tau[k]) / (2 - tau[k]);
        theta[s,k,2] = tau[k] / (2 - tau[k]);
      }
    }
  }
}

model {
  // --- Sparse Dirichlet prior on mixture weights ----------------------
  // alpha << 1 creates sparsity: most pi[k] shrink to ~0
  // The effective number of components is inferred from the data
  pi ~ dirichlet(rep_vector(alpha, K));

  // --- Prior on tau ---------------------------------------------------
  tau ~ beta(0.5, 0.5);

  // --- Likelihood: marginalise over cluster assignments ---------------
  for (s in 1:S) {
    vector[K] seg_lp;
    for (k in 1:K) {
      seg_lp[k] = log(pi[k]);
      for (i in 1:N) {
        if (seg_assignment[i] == s) {
          seg_lp[k] += log_sum_exp(
            log(theta[s,k,1]) + binomial_lpmf(NV[i] | DP[i], peaks[s,1]),
            log(theta[s,k,2]) + binomial_lpmf(NV[i] | DP[i], peaks[s,2])
          );
        }
      }
    }
    target += log_sum_exp(seg_lp);
  }
}

generated quantities {
  array[S] vector[K] seg_probs;
  array[S] int seg_assignment_hard;
  vector[S] log_lik;

  // --- Effective number of components ---------------------------------
  // Count components with non-negligible weight (pi[k] > threshold)
  real pi_threshold = 1.0 / K;   // component gets less than equal share
  int  K_eff = 0;
  for (k in 1:K)
    if (pi[k] > pi_threshold) K_eff += 1;

  // --- Soft and hard assignments --------------------------------------
  for (s in 1:S) {
    vector[K] seg_lp;
    for (k in 1:K) {
      seg_lp[k] = log(pi[k]);
      for (i in 1:N) {
        if (seg_assignment[i] == s) {
          seg_lp[k] += log_sum_exp(
            log(theta[s,k,1]) + binomial_lpmf(NV[i] | DP[i], peaks[s,1]),
            log(theta[s,k,2]) + binomial_lpmf(NV[i] | DP[i], peaks[s,2])
          );
        }
      }
    }

    log_lik[s]             = log_sum_exp(seg_lp);
    seg_probs[s]           = softmax(seg_lp);
    seg_assignment_hard[s] = 1;
    for (k in 2:K) {
      if (seg_probs[s][k] > seg_probs[s][seg_assignment_hard[s]])
        seg_assignment_hard[s] = k;
    }
  }
}
