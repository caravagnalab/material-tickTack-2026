
data {
  int S; int K; int N;
  array[S] int karyotype;
  array[N] int seg_assignment;
  array[S,2] real peaks;
  array[N] int NV;
  array[N] int DP;
}

parameters {
  simplex[K] pi;                          // global cluster weights (ONE, not per-segment)
  vector<lower=0.001, upper=0.999>[K] tau; // one tau per clock
}

transformed parameters {
  array[S,K,2] real<lower=0,upper=1> theta;
  for (s in 1:S) for (k in 1:K) {
    if (karyotype[s] == 1) {
      theta[s,k,1] = (3 - 2*tau[k]) / (3 - tau[k]);
      theta[s,k,2] = tau[k] / (3 - tau[k]);
    } else {
      theta[s,k,1] = (2 - 2*tau[k]) / (2 - tau[k]);
      theta[s,k,2] = tau[k] / (2 - tau[k]);
    }
  }
}

model {
  pi ~ dirichlet(rep_vector(1.0, K));
  // pi ~ dirichlet(rep_vector(0.01, K));
  tau ~ beta(0.5, 0.5);

  // Segment-level soft assignment: marginalise over K
  for (s in 1:S) {
    vector[K] seg_lp;
    // collect all mutations in segment s
    for (k in 1:K) {
      seg_lp[k] = log(pi[k]);
      for (i in 1:N) {
        if (seg_assignment[i] == s) {
          vector[2] allele_lp;
          for (j in 1:2)
            allele_lp[j] = log(0.5) +   // or learned allele weight
                            binomial_lpmf(NV[i] | DP[i], peaks[s,j]);
          // theta weighting on allele states
          seg_lp[k] += log_sum_exp(
            log(theta[s,k,1]) + binomial_lpmf(NV[i] | DP[i], peaks[s,1]),
            log(theta[s,k,2]) + binomial_lpmf(NV[i] | DP[i], peaks[s,2])
          );
        }
      }
    }
    target += log_sum_exp(seg_lp);  // marginalise z_s
  }
}

generated quantities {
  // Soft assignment: posterior probability that segment s belongs to clock k
  array[S] vector[K] seg_probs;
  array[S] int seg_assignment_hard;  // hard assignment: argmax
  vector[S] log_lik;

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
    
    log_lik[s] = log_sum_exp(seg_lp);  // marginal log-lik for segment s

    // Normalise to get probabilities: this is the E-step
    seg_probs[s] = softmax(seg_lp);

    // Hard assignment: which clock has highest posterior probability
    seg_assignment_hard[s] = 1;
    for (k in 2:K) {
      if (seg_probs[s][k] > seg_probs[s][seg_assignment_hard[s]])
        seg_assignment_hard[s] = k;
    }
  }
}
