
data {
  int S; int N;
  array[S] int karyotype;
  array[N] int seg_assignment;
  array[S,2] real peaks;
  array[N] int NV;
  array[N] int DP;
}

parameters {
  array[S] simplex[2] w;
}

transformed parameters {
  vector<lower=0.001, upper=0.999>[S] tau;
  for (s in 1:S) {
    real w1 = w[s][1];
    real numerator = (karyotype[s] == 1) ? 3.0 * (1.0 - w1) : 2.0 * (1.0 - w1);
    tau[s] = fmax(0.001, fmin(0.999, numerator / (2.0 - w1)));
  }
}

model {
  for (s in 1:S)
    w[s] ~ dirichlet(rep_vector(2.0, 2));

  for (s in 1:S) {
    for (i in 1:N) {
      if (seg_assignment[i] == s) {
        vector[2] contributions;
        for (k in 1:2)
          contributions[k] = log(w[s][k]) + binomial_lpmf(NV[i] | DP[i], peaks[s, k]);
        target += log_sum_exp(contributions);
      }
    }
  }
}

generated quantities {
  vector[S] log_lik;
  for (s in 1:S) {
    log_lik[s] = 0;
    for (i in 1:N) {
      if (seg_assignment[i] == s) {
        vector[2] contributions;
        for (k in 1:2)
          contributions[k] = log(w[s][k]) + binomial_lpmf(NV[i] | DP[i], peaks[s, k]);
        log_lik[s] += log_sum_exp(contributions);
      }
    }
  }
}
