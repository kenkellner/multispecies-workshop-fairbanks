data{
  int J;             // maximum number of replicate surveys
  int N;             // number of sites
  int K_site;        // number of occupancy covariates
  int K_dets;        // number of detection covariates
  array[N, J] int y; // detection / non-detection matrix
  array[N] row_vector[K_site + 1] obs_covs;    // occupancy design matrix
  array[N, J] row_vector[K_dets + 1] det_covs; // detection design matrix
}
parameters{
  vector[K_site + 1] beta_obs;  // occupancy covariates
  vector[K_dets + 1] alpha_det; // detection covariates
}
transformed parameters{
  vector[N] psi;       // site occupancy probability
  matrix[N, J] p;      // detection probability
  matrix[N, J] lpmf_y; // log probability of observed detection
  
  for(i in 1:N){                  // looping through all sites
    psi[i] = inv_logit(obs_covs[i] * beta_obs);
    for(j in 1:J){                // looping through all replicate observations
      p[i, j] = inv_logit(det_covs[i, j] * alpha_det);
      lpmf_y[i, j] = bernoulli_lpmf(y[i, j] | p[i, j]);
    }
  }
  
}
model{
  array[N] vector[2] log_elem; // log of each likelihood component
  for(i in 1:N){
    log_elem[i][1] = log(psi[i]) + sum(lpmf_y[i, 1:J]);
    log_elem[i][2] = log(1 - psi[i]) + log(max(y[i, 1:J]) == 0);
    target += log_sum_exp(log_elem[i]);
  }
}
