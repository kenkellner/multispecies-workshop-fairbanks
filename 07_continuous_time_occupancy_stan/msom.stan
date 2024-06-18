data{
  int J;                // number of replicate surveys
  int N;                // number of sites
  int S;                // number of species
  int K_site;           // number of occupancy covariates
  int K_dets;           // number of detection covariates
  array[S, N, J] int y; // detection / non-detection data
  array[N] row_vector [K_site + 1] x_s1; // species 1 occupancy covariates
  array[N] row_vector [K_site + 1] x_s2; // species 2 occupancy covariates
  array[N] row_vector [K_site + 1] x_in; // species interaction covariates
  array[S, N, J] row_vector[K_dets + 1] det_covs; // detection covariates
}
parameters{
  array[to_int(2 ^ S - 1)] vector[K_site + 1] beta_obs; // occupancy covariates
  array[S] vector[K_dets + 1] alpha_det;                // detection covariates
}
transformed parameters{
  
  array[N] vector [to_int(2 ^ S - 1)] f;   // natural parameters
  array[N] vector [to_int(2 ^ S)] log_psi; // log occupancy probability
  array[S] matrix[N, J] p;                 // detection probability 
  array[S] matrix[N, J] lpmf_y;            // log probability observed y
  
  for(i in 1:N){ // looping through all sites
    
    f[i][1] = x_s1[i] * beta_obs[1];
    f[i][2] = x_s2[i] * beta_obs[2];
    f[i][3] = x_in[i] * beta_obs[3];
    
    log_psi[i] = log_softmax([sum(f[i]), f[i][1], f[i][2], 0]');
    
    for(j in 1:J){ // looping through all replicate surveys
    
      p[1][i, j] = inv_logit(det_covs[1, i, j] * alpha_det[1]);
      lpmf_y[1][i, j] = bernoulli_lpmf(y[1, i, j] | p[1][i, j]);
      
      p[2][i, j] = inv_logit(det_covs[2, i, j] * alpha_det[2]);
      lpmf_y[2][i, j] = bernoulli_lpmf(y[2, i, j] | p[2][i, j]);
      
    }
  }
}
model{

  array[N] vector[to_int(2 ^ S)] log_elem; // log of each likelihood component
  
  for(i in 1:N){
  
    log_elem[i][1] = log_psi[i][1] + sum(lpmf_y[1][i, 1:J]) +
                      sum(lpmf_y[2][i, 1:J]);
    log_elem[i][2] = log_psi[i][2] + sum(lpmf_y[1][i, 1:J]) +
                      log(max(y[2, i, 1:J]) == 0);
    log_elem[i][3] = log_psi[i][3] + sum(lpmf_y[2][i, 1:J]) +
                      log(max(y[1, i, 1:J]) == 0);
    log_elem[i][4] = log_psi[i][4] + log(max(y[1, i, 1:J]) == 0) +
                      log(max(y[2, i, 1:J]) == 0);
    
    target += log_sum_exp(log_elem[i]);
 }
}
