data{
  int S;                // number of species
  int N;                // number of sites
  int K_site;           // number of occupancy covariates
  int K_dets;           // number of detection covariates
  int ints_wtd;         // number of wtd time intervals
  int ints_coy;         // number of coy time intervals
  array[N] row_vector [K_site + 1] x_s1; // species 1 occupancy covariates
  array[N] row_vector [K_site + 1] x_s2; // species 2 occupancy covariates
  array[N] row_vector [K_site + 1] x_in; // species interaction covariates
  array[N] int start_wtd; // index of 1st element of interval vector / site
  array[N] int start_coy; // index of 1st element of interval vector / site
  array[N] int end_wtd;   // index of last element of interval vector / site
  array[N] int end_coy;   // index of last element of interval vector / site
  vector[ints_wtd] delta_wtd; // time between detections
  vector[ints_coy] delta_coy; // time between detections
  vector[ints_wtd] y_wtd;     // quadrature / observation
  vector[ints_coy] y_coy;     // quadrature / observation
  array[ints_wtd] row_vector[K_dets + 1] x_wtd; // detection covariates
  array[ints_coy] row_vector[K_dets + 1] x_coy; // detection covariates
}
parameters{
  array[to_int(2 ^ S - 1)] vector[K_site + 1] beta_obs; // occupancy covariates
  array[S] vector[K_dets + 1] alpha_det; // detection covariates
}
transformed parameters{
  array[N] vector [to_int(2 ^ S - 1)] f;   // natural parameters
  array[N] vector [to_int(2 ^ S)] log_psi; // log occupancy probability
  vector[ints_wtd] lpmf_wtd;               // log prob vector, detections
  vector[ints_coy] lpmf_coy;               // log prob vector, detections
  
  for(i in 1:N){ // looping through all sites
    f[i][1] = x_s1[i] * beta_obs[1];
    f[i][2] = x_s2[i] * beta_obs[2];
    f[i][3] = x_in[i] * beta_obs[3];
    log_psi[i] = log_softmax([sum(f[i]), f[i][1], f[i][2], 0]');
    
    for(j in start_wtd[i]:end_wtd[i]){ // looping through all time intervals
      lpmf_wtd[j] = delta_wtd[j] *
        (y_wtd[j] * (x_wtd[j] * alpha_det[1]) - exp(x_wtd[j] * alpha_det[1]));
    }
    
    for(j in start_coy[i]:end_coy[i]){ // looping through all time intervals
      lpmf_coy[j] = delta_coy[j] *
        (y_coy[j] * (x_coy[j] * alpha_det[2]) - exp(x_coy[j] * alpha_det[2]));
    }
  }
}
model{

  array[N] vector[to_int(2 ^ S)] log_elem; // log of each likelihood component
  
  for(i in 1:N){
  
    log_elem[i][1] = log_psi[i][1] +
                      sum(lpmf_wtd[start_wtd[i]:end_wtd[i]]) +
                      sum(lpmf_coy[start_coy[i]:end_coy[i]]);
    log_elem[i][2] = log_psi[i][2] + sum(lpmf_wtd[start_wtd[i]:end_wtd[i]]) +
                      log(sum(y_coy[start_coy[i]:end_coy[i]]) == 0);
    log_elem[i][3] = log_psi[i][3] + sum(lpmf_coy[start_coy[i]:end_coy[i]]) +
                      log(sum(y_wtd[start_wtd[i]:end_wtd[i]]) == 0);
    log_elem[i][4] = log_psi[i][4] +
                      log(sum(y_wtd[start_wtd[i]:end_wtd[i]]) == 0) +
                      log(sum(y_coy[start_coy[i]:end_coy[i]]) == 0);
    
    target += log_sum_exp(log_elem[i]);
 }
}
