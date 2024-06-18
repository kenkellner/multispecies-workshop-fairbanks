data{

  int S;                // number of species
  int N;                // number of sites
  int K_site;           // number of occupancy covariates
  int K_dets;           // number of detection covariates
  int J;                // number of split-time intervals
  int W;                // number of wtd time intervals
  int C;                // number of coy time intervals
  
  array[N] row_vector [K_site + 1] x_s1; // species 1 occupancy covariates
  array[N] row_vector [K_site + 1] x_s2; // species 2 occupancy covariates
  array[N] row_vector [K_site + 1] x_in; // species interaction covariates
  array[N] int start_j;                  // start index of intervals
  array[N] int end_j;                    // end index of intervals
  
  array[J] int J_wtd;     // no. detections in each split time interval
  array[J] int J_coy;     // no. detections in each split time interval
  array[J] int start_wtd; // index of 1st element of interval vector / site
  array[J] int start_coy; // index of 1st element of interval vector / site
  array[J] int end_wtd;   // index of last element of interval vector / site
  array[J] int end_coy;   // index of last element of interval vector / site
  array[J] row_vector[K_dets + 1] x_det; // detection covariates
  
  vector[W] delta_wtd; // time between detections
  vector[C] delta_coy; // time between detections
  
}

parameters{

  array[S + 1] vector[K_site + 1] beta_obs; // occupancy covariates
  array[S] vector[K_dets + 1] alpha_det;    // detection covariates
  array[S + 1] vector[2] mu;                // parameters for Q matrix
  
}
transformed parameters{

  array[N] vector [to_int(2 ^ S - 1)] f;   // natural parameters
  array[N] vector [to_int(2 ^ S)] log_psi; // log occupancy probability
  
  array[S + 1] row_vector[2] p;        // probability of being in either state
  array[S + 1] matrix[2, 2] Q;         // generator matrix
  
  vector[J] pdf_wtd;               // prob density, wtd detections
  vector[J] pdf_wtd_coy;           // prob density, wtd and coy
  vector[J] pdf_coy;               // prob density, coy detections
  
  p[1] = [ exp(mu[1][2]) / sum(exp(mu[1])), exp(mu[1][1]) / sum(exp(mu[1])) ];
  p[2] = [ exp(mu[2][2]) / sum(exp(mu[2])), exp(mu[2][1]) / sum(exp(mu[2])) ];
  p[3] = [ exp(mu[3][2]) / sum(exp(mu[3])), exp(mu[3][1]) / sum(exp(mu[3])) ];
  
  Q[1] = [ [ -exp(mu[1][1]), exp(mu[1][1]) ],
           [ exp(mu[1][2]), -exp(mu[1][2]) ] ]; 
  Q[2] = [ [ -exp(mu[2][1]), exp(mu[2][1]) ],
           [ exp(mu[2][2]), -exp(mu[2][2]) ] ];
  Q[3] = [ [ -exp(mu[3][1]), exp(mu[3][1]) ],
           [ exp(mu[3][2]), -exp(mu[3][2]) ] ];
  
  for(i in 1:N){
    
    f[i][1] = x_s1[i] * beta_obs[1];
    f[i][2] = x_s2[i] * beta_obs[2];
    f[i][3] = x_in[i] * beta_obs[3];
    log_psi[i] = log_softmax([sum(f[i]), f[i][1], f[i][2], 0]');
    
  }
  
  for(j in 1:J){ // looping through all split times
      
    if(J_wtd[j] == 0){
    
      pdf_wtd[j] = p[1] *
        matrix_exp((Q[1] - [[0, 0], [0, exp(x_det[j] * alpha_det[1])]]) *
        delta_wtd[start_wtd[j]]) * [1, 1]';
      pdf_wtd_coy[j] = p[3] *
        matrix_exp((Q[3] - [[0, 0], [0, exp(x_det[j] * alpha_det[1])]]) *
        delta_wtd[start_wtd[j]]) * [1, 1]';
        
    }
      
    else if(J_wtd[j] == 1){
    
      pdf_wtd[j] = p[1] *
        matrix_exp((Q[1] - [[0, 0], [0, exp(x_det[j] * alpha_det[1])]]) *
        delta_wtd[start_wtd[j]]) *
        [[0, 0], [0, exp(x_det[j] * alpha_det[1])]] *
        matrix_exp((Q[1] - [[0, 0], [0, exp(x_det[j] * alpha_det[1])]]) *
        delta_wtd[start_wtd[j] + 1]) * [1, 1]';
      pdf_wtd_coy[j] = p[3] *
        matrix_exp((Q[3] - [[0, 0], [0, exp(x_det[j] * alpha_det[1])]]) *
        delta_wtd[start_wtd[j]]) *
        [[0, 0], [0, exp(x_det[j] * alpha_det[1])]] *
        matrix_exp((Q[3] - [[0, 0], [0, exp(x_det[j] * alpha_det[1])]]) *
        delta_wtd[start_wtd[j] + 1]) * [1, 1]';
        
    }
      
    else {
    
      array[J_wtd[j]] row_vector[2] pdf_wtd_tmp;
      array[J_wtd[j]] row_vector[2] pdf_wtd_coy_tmp;
        
      pdf_wtd_tmp[1] = p[1] *
        matrix_exp((Q[1] - [[0, 0], [0, exp(x_det[j] * alpha_det[1])]]) *
        delta_wtd[start_wtd[j]]) *
        [[0, 0], [0, exp(x_det[j] * alpha_det[1])]];
      pdf_wtd_coy_tmp[1] = p[3] *
        matrix_exp((Q[3] - [[0, 0], [0, exp(x_det[j] * alpha_det[1])]]) *
        delta_wtd[start_wtd[j]]) *
        [[0, 0], [0, exp(x_det[j] * alpha_det[1])]];
        
      for(s in 2:J_wtd[j]){
      
        pdf_wtd_tmp[s] = pdf_wtd_tmp[s - 1] *
          matrix_exp((Q[1] - [[0, 0], [0, exp(x_det[j] * alpha_det[1])]]) *
          delta_wtd[start_wtd[j] + (s - 1)]) *
          [[0, 0], [0, exp(x_det[j] * alpha_det[1])]];
        pdf_wtd_coy_tmp[s] = pdf_wtd_coy_tmp[s - 1] *
          matrix_exp((Q[3] - [[0, 0], [0, exp(x_det[j] * alpha_det[1])]]) *
          delta_wtd[start_wtd[j] + (s - 1)]) *
          [[0, 0], [0, exp(x_det[j] * alpha_det[1])]];
          
      }
        
      pdf_wtd[j] = pdf_wtd_tmp[J_wtd[j]] *
        matrix_exp((Q[1] - [[0, 0], [0, exp(x_det[j] * alpha_det[1])]]) *
        delta_wtd[end_wtd[j]]) * [1, 1]';
      pdf_wtd_coy[j] = pdf_wtd_coy_tmp[J_wtd[j]] *
        matrix_exp((Q[3] - [[0, 0], [0, exp(x_det[j] * alpha_det[1])]]) *
        delta_wtd[end_wtd[j]]) * [1, 1]';
      
    }
    
    // coyote detection intensity
    if(J_coy[j] == 0){
    
      pdf_coy[j] = p[2] *
        matrix_exp((Q[2] - [[0, 0], [0, exp(x_det[j] * alpha_det[2])]]) *
        delta_coy[start_coy[j]]) *
        [1, 1]';
        
    }
      
    else if(J_coy[j] == 1){
    
      pdf_coy[j] = p[2] *
        matrix_exp((Q[2] - [[0, 0], [0, exp(x_det[j] * alpha_det[2])]]) *
        delta_coy[start_coy[j]]) *
        [[0, 0], [0, exp(x_det[j] * alpha_det[2])]] *
        matrix_exp((Q[2] - [[0, 0], [0, exp(x_det[j] * alpha_det[2])]]) *
        delta_coy[start_coy[j] + 1]) * [1, 1]';
        
    }
      
    else {
    
      array[J_coy[j]] row_vector[2] pdf_coy_tmp;
        
      pdf_coy_tmp[1] = p[2] *
        matrix_exp((Q[2] - [[0, 0], [0, exp(x_det[j] * alpha_det[2])]]) *
        delta_coy[start_coy[j]]) *
        [[0, 0], [0, exp(x_det[j] * alpha_det[2])]];
        
      for(s in 2:J_coy[j]){
      
        pdf_coy_tmp[s] = pdf_coy_tmp[s - 1] *
          matrix_exp((Q[2] - [[0, 0], [0, exp(x_det[j] * alpha_det[2])]]) *
          delta_coy[start_coy[j] + (s - 1)]) *
          [[0, 0], [0, exp(x_det[j] * alpha_det[2])]];
          
      }
        
      pdf_coy[j] = pdf_coy_tmp[J_coy[j]] *
        matrix_exp((Q[2] - [[0, 0], [0, exp(x_det[j] * alpha_det[2])]]) *
        delta_coy[end_coy[j]]) *
        [1, 1]';
      
    }
  }
}
model{
  
  array[N] vector[to_int(2 ^ S)] log_elem; // log of each likelihood component
  
  for(i in 1:N){ // looping through all sites
    
    log_elem[i][1] = log_psi[i][1] +
                      sum(log(pdf_wtd_coy[start_j[i]:end_j[i]])) +
                      sum(log(pdf_coy[start_j[i]:end_j[i]]));
                      
    log_elem[i][2] = log_psi[i][2] +
                      sum(log(pdf_wtd[start_j[i]:end_j[i]])) +
                      log(sum(J_coy[start_j[i]:end_j[i]]) == 0);
                      
    log_elem[i][3] = log_psi[i][3] +
                      sum(log(pdf_coy[start_j[i]:end_j[i]])) +
                      log(sum(J_wtd[start_j[i]:end_j[i]]) == 0);
                      
    log_elem[i][4] = log_psi[i][4] +
                      log(sum(J_wtd[start_j[i]:end_j[i]]) == 0) +
                      log(sum(J_coy[start_j[i]:end_j[i]]) == 0);
    
    target += log_sum_exp(log_elem[i]);
    
  }
  
}
