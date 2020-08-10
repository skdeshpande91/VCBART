//
//  update_alpha.cpp
//    Update the alpha_z parameter
//  Created by Sameer Deshpande on 4/7/20.


#include "update_alpha_z.h"

void update_alpha_z(double &alpha_z, double &rho_alpha, std::vector<double> &theta, size_t &R, size_t &N_u, double &a, double &b, RNG &gen){
  
  double sum_log_theta = 0.0;
  for(size_t r = 0; r < R; r++) sum_log_theta += log(theta[r]);
  
  std::vector<double> alpha_log_prob(N_u - 1);
  std::vector<double> alpha_prob(N_u - 1);
  double max_log_prob = 0.0;
  
  double tmp_u = 0.0;
  double tmp_alpha = 0.0;
  double tmp_sum = 0.0;
  
  for(size_t u_ix = 0; u_ix < N_u - 1; u_ix++){
    tmp_u = ( (double) u_ix + 1.0)/( (double) N_u);
    tmp_alpha = rho_alpha * tmp_u/(1.0 - tmp_u);
    
    alpha_log_prob[u_ix] = lgamma(tmp_alpha) - ( (double) R) * lgamma(tmp_alpha/( (double) R));
    alpha_log_prob[u_ix] += tmp_alpha/( (double) R) * sum_log_theta;
    alpha_log_prob[u_ix] += (a - 1.0) * log(tmp_u) + (b + 1) * log(1.0 - tmp_u);
    
    if( (u_ix == 0) || (max_log_prob < alpha_log_prob[u_ix]) ) max_log_prob = alpha_log_prob[u_ix];
  }
  
  for(size_t u_ix = 0; u_ix < N_u - 1; u_ix++){
    alpha_prob[u_ix] = exp(alpha_log_prob[u_ix] - max_log_prob);
    tmp_sum += alpha_prob[u_ix];
  }
  
  for(size_t u_ix = 0; u_ix < N_u - 1; u_ix++) alpha_prob[u_ix] /= tmp_sum;
  
  size_t num_alphas = N_u - 1;
  size_t u_ix_new = draw_multinomial(num_alphas, alpha_prob, gen);
  alpha_z = rho_alpha * ( (double) u_ix_new + 1.0)/( (double) N_u);
}

