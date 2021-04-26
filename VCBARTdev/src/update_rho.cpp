//
//  update_rho.cpp
//  
//
//  Created by Cecilia Balocchi on 2/22/21.
//

#include "update_rho.h"

void update_rho_cs(double &rho_eps, double &sigma, data_info &di, RNG &gen)
{
  
  // we have a current value of rho, we propose a new one, we compute the MH ratio
  // work with xi = log(rho/(1-rho))

  double xi_curr = log(rho_eps/(1-rho_eps)); 
  double xi_prop = xi_curr + gen.normal(0, 0.1);
  double rho_prop = 1.0/(1+exp(-xi_prop));


  double lprior_curr = 0.0; // assume uniform prior for now
  double lprior_prop = 0.0; // assume uniform prior for now
  double lpost_curr = log(rho_eps) +log(1-rho_eps) + lprior_curr;
  double lpost_prop = log(rho_prop) +log(1-rho_prop) + lprior_prop;

  double tmp_sum, tmp_sum2; // hold running sum (and sum of squares) of full residuals
  double ss; // old sum of squares of the residuals and precision matrix (Ri^T Omega Ri)
  for(size_t i = 0; i < di.N; i++){
    // di.n[i]
    lpost_curr += -0.5 * ((double) di.n[i] * log(1-rho_eps) + log(1 + (double) di.n[i] * rho_eps/(1-rho_eps) ));
    lpost_prop += -0.5 * ((double) di.n[i] * log(1-rho_prop) + log(1 + (double) di.n[i] * rho_prop/(1-rho_prop) ));
    
    tmp_sum = 0.0;
    tmp_sum2 = 0.0;
    for(size_t j = 0; j < di.n[i]; j++){
      tmp_sum += di.rf[j + di.start_index[i]];
      tmp_sum2 += di.rf[j + di.start_index[i]] * di.rf[j + di.start_index[i]];
    }
    ss = tmp_sum2 / (1.0 - rho_eps) - tmp_sum * tmp_sum * rho_eps/( (1.0 - rho_eps) * (1.0 + ((double) di.n[i] - 1) * rho_eps));
    lpost_curr += -0.5 * ss / (sigma * sigma);
    ss = tmp_sum2 / (1.0 - rho_prop) - tmp_sum * tmp_sum * rho_prop/( (1.0 - rho_prop) * (1.0 + ((double) di.n[i] - 1) * rho_prop));
    lpost_prop += -0.5 * ss / (sigma * sigma);
  }
    
  double laccMH = lpost_prop - lpost_curr;
  if(log(gen.uniform()) < laccMH){
    rho_eps = rho_prop;
  }
  
}

void update_rho_cs(double &rho_eps, double &sigma, double &xi_sd, double &xi_mean, int iter, data_info &di, RNG &gen)
{
  // before running the loop, initialize: xi_sd = 0.1, xi_mean = rho_eps
  
  // work with xi = log(rho/(1-rho))

  double xi_curr = log(rho_eps/(1-rho_eps)); 
  // double xi_prop = xi_curr + gen.normal(0, 0.1);
  double xi_prop = xi_curr;
  if(gen.uniform() < 0.95){
    xi_prop += gen.normal(0, 2.38 * xi_sd);
  } else {
    xi_prop += gen.normal(0, 0.1);
  }

  double rho_prop = 1.0/(1+exp(-xi_prop));


  double lprior_curr = 0.0; // assume uniform prior for now
  double lprior_prop = 0.0; // assume uniform prior for now
  double lpost_curr = log(rho_eps) +log(1-rho_eps) + lprior_curr;
  double lpost_prop = log(rho_prop) +log(1-rho_prop) + lprior_prop;

  double tmp_sum, tmp_sum2; // hold running sum (and sum of squares) of full residuals
  double ss; // old sum of squares of the residuals and precision matrix (Ri^T Omega Ri)
  for(size_t i = 0; i < di.N; i++){
    // di.n[i]
    lpost_curr += -0.5 * ((double) di.n[i] * log(1-rho_eps) + log(1 + (double) di.n[i] * rho_eps/(1-rho_eps) ));
    lpost_prop += -0.5 * ((double) di.n[i] * log(1-rho_prop) + log(1 + (double) di.n[i] * rho_prop/(1-rho_prop) ));
    
    tmp_sum = 0.0;
    tmp_sum2 = 0.0;
    for(size_t j = 0; j < di.n[i]; j++){
      tmp_sum += di.rf[j + di.start_index[i]];
      tmp_sum2 += di.rf[j + di.start_index[i]] * di.rf[j + di.start_index[i]];
    }
    ss = tmp_sum2 / (1.0 - rho_eps) - tmp_sum * tmp_sum * rho_eps/( (1.0 - rho_eps) * (1.0 + ((double) di.n[i] - 1) * rho_eps));
    lpost_curr += -0.5 * ss / (sigma * sigma);
    ss = tmp_sum2 / (1.0 - rho_prop) - tmp_sum * tmp_sum * rho_prop/( (1.0 - rho_prop) * (1.0 + ((double) di.n[i] - 1) * rho_prop));
    lpost_prop += -0.5 * ss / (sigma * sigma);
  }
    
  double laccMH = lpost_prop - lpost_curr;
  if(log(gen.uniform()) < laccMH){
    rho_eps = rho_prop;
    xi_curr = xi_prop;
  }
  
  // update xi_sd
  // s2n = (n-2)/(n-1) s2n-1 + (1/n) * (Xn - Xbar_n-1)^2
  
  int n = iter + 2; // in iter = 0 we sample 2nd
  double xi_sd2 = (n-2)/(n-1) * xi_sd * xi_sd + (xi_curr - xi_mean)*(xi_curr - xi_mean)/n;
  xi_sd = sqrt(xi_sd2);
  xi_mean = (xi_mean * (n-1) + xi_curr)/n;

}