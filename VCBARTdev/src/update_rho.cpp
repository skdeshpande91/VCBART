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