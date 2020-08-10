//
//  draw_dirichlet.cpp
//  
//
//  Created by Sameer Deshpande on 4/7/20.
//

#include "draw_dirichlet.h"

//a symmetric Dirichlet(alpha/R, ..., alpha/R) distribution
void draw_dirichlet(std::vector<double> &theta, double &alpha, size_t &R, RNG &gen)
{
  theta.clear();
  theta.resize(R);
  std::vector<double> tmp_gamma(R);
  double tmp_sum = 0.0;
  for(size_t r = 0; r < R; r++){
    tmp_gamma[r] = gen.gamma(alpha/( (double) R), 1.0);
    tmp_sum += tmp_gamma[r];
  }
  for(size_t r = 0; r < R; r++) theta[r] = tmp_gamma[r]/tmp_sum;
}


void draw_dirichlet(std::vector<double> &theta, std::vector<double> &alpha, size_t &R, RNG &gen)
{
  std::vector<double> tmp_gamma(R);
  double tmp_sum = 0.0;
  for(size_t r = 0; r < R; r++){
    tmp_gamma[r] = gen.gamma(alpha[r],1.0);
    tmp_sum += tmp_gamma[r];
  }
  for(size_t r = 0; r < R; r++) theta[r] = tmp_gamma[r]/tmp_sum;
}
