//
//  update_split_probs.cpp
//    Update the split probabilities, theta
//  Created by Sameer Deshpande on 4/7/20.


#include "update_split_probs.h"

void update_split_probs(std::vector<double> &theta, std::vector<size_t> &var_counts, double &alpha_z, size_t &R, RNG &gen)
{
  std::vector<double> alpha_new(R, alpha_z/( (double) R)); // this is basically the prior value
  for(size_t r = 0; r < R; r++) alpha_new[r] += (double) var_counts[r];
  draw_dirichlet(theta, alpha_new, R, gen);
}
