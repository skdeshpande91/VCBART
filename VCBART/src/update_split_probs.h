//
//  update_split_probs.h
//  
//
//  Created by Sameer Deshpande on 4/7/20.
//

#ifndef _GUARD_update_split_probs_h
#define _GUARD_update_split_probs_h
#include "rng.h"
#include "draw_dirichlet.h"
#include <stdio.h>


void update_split_probs(std::vector<double> &theta, std::vector<size_t> &var_counts, double &alpha_z, size_t &R, RNG &gen);



#endif /* update_split_probs_h */
