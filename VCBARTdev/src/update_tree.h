//
//  update_tree.h
//  
//
//  Created by Sameer Deshpande on 4/7/20.
//

#ifndef _GUARD_update_tree_h
#define _GUARD_update_tree_h

#include<RcppArmadillo.h>
#include "rng.h"
#include "info.h"
#include "tree.h"
#include "funs.h"
#include "mu_posterior.h"
#include "draw_dirichlet.h"
#include "draw_multinomial.h"
#include <stdio.h>


void update_tree_ind(tree &x, double &sigma, std::vector<double> &theta, std::vector<size_t> &var_counts, xinfo &xi, data_info &di, tree_prior_info &tree_pi, RNG &gen);

void update_tree_cs(tree &x, double &sigma, double &rho, std::vector<double> &theta, std::vector<size_t> &var_counts, xinfo &xi, data_info &di, tree_prior_info &tree_pi, RNG &gen);


#endif /* update_tree_h */
