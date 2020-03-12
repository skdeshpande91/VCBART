#ifndef _GUARD_BD_H_
#define _GUARD_BD_H_

#include<RcppArmadillo.h>
#include "rng.h"
#include "info.h"
#include "tree.h"
#include "funs.h"
#include "mu_posterior_new.h"
#include <stdio.h>

void bd_ind(tree &x, double &sigma, xinfo &xi, data_info &di, tree_prior_info &tree_pi, RNG &gen); // for independent error correlation structures
void bd_cs(tree &x, double &rho, double &sigma, xinfo &xi, data_info &di, tree_prior_info &tree_pi, RNG &gen); // for CS error correlation structure
void bd_ar(tree &x, double &rho, double &sigma, xinfo &xi, data_info &di, tree_prior_info &tree_pi, RNG &gen); // for AR error correlation structure




/*
// birth/death step for when we have residual AR correlation structure
//     birth/death step if we have global sigma
double bd_ar(tree &x, std::vector<arma::mat> Omega_vec, double &sigma, xinfo &xi, data_info &di, tree_prior_info &tree_pi, RNG &gen);
//    birth/death step if we have individual sigma
double bd_ar(tree &x, std::vector<arma::mat> Omega_vec, std::vector<double> &sigma, xinfo &xi, data_info &di, tree_prior_info &tree_pi, RNG &gen);

// birth/death step for when we there is no residual correlation structure

double bd_ind(tree &x, double &sigma, xinfo &xi, data_info &di, tree_prior_info &tree_pi, RNG &gen);


// birth/death when we marginalize out the random effect
void bd_cs(tree &x, double &rho, double &sigma, xinfo &xi, data_info &di, tree_prior_info &tree_pi, RNG &gen);
*/

#endif
