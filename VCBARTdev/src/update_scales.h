//
//  update_scales.hpp
//  
//
//  Created by Sameer Deshpande on 2/4/20.
//

#ifndef GUARD_update_scales_h
#define GUARD_update_scales_h

#include <RcppArmadillo.h>
#include <cmath>
#include "tree.h"
#include "info.h"
#include "funs.h"
#include <stdio.h>

#endif /* update_scales_hpp */


// update a global sigma with inverse gamma prior
void update_sigma_ig_ind(double &sigma, sigma_prior_info &sigma_pi, data_info &di, RNG &gen);
void update_sigma_ig_cs(double &sigma, double &rho_eps, sigma_prior_info &sigma_pi, data_info &di, RNG &gen);
void update_sigma_ig_ar(double &sigma, double &rho_eps, sigma_prior_info &sigma_pi, data_info &di, RNG &gen);


// update a global sigma with half-t prior
void update_sigma_ht_ind(double &sigma, sigma_prior_info &sigma_pi, data_info &di, RNG &gen);
void update_sigma_ht_cs(double &sigma, double &rho_eps, sigma_prior_info &sigma_pi, data_info &di, RNG &gen);
void update_sigma_ht_ar(double &sigma, double &rho_eps, sigma_prior_info &sigma_pi, data_info &di, RNG &gen);


// the matrix versions -- this is mainly for testing purposes
void update_sigma_ht_ind_mat(double &sigma, sigma_prior_info &sigma_pi, data_info &di, RNG &gen);
void update_sigma_ht_cs_mat(double &sigma, double &rho_eps, sigma_prior_info &sigma_pi, data_info &di, RNG &gen);
void update_sigma_ht_ar_mat(double &sigma, double &rho_eps, sigma_prior_info &sigma_pi, data_info &di, RNG &gen);


// update sigma_mu using a half-t prior
void update_sigma_mu(double* sigma_mu, std::vector<std::vector<tree > > t_vec, tree_prior_info & tree_pi, RNG &gen);
