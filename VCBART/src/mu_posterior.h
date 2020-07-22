//
//  mu_posterior.h
//  
//
//  Created by Sameer Deshpande on 2/3/20.
//

#ifndef GUARD_mu_posterior_h
#define GUARD_mu_posterior_h

#include<RcppArmadillo.h>
#include<cmath>
#include "tree.h"
#include "info.h"
#include <stdio.h>

#endif

void mu_posterior_cs(arma::vec &post_mean, arma::mat &post_cov_chol, double &ll, double &sigma, double &rho, std::vector<sinfo> &sv, data_info &di, tree_prior_info &tree_pi);
void mu_posterior_cs_mat(arma::vec &post_mean, arma::mat &post_cov_chol, double &ll, double &sigma, double &rho, std::vector<sinfo> &sv, data_info &di, tree_prior_info &tree_pi);
void mu_posterior_ar(arma::vec &post_mean, arma::mat &post_cov_chol, double &ll, double &sigma, double &rho, std::vector<sinfo> &sv, data_info &di, tree_prior_info &tree_pi);
void mu_posterior_ar_mat(arma::vec &post_mean, arma::mat &post_cov_chol, double &ll, double &sigma, double &rho, std::vector<sinfo> &sv, data_info &di, tree_prior_info &tree_pi);
void mu_posterior_ind(arma::vec &post_mean, arma::mat &post_cov_chol, double &ll, double &sigma, std::vector<sinfo> &sv, data_info &di, tree_prior_info &tree_pi);
void mu_posterior_ind_mat(arma::vec &post_mean, arma::mat & post_cov_chol, double &ll, double &sigma, std::vector<sinfo> &sv, data_info &di, tree_prior_info &tree_pi);
