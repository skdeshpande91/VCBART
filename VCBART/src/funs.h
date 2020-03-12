//
//  funs.h
//  
//
//  Created by Sameer Deshpande on 8/21/19.
//

#ifndef GUARD_funs_h
#define GUARD_funs_h
#include <RcppArmadillo.h>
#include <cmath>
#include "tree.h"
#include "info.h"
#include "kernels.h"
#include <stdio.h>

#endif 


//--------------------------------------------------
// center and scale each column of the X matrix and outcomes Y
void prepare_x(arma::mat &X_all, std::vector<double> &x_col_mean, std::vector<double> &x_col_sd);
void prepare_x(arma::mat &X_train, arma::mat &X_test, std::vector<double> &x_col_mean, std::vector<double> &x_col_sd);
void prepare_y(arma::vec &Y, double &y_mean, double &y_sd, double &y_max, double &y_min);
//void prepare_precision_ar(std::vector<std::vector<arma::mat> > &Omega_all, std::vector<std::vector<double> > &log_det_all, std::vector<double> rho, data_info &di);
//void prepare_precision_cs(std::vector<std::vector<arma::mat> > &Omega_all, std::vector<std::vector<double> > &log_det_all, std::vector<double> rho, data_info &di);
//void prepare_precision_ind(std::vector<std::vector<arma::mat> > &Omega_all, std::vector<std::vector<double> > &log_det_all, data_info &di);
//--------------------------------------------------

//does a (bottom) node have variables you can split on?
bool cansplit(tree::tree_p n, xinfo& xi);
//--------------------------------------------------
//compute prob of a birth, goodbots will contain all the good bottom nodes
double getpb(tree &t, xinfo &xi, tree_prior_info &tree_pi, tree::npv &goodbots);
//--------------------------------------------------
//find variables n can split on, put their indices in goodvars
void getgoodvars(tree::tree_p n, xinfo& xi,  std::vector<size_t>& goodvars);
//get prob a node grows, 0 if no good vars, else alpha/(1+d)^beta
//double pgrow(tree::tree_p n, xinfo &xi, tree_prior_info &tree_pi);
double pgrow(tree::tree_p n, xinfo &xi, tree_prior_info &tree_pi, size_t k); // 4 Nov 2019 -- tells us which beta function we're updating
//--------------------------------------------------
//get sufficients stats for all bottom nodes
void allsuff(tree &x, xinfo &xi, data_info &di, tree::npv &bnv, std::vector<sinfo> &sv);
//--------------------------------------------------
//get sufficient stats for children (v,c) of node nx in tree x
//[SKD]: used in the birth proposals
void getsuff(tree &x, tree::tree_cp nx, size_t v, size_t c, xinfo &xi, data_info &di, sinfo &sl, sinfo &sr, sinfo &st);
//--------------------------------------------------
//get sufficient stats for pair of bottom children nl(left) and nr(right) in tree x
//[SKD]: used in the death proposals
// st contains all of the observations in both nl and nr
void getsuff(tree &x, tree::tree_cp nl, tree::tree_cp nr, xinfo &xi, data_info &di, sinfo &sl, sinfo &sr, sinfo &st);

//--------------------------------------------------
// new functions for getting all sufficient statistics for use with bd_fast

// for birth proposal
void getsuff(tree &x, tree::tree_cp nx, size_t v, size_t c, xinfo &xi, data_info &di, std::vector<sinfo> &sv_x, std::vector<sinfo> &sv_y);
// for death proposal
void getsuff(tree &x, tree::tree_cp nl, tree::tree_cp nr, xinfo &xi, data_info &di, std::vector<sinfo> &sv_x, std::vector<sinfo> &sv_y);



//--------------------------------------------------
// fit
void fit(tree& t, xinfo& xi, data_info& di, double* fv);

//--------------------------------------------------
// update_sigma -- using inverse gamma prior
//void update_sigma_ig(double &sigma, std::vector<arma::mat> &Omega_vec, sigma_prior_info &sigma_pi, data_info &di, RNG &gen); // global sigma
//void update_sigma_ig(std::vector<double> &sigma, std::vector<arma::mat> &Omega_vec, std::vector<sigma_prior_info> &sigma_pi, data_info &di, RNG &gen); // individual sigma

//void update_sigma_ig(double &sigma, double &rho, sigma_prior_info &sigma_pi, data_info &di, RNG &gen);



// for when there is no residual correlation
//void update_sigma_ig(double &sigma, sigma_prior_info &sigma_pi, data_info &di, RNG &gen);

//--------------------------------------------------
// update_sigma -- using a half-cauchy prior
//void update_sigma_hc(double&sigma, std::vector<arma::mat> &Omega_vec, sigma_prior_info &sigma_pi, data_info &di, RNG &gen);
//void update_sigma_hc(std::vector<double> &sigma, std::vector<arma::mat> &Omega_vec, std::vector<sigma_prior_info> &sigma_pi, data_info &di, RNG &gen);

// half-cauchy prior on
//void update_sigma_hc(double &sigma, double &rho, sigma_prior_info &sigma_pi, data_info &di, RNG &gen);


// a more general half-t (we'll overload nu)
//void update_sigma_ht(double &sigma, double &rho, sigma_prior_info &sigma_pi, data_info &di, RNG &gen);

// update sigma w/ HC prior in the setting where there is no residual correlation
//void update_sigma_hc(double &sigma, sigma_prior_info &sigma_pi, data_info &di, RNG &gen);

//---------------------------------------------------
// update_sigma_mu
//void update_sigma_mu(double* sigma_mu, std::vector<std::vector<tree > > t_vec, tree_prior_info &tree_pi, RNG &gen);


//--------------------------------------------------
// update_rho

// update rho when we have global sigma
//void update_rho(size_t &rho_ix, std::vector<std::vector<arma::mat > > &Omega_all, std::vector<std::vector<double> > &log_det_all, double &sigma, data_info &di, RNG &gen);

// update rho when we have individual sigma's
//void update_rho(size_t &rho_ix, std::vector<std::vector<arma::mat> > &Omega_all, std::vector<std::vector<double> > &log_det_all, std::vector<double> &sigma, data_info &di, RNG &gen);


//---------------------------------------------------
// get the fit of a tree
void fit(tree& t, xinfo& xi, data_info& di, double* fv);



/*
//--------------------------------------------------
// compute the posterior mean and variance for ALL terminal node parameters
// this is only useful when there is residual correlation that induces dependence between leaves
void mu_posterior(arma::vec &post_mean, arma::mat &post_cov_chol, double &ll, double &sigma, const std::vector<arma::mat> &Omega_vec, std::vector<sinfo> &sv, data_info &di, tree_prior_info &tree_pi);

void mu_posterior(arma::vec &post_mean, arma::mat &post_cov_chol, double &ll, std::vector<double> sigma, const std::vector<arma::mat> &Omega_vec, std::vector<sinfo> &sv, data_info &di, tree_prior_info &tree_pi);
// sv is a vector of length # of bottom nodes.


// new posterior function
void mu_posterior(arma::vec &post_mean, arma::mat &post_cov_chol, double &ll, double &sigma, double &rho, std::vector<sinfo> &sv, data_info &di, tree_prior_info &tree_pi);


// compute posterior mean and variance for a single leaf --
// this is only useful when leafs are independent a posteriori. This happens eg when there is no residual correlation
void mu_posterior(double &post_mean, double &post_var, double &ll, double &sigma, sinfo &si, data_info &di, tree_prior_info &tree_pi);
*/
//void drmu(tree &t, double &sigma, xinfo &xi, data_info &di, tree_prior_info &tree_pi, RNG &gen);
//void update_alpha(double* alpha, double &sigma, double &sigma_alpha, data_info &di, RNG &gen);


