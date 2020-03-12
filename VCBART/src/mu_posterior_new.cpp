//
//  mu_posterior.cpp
//  
//
//  Created by Sameer Deshpande on 2/3/20.
//

#include "mu_posterior_new.h"


// Compute posterior mean and covariance when errors have compound symmetry structure
void mu_posterior_cs(arma::vec &post_mean, arma::mat &post_cov_chol, double &ll, double &sigma, double &rho, std::vector<sinfo> &sv, data_info &di, tree_prior_info &tree_pi)
{
  size_t L = sv.size() ; // number of bottom nodes
  size_t j = 0;
  size_t jj = 0;
  post_mean.set_size(L); // set appropriate size for the mean
  post_cov_chol.set_size(L,L);
  ll = 0.0;
  arma::mat Lambda_inv = arma::zeros<arma::mat>(L,L);
  arma::vec theta = arma::zeros<arma::mat>(L);
  
  std::vector<double> xx_sum(L); // holds sum of x_ijk^2 within each leaf
  std::vector<double> xr_sum(L); // holds sum of x_ijk * rf_ij within each leaf
  std::vector<double> x_sum(L); // holds sum of x_ijk within each leaf
  //std::vector<double> r_sum(L); // holds sum of r_ij within each leaf
  double r_sum = 0.0;// holds the total sum of r_partial
  
  // note that each of these terms needs to be updated for every individual i
  
  for(size_t i = 0; i < di.N; i++){
    r_sum = 0.0;
    for(size_t l = 0; l < L; l++){
      xx_sum[l] = 0.0;
      xr_sum[l] = 0.0;
      x_sum[l] = 0.0;
      
      if(sv[l].node_count[i] > 0){
        for(size_t j_ix = 0; j_ix < sv[l].node_count[i]; j_ix++){
          j = sv[l].I[i][j_ix];
          xx_sum[l] += di.x[di.k + di.p * (di.start_index[i] + j)] * di.x[di.k + di.p * (di.start_index[i] + j)];
          xr_sum[l] += di.x[di.k + di.p * (di.start_index[i] + j)] * tree_pi.rp[di.start_index[i] + j];
          x_sum[l] += di.x[di.k + di.p * (di.start_index[i] + j)];
          r_sum += tree_pi.rp[di.start_index[i] + j];
        }
      }
    }
    // we now have enough to compute X_i' Omega X_i and X_i' Omega r_i for individual i.
    // we can therefore update elements of Lambda_inv now
    for(size_t l = 0; l < L; l++){
      Lambda_inv(l,l) += 1.0/(sigma * sigma) * xx_sum[l] * 1.0/ (1.0 - rho);
      for(size_t ll = 0; ll < L; ll++){
        Lambda_inv(l,ll) -= 1.0/(sigma * sigma) * rho / ((1.0 - rho) * (1.0 + ( (double) di.n[i] - 1.0) * rho)) * x_sum[l] * x_sum[ll];
      }
      theta(l) += 1.0/(sigma * sigma) * 1.0 / (1.0 - rho) * xr_sum[l] - 1.0/(sigma * sigma) * rho / ((1.0 - rho) * (1.0 + ( (double) di.n[i] - 1.0) * rho)) * x_sum[l] * r_sum;
    }
  }
  Lambda_inv.diag() += 1.0/(tree_pi.sigma_mu[di.k] * tree_pi.sigma_mu[di.k]);
  
  if(Lambda_inv.has_nan()) Rcpp::Rcout << "Lambda_inv has an nan" << std::endl;
  
  //Rcpp::Rcout << "theta = " << std::endl;
  //theta.print();
  
  //Rcpp::Rcout << "Lambda_inv = " << std::endl;
  //Lambda_inv.print();
  
  arma::mat Lambda = arma::inv_sympd(Lambda_inv);
  post_mean = Lambda * theta;
  
  //post_cov_chol = arma::chol(Lambda, "lower");
  bool chol_success = arma::chol(post_cov_chol, Lambda, "lower");
  if(chol_success == false){
    Rcpp::Rcout << "cholesky decomposition failed!" << std::endl;
    Lambda.print();
  }
  
  double tmp_log_det = 0.0;
  double tmp_log_det_sign = 0.0;
  arma::log_det(tmp_log_det, tmp_log_det_sign, Lambda);
  ll = 0.5 * tmp_log_det + 0.5 * arma::as_scalar(theta.t() * Lambda * theta);
}


// for testing purposes: let's compute the posterior mean and covariance matrix in two different ways
void mu_posterior_cs_mat(arma::vec &post_mean, arma::mat &post_cov_chol, double &ll, double &sigma, double &rho, std::vector<sinfo> &sv, data_info &di, tree_prior_info &tree_pi)
{
  size_t L = sv.size();
  
  
  arma::mat Lambda_inv = arma::zeros<arma::mat>(L,L);
  arma::vec theta = arma::zeros<arma::vec>(L);
  
  arma::mat X; // the n_i x L matrix
  arma::vec rp; // the vector of partial residuals for individual i
  arma::mat Sigma;
  arma::mat Omega;
  
  for(size_t i = 0; i < di.N; i++){
    Sigma.set_size(di.n[i], di.n[i]);
    Sigma.fill(rho);
    Sigma.diag() += 1.0 - rho;
    
    rp.zeros(di.n[i]);
    X.zeros(di.n[i], L);
    
    for(size_t j = 0; j < di.n[i]; j++){
      rp(j) = tree_pi.rp[di.start_index[i] + j];
    }
    
    Omega = arma::inv_sympd(Sigma);
    for(size_t l = 0; l < L; l++){
      if(sv[l].node_count[i] > 0){
        for(size_t j_ix = 0; j_ix < sv[l].node_count[i]; j_ix++){
          X(sv[l].I[i][j_ix],l) = di.x[di.k + di.p*(di.start_index[i] + sv[l].I[i][j_ix])];
        }
      }
    }
    
    Lambda_inv += X.t() * Omega * X * 1.0/(sigma * sigma);
    theta += X.t() * Omega * rp * 1.0/(sigma * sigma);
    
  }
  Lambda_inv.diag() += 1.0/(tree_pi.sigma_mu[di.k] * tree_pi.sigma_mu[di.k]);
  if(Lambda_inv.has_nan()) Rcpp::Rcout << "Lambda_inv has an nan" << std::endl;
  arma::mat Lambda = arma::inv_sympd(Lambda_inv);
  post_mean = Lambda * theta;
  
  //post_cov_chol = arma::chol(Lambda, "lower");
  bool chol_success = arma::chol(post_cov_chol, Lambda, "lower");
  if(chol_success == false){
    Rcpp::Rcout << "cholesky decomposition failed!" << std::endl;
    Lambda.print();
  }
  
  double tmp_log_det = 0.0;
  double tmp_log_det_sign = 0.0;
  arma::log_det(tmp_log_det, tmp_log_det_sign, Lambda);
  ll = 0.5 * tmp_log_det + 0.5 * arma::as_scalar(theta.t() * Lambda * theta);
}



// [SKD]: Warning (3 February 2020): There is something wrong in this function.
// For non-zero rho, we have non-zero differences between this and mu_posterior_ar_mat
// Eventually we will need to fix this but in the mean time, let us stand up the new implementation using the *_mat version

void mu_posterior_ar(arma::vec &post_mean, arma::mat &post_cov_chol, double &ll, double &sigma, double &rho, std::vector<sinfo> &sv, data_info &di, tree_prior_info &tree_pi)
{
  size_t L = sv.size();
  
  
  arma::mat Lambda_inv = arma::zeros<arma::mat>(L,L);
  arma::vec theta = arma::zeros<arma::vec>(L);
  
  //std::vector<std::vector<size_t> > leaf_obs_map(di.N);
  
  std::vector<size_t> leaf_obs_map;

  for(size_t i = 0; i < di.N; i++){
    leaf_obs_map.clear();
    leaf_obs_map.resize(di.n[i]);
    for(size_t l = 0; l < L; l++){
      if(sv[l].node_count[i] > 0){
        for(size_t j_ix = 0; j_ix < sv[l].I[i].size(); j_ix++){
          leaf_obs_map[sv[l].I[i][j_ix]] = l;
        }
      }
    }
    for(size_t j = 0; j < di.n[i]; j++){
      if(j == 0){
        Lambda_inv(leaf_obs_map[j], leaf_obs_map[j]) += di.x[di.k + di.p*(di.start_index[i] + j)] * di.x[di.k + di.p*(di.start_index[i] + j)] / (1 - rho * rho) * 1.0/(sigma * sigma);
        Lambda_inv(leaf_obs_map[j], leaf_obs_map[j+1]) += di.x[di.k + di.p*(di.start_index[i] + j)] * di.x[di.k + di.p*(di.start_index[i] + j+1)] * -1.0 * rho/(1 - rho * rho) * 1.0/(sigma * sigma);
        //Lambda_inv(leaf_obs_map[j+1], leaf_obs_map[j]) += di.x[di.k + di.p*(di.start_index[i] + j)] * di.x[di.k + di.p*(di.start_index[i] + j+1)] * -1.0 * rho/(1 - rho * rho) * 1.0/(sigma * sigma);
        
        theta(leaf_obs_map[j]) += di.x[di.k + di.p * di.start_index[i] + j] * tree_pi.rp[di.start_index[i] + j] * 1.0/(1.0 - rho * rho) * 1.0/(sigma * sigma);

        theta(leaf_obs_map[j]) += di.x[di.k + di.p * di.start_index[i] + j] * tree_pi.rp[di.start_index[i] + j + 1] * -1.0 * rho/(1 - rho * rho) * 1.0/(sigma * sigma);
        
      } else if(j == di.n[i] - 1){
        Lambda_inv(leaf_obs_map[j], leaf_obs_map[j]) += di.x[di.k + di.p*(di.start_index[i] + j)] * di.x[di.k + di.p*(di.start_index[i] + j)] / (1 - rho * rho) * 1.0/(sigma * sigma);
        Lambda_inv(leaf_obs_map[j], leaf_obs_map[j-1]) += di.x[di.k + di.p*(di.start_index[i] + j)] * di.x[di.k + di.p*(di.start_index[i] + j-1)] * -1.0 * rho/(1 - rho * rho) * 1.0/(sigma * sigma);
        //Lambda_inv(leaf_obs_map[j-1], leaf_obs_map[j]) += di.x[di.k + di.p*(di.start_index[i] + j)] * di.x[di.k + di.p*(di.start_index[i] + j-1)] * -1.0 * rho/(1 - rho * rho) * 1.0/(sigma * sigma);
        
        theta(leaf_obs_map[j]) += di.x[di.k + di.p * di.start_index[i] + j] * tree_pi.rp[di.start_index[i] + j] * 1.0/(1.0 - rho * rho) * 1.0/(sigma * sigma);
        theta(leaf_obs_map[j]) += di.x[di.k + di.p * di.start_index[i] + j] * tree_pi.rp[di.start_index[i] + j-1] * -1.0 * rho/(1.0 - rho * rho) * 1.0/(sigma * sigma);

        
      } else{
        Lambda_inv(leaf_obs_map[j], leaf_obs_map[j]) += di.x[di.k + di.p*(di.start_index[i] + j)] * di.x[di.k + di.p*(di.start_index[i] + j)] * (1 + rho * rho)/(1 - rho * rho) * 1.0/(sigma * sigma);
        
        Lambda_inv(leaf_obs_map[j], leaf_obs_map[j-1]) += di.x[di.k + di.p*(di.start_index[i] + j)] * di.x[di.k + di.p*(di.start_index[i] + j-1)] * -1.0 * rho/(1 - rho * rho) * 1.0/(sigma * sigma);
        Lambda_inv(leaf_obs_map[j], leaf_obs_map[j+1]) += di.x[di.k + di.p*(di.start_index[i] + j)] * di.x[di.k + di.p*(di.start_index[i] + j+ 1)] * -1.0 * rho/(1- rho * rho) * 1.0/(sigma * sigma);
        
        //Lambda_inv(leaf_obs_map[j-1], leaf_obs_map[j]) += di.x[di.k + di.p*(di.start_index[i] + j)] * di.x[di.k + di.p*(di.start_index[i] + j-1)] * -1.0 * rho/(1 - rho * rho) * 1.0/(sigma * sigma);
        //Lambda_inv(leaf_obs_map[j+1], leaf_obs_map[j]) += di.x[di.k + di.p*(di.start_index[i] + j)] * di.x[di.k + di.p*(di.start_index[i] + j+ 1)] * -1.0 * rho/(1- rho * rho) * 1.0/(sigma * sigma);
        
        theta(leaf_obs_map[j]) += di.x[di.k + di.p * di.start_index[i] + j] * tree_pi.rp[di.start_index[i] + j] * (1 + rho * rho)/(1 - rho * rho) * 1.0/(sigma * sigma);
        theta(leaf_obs_map[j]) += di.x[di.k + di.p * di.start_index[i] + j] * tree_pi.rp[di.start_index[i] + j-1] * -1.0 * rho/(1 - rho * rho) * 1.0/(sigma * sigma);
        theta(leaf_obs_map[j]) += di.x[di.k + di.p * di.start_index[i] + j] * tree_pi.rp[di.start_index[i] + j+1] * -1.0 * rho/(1 - rho * rho) * 1.0/(sigma * sigma);
      }
    }
  }

  Lambda_inv.diag() += 1.0/(tree_pi.sigma_mu[di.k] * tree_pi.sigma_mu[di.k]);
  if(Lambda_inv.has_nan()) Rcpp::Rcout << "Lambda_inv has an nan" << std::endl;
  arma::mat Lambda = arma::inv_sympd(Lambda_inv);
  post_mean = Lambda * theta;
  
  //post_cov_chol = arma::chol(Lambda, "lower");
  bool chol_success = arma::chol(post_cov_chol, Lambda, "lower");
  if(chol_success == false){
    Rcpp::Rcout << "cholesky decomposition failed!" << std::endl;
    Lambda.print();
  }
  
  double tmp_log_det = 0.0;
  double tmp_log_det_sign = 0.0;
  arma::log_det(tmp_log_det, tmp_log_det_sign, Lambda);
  ll = 0.5 * tmp_log_det + 0.5 * arma::as_scalar(theta.t() * Lambda * theta);
  
  // in mu_posterior_cs we have a nested loop: loop over individuals and then loop over leafs, and compute X'X and X'r
  // in this new setting, we will loop over leafs and update a lookup table that tells us which leafs each observation is in
}

void mu_posterior_ar_mat(arma::vec &post_mean, arma::mat &post_cov_chol, double &ll, double &sigma, double &rho, std::vector<sinfo> &sv, data_info &di, tree_prior_info &tree_pi)
{
  size_t L = sv.size();
  
  
  arma::mat Lambda_inv = arma::zeros<arma::mat>(L,L);
  arma::vec theta = arma::zeros<arma::vec>(L);
  
  arma::mat X; // the n_i x L matrix
  arma::vec rp; // the vector of partial residuals for individual i
  arma::mat Sigma;
  arma::mat Omega;
  
  for(size_t i = 0; i < di.N; i++){
    Sigma.set_size(di.n[i], di.n[i]);
    rp.zeros(di.n[i]);
    X.zeros(di.n[i], L);
    
    for(size_t j = 0; j < di.n[i]; j++){
      for(size_t jj = j; jj < di.n[i]; jj++){
        Sigma(j,jj) = pow(rho, jj-j); // remember that jj > j so this is always >= 0
        Sigma(jj,j) = pow(rho, jj - j); // remember that jj > j so this is always >= 0
      }
      rp(j) = tree_pi.rp[di.start_index[i] + j];
    }
    Omega = arma::inv_sympd(Sigma);
    for(size_t l = 0; l < L; l++){
      if(sv[l].node_count[i] > 0){
        for(size_t j_ix = 0; j_ix < sv[l].node_count[i]; j_ix++){
          X(sv[l].I[i][j_ix],l) = di.x[di.k + di.p*(di.start_index[i] + sv[l].I[i][j_ix])];
        }
      }
    }
    
    Lambda_inv += X.t() * Omega * X * 1.0/(sigma * sigma);
    theta += X.t() * Omega * rp * 1.0/(sigma * sigma);
  }
  Lambda_inv.diag() += 1.0/(tree_pi.sigma_mu[di.k] * tree_pi.sigma_mu[di.k]);
  if(Lambda_inv.has_nan()) Rcpp::Rcout << "Lambda_inv has an nan" << std::endl;
  arma::mat Lambda = arma::inv_sympd(Lambda_inv);
  post_mean = Lambda * theta;
  
  //post_cov_chol = arma::chol(Lambda, "lower");
  bool chol_success = arma::chol(post_cov_chol, Lambda, "lower");
  if(chol_success == false){
    Rcpp::Rcout << "cholesky decomposition failed!" << std::endl;
    Lambda.print();
  }
  
  double tmp_log_det = 0.0;
  double tmp_log_det_sign = 0.0;
  arma::log_det(tmp_log_det, tmp_log_det_sign, Lambda);
  ll = 0.5 * tmp_log_det + 0.5 * arma::as_scalar(theta.t() * Lambda * theta);
}


void mu_posterior_ind(arma::vec &post_mean, arma::mat &post_cov_chol, double &ll, double &sigma, std::vector<sinfo> &sv, data_info &di, tree_prior_info &tree_pi)
{
  size_t L = sv.size(); // number of nodes
  post_mean.set_size(L); // set appropriate size for the mean
  post_cov_chol.set_size(L,L);
  arma::mat Lambda_inv = arma::zeros<arma::mat>(L,L);
  arma::vec theta = arma::zeros<arma::vec>(L);
  
  size_t j = 0;
  
  for(size_t l = 0; l < L; l++){
    for(size_t i = 0; i < di.N; i++){
      if(sv[l].node_count[i] > 0){
        for(size_t j_ix = 0; j_ix < sv[l].node_count[i]; j_ix++){
          j = sv[l].I[i][j_ix];
          Lambda_inv(l,l) += 1.0/(sigma * sigma) * di.x[di.k + di.p * (di.start_index[i] + j)] * di.x[di.k + di.p * (di.start_index[i] + j)];
          theta(l) += 1.0/(sigma * sigma) * di.x[di.k + di.p * (di.start_index[i] + j)] * tree_pi.rp[di.start_index[i] + j];
        }
      }
    }
  }
  Lambda_inv.diag() += 1.0/(tree_pi.sigma_mu[di.k] * tree_pi.sigma_mu[di.k]);
  if(Lambda_inv.has_nan()) Rcpp::Rcout << "Lambda_inv has an nan" << std::endl;
  arma::mat Lambda = arma::inv_sympd(Lambda_inv);
  post_mean = Lambda * theta;
  bool chol_success = arma::chol(post_cov_chol, Lambda, "lower");
  if(chol_success == false){
    Rcpp::Rcout << "cholesky decomposition failed!" << std::endl;
    Lambda.print();
  }
  
  double tmp_log_det = 0.0;
  double tmp_log_det_sign = 0.0;
  arma::log_det(tmp_log_det, tmp_log_det_sign, Lambda);
  ll = 0.5 * tmp_log_det + 0.5 * arma::as_scalar(theta.t() * Lambda * theta);
}

void mu_posterior_ind_mat(arma::vec &post_mean, arma::mat &post_cov_chol, double &ll, double &sigma, std::vector<sinfo> &sv, data_info &di, tree_prior_info &tree_pi)
{
  size_t L = sv.size();
  arma::mat Lambda_inv = arma::zeros<arma::mat>(L,L);
  arma::vec theta = arma::zeros<arma::vec>(L);
  
  arma::mat X; // the n_i x L matrix
  arma::vec rp; // the vector of partial residuals for individual i
  arma::mat Omega;
  
  for(size_t i = 0; i < di.N; i++){
    Omega.set_size(di.n[i], di.n[i]);
    Omega.eye();
    
    rp.zeros(di.n[i]);
    X.zeros(di.n[i], L);
    
    for(size_t j = 0; j < di.n[i]; j++){
      rp(j) = tree_pi.rp[di.start_index[i] + j];
    }

    for(size_t l = 0; l < L; l++){
      if(sv[l].node_count[i] > 0){
        for(size_t j_ix = 0; j_ix < sv[l].node_count[i]; j_ix++){
          X(sv[l].I[i][j_ix],l) = di.x[di.k + di.p*(di.start_index[i] + sv[l].I[i][j_ix])];
        }
      }
    }
    
    Lambda_inv += X.t() * Omega * X * 1.0/(sigma * sigma);
    theta += X.t() * Omega * rp * 1.0/(sigma * sigma);
    
  }
  Lambda_inv.diag() += 1.0/(tree_pi.sigma_mu[di.k] * tree_pi.sigma_mu[di.k]);
  if(Lambda_inv.has_nan()) Rcpp::Rcout << "Lambda_inv has an nan" << std::endl;
  arma::mat Lambda = arma::inv_sympd(Lambda_inv);
  post_mean = Lambda * theta;
  
  //post_cov_chol = arma::chol(Lambda, "lower");
  bool chol_success = arma::chol(post_cov_chol, Lambda, "lower");
  if(chol_success == false){
    Rcpp::Rcout << "cholesky decomposition failed!" << std::endl;
    Lambda.print();
  }
  
  double tmp_log_det = 0.0;
  double tmp_log_det_sign = 0.0;
  arma::log_det(tmp_log_det, tmp_log_det_sign, Lambda);
  ll = 0.5 * tmp_log_det + 0.5 * arma::as_scalar(theta.t() * Lambda * theta);
}
