//
//  update_scales.cpp
//  
//
//  Created by Sameer Deshpande on 2/4/20.
//

#include "update_scales.h"

void update_sigma_ig_ind(double &sigma, sigma_prior_info &sigma_pi, data_info &di, RNG &gen)
{
  double scale_post = sigma_pi.lambda * sigma_pi.nu;
  double nu_post = sigma_pi.nu;
  for(size_t i = 0; i < di.N; i++){
    nu_post += di.n[i];
    for(size_t j = 0; j < di.n[i]; j++){
      scale_post += di.rf[j + di.start_index[i]] * di.rf[j + di.start_index[i]];
    }
  }
  sigma = sqrt( (scale_post)/gen.chi_square(nu_post));
}

void update_sigma_ht_ind(double &sigma, sigma_prior_info &sigma_pi, data_info &di, RNG &gen)
{
  double scale_post = 0.0;
  double nu_post = 0.0;
  
  for(size_t i = 0; i < di.N; i++){
    nu_post += di.n[i]; // add number of observations from individual i to the overall degrees of freedom
    for(size_t j = 0; j < di.n[i]; j++){
      scale_post += di.rf[j + di.start_index[i]] * di.rf[j + di.start_index[i]];
    }
  }
  //Rcpp::Rcout << "[ind]: nu_post = " << nu_post << "  scale_post = " << scale_post << std::endl;

  double sigma_prop = sqrt( (scale_post)/gen.chi_square(nu_post)); // proposed sigma'
  
  double log_prior_prop = -0.5 * (1.0 + sigma_pi.nu) * log(1.0 + 1.0/sigma_pi.nu * sigma_prop * sigma_prop/(sigma_pi.A * sigma_pi.A));
  double log_prior = -0.5 * (1.0 + sigma_pi.nu) * log(1.0 + 1.0/sigma_pi.nu * sigma * sigma/(sigma_pi.A * sigma_pi.A));
  double alpha = exp(3.0 * log(sigma_prop) + log_prior_prop - 3.0 * log(sigma) - log_prior);
  if(gen.uniform() <= alpha) sigma = sigma_prop;
}

void update_sigma_ht_ind_mat(double &sigma, sigma_prior_info &sigma_pi, data_info &di, RNG &gen)
{
  double scale_post = 0.0;
  double nu_post = 0.0;
  arma::mat Omega;
  arma::vec rf;
  
  for(size_t i = 0; i < di.N; i++){
    nu_post += di.n[i];
    rf.set_size(di.n[i]);
    Omega.set_size(di.n[i], di.n[i]);
    Omega.eye();
    for(size_t j = 0; j < di.n[i]; j++) rf(j) = di.rf[j + di.start_index[i]];
    
    scale_post += arma::as_scalar(rf.t() * Omega * rf);
  }
  //Rcpp::Rcout << "[ind_mat]: nu_post = " << nu_post << "  scale_post = " << scale_post << std::endl;
  
  double sigma_prop = sqrt( (scale_post)/gen.chi_square(nu_post)); // this is the proposed sigma'
  
  double log_prior_prop = -0.5 * (1.0 + sigma_pi.nu) * log(1.0 + 1.0/sigma_pi.nu * sigma_prop * sigma_prop/(sigma_pi.A * sigma_pi.A));
  double log_prior = -0.5 * (1.0 + sigma_pi.nu) * log(1.0 + 1.0/sigma_pi.nu * sigma * sigma/(sigma_pi.A * sigma_pi.A));
  double alpha = exp(3.0 * log(sigma_prop) + log_prior_prop - 3.0 * log(sigma) - log_prior);
  if(gen.uniform() <= alpha) sigma = sigma_prop;
  
}


void update_sigma_ig_cs(double &sigma, double &rho, sigma_prior_info &sigma_pi, data_info &di, RNG &gen)
{
  double scale_post = sigma_pi.lambda * sigma_pi.nu;
  double nu_post = sigma_pi.nu;
  
  double tmp_sum = 0.0; // hold running sum of full residuals
  double tmp_sum2 = 0.0; // hold running sum of squares of full residuals
  for(size_t i = 0; i < di.N; i++){
    nu_post += di.n[i];
    tmp_sum = 0.0;
    tmp_sum2 = 0.0;
    for(size_t j = 0; j < di.n[i]; j++){
      tmp_sum += di.rf[j + di.start_index[i]];
      tmp_sum2 += di.rf[j + di.start_index[i]] * di.rf[j + di.start_index[i]];
    }
    scale_post += tmp_sum2 / (1.0 - rho) - tmp_sum * tmp_sum * rho/( (1.0 - rho) * (1.0 + ((double) di.n[i] - 1) * rho));
  }
  sigma = sqrt( (scale_post)/gen.chi_square(nu_post));
}

void update_sigma_ht_cs(double &sigma, double &rho,sigma_prior_info &sigma_pi, data_info &di, RNG &gen)
{
  double scale_post = 0.0;
  double nu_post = 0.0;
  
  double tmp_sum = 0.0; // holds running sum of full residuals
  double tmp_sum2 = 0.0; // holds running sum of squares of full residuals
  
  for(size_t i = 0; i < di.N; i++){
    nu_post += di.n[i]; // add number of observations from individual i to the degrees of freedom
    
    tmp_sum = 0.0;
    tmp_sum2 = 0.0;
    for(size_t j = 0; j < di.n[i]; j++){
      tmp_sum += di.rf[j + di.start_index[i]];
      tmp_sum2 += di.rf[j + di.start_index[i]] * di.rf[j + di.start_index[i]];
    }
    scale_post += tmp_sum2 / (1.0 - rho) - tmp_sum * tmp_sum * rho/( (1.0 - rho) * (1.0 + ((double) di.n[i] - 1) * rho));
    if(scale_post != scale_post){
      Rcpp::Rcout <<"[update_sigma]: i = " << i << "  tmp_sum = " << tmp_sum << "  tmp_sum2 = " << tmp_sum2 << std::endl;
      Rcpp::stop("[update_sigma]: scale_post is nan. quitting");
    }
  }
  //Rcpp::Rcout << "[cs]: nu_post = " << nu_post << "  scale_post = " << scale_post << std::endl;

  double sigma_prop = sqrt( (scale_post)/gen.chi_square(nu_post)); // this is the proposed sigma'
  
  double log_prior_prop = -0.5 * (1.0 + sigma_pi.nu) * log(1.0 + 1.0/sigma_pi.nu * sigma_prop * sigma_prop/(sigma_pi.A * sigma_pi.A));
  double log_prior = -0.5 * (1.0 + sigma_pi.nu) * log(1.0 + 1.0/sigma_pi.nu * sigma * sigma/(sigma_pi.A * sigma_pi.A));
  double alpha = exp(3.0 * log(sigma_prop) + log_prior_prop - 3.0 * log(sigma) - log_prior);
  if(gen.uniform() <= alpha) sigma = sigma_prop;
}

void update_sigma_ht_cs_mat(double &sigma, double &rho, sigma_prior_info &sigma_pi, data_info &di, RNG &gen)
{
  double scale_post = 0.0;
  double nu_post = 0.0;
  
  arma::mat Sigma;
  arma::mat Omega;
  arma::vec rf;
  
  
  for(size_t i = 0; i < di.N; i++){
    nu_post += di.n[i];
    rf.set_size(di.n[i]);
    Sigma.set_size(di.n[i], di.n[i]);
    Sigma.fill(rho);
    Sigma.diag() += 1.0 - rho;
    Omega = arma::inv_sympd(Sigma);
    for(size_t j = 0; j < di.n[i]; j++) rf(j) = di.rf[j + di.start_index[i]];
    
    scale_post += arma::as_scalar(rf.t() * Omega * rf);
  }
  //Rcpp::Rcout << "[cs_mat]: nu_post = " << nu_post << "  scale_post = " << scale_post << std::endl;
  
  
  double sigma_prop = sqrt( (scale_post)/gen.chi_square(nu_post)); // this is the proposed sigma'
  
  double log_prior_prop = -0.5 * (1.0 + sigma_pi.nu) * log(1.0 + 1.0/sigma_pi.nu * sigma_prop * sigma_prop/(sigma_pi.A * sigma_pi.A));
  double log_prior = -0.5 * (1.0 + sigma_pi.nu) * log(1.0 + 1.0/sigma_pi.nu * sigma * sigma/(sigma_pi.A * sigma_pi.A));
  double alpha = exp(3.0 * log(sigma_prop) + log_prior_prop - 3.0 * log(sigma) - log_prior);
  if(gen.uniform() <= alpha) sigma = sigma_prop;
}


void update_sigma_ig_ar(double &sigma, double &rho, sigma_prior_info &sigma_pi, data_info & di, RNG &gen)
{
  double scale_post = sigma_pi.lambda * sigma_pi.nu;
  double nu_post = sigma_pi.nu;
  
  for(size_t i = 0; i < di.N; i++){
    nu_post += di.n[i];
    for(size_t j = 0; j < di.n[i]; j++){
      if(j == 0){
        scale_post += 1.0/(1 - rho * rho) * di.rf[j + di.start_index[i]] * di.rf[j + di.start_index[i]];
        scale_post += -1.0 * rho/(1 - rho * rho) * di.rf[j + di.start_index[i]] * di.rf[j + 1 + di.start_index[i]];
      } else if(j == di.n[i] - 1){
        scale_post += 1.0/(1.0 - rho * rho) * di.rf[j + di.start_index[i]] * di.rf[j + di.start_index[i]];
        scale_post += -1.0 * rho/(1.0 - rho * rho) * di.rf[j + di.start_index[i]] * di.rf[j - 1 + di.start_index[i]];
      } else{
        scale_post += (1 + rho * rho)/(1.0 - rho * rho) * di.rf[j + di.start_index[i]] * di.rf[j + di.start_index[i]];
        scale_post += -1.0 * rho /(1.0 - rho * rho) * di.rf[j + di.start_index[i]] * di.rf[j - 1 + di.start_index[i]];
        scale_post += -1.0 * rho/(1.0 - rho * rho) * di.rf[j + di.start_index[i]] * di.rf[j + 1 + di.start_index[i]];
      }
    }
  }
  sigma = sqrt( (scale_post)/gen.chi_square(nu_post));
}


void update_sigma_ht_ar(double &sigma, double &rho, sigma_prior_info &sigma_pi, data_info &di, RNG &gen)
{
  double scale_post = 0.0;
  double nu_post = 0.0;
  
  for(size_t i = 0; i < di.N; i++){
    nu_post += di.n[i];
    for(size_t j = 0; j < di.n[i]; j++){
      if(j == 0){
        scale_post += 1.0/(1 - rho * rho) * di.rf[j + di.start_index[i]] * di.rf[j + di.start_index[i]];
        scale_post += -1.0 * rho/(1 - rho * rho) * di.rf[j + di.start_index[i]] * di.rf[j + 1 + di.start_index[i]];
      } else if(j == di.n[i] - 1){
        scale_post += 1.0/(1.0 - rho * rho) * di.rf[j + di.start_index[i]] * di.rf[j + di.start_index[i]];
        scale_post += -1.0 * rho/(1.0 - rho * rho) * di.rf[j + di.start_index[i]] * di.rf[j - 1 + di.start_index[i]];
      } else{
        scale_post += (1 + rho * rho)/(1.0 - rho * rho) * di.rf[j + di.start_index[i]] * di.rf[j + di.start_index[i]];
        scale_post += -1.0 * rho /(1.0 - rho * rho) * di.rf[j + di.start_index[i]] * di.rf[j - 1 + di.start_index[i]];
        scale_post += -1.0 * rho/(1.0 - rho * rho) * di.rf[j + di.start_index[i]] * di.rf[j + 1 + di.start_index[i]];
      }
    }
  }
  //Rcpp::Rcout << "[ar]: nu_post = " << nu_post << "  scale_post = " << scale_post << std::endl;

  double sigma_prop = sqrt( (scale_post)/gen.chi_square(nu_post)); // this is the proposed sigma'
  
  double log_prior_prop = -0.5 * (1.0 + sigma_pi.nu) * log(1.0 + 1.0/sigma_pi.nu * sigma_prop * sigma_prop/(sigma_pi.A * sigma_pi.A));
  double log_prior = -0.5 * (1.0 + sigma_pi.nu) * log(1.0 + 1.0/sigma_pi.nu * sigma * sigma/(sigma_pi.A * sigma_pi.A));
  double alpha = exp(3.0 * log(sigma_prop) + log_prior_prop - 3.0 * log(sigma) - log_prior);
  if(gen.uniform() <= alpha) sigma = sigma_prop;
  
  
}

void update_sigma_ht_ar_mat(double &sigma, double &rho, sigma_prior_info &sigma_pi, data_info &di, RNG &gen)
{
  arma::mat Sigma;
  arma::mat Omega;
  arma::vec rf;
  
  double scale_post = 0.0;
  double nu_post = 0.0;
  
  for(size_t i = 0; i < di.N; i++){
    nu_post += di.n[i];
    rf.set_size(di.n[i]);
    Sigma.set_size(di.n[i], di.n[i]);
    
    for(size_t j = 0; j < di.n[i]; j++){
      for(size_t jj = j; jj < di.n[i]; jj++){
        Sigma(j,jj) = pow(rho, jj-j);
        Sigma(jj,j) = pow(rho, jj-j);
      }
      rf(j) = di.rf[j + di.start_index[i]];
    }
    Omega = arma::inv_sympd(Sigma);
    scale_post += arma::as_scalar(rf.t() * Omega * rf);
  }
  
  //Rcpp::Rcout << "[ar_mat]: nu_post = " << nu_post << "  scale_post = " << scale_post << std::endl;
  
  double sigma_prop = sqrt( (scale_post)/gen.chi_square(nu_post)); // this is the proposed sigma'
  
  double log_prior_prop = -0.5 * (1.0 + sigma_pi.nu) * log(1.0 + 1.0/sigma_pi.nu * sigma_prop * sigma_prop/(sigma_pi.A * sigma_pi.A));
  double log_prior = -0.5 * (1.0 + sigma_pi.nu) * log(1.0 + 1.0/sigma_pi.nu * sigma * sigma/(sigma_pi.A * sigma_pi.A));
  double alpha = exp(3.0 * log(sigma_prop) + log_prior_prop - 3.0 * log(sigma) - log_prior);
  if(gen.uniform() <= alpha) sigma = sigma_prop;
  
}



void update_sigma_mu(double* sigma_mu, std::vector<std::vector<tree > > t_vec, tree_prior_info &tree_pi, RNG &gen)
{
  tree::npv bnv; // holds pointers to bottom nodes of the trees
  size_t p = t_vec.size();
  size_t L_tot = 0;
  
  double M_tot = 0.0;
  double alpha = 0.0;
  double sigma_prop = 0.0;
  
  double log_prior_prop = 0.0;
  double log_prior = 0.0;
  
  for(size_t k = 0; k < p; k++){
    L_tot = 0;
    M_tot = 0.0;
    
    for(size_t m = 0; m < t_vec[k].size(); m++){
      bnv.clear();
      t_vec[k][m].getbots(bnv); // get bottom nodes for tree m for beta_k
      for(size_t l = 0; l < bnv.size(); l++){
        L_tot++;
        M_tot += pow(bnv[l]->getm(),2);
      }
      sigma_prop = sqrt( (M_tot)/gen.chi_square(L_tot));
      // this is for half-t
      log_prior_prop = -0.5 * (1.0 + tree_pi.nu[k]) * log(1.0 + 1.0/tree_pi.nu[k] * sigma_prop * sigma_prop/(tree_pi.A[k] * tree_pi.A[k]));
      log_prior = -0.5 * (1.0 + tree_pi.nu[k]) * log(1.0 + 1.0/tree_pi.nu[k] * sigma_mu[k] * sigma_mu[k]/(tree_pi.A[k] * tree_pi.A[k]));
      alpha = exp(3.0 * log(sigma_prop) + log_prior_prop - 3.0 * log(sigma_mu[k]) - log_prior);
      if(gen.uniform() <= alpha) sigma_mu[k] = sigma_prop;
    }
  }
}
