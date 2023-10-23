#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat rescale_beta_mean(arma::mat beta_input,
                            double y_mean,
                            double y_sd,
                            Rcpp::NumericVector x_mean,
                            Rcpp::NumericVector x_sd)
{
  int p = beta_input.n_cols;
  int nd = beta_input.n_rows;
  arma::mat beta_out = arma::zeros<arma::mat>(nd,p);
  
  
  beta_out.col(0) = beta_input.col(0) * y_sd;
  beta_out.col(0) += y_mean;
  
  for(int j = 1; j < p; j++){
    beta_out.col(j) = y_sd/x_sd[j] * beta_input.col(j);
    beta_out.col(0) -= y_sd/x_sd[j] * x_mean[j] * beta_input.col(j);
  }
  
  return beta_out;
}

// [[Rcpp::export]]
arma::cube rescale_beta(arma::cube beta_input,
                        double y_mean,
                        double y_sd,
                        Rcpp::NumericVector x_mean,
                        Rcpp::NumericVector x_sd)
{
  // beta_input is nd x n x p
  int p = beta_input.n_slices; // number of functions
  int n = beta_input.n_cols; // number of individuals
  int nd = beta_input.n_rows; // number of samples
  arma::cube beta_out = arma::zeros<arma::cube>(nd,n,p);
  
  
  beta_out.slice(0) = beta_input.slice(0) * y_sd;
  beta_out.slice(0) += y_mean;
  
  for(int j = 1; j < p; j++){
    beta_out.slice(j) = y_sd/x_sd[j] * beta_input.slice(j);
    beta_out.slice(0) -= y_sd/x_sd[j] * x_mean[j] * beta_input.slice(j);
  }
  return beta_out;
}
