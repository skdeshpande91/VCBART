//
//  sparse_vcbart_ind:
//    Indepednent errors and adaptive split probabilities
//

# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <vector>
#include <ctime>
#include <algorithm>
#include <stdio.h>

#include "rng.h"
#include "tree.h"
#include "funs.h"
#include "mu_posterior.h"
#include "update_scales.h"
#include "update_tree.h"
#include "update_split_probs.h"
#include "update_alpha_z.h"

// [[Rcpp::export(name = ".vcbart_ind_adapt_split")]]
Rcpp::List vcbart_ind_adapt_split(arma::vec Y, // n_train x 1 ... concatenation of all observed Y's
                                  arma::mat X_train, // n_obs x p  ... stack of all observed covariates x for training
                                  arma::mat Z_train, // n_obs x R ... stack of all observed modifiers z for training
                                  arma::vec n_vec_train, // number of observations per individual (training)
                                  arma::vec start_index_vec_train, // start_index_vec_train(i) tells us where individual i's observations start (training)
                                  arma::mat X_test, // n_test x p ... stack of all covariates x for testing
                                  arma::mat Z_test, // n_test x R ... stack of all modifiers z for testing
                                  arma::vec n_vec_test, // number of observations per individual (testing)
                                  arma::vec start_index_vec_test, // start_index_vec_test(i) tells us where individual i's observation start (testing)
                                  Rcpp::List xinfo_list, // cutpoints for z
                                  size_t M, // number of trees
                                  bool ht_sigma_y, bool ht_tau, // whether to use half-t priors for sigma & tau
                                  size_t burn, size_t nd, // number of iterations to burn-in and save
                                  bool verbose, size_t print_every, // print progress?
                                  double a, double b, size_t N_u, double rho_alpha, // hyper-parameters for split_probs
                                  arma::vec tau_vec, // leaf variances
                                  arma::vec alpha_vec, // alpha for tree prior
                                  arma::vec beta_vec, // beta for tree prior
                                  double sigma_hat, double nu_sigma, double nu_tau, double variance_prob) // arguments for CGM prior elicitation

{
  Rcpp::RNGScope scope;
  RNG gen;
  
  
  size_t n_obs_train = X_train.n_rows;
  size_t n_obs_test = X_test.n_rows;
  
  size_t N_train = n_vec_train.size(); // total number of individuals in training set
  size_t N_test = n_vec_test.size(); // total number of individuals in testing set
  size_t p = X_train.n_cols; // number of columns in X... for intercepts, assume first column of X is all 1's
  size_t R = Z_train.n_cols; // number of modifiers
  
  if( (arma::all(X_train.col(0) == 1.0)) != (arma::all(X_test.col(0) == 1.0)) ) Rcpp::stop("Only one of training and testing X includes intercept");
  
  bool intercept = arma::all(X_train.col(0) == 1.0);
  
  /*
  if(verbose == true){
    Rcpp::Rcout << " n_train  = " << N_train << "  n_test = " << N_test;
    Rcpp::Rcout << " n_obs_train = " << n_obs_train << "  n_obs_test = " << n_obs_test;
    Rcpp::Rcout << " p = " << p; // includes intercept
    Rcpp::Rcout << " R = " << R;
    Rcpp::Rcout << std::endl;
  }
  */
  std::vector<size_t> n_train(N_train);
  std::vector<size_t> n_test(N_test);
  
  std::vector<size_t> start_index_train(N_train);
  std::vector<size_t> start_index_test(N_test);
  
  
  for(size_t i = 0; i < N_train; i++){
    n_train[i] = n_vec_train(i);
    start_index_train[i] = start_index_vec_train(i)-1; // remember R is 1-indexed!
  }
  for(size_t i = 0; i < N_test; i++){
    n_test[i] = n_vec_test(i);
    start_index_test[i] = start_index_vec_test(i) - 1; // remember R is 1-indexed!
  }
  // Rcpp::Rcout << "  Got indices for each individual's start and end" << std::endl;
  
  // need to center and scale Y.
  double y_mean = 0.0;
  double y_sd = 0.0;
  double y_max = 0.0;
  double y_min = 0.0;
  
  prepare_y(Y, y_mean, y_sd, y_max, y_min);
  //Rcpp::Rcout << "Centered and scaled Y!" << endl;
  
  // we also need to center and scale columns of X_train & X_test!
  std::vector<double> x_col_mean(p);
  std::vector<double> x_col_sd(p);
  prepare_x(X_train, X_test, x_col_mean, x_col_sd);
  
  //if(verbose == true) Rcpp::Rcout << "  intercept = " << intercept << std::endl;
  //Rcpp::Rcout << "Centered and scaled X!" << endl;
  
  // need to create pointers for the data
  double* y_ptr = new double[n_obs_train]; // pointer to training set Y's
  double* x_ptr = new double[n_obs_train * p]; // pointer to training set X's
  double* z_ptr = new double[n_obs_train * R]; // pointer to training set Z's
  
  double* x_pred_ptr = new double[n_obs_test * p]; // pointer to test set X's
  double* z_pred_ptr = new double[n_obs_test * R]; // pointer to test set Z's
  
  for(size_t i = 0; i < N_train; i++){
    for(size_t j = 0; j < n_train[i]; j++){
      y_ptr[start_index_train[i] + j] = Y(start_index_train[i] + j); // Y is n_obs x 1
      for(size_t k = 0; k < p; k++) x_ptr[k + (start_index_train[i] + j)*p] = X_train(start_index_train[i] + j, k); // vectorize X_train (n_obs_train x p) by row
      for(size_t r = 0; r < R; r++) z_ptr[r + (start_index_train[i] + j)*R] = Z_train(start_index_train[i] + j, r); // vectorize Z_train (n_obs_train x R) by row
    }
  }
  
  for(size_t i = 0; i < N_test; i++){
    for(size_t j = 0; j < n_test[i]; j++){
      for(size_t k = 0; k < p; k++) x_pred_ptr[k + (start_index_test[i] + j) * p] = X_test(start_index_test[i] + j, k); // vectorize X_test (n_obs_test x p) by row
      for(size_t r = 0; r < R; r++) z_pred_ptr[r + (start_index_test[i] + j) * R] = Z_test(start_index_test[i] + j, r); // vectorize Z_test(n_obs_test x R) by row
    }
  }
  
  // in fit, allsuff, bn, etc. we point to the start of an observed vector of covariates and then increment by p to get to the next individual
  
  
  //Rcpp::Rcout << "Created y_ptr, x_ptr, z_ptr!" << endl;
  // Read in and format the cutpoints
  xinfo xi;
  xi.resize(R);
  for(size_t r = 0; r < R; r++){
    Rcpp::NumericVector tmp = xinfo_list[r];
    std::vector<double> tmp2;
    for(int rr = 0; rr < tmp.size(); rr++) tmp2.push_back(tmp[rr]);
    xi[r] = tmp2;
  }
  //Rcpp::Rcout << "Created xi!" << endl;
  
  
  // trees and points for fits and residuals
  std::vector<std::vector<tree > > t_vec(p, std::vector<tree>(M)); // a vector of vector of trees!
  std::vector<std::vector<double> > theta_vec(p, std::vector<double>(R, 1.0/R));
  std::vector<std::vector<size_t> > var_counts(p, std::vector<size_t>(R,0)); // var_counts[k][r] counts #times we split on Z_r in beta_k ensemble
  
  // some prior stuff for alpha
  std::vector<double> alpha_z(p);
  double tmp_alpha = 0.0;
  double u_init = 0.0;
  for(size_t k = 0; k < p; k++){
    u_init = gen.beta(a,b);
    alpha_z[k] = rho_alpha * u_init/(1.0 - u_init);
  }

  double* allfit = new double[n_obs_train]; //
  double* ftemp = new double[n_obs_train]; // holds the temporary fit from single tree on training set
  double* ftemp_pred = new double[n_obs_test]; // holds temporary fit from single tree on testing set
  double* betafit = new double[n_obs_train * p]; // holds the values of the beta functions evaluated on training set
  double* betafit_pred = new double[n_obs_test * p]; // holds the values of the beta functions evaluate on testing set
  double* r_full = new double[n_obs_train]; // r_full = y_ptr - allfit
  double* r_partial = new double[n_obs_train]; // this will hold partial residual and is ultimately what is passed to bd()
  
  // set the initial value of the tree fits to 0
  for(size_t k = 0; k < p; k++){
    for(size_t t = 0; t < M; t++) t_vec[k][t].setm(0.0);
    for(size_t i = 0; i < N_train; i++){
      for(size_t j = 0; j < n_train[i]; j++) betafit[k + p*(start_index_train[i] + j)] = 0.0;
    }
    
    for(size_t i = 0; i < N_test; i++){
      for(size_t j = 0; j < n_test[i]; j++) betafit_pred[k + p*(start_index_test[i] + j)] = 0.0;
    }
  }
  
  for(size_t i = 0; i < N_train; i++){
    for(size_t j = 0; j < n_train[i]; j++){
      ftemp[start_index_train[i] + j] = 0.0;
      allfit[start_index_train[i] + j] = 0.0;
      r_full[start_index_train[i] + j] = y_ptr[start_index_train[i] + j] - allfit[start_index_train[i] + j];
      r_partial[start_index_train[i] + j] = 0.0; // initial value doesn't matter much
    }
  }
  
  for(size_t i = 0; i < N_test; i++){
    for(size_t j = 0; j < n_test[i]; j++){
      ftemp_pred[start_index_test[i] + j] = 0.0;
    }
  }
  
  //Rcpp::Rcout << "Created allfit, ftemp, betafit, r_full!" << endl;
  
  // make the data info objects
  data_info di;
  di.N = N_train;
  di.p = p;
  di.R = R;
  di.n = n_train;
  di.start_index = start_index_train;
  di.y = &y_ptr[0];
  di.x = &x_ptr[0];
  di.z = &z_ptr[0];
  di.af = &allfit[0];
  di.rf = &r_full[0];
  di.k = 0;
  
  data_info dip;
  dip.N = N_test;
  dip.p = p;
  dip.R = R;
  dip.n = n_test;
  dip.start_index = start_index_test;
  dip.x = &x_pred_ptr[0];
  dip.z = &z_pred_ptr[0];
  dip.k = 0;
  
  
  //Rcpp::Rcout << "Created the data info object" << endl;
  
  tree_prior_info b_tree_pi;
  b_tree_pi.pbd = 1.0;
  b_tree_pi.pb = 0.5;
  b_tree_pi.alpha.resize(p, 0.95);
  b_tree_pi.beta.resize(p, 2.0);
  
  for(size_t k = 0; k < p; k++){
    b_tree_pi.alpha[k] = alpha_vec(k);
    b_tree_pi.beta[k] = beta_vec(k);
  }
  
  b_tree_pi.rp = &r_partial[0];
  b_tree_pi.A.clear();
  b_tree_pi.A.resize(p, 0.25/sqrt( (double) M)); // this is what Linero & Yang use.
  b_tree_pi.nu.clear();
  b_tree_pi.nu.resize(p, nu_tau) ; // use a half-t_7 distribution centered at A
  
  double* tau = new double[p]; //
  for(size_t k = 0; k < p; k++) tau[k] = tau_vec(k);
  b_tree_pi.tau = &tau[0];
  
  //Rcpp::Rcout << "Created tree_prior objects" << endl;
  
  // create sigma_prior objects
  // Initialize sigma and the precision matrices...
  
  double chisq_quantile = 0.0;
  Rcpp::Function qchisq("qchisq"); // uses R's built-in quantile function for chi-square
  Rcpp::NumericVector tmp_quantile = qchisq(Rcpp::Named("p") = 1.0 - variance_prob, Rcpp::Named("df") = nu_sigma);
  chisq_quantile = tmp_quantile[0];
  
  double sigma = 1.0; // global sigma
  sigma_prior_info sigma_pi; //
  
  sigma_pi.sigma_hat = sigma_hat; // default argument sigma_hat = 1.0
  sigma_pi.nu = nu_sigma;
  sigma_pi.lambda = chisq_quantile/nu_sigma * sigma_hat;
  sigma_pi.A = sigma_hat; // I think this should be equal to sigma_hat.
  
  // containers for output
  arma::mat f_train_samples(n_obs_train, nd); // training fits
  arma::mat f_test_samples(n_obs_test, nd); // testing fits
  arma::cube beta_train_samples(n_obs_train, p, nd); // training fits
  arma::cube beta_test_samples(n_obs_test, p, nd); // testing fits
  
  arma::vec sigma_samples(nd+burn);
  arma::cube theta_samples(R, p, nd+burn);
  arma::mat alpha_samples(p, nd+burn);
  arma::cube var_counts_samples(R, p, nd+burn);
  double tmp_int = 0.0; // holds value of the offset
  double tmp_fit = 0.0; // holds individual prediction for E[y_ij]
  
  if(verbose == true) Rcpp::Rcout << "[VCBART]: Entering MCMC" << std::endl;
  clock_t start_time = clock();
  
  for(size_t iter = 0; iter < nd+burn; iter++){
    
    if(verbose == true){
      if( (iter < burn) && (iter%print_every == 0)) Rcpp::Rcout << "  MCMC Iteration: " << iter << " of " << nd + burn << "; Burn-in" << std::endl;
      else if( ( (iter > burn) && (iter%print_every == 0) ) || (iter == burn)) Rcpp::Rcout << "  MCMC Iteration: " << iter << " of " << nd + burn << "; Sampling" << std::endl;
    }
    
    // update the trees
    //Rcpp::Rcout << "Before tree update r_full[0] = " << r_full[0] << std::endl;
    for(size_t k = 0; k < p; k++){
      di.k = k; // keep track of which set of beta functions we're updating
      
      for(size_t m = 0; m < M; m++){
        fit(t_vec[k][m], xi, di, ftemp); // get fit of tree
        for(size_t i = 0; i < N_train; i++){
          for(size_t j = 0; j < n_train[i]; j++){
            allfit[start_index_train[i] + j] -= x_ptr[k + p*(start_index_train[i] + j)] * ftemp[start_index_train[i] + j]; // remove contribution of single tree
            betafit[k + p*(start_index_train[i] + j)] -= ftemp[start_index_train[i] + j]; // remove fit of single  tree from beta estimate
            r_partial[start_index_train[i] + j] = y_ptr[start_index_train[i] + j] - allfit[start_index_train[i] + j]; // update r_partial
          }
        } // closes loop that removes fit of single tre from allfit and r_partial
        update_tree_ind(t_vec[k][m], sigma, theta_vec[k], var_counts[k], xi, di, b_tree_pi, gen); // update the tree
        fit(t_vec[k][m], xi, di, ftemp);
        for(size_t i = 0; i < N_train; i++){
          for(size_t j = 0; j < n_train[i]; j++){
            allfit[start_index_train[i] + j] += x_ptr[k + p*(start_index_train[i] + j)] * ftemp[start_index_train[i] + j];
            betafit[k + p*(start_index_train[i] + j)] += ftemp[start_index_train[i] + j];
          }
        }
      } // closes loop over trees for beta_k
      
      // update the split probabilities for beta_k.
      update_split_probs(theta_vec[k], var_counts[k], alpha_z[k], R, gen);
      
      // update the alpha parameter for beta_k:
      update_alpha_z(tmp_alpha, rho_alpha, theta_vec[k], R, N_u, a, b, gen);
      alpha_z[k] = tmp_alpha;
      alpha_samples(k, iter) = alpha_z[k];
      // save theta_samples and also save var_counts
      for(size_t r = 0; r < R; r++){
        theta_samples(r, k, iter) = theta_vec[k][r];
        var_counts_samples(r,k,iter) = var_counts[k][r];
      }
    } // closes loop over k
    
    // now that all trees have been updated, update the full residuals
    for(size_t i = 0; i < N_train; i++){
      for(size_t j = 0; j < n_train[i]; j++){
        r_full[start_index_train[i] + j] = y_ptr[start_index_train[i] + j] - allfit[start_index_train[i] + j];
      }
    }

    // update sigma
    if(ht_sigma_y == true) update_sigma_ht_ind(sigma, sigma_pi, di, gen);
    else update_sigma_ig_ind(sigma, sigma_pi, di, gen);
    sigma_samples(iter) = sigma * y_sd;
    
    
    
    
    if(iter >= burn){
      for(size_t i = 0; i < N_train; i++){
        for(size_t j = 0; j < n_train[i]; j++){
          tmp_fit = 0.0;
          tmp_int = y_mean;
          for(size_t k = 0; k < p; k++){
            beta_train_samples(start_index_train[i] + j, k, iter-burn) = betafit[k + p*(start_index_train[i] + j)] * y_sd/x_col_sd[k];
            tmp_fit += x_ptr[k + p*(start_index_train[i]  + j)] * betafit[k + p*(start_index_train[i] + j)]; // build up value of E[y | x,z]
            tmp_int -= x_col_mean[k] * betafit[k + p*(start_index_train[i] + j)] * y_sd/x_col_sd[k]; // this is the re-centering we need to do to all predictions
          }
          if(intercept == true) beta_train_samples(start_index_train[i] + j, 0, iter-burn) += tmp_int;
          f_train_samples(start_index_train[i] + j, iter-burn) = y_sd * tmp_fit + y_mean;
        }
      }
      
      for(size_t k = 0; k < p; k++){
        dip.k = k;
        for(size_t i = 0; i < N_test; i++){
          for(size_t j = 0; j < n_test[i]; j++) betafit_pred[k + p*(start_index_test[i] + j)] = 0.0; // reset value of betafit_pred
        }
        
        for(size_t m = 0; m < M; m++){
          fit(t_vec[k][m], xi, dip, ftemp_pred);
          for(size_t i = 0; i < N_test; i++){
            for(size_t j = 0; j < n_test[i]; j++) betafit_pred[k + p*(start_index_test[i] + j)] += ftemp_pred[start_index_test[i] + j];
          }
        }
      }
      for(size_t i = 0; i < N_test; i++){
        for(size_t j = 0; j < n_test[i]; j++){
          tmp_fit = 0.0;
          tmp_int = y_mean;
          for(size_t k = 0; k < p; k++){
            beta_test_samples(start_index_test[i] + j, k, iter-burn) = betafit_pred[k + p*(start_index_test[i] + j)] * y_sd/x_col_sd[k];
            tmp_fit += x_pred_ptr[k + p*(start_index_test[i] + j)] * betafit_pred[k + p*(start_index_test[i] + j)];
            tmp_int -= x_col_mean[k] * betafit_pred[k + p*(start_index_test[i] + j)] * y_sd/x_col_sd[k];
          }
          if(intercept == true) beta_test_samples(start_index_test[i] + j, 0, iter-burn) += tmp_int;
          f_test_samples(start_index_test[i] + j, iter-burn) = y_sd * tmp_fit + y_mean;
        }
      }
    } // closes if checking that we need to save samples
    
  } // closes main MCMC loop
  clock_t end_time = clock();
  long double sampler_time = ( (long double) (end_time - start_time))/CLOCKS_PER_SEC;
  
  if(verbose == true) Rcpp::Rcout << "[VCBART]: Finished MCMC! time = " << sampler_time << std::endl;;
  
  
  delete[] y_ptr;
  delete[] x_ptr;
  delete[] x_pred_ptr;
  delete[] z_ptr;
  delete[] z_pred_ptr;
  delete[] allfit;
  delete[] ftemp;
  delete[] ftemp_pred;
  delete[] betafit;
  delete[] betafit_pred;
  delete[] r_full;
  delete[] r_partial;
  
  Rcpp::List results;
  results["f_train_samples"] = f_train_samples;
  results["f_test_samples"] = f_test_samples;
  results["beta_train_samples"] = beta_train_samples;
  results["beta_test_samples"] = beta_test_samples;
  results["sigma_samples"] = sigma_samples;
  results["theta_samples"] = theta_samples;
  //results["alpha_samples"] = alpha_samples;
  results["var_counts_samples"] = var_counts_samples;
  results["time"] = sampler_time;
  return(results);
}

