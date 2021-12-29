#ifndef GUARD_helper_funs_h
#define GUARD_helper_funs_h

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstddef>
#include <vector>


typedef std::vector<int>::iterator int_it; // iterator type for vectors of integers
typedef std::vector<double>::iterator dbl_it; // iterator type vectors of doubles
typedef std::vector<std::vector<int> > vv_int; // vector of vector of integers
typedef std::map<int, vvint> suff_stat;
typedef suff_stat::itetator suff_stat_it; // iterator for sufficient statistic map
// key: node id
// value is a vector of length n; element i corresponds to subject i and contains vector of all the indices that land in that bottom node

// splitting rule: left child if x[v] < c and right child if x[v] >= c


class rule_t{
public:
  bool binary; // if we're splitting on binary value, we just set c = 1 and move on
  int v; // index of the variable that we split on
  double c; // cutpoint
  rule_t(){binary = false; v = 0; c = 0;}
  void clear(){binary = false; v = 0; c = 0;}
}

// class holding data dimensions, pointers to the covariate data, etc
class data_info{
public:
  int N; // total number of observations
  int n; // total number of subjects
  int p; // total number of covariates
  int R; // total number of modifiers
    
  int* n_vec; // pointer to first element of n_vec
  int* start_index; // pointer to first element of start_index
  
  double* y; // outcome for subject i, observation j is *(y + (start_index[i] + j));
  double* x; // j-th observation of X_k individual i is *(x + (k  + (start_index[i] + j)*p));
  double* z; // j-th observation of Z_r individual i is *(z + (r + (start_index[i] + j)*R));
  double* rp; // partial residual
  int k; // keeps track of which beta function we are updating at any time
  data_info(){N = 0; n = 0; p = 0; R = 0; n_vec = 0; start_index = 0; rp = 0; k = 0;}
}

class tree_prior_info{
public:
  double pbd; // prob. of doing a grow or prune. set to 1 typically
  double pb; // prob. of a grow (usually 0.5)
  double alpha;
  double beta;
  
  std::vector<std::vector<int> >* var_count; // counts how many time variable used in each ensemble
  std::vector<std::vector<double>* theta; // split probabilities for each ensemble
  double* tau; // points to the NumericVector holding the values of tau
  
  tree_prior_into(){
    pbd = 1.0;
    pb = 0.5;
    alpha = 0.95;
    beta = 2.0;
    var_count = 0;
    theta = 0;
    tau = 0;
  }
}

// just to improve readability of main function

inline void parse_input(int &N_train, int &n_train, int &p, int &R, int &N_test, int &n_test,
                         Rcpp::NumericVector &Y_train,
                         Rcpp::NumericMatrix &tX_train, Rcpp::NumericMatrix &tZ_train, Rcpp::IntegerVector n_vec_train,
                         Rcpp::NumericMatrix &tX_test, Rcpp::NumericMatrix &tZ_test, Rcpp::IntegerVector n_vec_test)
{
  N_train = Y_train.length();
  n_train = n_vec_train.length();
  p = tX_train.rows(); // nrow(t(X_train)) = ncol(X_train)
  R = tZ_train.rows(); // nrow(t(Z_train)) = ncol(Z_train)
  
  if(tX_test.size() == 1){
    N_test = 0;
    n_test = 0;
  } else{
    N_test = tX_test.cols(); // ncol(t(X_test)) = nrow(X_test)
    n_test = n_vec_test.length();
  }
}

#endif






