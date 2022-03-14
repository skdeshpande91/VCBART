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


// sufficient statistic map:
//   key: node id
//   value: vector of length n; element i corresponds to subject i and contains vector of the indices for that subject that land in the corresponding bottom node
typedef std::map<int, std::vector<std::vector<int>> suff_stat; //
typedef suff_stat::iterator suff_stat_it; // iterator for sufficient statistic map

// splitting rule: left child if x[v] < c and right child if x[v] >= c


class rule_t{
public:
  bool is_aa; // is it an axis-aligned split?
  bool is_cat; // is it a categorical split?
  int v_aa; // index of continuous variable on which we split
  double c; // cutpoint
  
  int v_cat; // index of categorical variable on which we split (always between 0 and R_cat)
  std::set<int> l_vals; // holds unique levels of v_cat associated with left child
  std::set<int> r_vals; // holds unique levels of v_cat associated with right child
  
  rule_t(){
    is_aa = true;
    is_cat = false;
    v_aa = 0;
    c = 0.0;
    v_cat = 0;
    l_vals = std::set<int>();
    r_vals = std::set<int>();
  }
  void clear(){
    is_aa = true;
    is_cat = false;
    v_aa = 0;
    c = 0.0;
    v_cat = 0;
    l_vals = std::set<int>();
    r_vals = std::set<int>();
  }
  void copy(rule_t &old_rule){
    is_aa = old_rule.is_aa;
    is_cat = old_rule.is_cat;
    v_aa = old_rule.v_aa;
    c = old_rule.c;
    v_cat = old_rule.v_cat;
    l_vals = old_rule.l_vals;
    r_vals = old_rule.r_vals;
  };
};

// class holding data dimensions, pointers to the covariate data, etc
class data_info{
public:
  int N; // total number of observations
  int n; // total number of subjects
  int p; // total number of covariates
  int R; // total number of modifiers
  int R_cont; // total number of continuous modifiers
  int R_cat; // total number of categorical modifiers
    
  int* obs_id; // id[i] tells us to which subject observation i is associated
  
  double* x; // pointer to matrix of covariates
  
  double* z_cont; // pointer to matrix of continuous effect modifiers
  bool* unif_cuts; // pointer to collection of booleans that dictate whether we do uniform cuts or if we use pre-specified cutpoints
  std::vector<std::set<double> >* cutpoints;
  
  int* z_cat; // pointer to matrix of categorical effect modifiers
  std::vector<int>* K; // number of levels of each categorical variable
  std::vector<std::set<int> >* cat_levels; // holds unique values of each categorical variable (just 0,...., K-1)
  std::vector<std::vector<unsigned int> >* adj_support; // support of lower triangle of adjacency matrix for categorical levels
  double* rp; // partial residual
 
  int k; // keeps track of which beta function we are updating at any time
  data_info(){
    N = 0;
    n = 0;
    p = 0;
    R = 0;
    R_cont = 0;
    R_cat = 0;
    obs_id = 0; // 0 pointer
    x = 0; // 0 pointer
    z_cont = 0; // 0 pointer
    unif_cuts = 0; // 0 pointer
    cutpoints = 0; // 0 pointer
    z_cat = 0; // 0 pointer
    K = 0; // 0 pointer
    cat_levels = 0; // 0 pointer
    adj_support = 0; // 0 pointer
    rp = 0; // 0 pointer
  }
};

class tree_prior_info{
public:
  double alpha; // 1st parameter of branching process prior
  double beta; // 2nd parameter of branching process prior
  
  double prob_bd; // prob of proposing a grow (birth) or prune (death) move. almost always set to 1
  double prob_b; // prob of proposing a grow (birth) move, conditional on doing grow or prune. almost always set  to 0.5
  
  bool mst_split; // do we split categorical variable using a random MST?
  
  std::vector<int>* var_count; // counts how many times we have used each variable in a splitting rule
  
  std::vector<double>* theta; // probabilities of selecting a predictor in the decision rule
  
  double tau; // prior sd for leaf
  double mu0; // prior mean for leaf
  
  // constructor
  tree_prior_info(){
    alpha = 0.95;
    beta = 2.0;
    prob_bd = 1.0;
    prob_b = 0.5;
    mst_split = false;
    var_count = 0; // 0 pointer
    theta = 0; // 0 pointer
    tau = 1.0;
    mu = 0.0;
  }
  
};


// processes the inputted information about categorical predictors (if they exist)

inline void parse_cat_levels(std::vector<std::set<int> > &cat_levels, std::vector<int> &K, int R_cat, Rcpp::List &tmp_cat_levels)
{
  cat_levels.clear();
  K.clear();
  if(tmp_cat_levels.size() == R_cat){
    for(int j = 0; j < R_cat; j++){
      Rcpp::IntegerVector levels_vec = Rcpp::as<Rcpp::IntegerVector>(tmp_cat_levels[j]);
      std::set<int> levels_set;
      for(int l = 0; l < levels_vec.size(); l++) levels_set.insert(levels_vec[l]);
      cat_levels.push_back(levels_set);
      K.push_back(levels_set.size());
    }
  } else{
    Rcpp::Rcout << "R_cat = " << R_cat;
    Rcpp::Rcout << "cat_levels_list.size() = " << tmp_cat_levels.size();
    Rcpp::stop("cat_levels_list must have size equal to R_cat!");
  }
}

inline void parse_cat_adj(std::vector<std::vector<unsigned int>> &adj_support, int R_cat, Rcpp::List &tmp_adj_support)
{
  adj_support.clear();
  if(tmp_adj_support.size() == R_cat){
    Rcpp::IntegerVector adj_rvec = Rcpp::as<Rcpp::IntegerVector>(tmp_adj_support[j]);
    std::vector<unsisnged int> adj_uvec;
    for(int l = 0; l < adj_rvec.size(); l++) adj_uvec.push_back( (unsigned int) adj_rvec[l]);
    adj_support.push_back(adj_uvec);
  } else{
    Rcpp::Rcout << "R_cat = " << R_cat;
    Rcpp::Rcout << "adj_levels_list.size() = " << tmp_adj_support.size() << std::endl;
    Rcpp::stop("adj_levels_list must have size equal to R_cat!");
  }
}

inline void parse_cutpoints(std::vector<std::set<double>> &cutpoints, int R_cont, Rcpp::List &tmp_cutpoints)
{
  cutpoints.clear();
  if(tmp_cutpoints.size() == R_cont){
    for(int j = 0; j < R_cont; j++){
      Rcpp::NumericVector cutpoints_vec = Rcpp::as<Rcpp::NumericVector>(tmp_cutpoints[j]);
      std::set<double> xi_set;
      for(int l = 0; l < cutpoints_vec.size(); l++) xi_set.insert(cutpoints_vec[l]);
      cutpoints.push_back(xi_set);
    }
  } else{
    Rcpp::Rcout << "R_cont = " << R_cont;
    Rcpp::Rcout << "  cutpoints_list.size() = " << tmp_cutpoints.size() << std::endl;
    Rcpp::stop("cutpoints_list needs to have length R_cont!");
  }
}

inline void parse_training_data(int &N_train, int &p, int &R, int &R_cont, int &R_cat, Rcpp::NumericMatrix &tX_train, Rcpp::NumericMatrix &tZ_cont_train, Rcpp::IntegerMatrix &tZ_cat_train)
{
  N_train = tX_train.cols(); // N_train is nrow(X_train) so it is ncol(t(X_train)). this is total number of observations
  p = tX_train.rows();
  if(tZ_cont_train.size() > 1 && tZ_cat_train.size() == 1){
    // only continuous modifiers are available
    if(tZ_cont_train.cols() == N_train){
      R_cont = tZ_cont_train.rows();
      R_cat = 0;
    } else{
      Rcpp::Rcout << "Z_cont_train has " << tZ_cont_train.cols() << " rows" << std::endl;
      Rcpp::Rcout << "X_train has " << N_train << " rows" << std::endl;
      Rcpp::stop("Z_cont_train and X_train must have same number of rows!");
    }
  } else if(tZ_cont_train.size() == 1 && tZ_cat_train.size() > 1){
    // only categorical predictors are available
    if(tZ_cat_train.cols() == N_train){
      R_cont = 0;
      R_cat = tZ_cat_train.rows();
    } else{
      Rcpp::Rcout << "Z_cat_train has " << tZ_cat_train.cols() << " rows" << std::endl;
      Rcpp::Rcout << "X_train has " << N_train << " rows" << std::endl;
      Rcpp::stop("Z_cat_train and X_train must have same number of rows!");
    }
  } else if(tZ_cont_train.size() > 1 && tZ_cat_train.size() > 1){
    if(tZ_cont_train.cols() == N_train && tZ_cat_train.cols() == N_train){
      R_cont = tZ_cont_train.rows();
      R_cat = tZ_cat_train.rows();
    } else{
      Rcpp::Rcout << "Z_cont_train has " << tZ_cont_train.cols() << " rows" << std::endl;
      Rcpp::Rcout << "Z_cat_train has " << tZ_cat_train.cols() << " rows" << std::endl;
      Rcpp::Rcout << "X_train has " << N_train << " rows" << std::endl;
      Rcpp::stop("Z_cont_train, Z_cat_train, and X_train must have same number of rows!");
    }
  } else{
    Rcpp::stop("No continuous or categorical modifiers are available!");
  }
  R = R_cont + R_cat;

}

inline void parse_testing_data(int &N_test, int p, int R_cont, int R_cat, Rcpp::NumericMatrix &tX_test, Rcpp::NumericMatrix &tZ_cont_test, Rcpp::IntegerMatrix &tZ_cat_test)
{
  N_test = tX_test.cols(); // N_train is nrow(X_train) so it is ncol(t(X_train)). this is total number of observations
  if(tX_test.rows() != p){
    Rcpp::Rcout << "X_train has " << p << " columns";
    Rcpp::Rcout << "  X_test has " << tX_test.rows() << " columns" << std::endl;
    Rcpp::stop("must have same number of covariates in training and testing data!");
  }
  if(tZ_cont_test.size() > 1){
    if(tZ_cont_test.rows() != R_cont){
      Rcpp::Rcout << "Z_cont_train has " << R_cont << " columns";
      Rcpp::Rcout << "Z_cont_test has " << tZ_cont_test.rows() << " columns" << std::endl;
      Rcpp::Rcout << "must have same number of continuous modifiers in training and testing data!");
    } else if(tZ_cont_test.cols() != N_test){
      Rcpp::Rcout << "Z_cont_test has " << tZ_cont_test.cols() << " rows" << std::endl;
      Rcpp::Rcout << "X_test has " << N_test << " rows" << std::endl;
      Rcpp::stop("Z_cont_test and X_test must have same number of rows!");
    }
  }
  
  if(tZ_cat_test.size() > 1){
    if(tZ_cat_test.rows() != R_cat){
      Rcpp::Rcout << "Z_cat_train has " << R_cat << " columns";
      Rcpp::Rcout << "Z_cat_test has " << tZ_cat_test.rows() << " columns" << std::endl;
      Rcpp::Rcout << "must have same number of categorical modifiers in training and testing data!");
    } else if(tZ_cat_test.cols() != N_test){
      Rcpp::Rcout << "Z_cat_test has " << tZ_cat_test.cols() << " rows" << std::endl;
      Rcpp::Rcout << "X_test has " << N_test << " rows" << std::endl;
      Rcpp::stop("Z_cat_test and X_test must have same number of rows!");
    }
  }
}




#endif






