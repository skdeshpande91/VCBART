#ifndef GUARD_info_h
#define GUARD_info_h

#include <RcppArmadillo.h>
using namespace arma;

//============================================================
//data
//============================================================

class data_info{
public:
  size_t N; // number of individuals
  size_t p; // number of predictors
  size_t R; // number of modifiers --> everywhere that we used to use p, we need to use r
  std::vector<size_t> n; // vector of length N, tells us how many observations for individual i
  std::vector<size_t> start_index;
  std::vector<double> rho_log_prior;
  
  double *y; // outcome for individual i, observation j is *(y + (start_index[i] + j)) or y[start_index[i] + j]
  double *x; // j^th observation of X_k, individual i is: *(x + (k + (start_index[i] + j)*p)
  double *z; // j^th observation of Z_r, individuali is: *(z + (r + (start_index[i] + j)*R)
  
  
  double* af; // allfit
  double* rf; // full residual
  
  size_t k; // keeps track of which function we are updating at any one time
  // constructor
  data_info(){N = 0; p = 0; R = 0; n = std::vector<size_t>(1); start_index = std::vector<size_t>(1); rho_log_prior = std::vector<double>(1); y=0; x=0;z=0;af=0;rf=0;k = 0;}
};
// possibly we will have to add some attributes to handle the sparse BART parameters


// We will instead track tau with a pointer -- this will be useful when we have half-Cauchy prior on the tau
class tree_prior_info{
public:
  double pbd; // probability of birth/death move
  double pb; // probability of birth given birth/death move
  std::vector<double> alpha;
  std::vector<double> beta;
  std::vector<double> A; // used to center the half-t prior on tau
  std::vector<double> nu; // degrees of freedome for the half-t prior on tau
  double* tau;
  double* rp; // will track partial residuals
  // constructor
  tree_prior_info(){pbd = 1.0; pb = 0.5; alpha = std::vector<double>(1); beta = std::vector<double>(1); A = std::vector<double>(1); nu = std::vector<double>(1); tau = 0; rp = 0;}

};

class sigma_prior_info{
public:
  double sigma_hat; // honestly this is never really used
  double lambda; // only used if we use an inverse Gamma
  double nu; // only used if we use an inverse Gamma
  double A; // used to center the half-cauchy prior
  // constructor
  sigma_prior_info(){sigma_hat = 1.0; lambda = 1.0; nu = 3; A = 0.0;}
};


//============================================================
//sufficient statistics for 1 node
//============================================================

// for VC-BART, we need to keep track of which observations for each individual
// belng to which node. so sv.I needs to be a vector of length N
class sinfo
{
public:
  size_t n;
  std::vector<size_t> node_count; // replaces sv.n. sv.node_count[i] counts the number of observations j for individual i that land in this node
  std::vector<std::vector<size_t > > I; // sv.I is vector of length N. sv.I[i] is a vector containing the sv.node_count[i] indices that land in this node
  sinfo(){n = 0;node_count = std::vector<size_t>(1); I = std::vector<std::vector<size_t> > (1, std::vector<size_t>(1));}
};






#endif

// Below are the old versions





