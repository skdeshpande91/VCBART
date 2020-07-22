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
double pgrow(tree::tree_p n, xinfo &xi, tree_prior_info &tree_pi, size_t k); // k tells us which beta function we're updating
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
