

#ifndef GUARD_funs_h
#define GUARD_funs_h

#include "tree.h"


//void fit(tree &t, data_info &di, double* ftemp); // vectorized fit of individual tree
//void fit_vec(std::vector<tree> &t_vec, data_info &di, double* ftemp); // vectorized fit of entire ensemble

//void fit(double* ftemp, tree &t, data_info &di);
//void fit_vec(double* ftemp, std::vector<tree> &t_vec, data_info &di);

//void fit(arma::vec &ftemp, tree &t, data_info &di);
//void fit_vec(arma::vec &ftemp, std::vector<tree> &t_vec, data_info &di);

//void compute_suff_stat_grow(suff_stat &sl, suff_stat &sr, suff_stat &sp, tree &t, tree::tree_cp nx, rule_t &rule, data_info &di);
//void compute_suff_stat_prune(suff_stat &sl, suff_stat &sr, suff_stat &sp, tree &t, tree::tree_cp nl, tree::tree_cp nr, data_info &di);
//void compute_suff_stat_vec(std::vector<suff_stat> &sv, tree::npv &bnv, tree &t, data_info &di);

//double compute_lil(suff_stat &ss, tree_prior_info &tree_pi); // compute single leaf contribution to log integrated likelihood
//double compute_lil_ratio_grow(tree &t, tree::tree_cp nx, rule_t &rule, data_info &di, tree_prior_info &tree_pi, bool debug);
//double compute_lil_ratio_prune(tree &t, tree::tree_cp nxl, tree::tree_cp nxr, data_info &di, tree_prior_info &tree_pi, bool debug);

//void compute_mu_post(double& post_mean, double& post_sd, suff_stat &ss, tree_prior_info &tree_pi);
//void draw_mu(tree &t, data_info &di, tree_prior_info &tree_pi, RNG &gen);

//void compute_cf_probs(size_t v_target, double* cf_probs, std::vector<tree> &t_vec, data_info &di);
//void compute_cf_probs(size_t v_target, arma::vec &cf_probs, std::vector<tree> &t_vec, data_info &di);

#endif /* funs_h */
