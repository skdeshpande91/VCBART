#ifndef GUARD_update_trees_h
#define GUARD_update_trees_h

#include "funs.h"

void grow_tree_ind(tree &t, int &accept, int &j, double &sigma, suff_stat &ss_train, suff_stat &ss_test, data_info &di_train, data_info &di_test, tree_prior_info &tree_pi, RNG &gen);
void prune_tree_ind(tree &t, int &accept, int &j, double &sigma, suff_stat &ss_train, suff_stat &ss_test, data_info &di_train, data_info &di_test, tree_prior_info &tree_pi, RNG &gen);
void update_tree_ind(tree &t, int &accept, int &j, double &sigma, suff_stat &ss_train, suff_stat &ss_test, data_info &di_train, data_info &di_test, tree_prior_info &tree_pi, RNG &gen);

void grow_tree_cs(tree &t, int &accept, int &j, double &rho, double &sigma, suff_stat &ss_train, suff_stat &ss_test, data_info &di_train, data_info &di_test, tree_prior_info &tree_pi, RNG &gen);
void prune_tree_cs(tree &t, int &accept, int &j, double &rho, double &sigma, suff_stat &ss_train, suff_stat &ss_test, data_info &di_train, data_info &di_test, tree_prior_info &tree_pi, RNG &gen);
void update_tree_cs(tree &t, int &accept, int &j, double &rho, double &sigma, suff_stat &ss_train, suff_stat &ss_test, data_info &di_train, data_info &di_test, tree_prior_info &tree_pi, RNG &gen);

#endif
