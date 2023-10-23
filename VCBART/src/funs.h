

#ifndef GUARD_funs_h
#define GUARD_funs_h

#include "tree.h"


void tree_traversal(suff_stat &ss, tree &t, data_info &di);


//void fit_single_tree(double* ftemp, tree &t, data_info &di);
void fit_ensemble(std::vector<double> &fit, std::vector<tree> &t_vec, data_info &di);

void compute_suff_stat_grow(suff_stat &orig_suff_stat, suff_stat &new_suff_stat, int &nx_nid, rule_t &rule, tree &t, data_info &di);
void compute_suff_stat_prune(suff_stat &orig_suff_stat, suff_stat &new_suff_stat, int &nl_nid, int &nr_nid, int &np_nid, tree &t, data_info &di);

void compute_p_theta_ind(int &j, arma::mat &P, arma::vec &Theta, std::map<int,int> &leaf_map, suff_stat &ss, double &sigma, data_info &di, tree_prior_info &tree_pi);

void compute_p_theta_cs(int &j, arma::mat &P, arma::vec &Theta, std::map<int,int> &leaf_map, suff_stat &ss, double &rho, double &sigma, data_info &di, tree_prior_info &tree_pi);


double compute_lil(arma::mat &P, arma::vec &Theta, tree_prior_info &tree_pi);

void draw_mu(tree &t, arma::mat &P, arma::vec &Theta, std::map<int,int> &leaf_map, RNG &gen);

void draw_rule(rule_t &rule, tree &t, int &nid, data_info &di, tree_prior_info &tree_pi, RNG &gen);
std::string write_tree(tree &t, tree_prior_info &tree_pi, set_str_conversion &set_str);
void read_tree(tree &t, std::string &tree_string, set_str_conversion &set_str);

void build_symmetric_edge_map(edge_map &emap, std::vector<edge> &edges, std::set<int> &vertices);
std::vector<edge> get_induced_edges(std::vector<edge> &edges, std::set<int> &vertex_subset);
void dfs(int v, std::map<int, bool> &visited, std::vector<int> &comp, edge_map &emap);
void find_components(std::vector<std::vector<int> > &components, std::vector<edge> &edges, std::set<int> &vertices);
void get_unique_edges(std::vector<edge> &edges);
void boruvka(std::vector<edge> &mst_edges, std::vector<edge> &edges, std::set<int> &vertices);
void wilson(std::vector<edge> &mst_edges, std::vector<edge> &edges, std::set<int> &vertices, RNG &gen);
void graph_partition(std::set<int> &avail_levels, std::set<int> &l_vals, std::set<int> &r_vals, std::vector<edge> &orig_edges, int &K, int &cut_type, RNG &gen);
void update_theta_u(std::vector<double> &theta, double &u, std::vector<int> &var_count, int &R, double &a_u, double &b_u, RNG &gen);
void update_theta_rc(double& theta_rc, int &rc_var_count, int &rc_rule_count, double &a_rc, double &b_rc, int &R_cont, RNG &gen);
void update_rho(double &rho, double &sigma, data_info &di, RNG &gen);

void update_sigma_ind(double &sigma, double &nu, double &lambda, data_info &di, RNG &gen);
void update_sigma_cs(double &sigma, double &rho, double &nu, double &lambda, data_info &di, RNG &gen);
#endif /* funs_h */
