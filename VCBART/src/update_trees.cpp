#include "update_trees.h"
void grow_tree_ind(tree &t, int &accept, int &j, double &sigma, suff_stat &ss_train, suff_stat &ss_test, data_info &di_train, data_info &di_test, tree_prior_info &tree_pi, RNG &gen)
{
  std::vector<int> bn_nid_vec; // vector to hold node id's for all bottom nodes/leafs in the tree
  for(suff_stat_it l_it = ss_train.begin(); l_it != ss_train.end(); ++l_it) bn_nid_vec.push_back(l_it->first);
  int ni = floor(gen.uniform() * bn_nid_vec.size()); // randomly pick the index of the leaf from which we propose the GROW
  int nx_nid = bn_nid_vec[ni];// id of the leaf from which we propose the GROW
  int nxl_nid = 2*nx_nid;
  int nxr_nid = 2*nx_nid+1;
  tree::tree_p nx = t.get_ptr(nx_nid); // pointer to nx
  tree::tree_cp nxp = nx->get_p(); // constant pointer to the PARENT of nx.
  
  // we are ready to compute the log transition ratio:
  double q_grow_old = tree_pi.prob_b; // transition prob. of GROWing old tree into new tree
  double q_prune_new = 1.0 - tree_pi.prob_b; // transition prob. of pruning new tree into old tree
  
  int nleaf_old = t.get_nbots(); // number of leaves in old tree
  int nnog_old = t.get_nnogs(); // number of nodes in old tree with no grandchildren (nog node)
  int nnog_new = nnog_old; // number of nog nods in new tree
  
  if(nxp == 0){
    // nx is the root so it has no parent. we ALWAYS propose growing here
    q_grow_old = 1.0; //
  } else if(nxp->is_nog()){
    // parent of nx has no grandchildren in old tree
    // in new tree, nxp has grandchildren but nx does not
    // hence nnog_new = nnog_old
    nnog_new = nnog_old;
  } else{
    // parent of nx has grandchildren in old tree and will continue to have grandchildren in new tree
    // nx has no grandchildren in the new tree
    nnog_new = 1 + nnog_old;
  }
  
  // numerator of transition ratio: P(uniformly pick a nod node in new tree) * P(decide to prune new tree)
  // denominator of transition ratio: P(uniformly pick a leaf node in old tree) * P(decide to grow old tree)
  
  double log_trans_ratio = (log(q_prune_new) - log( (double) nnog_new)) - (log(q_grow_old) - log( (double) nleaf_old));
  
  // for prior ratio:
  // numerator: P(grow at nx in new tree) * (1 - P(grow at nxl in new tree)) * (1 - P(grow at nxr in new tree))
  // denominator 1 - P(grow at nx in old tree)
  // Note: P(grow at nx in old tree) = alpha(1 + depth(nx))^(-beta)
  // P(grow at nxl in new tree) = alpha(1 + depth(nx) + 1)^(-beta)
  
  double p_grow_nx = tree_pi.alpha/pow(1.0 + (double) nx->get_depth(), tree_pi.beta);
  double p_grow_nxl = tree_pi.alpha/pow(2.0 + (double) nx->get_depth(), tree_pi.beta);
  double p_grow_nxr = tree_pi.alpha/pow(2.0 + (double) nx->get_depth(), tree_pi.beta);
  double log_prior_ratio = log(p_grow_nx) + log(1.0 - p_grow_nxl) + log(1.0 - p_grow_nxr) - log(1.0 - p_grow_nx);
  
  // we are now ready to draw a decision rule
  rule_t rule;
  draw_rule(rule, t, nx_nid, di_train, tree_pi, gen); // draw the rule
  
  // ready to compute the sufficient statistics map for the new tree
  suff_stat prop_ss_train; //
  compute_suff_stat_grow(ss_train, prop_ss_train, nx_nid, rule, t, di_train);
  
  suff_stat prop_ss_test;
  if(di_test.N > 0) compute_suff_stat_grow(ss_test, prop_ss_test, nx_nid, rule, t, di_test);
  
  // compute P & Theta for the old tree
  std::map<int,int> orig_leaf_map;
  arma::mat orig_P = arma::zeros<arma::mat>(ss_train.size(), ss_train.size());
  arma::vec orig_Theta = arma::zeros<arma::mat>(ss_train.size());
  compute_p_theta_ind(j, orig_P, orig_Theta, orig_leaf_map, ss_train, sigma, di_train, tree_pi);
  
  // compute P & Theta for the new tree
  std::map<int, int> prop_leaf_map;
  arma::mat prop_P = arma::zeros<arma::mat>(prop_ss_train.size(), prop_ss_train.size());
  arma::vec prop_Theta = arma::zeros<arma::vec>(prop_ss_train.size());
  compute_p_theta_ind(j, prop_P, prop_Theta, prop_leaf_map, prop_ss_train, sigma, di_train, tree_pi);
  
  double orig_lil = compute_lil(orig_P, orig_Theta, tree_pi);
  double prop_lil = compute_lil(prop_P, prop_Theta, tree_pi);
  
  double log_like_ratio = prop_lil - orig_lil;
  double log_alpha = log_like_ratio + log_prior_ratio + log_trans_ratio; // MH ratio
  if(gen.log_uniform() <= log_alpha){
    // accept the transition!!
    
    // bookkeping to update number of rules and count # of times variable is split upon
    ++(*tree_pi.rule_count); // increment running count of total number of ruels in j-th ensemble
    if(rule.is_aa && !rule.is_cat){
      ++(tree_pi.var_count->at(rule.v_aa));
    } else if(!rule.is_aa && rule.is_cat){
      int v_raw = rule.v_cat + di_train.R_cat;
      ++(tree_pi.var_count->at(v_raw));
    } else if(!rule.is_aa && !rule.is_cat){
      ++(*tree_pi.rc_rule_count);
      for(rc_it it = rule.rc_weight.begin(); it != rule.rc_weight.end(); ++it){
        ++(*tree_pi.rc_var_count); // update the *total* number of variables involved in random combination rule
      }
    } else{
      // we really should never hit this
      Rcpp::stop("[grow_tree_ind]: after accepting a GROW move, we cannot resolve rule type!");
    } // closes if/else about what type of rule we accepted
    
    // we now are ready to update ss, the sufficiant statistic map
    suff_stat_it nxl_it = prop_ss_train.find(nxl_nid); // iterator at element for nxl in the proposed suff_stat_map
    suff_stat_it nxr_it = prop_ss_train.find(nxr_nid); // iterator at element for nxr in the proposed suff_stat_map
    
    if(nxl_it == prop_ss_train.end() || nxr_it == prop_ss_train.end()){
      // could not find the key in prop_ss_train equal to nxl_nid or nxr_nid
      Rcpp::Rcout << "[grow_tree_ind]: sufficient stat map for training data not updated correctly in grow move!" << std::endl;
      Rcpp::Rcout << "  left child id = " << nxl_nid << " right child id " << nxr_nid << std::endl;
      Rcpp::Rcout << "  available ids in map:";
      for(suff_stat_it l_it = prop_ss_train.begin(); l_it != prop_ss_train.end(); ++l_it) Rcpp::Rcout << " " << l_it->first;
      Rcpp::Rcout << std::endl;
      Rcpp::stop("missing id for either left or right child in the proposed suff_stat_map!");
    } else{
      ss_train.insert(std::pair<int, std::vector<int>>(nxl_nid, nxl_it->second));
      ss_train.insert(std::pair<int, std::vector<int>>(nxr_nid, nxr_it->second));
      ss_train.erase(nx_nid); // remove element for nx in sufficient stat map.
      // note: erase's argument is the key, which in our case is just the node id for the leaf
    }
    
    if(di_test.N > 0){
      nxl_it = prop_ss_test.find(nxl_nid); // iterator at element for nxl in the proposed suff_stat_map
      nxr_it = prop_ss_test.find(nxr_nid); // iterator at element for nxr in the proposed suff_stat_map
      
      if(nxl_it == prop_ss_test.end() || nxr_it == prop_ss_test.end()){
        // could not find the key in prop_ss_test equal to nxl_nid or nxr_nid
        Rcpp::Rcout << "[grow_tree_ind]: sufficient stat map for testing data not updated correctly in grow move!" << std::endl;
        Rcpp::Rcout << "  left child id = " << nxl_nid << " right child id " << nxr_nid << std::endl;
        Rcpp::Rcout << "  available ids in map:";
        for(suff_stat_it l_it = prop_ss_test.begin(); l_it != prop_ss_test.end(); ++l_it) Rcpp::Rcout << " " << l_it->first;
        Rcpp::Rcout << std::endl;
        Rcpp::stop("missing id for either left or right child in the proposed suff_stat_map!");
      } else{
        ss_test.insert(std::pair<int, std::vector<int>>(nxl_nid, nxl_it->second));
        ss_test.insert(std::pair<int, std::vector<int>>(nxr_nid, nxr_it->second));
        ss_test.erase(nx_nid); // remove element for nx in sufficient stat map.
        // note: erase's argument is the key, which in our case is just the node id for the leaf
      }
    }
    t.birth(nx_nid, rule); // actually perform the birth
    draw_mu(t, prop_P, prop_Theta, prop_leaf_map, gen);
    accept = 1;
    
  } else{
    // reject the GROW move
    // just update the mu's
    draw_mu(t, orig_P, orig_Theta, orig_leaf_map, gen);
    accept = 0;
  }
}
void prune_tree_ind(tree &t, int &accept, int &j, double &sigma, suff_stat &ss_train, suff_stat &ss_test, data_info &di_train, data_info &di_test, tree_prior_info &tree_pi, RNG &gen)
{
  // first randomly select a nog node
  tree::npv nogs_vec;
  t.get_nogs(nogs_vec);
  
  int ni = floor(gen.uniform() * nogs_vec.size());
  tree::tree_p nx = nogs_vec[ni]; // pointer to nx, the randomly select nog node that will eventually be a leaf in new tree
  tree::tree_p nxl = nx->get_l(); // pointer to left child of nx that will be pruned
  tree::tree_p nxr = nx->get_r(); // pointer to right child of nx that will be pruned
  
  // transition ratio
  double q_prune_old = 1.0 - tree_pi.prob_b; // transition prob that we prune old tree
  double q_grow_new = tree_pi.prob_b; // transition prob that we grow new tree
  tree::tree_p nxp = nx->get_p(); // pointer to the parent node of nx in old tree
  if(nxp == 0) q_grow_new = 1.0; // nx is the top node so new tree is just the root & we'd always propose GROW
  else q_grow_new = tree_pi.prob_b;
  
  int nleaf_new = t.get_nbots() - 1; // new tree has one less leaf node than old tree
  int nnog_old = t.get_nnogs(); // number of nog nodes in old tree
  
  // numerator of transition ratio: P(uniformly pick a leaf node in new tree) * P(decide to grow new tree)
  // denominator of transition ratio: P(uniformly pick a nog node in old tree) * P(decide to prune old tree)
  
  double log_trans_ratio = (log(q_grow_new) - log(nleaf_new)) - (log(q_prune_old) - log(nnog_old));
  
  // prior ratio:
  // numerator: 1 - P(grow at nx in new tree)
  // denominator: P(grow at nx in old tree) * (1 - P(grow at nxl in old tree)) * (1 - P(grow at nxr in old tree)
  double p_grow_nx = tree_pi.alpha/pow(1.0 + (double) nx->get_depth(), tree_pi.beta);
  double p_grow_nxl = tree_pi.alpha/pow(2.0 + (double) nx->get_depth(), tree_pi.beta);
  double p_grow_nxr = tree_pi.alpha/pow(2.0 + (double) nx->get_depth(), tree_pi.beta);
  double log_prior_ratio = log(1.0 - p_grow_nx) - (log(1.0 - p_grow_nxl) + log(1.0 - p_grow_nxr) + log(p_grow_nx));
  
  // likelihood ratio
  suff_stat prop_ss_train;
  int nx_nid = nx->get_nid();
  int nxl_nid = nxl->get_nid();
  int nxr_nid = nxr->get_nid();
  compute_suff_stat_prune(ss_train, prop_ss_train, nxl_nid, nxr_nid, nx_nid, t, di_train);
  
  suff_stat prop_ss_test;
  if(di_test.N > 0) compute_suff_stat_prune(ss_test, prop_ss_test, nxl_nid, nxr_nid, nx_nid, t, di_test);
  
  // compute P & Theta for the old tree
  std::map<int,int> orig_leaf_map;
  arma::mat orig_P = arma::zeros<arma::mat>(ss_train.size(), ss_train.size());
  arma::vec orig_Theta = arma::zeros<arma::mat>(ss_train.size());
  compute_p_theta_ind(j, orig_P, orig_Theta, orig_leaf_map, ss_train, sigma, di_train, tree_pi);
  
  // compute P & Theta for the new tree
  std::map<int, int> prop_leaf_map;
  arma::mat prop_P = arma::zeros<arma::mat>(prop_ss_train.size(), prop_ss_train.size());
  arma::vec prop_Theta = arma::zeros<arma::vec>(prop_ss_train.size());
  compute_p_theta_ind(j, prop_P, prop_Theta, prop_leaf_map, prop_ss_train, sigma, di_train, tree_pi);
  
  double orig_lil = compute_lil(orig_P, orig_Theta, tree_pi);
  double prop_lil = compute_lil(prop_P, prop_Theta, tree_pi);
  
  double log_like_ratio = prop_lil - orig_lil;
  double log_alpha = log_like_ratio + log_prior_ratio + log_trans_ratio; // MH ratio
  
  if(gen.log_uniform() <= log_alpha){
    // accept PRUNE proposal
    
    --(*tree_pi.rule_count);
    if(nx->get_is_aa() && !nx->get_is_cat()){
      --(tree_pi.var_count->at(nx->get_v_aa()));
    } else if(!nx->get_is_aa() && nx->get_is_cat()){
      int v_raw = di_train.R_cont + nx->get_v_cat();
      --(tree_pi.var_count->at(v_raw));
    } else if(!nx->get_is_aa() && !nx->get_is_cat()){
      std::map<int,double> rc_weight = nx->get_rc_weight();
      --(*tree_pi.rc_rule_count);
      for(rc_it it = rc_weight.begin(); it != rc_weight.end(); ++it){
        --(*tree_pi.rc_var_count);
      }
    } else{
      Rcpp::Rcout << "[prune tree]: accepted a prune at nog node " << nx_nid << " but unable to figure out its rule type" << std::endl;
      t.print();
      Rcpp::stop("Cannot resolve rule type!");
    }
    
    // adjust suff_stat map
    suff_stat_it nx_it = prop_ss_train.find(nx_nid); // iterator at element for nx in suff_stat map for new tree
    if(nx_it == prop_ss_train.end()){
      // did not find nx_nid among the keys of prop_ss_train
      Rcpp::Rcout << "[prune_tree]: did not find id of new leaf nod in keys of training sufficient stat map" << std::endl;
      Rcpp::Rcout << "  id of new leaf: " << nx_nid << std::endl;
      Rcpp::Rcout << "  ids in map:";
      for(suff_stat_it l_it = prop_ss_train.begin(); l_it != prop_ss_train.end(); ++l_it) Rcpp::Rcout << " " << l_it->first;
      Rcpp::Rcout << std::endl;
    } else{
      ss_train.erase(nxl_nid); // delete entry for nxl in suff_stat map
      ss_train.erase(nxr_nid); // delete entry for nxr in suff_stat map
      ss_train.insert(std::pair<int, std::vector<int>>(nx_nid, nx_it->second)); // add an entry for nx in suff stat map
    }
    
    if(di_test.N > 0){
      nx_it = prop_ss_test.find(nx_nid); // iterator at element for nx in suff_stat map for new tree
      if(nx_it == prop_ss_test.end()){
        // did not find nx_nid among the keys of prop_ss_train
        Rcpp::Rcout << "[prune_tree]: did not find id of new leaf nod in keys of testing sufficient stat map" << std::endl;
        Rcpp::Rcout << "  id of new leaf: " << nx_nid << std::endl;
        Rcpp::Rcout << "  ids in map:";
        for(suff_stat_it l_it = prop_ss_test.begin(); l_it != prop_ss_test.end(); ++l_it) Rcpp::Rcout << " " << l_it->first;
        Rcpp::Rcout << std::endl;
      } else{
        ss_test.erase(nxl_nid); // delete entry for nxl in suff_stat map
        ss_test.erase(nxr_nid); // delete entry for nxr in suff_stat map
        ss_test.insert(std::pair<int, std::vector<int>>(nx_nid, nx_it->second)); // add an entry for nx in suff stat map
      }
    }
    t.death(nx_nid);
    accept = 1;
    draw_mu(t, prop_P, prop_Theta, prop_leaf_map, gen);
  } else{
    // reject proposal
    draw_mu(t, orig_P, orig_Theta, orig_leaf_map, gen);
    accept = 0;
  }


  
}
void update_tree_ind(tree &t, int &accept, int &j, double &sigma, suff_stat &ss_train, suff_stat &ss_test, data_info &di_train, data_info &di_test, tree_prior_info &tree_pi, RNG &gen)
{
  accept = 0;
  double PBx = tree_pi.prob_b; // prob of proposing a GROW move (typically 0.5)
  if(t.get_treesize() == 1) PBx = 1.0; // if tree is the root, we always propose a GROW move
  if(gen.uniform() < PBx) grow_tree_ind(t, accept, j, sigma, ss_train, ss_test, di_train, di_test, tree_pi, gen);
  else prune_tree_ind(t, accept, j, sigma, ss_train, ss_test, di_train, di_test, tree_pi, gen);
  // NOTE: grow_tree and prune_tree call draw_mu internally so no need to call it again here.
}
void grow_tree_cs(tree &t, int &accept, int &j, double &rho, double &sigma, suff_stat &ss_train, suff_stat &ss_test, data_info &di_train, data_info &di_test, tree_prior_info &tree_pi, RNG &gen)
{
  std::vector<int> bn_nid_vec; // vector to hold node id's for all bottom nodes/leafs in the tree
  for(suff_stat_it l_it = ss_train.begin(); l_it != ss_train.end(); ++l_it) bn_nid_vec.push_back(l_it->first);
  int ni = floor(gen.uniform() * bn_nid_vec.size()); // randomly pick the index of the leaf from which we propose the GROW
  int nx_nid = bn_nid_vec[ni];// id of the leaf from which we propose the GROW
  int nxl_nid = 2*nx_nid;
  int nxr_nid = 2*nx_nid+1;
  
  tree::tree_p nx = t.get_ptr(nx_nid); // pointer to nx
  tree::tree_cp nxp = nx->get_p(); // constant pointer to the PARENT of nx.
  
  // we are ready to compute the log transition ratio:
  double q_grow_old = tree_pi.prob_b; // transition prob. of GROWing old tree into new tree
  double q_prune_new = 1.0 - tree_pi.prob_b; // transition prob. of pruning new tree into old tree
  
  int nleaf_old = t.get_nbots(); // number of leaves in old tree
  int nnog_old = t.get_nnogs(); // number of nodes in old tree with no grandchildren (nog node)
  int nnog_new = nnog_old; // number of nog nods in new tree
  
  if(nxp == 0){
    // nx is the root so it has no parent. we ALWAYS propose growing here
    q_grow_old = 1.0; //
  } else if(nxp->is_nog()){
    // parent of nx has no grandchildren in old tree
    // in new tree, nxp has grandchildren but nx does not
    // hence nnog_new = nnog_old
    nnog_new = nnog_old;
  } else{
    // parent of nx has grandchildren in old tree and will continue to have grandchildren in new tree
    // nx has no grandchildren in the new tree
    nnog_new = 1 + nnog_old;
  }
  
  // numerator of transition ratio: P(uniformly pick a nod node in new tree) * P(decide to prune new tree)
  // denominator of transition ratio: P(uniformly pick a leaf node in old tree) * P(decide to grow old tree)
  
  double log_trans_ratio = (log(q_prune_new) - log( (double) nnog_new)) - (log(q_grow_old) - log( (double) nleaf_old));
  
  // for prior ratio:
  // numerator: P(grow at nx in new tree) * (1 - P(grow at nxl in new tree)) * (1 - P(grow at nxr in new tree))
  // denominator 1 - P(grow at nx in old tree)
  // Note: P(grow at nx in old tree) = alpha(1 + depth(nx))^(-beta)
  // P(grow at nxl in new tree) = alpha(1 + depth(nx) + 1)^(-beta)
  
  double p_grow_nx = tree_pi.alpha/pow(1.0 + (double) nx->get_depth(), tree_pi.beta);
  double p_grow_nxl = tree_pi.alpha/pow(2.0 + (double) nx->get_depth(), tree_pi.beta);
  double p_grow_nxr = tree_pi.alpha/pow(2.0 + (double) nx->get_depth(), tree_pi.beta);
  double log_prior_ratio = log(p_grow_nx) + log(1.0 - p_grow_nxl) + log(1.0 - p_grow_nxr) - log(1.0 - p_grow_nx);
  
  // we are now ready to draw a decision rule
  rule_t rule;
  draw_rule(rule, t, nx_nid, di_train, tree_pi, gen); // draw the rule
  
  // ready to compute the sufficient statistics map for the new tree
  suff_stat prop_ss_train; //
  compute_suff_stat_grow(ss_train, prop_ss_train, nx_nid, rule, t, di_train);
  
  suff_stat prop_ss_test;
  if(di_test.N > 0) compute_suff_stat_grow(ss_test, prop_ss_test, nx_nid, rule, t, di_test);
  
  // compute P & Theta for the old tree
  std::map<int,int> orig_leaf_map;
  arma::mat orig_P = arma::zeros<arma::mat>(ss_train.size(), ss_train.size());
  arma::vec orig_Theta = arma::zeros<arma::mat>(ss_train.size());
  compute_p_theta_cs(j, orig_P, orig_Theta, orig_leaf_map, ss_train, rho, sigma, di_train, tree_pi);
  
  // compute P & Theta for the new tree
  std::map<int, int> prop_leaf_map;
  arma::mat prop_P = arma::zeros<arma::mat>(prop_ss_train.size(), prop_ss_train.size());
  arma::vec prop_Theta = arma::zeros<arma::vec>(prop_ss_train.size());
  compute_p_theta_cs(j, prop_P, prop_Theta, prop_leaf_map, prop_ss_train, rho, sigma, di_train, tree_pi);
  
  double orig_lil = compute_lil(orig_P, orig_Theta, tree_pi);
  double prop_lil = compute_lil(prop_P, prop_Theta, tree_pi);
  
  double log_like_ratio = prop_lil - orig_lil;
  double log_alpha = log_like_ratio + log_prior_ratio + log_trans_ratio; // MH ratio
  if(gen.uniform() <= log_alpha){
    // accept the transition!!
    
    // bookkeping to update number of rules and count # of times variable is split upon
    ++(*tree_pi.rule_count); // increment running count of total number of ruels in j-th ensemble
    if(rule.is_aa && !rule.is_cat){
      ++(tree_pi.var_count->at(rule.v_aa));
    } else if(!rule.is_aa && rule.is_cat){
      int v_raw = rule.v_cat + di_train.R_cat;
      ++(tree_pi.var_count->at(v_raw));
    } else if(!rule.is_aa && !rule.is_cat){
      ++(*tree_pi.rc_rule_count);
      for(rc_it it = rule.rc_weight.begin(); it != rule.rc_weight.end(); ++it){
        ++(*tree_pi.rc_var_count); // update the *total* number of variables involved in random combination rule
      }
    } else{
      // we really should never hit this
      Rcpp::stop("[grow_tree_ind]: after accepting a GROW move, we cannot resolve rule type!");
    } // closes if/else about what type of rule we accepted
    
    // we now are ready to update ss, the sufficiant statistic map
    suff_stat_it nxl_it = prop_ss_train.find(nxl_nid); // iterator at element for nxl in the proposed suff_stat_map
    suff_stat_it nxr_it = prop_ss_train.find(nxr_nid); // iterator at element for nxr in the proposed suff_stat_map
    
    if(nxl_it == prop_ss_train.end() || nxr_it == prop_ss_train.end()){
      // could not find the key in prop_ss_train equal to nxl_nid or nxr_nid
      Rcpp::Rcout << "[grow_tree_ind]: sufficient stat map for training data not updated correctly in grow move!" << std::endl;
      Rcpp::Rcout << "  left child id = " << nxl_nid << " right child id " << nxr_nid << std::endl;
      Rcpp::Rcout << "  available ids in map:";
      for(suff_stat_it l_it = prop_ss_train.begin(); l_it != prop_ss_train.end(); ++l_it) Rcpp::Rcout << " " << l_it->first;
      Rcpp::Rcout << std::endl;
      Rcpp::stop("missing id for either left or right child in the proposed suff_stat_map!");
    } else{
      ss_train.insert(std::pair<int, std::vector<int>>(nxl_nid, nxl_it->second));
      ss_train.insert(std::pair<int, std::vector<int>>(nxr_nid, nxr_it->second));
      ss_train.erase(nx_nid); // remove element for nx in sufficient stat map.
      // note: erase's argument is the key, which in our case is just the node id for the leaf
    }
    
    if(di_test.N > 0){
      nxl_it = prop_ss_test.find(nxl_nid); // iterator at element for nxl in the proposed suff_stat_map
      nxr_it = prop_ss_test.find(nxr_nid); // iterator at element for nxr in the proposed suff_stat_map
      
      if(nxl_it == prop_ss_test.end() || nxr_it == prop_ss_test.end()){
        // could not find the key in prop_ss_test equal to nxl_nid or nxr_nid
        Rcpp::Rcout << "[grow_tree_ind]: sufficient stat map for testing data not updated correctly in grow move!" << std::endl;
        Rcpp::Rcout << "  left child id = " << nxl_nid << " right child id " << nxr_nid << std::endl;
        Rcpp::Rcout << "  available ids in map:";
        for(suff_stat_it l_it = prop_ss_test.begin(); l_it != prop_ss_test.end(); ++l_it) Rcpp::Rcout << " " << l_it->first;
        Rcpp::Rcout << std::endl;
        Rcpp::stop("missing id for either left or right child in the proposed suff_stat_map!");
      } else{
        ss_test.insert(std::pair<int, std::vector<int>>(nxl_nid, nxl_it->second));
        ss_test.insert(std::pair<int, std::vector<int>>(nxr_nid, nxr_it->second));
        ss_test.erase(nx_nid); // remove element for nx in sufficient stat map.
        // note: erase's argument is the key, which in our case is just the node id for the leaf
      }
    }
    
    t.birth(nx_nid, rule); // actually perform the birth
    draw_mu(t, prop_P, prop_Theta, prop_leaf_map, gen);
    accept = 1;
    
  } else{
    // reject the GROW move
    // just update the mu's
    draw_mu(t, orig_P, orig_Theta, orig_leaf_map, gen);
    accept = 0;
  }
}
void prune_tree_cs(tree &t, int &accept, int &j, double &rho, double &sigma, suff_stat &ss_train, suff_stat &ss_test, data_info &di_train, data_info &di_test, tree_prior_info &tree_pi, RNG &gen)
{
  // first randomly select a nog node
  tree::npv nogs_vec;
  t.get_nogs(nogs_vec);
  
  int ni = floor(gen.uniform() * nogs_vec.size());
  tree::tree_p nx = nogs_vec[ni]; // pointer to nx, the randomly select nog node that will eventually be a leaf in new tree
  tree::tree_p nxl = nx->get_l(); // pointer to left child of nx that will be pruned
  tree::tree_p nxr = nx->get_r(); // pointer to right child of nx that will be pruned
  
  // transition ratio
  double q_prune_old = 1.0 - tree_pi.prob_b; // transition prob that we prune old tree
  double q_grow_new = tree_pi.prob_b; // transition prob that we grow new tree
  tree::tree_p nxp = nx->get_p(); // pointer to the parent node of nx in old tree
  if(nxp == 0) q_grow_new = 1.0; // nx is the top node so new tree is just the root & we'd always propose GROW
  else q_grow_new = tree_pi.prob_b;
  
  int nleaf_new = t.get_nbots() - 1; // new tree has one less leaf node than old tree
  int nnog_old = t.get_nnogs(); // number of nog nodes in old tree
  
  // numerator of transition ratio: P(uniformly pick a leaf node in new tree) * P(decide to grow new tree)
  // denominator of transition ratio: P(uniformly pick a nog node in old tree) * P(decide to prune old tree)
  
  double log_trans_ratio = (log(q_grow_new) - log(nleaf_new)) - (log(q_prune_old) - log(nnog_old));
  
  // prior ratio:
  // numerator: 1 - P(grow at nx in new tree)
  // denominator: P(grow at nx in old tree) * (1 - P(grow at nxl in old tree)) * (1 - P(grow at nxr in old tree)
  double p_grow_nx = tree_pi.alpha/pow(1.0 + (double) nx->get_depth(), tree_pi.beta);
  double p_grow_nxl = tree_pi.alpha/pow(2.0 + (double) nx->get_depth(), tree_pi.beta);
  double p_grow_nxr = tree_pi.alpha/pow(2.0 + (double) nx->get_depth(), tree_pi.beta);
  double log_prior_ratio = log(1.0 - p_grow_nx) - (log(1.0 - p_grow_nxl) + log(1.0 - p_grow_nxr) + log(p_grow_nx));
  
  // likelihood ratio
  suff_stat prop_ss_train;
  int nx_nid = nx->get_nid();
  int nxl_nid = nxl->get_nid();
  int nxr_nid = nxr->get_nid();
  compute_suff_stat_prune(ss_train, prop_ss_train, nxl_nid, nxr_nid, nx_nid, t, di_train);
  
  suff_stat prop_ss_test;
  if(di_test.N > 0) compute_suff_stat_prune(ss_test, prop_ss_test, nxl_nid, nxr_nid, nx_nid, t, di_test);
  
  // compute P & Theta for the old tree
  std::map<int,int> orig_leaf_map;
  arma::mat orig_P = arma::zeros<arma::mat>(ss_train.size(), ss_train.size());
  arma::vec orig_Theta = arma::zeros<arma::mat>(ss_train.size());
  compute_p_theta_cs(j, orig_P, orig_Theta, orig_leaf_map, ss_train, rho, sigma, di_train, tree_pi);
  
  // compute P & Theta for the new tree
  std::map<int, int> prop_leaf_map;
  arma::mat prop_P = arma::zeros<arma::mat>(prop_ss_train.size(), prop_ss_train.size());
  arma::vec prop_Theta = arma::zeros<arma::vec>(prop_ss_train.size());
  compute_p_theta_cs(j, prop_P, prop_Theta, prop_leaf_map, prop_ss_train, rho, sigma, di_train, tree_pi);
  
  double orig_lil = compute_lil(orig_P, orig_Theta, tree_pi);
  double prop_lil = compute_lil(prop_P, prop_Theta, tree_pi);
  
  double log_like_ratio = prop_lil - orig_lil;
  double log_alpha = log_like_ratio + log_prior_ratio + log_trans_ratio; // MH ratio
  
  if(gen.log_uniform() <= log_alpha){
    // accept PRUNE proposal
    
    --(*tree_pi.rule_count);
    if(nx->get_is_aa() && !nx->get_is_cat()){
      --(tree_pi.var_count->at(nx->get_v_aa()));
    } else if(!nx->get_is_aa() && nx->get_is_cat()){
      int v_raw = di_train.R_cont + nx->get_v_cat();
      --(tree_pi.var_count->at(v_raw));
    } else if(!nx->get_is_aa() && !nx->get_is_cat()){
      std::map<int,double> rc_weight = nx->get_rc_weight();
      --(*tree_pi.rc_rule_count);
      for(rc_it it = rc_weight.begin(); it != rc_weight.end(); ++it){
        --(*tree_pi.rc_var_count);
      }
    } else{
      Rcpp::Rcout << "[prune tree]: accepted a prune at nog node " << nx_nid << " but unable to figure out its rule type" << std::endl;
      t.print();
      Rcpp::stop("Cannot resolve rule type!");
    }
    
    // adjust suff_stat map
    suff_stat_it nx_it = prop_ss_train.find(nx_nid); // iterator at element for nx in suff_stat map for new tree
    if(nx_it == prop_ss_train.end()){
      // did not find nx_nid among the keys of prop_ss_train
      Rcpp::Rcout << "[prune_tree]: did not find id of new leaf nod in keys of training sufficient stat map" << std::endl;
      Rcpp::Rcout << "  id of new leaf: " << nx_nid << std::endl;
      Rcpp::Rcout << "  ids in map:";
      for(suff_stat_it l_it = prop_ss_train.begin(); l_it != prop_ss_train.end(); ++l_it) Rcpp::Rcout << " " << l_it->first;
      Rcpp::Rcout << std::endl;
    } else{
      ss_train.erase(nxl_nid); // delete entry for nxl in suff_stat map
      ss_train.erase(nxr_nid); // delete entry for nxr in suff_stat map
      ss_train.insert(std::pair<int, std::vector<int>>(nx_nid, nx_it->second)); // add an entry for nx in suff stat map
    }
    
    if(di_test.N > 0){
      nx_it = prop_ss_test.find(nx_nid); // iterator at element for nx in suff_stat map for new tree
      if(nx_it == prop_ss_test.end()){
        // did not find nx_nid among the keys of prop_ss_train
        Rcpp::Rcout << "[prune_tree]: did not find id of new leaf nod in keys of testing sufficient stat map" << std::endl;
        Rcpp::Rcout << "  id of new leaf: " << nx_nid << std::endl;
        Rcpp::Rcout << "  ids in map:";
        for(suff_stat_it l_it = prop_ss_test.begin(); l_it != prop_ss_test.end(); ++l_it) Rcpp::Rcout << " " << l_it->first;
        Rcpp::Rcout << std::endl;
      } else{
        ss_test.erase(nxl_nid); // delete entry for nxl in suff_stat map
        ss_test.erase(nxr_nid); // delete entry for nxr in suff_stat map
        ss_test.insert(std::pair<int, std::vector<int>>(nx_nid, nx_it->second)); // add an entry for nx in suff stat map
      }
    }
    t.death(nx_nid);
    accept = 1;
    draw_mu(t, prop_P, prop_Theta, prop_leaf_map, gen);
  } else{
    // reject proposal
    draw_mu(t, orig_P, orig_Theta, orig_leaf_map, gen);
    accept = 0;
  }
}
void update_tree_cs(tree &t, int &accept, int &j, double &rho, double &sigma, suff_stat &ss_train, suff_stat &ss_test, data_info &di_train, data_info &di_test, tree_prior_info &tree_pi, RNG &gen)
{
  accept = 0; // initialize MH acceptance indicator to 0
  double PBx = tree_pi.prob_b; // prob of proposing a GROW move (typically 0.5)
  if(t.get_treesize() == 1) PBx = 1.0; // if tree is just the root, we always propose a GROW move
  
  if(gen.uniform() < PBx) grow_tree_cs(t, accept, j, rho, sigma, ss_train, ss_test, di_train, di_test, tree_pi, gen);
  else prune_tree_cs(t, accept, j, rho, sigma, ss_train, ss_test, di_train, di_test, tree_pi, gen);
  // NOTE: grow_tree and prune_tree call draw_mu internally so there's no need to call it here
}
