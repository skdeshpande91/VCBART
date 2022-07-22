#include "funs.h"

void tree_traversal(suff_stat &ss, tree &t, data_info &di)
{
  double* zz_cont = 0;
  int* zz_cat = 0;
  tree::tree_cp bn;
  int nid;
  ss.clear(); // clear out the sufficient statistic map
  
  tree::npv bnv;
  t.get_bots(bnv);
  suff_stat_it leaf_it;
  
  // add an element to suff stat map for each bottom node
  for(tree::npv_it it = bnv.begin(); it != bnv.end(); ++it){
    nid = (*it)->get_nid(); // bnv is a vector of pointers. it points to elements of bnv so we need (*it) to access the members of these elements
    // add an element to ss with key = nid and value = vector
    ss.insert(std::pair<int, std::vector<int>>(nid, std::vector<int>()));
  }
  
  // ready to loop over all observations
  for(int i = 0; i < di.N; i++){
    if(di.z_cont != 0) zz_cont = di.z_cont + i * di.R_cont;
    if(di.z_cat != 0) zz_cat = di.z_cat + i * di.R_cat;
    bn = t.get_bn(zz_cont, zz_cat);
    if(bn == 0){
      Rcpp::Rcout << "i = " << i << std::endl;
      t.print();
      Rcpp::stop("[tree_traversal]: could not find bottom node!");
    }
    else{
      nid = bn->get_nid();
      if(ss.count(nid) != 1) Rcpp::stop("[tree_traversal]: bottom node not included in sufficient statistic map!"); // we should *NEVER* hit this error
      else{
        // ss.find(nid)->second.push_back(i); // this might be getting too cute...
        leaf_it = ss.find(nid);
        leaf_it->second.push_back(i);
      } // closes if/else checking that ss has an element corresponding to the leaf node to which i was mapped.
    } // closes if/else checking that observation i maps to a valid leaf node
  } // closes loop over all observation
}


void fit_ensemble(std::vector<double> &fit, std::vector<tree> &t_vec, data_info &di){
  if(fit.size() != di.N) Rcpp::stop("[fit_ensemble]: size of fit must be equal to di.N!"); // honestly should never get triggered
  double* zz_cont = 0;
  int* zz_cat = 0;
  for(int i = 0; i < di.N; i++){
    if(di.z_cont != 0) zz_cont = di.z_cont + i * di.R_cont;
    if(di.z_cat != 0) zz_cat = di.z_cat + i * di.R_cat;
    fit[i] = 0.0;
    for(int m = 0; m < t_vec.size(); m++) fit[i] += t_vec[m].evaluate(zz_cont, zz_cat); // involves a tree traversal
  }
}

void compute_suff_stat_grow(suff_stat &orig_suff_stat, suff_stat &new_suff_stat, int &nx_nid, rule_t &rule, tree &t, data_info &di)
{
  double* zz_cont;
  int* zz_cat;
  double tmp_z = 0.0;
  int i;
  int l_count;
  int r_count;
  
  // we are growing tree from node nx, which has id of nx_nid
  
  int nxl_nid = 2*nx_nid; // id of proposed left child of nx
  int nxr_nid = 2*nx_nid+1; // id of proposed right child of nx
  
  suff_stat_it nx_it = orig_suff_stat.find(nx_nid); // iterator at element for nx in original sufficient statistic map
  new_suff_stat.clear();
  
  // copy orig_suff_stat into new_suff_stat
  for(suff_stat_it it = orig_suff_stat.begin(); it != orig_suff_stat.end(); ++it){
    new_suff_stat.insert(std::pair<int,std::vector<int>>(it->first, it->second));
  }
  
  // now we manipulate new_suff_stat to drop nx and add nxl and nxr
  new_suff_stat.insert(std::pair<int,std::vector<int>>(nxl_nid, std::vector<int>())); // create map element for left child of nx
  new_suff_stat.insert(std::pair<int,std::vector<int>>(nxr_nid, std::vector<int>())); // create map element for right child of nx
  new_suff_stat.erase(nx_nid); // remove map element for nx as it is not a bottom leaf node in new tree
  
  suff_stat_it nxl_it = new_suff_stat.find(nxl_nid); // iterator at element for nxl in new sufficient stat map
  suff_stat_it nxr_it = new_suff_stat.find(nxr_nid); // iterator at element for nxr in new sufficient stat map
  
  // loop over all of the observations that were assigned to nx in the original tree
  // nx_it->first is just the node id for nx (nx_nid)
  // nx_it->second is a vector that contains the indices of observations that land in nx
  
  for(int_it it = nx_it->second.begin(); it != nx_it->second.end(); ++it){
    i = *it;
    if(di.z_cont != 0) zz_cont = di.z_cont + i * di.R_cont;
    if(di.z_cat != 0) zz_cat = di.z_cat + i * di.R_cat;
    tmp_z = 0.0;
    
    if(rule.is_aa && !rule.is_cat){
      // axis-aligned rule
      if(zz_cont[rule.v_aa] < rule.c) nxl_it->second.push_back(i);
      else if(zz_cont[rule.v_aa] >= rule.c) nxr_it->second.push_back(i);
      else{
        Rcpp::Rcout << "  i = " << i << " v = " << rule.v_aa+1 << "  value = " << zz_cont[rule.v_aa] << " cutpoint = " << rule.c << std::endl;
        Rcpp::stop("[compute_ss_grow]: could not assign observation to left or right child in axis-aligned decision!");
      }
    } else if(!rule.is_aa && rule.is_cat){
      // categorical rule
      // we need to see whether i-th observation's value of the categorical pred goes to left or right
      // std::set.count returns 1 if the value is in the set and 0 otherwise
      l_count = rule.l_vals.count(zz_cat[rule.v_cat]);
      r_count = rule.r_vals.count(zz_cat[rule.v_cat]);
      
      if(l_count == 1 && r_count == 0) nxl_it->second.push_back(i);
      else if(l_count == 0 && r_count == 1) nxr_it->second.push_back(i);
      else if(l_count == 1 && r_count == 1) Rcpp::stop("[compute_ss_grow]: observation goes to both left and right child!");
      else{
        Rcpp::Rcout << "  i = " << i << " v = " << rule.v_cat+1 << "  value = " << zz_cat[rule.v_aa] << std::endl;
        Rcpp::Rcout << "left values:";
        for(set_it levels_it = rule.l_vals.begin(); levels_it != rule.l_vals.end(); ++levels_it) Rcpp::Rcout << " " << *levels_it;
        Rcpp::Rcout << std::endl;
        
        Rcpp::Rcout << "right values:";
        for(set_it levels_it = rule.r_vals.begin(); levels_it != rule.r_vals.end(); ++levels_it) Rcpp::Rcout << " " << *levels_it;
        Rcpp::Rcout << std::endl;
        
        Rcpp::stop("[compute_ss_grow]: could not assign observation to left or right child in categorical rule!");
      }
    } else if(!rule.is_aa && !rule.is_cat){
      // random combination rule
      tmp_z = 0.0;
      for(rc_it rcit = rule.rc_weight.begin(); rcit != rule.rc_weight.end(); ++rcit) tmp_z += (rcit->second) * zz_cont[rcit->first];
      if(tmp_z < rule.c) nxl_it->second.push_back(i);
      else if(tmp_z >= rule.c) nxr_it->second.push_back(i);
      else{
        Rcpp::Rcout << "  i = " << i << " tmp_z = " << tmp_z << " cutpoint = " << rule.c << std::endl;
        Rcpp::stop("[compute_ss_grow]: could not resolve random combination rule!");
      }
    } else{
      // we should never hit this part of the code
      Rcpp::stop("[compute_ss_grow]: cannot resolve the type of decision rule!");
    }
  } // closes loop over observations in the leaf
}

void compute_suff_stat_prune(suff_stat &orig_suff_stat, suff_stat &new_suff_stat, int &nl_nid, int &nr_nid, int &np_nid, tree &t, data_info &di)
{
  if(orig_suff_stat.count(nl_nid) != 1) Rcpp::stop("[compute_ss_prune]: did not find left node in suff stat map");
  if(orig_suff_stat.count(nr_nid) != 1) Rcpp::stop("[compute_ss_prune]: did not find right node in suff stat map");
  
  suff_stat_it nl_it = orig_suff_stat.find(nl_nid); // iterator at element for nl in original suff stat map
  suff_stat_it nr_it = orig_suff_stat.find(nr_nid); // iterator at element for nr in original suff stat map
  
  new_suff_stat.clear();
  // this makes a completely new copy of orig_suff_stat
  for(suff_stat_it ss_it = orig_suff_stat.begin(); ss_it != orig_suff_stat.end(); ++ss_it){
    new_suff_stat.insert(std::pair<int,std::vector<int>>(ss_it->first, ss_it->second));
  }
  new_suff_stat.insert(std::pair<int,std::vector<int>>(np_nid, std::vector<int>())); // add element for np in new suff stat map
  new_suff_stat.erase(nl_nid); // delete element for nl in new suff stat map since nl has been pruned
  new_suff_stat.erase(nr_nid); // delete element for nr in new suff stat map since nr has been pruned
  
  if(new_suff_stat.count(np_nid) != 1) Rcpp::stop("[compute_ss_prune]: didn't create element in new suff stat map for np correctly");
  suff_stat_it np_it = new_suff_stat.find(np_nid); // iterator at element for np in new suff stat map
  
  // time to populate np_it
  // first let's add the elements from nl_it
  for(int_it it = nl_it->second.begin(); it != nl_it->second.end(); ++it) np_it->second.push_back( *it );
  for(int_it it = nr_it->second.begin(); it != nr_it->second.end(); ++it) np_it->second.push_back( *it );
}


void compute_p_theta_ind(int &j, arma::mat &P, arma::vec &Theta, std::map<int,int> &leaf_map, suff_stat &ss, double &sigma, data_info &di, tree_prior_info &tree_pi)
{
  int L = ss.size(); // number of leafs
  P.set_size(L,L);
  P.eye();
  P.diag() *= 1.0/pow(tree_pi.tau, 2.0); // at this point P is just the diagonal matrix
  
  Theta.set_size(L);
  Theta.ones();
  Theta *= tree_pi.mu0/pow(tree_pi.tau, 2.0);
  // mu | all else ~ N(P^-1 Theta, P^-1)
  
  // the keys in ss are the node ids (nid) of the leaf node in each tree
  // we need to map those nid's to the set 0,..., L-1 so that we can populate the matrix P correctly
  // we need to save sum of x_itj^2 in each leaf
  leaf_map.clear(); // maps each leaf's nid to the correct row in P
  int l = 0;
  int i = 0;
  for(suff_stat_it l_it = ss.begin(); l_it != ss.end(); ++l_it){
    leaf_map.insert(std::pair<int, int>(l_it->first, l));
    ++l;
  }
  for(suff_stat_it l_it = ss.begin(); l_it != ss.end(); ++l_it){
    // l_it loops over all of the leafs in the tree
    l = leaf_map.find(l_it->first)->second; // gets the index in P and theta for the corresponding leaf
    for(int_it it = l_it->second.begin(); it != l_it->second.end(); ++it){
      // we need the sum of all x_itj's for this leaf
      i = *it;
      P(l,l) += 1.0/pow(sigma, 2.0) * pow(di.x[j + i * di.p], 2.0);
      Theta(l) += 1.0/pow(sigma, 2.0) * di.x[j + i * di.p] * di.rp[i];
    }
  }
}


void compute_p_theta_cs(int &j, arma::mat &P, arma::vec &Theta, std::map<int,int> &leaf_map, suff_stat &ss, double &rho, double &sigma, data_info &di, tree_prior_info &tree_pi)
{
  int L = ss.size();
  P.set_size(L,L);
  P.eye();
  P.diag() *= 1.0/pow(tree_pi.tau, 2.0); // at this point P is just the prior precision matrix of the jumps mu
  
  //arma::vec Theta = arma::ones<arma::vec>(L);
  //Theta *= tree.pi[j].mu0/pow(tree_pi[j].tau, 2.0); // at this point Theta just captures the contribution from the prior to the posterior mean of jumps mu
  
  Theta.set_size(L);
  Theta.zeros();
  
  
  // the keys in ss are the node ids (nid) of the leaf node in each tree
  // we need to map those nid's to the set 0,..., L-1 so that we can populate the matrix P correctly
  
  // we need to compute the following things:
  // 1. Sum of partial residuals (across all leafs) for each individual
  // 2. Sum of x_itj for each individual in each leaf
  // 3. Sum of x_itj * rp_i across all individuals in each leaf
  std::map<int, std::vector<double>> x_sums; // holds the sum of x_itj for each individual in each leaf
  std::map<int, double> xx_sums; // holds the sum of (x_itj)^2 across all individuals in each leaf
  std::map<int, double> xr_sums; // holds sum of (x_itj)*r_it across all individuals in each leaf
  std::vector<double> r_sums(di.n,0.0); // essentially obsolete since we now track each individuals sum of squared residuals in di.
  
  std::map<int,std::vector<double>>::iterator xs_it = x_sums.begin();
  std::map<int,std::vector<double>>::iterator xs_it_ll = x_sums.begin();
  std::map<int,double>::iterator xxs_it = xx_sums.begin();
  std::map<int, double>::iterator xrs_it = xr_sums.begin();
  
  leaf_map.clear();
  int l = 0;
  int ll = 0;
  int i;
  int subj_id;
  for(suff_stat_it l_it = ss.begin(); l_it != ss.end(); ++l_it){
    leaf_map.insert(std::pair<int, int>(l_it->first, l)); // we can now map leaf with node id l_it->first to the l-th row of P
    x_sums.insert(std::pair<int, std::vector<double>>(l_it->first, std::vector<double>(di.n, 0.0)));
    xx_sums.insert(std::pair<int, double>(l_it->first, 0.0)); // initialize total sum of x_ijt^2 in each leaf
    xr_sums.insert(std::pair<int, double>(l_it->first, 0.0)); // initialze total sum of x_itj*r_it in each leaf
    
    xs_it = x_sums.find(l_it->first);
    xxs_it = xx_sums.find(l_it->first);
    xrs_it = xr_sums.find(l_it->first);
  
    for(int_it it = l_it->second.begin(); it != l_it->second.end(); ++it){
      // it iterates over the indices of all observations that landed in leaf l_it
      i = *it;
      subj_id = di.subj_id[i]; // which subject contributed observation i to leaf l_it
      xs_it->second[subj_id] += di.x[j + i*di.p]; // increment running sum of x_ijt for this subject in this leaf
      xxs_it->second += pow(di.x[j + i*di.p], 2.0); // increment running total of x_ijt*x_ijt for ALL subjects in this leaf
      xrs_it->second += di.x[j + i*di.p] * di.rp[i]; // increment running total of x_ijt * r_it for ALL subjects in this leaf
      r_sums[subj_id] += di.rp[i]; // increment running sum of partial residual for this subject.
    }
    ++l; // incremember our counter. safer to write for loop with a paired iterator but that's messy...
  }
  
  // debugging time
  // for now these seem to be have been computed correctly
  /*
  for(suff_stat_it l_it = ss.begin(); l_it != ss.end(); ++l_it){
    Rcpp::Rcout << "leaf node " << l_it->first;
    Rcpp::Rcout << " xxs = " << xx_sums.find(l_it->first)->second;
    Rcpp::Rcout << "  xrs = " << xr_sums.find(l_it->first)->second << std::endl;
  }
  */
  
  // now that we have all of the leaf-wise summaries computed, we are ready to asemble P and Theta
  
  // we start with the diagonal elements of P
  for(suff_stat_it l_it = ss.begin(); l_it != ss.end(); ++l_it){
    l = leaf_map.find(l_it->first)->second;
    P(l,l) += 1.0/(1.0 - rho) * 1.0/pow(sigma, 2.0) * xx_sums.find(l_it->first)->second;
    Theta(l) += xr_sums.find(l_it->first)->second;
  }
  
  // we now need to loop over *all* subjects
  for(int subj_ix = 0; subj_ix < di.n; subj_ix++){
    
    for(suff_stat_it l_it = ss.begin(); l_it != ss.end(); ++l_it){
      l = leaf_map.find(l_it->first)->second;
      xs_it = x_sums.find(l_it->first);
      P(l,l) -= 1.0/pow(sigma,2.0) * rho/(1.0 - rho) * 1.0/(1 + rho * ( (double) di.ni[subj_ix] - 1.0)) * pow(xs_it->second[subj_ix],2.0);
      Theta(l) -= rho * r_sums[subj_ix]/(1 + rho * ( (double) di.ni[subj_ix] - 1.0)) * xs_it->second[subj_ix];
      
      // now we're ready to do the off-diagonal elements of P
      for(suff_stat_it ll_it = ss.begin(); ll_it != l_it; ++ll_it){
        ll = leaf_map.find(ll_it->first)->second;
        if(l ==ll) Rcpp::Rcout << "l = ll = " << l << std::endl;
        xs_it_ll = x_sums.find(ll_it->first);
        
        P(l,ll) -= 1.0/pow(sigma,2.0) * rho/(1.0 - rho) * 1.0/(1 + rho * ( (double) di.ni[subj_ix] - 1.0)) * xs_it->second[subj_ix] * xs_it_ll->second[subj_ix];
        P(ll,l) -= 1.0/pow(sigma,2.0) * rho/(1.0 - rho) * 1.0/(1 + rho * ( (double) di.ni[subj_ix] - 1.0)) * xs_it->second[subj_ix] * xs_it_ll->second[subj_ix];
      }
    }
  }
  
  // remember to rescale and add in the prior stuff for theta
  Theta *= 1.0/pow(sigma, 2.0) * 1.0/(1.0 - rho);
  Theta += tree_pi.mu0/pow(tree_pi.tau, 2.0);
  
}

double compute_lil(arma::mat &P, arma::vec &Theta, tree_prior_info &tree_pi)
{
  // -0.5 * log det P  - L(T) * log(tau) + 0.5 * Theta' P^-1 Theta - 0.5 * L * mu0^2/tau^2
  int L = Theta.size(); // number of leafs
  arma::vec Pinv_theta;
  bool solve_flag = arma::solve(Pinv_theta, P, Theta); // computes p^-1 * Theta;
  if(!solve_flag){
    Rcpp::Rcout << "[[compute_lil]]: could not compute P_inv" << std::endl;
    P.print();
    Rcpp::stop("Invalid posterior precision matrix for this tree.");
  }
  double lil = -1.0 * ( (double) L) * log(tree_pi.tau) - 0.5 * pow(tree_pi.mu0, 2.0)/pow(tree_pi.tau, 2.0) * ( (double) L);
  lil += 0.5 * arma::dot(Theta, Pinv_theta);
  lil -= 0.5 * arma::log_det_sympd(P);
  return lil;
}


void draw_mu(tree &t, arma::mat &P, arma::vec &Theta, std::map<int,int> &leaf_map, RNG &gen)
{
  tree::tree_p bn;
  arma::vec mu = gen.mvnormal(Theta, P); // remember mu ~ mvnormal(P^-1Theta, P^-1)
  for(std::map<int,int>::iterator it = leaf_map.begin(); it != leaf_map.end(); ++it){
    // it->first is the id of the leaf
    // we need to get the pointer to the leaf
    bn = t.get_ptr(it->first);
    if(bn == 0){
      Rcpp::Rcout << "Trying to assign mu for leaf with id " << it->first << " but could not find pointer to this node" << std::endl;
      t.print();
      Rcpp::stop("[draw_mu]: could not find node in the tree");
    } else{
      // it->second tells us which element of mu corresponds to this particular leaf.
      bn->set_mu(mu(it->second));
    }
  }
}

void draw_rule(rule_t &rule, tree &t, int &nid, data_info &di, tree_prior_info &tree_pi, RNG &gen){
  rule.clear(); // clear out the rule
  // we now are ready to draw a decision rule

  int rule_counter = 0; // we are allowed multiple tries to draw a valid random combination or categorical rule
  double c_upper = 1.0; // upper bound for range of cutpoints in axis aligned split
  double c_lower = -1.0; // lower bound for range of cutpoints in axis aligned split
  double tmp_weight = 0.0; // weights of random combination
  double c_max = 1.0; // upper bound for absolute value of cutpoint in random combination split
  tree::tree_p nx = t.get_ptr(nid); // at what node are we proposing this rule.
  
  double unif = gen.uniform();
  if( (!tree_pi.rc_split) || (tree_pi.rc_split && unif > *tree_pi.prob_rc) ){
    // either we were never allowed to try a rc split OR (1) we are allowed to try rc splits but (2) randomly choose a different type of rule
    //int v_raw = gen.multinomial(di.p, tree_pi.theta);
    int v_raw = gen.categorical(tree_pi.theta);
    if(v_raw < di.R_cont){
      // continuous variable so it's an axis-aligned splitting rule
      rule.is_aa = true;
      rule.is_cat = false;
      rule.v_aa = v_raw;
      if(tree_pi.unif_cuts[rule.v_aa] == 0){
        // draw the cutpoint from the supplied cutpoints
        c_lower = *(tree_pi.cutpoints->at(rule.v_aa).begin()); // returns smallest element in set
        c_upper = *(tree_pi.cutpoints->at(rule.v_aa).rbegin()); // reverse iterator, returns largest value in set
        nx->get_rg_aa(rule.v_aa, c_lower, c_upper);
        if(c_lower >= c_upper){
          // this is a weird tree and we'll just propose a trivial split
          c_lower = *(tree_pi.cutpoints->at(rule.v_aa).begin());
          c_upper = *(tree_pi.cutpoints->at(rule.v_aa).rbegin());
        }
        std::vector<double> valid_cutpoints;
        if(tree_pi.cutpoints->at(rule.v_aa).count(c_lower) != 1 || tree_pi.cutpoints->at(rule.v_aa).count(c_upper) != 1){
          // c_lower and c_upper were not found in the set of available cutpoints
          Rcpp::Rcout << "[grow tree]: attempting to select a cutpoint from given set" << std::endl;
          Rcpp::Rcout << "  lower bound is: " << c_lower << " count in set is " << tree_pi.cutpoints->at(rule.v_aa).count(c_lower) << std::endl;
          Rcpp::Rcout << "  upper bound is: " << c_upper << " count in set is " << tree_pi.cutpoints->at(rule.v_aa).count(c_upper) << std::endl;
          //Rcpp::Rcout << "  cutpoints are:";
          //for(std::set<double>::iterator it = tree_pi.cutpoints->at(rule.v_aa).begin(); it != tree_pi.cutpoints->at(rule.v_aa).end(); ++it) Rcpp::Rcout << " " << *it;
          //Rcpp::Rcout << std::endl;
          Rcpp::stop("we should never have a c that is outside the pre-defined set of cutpoints!");
        }
        // we want to draw from the cutpoints exclusive of c_lower & c_upper;
        // i.e. we want to start with the one just after c_lower and just before c_upper
        // std::set::lower_bound: iterator at first element that is not considered to come before
        // std::set::upper_bound: iterator at first element considered to come after
        // if value is not in set, lower_bound and upper_bound give same result
        // if value is in set: lower bound returns the value, upper bound returns the next value
        for(std::set<double>::iterator it = tree_pi.cutpoints->at(rule.v_aa).upper_bound(c_lower); it != tree_pi.cutpoints->at(rule.v_aa).lower_bound(c_upper); ++it){
          valid_cutpoints.push_back(*it);
        }
        int num_cutpoints = valid_cutpoints.size();
        if(num_cutpoints < 1){
          // no valid splits are available; we will just pick something, all of the observations will go to one child anyway...
          valid_cutpoints.clear();
          for(std::set<double>::iterator it = tree_pi.cutpoints->at(rule.v_aa).begin(); it != tree_pi.cutpoints->at(rule.v_aa).end(); ++it){
            valid_cutpoints.push_back(*it);
          }
          num_cutpoints = valid_cutpoints.size();
        }
        // at this point, valid cutpoints is a vector containing the available cutpoints at this node. we pick one uniformly.
        rule.c = valid_cutpoints[floor(gen.uniform() * num_cutpoints)];
        
      } else{
        // draw cutpoints uniformly
        c_upper = 1.0;
        c_lower = -1.0;
        nx->get_rg_aa(rule.v_aa, c_lower, c_upper);
        if(c_lower >= c_upper){
          c_lower = -1.0;
          c_upper = 1.0;
        }
        rule.c = gen.uniform(c_lower, c_upper);
      }
    } else{
      // categorical decision rule
      rule.is_aa = false;
      rule.is_cat = true;
      rule.v_cat = v_raw - di.R_cont;

      std::set<int> avail_levels = tree_pi.cat_levels->at(rule.v_cat); // get the full set of levels for this variable
      nx->get_rg_cat(rule.v_cat, avail_levels); // determine the set of levels available at nx.
      // if there is only one level left for this variable at nx, we will just propose a trivial split
      // and will reset the value of avail_levels to be the full set of all levels for the variable
      if(avail_levels.size() <= 1) avail_levels = tree_pi.cat_levels->at(rule.v_cat);
      
      rule.l_vals.clear();
      rule.r_vals.clear();
      
      if(tree_pi.graph_split[rule.v_cat] == 1 && tree_pi.edges->at(rule.v_cat).size() > 0){
        // if we explicitly say to use the graph to split the variables
        //graph_partition(avail_levels, rule.l_vals, rule.r_vals, di.adj_support->at(rule.v_cat), di.K->at(rule.v_cat), tree_pi.graph_cut_type, gen);
        graph_partition(avail_levels, rule.l_vals, rule.r_vals, tree_pi.edges->at(rule.v_cat), tree_pi.K->at(rule.v_cat), tree_pi.graph_cut_type, gen);
      } else{
        // otherwise we default to splitting the available levels uniformly at random: prob 0.5 to go to each child
        rule_counter = 0;
        while( ((rule.l_vals.size() == 0) || (rule.r_vals.size() == 0)) && rule_counter < 1000 ){
          rule.l_vals.clear();
          rule.r_vals.clear();
          for(set_it it = avail_levels.begin(); it != avail_levels.end(); ++it){
            if(gen.uniform() <= 0.5) rule.l_vals.insert(*it);
            else rule.r_vals.insert(*it);
          }
          ++(rule_counter);
        }
        if(rule_counter == 1000){
          Rcpp::stop("failed to generate valid categorical split in 1000 attempts"); // this should almost surely not get triggered.
        }
      }
      if( (rule.l_vals.size() == 0) || (rule.r_vals.size() == 0) ){
        Rcpp::Rcout << "[draw_rule]: proposed an invalid categorical rule" << std::endl;
        Rcpp::stop("proposed an invalid categorical rule!");
      }
    } // closes if/else determining whether we do an axis-aligned or categorical decision rule
  } else if(tree_pi.rc_split && unif <= *tree_pi.prob_rc){
    // random combination rule
    rule.is_aa = false;
    rule.is_cat = false;

    while( (rule.rc_weight.size() < 2) && (rule_counter < 1000) ){
      rule.rc_weight.clear();
      c_max = 0.0;
      for(int j = 0; j < di.R_cont; j++){
        if(gen.uniform() < (*tree_pi.theta_rc)){
          tmp_weight = gen.uniform(-1.0,1.0); // Breiman used Uniform(-1,1) weights and so shall we
          rule.rc_weight.insert(std::pair<int,double>(j,tmp_weight));
          c_max += fabs(tmp_weight);
        }
      }
      ++(rule_counter);
    }
    if(rule.rc_weight.size() < 2) Rcpp::stop("[propose_rule]: failed to generate a valid random combination rule in 1000 attempts!");
    else{
      rule.c = gen.uniform(-1.0,1.0) * c_max;
    }
  } else{
    Rcpp::stop("[draw_rule]: unable to draw a rule. check values of rc_split & prob_rc");
  } // closes all if/else determining the type of rule (axis-aligned or categorical) OR random combination
}


std::string write_tree(tree &t, tree_prior_info &tree_pi, set_str_conversion &set_str)
{
  std::ostringstream os;
  os.precision(32);
  
  tree::cnpv nds;
  rule_t rule;
  t.get_nodes(nds);
  
  for(tree::cnpv_it nd_it = nds.begin(); nd_it != nds.end(); ++nd_it){
    os << (*nd_it)->get_nid() << " ";
    if( (*nd_it)->get_ntype() == 'b'){
      // it's a bottom node
      os << "m " << (*nd_it)->get_mu() << " " << std::endl; // m tells us to expect mu next
    } else if( ((*nd_it)->get_ntype() == 't') && ( !(*nd_it)->l) ){ // because we need to look at left child of a node, make write_tree a friend in tree class
      // it's the top node and it has no children
      os << "m " << (*nd_it)->get_mu() << " " << std::endl; // tree is a stump, m tells us to expect mu next
    } else{
      // we need to print out the rule
      //os << "make rule " << std::endl;
      os << "r "; // node has a decision rule. r tells us to expect a rule next
      rule.clear();
      rule = (*nd_it)->get_rule();
      
      os << rule.is_aa << " " << rule.is_cat << " ";

      if(rule.is_aa && !rule.is_cat){
        // axis-aligned rule
        os << rule.c << " " << rule.v_aa;
      } else if(!rule.is_aa && rule.is_cat){
        // categorical rule
        int K = tree_pi.K->at(rule.v_cat); // how many levels
        os << rule.v_cat << " " << K << " ";
        os << set_str.set_to_hex(K, rule.l_vals) << " ";
        os << set_str.set_to_hex(K, rule.r_vals) << " ";
      } else if(!rule.is_aa && !rule.is_cat){
        // random combination
        os << rule.c << " ";
        for(rc_it rcit = rule.rc_weight.begin(); rcit != rule.rc_weight.end(); ++rcit){
          os << rcit->first << " " << rcit->second << " ";
        }
      } else{
        Rcpp::stop("[write tree]: rule cannot be both axis-aligned and categorical!");
      }
      os << std::endl;
    } // closes if/else checking what type of node we are writing
  } // closes loop over the nodes in t
  
  return os.str();
  
}

void read_tree(tree &t, std::string &tree_string, set_str_conversion &set_str)
{
  std::istringstream tree_ss(tree_string); // an in stringstream of the tree's string representation
  std::string node_string; // string for each individual node in the tree
  
  int nid;
  char stream_type; // either 'm' to indicate that next element in stream is mu or 'r' to indicate a rule follows
  
  double tmp_mu; // holds the value of mu for a leaf node

  char aa; // '0' or '1' for rule.is_aa
  char cat; // '0' or '1' for rule.is_cat
  rule_t tmp_rule; // temporary rule that gets populated as we read the tree's string/stream
  
  int tmp_v; // for reading in rc weights
  double tmp_phi; // for reading in rc weights
  
  int K; // tells us how many levels there were to the categorical variable
  std::string l_hex; // string representation of the l_vals in a categorical rule
  std::string r_hex; // string representation of the r_vals in a categorical rule
  
  std::map<int, rule_t> decision_nodes;
  std::map<int, double> leaf_nodes;
  
  while(tree_ss){
    std::getline(tree_ss, node_string, '\n');
    if(node_string.size() > 0){
      std::istringstream node_ss(node_string); // in stream for the single node
      node_ss >> nid; // get the node nid
      node_ss >> stream_type;
      
      if(stream_type == 'm'){
        node_ss >> tmp_mu;
        leaf_nodes.insert(std::pair<int,double>(nid, tmp_mu));
      } else if(stream_type == 'r'){
        tmp_rule.clear();
        node_ss >> aa;
        node_ss >> cat;
        
        if(aa == '0') tmp_rule.is_aa = false;
        else tmp_rule.is_aa = true;
        
        if(cat == '0') tmp_rule.is_cat = false;
        else tmp_rule.is_cat = true;

        
        if(tmp_rule.is_aa && !tmp_rule.is_cat){
          // axis-aligned
          node_ss >> tmp_rule.c;
          node_ss >> tmp_rule.v_aa;
        } else if(!tmp_rule.is_aa && tmp_rule.is_cat){
          // categorical rule
          node_ss >> tmp_rule.v_cat; // get the variable index
          node_ss >> K; // we now know how many levels of the categorical variable there were
          node_ss >> l_hex;
          node_ss >> r_hex;
          
          tmp_rule.l_vals = set_str.hex_to_set(K, l_hex);
          tmp_rule.r_vals = set_str.hex_to_set(K, r_hex);
        } else if(!tmp_rule.is_aa && !tmp_rule.is_cat){
          // random combination
          node_ss >> tmp_rule.c; // get the cutpoint first
          while(node_ss){
            node_ss >> tmp_v;
            node_ss >> tmp_phi;
            tmp_rule.rc_weight.insert(std::pair<int,double>(tmp_v, tmp_phi));
          }
        } else{
          Rcpp::stop("[read tree]: rule cannot be axis-aligned and categorical");
        }
        decision_nodes.insert(std::pair<int, rule_t>(nid, tmp_rule));
      } // closes if/else checking what type of node we're parsing
    } // closes if checking that we found a valid node in the stream
  } // closes while that parses stream for the tree
  
  // we now have decision_nodes and leaf_nodes and are ready to build up our tree
  t.to_null(); // clear out the tree if there was anything there
  
  // remember std::map is sorted by key.
  // we have always used node id as the key so by iterating over our map, we will *never*
  // attempt to birth from a node that has not already been created.

  for(std::map<int,rule_t>::iterator it = decision_nodes.begin(); it != decision_nodes.end(); ++it){
    t.birth(it->first, it->second); // do the birth.
  }

  tree::npv bnv;
  t.get_bots(bnv); // get the bottom nodes
  std::map<int,double>::iterator leaf_it;

  for(tree::npv_it it = bnv.begin(); it != bnv.end(); ++it){
    leaf_it = leaf_nodes.find( (*it)->get_nid() );
    if(leaf_it == leaf_nodes.end()){
      Rcpp::Rcout << "[read_tree]: we didn't read a leaf node with nid" << (*it)->get_nid() << std::endl;
      Rcpp::stop("mistake in reading in tree");
    } else{
      (*it)->set_mu(leaf_it->second);
    }
  }
  
}

// our graph algorithms defined here

// sometimes it's convenient to have a redundant (i.e. symmetric) representation of the edges
// each edge represented in two elements of the map with keys correspond to the "source" and "sink"
// will call this specifically in boruvka's and in the LERW style code
void build_symmetric_edge_map(edge_map &emap, std::vector<edge> &edges, std::set<int> &vertices)
{
  emap.clear();
  
  // create an element with a key for each vertex in the specified set
  for(std::set<int>::iterator v_it = vertices.begin(); v_it != vertices.end(); ++v_it) emap.insert(std::pair<int, std::vector<edge>>(*v_it, std::vector<edge>()));
  
  for(edge_vec_it it = edges.begin(); it != edges.end(); ++it){
    
    if(emap.count(it->source) != 1 || emap.count(it->sink) != 1){
      Rcpp::Rcout << "[build_symmetric_edge_map]: Supplied set includes an edge with one vertex outside set of supplied vertices!" << std::endl;
      Rcpp::Rcout << "  edge is from " << it->source << " to " << it->sink << std::endl;
      Rcpp::Rcout << "  Supplied vertices:";
      for(std::set<int>::iterator v_it = vertices.begin(); v_it != vertices.end(); ++v_it) Rcpp::Rcout << " " << *v_it;
      Rcpp::Rcout << std::endl;
    } else{
      emap.find(it->source)->second.push_back(edge(it->source, it->sink, it->weight));
      emap.find(it->sink)->second.push_back(edge(it->sink, it->source, it->weight));
    }
  } // finish looping over all of the edges
}


// in the main BART loop, we will often have a set of categories that form a subset of the vertices of a network
// we will have to partition the subgraph induced by this set of vertices
// get_induced_edges loops over all edges and pulls out those whose source & sink are in the vertex_subset
// we need to pull out the edges
std::vector<edge> get_induced_edges(std::vector<edge> &edges, std::set<int> &vertex_subset)
{
  std::vector<edge> subset_edges;
  // loop over the edges
  for(edge_vec_it it = edges.begin(); it != edges.end(); ++it){
    // check that both the source and the sink belong to the vertex_set
    if(vertex_subset.count(it->source) == 1 && vertex_subset.count(it->sink) == 1){
      subset_edges.push_back(edge(it->source, it->sink, it->weight));
    }
  }
  return subset_edges;
}

// stuff for finding connected components:
//  depth-first search: we maintain a map with keys corresponding to vertices and boolean values; this records which vertices our DFS has visited
//  dfs should only ever be called using a *symmetric* edge_map; we need a key for every vertex
void dfs(int v, std::map<int, bool> &visited, std::vector<int> &comp, edge_map &emap)
{
  // this is the first time we have visited vertex labelled v
  std::map<int,bool>::iterator v_it = visited.find(v);
  v_it->second = true; // mark that we have visited vertex v
  comp.push_back(v); // now that v has been marked, we can add it to the connected component

  // now find the element of our map containing all verties whose source is v
  if(emap.count(v) != 1){
    // if we are calling dfs correctly this should *never* be hit
    Rcpp::Rcout << "[dfs]: at vertex" << v << std::endl;
    Rcpp::Rcout << "keys in emap are:";
    for(edge_map_it tmp_it = emap.begin(); tmp_it != emap.end(); ++tmp_it) Rcpp::Rcout << " " << tmp_it->first << std::endl;
    Rcpp::stop("something is wrong with our edge map");
  } else{
    edge_map_it v_edges_it = emap.find(v);
    for(edge_vec_it it = v_edges_it->second.begin(); it != v_edges_it->second.end(); ++it){
      // it points to an edge leaving v, let us see if the sink has been visited
      int vv = it->sink;
      std::map<int,bool>::iterator vv_it = visited.find(vv);
      if(vv_it == visited.end()){
        // we were unable to find an entry for vv in visited.
        // if we are calling dfs correctly, this should *never* be hit
        Rcpp::Rcout<< "[dfs]: Traversing edges from " << v << " to " << vv << std::endl;
        Rcpp::Rcout << "  did not find entry for " << vv << " in visited" << std::endl;
        Rcpp::Rcout << "  keys of visited:";
        for(std::map<int,bool>::iterator visit_it = visited.begin(); visit_it != visited.end(); ++visit_it) Rcpp::Rcout << " " << visit_it->first;
        Rcpp::Rcout << std::endl;
        Rcpp::stop("something is wrong with our edge map or visited!");
      } else{
        if(!vv_it->second){
          // we have not yet visited the vertex vv and so we continue our dfs from there
          dfs(vv, visited, comp, emap);
        }
      }
    } // closes loop over edges incident to v
  } // closes if/else checking that we have an entry for v in our edge_map
}


void find_components(std::vector<std::vector<int> > &components, std::vector<edge> &edges, std::set<int> &vertices)
{
  components.clear();
  
  if(edges.size() == 0){
    // no edges: every vertex is is its own component
    //Rcpp::Rcout << "[find_components]: no edges" << std::endl;
    for(std::set<int>::iterator v_it = vertices.begin(); v_it != vertices.end(); ++v_it) components.push_back(std::vector<int>(1, *v_it));
  } else{
    std::map<int, bool> visited;
    // initially mark every vertex as unvisited
    for(std::set<int>::iterator v_it = vertices.begin(); v_it != vertices.end(); ++v_it) visited.insert(std::pair<int,bool>(*v_it, false));
    edge_map emap;
    build_symmetric_edge_map(emap, edges, vertices);
    
    // loop over the vertices and run dfs if we haven't already visited that vertex
    for(std::set<int>::iterator v_it = vertices.begin(); v_it != vertices.end(); ++v_it){
      //Rcpp::Rcout << " starting from " << *v_it << std::endl;
      if(!visited.find(*v_it)->second){
        // we have not yet visited vertex labelled *v_it, so it represents a new component
        std::vector<int> new_comp;
        dfs(*v_it, visited, new_comp, emap);
        // by the time the dfs finishes, we have visited everything in the current connected component
        components.push_back(new_comp); // add current connected component to our vector of components
      }
    }
  }
}

void get_unique_edges(std::vector<edge> &edges)
{
  std::vector<edge> unik_edges;
  if(edges.size() > 0){
    unik_edges.push_back(edges[0]);
    for(std::vector<edge>::iterator it = edges.begin(); it != edges.end(); ++it){
      bool add_edge = true;
      for(std::vector<edge>::iterator uit = unik_edges.begin(); uit != unik_edges.end(); ++uit){
        if( (it->source == uit->source && it->sink == uit->sink) || (it->sink == uit->source && it->source == uit->sink)){
          // *it (or its reverse) is already in our set of unique edges. DO NOT keep comparing *it to elements of unik_edges & DO NOT add it to unik_edges
          //Rcpp::Rcout << "duplicate edge found!" << std::endl;
          add_edge = false;
          break;
        }
      }
      if(add_edge) unik_edges.push_back(*it);
    }
  }
  edges.clear();
  for(std::vector<edge>::iterator it = unik_edges.begin(); it != unik_edges.end(); ++it) edges.push_back(edge(it->source, it->sink, it->weight));
}

void boruvka(std::vector<edge> &mst_edges, std::vector<edge> &edges, std::set<int> &vertices)
{
  
  mst_edges.clear();
  
  edge_map emap;
  build_symmetric_edge_map(emap, edges, vertices);

  std::vector<std::vector<int> > mst_components;
  find_components(mst_components, mst_edges, vertices);
  int n_vertex = vertices.size();
  int counter = 0;
  // Boruvka has worst-case complexity of O(log(n)), so we add a considerable buffer
  while( mst_components.size() > 1 && counter < 2*n_vertex ){
    std::vector<edge> new_edges; // the new edges that will get added to the MST
    //Rcpp::Rcout << " Starting Round " << counter << " of Boruvka. Edges are:" << std::endl;
    //for(std::vector<edge>::iterator mst_eit = mst_edges.begin(); mst_eit != mst_edges.end(); ++mst_eit){
    //  Rcpp::Rcout << mst_eit->source << " to " << mst_eit->sink << std::endl;
    //}
    for(std::vector<std::vector<int> >::iterator comp_it = mst_components.begin(); comp_it != mst_components.end(); ++comp_it){
      // looping over all of the existing components in the current MST
      int min_source = 0;
      int min_sink = 0;
      double min_weight = 1.0;
      for(std::vector<int>::iterator v_it = comp_it->begin(); v_it != comp_it->end(); ++v_it){
        // v_it points to a particular vertex in the component *comp_it.
        // We need to loop over all incident edges to *v_it
        // we check whether (A) the "sink" is outside the component and (B) whether it has smallest weight
        
        //Rcpp::Rcout << "    Visiting vertex " << *v_it << std::endl;
        
        edge_map_it ve_it = emap.find(*v_it); // ve_it points to the vector of edges incident to *v_it
        for(std::vector<edge>::iterator e_it = ve_it->second.begin(); e_it != ve_it->second.end(); ++e_it){
          // e_it points to a specific edge
          // we check whether (A) the "sink" is outside the component and (B) whether it has smallest weight
          //Rcpp::Rcout << "   checking edge from " << e_it->source << " to " << e_it->sink << " count = " << std::count(comp_it->begin(), comp_it->end(), e_it->sink) << std::endl;
          if(std::count(comp_it->begin(), comp_it->end(), e_it->sink) == 0 && e_it->weight < min_weight){
            // found a new minimum edge weight leaving the component!
            min_source = e_it->source; // this had better be *v_it
            min_sink = e_it->sink;
            min_weight = e_it->weight;
          }
        } // closes loop over edge incident to vertex *v_it
      } // closes loop over vertices in component *comp_it
      new_edges.push_back(edge(min_source, min_sink, min_weight));
    } // closes loop over components
    
    // at this point, we've finished looping over all of the vertices in a particular component and we have found the edge
    // that leaves the component and has minimum weight. we dump that edge into our running collection of edges in the MSt
    
    for(std::vector<edge>::iterator it = new_edges.begin(); it != new_edges.end(); ++it){
      mst_edges.push_back(*it);
    }
    get_unique_edges(mst_edges); // we have may duplicate edges so kill them off here.
    find_components(mst_components, mst_edges, vertices); // re-compute the number of components
    counter++;
  } // closes main while loop
  
  if(mst_components.size() > 1){
    Rcpp::Rcout << "[boruvka]: after " << counter << " rounds, we have not yet formed a single connected MST" << std::endl;
    Rcpp::Rcout << "  returning an empty set of edges" << std::endl;
    mst_edges.clear();
  }
}
void wilson(std::vector<edge> &mst_edges, std::vector<edge> &edges, std::set<int> &vertices, RNG &gen)
{
  int n_vertex = vertices.size();
  edge_map emap;
  build_symmetric_edge_map(emap, edges, vertices);
  
  
  std::map<int, bool> in_tree; // key is vertex label, value is boolean of whether vertex is in tree
  std::vector<int> possible_roots; // we have to start the tree from somewhere

  for(std::set<int>::iterator v_it = vertices.begin(); v_it != vertices.end(); ++v_it){
    in_tree.insert(std::pair<int, bool>(*v_it, false)); // each vertex is not tree
    possible_roots.push_back(*v_it);
  }
  int root = possible_roots[floor(gen.uniform() * n_vertex)];
  in_tree.find(root)->second = true; // mark the root as a member of the tree
  
  
  bool tree_full = true; // true when all vertices are members of the tree
  
  // when we start, we know most values in in_tree are false so this loop isn't super expensive
  for(std::set<int>::iterator v_it = vertices.begin(); v_it != vertices.end(); ++v_it){
    if(!in_tree.find(*v_it)->second){
      tree_full = false;
      break;
    }
  }
  
  possible_roots.clear();
  // from now on, possible_roots will store the potential roots for the LERW to the tree
  int outer_counter = 0; // number of rounds in which we try to walk to the tree
  
  while(!tree_full && outer_counter <= n_vertex){
    // each time through the loop, we choose a new vertex and walk from it to the tree
    // so long as the inner loop terminates (i.e. we successfully walk to the tree)
    // we should never have to run the outer loop more than n_vertex times
    
    possible_roots.clear();
    for(std::set<int>::iterator v_it = vertices.begin(); v_it != vertices.end(); ++v_it){
      if(!in_tree.find(*v_it)->second) possible_roots.push_back(*v_it);
    }
    int next_state = possible_roots[floor(gen.uniform() * possible_roots.size())]; // next_state is start of the walk
    int old_state;
    
    std::vector<int> lew;
    std::vector<edge> lew_edges;
    // as we do the loop erasure, we have to keep track of which vertices have been visited
    // key is vertex label and value is boolean for whether we visited or not
    std::map<int, bool> visited;
    for(std::set<int>::iterator v_it = vertices.begin(); v_it != vertices.end(); ++v_it){
      visited.insert(std::pair<int,bool>(*v_it, false));
    }
    
    lew.push_back(next_state); // only used for testing
    visited.find(next_state)->second = true; // mark the root as having been visited
    
    
    int inner_counter = 0;
    int edge_index;
    
    // next few lines only for testing purposes
    //Rcpp::Rcout << "Starting Round " << outer_counter << "! Tree contains:";
    //for(std::map<int,bool>::iterator tr_it = in_tree.begin(); tr_it != in_tree.end(); ++tr_it){
    //  if(tr_it->second) Rcpp::Rcout << " " << tr_it->first;
    //}
    //Rcpp::Rcout << std::endl;
    // done with the printing used to test
    
    while(!in_tree.find(next_state)->second && inner_counter < (int) 10 * pow(n_vertex,3.0)){
      old_state = next_state;
      edge_map_it em_it = emap.find(old_state); // em_it now points to vector of edges incident to old_state
      if(em_it->second.size() == 0){
        //Rcpp::Rcout << "Random walk is a vertex that has no incident edges. Stopping now" << std::endl;
        mst_edges.clear(); // clear out mst_edges
        tree_full = true;
        break; // break out of inner loop
      } else{
        std::vector<double> weights;
        for(std::vector<edge>::iterator e_it = em_it->second.begin(); e_it != em_it->second.end(); ++e_it){
          weights.push_back(e_it->weight);
        }
        edge_index = gen.categorical(weights);
        
        next_state = em_it->second[edge_index].sink;
        if(!visited.find(next_state)->second){
          // we have not yet visited the vertex next_state
          lew.push_back(next_state); // only for testing purposes
          lew_edges.push_back( em_it->second[edge_index] ); // add the edge to next_state
          visited.find(next_state)->second = true; // mark next_state has having been visited
        } else{
          // we have already visited next_state
          std::vector<int>::iterator prev_it = lew.end()-1; // prev_it points to the last current state in the LERW
          std::vector<edge>::iterator prev_eit = lew_edges.end() - 1; // points to last edge in the LERW
          for(;prev_it != lew.begin(); --prev_it){
            if(*prev_it == next_state) break;
          }
          // prev_it should point to the last instance of next_state in lew
          for(;prev_eit != lew_edges.begin(); --prev_eit){
            if(prev_eit->sink == next_state) break;
          }
          // prev_eit should now point to the last edge in lew that goes into next_state
          
          if(prev_it == lew.begin() && lew[0] != next_state){
            // if we've reached the beginning of lew but haven't seen next_state something has gone quite wrong
            Rcpp::Rcout << "[wilson]: LERW trying to re-visit " << next_state << " but could not find last visit" << std::endl;
            Rcpp::Rcout << "  LERW currently contains:";
            for(std::vector<int>::iterator lew_it = lew.begin(); lew_it != lew.end(); ++lew_it) Rcpp::Rcout << " " << *lew_it;
            Rcpp::Rcout << std::endl;
            Rcpp::stop("could not perform loop-erasure");
            
          } else if(prev_it == lew.begin() && lew[0] == next_state){
            // we are attempting to revisit the the starting point of the LERW
            // we reset lew to contain just next_state and empty lew_edges
            lew.clear();
            lew.push_back(next_state);
            lew_edges.clear();
            
            // we need to update visited: in this case, only next_state should be marked as visited
            for(std::map<int, bool>::iterator visit_it = visited.begin(); visit_it != visited.end(); ++visit_it){
              if(visit_it->first != next_state) visit_it->second = false;
              else visit_it->second = true;
            }
          } else{
            // we have found the last time the LERW has visited next state and it wasn't at the beginning
            lew.erase(prev_it, lew.end()); // erase also kills the last visit to next_state...
            lew.push_back(next_state); // ...so we add that back in!
            
            // prev_eit points to the last edge that goes INTO next_state
            // this is not an edge in the cycle being popped
            int prev_source = prev_eit->source;
            int prev_sink = prev_eit->sink; // had better be equal to next_state!
            double prev_weight = prev_eit->weight;
            
            lew_edges.erase(prev_eit, lew_edges.end()); // not really sure what prev_eit will point to after the erase...
            lew_edges.push_back(edge(prev_source, prev_sink, prev_weight)); // ... so we make a new copy of the edge to be safe
            
            // we must update visited.
            // looping over lew is not sufficient because we may have removed a visited vertex while popping a cycle
            for(std::map<int, bool>::iterator visit_it = visited.begin(); visit_it != visited.end(); ++visit_it){
              if(std::count(lew.begin(), lew.end(), visit_it->first) == 1) visit_it->second = true;
              else if(std::count(lew.begin(), lew.end(), visit_it->first) == 0) visit_it->second = false;
              else{
                Rcpp::Rcout << "[wilson]: After popping cycle, our LERW visits a vertex twice!" << std::endl;
                Rcpp::Rcout << "  LERW currently contains:";
                for(std::vector<int>::iterator lew_it = lew.begin(); lew_it != lew.end(); ++lew_it) Rcpp::Rcout << " " << *lew_it;
                Rcpp::Rcout << std::endl;
                Rcpp::stop("LERW somehome contains a loop...");
              }
            } // closes loop updating visited after cycle popping
          } // closes if/else checking when we visited next_state
        } // closes if/else checking whether we have visited next_state (and popping cycle if so)
      } // closes if/else checking that old_state has incident edges
      inner_counter++;
    } // closes inner while loop
    
    if(inner_counter == (int) 10*pow(n_vertex,3.0)){
      Rcpp::Rcout << "[wilson]: inner loop did not hit tree in " << 10*pow(n_vertex,3.0) <<  " steps. stopping now" << std::endl;
      Rcpp::Rcout << "  LERW currently contains :";
      for(std::vector<int>::iterator lew_it = lew.begin(); lew_it != lew.end(); ++lew_it) Rcpp::Rcout << " " << *lew_it;
      Rcpp::Rcout << std::endl;
      Rcpp::stop("may need to increase number of steps allowed for inner loop!");
    } else{
      // dump edges from lew into mst_edges
      for(std::vector<edge>::iterator e_it = lew_edges.begin(); e_it != lew_edges.end(); ++e_it) mst_edges.push_back(*e_it);
      // update in_tree
      for(std::vector<int>::iterator lew_it = lew.begin(); lew_it != lew.end(); ++lew_it) in_tree.find(*lew_it)->second = true;
      
      tree_full = true;
      for(std::map<int, bool>::iterator tr_it = in_tree.begin(); tr_it != in_tree.end(); ++tr_it){
        if(!tr_it->second){
          tree_full = false;
          break;
        }
      }
      
      // some printing only used for testing
      //Rcpp::Rcout << "Ending Round " << outer_counter << "! Tree contains:";
      //for(std::map<int,bool>::iterator tr_it = in_tree.begin(); tr_it != in_tree.end(); ++tr_it){
      //  if(tr_it->second) Rcpp::Rcout << " " << tr_it->first;
      //}
      //Rcpp::Rcout << std::endl;
    }
    outer_counter++;
  } // closes outer loop
  
  if(!tree_full){
    Rcpp::Rcout << "[wilson]: unable to reach every vertex in the graph. returning an empty set of edges" << std::endl;
    mst_edges.clear();
  }
}
void graph_partition(std::set<int> &avail_levels, std::set<int> &l_vals, std::set<int> &r_vals, std::vector<edge> &orig_edges, int &K, int &cut_type, RNG &gen)
{
  // get the edges for the induced subgraph
  std::vector<edge> edges = get_induced_edges(orig_edges, avail_levels);
  
  // check if the induced subgraph is connected
  std::vector<std::vector<int> > components;
  find_components(components, edges, avail_levels);
  if(components.size() == 0){
    // this should never be reached
    Rcpp::stop("[graph_partition]: graph has no components...");
  } else if(components.size() > 1){
    // graph is not connected. assign whole components to left and right children uniformly at random
    l_vals.clear();
    r_vals.clear();
    int rule_counter = 0;
    while( ((l_vals.size() == 0 || r_vals.size() == 0)) && rule_counter < 1000 ){
      l_vals.clear();
      r_vals.clear();
      for(int comp_ix = 0; comp_ix < components.size(); comp_ix++){
        if(gen.uniform() <= 0.5){
          // send everything in this component to the left child
          for(int_it it = components[comp_ix].begin(); it != components[comp_ix].end(); ++it) l_vals.insert(*it);
        } else{
          // send everything in this component to the right child
          for(int_it it = components[comp_ix].begin(); it != components[comp_ix].end(); ++it) r_vals.insert(*it);
        }
      }
      ++rule_counter;
    }
    if(rule_counter == 1000){
      Rcpp::stop("[graph partition]: graph disconnected. failed to generate a non-trivial partiton of components in 1000 attempts!");
    }
  } else{
    // induced subgraph is connected and we can form a partition
    if(cut_type <= 1){
      std::vector<edge> mst_edges;
      if(cut_type == 0) wilson(mst_edges, edges, avail_levels, gen);
      else{
        // need to assign random weights to each edge
        for(std::vector<edge>::iterator e_it = edges.begin(); e_it != edges.end(); ++e_it) e_it->weight = gen.uniform();
        boruvka(mst_edges, edges, avail_levels);
      }
      
      if(mst_edges.size() == 0){
        // something failed in Wilson or Boruvka algorithm, we should stop
        Rcpp::stop("[graph_partition]: unable to compute the random spanning tree. Quitting now!");
      } else{
        int cut_index = 0;
        std::vector<double> cut_probs;
        for(std::vector<edge>::iterator e_it = mst_edges.begin(); e_it != mst_edges.end(); ++e_it){
          cut_probs.push_back(e_it->weight);
        }
        cut_index = gen.categorical(cut_probs);
        std::vector<edge> cut_mst_edges;
        for(int e = 0; e < mst_edges.size(); e++){
          if(e != cut_index) cut_mst_edges.push_back(mst_edges[e]);
        }
        std::vector<std::vector<int> > cut_components;
        find_components(cut_components, cut_mst_edges, avail_levels);
        
        if(cut_components.size() != 2){
          Rcpp::Rcout << "[graph_partition]: Attempted to partition connected subgraph by cutting a random spanning tree" << std::endl;
          Rcpp::Rcout << "   but after deleting edge from the spanning tree, resulting graph does not have 2 components. something is wrong" << std::endl;
          Rcpp::stop("error in cutting edge from spanning tree");
        } else{
          l_vals.clear();
          r_vals.clear();
          for(std::vector<int>::iterator it = cut_components[0].begin(); it != cut_components[0].end(); ++it) l_vals.insert(*it);
          for(std::vector<int>::iterator it = cut_components[1].begin(); it != cut_components[1].end(); ++it) r_vals.insert(*it);
        } // closes if/else checking that we have 2 connected components after deleting edge from spanning tree
      } // closes if/else checking that our spanning tree search was successful
    } else{
      Rcpp::Rcout << "[graph_partition]: cut_type argument is " << cut_type << std::endl;
      Rcpp::stop("cut_type argument must be 0 or 1");
    } // cloeses if/else checking the cut type
  } // closes if/else checking whether the graph induced by avail_levels is connected or not
}

void update_theta_u(std::vector<double> &theta, double &u, std::vector<int> &var_count, int &R, double &a_u, double &b_u, RNG &gen)
{
  if(theta.size() != R){
    Rcpp::Rcout << "theta has size " << theta.size() << "  R = " << R << std::endl;
    Rcpp::stop("theta must have size R!");
  } else{
    double tmp_sum = 0.0;
    double tmp_concentration = 0.0;
    double sum_log_theta = 0.0;
    int v_count;
    std::vector<double> tmp_gamma(R, 0.0);
    
    // update theta first
    double u_orig = u;
    for(int j = 0; j < R; j++){
      v_count = var_count[j];
      tmp_concentration = u_orig/(1.0 - u_orig) + (double) v_count;
      tmp_gamma[j] = gen.gamma(tmp_concentration, 1.0);
      tmp_sum += tmp_gamma[j];
    }
    for(int j = 0; j < R; j++){
      theta[j] = tmp_gamma[j]/tmp_sum;
      sum_log_theta += log(theta[j]);
    }
    
    // we're now ready to update u
    double u_prop = gen.beta(a_u,b_u);
    double log_like_prop = (u_prop)/(1.0 - u_prop) * sum_log_theta;
    double log_like_orig = (u_orig)/(1.0 - u_orig) * sum_log_theta;
    
    log_like_prop += lgamma( (double) R * u_prop/(1.0 - u_prop)) - ((double) R) * lgamma(u_prop/(1.0 - u_prop));
    log_like_orig += lgamma( (double) R * u_orig/(1.0 - u_orig)) - ((double) R) * lgamma(u_orig/(1.0 - u_orig));
    double log_accept = log_like_prop - log_like_orig;
    if(gen.log_uniform() <= log_accept) u = u_prop;
    else u = u_orig;
  }
}
void update_theta_rc(double& theta_rc, int &rc_var_count, int &rc_rule_count, double &a_rc, double &b_rc, int &R_cont, RNG &gen)
{
  // since we disallow rc rules with 0 or 1 variable, the posterior is proportional to
  // theta^(a + rc_var_count - 1) * (1 - theta)^(b + R_cont * rc_rule_count - rc_var_count - 1)/(1 - (1- theta)^R_cont - R_cont * theta * 1 - theta)^(R_cont - 1))
  // note that the denominator here is just P(Bin(R_cont, theta) >= 2) as a function of theta
  
  // use independence MH: transition proposal is Beta(a  + rc_var_count, b + R_cont * rc_rule_count)
  // acceptance ratio turns out to:
  // P( Bin(R_cont, theta_orig) >= 2) / P(Bin(R_cont, theta_prop))
  
  double a_post = a_rc + (double) rc_var_count;
  double b_post = b_rc + ((double) R_cont) * rc_rule_count - (double) rc_var_count;
  
  double theta_orig = theta_rc;
  double theta_prop = gen.beta(a_post, b_post);
  
  double log_post_prop = log(1.0 - pow(1.0 - theta_prop, R_cont) - (double) R_cont * theta_prop * pow(1.0 - theta_prop, R_cont-1));
  double log_post_orig = log(1.0 - pow(1.0 - theta_orig, R_cont) - (double) R_cont * theta_orig * pow(1.0 - theta_orig, R_cont-1));
 
  double log_alpha = log_post_orig - log_post_prop;
  if(gen.log_uniform() < log_alpha) theta_rc = theta_prop;
  else theta_rc = theta_orig;
}

void update_rho(double &rho, double &sigma, data_info &di, RNG &gen)
{
  // we have current value of rho, we propose a new one, and then take an MH step
  // easier to work with xi = log(rho/(1-rho))
  
  double rho_curr = rho;
  double xi_curr = log(rho_curr/(1.0 - rho_curr));
  double xi_prop = xi_curr + gen.normal(0.0,0.1); // add N(0,0.1^2) noise to xi
  double rho_prop = 1.0/(1.0 + exp(-1.0 * xi_prop)); // new value of rho
  
  double log_prior_curr = 0.0; // assumes a uniform prior on rho
  double log_prior_prop = 0.0;
  double log_post_curr = log(rho_curr) + log(1.0 - rho_curr) + log_prior_curr;
  double log_post_prop = log(rho_prop) + log(1.0 - rho_prop) + log_prior_prop;
  
  double ss_curr = 0.0;
  double ss_prop = 0.0;
  double n;
  for(int subj_ix = 0; subj_ix < di.n; subj_ix ++){
    n = (double) di.ni[subj_ix];
    ss_curr = di.r2_sum[subj_ix]/(1.0 - rho_curr) - pow(di.r_sum[subj_ix], 2) * rho_curr/( (1.0 - rho_curr) * (1.0 + rho_curr * (n - 1.0)));
    ss_prop = di.r2_sum[subj_ix]/(1.0 - rho_prop) - pow(di.r_sum[subj_ix], 2) * rho_prop/( (1.0 - rho_prop) * (1.0 + rho_prop * (n - 1.0)));
    log_post_curr += -0.5 * ss_curr/pow(sigma, 2.0);
    log_post_prop += -0.5 * ss_prop/pow(sigma, 2.0);
  }
  /*
   no need to recompute each individual's sum of residuals or sum of squared residuals, as these are being tracked by di.
  std::vector<double> tmp_sum(di.N, 0.0); // hold sum of each individual's full residual
  std::vector<double> tmp_sum2(di.N, 0.0);// hold sum of each individual's squared residual
  int subj_id = 0;
  for(int i = 0; i < di.N; i++){
    subj_id = di.subj_id[i]; // subject for the i-th observation
    tmp_sum[subj_id] += di.rp[i];
    tmp_sum2[subj_id] += pow(di.rp[i], 2.0);
  }
  
  double ss_curr = 0.0;
  double ss_prior;
  for(int subj_ix = 0; subj_ix < di.n; subj_ix++){
    double n = (double) di.ni[subj_ix];
    
    ss_curr = tmp_sum2[subj_ix]/(1.0 - rho_curr) - pow(tmp_sum[subj_ix],2) * rho_curr/( (1.0 - rho_curr) * (1 + rho_curr * (n-1.0)));
    ss_prop = tmp_sum2[subj_ix]/(1.0 - rho_prop) - pow(tmp_sum[subj_ix],2) * rho_prop/( (1.0 - rho_prop) * (1 + rho_prop * (n-1.0)));
    
    log_post_curr += -0.5 * ss/pow(sigma, 2);
    log_post_prop += -0.5 * ss_prop.pow(sigma,2);
  }
  */
  double log_alpha = log_post_prop - log_post_curr;
  if(gen.log_uniform() <= log_alpha) rho = rho_prop;
  else rho = rho_curr;

}

void update_sigma_ind(double &sigma, double &nu, double &lambda, data_info &di, RNG &gen)
{
  double scale_post = lambda * nu;
  double nu_post = nu + di.N; // total number of observations
  for(int subj_ix = 0; subj_ix < di.n; subj_ix++) scale_post += di.r2_sum[subj_ix]; // add sum of squared residual for each subject
  sigma = sqrt( scale_post/gen.chi_square(nu_post));
}
void update_sigma_cs(double &sigma, double &rho, double &nu, double &lambda, data_info &di, RNG &gen)
{
  double scale_post = lambda * nu;
  double nu_post = nu + di.N;
  double n;
  for(int subj_ix = 0; subj_ix < di.n; subj_ix++){
    n = (double) di.ni[subj_ix];
    if(di.r2_sum[subj_ix] <= rho * pow(di.r_sum[subj_ix],2)/(1 + rho * (n - 1.0))){
      Rcpp::Rcout << "[update_sigma_cs]: something is fishy" << std::endl;
      Rcpp::Rcout << "  subj_ix = " << subj_ix << " r_sum = " << di.r_sum[subj_ix] << " r2_sum = " << di.r2_sum[subj_ix] << std::endl;
      Rcpp::Rcout << "  first term = " << di.r2_sum[subj_ix] << std::endl;
      Rcpp::Rcout << "  second term = " << rho * pow(di.r_sum[subj_ix],2)/(1 + rho * (n - 1.0)) << std::endl;
    }
    scale_post += di.r2_sum[subj_ix]/(1.0 - rho) - rho * pow(di.r_sum[subj_ix],2)/( (1.0 - rho) * (1.0 + rho * (n - 1.0)));
  }
  //Rcpp::Rcout << "[update_sigma_cs]: scale_post = " << scale_post << " nu_post = " << nu_post << std::endl;
  sigma = sqrt( scale_post/gen.chi_square(nu_post));
}
