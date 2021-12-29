#include "funs.h"

void tree_traversal(suff_stat &ss, tree &t, data_info &di){
  double* zz;
  tree::tree_cp bn;
  int nid;
  ss.clear(); // clear out the sufficient stat
  tree::npv; // vector of all bottom node vectors
  t.get_bots(npv);
  suff_stat_it ss_it; // used to look up which element of the map we update
  
  for(npv_it it = bnv.begin(); it != bnv.end(); ++it){
    nid = it->get_nid(); // get id of the bottom node
    // add element to map with key = id of the bottom node referenced by it
    // value is a vector of length di.n, each element is a vector of length 0 by default
    // save to start pushing back now
    ss.insert(std::pair<int, vv_int>(nid, vv_int(di.n))); //
  }
  
  // ready to loop over all observations
  for(size_t i = 0; i < di.n; i++){
    for(size_t j = 0; j < di.n_vec[i]; j++){
      zz = di.z + (j + di.start_index[i])*di.R; // zz now points to 1st modifier of j-th observation of subject i
      bn = t.get_bn(zz);
      if(bn == 0) Rcpp::stop("[tree_traversal]: bottom node not found!");
      else{
        nid = bn->get_nid(); // now we have the id of the bottom node containing j-th observation of subject i
        if(ss.count(nid) != 1) Rcpp::stop("[tree_traversal]: bottom node id is not a key in suff_stat map");
        else{
          ss_it = ss.find(nid); // iterator now set at the correct element of our suff_stat map
          // ss_it->second is a vector of vector of integers, telling us which observations are in this bottom node
          // ss_it->second[i] is a vector of integers, meant to hold the indices of subject i's observations landing in this bottom node
          ss_it->second[i].push_back(j); // ss_it->second[i] is vector
        }
      } // closes if/else checking that j-th obs of subject i has a valid bottom node
    } // closes loop over subject i's observations
  } // closes loop over subjects

}



/*
void fit(double* ftemp, tree &t, data_info &di)
{
  double* xx_cont;
  int* xx_cat;
  tree::tree_cp bn;
  for(size_t i = 0; i < di.n; i++){
    xx_cont = di.x_cont + i * di.p_cont;
    xx_cat = di.x_cat + i * di.p_cat;
    bn = t.get_bn(xx_cont, xx_cat);
    
    if(bn == 0) Rcpp::stop("[fit]: bottom node not found!");
    else ftemp[i] = bn->get_mu();
  }
}
*/
/*
void fit_vec(double* ftemp, std::vector<tree> &t_vec, data_info &di){
  double* xx_cont;
  int* xx_cat;
  tree::tree_cp bn;
  for(size_t i = 0; i < di.n; i++) ftemp[i] = 0.0; // reset value of ftemp
  for(size_t m = 0; m < t_vec.size(); m++){
    for(size_t i = 0; i < di.n; i++){
      xx_cont = di.x_cont + i * di.p_cont;
      xx_cat = di.x_cat + i * di.p_cat;
      bn = t_vec[m].get_bn(xx_cont, xx_cat);
      if(bn == 0) Rcpp::stop("[fit_vec]: bottom node not found!");
      else ftemp[i] += bn->get_mu();
    }
  }
}
*/
/*
void fit(arma::vec &ftemp, tree &t, data_info &di)
{
  double* xx_cont;
  int* xx_cat;
  tree::tree_cp bn;
  for(size_t i = 0; i < di.n; i++){
    xx_cont = di.x_cont + i * di.p_cont;
    xx_cat = di.x_cat + i * di.p_cat;
    bn = t.get_bn(xx_cont, xx_cat);
    
    if(bn == 0) Rcpp::stop("[fit]: bottom node not found!");
    else ftemp(i) = bn->get_mu();
  }
}
void fit_vec(arma::vec &ftemp, std::vector<tree> &t_vec, data_info &di)
{
  double* xx_cont;
  int* xx_cat;
  tree::tree_cp bn;
  ftemp.zeros(); // reset every entry in ftemp to 0
  
  for(size_t m = 0; m < t_vec.size(); m++){
    for(size_t i = 0; i < di.n; i++){
      xx_cont = di.x_cont + i * di.p_cont;
      xx_cat = di.x_cat + i * di.p_cat;
      bn = t_vec[m].get_bn(xx_cont, xx_cat);
      if(bn == 0) Rcpp::stop("[fit_vec]: bottom node not found!");
      else ftemp[i] += bn->get_mu();
    }
  }
}

*/


void compute_suff_stat_grow(suff_stat &orig_suff_stat, suff_stat &new_suff_stat, int &nx_nid, int &v, double &c, tree &t, data_info &di)
{
  double* zz;

  int nxl_nid = 2*nx_nid; // id of proposed left child of nx
  int nxr_nid = 2*nx_nid+1; // id of propoes right child of nx
  int j;
  suff_stat_it nx_it = orig_suff_stat.find(nx_nid); // iterator corresponding to nx in suff_stat map
  
  new_suff_stat.clear();
  // new_suff_stat should be a copy of orig_suff_stat
  for(suff_stat_it it = orig_suff_stat.begin(); it != orig_suff_stat.end(); ++it){
    new_suff_stat.insert(std::pair<int, vv_int>(it->first, it->second));
  }
  new_suff_stat.insert(std::pair<int, vv_int>(nxl_nid, vv_int(di.n))); // new_suff_stat now has an element for the proposed left child of nx
  new_suff_stat.insert(std::pair<int, vv_int>(nxr_nid, vv_int(di.n))); // new_suff_stat now has an element for the proposed right child of nx
  new_suff_stat.erase(nx_id); // new_suff_stat no longer has an element for nx!
  
  suff_stat_it nxl_it = new_suff_stat.find(nxl_nid); // iterator telling us to where in new_suff_stat node nxl corresponds
  suff_stat_it nxr_it = new_suff_stat.find(nxr_nid); // iterator telling us to where in new_suff_stat node nxr corresponds

  for(int i = 0; i < di.n; i++){
    // now we loop over all the elements contained in the i-th element of orig_ss_it->second
    for(int_it j_ix = nx_it->second[i].begin(); j_ix != nx_it->second[i].end(); ++j_ix){
      j = *(j_ix); // get the actual index of the observation
      zz = di.z + (j + di.start_index[i]) * di.R; // zz now points to first modifier of j-th observation of subject i
      if(zz[v] < c){
        // j-th observation of i-th subject was at nx and now will go to nxl
        nxl_it->second[i].push_back(j);
      } else if(zz[v] >= c){
        // j-th observation of i-th subject was at nx and now will go to nxr
        nxr_it->second[i].push_back(j);
      } else Rcpp::stop("[compute_ss_grow]: observation doesn't go to left or right child!"); // should absolutely never be triggered
    } // closes loop over observations of subject i wh
  } // closes loop over all subjects
  // at this point, we have mapped all observations to their bottom nodes in the original tree and the tree proposed by our GROW move
}


void compute_ss_prune(stuff_stat &orig_suff_stat, suff_stat &new_suff_stat, int &nl_nid, int &nr_nid, int &np_nid, tree &t, data_info &di)
{
  
  suff_stat_it nl_it = orig_suff_stat.find(nl_nid); // iterator in orig_suff_stat to which nl corresponds
  suff_stat_it nr_it = orig_suff_stat.find(nr_nid); // iterator in orig_suff_stat to which nr_corresponds
  int j;
  
  new_suff_stat.clear();
  // new_suff_stat should be a copy of orig_suff_stat
  for(suff_stat_it it = orig_suff_stat.begin(); it != orig_suff_stat.end(); ++it){
    new_suff_stat.insert(std::pair<int, vv_int>(it->first, it->second));
  }
  new_suff_stat.insert(std::pair<int, vv_int>(np_nid, vv_int(di.n))); // new_suff_stat now has an element for the parent node of nl and nr
  new_suff_stat.erase(nl_nid); // delete nl from new_suff_stat
  new_suff_stat.erase(nr_nid); // delete nr from new_suff_stat
  
  suff_stat_it np_it = new_suff_stat.find(np_nid);// iterator in new_suff_stat corresponding to leaf node np in proposed tree
  for(int i = 0; i < di.n; i++){
    // add all observations from subject i that land in nl to the new suff stat element corresponding to np (parent of nl)
    for(int_it j_ix = nl_it->second[i].begin(); j_ix != nl_it->second[i].end(); ++j_ix){
      j = *(j_ix);
      np_it->second[i].push_back(j);
    }
    // add all observations from subject i that land in nr to the new suff stat element corresponding to np (parent of nr)
    for(int_it j_ix = nr_it->second[i].begin(); j_ix != nr_it->second[i].end(); ++j_ix){
      j = *(j_ix);
      np_it->second[i].push_back(j);
    }
  }
  
  
  
}




// no longer need this

void compute_suff_stat_prune(suff_stat &sl, suff_stat &sr, suff_stat &sp, tree &t,
                             tree::tree_cp nl, tree::tree_cp nr, data_info &di)
{
  double* xx_cont; // pointer to the continuous predictors of current observation
  int* xx_cat; // pointer to categorical predictors of current observation
  
  sl.n = 0; // resent sufficient statistics for Left(nx) in new tree
  sl.val = 0.0;
  
  sr.n = 0; // reset sufficient statistics for Right(nx) in new tree
  sr.val = 0.0;
  
  sp.n = 0; // reset sufficient statistics for nx in old tree
  sp.val = 0.0;
  
  tree::tree_cp bn;
  
  for(size_t i = 0; i < di.n; i++){
    xx_cont = di.x_cont + i * di.p_cont;
    xx_cat = di.x_cat + i * di.p_cat;
    bn = t.get_bn(xx_cont, xx_cat);
    
    if(bn == nl){
      // i-th obs goes to nl in old tree
      ++(sl.n);
      sl.val += di.rp[i];
      
      ++(sp.n);
      sp.val += di.rp[i];
    } else if(bn == nr){
      // i-th obs goes to nr in old tree
      ++(sr.n);
      sr.val += di.rp[i];
      
      ++(sp.n);
      sp.val += di.rp[i];
    } else if(bn == 0) Rcpp::stop("[compute_suff_stat_prune]: bottom node not found!");
  }
  
}


double compute_lil(suff_stat &ss, tree_prior_info &tree_pi)
{
  // reminder: posterior of jump mu is N(P^-1 Theta, P^-1)
  double P = 1.0/pow(tree_pi.tau, 2.0) + ss.n;
  double Theta = tree_pi.mu0/pow(tree_pi.tau, 2.0) + ss.val;
  return(-0.5 * log(P) + 0.5 * Theta * Theta / P);
}

void compute_mu_post(double& post_mean, double& post_sd, suff_stat &ss, tree_prior_info &tree_pi)
{
  double P = 1.0/pow(tree_pi.tau, 2.0) + ss.n;
  double Theta = tree_pi.mu0/pow(tree_pi.tau, 2.0) + ss.val;
  post_sd = sqrt(1.0/P);
  post_mean = Theta/P;
}

double compute_lil_ratio_grow(tree &t, tree::tree_cp nx, rule_t &rule, data_info &di, tree_prior_info &tree_pi, bool debug)
{
  suff_stat sp;
  suff_stat sl;
  suff_stat sr;
  compute_suff_stat_grow(sl, sr, sp, t, nx, rule, di);
  
  double lill = compute_lil(sl, tree_pi);
  double lilr = compute_lil(sr, tree_pi);
  double lilp = compute_lil(sp, tree_pi);
  
  if(debug){
    Rcpp::Rcout << "  sl.n = " << sl.n << "  sl.val = " << sl.val << " lill = " << lill << std::endl;
    Rcpp::Rcout << "  sr.n = " << sr.n << "  sr.val = " << sr.val << " lilr = " << lilr << std::endl;
    Rcpp::Rcout << "  sp.n = " << sp.n << "  sp.val = " << sp.val << " lilp = " << lilp << std::endl;
  }
  double lil_ratio = lill + lilr - lilp - 1.0 * log(tree_pi.tau) - 0.5 * pow(tree_pi.mu0/tree_pi.tau, 2.0);
  return(lil_ratio);
  
  
  
}

double compute_lil_ratio_prune(tree &t, tree::tree_cp nxl, tree::tree_cp nxr, data_info &di, tree_prior_info &tree_pi, bool debug)
{
  suff_stat sp;
  suff_stat sl;
  suff_stat sr;
  compute_suff_stat_prune(sl, sr, sp, t, nxl, nxr, di);
  double lill = compute_lil(sl, tree_pi);
  double lilr = compute_lil(sr, tree_pi);
  double lilp = compute_lil(sp, tree_pi);
  
  if(debug){
    Rcpp::Rcout << "  sl.n = " << sl.n << "  sl.val = " << sl.val << " lill = " << lill << std::endl;
    Rcpp::Rcout << "  sr.n = " << sr.n << "  sr.val = " << sr.val << " lilr = " << lilr << std::endl;
    Rcpp::Rcout << "  sp.n = " << sp.n << "  sp.val = " << sp.val << " lilp = " << lilp << std::endl;
  }
  
  double lil_ratio = lilp - lill - lilr + log(tree_pi.tau) + 0.5 * pow(tree_pi.mu0/tree_pi.tau, 2.0);
  return(lil_ratio);
  

}

void draw_mu(tree &t, data_info &di, tree_prior_info &tree_pi, RNG &gen)
{
  std::vector<suff_stat> sv;
  tree::npv bnv;
  
  compute_suff_stat_vec(sv, bnv, t, di); // sv and bnv are now aligned
  
  double post_mean = 0.0;
  double post_sd = 0.0;
  
  for(size_t l = 0; l < bnv.size(); l++){
    compute_mu_post(post_mean, post_sd, sv[l], tree_pi);
    bnv[l]->set_mu(gen.normal(post_mean, post_sd));
  }
}
