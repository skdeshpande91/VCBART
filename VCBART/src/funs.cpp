//
//  funs.cpp
//  
//
//  Created by Sameer Deshpande on 8/21/19.
//
#include "funs.h"

//--------------------------------------------------
// center and scale the columns of X and the outcome Y

void prepare_x(arma::mat &X_all, std::vector<double> &x_col_mean, std::vector<double> &x_col_sd)
{
  
  bool intercept = arma::all(X_all.col(0) == 1); // check if first column of X is all ones.
  
  size_t n_obs = X_all.n_rows; // total number of observations
  size_t p = X_all.n_cols;
  double tmp_sum = 0.0; // hold running sum of X
  double tmp_sum2 = 0.0; // hold running sum of X^2
  size_t tmp_count = 0; // counts number of non-missing observations for each X. This will almost always be 0
  
  // re-size containers
  x_col_mean.clear();
  x_col_mean.resize(p);
  x_col_sd.clear();
  x_col_sd.resize(p);
  
  for(size_t k = 0; k < p; k++){
    tmp_sum = 0.0;
    tmp_sum2 = 0.0;
    tmp_count = 0;
    for(size_t i = 0; i < n_obs; i++){
      if(X_all(i,k) == X_all(i,k)){
        tmp_count++;
        tmp_sum += X_all(i,k);
        tmp_sum2 += X_all(i,k) * X_all(i,k);
      } // closes loop checking that X_all(i,k) == X_all(i,k)
    } // closes loop over rows of X_all
    
    if(tmp_count < 2) Rcpp::stop("Must have at least two observations per task!");
    else{
      x_col_mean[k] = tmp_sum/tmp_count;
      x_col_sd[k] = sqrt(1.0/( (double) tmp_count - 1.0) * (tmp_sum2 - tmp_sum * tmp_sum/( (double) tmp_count)));
      if( (k != 0) || (intercept == false) ){
        for(size_t i = 0; i < n_obs; i++){
          if(X_all(i,k) == X_all(i,k)){
            X_all(i,k) -= x_col_mean[k];
            X_all(i,k) /= x_col_sd[k];
          } // inner condition verifies that there's not missing values in X
        } // closes loop that performs centering and re-scaling
      } // outer condition ensures we standardize X[,k] only if (i) k != 0 or (ii) k == 0 & intercept = false
    } // closes else checking that we have at least 2 valid observations of X_k
  } // closes loop over columns of X_all
  if(intercept == true){
    // if we have an intercept, when we do not want to mess around with the scaling
    x_col_mean[0] = 0.0;
    x_col_sd[0] = 1.0;
  }
  
  
  
}


void prepare_x(arma::mat &X_train, arma::mat &X_test, std::vector<double> &x_col_mean, std::vector<double> &x_col_sd)
{
  
  bool intercept = arma::all(X_train.col(0) == 1); // check if first column of X is all ones.
  
  size_t n_obs = X_train.n_rows; // total number of observations
  size_t p = X_train.n_cols;
  double tmp_sum = 0.0; // hold running sum of X
  double tmp_sum2 = 0.0; // hold running sum of X^2
  size_t tmp_count = 0; // counts number of non-missing observations for each X. This will almost always be 0
  
  // re-size containers
  x_col_mean.clear();
  x_col_mean.resize(p);
  x_col_sd.clear();
  x_col_sd.resize(p);
  
  for(size_t k = 0; k < p; k++){
    tmp_sum = 0.0;
    tmp_sum2 = 0.0;
    tmp_count = 0;
    for(size_t i = 0; i < n_obs; i++){
      if(X_train(i,k) == X_train(i,k)){
        tmp_count++;
        tmp_sum += X_train(i,k);
        tmp_sum2 += X_train(i,k) * X_train(i,k);
      } // closes loop checking that X_all(i,k) == X_all(i,k)
    } // closes loop over rows of X_all
    
    if(tmp_count < 2) Rcpp::stop("Must have at least two observations per task!");
    else{
      x_col_mean[k] = tmp_sum/tmp_count;
      x_col_sd[k] = sqrt(1.0/( (double) tmp_count - 1.0) * (tmp_sum2 - tmp_sum * tmp_sum/( (double) tmp_count)));
      if( (k != 0) || (intercept == false) ){
        for(size_t i = 0; i < n_obs; i++){
          if(X_train(i,k) == X_train(i,k)){
            X_train(i,k) -= x_col_mean[k];
            X_train(i,k) /= x_col_sd[k];
          } // inner condition verifies that there's not missing values in X
        } // closes loop that performs centering and re-scaling
        
        for(size_t i = 0; i < X_test.n_rows; i++){
          if(X_test(i,k) == X_test(i,k)){
            X_test(i,k) -= x_col_mean[k];
            X_test(i,k) /= x_col_sd[k];
          }
        } // closes loop that centers and re-scales X_test
        
      } // outer condition ensures we standardize X[,k] only if (i) k != 0 or (ii) k == 0 & intercept = false
    } // closes else checking that we have at least 2 valid observations of X_k
  } // closes loop over columns of X_all
  if(intercept == true){
    // if we have an intercept, when we do not want to mess around with the scaling
    x_col_mean[0] = 0.0;
    x_col_sd[0] = 1.0;
  }
  
  
  
}


void prepare_y(arma::vec &Y, double &y_mean, double &y_sd, double &y_max, double &y_min)
{
  size_t n_obs = Y.size();
  double tmp_sum = 0.0; // hold running sum of Y
  double tmp_sum2 = 0.0; // holds running sum of Y^2
  size_t tmp_count = 0; // counts number of non-missing observations
  size_t min_index = -1; // holds index of min element of Y
  size_t max_index = -1; // holds index of max element of Y
  for(size_t i = 0; i < n_obs; i++){
    if(Y(i) == Y(i)){
      tmp_count++;
      tmp_sum += Y(i);
      tmp_sum2 += Y(i) * Y(i);
      if(max_index == -1) max_index = i;
      else if(Y(i) > Y(max_index)) max_index = i;
      
      if(min_index == -1) min_index = i;
      else if(Y(i) < Y(min_index)) min_index = i;
    } // closes loop checking that Y(i) == Y(i)
  } // closes loop over all observations
  if(tmp_count < 2) Rcpp::stop("Must have at least two observations");
  else{
    y_mean = tmp_sum/tmp_count;
    y_sd = sqrt( 1.0/( (double) tmp_count -1.0) * (tmp_sum2 - tmp_sum * tmp_sum/( (double) tmp_count)));
    for(size_t i = 0; i < n_obs; i++){
      if(Y(i) == Y(i)){
        Y(i) -= y_mean;
        Y(i) /= y_sd;
      }
    }
    y_min = Y(min_index);
    y_max = Y(max_index);
  } // closes else checking that we have at least 2 observed values of Y
}

/*
//--------------------------------------------------
void prepare_precision_ar(std::vector<std::vector<arma::mat> > &Omega_all, std::vector<std::vector<double> > &log_det_all, std::vector<double> rho, data_info &di)
{
  size_t W = rho.size();
  Omega_all.clear();
  log_det_all.clear();
  
  Omega_all.resize(W);
  log_det_all.resize(W);
  arma::mat tmp_Sigma;
  arma::mat tmp_Omega;
  double tmp_log_det = 0.0;
  double tmp_log_det_sign = 0.0;

  for(size_t w = 0; w < W; w++){
    //Rcpp::Rcout << "  w = " << w + 1 << "log_det:" << std::endl; // remember R is 1-indexed!
    for(size_t i = 0; i < di.N; i++){
      tmp_Sigma.set_size(di.n[i], di.n[i]);
      tmp_Omega.set_size(di.n[i], di.n[i]);
      for(size_t j = 0; j < di.n[i]; j++){
        for(size_t jj = 0; jj < di.n[i]; jj++){
          // for AR kernel, we treat the first column of Z as the "time" index.
          tmp_Sigma(j,jj) = ar_kernel(di.z[di.R * (di.start_index[i] + j)], di.z[di.R * (di.start_index[i] + jj)], rho[w]);
          tmp_Sigma(jj,j) = tmp_Sigma(j,jj);
        } // closes inner loop over observations jj
      } // closes outer loop over observations j
      tmp_Omega = arma::inv_sympd(tmp_Sigma);
      arma::log_det(tmp_log_det, tmp_log_det_sign, tmp_Omega);
      Omega_all[w].push_back(tmp_Omega);
      log_det_all[w].push_back(tmp_log_det);
      //Rcpp::Rcout << " " << tmp_log_det ;
    } // closes loop over individual's i
    //Rcpp::Rcout << std::endl;
  } // closes loop over w
}



void prepare_precision_cs(std::vector<std::vector<arma::mat> > &Omega_all, std::vector<std::vector<double> > &log_det_all, std::vector<double> rho, data_info &di)
{
  size_t W = rho.size();
  Omega_all.clear();
  log_det_all.clear();
  
  Omega_all.resize(W);
  log_det_all.resize(W);
  arma::mat tmp_Sigma;
  arma::mat tmp_Omega;
  double tmp_log_det = 0.0;
  double tmp_log_det_sign = 0.0;
  
  for(size_t w = 0; w < W; w++){
    for(size_t i = 0; i < di.N; i++){
      tmp_Sigma.set_size(di.n[i], di.n[i]);
      tmp_Omega.set_size(di.n[i], di.n[i]);
      for(size_t j = 0; j < di.n[i]; j++){
        for(size_t jj = 0; jj < di.n[i]; jj++){
          // for AR & CS kernel, only the first modifier Z_1 matters. This is indexed by di.R*(di.start_index[i] + j)
          tmp_Sigma(j,jj) = cs_kernel(di.z[di.R * (di.start_index[i] + j)], di.z[di.R * (di.start_index[i] + jj)], rho[w]);
          tmp_Sigma(jj,j) = tmp_Sigma(j,jj);
        } // closes inner loop over observations jj
      } // closes outer loop over observations j
      tmp_Omega = arma::inv_sympd(tmp_Sigma);
      arma::log_det(tmp_log_det, tmp_log_det_sign, tmp_Omega);
      Omega_all[w].push_back(tmp_Omega);
      log_det_all[w].push_back(tmp_log_det);
    } // closes loop over individual's i
  } // closes loop over w}

}

void prepare_precision_ind(std::vector<std::vector<arma::mat> > &Omega_all, std::vector<std::vector<double> > &log_det_all, data_info &di)
{
  size_t W = 1;
  Omega_all.clear();
  log_det_all.clear();
  
  Omega_all.resize(W);
  log_det_all.resize(W);
  arma::mat tmp_Sigma;
  arma::mat tmp_Omega;
  
  double tmp_log_det = 0.0;
  double tmp_log_det_sign = 0.0;
  
  for(size_t w = 0; w < W; w++){
    for(size_t i = 0; i < di.N; i++){
      tmp_Sigma.set_size(di.n[i], di.n[i]);
      tmp_Omega.set_size(di.n[i], di.n[i]);
      tmp_Omega.eye(); // Omega is the identity matrix
      Omega_all[w].push_back(tmp_Omega);
      log_det_all[w].push_back(0.0);
    }
  }
}
*/
//--------------------------------------------------
//does a (bottom) node have variables you can split on?
bool cansplit(tree::tree_p n, xinfo& xi)
{
  int L,U;
  bool v_found = false; //have you found a variable you can split on
  size_t v=0;
  while(!v_found && (v < xi.size())) { //invar: splitvar not found, vars left
    L=0; U = xi[v].size()-1;
    n->rg(v,&L,&U);
    if(U>=L) v_found=true;
    v++;
  }
  return v_found;
}
//--------------------------------------------------
//compute prob of a birth, goodbots will contain all the good bottom nodes
double getpb(tree &t, xinfo &xi, tree_prior_info &tree_pi, tree::npv &goodbots){
  double pb; // prob of birth to be returned
  tree::npv bnv; // all the bottom nodes
  t.getbots(bnv); // actually find all of the bottom nodes
  for(size_t i = 0; i != bnv.size(); i++){
    if(cansplit(bnv[i], xi)) goodbots.push_back(bnv[i]);
  }
  if(goodbots.size() == 0) pb = 0.0; // there are no bottom nodes you can split on
  else{
    if(t.treesize() == 1) pb = 1.0; // tree only has one node
    else pb = tree_pi.pb;
  }
  return pb;
}
//--------------------------------------------------
//find variables n can split on, put their indices in goodvars
void getgoodvars(tree::tree_p n, xinfo& xi,  std::vector<size_t>& goodvars)
{
  int L,U;
  for(size_t v=0;v!=xi.size();v++) {//try each variable
    L=0; U = xi[v].size()-1;
    n->rg(v,&L,&U);
    if(U>=L) goodvars.push_back(v);
  }
}
//--------------------------------------------------
//get prob a node grows, 0 if no good vars, else alpha/(1+d)^beta
/*
double pgrow(tree::tree_p n, xinfo &xi, tree_prior_info &tree_pi)
{
  if(cansplit(n,xi)) return tree_pi.alpha/pow(1.0 + n->depth(), tree_pi.beta);
  else return 0.0;
}
*/
double pgrow(tree::tree_p n, xinfo &xi, tree_prior_info &tree_pi, size_t k)
{
  if(cansplit(n,xi)) return tree_pi.alpha[k]/pow(1.0 + n->depth(), tree_pi.beta[k]);
  else return 0.0;
}
//--------------------------------------------------
//get sufficients stats for all bottom nodes
void allsuff(tree &x, xinfo &xi, data_info &di, tree::npv &bnv, std::vector<sinfo> &sv)
{
  // Bottom nodes are written to bnv.
  // Suff stats for each bottom node are written to elements (each of class sinfo) of sv.
  // Initialize data structures
  tree::tree_cp tbn; //the pointer to the bottom node for the current observations.  tree_cp bc not modifying tree directly.
  size_t ni; //the  index into vector of the current bottom node
  double *zz; //current z
  bnv.clear(); // Clear the bnv variable if any value is already saved there.
  x.getbots(bnv); // Save bottom nodes for x to bnv variable.
  
  
  typedef tree::npv::size_type bvsz;  // Is a better C way to set type.  (tree::npv::size_type) will resolve to an integer,
  // or long int, etc.  We don't have to know that ahead of time by using this notation.
  bvsz nb = bnv.size();   // Initialize new var nb of type bvsz for number of bottom nodes, then...
  
  sv.clear();
  sv.resize(nb);
  
  // need to re-size the members of sv
  // use l to index the bottom nodes
  for(size_t l = 0; l != bnv.size(); l++){
    sv[l].n = 0; // total number of observations in this node
    sv[l].node_count.clear();
    sv[l].node_count.resize(di.N, 0);
    sv[l].I.clear();
    sv[l].I.resize(di.N, std::vector<size_t>(1)); // this makes sv.I have length N and each element have length 1
    for(size_t i = 0; i < di.N; i++) sv[l].I[i].clear(); // this clears out each element of N to now have length 0. We can safely push_back now.
  }

  // bnmap is a tuple (lookups, like in Python).  Want to index by bottom nodes.
  // [SKD] 27 Aug 2019 -- I changed the loop to iterate over l for consistenct with what's above
  std::map<tree::tree_cp,size_t> bnmap;
  for(bvsz l=0;l!=bnv.size();l++) bnmap[bnv[l]]=l;  // bnv[l]
  //map looks like
  // bottom node 1 ------ 1
  // bottom node 2 ------ 2
  // Loop through each observation.  Push each obs x down the tree and find its bottom node,
  // then index into the suff stat for the bottom node corresponding to that obs.
  
  for(size_t i = 0; i < di.N; i++){
    for(size_t j = 0; j < di.n[i]; j++){
      zz = di.z + (di.start_index[i] + j)*di.R; // Now points to beginning of j^th set of modifiers for individual i.
      tbn = x.bn(zz,xi); // Finds bottom node for observation j, individual i
      ni = bnmap[tbn]; // Maps bottom node to integer index
      ++(sv[ni].n); // increment total count (across individuals)
      ++(sv[ni].node_count[i]); // increment the count at the node for individual i
      sv[ni].I[i].push_back(j); // remember I only holds the index for individual i.
    }
  }
}
//get sufficient stats for children (v,c) of node nx in tree x
//[SKD]: used in the birth proposals.
// in node nx we split on variable v at cutpoint c
// sl will contain the observations sent to the left child, sr contains the observations sent to the right child
void getsuff(tree &x, tree::tree_cp nx, size_t v, size_t c, xinfo &xi, data_info &di, sinfo &sl, sinfo &sr, sinfo &st)
{
  double* zz;
  sl.n = 0;
  sl.node_count.resize(di.N, 0);
  sl.I.resize(di.N, std::vector<size_t>(1));
  
  sr.n = 0;
  sr.node_count.resize(di.N, 0);
  sr.I.resize(di.N, std::vector<size_t>(1));
  
  st.n = 0;
  st.node_count.resize(di.N, 0);
  st.I.resize(di.N, std::vector<size_t>(1));
  
  for(size_t i = 0; i < di.N; i++){
    sl.I[i].clear();
    sr.I[i].clear();
    st.I[i].clear();
  }
  
  for(size_t i =0; i < di.N; i++){
    for(size_t j = 0; j < di.n[i]; j++){
      zz = di.z + (di.start_index[i] + j)*di.R; // zz now points to start of covariates for observation j, individual i
      if(nx == x.bn(zz,xi)){ // in original tree x, does this observation land in the node being split
        if(zz[v] < xi[v][c]){ // it goes to the left child!
          ++(sl.n); // increment total count
          ++(sl.node_count[i]); // increment the count for this node
          sl.I[i].push_back(j); // add the index j to the appropriate leaf index vector
        } else{
          ++(sr.n); // increment total count
          ++(sr.node_count[i]); // increment the count for this node
          sr.I[i].push_back(j); // add the index j to the appropriate leaf index vect
        }
        ++(st.n); // increment total count
        ++(st.node_count[i]);
        st.I[i].push_back(j);
        
        
      } // closes if checking whether observation lands in node nx, which is being split
    } // closes loop over observations for individual i
  } // closes loop over individuals i
}

//get sufficient stats for children (v,c) of node nx in tree x
//[SKD]: used in the death
// we are trying to collapse nodes nl and nr into their parent
// this function will figure out which observations go to nl and nr
// sl will contain the observations sent to the left child, sr contains the observations sent to the right child
// st is basically the set of sufficient statistics for the proposed combined node
void getsuff(tree &x, tree::tree_cp nl, tree::tree_cp nr, xinfo &xi, data_info &di, sinfo &sl, sinfo &sr, sinfo &st)
{
  double* zz;
  sl.n = 0;
  sl.node_count.resize(di.N, 0);
  sl.I.resize(di.N, std::vector<size_t>(1));
  
  sr.n = 0;
  sr.node_count.resize(di.N, 0);
  sr.I.resize(di.N, std::vector<size_t>(1));
  
  st.n = 0;
  st.node_count.resize(di.N,0);
  st.I.resize(di.N, std::vector<size_t>(1));
  
  for(size_t i = 0; i < di.N; i++){
    sl.I[i].clear();
    sr.I[i].clear();
    st.I[i].clear();
  }
  
  for(size_t i = 0; i < di.N; i++){
    for(size_t j = 0; j < di.n[i]; j++){
      zz = di.z + (di.start_index[i] + j)*di.R;
      tree::tree_cp bn = x.bn(zz,xi); // gets the bottom node into which observation falls
      if(bn == nl){
        ++(sl.n);
        ++(sl.node_count[i]);
        sl.I[i].push_back(j);
        
        ++(st.n);
        ++(st.node_count[i]);
        st.I[i].push_back(j);
      } else if(bn == nr){
        ++(sr.n);
        ++(sr.node_count[i]);
        sr.I[i].push_back(j);
        
        ++(st.n);
        ++(st.node_count[i]);
        st.I[i].push_back(j);
        
      }
    } // closes loop over observations for individual i
  } // closes loop over individuals
}
//--------------------------------------------------
// new getsuff functions for use in bd_fast
/*
 for growth proposal we do the following:
 1. Make a mapping from the bottom nodes of original tree x to indices of sv_x
 2. Loop over observations
     if bottom node for observation j, individual i is not nx, then the index of this bottom node in sv_y is the same as index in sv_x
     otherwise, if observation j, individual i lands in node nx that is being split we have
         if it goes to left child, index of the new bottom node in y is going to be bnmap[nx] (recycling the index)
         if it goes to the right child, index of new bottom node in y is going to bnv.size() -- remember max index of bottom node in x is bnv.size()-1
 3. Function returns two sets of sufficient statistics, one for original tree x and one for proposed tree y
 */
void getsuff(tree &x, tree::tree_cp nx, size_t v, size_t c, xinfo &xi, data_info &di, std::vector<sinfo> &sv_x, std::vector<sinfo> &sv_y)
{
  double* zz;
  
  tree::npv bnv;
  bnv.clear(); // vector of bottom nodes for original tree x
  x.getbots(bnv); // actually get the bottom nodes of x
  
  typedef tree::npv::size_type bvsz;
  bvsz nb = bnv.size(); // number of bottom nodes for x
  
  size_t ni_x; // leaf index for original tree x
  size_t ni_y; // leaf index for proposed tree y
  
  tree::tree_cp tbn; // pointer for bottom node
  
  sv_x.clear();
  sv_y.clear();
  
  sv_x.resize(nb);
  sv_y.resize(nb+1); // proposal tree will have 1 more bottom node
  
  for(bvsz l = 0; l != nb; l++){
    sv_x[l].n = 0;
    sv_x[l].node_count.resize(di.N, 0);
    sv_x[l].I.resize(di.N, std::vector<size_t>(1));
    for(size_t i = 0; i < di.N; i++) sv_x[l].I[i].clear();
  }
  for(bvsz l = 0; l != nb + 1; l++){
    sv_y[l].n = 0;
    sv_y[l].node_count.resize(di.N, 0);
    sv_y[l].I.resize(di.N, std::vector<size_t>(1));
    for(size_t i = 0; i < di.N; i++) sv_y[l].I[i].clear();
  }
  
  // make a map of bottom nodes for the original tree x
  std::map<tree::tree_cp, size_t> bnmap;
  for(bvsz l = 0; l != nb; l++) bnmap[bnv[l]] = l;
  // map looks like:
  // bottom node 0 --- 0
  // bottom node 1 --- 1
  
  for(size_t i = 0; i < di.N; i++){
    for(size_t j = 0; j < di.n[i]; j++){
      zz = di.z + (di.start_index[i] + j) * di.R; // now points to beginning of modifiers for observation j, individual i
      tbn = x.bn(zz, xi); // bottom node
      
      ni_x = bnmap[tbn];
      ++(sv_x[ni_x].n); // increment total count of observations (across individuals) in the leaf
      ++(sv_x[ni_x].node_count[i]); // increment count of observations for individual i in this leaf
      sv_x[ni_x].I[i].push_back(j); // remember I only holds index for individual i
      
      // now get the right index for *proposed tree* y
      if(tbn != nx) ni_y = bnmap[tbn];
      else{
        // observation j, individual i lands in node getting split
        if(zz[v] < xi[v][c]) ni_y = bnmap[tbn];
        else ni_y = bnv.size();
        // if node l (out of L) is split in original tree, for proposed tree we let bottom node index l be used for left child
        // and L be used for right child
      }
      ++(sv_y[ni_y].n); // increment total count of observations (across individuals) in the leaf
      ++(sv_y[ni_y].node_count[i]); // increment count of observations for individual i in this leaf
      sv_y[ni_y].I[i].push_back(j); // remember I only holds index for individual i
    } // closes inner loop over observations per individual
  } // closes outer loop over individuals
}

// new getsuff for death proposal in bd_fast
/*
 This function works slightly differently than the one in the growth proposal.
 1. We first make a map bnv_x from bottom nodes of original tree x to the index in sv_x, the vector of suff stats for original tree
 2. We then make a map bnv_y from bottom nodes of original tree x to the index in sv_y, the vector of suff stats for proposed tree
 3. We then loop over observations and update sv_x and sv_y accordingly
 */
void getsuff(tree &x, tree::tree_cp nl, tree::tree_cp nr, xinfo &xi, data_info &di, std::vector<sinfo> &sv_x, std::vector<sinfo> &sv_y)
{
  double* zz;
  tree::npv bnv;
  bnv.clear(); // vector of bottom nodes for original tree x
  x.getbots(bnv); // actually get the bottom nodes of x
  
  typedef tree::npv::size_type bvsz;
  bvsz nb = bnv.size(); // number of bottom nodes for x
  
  size_t ni_x; // leaf index for original tree x
  size_t ni_y; // leaf index for proposed tree y
  
  tree::tree_cp tbn; // pointer for bottom node
  
  sv_x.clear();
  sv_y.clear();
  
  sv_x.resize(nb);
  sv_y.resize(nb-1); // proposal tree will have 1 fewer bottom node
  
  for(bvsz l = 0; l != nb; l++){
    sv_x[l].n = 0;
    sv_x[l].node_count.resize(di.N, 0);
    sv_x[l].I.resize(di.N, std::vector<size_t>(1));
    for(size_t i = 0; i < di.N; i++) sv_x[l].I[i].clear();
  }
  for(bvsz l = 0; l != nb - 1; l++){
    sv_y[l].n = 0;
    sv_y[l].node_count.resize(di.N, 0);
    sv_y[l].I.resize(di.N, std::vector<size_t>(1));
    for(size_t i = 0; i < di.N; i++) sv_y[l].I[i].clear();
  }
  
  // make a map of bottom nodes for original tree x
  size_t ni_x_l; // index in original tree x of left child that is being collapsed
  size_t ni_x_r; // index in original tree x of right child that is being collapsed
  std::map<tree::tree_cp, size_t> bnmap;
  for(bvsz l = 0; l != nb; l++){
    bnmap[bnv[l]] = l;
    if(bnv[l] == nl) ni_x_l = l;
    else if(bnv[l] == nr) ni_x_r = l;
  }
  
  
  size_t min_index;
  size_t max_index;
  if(ni_x_l < ni_x_r){
    min_index = ni_x_l;
    max_index = ni_x_r;
  } else{
    min_index = ni_x_r;
    max_index = ni_x_l;
  }
  
  std::map<tree::tree_cp, size_t> bnmap_y; // make a map for the new tree y.
  // this maps bottom nodes of original tree x to leaf indices of new tree y
  for(bvsz l = 0; l != nb; l++){
    if(bnmap[bnv[l]] < max_index) bnmap_y[bnv[l]] = bnmap[bnv[l]];
    else if(bnmap[bnv[l]] == max_index) bnmap_y[bnv[l]] = min_index;
    else if(bnmap[bnv[l]] > max_index) bnmap_y[bnv[l]] = bnmap[bnv[l]] - 1;
  }
  
  
  for(size_t i = 0; i < di.N; i++){
    for(size_t j = 0; j < di.n[i]; j++){
      zz = di.z + (di.start_index[i] + j) * di.R; // now points to beginning of modifiers for observation j, individual i
      tbn = x.bn(zz, xi); // bottom node
      
      ni_x = bnmap[tbn]; // get the correct index for original tree x
      ni_y = bnmap_y[tbn]; // get the correct index for proposed tree y
      
      ++(sv_x[ni_x].n); // increment total count of observations (across individuals) in the leaf
      ++(sv_x[ni_x].node_count[i]); // increment count of observations for individual i in this leaf
      sv_x[ni_x].I[i].push_back(j); // remember I only holds index for individual i
      
      ++(sv_y[ni_y].n); // increment total count of observations (across individuals) in the leaf
      ++(sv_y[ni_y].node_count[i]); // increment count of observations for individual i in this leaf
      sv_y[ni_y].I[i].push_back(j); // remember I only holds index for individual i
      
    }
  }
}

//--------------------------------------------------
// fit
void fit(tree& t, xinfo& xi, data_info& di, double* fv)
{
  double *zz;
  tree::tree_cp bn;
  
  for(size_t i = 0; i < di.N; i++){
    for(size_t j = 0; j < di.n[i]; j++){
      zz = di.z + (di.start_index[i] + j) * di.R;
      bn = t.bn(zz, xi);
      fv[di.start_index[i] + j] = bn->getm();
    }
  }
}


/*
//--------------------------------------------------
// update_sigma -- using inverse gamma prior

// update individual sigma's
void update_sigma_ig(std::vector<double> &sigma, std::vector<arma::mat> &Omega_vec, std::vector<sigma_prior_info> &sigma_pi, data_info &di, RNG &gen)
{
  // a posteriori sigma_2^2 ~ IG( (nu_post/2, scale_post/2) )
  //generated by scale_post/(chisq(nu))
  
  double scale_post = 0.0; //
  double nu_post = 0.0; //
  double quad_form = 0.0;
  for(size_t i = 0; i < di.N; i++){
    scale_post = sigma_pi[i].lambda * sigma_pi[i].nu; // this is from the prior scaling
    nu_post = sigma_pi[i].nu + di.n[i]; // posterior degrees of freedom
    quad_form = 0.0;
    for(size_t j = 0; j < di.n[i]; j++){
      for(size_t jj = 0; jj < di.n[i]; jj++){
        scale_post += di.rf[j + di.start_index[i]] * Omega_vec[i](j,jj) * di.rf[jj + di.start_index[i]];
      }
    }
    //Rcpp::Rcout << "i = " << i << "  scale_post = " << scale_post << " nu_post = " << nu_post << endl;
    sigma[i] = sqrt( (scale_post)/gen.chi_square(nu_post));
  } // closes loop over individuals i
}

void update_sigma_ig(double &sigma, std::vector<arma::mat> &Omega_vec, sigma_prior_info &sigma_pi, data_info &di, RNG &gen)
{
  double scale_post = sigma_pi.lambda * sigma_pi.nu;
  double nu_post = sigma_pi.nu;
  for(size_t i = 0; i < di.N; i++){
    nu_post += di.n[i]; // add the number of observations from individual i to the degrees of freedom
    for(size_t j = 0; j < di.n[i]; j++){
      for(size_t jj = 0; jj < di.n[i]; jj++){
        scale_post += di.rf[j + di.start_index[i]] * Omega_vec[i](j,jj) * di.rf[jj + di.start_index[i]]; // add the quadratic form (rf_i)'Omega_i (rf_i)
      }
    }
  }
  sigma = sqrt( (scale_post)/gen.chi_square(nu_post));
}

// overloaded version for margin_re
void update_sigma_ig(double &sigma, double &rho, sigma_prior_info &sigma_pi, data_info &di, RNG &gen)
{
  double scale_post = sigma_pi.lambda * sigma_pi.nu;
  double nu_post = sigma_pi.nu;
  
  //double omega_offdiag = 0.0;
 // double omega_diag = 0.0;
  
  double tmp_sum = 0.0; // holds running sum of partial residuals
  double tmp_sum2 = 0.0; // holds running sum of squares of partial residuals
  
  for(size_t i = 0; i < di.N; i++){
    nu_post += di.n[i]; // add number of observations from individual i to the degrees of freedom
    // new way to compute the relevant quadratic forms is linear O(n_i).
    
    tmp_sum = 0.0;
    tmp_sum2 = 0.0;
    for(size_t j = 0; j < di.n[i]; j++){
      tmp_sum += di.rf[j + di.start_index[i]];
      tmp_sum2 += di.rf[j + di.start_index[i]] * di.rf[j + di.start_index[i]];
      if(tmp_sum != tmp_sum){
        Rcpp::Rcout << "[update_sigma]: i = " << i << " j = " << j << "  di.rf[j + di.start_index[i]] = " << di.rf[j + di.start_index[i]] << std::endl;
        Rcpp::stop("quitting");
      }
      if(tmp_sum2 != tmp_sum2){
        Rcpp::Rcout << "[update_sigma]: i = " << i << " j = " << j << "  di.rf[j + di.start_index[i]] = " << di.rf[j + di.start_index[i]] << std::endl;
        Rcpp::stop("quitting because of tmp_sum2");
      }
    }
    scale_post += tmp_sum2 / (1.0 - rho) - tmp_sum * tmp_sum * rho/( (1.0 - rho) * (1.0 + ((double) di.n[i] - 1) * rho));
    if(scale_post != scale_post){
      Rcpp::Rcout <<"[update_sigma]: i = " << i << "  tmp_sum = " << tmp_sum << "  tmp_sum2 = " << tmp_sum2 << std::endl;
      Rcpp::stop("[update_sigma]: scale_post is nan. quitting");
    }
  }
  sigma = sqrt( (scale_post)/gen.chi_square(nu_post));
  //Rcpp::Rcout << "[update_sigma]: scale_post = " << scale_post << "  nu_post = " << nu_post << "  sigma = " << sigma << std::endl;
}


void update_sigma_ig(double &sigma, sigma_prior_info &sigma_pi, data_info &di, RNG &gen)
{
  double scale_post = sigma_pi.lambda * sigma_pi.nu;
  double nu_post = sigma_pi.nu;
  
  for(size_t i = 0; i < di.N; i++){
    nu_post += di.n[i];
    for(size_t j = 0; j < di.n[i]; j++){
      scale_post += di.rf[j + di.start_index[i]] * di.rf[j + di.start_index[i]];
    }
  }
  sigma = sqrt( (scale_post)/gen.chi_square(nu_post));
}

//--------------------------------------------------
// update_sigma -- using half-Cauchy prior
void update_sigma_hc(double&sigma, std::vector<arma::mat> &Omega_vec, sigma_prior_info &sigma_pi, data_info &di, RNG &gen){
  double scale_post = 0.0; // this will be used to draw the proposal
  double nu_post = 0.0; // this will be used to draw the proposal
  for(size_t i = 0; i < di.N; i++){
    nu_post += di.n[i]; // add number of observations from individual i to the degrees of freedome
    for(size_t j = 0; j < di.n[i]; j++){
      for(size_t jj = 0; jj < di.n[i]; jj++){
        scale_post += di.rf[j + di.start_index[i]] * Omega_vec[i](j,jj) * di.rf[jj + di.start_index[i]]; // add the quadratic form (rf_i)'Omega_i(rf_i)
      }
    }
  }
  double sigma_prop = sqrt((scale_post)/gen.chi_square(nu_post)); // this is the proposed sigma'
  
  double alpha = std::min(1.0, pow(sigma_prop, 3)/pow(sigma,3) * (pow(sigma_pi.A, 2) + pow(sigma, 2))/(pow(sigma_pi.A, 2) + pow(sigma_prop,2))); // MH ratio
  if(gen.uniform() <= alpha) sigma = sigma_prop;
}

void update_sigma_hc(double &sigma, double &rho, sigma_prior_info &sigma_pi, data_info &di, RNG &gen)
{
  double scale_post = 0.0;
  double nu_post = 0.0;

  double tmp_sum = 0.0; // holds running sum of partial residuals
  double tmp_sum2 = 0.0; // holds running sum of squares of partial residuals
  
  for(size_t i = 0; i < di.N; i++){
    nu_post += di.n[i]; // add number of observations from individual i to the degrees of freedom
    
    tmp_sum = 0.0;
    tmp_sum2 = 0.0;
    for(size_t j = 0; j < di.n[i]; j++){
      tmp_sum += di.rf[j + di.start_index[i]];
      tmp_sum2 += di.rf[j + di.start_index[i]] * di.rf[j + di.start_index[i]];
    }
    scale_post += tmp_sum2 / (1.0 - rho) - tmp_sum * tmp_sum * rho/( (1.0 - rho) * (1.0 + ((double) di.n[i] - 1) * rho));
    if(scale_post != scale_post){
      Rcpp::Rcout <<"[update_sigma]: i = " << i << "  tmp_sum = " << tmp_sum << "  tmp_sum2 = " << tmp_sum2 << std::endl;
      Rcpp::stop("[update_sigma]: scale_post is nan. quitting");
    }
  }
  
  double sigma_prop = sqrt( (scale_post)/gen.chi_square(nu_post)); // this is the proposed sigma'
  double alpha = std::min(1.0, pow(sigma_prop, 3)/pow(sigma,3) * (pow(sigma_pi.A, 2) + pow(sigma, 2))/(pow(sigma_pi.A, 2) + pow(sigma_prop,2))); // MH ratio
  if(gen.uniform() <= alpha) sigma = sigma_prop;
}

void update_sigma_ht(double &sigma, double &rho, sigma_prior_info &sigma_pi, data_info &di, RNG &gen)
{
  double scale_post = 0.0;
  double nu_post = 0.0;
  
  double tmp_sum = 0.0; // holds running sum of partial residuals
  double tmp_sum2 = 0.0; // holds running sum of squares of partial residuals
  
  for(size_t i = 0; i < di.N; i++){
    nu_post += di.n[i]; // add number of observations from individual i to the degrees of freedom
    
    tmp_sum = 0.0;
    tmp_sum2 = 0.0;
    for(size_t j = 0; j < di.n[i]; j++){
      tmp_sum += di.rf[j + di.start_index[i]];
      tmp_sum2 += di.rf[j + di.start_index[i]] * di.rf[j + di.start_index[i]];
    }
    scale_post += tmp_sum2 / (1.0 - rho) - tmp_sum * tmp_sum * rho/( (1.0 - rho) * (1.0 + ((double) di.n[i] - 1) * rho));
    if(scale_post != scale_post){
      Rcpp::Rcout <<"[update_sigma]: i = " << i << "  tmp_sum = " << tmp_sum << "  tmp_sum2 = " << tmp_sum2 << std::endl;
      Rcpp::stop("[update_sigma]: scale_post is nan. quitting");
    }
  }
  
  double sigma_prop = sqrt( (scale_post)/gen.chi_square(nu_post)); // this is the proposed sigma'
  
  double log_prior_prop = -0.5 * (1.0 + sigma_pi.nu) * log(1.0 + 1.0/sigma_pi.nu * sigma_prop * sigma_prop/(sigma_pi.A * sigma_pi.A));
  double log_prior = -0.5 * (1.0 + sigma_pi.nu) * log(1.0 + 1.0/sigma_pi.nu * sigma * sigma/(sigma_pi.A * sigma_pi.A));
  double alpha = exp(3.0 * log(sigma_prop) + log_prior_prop - 3.0 * log(sigma) - log_prior);
  if(gen.uniform() <= alpha) sigma = sigma_prop;
}



void update_sigma_hc(std::vector<double> &sigma, std::vector<arma::mat> &Omega_vec, std::vector<sigma_prior_info> &sigma_pi, data_info &di, RNG &gen)
{
  double scale_post = 0.0;
  double nu_post = 0.0;
  double alpha = 0.0;
  double sigma_prop = 0.0;
  for(size_t i = 0; i < di.N; i++){
    nu_post = di.n[i];
    scale_post = 0.0;
    for(size_t j = 0; j < di.n[i]; j++){
      for(size_t jj = 0; j < di.n[i];jj++){
        scale_post += di.rf[j + di.start_index[i]] * Omega_vec[i](j,jj) * di.rf[jj + di.start_index[i]]; // add the quadratic form (rf_i)' Omega_i (rf_i)
      }
    }
    sigma_prop = sqrt( (scale_post)/gen.chi_square(nu_post)); // this is the proposed sigma'
    alpha = std::min(1.0, pow(sigma_prop,3)/pow(sigma[i], 3) * (pow(sigma_pi[i].A, 2) + pow(sigma[i],2))/(pow(sigma_pi[i].A, 2) + pow(sigma_prop,2)));
    if(gen.uniform() <= alpha) sigma[i] = sigma_prop;
  }
}




// update sigma (with HC prior) when there is no residual correlation
void update_sigma_hc(double &sigma, sigma_prior_info &sigma_pi, data_info &di, RNG &gen){
  double scale_post = 0.0 ; // used to draw the proposal
  double nu_post = 0.0; // used to draw the proposal
  for(size_t i = 0; i < di.N; i++){
    nu_post += di.n[i];
    for(size_t j = 0; j < di.n[i]; j++){
      scale_post += di.rf[j + di.start_index[i]] * di.rf[j + di.start_index[i]];
    }
  }
  double sigma_prop = sqrt( (scale_post)/gen.chi_square(nu_post)); // this is proposed sigma'
  
  double alpha = std::min(1.0, pow(sigma_prop,3)/pow(sigma,3) * (pow(sigma_pi.A,2) + pow(sigma,2))/(pow(sigma_pi.A,2) + pow(sigma_prop,2)));
  if(gen.uniform() <= alpha) sigma = sigma_prop;
}
*/

/*
//--------------------------------------------------
// update sigma_mu when we give it a half-Cauchy prior
void update_sigma_mu(double* sigma_mu, std::vector<std::vector<tree > > t_vec, tree_prior_info &tree_pi, RNG &gen)
{
  tree::npv bnv; // holds point to bottom nodes of the tree
  size_t p = t_vec.size();
  size_t L_tot = 0;
  double M_tot = 0.0;
  double alpha = 0.0;
  double sigma_prop = 0.0;
  for(size_t k = 0; k < p; k++){
    L_tot = 0;
    M_tot = 0.0;
    
    for(size_t m = 0; m < t_vec[k].size(); m++){
      bnv.clear();
      t_vec[k][m].getbots(bnv); // get the bottom node for tree m for beta_k
      for(size_t l = 0; l < bnv.size(); l++){
        L_tot++;
        M_tot += pow(bnv[l]->getm(),2);
      }
      sigma_prop = sqrt( (M_tot)/gen.chi_square(L_tot));
      alpha = std::min(1.0, pow(sigma_prop,3)/pow(sigma_mu[k],3) * (pow(tree_pi.A[k], 2) + pow(sigma_mu[k],2))/(pow(tree_pi.A[k],2) + pow(sigma_prop, 2)));
      if(gen.uniform() <= alpha) sigma_mu[k] = sigma_prop;
    } // closes loop over the trees in ensemble for beta_k
  } // closes loop over covariates.
}
*/
/*
//--------------------------------------------------
// update rho
// we update only rho_ix.

// the rho update for a global sigma
void update_rho(size_t &rho_ix, std::vector<std::vector<arma::mat > > &Omega_all, std::vector<std::vector<double> > &log_det_all, double &sigma, data_info &di, RNG &gen)
{
  size_t W = Omega_all.size(); // number of possible values
  if(W == 1){
    rho_ix = 0; // easy enough. if there's only one value of rho to consider...
  } else{
    std::vector<double> log_post_prob(W, 0.0); // holds the posterior probabilities
    double max_log_post = 0.0;
    for(size_t w = 0; w < W; w++){
      log_post_prob[w] += di.rho_log_prior[w]; // add in the log-prior probability for this value of rho
      for(size_t i = 0; i < di.N; i++){
        log_post_prob[w] += 0.5 * log_det_all[w][i];
        for(size_t j = 0; j < di.n[i]; j++){
          for(size_t jj = 0; jj < di.n[i]; jj++){
            log_post_prob[w] += (-0.5) * 1.0/(sigma * sigma) * di.rf[j + di.start_index[i]] * Omega_all[w][i](j,jj) * di.rf[jj + di.start_index[i]];
          } // closes inner loop for computing quadratic form (rf_i)' Omega_i (rf_i)
        } // closes outer loop for computing quadratic form (rf_i)' Omega_i (rf_i)
      } // closes loop over the individuals
      if( (w == 0) || (log_post_prob[w] > max_log_post)) max_log_post = log_post_prob[w];
    } // closes loop over w
    //Rcpp::Rcout << "max_log_post = " << max_log_post << std::endl;
    // subtract off the maximum log_post
    std::vector<double> post_prob(W);
    std::vector<double> cdf(W); // compute the posterior cdf
    double tmp_sum = 0.0;
    for(size_t w = 0; w < W; w++){
      log_post_prob[w] -= max_log_post;
      post_prob[w] = exp(log_post_prob[w]);
      tmp_sum += post_prob[w];
      cdf[w] = tmp_sum; // keeps a running sum
    }
    for(size_t w = 0; w < W; w++){
      post_prob[w] /= tmp_sum; // normalize the posterior probabilies
      cdf[w] /= tmp_sum; // normalizes the CDF
    }

    double unif = gen.uniform(); // draw a random uniform
    //Rcpp::Rcout << "unif = " << unif << std::endl;
    std::vector<size_t> possible_ix; // all of the indexes for which the cdf is less than unif
    
    // 17 September 2019 -- generate the uniform OUTSIDE of the loop and then check it
    
    for(size_t w = 0; w < W; w++){
      if(cdf[w] <= unif) possible_ix.push_back(w);
    }
    if(possible_ix.size() == 0) rho_ix = 0; // all non-zero values of cdf are greater than unif. goes to first value
    else if(possible_ix.size() == W) rho_ix = W-1; // unif = 1, which is never going to happen...
    else{
      rho_ix = possible_ix[possible_ix.size()-1] + 1;// cdf[rho_ix] < unif < cdf[rho_ix+1]
    }
    //Rcpp::Rcout << "rho_ix = " << rho_ix << std::endl;
  } // closes if/else checking that W = 1

}

// the rho update when we have individual sigma's
void update_rho(size_t &rho_ix, std::vector<std::vector<arma::mat> > &Omega_all, std::vector<std::vector<double> > &log_det_all, std::vector<double> &sigma, data_info &di, RNG &gen)
{
  size_t W = Omega_all.size(); // number of possible values
  if(W == 1){
    rho_ix = 0; // easy enough. if there's only one value of rho to consider...
  } else{
    std::vector<double> log_post_prob(W, 0.0); // holds the posterior probabilities
    double max_log_post = 0.0;
    for(size_t w = 0; w < W; w++){
      log_post_prob[w] += di.rho_log_prior[w]; // add in the log-prior probability for this value of rho
      for(size_t i = 0; i < di.N; i++){
        log_post_prob[w] += 0.5 * log_det_all[w][i];
        for(size_t j = 0; j < di.n[i]; j++){
          for(size_t jj = 0; jj < di.n[i]; jj++){
            log_post_prob[w] += (-0.5) * 1.0/(sigma[i] * sigma[i]) * di.rf[j + di.start_index[i]] * Omega_all[w][i](j,jj) * di.rf[jj + di.start_index[i]];
          } // closes inner loop for computing quadratic form (rf_i)' Omega_i (rf_i)
        } // closes outer loop for computing quadratic form (rf_i)' Omega_i (rf_i)
      } // closes loop over the individuals
      if( (w == 0) || (log_post_prob[w] > max_log_post)) max_log_post = log_post_prob[w];
    } // closes loop over w
    //Rcpp::Rcout << "max_log_post = " << max_log_post << std::endl;
    // subtract off the maximum log_post
    std::vector<double> post_prob(W);
    std::vector<double> cdf(W); // compute the posterior cdf
    double tmp_sum = 0.0;
    for(size_t w = 0; w < W; w++){
      log_post_prob[w] -= max_log_post;
      post_prob[w] = exp(log_post_prob[w]);
      tmp_sum += post_prob[w];
      cdf[w] = tmp_sum; // keeps a running sum
    }
    for(size_t w = 0; w < W; w++){
      post_prob[w] /= tmp_sum; // normalize the posterior probabilies
      cdf[w] /= tmp_sum; // normalizes the CDF
    }
    // printing for testing
    //Rcpp::Rcout << "log_post_prob: ";
    //for(size_t w = 0; w < W; w++) Rcpp::Rcout << " " << log_post_prob[w];
    //Rcpp::Rcout << std::endl;
    //Rcpp::Rcout << "post_prob: ";
    //for(size_t w = 0; w < W; w++) Rcpp::Rcout << " " << post_prob[w];
    //Rcpp::Rcout << std::endl;
    //Rcpp::Rcout << "cdf: ";
    //Rcpp::Rcout.precision(4);
    //for(size_t w = 0; w < W; w++) Rcpp::Rcout << " " << cdf[w];
    //Rcpp::Rcout << std::endl;
    
    double unif = gen.uniform();
    //Rcpp::Rcout << "unif = " << unif << std::endl;
    std::vector<size_t> possible_ix; // all of the indexes for which the cdf is less than unif
    for(size_t w = 0; w < W; w++){
      if(cdf[w] <= unif) possible_ix.push_back(w);
    }
    if(possible_ix.size() == 0) rho_ix = 0; // all non-zero values of cdf are greater than unif. goes to first value
    else if(possible_ix.size() == W) rho_ix = W-1; // unif = 1, which is never going to happen...
    else{
      rho_ix = possible_ix[possible_ix.size()-1] + 1;// cdf[rho_ix] < unif < cdf[rho_ix+1]
    }
    //Rcpp::Rcout << "rho_ix = " << rho_ix << std::endl;
  } // closes if/else checking that W = 1
}

//--------------------------------------------------
// compute the posterior mean and variance for terminal node parameters

void mu_posterior(arma::vec &post_mean, arma::mat &post_cov_chol, double &ll, std::vector<double> sigma, const std::vector<arma::mat> &Omega_vec, std::vector<sinfo> &sv, data_info &di, tree_prior_info &tree_pi)
{
  size_t L = sv.size() ; // number of bottom nodes
  size_t j = 0;
  size_t jj = 0;
  post_mean.set_size(L); // set appropriate size for the mean
  post_cov_chol.set_size(L,L);
  ll = 0.0;
  arma::mat Lambda_inv = arma::zeros<arma::mat>(L,L);
  arma::vec theta = arma::zeros<arma::mat>(L);
  
  // The indexes align with the output of bnv
  
  //Rcpp::Rcout << "[mu_posterior]: About to start loop" << std::endl;
  // populate elements of Lamba_inv and theta
  for(size_t l = 0; l < L; l++){
    for(size_t ll = l; ll < L; ll++){
      //Rcpp::Rcout << "l = " << l << " ll = " << ll << endl;
      if(ll == l) Lambda_inv(l,l) += 1.0/(tree_pi.sigma_mu[di.k] * tree_pi.sigma_mu[di.k]);
      
      // loop over elements j in sv[l].I[i] and j' in sv[ll].I[i]
      for(size_t i = 0; i < di.N; i++){
        if( (sv[l].node_count[i] > 0) && (sv[ll].node_count[i] > 0)){ // contributions arise only if >= 1 observation in leafs l and ll
          for(size_t j_ix = 0; j_ix < sv[l].node_count[i]; j_ix++){
            for(size_t jj_ix = 0; jj_ix < sv[ll].node_count[i]; jj_ix++){
              j = sv[l].I[i][j_ix];
              jj = sv[ll].I[i][jj_ix];
              Lambda_inv(l,ll) += 1.0/(sigma[i] * sigma[i]) * di.x[di.k + di.p*(di.start_index[i] + j)] * Omega_vec[i](j,jj) * di.x[di.k + di.p*(di.start_index[i] + jj)];
              if(Lambda_inv(l,ll) != Lambda_inv(l,ll)){
                Rcpp::Rcout << "l = " << l << " ll = " << ll << " i = " << i << 1.0/(sigma[i] * sigma[i]) * di.x[di.k + di.p*(di.start_index[i] + j)] * Omega_vec[i](j,jj) * di.x[di.k + di.p*(di.start_index[i] + jj)] << endl;
                Rcpp::stop("Something is wrong here");
              }
                
            } // closes inner loop over observations jj who land in leaf ll
          } // closes outer loop over observations j who land in leaf l
        } // closes if checking that there's at least one observation in leaf l and ll
      } // closes loop over all individuals
      Lambda_inv(ll,l) = Lambda_inv(l,ll); // symmetrize the matrix
    } // closes inner loop over ll > l
  } // closes outer loop over the leafs
  // for testing purposes we whould print out Lambda here
  //Rcpp::Rcout << "  Current value of Lambda_inv " << endl;
  //Lambda_inv.print();
  
  arma::mat Lambda = arma::inv_sympd(Lambda_inv); // invert the matrix to get the posterior covariance matrix
  
  // now it's time for theta
  for(size_t l = 0; l < L; l++){
    for(size_t i = 0; i < di.N; i++){
      if(sv[l].node_count[i] > 0){
        for(size_t j_ix = 0; j_ix < sv[l].node_count[i]; j_ix++){
          j = sv[l].I[i][j_ix];
          for(size_t jj = 0; jj < di.n[i]; jj++){
            theta(l) += 1.0/(sigma[i] * sigma[i]) * di.x[di.k + di.p*(di.start_index[i] + j)] * Omega_vec[i](j,jj) * tree_pi.rp[di.start_index[i] + jj];
          } // closes loop over jj indexing *ALL* observations for individual i
        } // closes loop over observations j in leaf l
      } // closes if checking that there are observation in this leaf
    } // closes loop over individuals i
  } // closes loop over the leafs l
  
  //Rcpp::Rcout << "  Current value of Theta " << endl;
  //theta.print();
  post_mean = Lambda * theta;
  //Rcpp::Rcout << "    Current value of post_mean" << endl;
  //post_mean.print();
  // compute the lower cholesky decomposition of Lambda ... this is used to generate the multivariate normal draws
  post_cov_chol = arma::chol(Lambda, "lower");
  // compute the log-determinant of Lambda
  
  double tmp_log_det = 0.0;
  double tmp_log_det_sign = 0.0;
  arma::log_det(tmp_log_det, tmp_log_det_sign, Lambda);
  ll = 0.5 * tmp_log_det + 0.5 * arma::as_scalar(theta.t() * Lambda * theta);
  
}



void mu_posterior(arma::vec &post_mean, arma::mat &post_cov_chol, double &ll, double &sigma, const std::vector<arma::mat> &Omega_vec, std::vector<sinfo> &sv, data_info &di, tree_prior_info &tree_pi)
{
  size_t L = sv.size() ; // number of bottom nodes
  size_t j = 0;
  size_t jj = 0;
  post_mean.set_size(L); // set appropriate size for the mean
  post_cov_chol.set_size(L,L);
  ll = 0.0;
  arma::mat Lambda_inv = arma::zeros<arma::mat>(L,L);
  arma::vec theta = arma::zeros<arma::mat>(L);
  
  // The indexes align with the output of bnv
  
  //Rcpp::Rcout << "[mu_posterior]: About to start loop" << std::endl;

  // populate elements of Lamba_inv and theta
  for(size_t l = 0; l < L; l++){
    for(size_t ll = l; ll < L; ll++){
      if(ll == l) Lambda_inv(l,l) += 1.0/(tree_pi.sigma_mu[di.k] * tree_pi.sigma_mu[di.k]);
        
      // loop over elements j in sv[l].I[i] and j' in sv[ll].I[i]
      for(size_t i = 0; i < di.N; i++){
        //Rcpp::Rcout << " " << i;
        //if(Omega_vec[i].n_rows != di.n[i]) Rcpp::Rcout << "  i = " << i << " n[i] = " << di.n[i] << "  Omega_vec[i] has " << Omega_vec[i].n_rows << "rows" << std::endl;
        if( (sv[l].node_count[i] > 0) && (sv[ll].node_count[i] > 0)){ // contributions arise only if >= 1 observation in leafs l and ll
          for(size_t j_ix = 0; j_ix < sv[l].node_count[i]; j_ix++){
            for(size_t jj_ix = 0; jj_ix < sv[ll].node_count[i]; jj_ix++){
              j = sv[l].I[i][j_ix];
              jj = sv[ll].I[i][jj_ix];
              if( (j >= Omega_vec[i].n_rows) || (jj >= Omega_vec[i].n_cols)){
                Rcpp::Rcout << "    l = " << l << " ll = " << ll << endl;
                Rcpp::Rcout << "i = " << i << " j = " << j << "jj_ix = " << jj_ix << " jj = " << jj << std::endl;
                Rcpp::Rcout << "sv[l].node_count[i] = " << sv[l].node_count[i] << std::endl;
                Rcpp::Rcout << "sv[ll].node_count[i] = " << sv[ll].node_count[i] << std::endl;
                Rcpp::stop("j or jj exceeds number of rows/columns of Omega_vec[i]");
              }
              Lambda_inv(l,ll) += 1.0/(sigma * sigma) * di.x[di.k + di.p*(di.start_index[i] + j)] * Omega_vec[i](j,jj) * di.x[di.k + di.p*(di.start_index[i] + jj)];
              if(Lambda_inv(l,ll) != Lambda_inv(l,ll)){
                Rcpp::Rcout << "l = " << l << " ll = " << ll << " i = " << i << 1.0/(sigma * sigma) * di.x[di.k + di.p*(di.start_index[i] + j)] * Omega_vec[i](j,jj) * di.x[di.k + di.p*(di.start_index[i] + jj)] << endl;
                Rcpp::stop("Something is wrong here");
              }
            } // closes inner loop over observations jj who land in leaf ll
          } // closes outer loop over observations j who land in leaf l
        } // closes if checking that there's at least one observation in leaf l and ll
      } // closes loop over all individuals
      //Rcpp::Rcout << std::endl;
      Lambda_inv(ll,l) = Lambda_inv(l,ll); // symmetrize the matrix
    } // closes inner loop over ll > l
  } // closes outer loop over the leafs
  // for testing purposes we whould print out Lambda here
  //Rcpp::Rcout << "  Current value of Lambda_inv " << endl;
  //Lambda_inv.print();
                  
  arma::mat Lambda = arma::inv_sympd(Lambda_inv); // invert the matrix to get the posterior covariance matrix
  //Rcpp::Rcout << "[mu_posterior]: Computed lambda" << std::endl;
                  
  // now it's time for theta
  for(size_t l = 0; l < L; l++){
    for(size_t i = 0; i < di.N; i++){
      if(sv[l].node_count[i] > 0){
        for(size_t j_ix = 0; j_ix < sv[l].node_count[i]; j_ix++){
          j = sv[l].I[i][j_ix];
          for(size_t jj = 0; jj < di.n[i]; jj++){
            theta(l) += 1.0/(sigma * sigma) * di.x[di.k + di.p*(di.start_index[i] + j)] * Omega_vec[i](j,jj) * tree_pi.rp[di.start_index[i] + jj];
          } // closes loop over jj indexing *ALL* observations for individual i
        } // closes loop over observations j in leaf l
      } // closes if checking that there are observation in this leaf
    } // closes loop over individuals i
  } // closes loop over the leafs l
                  
  //Rcpp::Rcout << "  Current value of Theta " << endl;
  //theta.print();
  //Rcpp::Rcout << "  Current value of Lambda_inv" << std::endl;
  //Lambda_inv.print();
  post_mean = Lambda * theta;
  //Rcpp::Rcout << "    Current value of post_mean" << endl;
  //post_mean.print();
  // compute the lower cholesky decomposition of Lambda ... this is used to generate the multivariate normal draws
  post_cov_chol = arma::chol(Lambda, "lower");
                  // compute the log-determinant of Lambda
                  
  double tmp_log_det = 0.0;
  double tmp_log_det_sign = 0.0;
  arma::log_det(tmp_log_det, tmp_log_det_sign, Lambda);
  ll = 0.5 * tmp_log_det + 0.5 * arma::as_scalar(theta.t() * Lambda * theta);
}

// overloaded for margin_re
void mu_posterior(arma::vec &post_mean, arma::mat &post_cov_chol, double &ll, double &sigma, double &rho, std::vector<sinfo> &sv, data_info &di, tree_prior_info &tree_pi)
{
  size_t L = sv.size() ; // number of bottom nodes
  size_t j = 0;
  size_t jj = 0;
  post_mean.set_size(L); // set appropriate size for the mean
  post_cov_chol.set_size(L,L);
  ll = 0.0;
  arma::mat Lambda_inv = arma::zeros<arma::mat>(L,L);
  arma::vec theta = arma::zeros<arma::mat>(L);
  
  
  std::vector<double> xx_sum(L); // holds sum of x_ijk^2 within each leaf
  std::vector<double> xr_sum(L); // holds sum of x_ijk * rf_ij within each leaf
  std::vector<double> x_sum(L); // holds sum of x_ijk within each leaf
  //std::vector<double> r_sum(L); // holds sum of r_ij within each leaf
  double r_sum = 0.0;// holds the total sum of r_partial
 
  // note that each of these terms needs to be updated for every individual i
  
  for(size_t i = 0; i < di.N; i++){
    r_sum = 0.0;
    for(size_t l = 0; l < L; l++){
      xx_sum[l] = 0.0;
      xr_sum[l] = 0.0;
      x_sum[l] = 0.0;
      
      if(sv[l].node_count[i] > 0){
        for(size_t j_ix = 0; j_ix < sv[l].node_count[i]; j_ix++){
          j = sv[l].I[i][j_ix];
          xx_sum[l] += di.x[di.k + di.p * (di.start_index[i] + j)] * di.x[di.k + di.p * (di.start_index[i] + j)];
          xr_sum[l] += di.x[di.k + di.p * (di.start_index[i] + j)] * tree_pi.rp[di.start_index[i] + j];
          x_sum[l] += di.x[di.k + di.p * (di.start_index[i] + j)];
          r_sum += tree_pi.rp[di.start_index[i] + j];
        }
      }
    }
    // we now have enough to compute X_i' Omega X_i and X_i' Omega r_i for individual i.
    // we can therefore update elements of Lambda_inv now
    for(size_t l = 0; l < L; l++){
      Lambda_inv(l,l) += 1.0/(sigma * sigma) * xx_sum[l] * 1.0/ (1.0 - rho);
      for(size_t ll = 0; ll < L; ll++){
        Lambda_inv(l,ll) -= 1.0/(sigma * sigma) * rho / ((1.0 - rho) * (1.0 + ( (double) di.n[i] - 1.0) * rho)) * x_sum[l] * x_sum[ll];
      }
      theta(l) += 1.0/(sigma * sigma) * 1.0 / (1.0 - rho) * xr_sum[l] - 1.0/(sigma * sigma) * rho / ((1.0 - rho) * (1.0 + ( (double) di.n[i] - 1.0) * rho)) * x_sum[l] * r_sum;
    }
  }
  Lambda_inv.diag() += 1.0/(tree_pi.sigma_mu[di.k] * tree_pi.sigma_mu[di.k]);
  
  if(Lambda_inv.has_nan()) Rcpp::Rcout << "Lambda_inv has an nan" << std::endl;
  
  //Rcpp::Rcout << "theta = " << std::endl;
  //theta.print();
  
  //Rcpp::Rcout << "Lambda_inv = " << std::endl;
  //Lambda_inv.print();
  
  arma::mat Lambda = arma::inv_sympd(Lambda_inv);
  post_mean = Lambda * theta;

  //post_cov_chol = arma::chol(Lambda, "lower");
  bool chol_success = arma::chol(post_cov_chol, Lambda, "lower");
  if(chol_success == false){
    Rcpp::Rcout << "cholesky decomposition failed!" << std::endl;
    Lambda.print();
  }
  
  double tmp_log_det = 0.0;
  double tmp_log_det_sign = 0.0;
  arma::log_det(tmp_log_det, tmp_log_det_sign, Lambda);
  ll = 0.5 * tmp_log_det + 0.5 * arma::as_scalar(theta.t() * Lambda * theta);
}






void mu_posterior(double &post_mean, double &post_var, double &ll, double &sigma, sinfo &si, data_info &di, tree_prior_info &tree_pi)
{
  double V_inv = 1.0/(tree_pi.sigma_mu[di.k] * tree_pi.sigma_mu[di.k]);

  double theta = 0.0;
  size_t j = 0;

  for(size_t i = 0; i < di.N; i++){
    if(si.node_count[i] > 0){
      // only need to do stuff if some observations from individual i land in the leaf we care about
      for(size_t j_ix = 0; j_ix < si.node_count[i]; j_ix++){
        j = si.I[i][j_ix]; // get the index of the observation from individual i that lands in this leaf
        V_inv += 1.0/(sigma * sigma) * di.x[di.k + di.p*(di.start_index[i] + j)] * di.x[di.k + di.p*(di.start_index[i] + j)];
        theta += 1.0/(sigma * sigma) * di.x[di.k + di.p*(di.start_index[i] + j)] * tree_pi.rp[j + di.start_index[i]];
      }
    }
  }
  post_var = 1.0/V_inv;
  post_mean = post_var * theta;
  ll = 0.5 * log(post_var)  + 0.5 * (post_mean * post_mean)/post_var;

}
*/

/*
void drmu(tree &t, double &sigma, xinfo &xi, data_info &di, tree_prior_info &tree_pi, RNG &gen)
{
  tree::npv bnv; // bottom nodes
  std::vector<sinfo> sv;
  allsuff(t, xi, di, bnv, sv);
  
  double mu_bar = 0.0;
  double V = 0.0;
  double ll = 0.0;
  for(tree::npv::size_type i = 0; i != bnv.size(); i++){
    mu_bar = 0.0;
    V = 0.0;
    mu_posterior(mu_bar, V, ll, sigma, sv[i],di, tree_pi);
    bnv[i]->setm(mu_bar + sqrt(V) * gen.normal());
    if(bnv[i]->getm() != bnv[i]->getm()){
      Rcpp::stop("nan in drmu");
    }
  }
}
*/
/*
void update_alpha(double* alpha, double &sigma, double &sigma_alpha, data_info &di, RNG &gen)
{
  
  //This is probably bad practice but I'm over-loading the full residuals
   //Immediately before the function is called, both allfits and r_full have been adjusted
   //to remove the previous value of alpha
  double alpha_mean = 0.0;
  double alpha_var_inv = 1.0/(sigma_alpha * sigma_alpha);
  double alpha_var = 1.0;
  
  for(size_t i = 0; i < di.N; i++){
    alpha_mean = 0.0;
    alpha_var_inv = 1.0/(sigma_alpha * sigma_alpha);
    for(size_t j = 0; j < di.n[i]; j++){
      alpha_mean += 1.0/(sigma * sigma) * di.rf[di.start_index[i] + j];
      alpha_var_inv += 1.0/(sigma * sigma);
    }
    alpha_var = 1.0/alpha_var_inv;
    alpha_mean *= alpha_var;
    alpha[i] = alpha_mean + sqrt(alpha_var) * gen.normal();
  }
}
*/
