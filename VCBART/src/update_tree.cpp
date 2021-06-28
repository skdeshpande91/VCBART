//
//  update_tree.cpp
//    We will use this instead of bd_new*
//
//    For now, we will make a copy of the tree (xx).
//    The proposal tree will be y.
//  Created by Sameer Deshpande on 4/7/20.
//

#include "update_tree.h"
void update_tree_ind(tree &x, double &sigma, std::vector<double> &theta, std::vector<size_t> &var_counts, xinfo &xi, data_info &di, tree_prior_info &tree_pi, RNG &gen)
{
  
  tree::npv goodbots_x; // nodes that we can split in x, original tree
  double PBx = getpb(x, xi, tree_pi, goodbots_x); // probability that we propose a birth given x
  
  

  // initialize a bunch of stuff that is common to GROW and PRUNE moves
  
  tree::npv bnv_x; // bottom nodes of x, original tree
  tree::npv bnv_xx; // bottom nodes of xx, copy of original tree
  tree::npv bnv_y; // bottom nodes of y, proposed tree
  
  std::vector<sinfo> sv_xx; // sufficient statistics for xx, copy of original tree
  std::vector<sinfo> sv_y; // sufficient statistics for y, proposed tree
  
  arma::mat post_cov_chol_xx; // lowere cholesky of post. covariance of mu, conditional on y the proposed tree
  arma::mat post_cov_chol_y; // lower cholesky of post. covariance of mu, conditional on y the proposed tree
  
  
  double ll_xx; // log marginal likelihood given xx, copy of original tree
  double ll_y; // log marginal likelihood given y, proposed tree
  
  arma::vec post_mean_xx; // posterior mean of mu given xx, copy of original tree
  arma::vec post_mean_y; // posterior mean of mu given y, proposed tree
  
  double alpha1;
  double alpha2;
  double alpha;
  
  arma::vec norm_samples; // holds standard normal samples
  arma::vec mu_samples; // holds the new sampled mus
  
  if(gen.uniform() < PBx){
    // try a GROW proposal
    size_t ni = floor(gen.uniform() * goodbots_x.size());
    tree::tree_p nx = goodbots_x[ni]; // randomly select a bottom node split
    
    std::vector<size_t> goodvars; // variables nx can non-trivially split on
    getgoodvars(nx, xi, goodvars); // collection of variables that yield non-trivial splits
    
    size_t v = draw_multinomial(di.R, theta, gen); // draw the split index
    
    int L,U;
    L = 0;
    U = xi[v].size() - 1;
    nx->rg(v, &L, &U);
    if(U < L){
      // we have a trivial split
      L = 0;
      U = xi[v].size() - 1;
    }
    size_t c = L + floor(gen.uniform() * (U - L + 1));
    
    double Pbotx = 1.0/goodbots_x.size();
    size_t dnx = nx->depth(); // depth of node nx
    double PGnx = tree_pi.alpha[di.k]/pow(1.0 + dnx, tree_pi.beta[di.k]);
    double PGly, PGry;
    if(goodvars.size() > 1){ // we know there are variables we could split l and r on
      PGly = tree_pi.alpha[di.k]/pow(1.0 + 1.0 + dnx, tree_pi.beta[di.k]); // [SKD]: 28 June 2021 -- we were missing one of the 1.0's
      PGry = PGly;
    } else{ // only have v to work with, if it is exhausted at either child need PG = 0
      if( (int)(c-1) < L){ // v exhausted in the new left child l, new upper limit would be c-1
        PGly = 0.0;
      } else{
        PGly = tree_pi.alpha[di.k]/pow(1.0 + 1.0 + dnx, tree_pi.beta[di.k]);
      }
      if(U < (int)(c+1)){ // v exhausted in new right child, new lower limit would be c+1
        PGry = 0.0;
      } else{
        PGry = tree_pi.alpha[di.k]/pow(1.0 + 1.0 + (double) dnx, tree_pi.beta[di.k]);
      }
    }
    double PDy; // prob of proposing death at y
    if(goodbots_x.size() > 1){ // can birth at y as there are splittable nodes left
      PDy = 1.0 - tree_pi.pb;
    } else{ // nx was the only node we could split on
      if( (PGry == 0) && (PGly == 0)){ // cannot birth at y
        PDy = 1.0;
      } else{ // y can birth at either l or r
        PDy = 1.0 - tree_pi.pb;
      }
    }
    
    double Pnogy; // death prob of choosing the nog node nx in y
    size_t nnogs = x.nnogs();
    tree::tree_cp nxp = nx->getp();
    if(nxp == 0){ // no parent, nx is thetop and only node
      Pnogy = 1.0;
    } else{
      if(nxp->isnog()){ // if parent is a nog, number of nogs same at x and y
        Pnogy = 1.0/nnogs;
      } else{ // if parent is not a nog, y has one more nog.
        Pnogy = 1.0/(nnogs + 1.0);
      }
    }
    
    // Time to compute log-likelihoods
    tree xx = tree(x); // a copy of x
    tree y = tree(x); // a copy of x, which we will grow.
    
    y.birth(nx->nid(), v, c, 0.0, 0.0); // perfor the birth
    
    allsuff(xx, xi, di, bnv_xx, sv_xx); // note bnv_xx gets updated within allsuff
    allsuff(y, xi, di, bnv_y, sv_y); // note bnv_y gets updated within allsuff
    
    // let's compute the log-likelihood for the original tree. also we will get the posterior mean and covariance
    
    mu_posterior_ind(post_mean_xx, post_cov_chol_xx, ll_xx, sigma,sv_xx, di, tree_pi);
    
    bool flag = true; // flag == true means we can continue with the birth proposal
    for(size_t l = 0; l != bnv_y.size(); l++){
      if(sv_y[l].n < 5){
        flag = false;
        break;
      }
    }
    if(flag == true){
      mu_posterior_ind(post_mean_y, post_cov_chol_y, ll_y, sigma, sv_y, di, tree_pi);
      alpha1 = (PGnx*(1.0-PGly)*(1.0-PGry)*PDy*Pnogy)/((1.0-PGnx)*PBx*Pbotx);
      alpha2 = alpha1 * exp(ll_y - ll_xx) * 1.0/(tree_pi.tau[di.k]);
      alpha = std::min(1.0, alpha2);
      if(gen.uniform() < alpha){
        // we accept the GROW proposal
        var_counts[v]++; // update the running counts of how many times v was selected
        x.birth(nx->nid(),v,c,0.0,0.0);
        bnv_x.clear();
        x.getbots(bnv_x);
        if(bnv_x.size() != bnv_y.size()) Rcpp::stop("[bd]: after accepting birth bnv_x and bnv_y have different sizes!");
        norm_samples.set_size(bnv_x.size());
        mu_samples.set_size(bnv_x.size());
        for(size_t l = 0; l < bnv_x.size(); l++) norm_samples(l) = gen.normal();
        mu_samples = post_mean_y + post_cov_chol_y * norm_samples;
        for(size_t l = 0; l < bnv_x.size(); l++) bnv_x[l]->setm(mu_samples(l));
      } else{
        // we reject GROW proposal
        bnv_x.clear();
        x.getbots(bnv_x);
        if(bnv_x.size() != bnv_xx.size()) Rcpp::stop("[bd]: After rejecting GROW, bnv_x and bnv_xx have different sizes!");
        norm_samples.set_size(bnv_x.size());
        mu_samples.set_size(bnv_x.size());
        for(size_t l = 0; l < bnv_x.size(); l++) norm_samples(l) = gen.normal();
        mu_samples = post_mean_xx + post_cov_chol_xx * norm_samples;
        for(size_t l = 0; l < bnv_x.size(); l++) bnv_x[l]->setm(mu_samples(l));
      } // closes if/else checking whether we accept/reject the GROW proposal
    } // closes if/else checking that all terminal nodes have at least 5 observations
  } else{
    // try a PRUNE proposal
    tree::npv nognds; // nog nods
    x.getnogs(nognds);
    size_t ni = floor(gen.uniform() * nognds.size());
    tree::tree_p nx = nognds[ni]; // the nog node we might kill children at
    
    // compute stuff for metropolis-hastings ratio
    double PGny; // prob the nog node grows
    size_t dny = nx->depth();
    PGny = tree_pi.alpha[di.k]/pow(1.0+dny,tree_pi.beta[di.k]);
    //Rcpp::Rcout << "[bd]:    PGny = " << PGny << endl;
    
    //better way to code these two?
    double PGlx = pgrow(nx->getl(),xi,tree_pi, di.k);
    double PGrx = pgrow(nx->getr(),xi,tree_pi, di.k);
    
    double PBy;  //prob of birth move at y
                 //if(nx->ntype()=='t') { //is the nog node nx the top node
    if(!(nx->p)) { //is the nog node nx the top node
      PBy = 1.0;
    } else {
      PBy = tree_pi.pb;
    }
    
    double Pboty;  //prob of choosing the nog as bot to split on when y
    int ngood = goodbots_x.size();
    if(cansplit(nx->getl(),xi)) --ngood; //if can split at left child, lose this one
    if(cansplit(nx->getr(),xi)) --ngood; //if can split at right child, lose this one
    ++ngood;  //know you can split at nx
    Pboty=1.0/ngood;
    
    double PDx = 1.0-PBx; //prob of a death step at x
    double Pnogx = 1.0/nognds.size();
    
    // what was the split index of nx
    size_t v = nx->v;
    
    
    // now actually perform the death
    tree xx = tree(x); // copy of the original tree
    tree y = tree(x); // will run the death on this copy
    
    y.death(nx->nid(),0.0); // this is the proposal
    
    allsuff(xx, xi, di, bnv_xx, sv_xx); // note bnv_xx gets updated within allsuff
    allsuff(y, xi, di, bnv_y, sv_y); // note bnv_y gets updated within allsuff
    
    // let's compute the log-likelihood for the original tree. also we will get the posterior mean and covariance
    mu_posterior_ind(post_mean_xx, post_cov_chol_xx, ll_xx, sigma, sv_xx, di, tree_pi);
    mu_posterior_ind(post_mean_y, post_cov_chol_y, ll_y, sigma, sv_y, di, tree_pi);
    
    alpha1 = ((1.0-PGny)*PBy*Pboty)/(PGny*(1.0-PGlx)*(1.0-PGrx)*PDx*Pnogx);
    alpha2 = alpha1 * exp(ll_y - ll_xx) * tree_pi.tau[di.k];
    alpha = std::min(1.0, alpha2);
    if(gen.uniform() < alpha){
      // accept the PRUNE proposal
      x.death(nx->nid(),0.0);
      
      var_counts[v]--; // decrement the counter for variable v.
      
      
      bnv_x.clear();
      x.getbots(bnv_x); // get all of the bottom nodes of x. these will be arranged in the same order as bnv_y!
      if(bnv_x.size() != bnv_y.size()) Rcpp::stop("[bd]: after accepting PRUNE, bnv_x and bnv_y have different sizes!");
      
      norm_samples.set_size(bnv_x.size());
      mu_samples.set_size(bnv_x.size());
      for(size_t l = 0; l < bnv_x.size(); l++) norm_samples(l) = gen.normal();
      mu_samples = post_mean_y + post_cov_chol_y * norm_samples;
      // assign the mu parameters for x
      for(size_t l = 0; l < bnv_x.size(); l++) bnv_x[l]->setm(mu_samples(l));
    } else{
      // reject the death proposal
      bnv_x.clear();
      x.getbots(bnv_x);
      if(bnv_x.size() != bnv_xx.size()) Rcpp::stop("[bd]: bnv_x and bnv_xx have different sizes!");
      
      norm_samples.set_size(bnv_x.size());
      mu_samples.set_size(bnv_x.size());
      for(size_t l = 0; l < bnv_x.size(); l++) norm_samples(l) = gen.normal();
      mu_samples = post_mean_xx + post_cov_chol_xx * norm_samples;
      for(size_t l = 0; l < bnv_x.size(); l++) bnv_x[l]->setm(mu_samples(l));
    } // closes if/else checking whether the death proposal is accepted
  }
}



// version of update tree for compound symmetry error structure
void update_tree_cs(tree &x, double &sigma, double &rho, std::vector<double> &theta, std::vector<size_t> &var_counts, xinfo &xi, data_info &di, tree_prior_info &tree_pi, RNG &gen)
{
  
  tree::npv goodbots_x; // nodes that we can split in x, original tree
  double PBx = getpb(x, xi, tree_pi, goodbots_x); // probability that we propose a birth given x
  
  
  
  // initialize a bunch of stuff that is common to GROW and PRUNE moves
  
  tree::npv bnv_x; // bottom nodes of x, original tree
  tree::npv bnv_xx; // bottom nodes of xx, copy of original tree
  tree::npv bnv_y; // bottom nodes of y, proposed tree
  
  std::vector<sinfo> sv_xx; // sufficient statistics for xx, copy of original tree
  std::vector<sinfo> sv_y; // sufficient statistics for y, proposed tree
  
  arma::mat post_cov_chol_xx; // lowere cholesky of post. covariance of mu, conditional on y the proposed tree
  arma::mat post_cov_chol_y; // lower cholesky of post. covariance of mu, conditional on y the proposed tree
  
  
  double ll_xx; // log marginal likelihood given xx, copy of original tree
  double ll_y; // log marginal likelihood given y, proposed tree
  
  arma::vec post_mean_xx; // posterior mean of mu given xx, copy of original tree
  arma::vec post_mean_y; // posterior mean of mu given y, proposed tree
  
  double alpha1;
  double alpha2;
  double alpha;
  
  arma::vec norm_samples; // holds standard normal samples
  arma::vec mu_samples; // holds the new sampled mus
  
  if(gen.uniform() < PBx){
    // try a GROW proposal
    size_t ni = floor(gen.uniform() * goodbots_x.size());
    tree::tree_p nx = goodbots_x[ni]; // randomly select a bottom node split
    
    std::vector<size_t> goodvars; // variables nx can non-trivially split on
    getgoodvars(nx, xi, goodvars); // collection of variables that yield non-trivial splits
    
    size_t v = draw_multinomial(di.R, theta, gen); // draw the split index
    
    int L,U;
    L = 0;
    U = xi[v].size() - 1;
    nx->rg(v, &L, &U);
    if(U < L){
      // we have a trivial split
      L = 0;
      U = xi[v].size() - 1;
    }
    size_t c = L + floor(gen.uniform() * (U - L + 1));
    
    double Pbotx = 1.0/goodbots_x.size();
    size_t dnx = nx->depth(); // depth of node nx
    double PGnx = tree_pi.alpha[di.k]/pow(1.0 + dnx, tree_pi.beta[di.k]);
    double PGly, PGry;
    if(goodvars.size() > 1){ // we know there are variables we could split l and r on
      PGly = tree_pi.alpha[di.k]/pow(1.0 + 1.0 + dnx, tree_pi.beta[di.k]); // [SKD]: 28 June 2021 -- we were missing one of the 1.0's
      PGry = PGly;
    } else{ // only have v to work with, if it is exhausted at either child need PG = 0
      if( (int)(c-1) < L){ // v exhausted in the new left child l, new upper limit would be c-1
        PGly = 0.0;
      } else{
        PGly = tree_pi.alpha[di.k]/pow(1.0 + 1.0 + dnx, tree_pi.beta[di.k]);
      }
      if(U < (int)(c+1)){ // v exhausted in new right child, new lower limit would be c+1
        PGry = 0.0;
      } else{
        PGry = tree_pi.alpha[di.k]/pow(1.0 + 1.0 + (double) dnx, tree_pi.beta[di.k]);
      }
    }
    double PDy; // prob of proposing death at y
    if(goodbots_x.size() > 1){ // can birth at y as there are splittable nodes left
      PDy = 1.0 - tree_pi.pb;
    } else{ // nx was the only node we could split on
      if( (PGry == 0) && (PGly == 0)){ // cannot birth at y
        PDy = 1.0;
      } else{ // y can birth at either l or r
        PDy = 1.0 - tree_pi.pb;
      }
    }
    
    double Pnogy; // death prob of choosing the nog node nx in y
    size_t nnogs = x.nnogs();
    tree::tree_cp nxp = nx->getp();
    if(nxp == 0){ // no parent, nx is thetop and only node
      Pnogy = 1.0;
    } else{
      if(nxp->isnog()){ // if parent is a nog, number of nogs same at x and y
        Pnogy = 1.0/nnogs;
      } else{ // if parent is not a nog, y has one more nog.
        Pnogy = 1.0/(nnogs + 1.0);
      }
    }
    
    // Time to compute log-likelihoods
    tree xx = tree(x); // a copy of x
    tree y = tree(x); // a copy of x, which we will grow.
    
    y.birth(nx->nid(), v, c, 0.0, 0.0); // perfor the birth
    
    allsuff(xx, xi, di, bnv_xx, sv_xx); // note bnv_xx gets updated within allsuff
    allsuff(y, xi, di, bnv_y, sv_y); // note bnv_y gets updated within allsuff
    
    // let's compute the log-likelihood for the original tree. also we will get the posterior mean and covariance
    
    mu_posterior_cs(post_mean_xx, post_cov_chol_xx, ll_xx, sigma, rho, sv_xx, di, tree_pi);
    bool flag = true; // flag == true means we can continue with the birth proposal
    for(size_t l = 0; l != bnv_y.size(); l++){
      if(sv_y[l].n < 5){
        flag = false;
        break;
      }
    }
    if(flag == true){
      mu_posterior_cs(post_mean_y, post_cov_chol_y, ll_y, sigma, rho, sv_y, di, tree_pi);
      alpha1 = (PGnx*(1.0-PGly)*(1.0-PGry)*PDy*Pnogy)/((1.0-PGnx)*PBx*Pbotx);
      alpha2 = alpha1 * exp(ll_y - ll_xx) * 1.0/(tree_pi.tau[di.k]);
      alpha = std::min(1.0, alpha2);
      if(gen.uniform() < alpha){
        // we accept the GROW proposal
        var_counts[v]++; // update the running counts of how many times v was selected
        x.birth(nx->nid(),v,c,0.0,0.0);
        bnv_x.clear();
        x.getbots(bnv_x);
        if(bnv_x.size() != bnv_y.size()) Rcpp::stop("[update_tree]: after accepting birth bnv_x and bnv_y have different sizes!");
        norm_samples.set_size(bnv_x.size());
        mu_samples.set_size(bnv_x.size());
        for(size_t l = 0; l < bnv_x.size(); l++) norm_samples(l) = gen.normal();
        mu_samples = post_mean_y + post_cov_chol_y * norm_samples;
        for(size_t l = 0; l < bnv_x.size(); l++) bnv_x[l]->setm(mu_samples(l));
      } else{
        // we reject GROW proposal
        bnv_x.clear();
        x.getbots(bnv_x);
        if(bnv_x.size() != bnv_xx.size()) Rcpp::stop("[update_tree]: After rejecting GROW, bnv_x and bnv_xx have different sizes!");
        norm_samples.set_size(bnv_x.size());
        mu_samples.set_size(bnv_x.size());
        for(size_t l = 0; l < bnv_x.size(); l++) norm_samples(l) = gen.normal();
        mu_samples = post_mean_xx + post_cov_chol_xx * norm_samples;
        for(size_t l = 0; l < bnv_x.size(); l++) bnv_x[l]->setm(mu_samples(l));
      } // closes if/else checking whether we accept/reject the GROW proposal
    } // closes if/else checking that all terminal nodes have at least 5 observations
  } else{
    // try a PRUNE proposal
    tree::npv nognds; // nog nods
    x.getnogs(nognds);
    size_t ni = floor(gen.uniform() * nognds.size());
    tree::tree_p nx = nognds[ni]; // the nog node we might kill children at
    
    // compute stuff for metropolis-hastings ratio
    double PGny; // prob the nog node grows
    size_t dny = nx->depth();
    PGny = tree_pi.alpha[di.k]/pow(1.0+dny,tree_pi.beta[di.k]);
    //Rcpp::Rcout << "[bd]:    PGny = " << PGny << endl;
    
    //better way to code these two?
    double PGlx = pgrow(nx->getl(),xi,tree_pi, di.k);
    double PGrx = pgrow(nx->getr(),xi,tree_pi, di.k);
    
    double PBy;  //prob of birth move at y
                 //if(nx->ntype()=='t') { //is the nog node nx the top node
    if(!(nx->p)) { //is the nog node nx the top node
      PBy = 1.0;
    } else {
      PBy = tree_pi.pb;
    }
    
    double Pboty;  //prob of choosing the nog as bot to split on when y
    int ngood = goodbots_x.size();
    if(cansplit(nx->getl(),xi)) --ngood; //if can split at left child, lose this one
    if(cansplit(nx->getr(),xi)) --ngood; //if can split at right child, lose this one
    ++ngood;  //know you can split at nx
    Pboty=1.0/ngood;
    
    double PDx = 1.0-PBx; //prob of a death step at x
    double Pnogx = 1.0/nognds.size();
    
    // what was the split index of nx
    size_t v = nx->v;
    
    
    // now actually perform the death
    tree xx = tree(x); // copy of the original tree
    tree y = tree(x); // will run the death on this copy
    
    y.death(nx->nid(),0.0); // this is the proposal
    
    allsuff(xx, xi, di, bnv_xx, sv_xx); // note bnv_xx gets updated within allsuff
    allsuff(y, xi, di, bnv_y, sv_y); // note bnv_y gets updated within allsuff
    
    // let's compute the log-likelihood for the original tree. also we will get the posterior mean and covariance
    mu_posterior_cs(post_mean_xx, post_cov_chol_xx, ll_xx, sigma, rho,  sv_xx, di, tree_pi);
    mu_posterior_cs(post_mean_y, post_cov_chol_y, ll_y, sigma, rho, sv_y, di, tree_pi);
    
    alpha1 = ((1.0-PGny)*PBy*Pboty)/(PGny*(1.0-PGlx)*(1.0-PGrx)*PDx*Pnogx);
    alpha2 = alpha1 * exp(ll_y - ll_xx) * tree_pi.tau[di.k];
    alpha = std::min(1.0, alpha2);
    if(gen.uniform() < alpha){
      // accept the PRUNE proposal
      x.death(nx->nid(),0.0);
      
      var_counts[v]--; // decrement the counter for variable v.
      
      
      bnv_x.clear();
      x.getbots(bnv_x); // get all of the bottom nodes of x. these will be arranged in the same order as bnv_y!
      if(bnv_x.size() != bnv_y.size()) Rcpp::stop("[update_tree]: after accepting PRUNE, bnv_x and bnv_y have different sizes!");
      
      norm_samples.set_size(bnv_x.size());
      mu_samples.set_size(bnv_x.size());
      for(size_t l = 0; l < bnv_x.size(); l++) norm_samples(l) = gen.normal();
      mu_samples = post_mean_y + post_cov_chol_y * norm_samples;
      // assign the mu parameters for x
      for(size_t l = 0; l < bnv_x.size(); l++) bnv_x[l]->setm(mu_samples(l));
    } else{
      // reject the death proposal
      bnv_x.clear();
      x.getbots(bnv_x);
      if(bnv_x.size() != bnv_xx.size()) Rcpp::stop("[update_tree]: bnv_x and bnv_xx have different sizes!");
      
      norm_samples.set_size(bnv_x.size());
      mu_samples.set_size(bnv_x.size());
      for(size_t l = 0; l < bnv_x.size(); l++) norm_samples(l) = gen.normal();
      mu_samples = post_mean_xx + post_cov_chol_xx * norm_samples;
      for(size_t l = 0; l < bnv_x.size(); l++) bnv_x[l]->setm(mu_samples(l));
    } // closes if/else checking whether the death proposal is accepted
  }
}

