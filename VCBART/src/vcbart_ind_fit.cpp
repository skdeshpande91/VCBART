#include "update_trees.h"
// [[Rcpp::export(".vcbart_ind_fit")]]
Rcpp::List vcbart_ind_fit(Rcpp::NumericVector Y_train,
                          Rcpp::IntegerVector subj_id_train,
                          Rcpp::IntegerVector ni_train,
                          Rcpp::NumericMatrix tX_train,
                          Rcpp::NumericMatrix tZ_cont_train,
                          Rcpp::IntegerMatrix tZ_cat_train,
                          Rcpp::IntegerVector subj_id_test,
                          Rcpp::IntegerVector ni_test,
                          Rcpp::NumericMatrix tX_test,
                          Rcpp::NumericMatrix tZ_cont_test,
                          Rcpp::IntegerMatrix tZ_cat_test,
                          Rcpp::LogicalVector unif_cuts,
                          Rcpp::Nullable<Rcpp::List> cutpoints_list,
                          Rcpp::Nullable<Rcpp::List> cat_levels_list,
                          Rcpp::Nullable<Rcpp::List> edge_mat_list,
                          Rcpp::LogicalVector graph_split, int graph_cut_type,
                          bool rc_split, double prob_rc, double a_rc, double b_rc,
                          bool sparse, double a_u, double b_u,
                          Rcpp::NumericVector mu0, Rcpp::NumericVector tau,
                          double lambda, double nu,
                          int M,
                          int nd, int burn, int thin,
                          bool save_samples,
                          bool save_trees,
                          bool verbose, int print_every)
{
  Rcpp::RNGScope scope;
  RNG gen;
  
  set_str_conversion set_str;
  
  // Preprocessing stuff
  int N_train = 0;
  int n_train = 0;
  int p = 0;
  int R_cont = 0;
  int R_cat = 0;
  int R = 0;
  
  int N_test = 0;
  int n_test = 0;
  
  parse_training_data(N_train, n_train, p, R, R_cont, R_cat, ni_train, tX_train, tZ_cont_train, tZ_cat_train);
  if(Y_train.size() != N_train) Rcpp::stop("Number of observations in Y_train does not match number of rows in matrix of training X's");
  parse_testing_data(N_test, n_test, p, R_cont, R_cat, ni_test, tX_test, tZ_cont_test, tZ_cat_test);
  
  
  Rcpp::Rcout << " p = " << p << " R_cont = " << R_cont << " R_cat = " << R_cat << std::endl;
  Rcpp::Rcout << "N_train = " << N_train << " n_train = " << n_train << std::endl;
  Rcpp::Rcout << "N_test = " << N_test << " n_test = " << n_test << std::endl;
  
  std::vector<std::set<double>> cutpoints;
  if(R_cont > 0){
    if(cutpoints_list.isNotNull()){
      Rcpp::List tmp_cutpoints = Rcpp::List(cutpoints_list);
      parse_cutpoints(cutpoints, R_cont, tmp_cutpoints, unif_cuts);
    }
  }
  
  std::vector<std::set<int>> cat_levels;
  std::vector<int> K;
  std::vector<std::vector<edge>> edges;
  
  if(R_cat > 0){
    if(cat_levels_list.isNotNull()){
      Rcpp::List tmp_cat_levels = Rcpp::List(cat_levels_list);
      parse_cat_levels(cat_levels,K, R_cat, tmp_cat_levels);
    } else{
      Rcpp::stop("Must provide list of categorical levels");
    }
    
    if(edge_mat_list.isNotNull()){
      Rcpp::List tmp_edge_mat = Rcpp::List(edge_mat_list);
      parse_graphs(edges, R_cat, K, tmp_edge_mat, graph_split);
    }
  }
  
  double* allfit_train = new double[N_train]; // holds fit of regression function
  double* beta_fit_train = new double[p*N_train]; // holds fit of each coefficient function. beta_fit_train[j + i * p] is for obs i coefficient j
  double* residual = new double[N_train];
  
  double* r_sum = new double[n_train];
  double* r2_sum = new double[n_train];
  
  /*
  double* tmp_r_sum = new double[n_train];
  double* tmp_r2_sum = new double[n_train];
  */
  for(int subj_ix = 0; subj_ix < n_train; subj_ix++){
    r_sum[subj_ix] = 0.0;
    r2_sum[subj_ix] = 0.0;
  }
  
  for(int i = 0; i < N_train; i++){
    allfit_train[i] = 0.0; // all trees start off with 0 prediction
    for(int j = 0; j < p; j++) beta_fit_train[j + i*p] = 0.0;
  }
  
  data_info di_train;
  di_train.N = N_train;
  di_train.n = n_train;
  di_train.p = p;
  di_train.R = R;
  di_train.R_cont = R_cont;
  di_train.R_cat = R_cat;
  di_train.subj_id = subj_id_train.begin();
  di_train.ni = ni_train.begin();
  di_train.x = tX_train.begin();
  if(R_cont > 0) di_train.z_cont = tZ_cont_train.begin();
  if(R_cat > 0) di_train.z_cat = tZ_cat_train.begin();
  di_train.rp = residual;
  di_train.r_sum = r_sum;
  di_train.r2_sum = r2_sum;
  
  int tmp_N_test = 0;
  if(n_test > 0) tmp_N_test = N_test;
  else tmp_N_test = 1;
  
  double* allfit_test = new double[tmp_N_test];
  double* beta_fit_test = new double[p * tmp_N_test];
  
  if(n_test > 0){
    for(int i = 0; i < N_test; i++){
      allfit_test[i] = 0.0;
      for(int j = 0; j < p; j++) beta_fit_test[j + i * p] = 0.0;
    }
  }
  
  data_info di_test;
  if(n_test > 0){
    di_test.N = N_test;
    di_test.n = n_test;
    di_test.p = p;
    di_test.R = R;
    di_test.R_cont = R_cont;
    di_test.R_cat = R_cat;
    di_test.subj_id = subj_id_test.begin();
    di_test.ni = ni_test.begin();
    di_test.x = tX_test.begin();
    if(R_cont > 0) di_test.z_cont = tZ_cont_test.begin();
    if(R_cat > 0) di_test.z_cat = tZ_cat_test.begin();
  }
  
  // declare stuff for variable selection
  std::vector<std::vector<double>> theta(p, std::vector<double>(R, 1.0/( (double) R)));
  double* u = new double[p];
  for(int j = 0; j < p; j++) u[j] = 1.0/(1.0 + (double) R); // initialize u
  std::vector<std::vector<int>> var_count(p, std::vector<int>(R, 0));
  // following objects are for random combination rules
  int* rule_count = new int[p]; // how many rules are there in each ensemble
  for(int j = 0; j < p; j++) rule_count[j] = 0;
  int* rc_rule_count = new int[p]; // how many random combination rules are there in each ensemble?
  int* rc_var_count = new int[p]; // how many total variables get used in the random combination rules in each ensemble
  double* theta_rc = new double[p];
  if(R_cont >= 2 && rc_split){
    for(int j = 0; j < p; j++) theta_rc[j] = 2.0/( (double) R_cont);
  }
  
  // now ready to initialize the tree_prior_info objects for each ensemble
  std::vector<tree_prior_info> tree_pi_vec(p);
  for(int j = 0; j < p; j++){
    tree_pi_vec[j].theta = &(theta[j]);
    tree_pi_vec[j].var_count = &(var_count[j]);
    tree_pi_vec[j].rule_count = &(rule_count[j]);
    tree_pi_vec[j].unif_cuts = unif_cuts.begin();
    
    if(R_cont > 0){
      tree_pi_vec[j].unif_cuts = unif_cuts.begin();
      tree_pi_vec[j].cutpoints = &cutpoints;
      if(rc_split){
        tree_pi_vec[j].rc_split = rc_split;
        tree_pi_vec[j].prob_rc = &prob_rc;
        tree_pi_vec[j].theta_rc = &(theta_rc[j]);
        tree_pi_vec[j].rc_var_count = &(rc_var_count[j]);
        tree_pi_vec[j].rc_rule_count = &(rc_rule_count[j]);
      }
    }
    if(R_cat > 0){
      tree_pi_vec[j].cat_levels = &cat_levels;
      tree_pi_vec[j].edges = &edges;
      tree_pi_vec[j].K = &K;
      tree_pi_vec[j].graph_split = graph_split.begin();
    }
    tree_pi_vec[j].mu0 = mu0(j);
    tree_pi_vec[j].tau = tau(j);
  }
  
  // sigma stuff
  double sigma = 1.0;
  
  // stuff for MCMC loop
  int total_draws = 1 + burn + (nd - 1) * thin;
  int i = 0;
  int sample_index = 0;
  int accept = 0;
  int* total_accept = new int[p];
  for(int j = 0; j < p; j++) total_accept[j] = 0;
  double tmp_mu; // for holding the value of mu that we're updating
  
  // initialize the trees
  std::vector<std::vector<tree>> tree_vec(p, std::vector<tree>(M));
  std::vector<std::vector<suff_stat>> ss_train_vec(p, std::vector<suff_stat>(M));
  std::vector<std::vector<suff_stat>> ss_test_vec(p, std::vector<suff_stat>(M));
  for(int j = 0; j < p; j++){
    for(int m = 0; m < M; m++){
      tree_traversal(ss_train_vec[j][m], tree_vec[j][m], di_train);
      
      for(suff_stat_it l_it = ss_train_vec[j][m].begin(); l_it != ss_train_vec[j][m].end(); ++l_it){
        tmp_mu = tree_vec[j][m].get_ptr(l_it->first)->get_mu();
        for(int_it it = l_it->second.begin(); it != l_it->second.end(); ++it){
          i = *it;
          allfit_train[i] += tX_train(j,i) * tmp_mu;
          beta_fit_train[j + i*p] += tmp_mu;
        }
      }
      
      if(n_test > 0){
        tree_traversal(ss_test_vec[j][m], tree_vec[j][m], di_test);
        for(suff_stat_it l_it = ss_test_vec[j][m].begin(); l_it != ss_test_vec[j][m].end(); ++l_it){
          tmp_mu = tree_vec[j][m].get_ptr(l_it->first)->get_mu();
          for(int_it it = l_it->second.begin(); it != l_it->second.end(); ++it){
            i = *it;
            allfit_test[i] += tX_test(j,i) * tmp_mu;
            beta_fit_test[j+i*p] += tmp_mu;
          }
        }
      }
    }
  }

  for(int i = 0; i < N_train; i++){
    residual[i] = Y_train[i] - allfit_train[i];
    r_sum[subj_id_train[i]] += residual[i];
    r2_sum[subj_id_train[i]] += pow(residual[i], 2.0);
  }
  
  // output containers
  arma::vec fit_train_mean = arma::zeros<arma::vec>(N_train);
  arma::vec fit_test_mean = arma::zeros<arma::vec>(tmp_N_test);
  
  arma::mat beta_train_mean = arma::zeros<arma::mat>(N_train,p);
  arma::mat beta_test_mean = arma::zeros<arma::mat>(tmp_N_test,p);


  arma::mat fit_train = arma::zeros<arma::mat>(1,1); // posterior samples for training data
  arma::mat fit_test = arma::zeros<arma::mat>(1,1); // posterior samples for testing data
  
  arma::cube beta_train = arma::zeros<arma::cube>(1,1,1);
  arma::cube beta_test = arma::zeros<arma::cube>(1,1,1);
  
  if(save_samples){
    fit_train.zeros(nd,N_train);
    beta_train.zeros(nd, N_train,p);
    if(n_test > 0){
      fit_test.zeros(nd,N_test);
      beta_test.zeros(nd, N_test,p);
    }
  }
  
  arma::vec sigma_samples(total_draws);
  arma::mat total_accept_samples(nd,p);
  arma::cube theta_samples(1,1,1);
  if(sparse) theta_samples.zeros(total_draws, R, p);
  arma::cube var_count_samples(total_draws, R, p);
  
  Rcpp::List tree_draws(nd);
  

  // main MCMC loop starts here!
  for(int iter = 0; iter < total_draws; iter++){
    if(verbose){
      if( (iter < burn) && (iter % print_every == 0)){
        Rcpp::Rcout << "  MCMC Iteration: " << iter << " of " << total_draws << "; Warmup" << std::endl;
        Rcpp::checkUserInterrupt();
      } else if(( (iter > burn) && (iter%print_every == 0)) || (iter == burn)){
        Rcpp::Rcout << "  MCMC Iteration: " << iter << " of " << total_draws << "; Sampling" << std::endl;
      }
    }
    
    if( (iter == burn) && (n_test > 0)){
      // in first sampling iteration, we should populate allfit_test and beta_fit_test
      for(int j = 0; j < p; j++){
        for(int m = 0; m < M; m++){
          for(suff_stat_it l_it = ss_test_vec[j][m].begin(); l_it != ss_test_vec[j][m].end(); ++l_it){
            tmp_mu = tree_vec[j][m].get_ptr(l_it->first)->get_mu(); // get value of mu in the leaf
            for(int_it it = l_it->second.begin(); it != l_it->second.end(); ++it){
              i = *it;
              allfit_test[i] += tX_test(j,i) * tmp_mu;
              beta_fit_test[j + i*p] += tmp_mu;
            }
          }
        }
      }
    } // finished computing allfit_test and beta_fit_test in first post-warmup sample
    
    // loop over ensembles
    for(int j = 0; j < p; j++){
      total_accept[j] = 0;
      for(int m = 0; m < M; m++){
        // remove fit of the m-th tree in the j-th ensemble.
        // we do this by looping over the leafs of the tree
        for(suff_stat_it l_it = ss_train_vec[j][m].begin(); l_it != ss_train_vec[j][m].end(); ++l_it){
          tmp_mu = tree_vec[j][m].get_ptr(l_it->first)->get_mu(); // get the value of mu in the leaf
          for(int_it it = l_it->second.begin(); it != l_it->second.end(); ++it){
            // remove fit of m-th tree from allfit: allfit[*it] -= X(*it,j) * tmp_mu
            // for partial residual: we need to compute Y - allfit[*it]
            // numerically this exactly equal to ADDING X(*it,j) * tmp_mu to residual[*it]
            i = *it;
            allfit_train[i] -= tX_train(j,i) * tmp_mu;
            residual[i] += tX_train(j,i) * tmp_mu;
            beta_fit_train[j + i*p] -= tmp_mu; // remove fit of m-th tree from overall fit of beta's
            r_sum[subj_id_train[i]] -= tX_train(j,i) * tmp_mu; // remove contribution to total sum of residuals for subject contributing i-th obs
            r2_sum[subj_id_train[i]] += 2.0 * tX_train(j,i) * tmp_mu * residual[i] - pow(tX_train(j,i) * tmp_mu, 2);
          } // close loop over obs in single leaf
        } // close loop over all leafs
        
        // if we have testing data & we're post-warmup, then we should remove the value of m-th tree of j-th ensemble from allfit_test and beta_fit_test
        if(iter >= burn && n_test > 0){
          for(suff_stat_it l_it = ss_test_vec[j][m].begin(); l_it != ss_test_vec[j][m].end(); ++l_it){
            tmp_mu = tree_vec[j][m].get_ptr(l_it->first)->get_mu();
            for(int_it it = l_it->second.begin(); it != l_it->second.end(); ++it){
              i = *it;
              allfit_test[i] -= tX_test(j,i) * tmp_mu;
              beta_fit_test[j + i*p] -= tmp_mu;
            }
          }
        } // close if checking if iter >= burn and n_test > 0
        
        
        // now we do the tree update
        update_tree_ind(tree_vec[j][m],accept, j,sigma, ss_train_vec[j][m], ss_test_vec[j][m], di_train, di_test, tree_pi_vec[j], gen);
        total_accept[j] += accept;
        
        // now we add the fit of the m-th tree in the j-th ensemble back to the relevant stuff
        for(suff_stat_it l_it = ss_train_vec[j][m].begin(); l_it != ss_train_vec[j][m].end(); ++l_it){
          tmp_mu = tree_vec[j][m].get_ptr(l_it->first)->get_mu(); // get the value of mu in the leaf
          for(int_it it = l_it->second.begin(); it != l_it->second.end(); ++it){
            // add fit of m-th tree from allfit: allfit[*it] += X(*it,j) * tmp_mu
            // for partial residual: we need to compute Y - allfit[*it]
            // numerically this exactly equal to SUBTRACTING X(*it,j) * tmp_mu to residual[*it]
            i = *it;
            r2_sum[subj_id_train[i]] += pow(tX_train(j,i) * tmp_mu, 2) - 2.0 * tX_train(j,i) * tmp_mu * residual[i]; // remember, residual[i] currently holds the PARTIAL residual.
            r_sum[subj_id_train[i]] += tX_train(j,i) * tmp_mu; // restore contribution to total sum of residuals for subject contributing i-th obs
            allfit_train[i] += tX_train(j,i) * tmp_mu;
            residual[i] -= tX_train(j,i) * tmp_mu;
            beta_fit_train[j + i*p] += tmp_mu; // restore fit of m-th tree from overall fit of beta's
            //r2_sum[subj_id_train[i]] += pow(tX_train(j,i) * tmp_mu,2); // remove contribution to total sum of squared residuals for subject contribution i-th obs WRONG!!!
          } // close loop over obs in single leaf
        } // close loop over all leafs
        
        // if we have testing data & we're post-warmup, then we should add back in the value of m-th tree of j-th ensemble to allfit_test and beta_fit_test
        if(iter >= burn && n_test > 0){
          for(suff_stat_it l_it = ss_test_vec[j][m].begin(); l_it != ss_test_vec[j][m].end(); ++l_it){
            tmp_mu = tree_vec[j][m].get_ptr(l_it->first)->get_mu();
            for(int_it it = l_it->second.begin(); it != l_it->second.end(); ++it){
              i = *it;
              allfit_test[i] += tX_test(j,i) * tmp_mu;
              beta_fit_test[j + i*p] += tmp_mu;
            }
          }
        }
      } // closes loop over trees in j-th ensemble
      
      // now we can update theta
      if(sparse){
        update_theta_u(theta[j], u[j], var_count[j], R, a_u, b_u, gen);
        for(int r = 0; r < R; r++){
          theta_samples(iter,r,j) = theta[j][r];
          var_count_samples(iter,r,j) = var_count[j][r];
        }
      }
    } // closes loop over j-th ensemble
    /*
     This was for checking that we were computing residuals correctly
    // we need to update the sum of squared residuals for each individual
    for(int subj_ix = 0; subj_ix < n_train; subj_ix++){
      tmp_r_sum[subj_ix] = 0.0;
      tmp_r2_sum[subj_ix] = 0.0;
    }

    for(int i = 0; i < N_train; i++){
      //if(residual[i] != Y_train[i] - allfit_train[i]){
      if(abs(residual[i] - (Y_train[i] - allfit_train[i])) > 1e-12){
        Rcpp::Rcout << " iter = " << iter << " i = " << i << " residual = " << residual[i] << std::endl;
        Rcpp::Rcout << "  Y - allfit = " << Y_train[i] - allfit_train[i] << std::endl;
        Rcpp::Rcout << "diff = " << abs(residual[i] - (Y_train[i] - allfit_train[i])) << std::endl;
        Rcpp::stop("something is wrong with the residuals!");
      }
      tmp_r_sum[subj_id_train[i]] += residual[i];
      tmp_r2_sum[subj_id_train[i]] += pow(residual[i], 2.0);
    }
    
    for(int subj_ix = 0; subj_ix < n_train; subj_ix++){
      if(abs(tmp_r_sum[subj_ix] - r_sum[subj_ix]) > 1e-12){
        Rcpp::Rcout << "  iter = " << iter << " subj_ix = " << subj_ix << std::endl;
        Rcpp::Rcout << "  r_sum = " << r_sum[subj_ix] << " tmp_r_sum = " << tmp_r_sum[subj_ix] << std::endl;
        Rcpp::Rcout << "  diff = " << abs(tmp_r_sum[subj_ix] - r_sum[subj_ix]) << std::endl;
        Rcpp::stop("something is wrong with updating sum of residuals per subject");
      }
      if(abs(tmp_r2_sum[subj_ix] - r2_sum[subj_ix]) > 1e-12){
        Rcpp::Rcout << "  iter = " << iter << " subj_ix = " << subj_ix << std::endl;
        Rcpp::Rcout << "  r2_sum = " << r2_sum[subj_ix] << " tmp_r2_sum = " << tmp_r2_sum[subj_ix] << std::endl;
        Rcpp::Rcout << "  diff = " << abs(tmp_r2_sum[subj_ix] - r2_sum[subj_ix]) << std::endl;
        Rcpp::stop("something is wrong with updating sum of squared residuals per subject");
      }
      
    }
    */
    update_sigma_ind(sigma, nu, lambda, di_train, gen);
    sigma_samples(iter) = sigma;
    
    if( (iter >= burn) && ( (iter - burn) % thin == 0)){
      sample_index = (int) ( (iter-burn)/thin);
      for(int j = 0; j < p; j++) total_accept_samples(sample_index,j) = total_accept[j];
      
      if(save_trees){
        Rcpp::List tmp_tree_draws(p);
        for(int j = 0; j < p; j++){
          Rcpp::CharacterVector tree_string_vec(M);
          for(int m = 0; m < M; m++){
            tree_string_vec[m] = write_tree(tree_vec[j][m], tree_pi_vec[j], set_str);
          }
          tmp_tree_draws[j] = tree_string_vec;
        }
        tree_draws[sample_index] = tmp_tree_draws;
      }
      
      if(save_samples){
        for(int i = 0; i < N_train; i++){
          fit_train(sample_index,i) = allfit_train[i];
          fit_train_mean(i) += allfit_train[i];
          for(int j = 0; j < p; j++){
            beta_train(sample_index,i,j) = beta_fit_train[j + i * p];
            beta_train_mean(i,j) += beta_fit_train[j + i*p];
          }
        }
        
        if(n_test > 0){
          for(int i = 0; i < N_test; i++){
            fit_test(sample_index,i) = allfit_test[i];
            fit_test_mean(i) += allfit_test[i];
            for(int j = 0; j < p; j++){
              beta_test(sample_index,i,j) = beta_fit_test[j + i*p];
              beta_test_mean(i,j) += beta_fit_test[j+i*p];
            }
          }
        }
      } else{
        for(int i = 0; i < N_train; i++){
          fit_train_mean(i) += allfit_train[i];
          for(int j = 0; j < p; j++){
            beta_train_mean(i,j) += beta_fit_train[j + i*p];
          }
        }
        if(n_test > 0){
          for(int i = 0; i < N_test; i++){
            fit_test_mean(i) += allfit_test[i];
            for(int j = 0; j < p; j++){
              beta_test_mean(i,j) += beta_fit_test[j + i*p];
            }
          }
        }
      }
    } // closes if checking if we should be saving samples or not
  } // closes main MCMC loop
  
  fit_train_mean /= ( (double) nd);
  beta_train_mean /= ( (double) nd);
  if(n_test > 0){
    fit_test_mean /= ( (double) nd);
    beta_test_mean /= ( (double) nd);
  }
  
  Rcpp::List results;
  results["fit_train_mean"] = fit_train_mean;
  results["beta_train_mean"] = beta_train_mean;
  if(save_samples){
    results["fit_train"] = fit_train;
    results["beta_train"] = beta_train;
  }
  if(n_test > 0){
    results["fit_test_mean"] = fit_test_mean;
    results["beta_test_mean"] = beta_test_mean;
    if(save_samples){
      results["fit_test"] = fit_test;
      results["beta_test"] = beta_test;
    }
  }
  results["sigma"] = sigma_samples;
  results["total_accept"] = total_accept_samples;
  results["var_count"] = var_count_samples;
  if(save_trees) results["trees"] = tree_draws;
  if(sparse) results["theta"] = theta_samples;

  
  delete[] allfit_train;
  delete[] allfit_test;
  delete[] beta_fit_train;
  delete[] beta_fit_test;
  delete[] residual;
  delete[] r_sum;
  delete[] r2_sum;
  
  return results;
}
