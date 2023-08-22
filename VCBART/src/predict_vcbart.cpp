#include "funs.h"

// [[Rcpp::export(".predict_vcbart")]]
arma::cube predict_vcbart(Rcpp::List tree_draws,
                          int p, // will include intercept
                          int M,
                          Rcpp::NumericMatrix tZ_cont,
                          Rcpp::IntegerMatrix tZ_cat,
                          bool verbose = true)
{
  set_str_conversion set_str; // for converting sets of integers to strings
  
  int N;
  int R_cont;
  int R_cat;
  
  if(tZ_cont.size() > 1 && tZ_cat.size() == 1){
    // only continuous modifiers supplied
    N = tZ_cont.cols();
    R_cont = tZ_cont.rows();
    R_cat = 0;
  } else if(tZ_cont.size() == 1 && tZ_cat.size() > 1){
    // only categorical modifiers supplied
    N = tZ_cat.cols();
    R_cont = 0;
    R_cat = tZ_cat.rows();
  } else if(tZ_cont.size() > 1 && tZ_cat.size() > 1){
    if(tZ_cont.cols() != tZ_cat.cols()){
      Rcpp::Rcout << "Z_cont has " << tZ_cont.cols() << " rows" << std::endl;
      Rcpp::Rcout << "Z_cat has " << tZ_cat.cols() << " rows" << std::endl;
      Rcpp::stop("[predict_vcbart]: Z_cont and Z_cat must have same number of rows!");
    } else{
      N = tZ_cont.cols();
      R_cont = tZ_cont.rows();
      R_cat = tZ_cat.rows();
    }
  } else{
    Rcpp::stop("No modifiers supplied. Nothing to predict!");
  }
  
  int R = R_cont + R_cat;
  
  data_info di;
  di.N = N;
  di.p = p;
  di.R_cont = R_cont;
  di.R_cat = R_cat;
  di.R = R;
  
  if(R_cont > 0) di.z_cont = tZ_cont.begin();
  if(R_cat > 0) di.z_cat = tZ_cat.begin();
  
  int nd = tree_draws.size();

  //Rcpp::CharacterVector first_tree_vec = tree_draws[0];
  //int M = first_tree_vec.size(); // how many trees in the ensemble
  
  if(verbose){
    Rcpp::Rcout << "nd = " << nd << "M = " << M << std::endl;
  }
  
  std::vector<double> allfit(N);

  arma::cube pred_out = arma::zeros<arma::cube>(nd, N, p);
  

  int print_every = floor(nd/10);
  
  // tree_draws is a list of list of character vectors
  // outer list indexed by iterations
  // inner list indexed by covariates
  
  for(int iter = 0; iter < nd; ++iter){
    if(verbose && (iter%print_every == 0)){
      Rcpp::Rcout << " Iteration: " << iter << " of " << nd << std::endl;
      Rcpp::checkUserInterrupt();
    }
    
    //Rcpp::CharacterVector tmp_string_vec = tree_draws[iter];
    Rcpp::List tmp_draw = tree_draws[iter];
    if(tmp_draw.size() != di.p){
      Rcpp::Rcout << " iter = " << iter << " found " << tmp_draw.size() << " ensembles." << std::endl;
      Rcpp::Rcout << " Expected " << di.p << " ensembles." << std::endl;
      Rcpp::stop("Unexpected number of tree ensembles detected.");
    } else{
      for(int j = 0; j < di.p; ++j){
        Rcpp::CharacterVector tmp_string_vec = tmp_draw[j];
        if(tmp_string_vec.size() != M){
          // did we not record enough tree strings??
          // we should never hit this part of the code
          Rcpp::Rcout << " iter = " << iter << " j = " << j;
          Rcpp::Rcout << " # tree strings = " << tmp_string_vec.size() << std::endl;
          Rcpp::stop("Unexpeted number of tree strings!");
        } else{
          std::vector<tree> t_vec(M);
          for(int m = 0; m < M; ++m){
            // tmp_string is an Rcpp::CharacterVector
            // we will extract a single element and turn it into a std::string
            // that will be passed to read_tree
            std::string tmp_string = Rcpp::as<std::string>(tmp_string_vec[m]);
            read_tree(t_vec[m], tmp_string, set_str);
          }// closes loop over trees in ensemble
          
          fit_ensemble(allfit, t_vec, di);
          for(int i = 0; i < N; ++i) pred_out(iter, i, j) = allfit[i];
        } // closes if/else checking that we had enough trees in this draw of the beta-j ensemble
      } // closes loop over all beta ensembles
    } // closes if/else checking that we saved the right number of tree ensembles in this iteration
  } // closes loop over iterations

  return pred_out;
  
}
