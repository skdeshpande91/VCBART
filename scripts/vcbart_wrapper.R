# wrapper for VCBART
# Used for the main simulation studies
#library(VCBARTdev)
library(VCBART)
vcbart_wrapper <- function(Y_train, 
                           X_train, Z_train, n_train,
                           X_test, Z_test, n_test,
                           cutpoints,
                           error_structure,
                           split_probs_type,
                           M = 50,
                           tau = 1/sqrt(M),
                           burn, nd, verbose, print_every)
{
  p <- ncol(X_train)
  tau_vec = rep(tau, times = p+1) # we're always going to include an intercept
  run_time <- system.time({
    chain1 <- VCBART(Y_train = Y_train,
                     X_train = X_train, Z_train = Z_train, n_train = n_train,
                     X_test = X_test, Z_test = Z_test, n_test = n_test, 
                     cutpoints = cutpoints, M = M,
                     error_structure = error_structure,
                     split_probs_type = split_probs_type, tau_vec = tau_vec,
                     burn = burn, nd = nd, verbose = verbose, print_every = print_every)
    chain2 <- VCBART(Y_train = Y_train,
                     X_train = X_train, Z_train = Z_train, n_train = n_train,
                     X_test = X_test, Z_test = Z_test, n_test = n_test, 
                     cutpoints = cutpoints, M = M,
                     error_structure = error_structure,
                     split_probs_type = split_probs_type, tau_vec = tau_vec,
                     burn = burn, nd = nd, verbose = verbose, print_every = print_every)
  })["elapsed"]
  
  
  beta_summary <- summarize_beta(chain1, chain2, burn = burn)
  ystar_summary <- summarize_posterior_predictive(chain1, chain2, burn = burn)
  beta_support <- get_beta_support(chain1, chain2, burn = burn, max_cutoff = 10)
  
  results <- list(train = list(beta = beta_summary$train, ystar = ystar_summary$train),
                  test = list(beta = beta_summary$test, ystar = ystar_summary$test),
                  beta_support = beta_support,
                  time = run_time,
                  vcbart_time = chain1$time + chain2$time)
  
  return(results)
  
}