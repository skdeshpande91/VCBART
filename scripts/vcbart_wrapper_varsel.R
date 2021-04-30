# wrapper for VCBART used in the simulation studies about variable selection

#library(VCBARTdev)
library(VCBART)
vcbart_varsel_wrapper <- function(Y_train, 
                           X_train, Z_train, n_train,
                           X_test, Z_test, n_test,
                           cutpoints,
                           error_structure,
                           split_probs_type,
                           tau = 1/sqrt(50),
                           burn, nd, verbose, print_every)
{
  p <- ncol(X_train)
  R <- ncol(Z_train)
  tau_vec = rep(tau, times = p+1) # we're always going to include an intercept
  run_time <- system.time({
    chain1 <- VCBART(Y_train = Y_train,
                     X_train = X_train, Z_train = Z_train, n_train = n_train,
                     X_test = X_test, Z_test = Z_test, n_test = n_test, 
                     cutpoints = cutpoints, 
                     error_structure = error_structure,
                     split_probs_type = split_probs_type, tau_vec = tau_vec,
                     burn = burn, nd = nd, verbose = verbose, print_every = print_every)
    chain2 <- VCBART(Y_train = Y_train,
                     X_train = X_train, Z_train = Z_train, n_train = n_train,
                     X_test = X_test, Z_test = Z_test, n_test = n_test, 
                     cutpoints = cutpoints, 
                     error_structure = error_structure,
                     split_probs_type = split_probs_type, tau_vec = tau_vec,
                     burn = burn, nd = nd, verbose = verbose, print_every = print_every)
  })["elapsed"]
  
  
  beta_summary <- summarize_beta(chain1, chain2, burn = burn)
  ystar_summary <- summarize_posterior_predictive(chain1, chain2, burn = burn)
  
  varcount_mean <- 
    0.5 * apply(chain1$var_counts_samples[,,-(1:burn)], MARGIN = c(1,2), FUN = mean) +
    0.5 * apply(chain2$var_counts_samples[,,-(1:burn)], MARGIN = c(1,2), FUN = mean)

  tmp_dim <- dim(chain1$var_counts_samples[,,-(1:burn)])
  tmp_dim[3] <- 2 * tmp_dim[3]
  if(tmp_dim[3] != 2 * nd) stop("Something's wrong with dimenions of var_count_samples")
  
  varcount_samples <- array(dim = tmp_dim)
  varcount_samples[,,1:nd] <- chain1$var_counts_samples[,,-(1:burn)]
  varcount_samples[,,(nd+1):(2*nd)] <- chain2$var_counts_samples[,,-(1:burn)]
  varcount_prob <- apply(varcount_samples >= 1, FUN = mean, MARGIN = c(1,2))
  
  # varcount_mean[j,r]: posterior mean # times Z_r is split on in ensemble appxing beta_j
  # varcount_prob[j,r]: posterior prob Z_r is split on at least once in ensemble appxing beta_j
  
  if(split_probs_type == "adaptive"){
    theta_mean <- 
      0.5 * apply(chain1$theta_samples[,,-(1:burn)], MARGIN = c(1,2), FUN = mean) +
      0.5 * apply(chain2$theta_samples[,,-(1:burn)], MARGIN = c(1,2), FUN = mean)
  } else{
    theta_mean <- matrix(1/R, nrow = R, ncol = p+1)
  }
  
  results <- list(train = list(beta = beta_summary$train, ystar = ystar_summary$train),
                  test = list(beta = beta_summary$test, ystar = ystar_summary$test),
                  varcount_mean = varcount_mean,
                  varcount_prob = varcount_prob,
                  varcount_samples = varcount_samples,
                  theta_mean = theta_mean,
                  time = run_time,
                  vcbart_time = chain1$time + chain2$time)
  
  
  return(results)
  
}

dirichlet_wrapper <- function(Y_train, 
                              X_train, Z_train, n_train,
                              X_test, Z_test, n_test,
                              cutpoints,
                              init_alpha_z,
                              tau = 1/sqrt(50),
                              burn, nd, verbose, print_every)
{
  p <- ncol(X_train)
  R <- ncol(Z_train)
  tau_vec = rep(tau, times = p+1) # we're always going to include an intercept
  run_time <- system.time({
    chain1 <- VCBART(Y_train = Y_train,
                     X_train = X_train, Z_train = Z_train, n_train = n_train,
                     X_test = X_test, Z_test = Z_test, n_test = n_test, 
                     cutpoints = cutpoints, 
                     error_structure = "ind",
                     split_probs_type = "adaptive", 
                     fixed_dirichlet = TRUE, init_alpha_z = init_alpha_z,
                     tau_vec = tau_vec,
                     burn = burn, nd = nd, verbose = verbose, print_every = print_every)
    chain2 <- VCBART(Y_train = Y_train,
                     X_train = X_train, Z_train = Z_train, n_train = n_train,
                     X_test = X_test, Z_test = Z_test, n_test = n_test, 
                     cutpoints = cutpoints, 
                     error_structure = "ind",
                     split_probs_type = "adaptive", 
                     fixed_dirichlet = TRUE, init_alpha_z = init_alpha_z,
                     tau_vec = tau_vec,
                     burn = burn, nd = nd, verbose = verbose, print_every = print_every)
  })["elapsed"]
  
  beta_summary <- summarize_beta(chain1, chain2, burn = burn)
  ystar_summary <- summarize_posterior_predictive(chain1, chain2, burn = burn)
  
  varcount_mean <- 
    0.5 * apply(chain1$var_counts_samples[,,-(1:burn)], MARGIN = c(1,2), FUN = mean) +
    0.5 * apply(chain2$var_counts_samples[,,-(1:burn)], MARGIN = c(1,2), FUN = mean)
  
  tmp_dim <- dim(chain1$var_counts_samples[,,-(1:burn)])
  tmp_dim[3] <- 2 * tmp_dim[3]
  if(tmp_dim[3] != 2 * nd) stop("Something's wrong with dimenions of var_count_samples")
  
  varcount_samples <- array(dim = tmp_dim)
  varcount_samples[,,1:nd] <- chain1$var_counts_samples[,,-(1:burn)]
  varcount_samples[,,(nd+1):(2*nd)] <- chain2$var_counts_samples[,,-(1:burn)]
  varcount_prob <- apply(varcount_samples >= 1, FUN = mean, MARGIN = c(1,2))
  
  theta_mean <- 
    0.5 * apply(chain1$theta_samples[,,-(1:burn)], MARGIN = c(1,2), FUN = mean) +
    0.5 * apply(chain2$theta_samples[,,-(1:burn)], MARGIN = c(1,2), FUN = mean)
  
  results <- list(train = list(beta = beta_summary$train, ystar = ystar_summary$train),
                  test = list(beta = beta_summary$test, ystar = ystar_summary$test),
                  varcount_mean = varcount_mean,
                  varcount_prob = varcount_prob,
                  varcount_samples = varcount_samples,
                  theta_mean = theta_mean,
                  time = run_time,
                  vcbart_time = chain1$time + chain2$time)
  
  
  return(results)
}


