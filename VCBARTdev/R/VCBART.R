VCBART <- function(Y_train, 
                   X_train, Z_train, n_train, 
                   X_test, Z_test, n_test,
                   cutpoints, 
                   intercept = TRUE,
                   M = 50,
                   error_structure = c("ind", "cs"),
                   split_probs_type = c("adaptive", "fixed"),
                   ht_sigma_y = TRUE,
                   ht_tau = TRUE,
                   burn = 500, nd = 1000,
                   verbose = TRUE, print_every = 50,
                   rho_eps = NULL,
                   split_probs = NULL,
                   a = NULL, 
                   b = NULL,
                   N_u = NULL,
                   rho_alpha = NULL, 
                   tau_vec = NULL, 
                   alpha_vec = NULL, 
                   beta_vec = NULL)
{
  
  # Some sanity checks: dimensions and missing values
  if(!is.matrix(X_train)) stop("X_train must be a matrix!")
  if(!is.matrix(Z_train)) stop("Z_train must be a matrix!")
  if(nrow(X_train) != nrow(Z_train)) stop("X_train and Z_train must have same number of rows!")
  if(nrow(X_train) != length(Y_train)) stop("nrow(X_train) must equal length(Y_train) must have same number")
  
  if(!is.matrix(X_test)) stop("X_test must be a matrix!")
  if(!is.matrix(Z_test)) stop("Z_test must be a matrix!")
  if(nrow(X_test) != nrow(Z_test)) stop("X_test and Z_test must have same number of rows!")
  
  if(ncol(X_train) != ncol(X_test)) stop("X_train and X_test must have same number of columns!")
  if(ncol(Z_train) != ncol(Z_test)) stop("Z_train and Z_test must have same number of columns!")
  
  N_train <- length(n_train) # number of subjects in training set
  N_test <- length(n_test) # number of subjects in test set 
  
  start_index_train <- 1 + c(0, cumsum(n_train)[-N_train])
  start_index_test <- 1 + c(0, cumsum(n_test)[-N_test])
  
  
  n_obs_train <- nrow(X_train)
  n_obs_test <- nrow(X_test)

  if(sum(n_train) != n_obs_train) stop("Total number of training observations is not equal to number of rows of X_train")
  if(sum(n_test) != n_obs_test) stop("Total number of testing observations is not equal to number of rows of X_train")
  
  if(intercept == TRUE){
    X_train <- cbind(rep(1, times = n_obs_train), X_train)
    X_test <- cbind(rep(1, times = n_obs_test), X_test)
  }
  
  p <- ncol(X_train)
  R <- ncol(Z_train)
  
  if(!is.null(tau_vec)){
    if(length(tau_vec) != p){
      if(intercept == TRUE) stop("length(tau_vec) must have length 1 + ncol(X_train)!")
      else if(intercept == FALSE) stop("length(tau_vec) must have length ncol(X_train!")
    }
  } else{
    tau_vec <- rep(1/sqrt(M), times = p)
  }
  
  if(!is.null(alpha_vec)){
    if(length(alpha_vec) != p){
      if(intercept == TRUE) stop("alpha_vec must have length 1 + ncol(X_train)!")
      else if(intercept == FALSE) stop("alpha_vec must have length ncol(X_train!")
    }
  } else{
    alpha_vec <- rep(0.95, times = p)
  }
  
  if(!is.null(beta_vec)){
    if(length(beta_vec) != p){
      if(intercept == TRUE) stop("beta_vec must have length 1 + ncol(X_train)!")
      else if(intercept == FALSE) stop("beta_vec must have length ncol(X_train!")
    }
  } else{
    beta_vec <- rep(0.95, times = p)
  }

  # eventually handle these arguments in a better way
  sigma_hat <- 1
  nu_sigma <- 7
  nu_tau <- 7
  variance_prob <- 0.9
  ht_tau <- FALSE
  
  if(split_probs_type == "fixed"){
    if(is.null(split_probs)){
      print("Fixed split probability not specified. Resorting to default of uniform split probabilities")
      split_probs <- list()
      for(k in 1:p) split_probs[[k]] <- rep(1/R, times = R)
    } else{
      if(!is.list(split_probs)) stop("If passed, split_probs must be a list")
      else if(length(split_probs) != p & intercept == FALSE) stop("length(split_probs) must equal ncol(X_train)")
      else if(length(split_probs) != p & intercept == TRUE) stop("length(split_probs) must equal 1 + ncol(X_train)")
      else{
        # if we get to this statement, user has given a list of the right length
        split_probs_lengths <- sapply(split_probs, FUN = length)
        split_probs_pos <- sapply(split_probs, FUN = function(x){return(all(x >= 0))})
        split_probs_sum <- sapply(split_probs, FUN = sum)
        
        if(any(split_probs_lengths != R) | any(!split_probs_pos) | any(split_probs_sum != 1)){
          stop("Each element of split_probs must vector of length ncol(Z_train) of non-negative probabilities summing to 1!")
        }
      }
    }
    
    # If we have reached here without stopping, we will be ready to run VC-BART with fixed splitting probabilities
    if(error_structure == "ind"){
      print("Entering VCBART w/ independent errors & fixed split probabilities")
      fit <- .vcbart_ind_fixed_split(Y_train, 
                                     X_train, Z_train, n_train, start_index_train,
                                     X_test, Z_test, n_test, start_index_test,
                                     cutpoints, M, ht_sigma_y, ht_tau,
                                     burn, nd, verbose, print_every,
                                     split_probs, tau_vec, alpha_vec, beta_vec,
                                     sigma_hat, nu_sigma, nu_tau, variance_prob)
    } else if(error_structure == "cs" & is.null(rho_eps)){
      print("Entering VCBART w/ compound symmetry errors, adaptive autocorrelation, and fixed split probabilities")
      rho_eps <- 0.5 # initial value of rho_eps
      fit <- .vcbart_cs_fixed_split_adapt_rho(Y_train,
                                              X_train, Z_train, n_train, start_index_train,
                                              X_test, Z_test, n_test, start_index_test,
                                              cutpoints, M, ht_sigma_y, ht_tau,
                                              burn, nd, verbose, print_every, 
                                              rho_eps, split_probs, tau_vec, alpha_vec, beta_vec,
                                              sigma_hat, nu_sigma, nu_tau, variance_prob)
    } else if(error_structure == "cs" & !is.null(rho_eps)){
      print("Entering VCBART w/ compound symmetry errors, fixed autocorrelation, and fixed split probabilities")
      if(abs(rho_eps - 0.5) >= 0.5) stop("rho_eps must be in the open interval (0,1)!")
      fit <- .vcbart_cs_fixed_split_fixed_rho(Y_train,
                                              X_train, Z_train, n_train, start_index_train,
                                              X_test, Z_test, n_test, start_index_test,
                                              cutpoints, M, ht_sigma_y, ht_tau,
                                              burn, nd, verbose, print_every,
                                              rho_eps, split_probs, tau_vec, alpha_vec, beta_vec,
                                              sigma_hat, nu_sigma, nu_tau, variance_prob)
    }
  } else if(split_probs_type == "adaptive"){
    
    if(is.null(a)) a <- 1
    if(is.null(b)) b <- R
    if(is.null(rho_alpha)) rho_alpha <- R
    if(is.null(N_u)) N_u <- 100
    
    if(error_structure == "ind"){
      print("Entering VCBART w/ independent errors and adaptive split probabilities")
      fit <- .vcbart_ind_adapt_split(Y_train, 
                                     X_train, Z_train, n_train, start_index_train,
                                     X_test, Z_test, n_test, start_index_test,
                                     cutpoints, M, ht_sigma_y, ht_tau,
                                     burn, nd, verbose, print_every,
                                     a, b, N_u, rho_alpha,
                                     tau_vec, alpha_vec, beta_vec, 
                                     sigma_hat, nu_sigma, nu_tau, variance_prob)
    } else if(error_structure == "cs" & is.null(rho_eps)){
      print("Entering VCBART w/ compound symmetry errors, adaptive autocorrelation, and adaptive split probabilities")
      rho_eps <- 0.5 # initialize rho_eps
      fit <- .vcbart_cs_adapt_split_fixed_rho(Y_train, 
                                              X_train, Z_train, n_train, start_index_train,
                                              X_test, Z_test, n_test, start_index_test,
                                              cutpoints, M, ht_sigma_y, ht_tau,
                                              burn, nd, verbose, print_every,
                                              rho_eps,
                                              a, b, N_u, rho_alpha,
                                              tau_vec, alpha_vec, beta_vec, 
                                              sigma_hat, nu_sigma, nu_tau, variance_prob)
    } else if(error_structure == "cs" & !is.null(rho_eps)){
      print("Entering VCBART w/ compound symmetry errors, fixed autocorrelation, and adaptive split probabilities")
      if(abs(rho_eps - 0.5) >= 0.5) stop("rho_eps must be in the open interval (0,1)!")
      fit <- .vcbart_cs_adapt_split_fixed_rho(Y_train, 
                                              X_train, Z_train, n_train, start_index_train,
                                              X_test, Z_test, n_test, start_index_test,
                                              cutpoints, M, ht_sigma_y, ht_tau,
                                              burn, nd, verbose, print_every,
                                              rho_eps,
                                              a, b, N_u, rho_alpha,
                                              tau_vec, alpha_vec, beta_vec, 
                                              sigma_hat, nu_sigma, nu_tau, variance_prob)
    }
  } else stop("split_probs_type must be \"adaptive\" or \"fixed\" ")
  return(fit)
}
