vcbart_pred_wrapper <- function(Y_train,
                           subj_id_train,
                           ni_train,
                           X_train, 
                           Z_cont_train = matrix(0, nrow = 1, ncol = 1),
                           Z_cat_train = matrix(0, nrow = 1, ncol = 1),
                           X_test = matrix(0, nrow = 1, ncol = 1),
                           Z_cont_test = matrix(0, nrow = 1, ncol = 1),
                           Z_cat_test = matrix(0, nrow = 1, ncol = 1),
                           unif_cuts = rep(TRUE, times = ncol(Z_cont_train)),
                           cutpoints_list = NULL,
                           cat_levels_list = NULL,
                           sparse = TRUE,
                           M = 200, 
                           mu0 = NULL, tau = NULL, nu = NULL, lambda = NULL,
                           nd = 1000, burn = 1000, thin = 1,
                           save_samples = FALSE, save_trees = FALSE,
                           verbose = TRUE, print_every = floor( (nd * thin + burn)/10),
                           n_chains = 4)
{
  if(is.null(tau)) tau <- rep(0.5/sqrt(M), times = p+1)
  p <- ncol(X_train)
  
  if(length(Z_cont_train) != 1) R_cont <- ncol(Z_cont_train)
  else R_cont <- 0
  
  if(length(Z_cat_train) != 1) R_cat <- ncol(Z_cat_train)
  else R_cat <- 0
  R <- R_cont + R_cat
  
  
  mu_train_mean <- rep(0, times = N_train)
  mu_test_mean <- rep(0, times = N_test)
  beta_train_mean <- matrix(0, nrow = N_train, ncol = p+1)
  beta_test_mean <- matrix(0, nrow = N_test, ncol = p+1)
  
  varcounts_samples <- array(dim = c(n_chains*nd, R, p+1))
  timing <- rep(NA, times = n_chains)
  
  
  for(chain_num in 1:n_chains){
    print(paste("Starting chain", chain_num, "at", Sys.time()))
    train_time <- 
      system.time(
        fit <- VCBART::VCBART_ind(Y_train = Y_train,
                                  subj_id_train = subj_id_train, 
                                  ni_train = ni_train,
                                  X_train = X_train,
                                  Z_cont_train = Z_cont_train,
                                  X_test = X_test,
                                  Z_cont_test = Z_cont_test,
                                  unif_cuts = unif_cuts,
                                  cutpoints_list = cutpoints_list,
                                  cat_levels_list = cat_levels_list,
                                  sparse = sparse,
                                  M = M,
                                  mu0 = mu0, tau = tau, nu = nu, lambda = lambda,
                                  nd = nd, burn = burn, thin = 1,
                                  save_samples = save_samples, save_trees = save_trees,
                                  verbose = verbose, print_every = print_every))
    start_index <- (chain_num-1)*nd + 1
    end_index <- chain_num*nd
    
    mu_train_mean <- mu_train_mean + fit$yhat.train.mean
    mu_test_mean <- mu_test_mean + fit$yhat.test.mean
    beta_train_mean <- beta_train_mean + fit$betahat.train.mean
    beta_test_mean <- beta_test_mean + fit$betahat.test.mean
    varcounts_samples[start_index:end_index,,] <- fit$varcounts[-(1:burn),,]
    timing[chain_num] <- train_time["elapsed"]
    
  }
  
  mu_train_mean <- mu_train_mean/n_chains
  mu_test_mean <- mu_test_mean/n_chains
  beta_train_mean <- beta_train_mean/n_chains
  beta_test_mean <- beta_test_mean/n_chains

  selection_prob <- apply(varcounts_samples >= 1, FUN = mean, MARGIN = c(2,3))
  varcounts_mean <- apply(varcounts_samples, FUN = mean, MARGIN = c(2,3))
  
  
  return(
    list(time = timing, train_time = timing,
         mu_train_mean = mu_train_mean,
         mu_test_mean = mu_test_mean,
         beta_train_mean = beta_train_mean,
         beta_test_mean = beta_test_mean,
         varsel = list(varcounts_samples = varcounts_samples,selection_prob = selection_prob,
                       varcounts_means = varcounts_mean)))
}
