# mean only
vcbart_cs_pred_wrapper <- function(Y_train,
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
                                   M = 50, mu0 = NULL, tau = NULL,
                                   nd = 1000, burn = 1000, thin = 1,
                                   verbose = TRUE, print_every = floor( (nd * thin + burn)/10),
                                   n_chains = 4)
{
  p <- ncol(X_train)
  N_train <- nrow(X_train)
  N_test <- nrow(X_test)

  if(length(Z_cont_train) != 1) R_cont <- ncol(Z_cont_train)
  else R_cont <- 0
  
  if(length(Z_cat_train) != 1) R_cat <- ncol(Z_cat_train)
  else R_cat <- 0
  R <- R_cont + R_cat

  betahat_train <- matrix(0, nrow = N_train, ncol = p+1)
  betahat_test <- matrix(0, nrow = N_test, ncol = p+1)
  
  yhat_train <- rep(0, times = N_train)
  yhat_test <- rep(0, times = N_test)
  
  varcounts_samples <- array(dim = c(n_chains*nd, R, p+1))
  theta_samples <- array(dim = c(n_chains*nd, R, p+1))
  timing <- rep(NA, times = n_chains)
  
  if(is.null(mu0)) mu0 <- rep(0, times = p+1)
  if(is.null(tau)) tau <- rep(0.5/sqrt(M), times = p+1)
  
  for(chain_num in 1:n_chains){
    print(paste("Starting chain", chain_num, "at", Sys.time()))
    train_time <- 
      system.time(
        fit <- VCBART::VCBART_cs(Y_train = Y_train,
                                 subj_id_train = subj_id_train,
                                 ni_train = ni_train,
                                 X_train = X_train,
                                 Z_cont_train = Z_cont_train,
                                 Z_cat_train = Z_cat_train,
                                 X_test = X_test,
                                 Z_cont_test = Z_cont_test,
                                 Z_cat_test = Z_cat_test,
                                 unif_cuts = unif_cuts,
                                 cutpoints_list = cutpoints_list,
                                 cat_levels_list = cat_levels_list,
                                 sparse = TRUE,
                                 M = 50, mu0 = mu0, tau = tau,
                                 nd = nd, burn = burn, thin = thin,
                                 save_trees = FALSE, save_samples = FALSE,
                                 verbose = verbose))
        
        

    start_index <- (chain_num-1)*nd + 1
    end_index <- chain_num*nd
    
    betahat_train <- betahat_train + fit$betahat.train.mean
    betahat_test <- betahat_test + fit$betahat.test.mean
    
    yhat_train <- yhat_train + fit$yhat.train.mean
    yhat_test <- yhat_test + fit$yhat.test.mean
    
    varcounts_samples[start_index:end_index,,] <- fit$varcounts[-(1:burn),,]
    timing[chain_num] <- train_time["elapsed"]
    
  }
 
  betahat_train <- betahat_train/n_chains
  betahat_test <- betahat_test/n_chains
  yhat_train <- yhat_train/n_chains
  yhat_test <- yhat_test/n_chains
  
  selection_prob <- apply(varcounts_samples >= 1, FUN = mean, MARGIN = c(2,3))
  varcounts_mean <- apply(varcounts_samples, FUN = mean, MARGIN = c(2,3))
  
  
  return(
    list(time = timing, train_time = timing,
         train = list(beta = betahat_train, ystar = yhat_train),
         test = list(beta = betahat_test, ystar = yhat_test),
         varsel = list(selection_prob = selection_prob,
                       varcounts_means = varcounts_mean)))
}
