vcbart_cs_wrapper <- function(Y_train,
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
                              edge_mat_list = NULL,
                              graph_split = rep(FALSE, times = ncol(Z_cat_train)),
                              sparse = TRUE,
                              M = 50, 
                              mu0 = NULL, tau = NULL, nu = NULL, lambda = NULL,
                              nd = 1000, burn = 1000, thin = 1,
                              save_samples = TRUE, save_trees = FALSE,
                              verbose = TRUE, print_every = floor( (nd * thin + burn)/10),
                              n_chains = 4)
{
  N_train <- nrow(X_train)
  
  if(is.null(tau)) tau <- rep(0.5/sqrt(M), times = p+1)
  p <- ncol(X_train)
  
  if(length(Z_cont_train) != 1) R_cont <- ncol(Z_cont_train)
  else R_cont <- 0
  
  if(length(Z_cat_train) != 1) R_cat <- ncol(Z_cat_train)
  else R_cat <- 0
  R <- R_cont + R_cat
  
  beta_train_samples <- array(dim = c(n_chains*nd, N_train, p+1))
  timing <- rep(NA, times = n_chains)
  
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
                                 edge_mat_list = edge_mat_list,
                                 graph_split = graph_split,
                                 sparse = sparse,
                                 M = M,
                                 mu0 = mu0, tau = tau, nu = nu, lambda = lambda,
                                 nd = nd, burn = burn, thin = 1,
                                 save_samples = save_samples, save_trees = save_trees,
                                 verbose = verbose, print_every = print_every))
    start_index <- (chain_num-1)*nd + 1
    end_index <- chain_num*nd
    beta_train_samples[start_index:end_index,,] <- fit$betahat.train
    timing[chain_num] <- train_time["elapsed"]
  }
  
  # Now build containers to summarize all of the posterior samples of beta
  beta_summary_train <- array(dim = c(N_train, 3, p+1), dimnames = list(c(), c("MEAN", "L95", "U95"), c()))
  
  beta_summary_train[,"MEAN",] <- apply(beta_train_samples, MARGIN = c(2,3), FUN = mean)
  beta_summary_train[,"L95",] <- apply(beta_train_samples, MARGIN = c(2,3), FUN = quantile, probs = 0.025)
  beta_summary_train[,"U95",] <- apply(beta_train_samples, MARGIN = c(2,3), FUN = quantile, probs = 0.975)
  
  return(list(time = timing, train_time = timing,
              train = list(beta = beta_summary_train)))
}
