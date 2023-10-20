source("true_betas_p5R20.R")
source("vcbart_wrapper.R")
source("vcbart_pred_wrapper.R")
source("assess_support_recovery.R")
p <- 5

n_train_list <- c(25, 50, 100, 125, 250, 500,
                  1000, 2500, 5000)

#for(n_train in c(50, 100, 125, 250, 500)){
#for(n_train in c(5000)){
for(n_train in n_train_list){
  n_test <- 250
  n_sim <- 25
  varsel_cols <- c("tp", "fp", "fn", "tn",
                   "sen", "spec", "prec", "acc",
                   "mcc", "f1")
  
  ystar_rmse_train <- rep(NA, times = n_sim)
  ystar_rmse_test <- rep(NA, times = n_sim)
  ystar_cov_train <- rep(NA, times = n_sim)
  ystar_cov_test <- rep(NA, times = n_sim)
  
  beta_mse_train <- matrix(nrow = n_sim, ncol = p+1)
  beta_mse_test <- matrix(nrow = n_sim, ncol = p+1)
  beta_cov_train <- matrix(nrow = n_sim, ncol = p+1)
  beta_cov_test <- matrix(nrow = n_sim, ncol = p+1)
  
  timing <- rep(NA, times = n_sim)
  
  varsel_results <- matrix(nrow = n_sim,
                           ncol = 10,
                           dimnames = list(c(), varsel_cols))
  
  for(sim_number in 1:n_sim){
    set.seed(6821 + sim_number)
    source("generate_data.R")
    M <- 50
    tau <- rep(0.5/sqrt(M), times = p+1)
    print(paste("Simulation", sim_number, "at", Sys.time()))
    fit <- vcbart_wrapper(Y_train = Y_train,
                          subj_id_train = subj_id_train,
                          ni_train = ni_train,
                          X_train = X_train,
                          Z_cont_train = Z_cont_train,
                          X_test = X_test,
                          Z_cont_test = Z_cont_test,
                          unif_cuts = unif_cuts,
                          sparse = TRUE,
                          M = M,tau = tau,
                          verbose = FALSE,
                          n_chains = 4)
    
    beta_mse_train[sim_number,] <- 
      colMeans( (beta_train - fit$train$beta[,"MEAN",])^2 )
    beta_cov_train[sim_number,] <- 
      colMeans( (beta_train >= fit$train$beta[,"L95",] &
                   beta_train <= fit$train$beta[,"U95",]))
    
    beta_mse_test[sim_number,] <- 
      colMeans( (beta_test - fit$test$beta[,"MEAN",])^2 )
    beta_cov_test[sim_number,] <- 
      colMeans( (beta_test >= fit$test$beta[,"L95",] &
                   beta_test <= fit$test$beta[,"U95",]))
    
    ystar_rmse_train[sim_number] <- sqrt(mean( (Y_train - fit$train$ystar[,"MEAN"])^2 ))
    ystar_cov_train[sim_number] <- mean( (Y_train >= fit$train$ystar[,"L95"] &
                                            Y_train <= fit$train$ystar[,"U95"]))
    
    ystar_rmse_test[sim_number] <- sqrt(mean( (Y_test - fit$test$ystar[,"MEAN"])^2 ))
    ystar_cov_test[sim_number] <- mean( (Y_test >= fit$test$ystar[,"L95"] &
                                           Y_test <= fit$test$ystar[,"U95"]))
    
    
    tmp <- assess_support_recovery(fit$varsel$varcounts_means,
                                   fit$varsel$selection_prob,
                                   true_support, R)
    varsel_results[sim_number,] <- tmp$varcount_prob[51,]
    timing[sim_number] <- sum(fit$train_time)
  }
  assign(paste0("results_n", n_train),
         list(ystar_rmse_train = ystar_rmse_train,
              ystar_rmse_test = ystar_rmse_test,
              ystar_cov_train = ystar_cov_train,
              ystar_cov_test = ystar_cov_test,
              beta_mse_train = beta_mse_train,
              beta_mse_test = beta_mse_test,
              beta_cov_train = beta_cov_train,
              beta_cov_test = beta_cov_test,
              varsel = varsel_results,
              timing = timing))
  save(list = paste0("results_n", n_train), file = paste0("p5R20_demo_n", n_train, ".RData"))
}



# for n_train > 10000 we will use a different wrapper that only tracks posterior means
# storing the samples for the training data takes up too much memory
for(n_train in c(10000, 12500)){
  n_test <- 250
  n_sim <- 25
  varsel_cols <- c("tp", "fp", "fn", "tn",
                   "sen", "spec", "prec", "acc",
                   "mcc", "f1")
  
  ystar_rmse_train <- rep(NA, times = n_sim)
  ystar_rmse_test <- rep(NA, times = n_sim)
  
  beta_mse_train <- matrix(nrow = n_sim, ncol = p+1)
  beta_mse_test <- matrix(nrow = n_sim, ncol = p+1)
  
  timing <- rep(NA, times = n_sim)
  
  varsel_results <- matrix(nrow = n_sim,
                           ncol = 10,
                           dimnames = list(c(), varsel_cols))
  
  for(sim_number in 1:n_sim){
    set.seed(6821 + sim_number)
    source("generate_data.R")
    M <- 50
    tau <- rep(0.5/sqrt(M), times = p+1)
    print(paste("Simulation", sim_number, "at", Sys.time()))
    fit <- vcbart_pred_wrapper(Y_train = Y_train,
                               subj_id_train = subj_id_train,
                               ni_train = ni_train,
                               X_train = X_train,
                               Z_cont_train = Z_cont_train,
                               X_test = X_test,
                               Z_cont_test = Z_cont_test,
                               unif_cuts = unif_cuts,
                               sparse = TRUE,
                               M = M,tau = tau,
                               verbose = FALSE,
                               n_chains = 4)
    
    beta_mse_train[sim_number,] <- colMeans( (beta_train - fit$beta_train_mean)^2)
    beta_mse_test[sim_number,] <- colMeans( (beta_test - fit$beta_test_mean)^2)
    
    
    ystar_rmse_train[sim_number] <- sqrt(mean( (Y_train - fit$mu_train_mean)^2 ))
    ystar_rmse_test[sim_number] <- sqrt(mean( (Y_test - fit$mu_test_mean)^2))
    
    tmp <- assess_support_recovery(fit$varsel$varcounts_means,
                                   fit$varsel$selection_prob,
                                   true_support, R)
    varsel_results[sim_number,] <- tmp$varcount_prob[51,]
    timing[sim_number] <- sum(fit$train_time)
  }
  assign(paste0("results_n", n_train),
         list(ystar_rmse_train = ystar_rmse_train,
              ystar_rmse_test = ystar_rmse_test,
              beta_mse_train = beta_mse_train,
              beta_mse_test = beta_mse_test,
              varsel = varsel_results,
              timing = timing))
  save(list = paste0("results_n", n_train), file = paste0("p5R20_demo_n", n_train, ".RData"))
}
