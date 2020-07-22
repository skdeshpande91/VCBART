p <- 5
R <- 20
N_sim <- 25

tau_names <- paste0("tau_", c("1_1", "1_4", "1_2", "2_3", "3_2", "2_1", "4_1"))
methods <- tau_names

beta_mse_train <- array(dim = c(length(methods), p+1, N_sim), dimnames = list(methods, c(), c()))
beta_mse_test <- beta_mse_train

beta_cov_train <- beta_mse_train
beta_cov_test <- beta_mse_train

beta_int_train <- beta_mse_train
beta_int_test <- beta_mse_train

ystar_smse_train <- matrix(nrow = length(methods), ncol = N_sim,
                           dimnames = list(methods, c()))
ystar_smse_test <- ystar_smse_train

ystar_rmse_train <- ystar_smse_train
ystar_rmse_test <- ystar_smse_test

ystar_cov_train <- ystar_smse_train
ystar_cov_test <- ystar_smse_train

ystar_int_train <- ystar_smse_train
ystar_int_test <- ystar_smse_train

timing <- matrix(nrow = length(methods), ncol = N_sim,
                 dimnames = list(methods, c()))

for(sim_number in 1:N_sim){
  data_filename <- paste0("data/sim_p5R20/data_p5R20_", sim_number, ".RData")
  results_filename <- paste0("results/sim_p5R20_tau/results_tau_", sim_number, ".RData")
  
  if(file.exists(data_filename) & file.exists(results_filename)){
    load(data_filename)
    load(results_filename)
    
    
    orig_fit <- get(paste0("vcbart_tau_1_1_", sim_number))
    true_beta_int_train <- orig_fit[["train"]][["beta"]][,"U95",] - orig_fit[["train"]][["beta"]][,"L95",]
    true_beta_int_test <- orig_fit[["test"]][["beta"]][,"U95",] - orig_fit[["test"]][["beta"]][,"L95",]
    
    true_ystar_int_train <- orig_fit[["train"]][["ystar"]][,"U95"] - orig_fit[["train"]][["ystar"]][,"L95"]
    true_ystar_int_test <- orig_fit[["test"]][["ystar"]][,"U95"] - orig_fit[["test"]][["ystar"]][,"L95"]
    
    
    for(ix in 1:length(tau_names)){
      tmp_fit <- get(paste0("vcbart_", tau_names[ix], "_", sim_number))
      
      tmp_beta_int_train <- (tmp_fit[["train"]][["beta"]][,"U95",] - tmp_fit[["train"]][["beta"]][,"L95",])
      tmp_beta_int_test <- (tmp_fit[["test"]][["beta"]][,"U95",] - tmp_fit[["test"]][["beta"]][,"L95",])
      
      
      beta_mse_train[ix,,sim_number] <- colMeans( (beta_train - tmp_fit[["train"]][["beta"]][,"MEAN",])^2, na.rm = TRUE) 
      beta_mse_test[ix,,sim_number] <- colMeans( (beta_test - tmp_fit[["test"]][["beta"]][,"MEAN",])^2, na.rm = TRUE) 
      
      beta_cov_train[ix,,sim_number] <- 
        apply(
          (beta_train >= tmp_fit[["train"]][["beta"]][,"L95",]) & (beta_train <= tmp_fit[["train"]][["beta"]][,"U95",]),
          FUN = mean, MARGIN = 2, na.rm = TRUE)
      beta_cov_test[ix,,sim_number] <- 
        apply(
          (beta_test >= tmp_fit[["test"]][["beta"]][,"L95",]) & (beta_test <= tmp_fit[["test"]][["beta"]][,"U95",]),
          FUN = mean, MARGIN = 2, na.rm = TRUE)
      
      beta_int_train[ix,,sim_number] <- 
        apply(tmp_beta_int_train/true_beta_int_train, FUN = mean, MARGIN = 2, na.rm = TRUE)
      beta_int_test[ix,,sim_number] <- 
        apply(tmp_beta_int_test/true_beta_int_test, FUN = mean, MARGIN = 2, na.rm = TRUE)
      
      
      tmp_mse_train <- mean( (Y_train - tmp_fit[["train"]][["ystar"]][,"MEAN"])^2, na.rm = TRUE)
      tmp_mse_test <- mean( (Y_test - tmp_fit[["test"]][["ystar"]][,"MEAN"])^2, na.rm = TRUE)
      
      tmp_ystar_int_train <- (tmp_fit[["train"]][["ystar"]][,"U95"] - tmp_fit[["train"]][["ystar"]][,"L95"])
      tmp_ystar_int_test <- (tmp_fit[["test"]][["ystar"]][,"U95"] - tmp_fit[["test"]][["ystar"]][,"L95"])
      
      
      ystar_rmse_train[ix, sim_number] <- sqrt(tmp_mse_train)
      ystar_rmse_test[ix, sim_number] <- sqrt(tmp_mse_test)
      
      ystar_smse_train[ix, sim_number] <- tmp_mse_train/mean( (Y_train - mean(Y_train))^2, na.rm = TRUE)
      ystar_smse_test[ix, sim_number] <- tmp_mse_test/mean( (Y_test - mean(Y_train))^2, na.rm = TRUE)
      
      ystar_cov_train[ix, sim_number] <- mean( (Y_train >= tmp_fit[["train"]][["ystar"]][,"L95"]) &
                                                     (Y_train <= tmp_fit[["train"]][["ystar"]][,"U95"]),
                                                   na.rm = TRUE)
      
      ystar_cov_test[ix, sim_number] <- mean( (Y_test >= tmp_fit[["test"]][["ystar"]][,"L95"]) &
                                                    (Y_test <= tmp_fit[["test"]][["ystar"]][,"U95"]),
                                                  na.rm = TRUE)
      
      ystar_int_train[ix, sim_number] <- mean( tmp_ystar_int_train/true_ystar_int_train, na.rm = TRUE)
      ystar_int_test[ix, sim_number] <- mean( tmp_ystar_int_test/true_ystar_int_test, na.rm = TRUE)
      
      timing[ix, sim_number] <- tmp_fit$time
    }
    
  }
}


save(beta_mse_train, beta_mse_test, beta_cov_train, beta_cov_test, beta_int_train, beta_int_test,
     ystar_rmse_train, ystar_rmse_test, ystar_smse_train, ystar_smse_test, 
     ystar_cov_train, ystar_cov_test, ystar_int_train, ystar_int_test, file = paste0("results/sim_p5R20_tau/results_p5R20_tau.RData"))

