# Tabulate the panel simulations
rho_true <- 0.75
rho_list <- c(rho_true, seq(0.0, 0.9, by = 0.1))
p <- 5
R <- 20

N_sim <- 25

methods <- paste0("vcbart_adapt_cs", 100*rho_list)


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
  data_filename <- paste0("data/sim_p5R20_rho75/data_p5R20_rho75_", sim_number, ".RData")
  if(file.exists(data_filename)){
    load(data_filename)
    load(paste0("results/sim_p5R20_rho75/results_sim", sim_number, ".RData"))
    
    true_fit <- get(paste0("vcbart_adapt_cs", 100*rho_true, "_", sim_number))
    true_beta_int_train <- true_fit[["train"]][["beta"]][,"U95",] - true_fit[["train"]][["beta"]][,"L95",]
    true_beta_int_test <- true_fit[["test"]][["beta"]][,"U95",] - true_fit[["test"]][["beta"]][,"L95",]

    true_ystar_int_train <- true_fit[["train"]][["ystar"]][,"U95"] - true_fit[["train"]][["ystar"]][,"L95"]
    true_ystar_int_test <- true_fit[["test"]][["ystar"]][,"U95"] - true_fit[["test"]][["ystar"]][,"L95"]
    
    for(rho_ix in 1:length(rho_list)){
      rho <- rho_list[rho_ix]
      tmp_fit <- get(paste0("vcbart_adapt_cs", 100*rho, "_", sim_number))
      
      tmp_beta_int_train <- (tmp_fit[["train"]][["beta"]][,"U95",] - tmp_fit[["train"]][["beta"]][,"L95",])
      tmp_beta_int_test <- (tmp_fit[["test"]][["beta"]][,"U95",] - tmp_fit[["test"]][["beta"]][,"L95",])
      
      
      beta_mse_train[rho_ix,,sim_number] <- colMeans( (beta_train - tmp_fit[["train"]][["beta"]][,"MEAN",])^2, na.rm = TRUE) 
      beta_mse_test[rho_ix,,sim_number] <- colMeans( (beta_test - tmp_fit[["test"]][["beta"]][,"MEAN",])^2, na.rm = TRUE) 
      
      beta_cov_train[rho_ix,,sim_number] <- 
        apply(
          (beta_train >= tmp_fit[["train"]][["beta"]][,"L95",]) & (beta_train <= tmp_fit[["train"]][["beta"]][,"U95",]),
          FUN = mean, MARGIN = 2, na.rm = TRUE)
      beta_cov_test[rho_ix,,sim_number] <- 
        apply(
          (beta_test >= tmp_fit[["test"]][["beta"]][,"L95",]) & (beta_test <= tmp_fit[["test"]][["beta"]][,"U95",]),
          FUN = mean, MARGIN = 2, na.rm = TRUE)
      
      beta_int_train[rho_ix,,sim_number] <- 
        apply(tmp_beta_int_train/true_beta_int_train, FUN = mean, MARGIN = 2, na.rm = TRUE)
      beta_int_test[rho_ix,,sim_number] <- 
        apply(tmp_beta_int_test/true_beta_int_test, FUN = mean, MARGIN = 2, na.rm = TRUE)
      
      
      tmp_mse_train <- mean( (Y_train - tmp_fit[["train"]][["ystar"]][,"MEAN"])^2, na.rm = TRUE)
      tmp_mse_test <- mean( (Y_test - tmp_fit[["test"]][["ystar"]][,"MEAN"])^2, na.rm = TRUE)
      
      tmp_ystar_int_train <- (tmp_fit[["train"]][["ystar"]][,"U95"] - tmp_fit[["train"]][["ystar"]][,"L95"])
      tmp_ystar_int_test <- (tmp_fit[["test"]][["ystar"]][,"U95"] - tmp_fit[["test"]][["ystar"]][,"L95"])
      
      
      ystar_rmse_train[rho_ix, sim_number] <- sqrt(tmp_mse_train)
      ystar_rmse_test[rho_ix, sim_number] <- sqrt(tmp_mse_test)
      
      ystar_smse_train[rho_ix, sim_number] <- tmp_mse_train/mean( (Y_train - mean(Y_train))^2, na.rm = TRUE)
      ystar_smse_test[rho_ix, sim_number] <- tmp_mse_test/mean( (Y_test - mean(Y_train))^2, na.rm = TRUE)
      
      ystar_cov_train[rho_ix, sim_number] <- mean( (Y_train >= tmp_fit[["train"]][["ystar"]][,"L95"]) &
                                                (Y_train <= tmp_fit[["train"]][["ystar"]][,"U95"]),
                                              na.rm = TRUE)
      
      ystar_cov_test[rho_ix, sim_number] <- mean( (Y_test >= tmp_fit[["test"]][["ystar"]][,"L95"]) &
                                               (Y_test <= tmp_fit[["test"]][["ystar"]][,"U95"]),
                                             na.rm = TRUE)
      
      ystar_int_train[rho_ix, sim_number] <- mean( tmp_ystar_int_train/true_ystar_int_train, na.rm = TRUE)
      ystar_int_test[rho_ix, sim_number] <- mean( tmp_ystar_int_test/true_ystar_int_test, na.rm = TRUE)
      
      timing[rho_ix, sim_number] <- tmp_fit$time
    }
  }
}

save(beta_mse_train, beta_mse_test, beta_cov_train, beta_cov_test, beta_int_train, beta_int_test,
     ystar_rmse_train, ystar_rmse_test, ystar_smse_train, ystar_smse_test, 
     ystar_cov_train, ystar_cov_test, ystar_int_train, ystar_int_test, file = paste0("results/sim_p5R20_rho75/results_p5R20_rho75.RData"))
