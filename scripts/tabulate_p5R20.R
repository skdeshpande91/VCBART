# Tabulate results

p <- 5
R <- 20
N_sim <- 25

lin_methods <- c("vcbart_adapt", "vcbart_fixed", "lm", 
                 "kernel_smoothing", "boosted_tvcm", "tvc")
other_methods <- c("bart", "extraTrees", "gbm")

methods <- c(lin_methods, other_methods)

beta_mse_train <- array(dim = c(length(lin_methods), p+1, N_sim), dimnames = list(lin_methods, c(), c()))
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
  if(file.exists(data_filename)){
    load(data_filename)
    load(paste0("results/sim_p5R20/vcbart_adapt/vcbart_adapt_", sim_number, ".RData"))
    vcbart_adapt <- get(paste0("vcbart_adapt_", sim_number))
    
    adapt_beta_int_train <- vcbart_adapt[["train"]][["beta"]][,"U95",] - vcbart_adapt[["train"]][["beta"]][,"L95",]
    adapt_beta_int_test <- vcbart_adapt[["test"]][["beta"]][,"U95",] - vcbart_adapt[["test"]][["beta"]][,"L95",]
    
    adapt_ystar_int_train <- vcbart_adapt[["train"]][["ystar"]][,"U95"] - vcbart_adapt[["train"]][["ystar"]][,"L95"]
    adapt_ystar_int_test <- vcbart_adapt[["test"]][["ystar"]][,"U95"] - vcbart_adapt[["test"]][["ystar"]][,"L95"]
    
    
    for(m in methods){
      
      
      file_name <- paste0("results/sim_p5R20/", m, "/", m, "_", sim_number, ".RData")
      
      if(file.exists(file_name)){
        load(file_name)
        tmp_fit <- get(paste0(m, "_", sim_number))
        rm(list = paste0(m, "_", sim_number))
        if(m %in% lin_methods){
          
          beta_mse_train[m,,sim_number] <- colMeans( (beta_train - tmp_fit[["train"]][["beta"]][,"MEAN",])^2, na.rm = TRUE)
          beta_mse_test[m,,sim_number] <- colMeans( (beta_test - tmp_fit[["test"]][["beta"]][,"MEAN",])^2, na.rm = TRUE)
          
          beta_cov_train[m,,sim_number] <- 
            apply(
              (beta_train >= tmp_fit[["train"]][["beta"]][,"L95",]) & (beta_train <= tmp_fit[["train"]][["beta"]][,"U95",]),
              FUN = mean, MARGIN = 2, na.rm = TRUE)
          beta_cov_test[m,,sim_number] <- 
            apply(
              (beta_test >= tmp_fit[["test"]][["beta"]][,"L95",]) & (beta_test <= tmp_fit[["test"]][["beta"]][,"U95",]),
              FUN = mean, MARGIN = 2, na.rm = TRUE)
          
          
          tmp_beta_int_train <- (tmp_fit[["train"]][["beta"]][,"U95",] - tmp_fit[["train"]][["beta"]][,"L95",])
          tmp_beta_int_test <- (tmp_fit[["test"]][["beta"]][,"U95",] - tmp_fit[["test"]][["beta"]][,"L95",])
          
          beta_int_train[m,,sim_number] <- 
            apply(tmp_beta_int_train/adapt_beta_int_train, FUN = mean, MARGIN = 2, na.rm = TRUE)
          beta_int_test[m,,sim_number] <- 
            apply(tmp_beta_int_test/adapt_beta_int_test, FUN = mean, MARGIN = 2, na.rm = TRUE)
        }
        
        
        tmp_mse_train <- mean( (Y_train - tmp_fit[["train"]][["ystar"]][,"MEAN"])^2, na.rm = TRUE)
        tmp_mse_test <- mean( (Y_test - tmp_fit[["test"]][["ystar"]][,"MEAN"])^2, na.rm = TRUE)
        
        tmp_ystar_int_train <- (tmp_fit[["train"]][["ystar"]][,"U95"] - tmp_fit[["train"]][["ystar"]][,"L95"])
        tmp_ystar_int_test <- (tmp_fit[["test"]][["ystar"]][,"U95"] - tmp_fit[["test"]][["ystar"]][,"L95"])
        
        
        ystar_rmse_train[m, sim_number] <- sqrt(tmp_mse_train)
        ystar_rmse_test[m, sim_number] <- sqrt(tmp_mse_test)
        
        ystar_smse_train[m, sim_number] <- tmp_mse_train/mean( (Y_train - mean(Y_train))^2, na.rm = TRUE)
        ystar_smse_test[m, sim_number] <- tmp_mse_test/mean( (Y_test - mean(Y_train))^2, na.rm = TRUE)
        
        ystar_cov_train[m, sim_number] <- mean( (Y_train >= tmp_fit[["train"]][["ystar"]][,"L95"]) &
                                                  (Y_train <= tmp_fit[["train"]][["ystar"]][,"U95"]),
                                                na.rm = TRUE)
        
        ystar_cov_test[m, sim_number] <- mean( (Y_test >= tmp_fit[["test"]][["ystar"]][,"L95"]) &
                                                 (Y_test <= tmp_fit[["test"]][["ystar"]][,"U95"]),
                                               na.rm = TRUE)
        
        ystar_int_train[m, sim_number] <- mean( tmp_ystar_int_train/adapt_ystar_int_train, na.rm = TRUE)
        ystar_int_test[m, sim_number] <- mean( tmp_ystar_int_test/adapt_ystar_int_test, na.rm = TRUE)
        
        if(length(tmp_fit$time) > 1) timing[m, sim_number] <- tmp_fit$time["elapsed"]
        else timing[m, sim_number] <- tmp_fit$time
      }
    }
    
    
    # no need to keep Y_train, X_train, Z_train, etc.
    rm(Y_train, X_train, Z_train, beta_train, X_test, Z_test, beta_test)
    
  }
}

save(beta_mse_train, beta_mse_test, beta_cov_train, beta_cov_test, beta_int_train, beta_int_test,
     ystar_rmse_train, ystar_rmse_test, ystar_smse_train, ystar_smse_test, 
     ystar_cov_train, ystar_cov_test, ystar_int_train, ystar_int_test, file = paste0("results/sim_p5R20/results_p5R20.RData"))
