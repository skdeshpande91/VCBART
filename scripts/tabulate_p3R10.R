# Tabulate results

load("data/p3R10_data.RData")
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
  
  # Get training/testing splits
  set.seed(129 + sim_number)
  index <- sample(1:N, size = 500, replace = FALSE)
  train_index <- sort(index[1:375])
  test_index <- sort(index[376:500])
  
  
  X_train <- as.matrix(X_all[train_index,], nrow = length(train_index), ncol = p)
  Y_train <- Y_all[train_index]
  Z_train <- as.matrix(Z_all[train_index,], ncol = R)
  
  X_test <- as.matrix(X_all[test_index,], nrow = length(test_index), ncol = p)
  Y_test <- Y_all[test_index]
  Z_test <- as.matrix(Z_all[test_index,], ncol = R)
  
  beta_train <- beta_all[train_index,]
  beta_test <- beta_all[test_index,]
  
  # comparing against vcbart_adapt
  # load in the vcbart adapt results
  load(paste0("results/sim_p3R10/vcbart_adapt/new_vcbart_adapt_", sim_number, ".RData"))
  vcbart_adapt <- get(paste0("vcbart_adapt_", sim_number))
  
  adapt_beta_int_train <- vcbart_adapt[["train"]][["beta"]][,"U95",] - vcbart_adapt[["train"]][["beta"]][,"L95",]
  adapt_beta_int_test <- vcbart_adapt[["test"]][["beta"]][,"U95",] - vcbart_adapt[["test"]][["beta"]][,"L95",]
  
  adapt_ystar_int_train <- vcbart_adapt[["train"]][["ystar"]][,"U95"] - vcbart_adapt[["train"]][["ystar"]][,"L95"]
  adapt_ystar_int_test <- vcbart_adapt[["test"]][["ystar"]][,"U95"] - vcbart_adapt[["test"]][["ystar"]][,"L95"]
  
  
  for(m in methods){

    
    if(m %in% c("vcbart_fixed", "vcbart_adapt")){
      file_name <- paste0("results/sim_p3R10/", m, "/new_", m, "_", sim_number, ".RData")
    } else{
      file_name <- paste0("results/sim_p3R10/", m, "/", m, "_", sim_number, ".RData")
    }
    
    if(file.exists(file_name)){
      load(file_name)
      tmp_fit <- get(paste0(m, "_", sim_number))
      
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
}

save(beta_mse_train, beta_mse_test, beta_cov_train, beta_cov_test, beta_int_train, beta_int_test,
     ystar_rmse_train, ystar_rmse_test, ystar_smse_train, ystar_smse_test, 
     ystar_cov_train, ystar_cov_test, ystar_int_train, ystar_int_test, file = paste0("results/sim_p3R10/results_p3R10.RData"))
