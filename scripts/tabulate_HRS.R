# Tabulate results
# for HRS

N_sim <- 25


adapt_cs_methods <- paste0("vcbart_adapt_cs", c(25, 50, 75))
fixed_cs_methods <- paste0("vcbart_fixed_cs", c(25, 50, 75))
methods <- c("vcbart_adapt", adapt_cs_methods, "vcbart_fixed", fixed_cs_methods, 
             "lm",  "boosted_tvcm", "bart", "extraTrees", "gbm")


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
  print(sim_number)
  # Load the data
  load(paste0("data/HRS/HRS_data_", sim_number, ".RData"))
  
  # Load the VCBART-adapt data
  if(file.exists(paste0("results/HRS/vcbart_adapt/vcbart_adapt_", sim_number, ".RData"))){
    load(paste0("results/HRS/vcbart_adapt/vcbart_adapt_", sim_number, ".RData"))
    
    vcbart_adapt <- get(paste0("vcbart_adapt_", sim_number))
    adapt_ystar_int_train <- vcbart_adapt[["train"]][["ystar"]][,"U95"] - vcbart_adapt[["train"]][["ystar"]][,"L95"]
    adapt_ystar_int_test <- vcbart_adapt[["test"]][["ystar"]][,"U95"] - vcbart_adapt[["test"]][["ystar"]][,"L95"]
    rm(vcbart_adapt) # delete the full vcbart_adapt object
    
    for(m in methods){
      
      if(m %in% adapt_cs_methods){
        file_name <- paste0("results/HRS/vcbart_adapt/", m, "_", sim_number, ".RData")
      } else if(m %in% fixed_cs_methods){
        file_name <- paste0("results/HRS/vcbart_fixed/", m, "_", sim_number, ".RData")
      } else{
        file_name <- paste0("results/HRS/",m, "/", m, "_", sim_number, ".RData")
      }

      
      if(file.exists(file_name)){
        load(file_name)
        tmp_fit <- get(paste0(m, "_", sim_number))
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
        rm(list = paste0(m, "_", sim_number)) # get rid of the temporary fits as we go
      }
    }
  }
}

save(ystar_rmse_train, ystar_rmse_test, ystar_smse_train, ystar_smse_test, 
     ystar_cov_train, ystar_cov_test, ystar_int_train, ystar_int_test, file = paste0("results/HRS/results_HRS.RData"))
