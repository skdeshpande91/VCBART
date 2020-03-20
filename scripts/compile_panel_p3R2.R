# Compile the results from the panel p3R2 simulations

p <- 3
R <- 2

sim_type <- "ar75"
methods <- c("true", "ind", paste0("ar_rho", 1:9), paste0("cs_rho", 1:9))

beta_mse_train <- array(dim = c(length(methods), p+1, 10), dimnames = list(methods, c(), c()))
beta_mse_test <- array(dim = c(length(methods), p+1, 10), dimnames = list(methods, c(), c()))
beta_coverage_train <- array(dim = c(length(methods), p+1, 10), dimnames = list(methods, c(), c()))
beta_coverage_test <- array(dim = c(length(methods), p+1, 10), dimnames = list(methods, c(), c()))
beta_int_train <- array(dim = c(length(methods), p+1, 10), dimnames = list(methods, c(), c()))
beta_int_test <- array(dim = c(length(methods), p+1, 10), dimnames = list(methods, c(), c()))

ystar_mse_train <- matrix(nrow = 10, ncol = length(methods), dimnames = list(c(), methods))
ystar_mse_test <- matrix(nrow = 10, ncol = length(methods), dimnames = list(c(), methods))
ystar_coverage_train <- matrix(nrow = 10, ncol = length(methods), dimnames = list(c(), methods))
ystar_coverage_test <- matrix(nrow = 10, ncol = length(methods), dimnames = list(c(), methods))
ystar_int_train <- matrix(nrow = 10, ncol = length(methods), dimnames = list(c(), methods))
ystar_int_test <- matrix(nrow = 10, ncol = length(methods), dimnames = list(c(), methods))

sigma_est <- matrix(nrow = 10, ncol = length(methods), dimnames = list(c(), methods))

timing <- matrix(nrow = 10, ncol = length(methods), dimnames = list(c(), methods))


for(rep in 1:10){
  
  if(file.exists(paste0("results/panel_p3R2", "_", sim_type, "_", rep, ".RData"))){
    load(paste0("results/panel_p3R2", "_", sim_type, "_", rep, ".RData"))
    
    beta_mse_train[,,rep] <- get(paste0("beta_mse_train_", sim_type, "_", rep))
    beta_coverage_train[,,rep] <- get(paste0("beta_coverage_train_", sim_type, "_", rep))
    beta_int_train[,,rep] <- get(paste0("beta_int_train_", sim_type, "_", rep))
    
    beta_mse_test[,,rep] <- get(paste0("beta_mse_test_", sim_type, "_", rep))
    beta_coverage_test[,,rep] <- get(paste0("beta_coverage_test_", sim_type, "_", rep))
    beta_int_test[,,rep] <- get(paste0("beta_int_test_", sim_type, "_", rep))
    
    ystar_mse_train[rep,] <- get(paste0("ystar_mse_train_", sim_type, "_",rep))
    ystar_coverage_train[rep,] <- get(paste0("ystar_coverage_train_", sim_type, "_", rep))
    ystar_int_train[rep,] <- get(paste0("ystar_int_train_", sim_type, "_", rep))
    
    ystar_mse_test[rep,] <- get(paste0("ystar_mse_test_", sim_type, "_",rep))
    ystar_coverage_test[rep,] <- get(paste0("ystar_coverage_test_", sim_type, "_", rep))
    ystar_int_test[rep,] <- get(paste0("ystar_int_test_", sim_type, "_", rep))
    
    sigma_est[rep,] <- get(paste0("sigma_est_", sim_type, "_", rep))
    timing[rep,] <- get(paste0("timing_", sim_type, "_", rep))
    
    rm_list <- paste0(c("beta_mse_train_", "beta_coverage_train_", "beta_int_train_",
                        "ystar_mse_train_", "ystar_coverage_train_", "ystar_int_train_",
                        "beta_mse_test_", "beta_coverage_test_", "beta_int_test_",
                        "ystar_mse_test_", "ystar_coverage_test_", "ystar_int_test_",
                        "sigma_est_", "timing_"), sim_type, "_", rep)
    
    rm(list = rm_list)
  }
}

assign(paste0("beta_mse_train_p3R2_", sim_type), beta_mse_train)
assign(paste0("beta_coverage_train_p3R2_", sim_type), beta_coverage_train)
assign(paste0("beta_int_train_p3R2_", sim_type), beta_int_train)

assign(paste0("beta_mse_test_p3R2_", sim_type), beta_mse_test)
assign(paste0("beta_coverage_test_p3R2_", sim_type), beta_coverage_test)
assign(paste0("beta_int_test_p3R2_", sim_type), beta_int_test)


assign(paste0("ystar_mse_train_p3R2_", sim_type), ystar_mse_train)
assign(paste0("ystar_coverage_train_p3R2_", sim_type), ystar_coverage_train)
assign(paste0("ystar_int_train_p3R2_", sim_type), ystar_int_train)

assign(paste0("ystar_mse_test_p3R2_", sim_type), ystar_mse_test)
assign(paste0("ystar_coverage_test_p3R2_", sim_type), ystar_coverage_test)
assign(paste0("ystar_int_test_p3R2_", sim_type), ystar_int_test)

assign(paste0("sigma_est_p3R2_", sim_type), sigma_est)
assign(paste0("timing_p3R2_", sim_type), timing)

save_list <- paste0(c("beta_mse_train_", "beta_coverage_train_", "beta_int_train_",
                      "ystar_mse_train_", "ystar_coverage_train_", "ystar_int_train_",
                      "beta_mse_test_", "beta_coverage_test_", "beta_int_test_",
                      "ystar_mse_test_", "ystar_coverage_test_", "ystar_int_test_",
                      "sigma_est_", "timing_"), "p3R2_", sim_type)
save(list = save_list, file = paste0("results/panel_p3R2_", sim_type, ".RData"))

