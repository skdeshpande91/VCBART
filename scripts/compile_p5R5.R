# Compile p5R5 simulation results

load("data/example_p5R5_data.RData")

methods <- c("train_mean", "lm", "np", "btvcm", "vc_bart", "bart", "rf", "gbm")
lin_methods <- c("lm", "np", "tvcm", "vc_bart")

########
# Prepare containers to hold all output
########

n_sim <- 50

beta_mse_train_p5R5 <- array(dim = c(length(lin_methods), p+1, n_sim), dimnames = list(lin_methods, c(), c()))
beta_mse_test_p5R5 <- array(dim = c(length(lin_methods), p+1, n_sim), dimnames = list(lin_methods, c(), c()))
beta_coverage_train_p5R5 <- array(dim = c(length(lin_methods), p+1, n_sim), dimnames = list(lin_methods, c(), c()))
beta_coverage_test_p5R5 <- array(dim = c(length(lin_methods), p+1, n_sim), dimnames = list(lin_methods, c(), c()))
beta_int_train_p5R5 <- array(dim = c(length(lin_methods), p+1, n_sim), dimnames = list(lin_methods, c(), c()))
beta_int_test_p5R5 <- array(dim = c(length(lin_methods), p+1, n_sim), dimnames = list(lin_methods, c(), c()))

ystar_mse_train_p5R5 <- matrix(nrow = n_sim, ncol = length(methods), dimnames = list(c(), methods))
ystar_mse_test_p5R5 <- matrix(nrow = n_sim, ncol = length(methods), dimnames = list(c(), methods))
ystar_coverage_train_p5R5 <- matrix(nrow = n_sim, ncol = length(methods), dimnames = list(c(), methods))
ystar_coverage_test_p5R5 <- matrix(nrow = n_sim, ncol = length(methods), dimnames = list(c(), methods))
ystar_int_train_p5R5 <- matrix(nrow = n_sim, ncol = length(methods), dimnames = list(c(), methods))
ystar_int_test_p5R5 <- matrix(nrow = n_sim, ncol = length(methods), dimnames = list(c(), methods))

timing_p5R5 <- matrix(nrow = n_sim, ncol = length(methods), dimnames = list(c(), methods))

for(rep in 1:n_sim){
  
  if(file.exists(paste0("results/results_p5R5_", rep, ".RData"))){
    load(paste0("results/results_p5R5_", rep, ".RData"))
    
    beta_mse_train_p5R5[,,rep] <- get(paste0("beta_mse_train_", rep))
    beta_mse_test_p5R5[,,rep] <- get(paste0("beta_mse_test_", rep))
    
    beta_coverage_train_p5R5[,,rep] <- get(paste0("beta_coverage_train_", rep))
    beta_coverage_test_p5R5[,,rep] <- get(paste0("beta_coverage_test_", rep))
    
    beta_int_train_p5R5[,,rep] <- get(paste0("beta_int_train_", rep))
    beta_int_test_p5R5[,,rep] <- get(paste0("beta_int_test_", rep))
    
    ystar_mse_train_p5R5[rep,] <- get(paste0("ystar_mse_train_", rep))
    ystar_mse_test_p5R5[rep,] <- get(paste0("ystar_mse_test_", rep))
    
    ystar_coverage_train_p5R5[rep,] <- get(paste0("ystar_coverage_train_", rep))
    ystar_coverage_test_p5R5[rep,] <- get(paste0("ystar_coverage_test_", rep))
    
    ystar_int_train_p5R5[rep,] <- get(paste0("ystar_int_train_", rep))
    ystar_int_test_p5R5[rep,] <- get(paste0("ystar_int_test_", rep))
    
    timing_p5R5[rep,] <- get(paste0("timing_", rep))
    
    rm_list <- paste0(c("beta_mse_train_", "beta_coverage_train_", "beta_int_train_",
                        "beta_mse_test_", "beta_coverage_test_", "beta_int_test_",
                        "ystar_mse_train_", "ystar_coverage_train_", "ystar_int_train_",
                        "ystar_mse_test_", "ystar_coverage_test_", "ystar_int_test_", "timing_"), rep)
    rm(list = rm_list)
  }
  
}
save_list <- paste0(c("beta_mse_train_", "beta_coverage_train_", "beta_int_train_",
                      "beta_mse_test_", "beta_coverage_test_", "beta_int_test_",
                      "ystar_mse_train_", "ystar_coverage_train_", "ystar_int_train_",
                      "ystar_mse_test_", "ystar_coverage_test_", "ystar_int_test_", "timing_"), "p5R5")

save(list = save_list, file = "results/results_p5R5.RData")
