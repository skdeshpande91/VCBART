load("sim_settings.RData")
source("true_betas_p5R20.R")

source("lm_wrapper.R")
source("vcbart_wrapper.R")
source("tvc_wrapper.R")
source("kernel_smoothing_wrapper.R")
source("boosted_tvcm_wrapper.R")
source("bart_wrapper.R")
source("gbm_wrapper.R")
source("extraTrees_wrapper.R")

##############
#args <- commandArgs(TRUE)
#job_id <- as.numeric(args[1])
#sim_number <- sim_settings[job_id,"sim_number"]
#m <- as.character(sim_settings[job_id, "method"])
###############

sim_number <- 1
m <- "vcbart"
# m can also be one of
# "btvcm", "tvc", "ks", "lm", "gbm", "ert"

#################
# Generate data
#################
n_train <- 250
n_test <- 25
set.seed(6821 + sim_number)
source("generate_data.R")

#####################

if(m == "vcbart"){
  #M <- 200
  M <- 50
  tau <- rep(0.5/sqrt(M), times = p+1)
  
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
} else if(m == "btvcm"){
  fit <- boosted_tvcm_wrapper(Y_train, X_train, Z_cont_train, 
                              X_test, Z_cont_test, B = 50)
} else if(m == "tvc"){
  fit <- tvc_wrapper(Y_train, X_train, Z_cont_train, 
                         X_test, Z_cont_test,B = 50)
  
} else if(m == "ks"){
  fit <- kernel_smoothing_wrapper(Y_train, cov_train, mod_train, 
                                  cov_test, mod_test, B = 50)
  
} else if(m == "lm"){
  fit <- lm_wrapper(Y_train, X_train, Z_cont_train, X_test, Z_cont_test)
  
} else if(m == "bart"){
  fit <- bart_wrapper(Y_train, X_train, Z_cont_train, X_test, Z_cont_test,
                      nd = 1000, burn = 1000, thin = 1)
} else if(m == "gbm"){
  fit <- gbm_wrapper(Y_train, X_train, Z_cont_train, X_test, Z_cont_test)
} else if(m == "ert"){
  fit <- extraTrees_wrapper(Y_train, X_train, Z_cont_train, X_test, Z_cont_test)
}


ystar_rmse_train <- sqrt(mean( (Y_train - fit$train$ystar[,"MEAN"])^2 ))
ystar_cov_train <- mean( (Y_train >= fit$train$ystar[,"L95"] &
                               Y_train <= fit$train$ystar[,"U95"]))

ystar_rmse_test <- sqrt(mean( (Y_test - fit$test$ystar[,"MEAN"])^2 ))
ystar_cov_test <- mean( (Y_test >= fit$test$ystar[,"L95"] &
                              Y_test <= fit$test$ystar[,"U95"]))
train_timing <- fit$train_time
timing <- fit$time

results <- list(method = m,
                sim_number = sim_number,
                sigma = sigma,
                ystar_rmse_train = ystar_rmse_train,
                ystar_cov_train = ystar_cov_train,
                ystar_rmse_test = ystar_rmse_test,
                ystar_cov_test = ystar_cov_test,
                train_timing = train_timing,
                timing = timing)

if(m %in% c("vcbart", "btvcm", "tvc", "ks", "lm")){
  beta_mse_train <- colMeans( (beta_train - fit$train$beta[,"MEAN",])^2 )
  beta_cov_train <- 
    colMeans( (beta_train >= fit$train$beta[,"L95",] &
                 beta_train <= fit$train$beta[,"U95",]))
  
  beta_mse_test <- colMeans( (beta_test - fit$test$beta[,"MEAN",])^2 )
  beta_cov_test <- colMeans( (beta_test >= fit$test$beta[,"L95",] &
                                beta_test <= fit$test$beta[,"U95",]))
  
  results[["beta_mse_train"]] <- beta_mse_train
  results[["beta_cov_train"]] <- beta_cov_train
  results[["beta_mse_test"]] <- beta_mse_test
  results[["beta_cov_test"]] <- beta_cov_test
}

#assign(paste0("p5R20_results_", job_id),
#       results)
#save(list = paste0("p5R20_results_", job_id),
#     file = paste0("p5R20_results_", job_id, ".RData"))



