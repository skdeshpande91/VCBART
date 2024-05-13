library(coda)
source("assess_support_recovery.R")
source("true_betas_p5R20.R")
source("vcbart_wrapper.R")

nd <- 10000
burn <- nd
n_train <- 250
n_test <- 25
set.seed(41524)
source("generate_data.R")


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
                      nd = nd, burn = burn,
                      n_chains = 4,
                      save_samples = FALSE)

pdf("../../figures/mixing_traceplot20k.pdf", width = 6, height = 6*9/16)
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0))
plot(1, type = "n", xlim = c(0, 20000), ylim = c(0.95, 1.5),
     xlab = "Iterations", ylab = "Sigma", main = "Traceplot (20000 iterations)")
lines(1:20000, fit$sigma[,1], col = rgb(my_rgb[1,2], my_rgb[2,2], my_rgb[3,2], 0.75))
lines(1:20000, fit$sigma[,2], col = rgb(my_rgb[1,3], my_rgb[2,3], my_rgb[3,3], 0.75))
lines(1:20000, fit$sigma[,3], col = rgb(my_rgb[1,4], my_rgb[2,4], my_rgb[3,4], 0.75))
lines(1:20000, fit$sigma[,4], col = rgb(my_rgb[1,7], my_rgb[2,7], my_rgb[3,7], 0.75))
dev.off()
pdf("../../figures/mixing_traceplot2k.pdf", width = 6, height = 6*9/16)
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0))
plot(1, type = "n", xlim = c(0, 2000), ylim = c(0.95, 1.5),
     xlab = "Iterations", ylab = "Sigma", main = "Traceplot (2000 iterations)")
lines(1:2000, fit$sigma[1:2000,1], col = rgb(my_rgb[1,2], my_rgb[2,2], my_rgb[3,2], 0.75))
lines(1:2000, fit$sigma[1:2000,2], col = rgb(my_rgb[1,3], my_rgb[2,3], my_rgb[3,3], 0.75))
lines(1:2000, fit$sigma[1:2000,3], col = rgb(my_rgb[1,4], my_rgb[2,4], my_rgb[3,4], 0.75))
lines(1:2000, fit$sigma[1:2000,4], col = rgb(my_rgb[1,7], my_rgb[2,7], my_rgb[3,7], 0.75))
dev.off()


ystar_rmse_train <- sqrt(mean( (Y_train - fit$train$ystar[,"MEAN"])^2 ))
ystar_cov_train <- mean( (Y_train >= fit$train$ystar[,"L95"] &
                            Y_train <= fit$train$ystar[,"U95"]))

ystar_rmse_test <- sqrt(mean( (Y_test - fit$test$ystar[,"MEAN"])^2 ))
ystar_cov_test <- mean( (Y_test >= fit$test$ystar[,"L95"] &
                           Y_test <= fit$test$ystar[,"U95"]))
train_timing <- fit$train_time
timing <- fit$time

results <- list(sim_number = sim_number,
                nd = nd,
                sigma = sigma,
                ystar_rmse_train = ystar_rmse_train,
                ystar_cov_train = ystar_cov_train,
                ystar_rmse_test = ystar_rmse_test,
                ystar_cov_test = ystar_cov_test,
                train_timing = train_timing,
                timing = timing)

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

# referee suggestion: correlation b/w true beta's and estimated beta's
beta_corr_train <- rep(NA, times = p+1)
beta_corr_test <- rep(NA, times = p+1)
beta_lm_train <- matrix(nrow = 2, ncol = p+1, dimnames = list(c("Slope", "StdError"), c()))
beta_lm_test <- matrix(nrow = 2, ncol = p+1, dimnames = list(c("Slope", "StdError"), c()))

for(j in 1:(p+1)){
  if( (var(beta_train[,j]) != 0) & (var(fit$train$beta[,"MEAN",j]) != 0) ){
    beta_corr_train[j] <- cor(beta_train[,j], fit$train$beta[,"MEAN",j])
    train_lm <- lm(beta_train[,j] ~ fit$train$beta[,"MEAN",j])
    beta_lm_train[, j] <- summary(train_lm)$coefficients[2,1:2]
  }
  if((var(beta_test[,j]) != 0) & (var(fit$test$beta[,"MEAN",j]) != 0)){
    beta_corr_test[j] <- cor(beta_test[,j], fit$test$beta[,"MEAN",j])
    test_lm <- lm(beta_test[,j] ~ fit$test$beta[,"MEAN",j])
    beta_lm_test[,j] <- summary(test_lm)$coefficients[2,1:2]
  }
}
results[["beta_corr_train"]] <- beta_corr_train
results[["beta_lm_train"]] <- beta_lm_train
results[["beta_corr_test"]] <- beta_corr_test
results[["beta_lm_test"]] <- beta_lm_test

tmp_varsel <- assess_support_recovery(fit$varsel$varcounts_means,
                                      fit$varsel$selection_prob,
                                      true_support, R)
results[["varsel"]] <- tmp_varsel$varcount_prob[51,]

sigma_chain <- 
  coda::mcmc.list(
    coda::mcmc(fit$sigma[,1]),
    coda::mcmc(fit$sigma[,2]),
    coda::mcmc(fit$sigma[,3]),
    coda::mcmc(fit$sigma[,4]))
results[["sigma_Rhat"]] <- coda::gelman.diag(sigma_chain)
results[["sigma_neff"]] <- coda::effectiveSize(sigma_chain)

assign(paste0("longer_run_results_", job_id),
       results)
save(list = paste0("longer_run_results_", job_id),
     file = paste0("longer_run_results_", job_id, ".RData"))
