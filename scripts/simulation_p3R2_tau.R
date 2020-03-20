# Re-run the p3R2 simulation but with different values of tau

library(VCBART)
load("data/p3R2_data.RData")

args <- commandArgs(TRUE)
rep <- as.numeric(args[1])
set.seed(129 + rep)

## Generate the test/train split
train_index <- sort(sample(1:N, size = floor(0.75 * N), replace = FALSE))
test_index <- (1:N)[which(! (1:N) %in% train_index)]

X_train <- X_all[train_index,]
Y_train <- Y_all[train_index]
Z_train <- as.matrix(Z_all[train_index,], ncol = R)
n_vec_train <- length(train_index)
start_index_train <- 1

X_test <- X_all[test_index,]
Y_test <- Y_all[test_index]
Z_test <- as.matrix(Z_all[test_index,], ncol = R)
n_vec_test <- length(test_index)
start_index_test <- 1

#########
# Collect the true betas together in a matrix
beta_all <- cbind(beta0_all, beta1_all, beta2_all, beta3_all)
beta_train <- beta_all[train_index,]
beta_test <- beta_all[test_index,]
#########

tau_seq <- c(1, 1/4, 1/2, 2/3, 3/2,2, 4)
tau_names <- paste0("tau_", c("1_1", "1_4", "1_2", "2_3", "3_2", "2_1", "4_1"))


beta_mse_train <- matrix(nrow = length(tau_names), ncol = p + 1, dimnames = list(tau_names, c()))
beta_mse_test <- matrix(nrow = length(tau_names), ncol = p + 1, dimnames = list(tau_names, c()))
beta_coverage_train <- matrix(nrow = length(tau_names), ncol = p + 1, dimnames = list(tau_names, c()))
beta_coverage_test <- matrix(nrow = length(tau_names), ncol = p + 1, dimnames = list(tau_names, c()))
beta_int_train <- matrix(nrow = length(tau_names), ncol = p + 1, dimnames = list(tau_names, c()))
beta_int_test <- matrix(nrow = length(tau_names), ncol = p + 1, dimnames = list(tau_names, c()))
# For predictions
ystar_mse_train <- rep(NA, times = length(tau_names))
names(ystar_mse_train) <- tau_names
ystar_mse_test <- ystar_mse_train
ystar_coverage_train <- ystar_mse_train
ystar_coverage_test <- ystar_mse_train
ystar_int_train <- ystar_mse_train
ystar_int_test <- ystar_mse_test

timing <- ystar_mse_train
#########
# Run the baseline method
#########

tau_ix <- 1
vc_chain1 <- vc_BART_ind(Y_train, X_train = X_train, Z_train = Z_train,
                         X_test = X_test, Z_test = Z_test, xinfo_list = cutpoints, burn = 250, nd = 1000, verbose = TRUE, print_every = 50)
vc_chain2 <- vc_BART_ind(Y_train, X_train = X_train, Z_train = Z_train,
                         X_test = X_test, Z_test = Z_test, xinfo_list = cutpoints, burn = 250, nd = 1000, verbose = TRUE, print_every = 50)
vc_sum <- get_summary(vc_chain1, vc_chain2)

beta_mse_train[tau_names[tau_ix],] <- colMeans( (beta_train - vc_sum$train$beta[,"MEAN",])^2 , na.rm = TRUE)
beta_mse_test[tau_names[tau_ix],] <- colMeans( (beta_test - vc_sum$test$beta[,"MEAN",])^2 , na.rm = TRUE)

beta_coverage_train[tau_names[tau_ix],] <- apply( (beta_train >= vc_sum$train$beta[,"L95",]) & (beta_train <= vc_sum$train$beta[,"U95",]), FUN = mean, MARGIN = 2)
beta_coverage_test[tau_names[tau_ix],] <- apply( (beta_test >= vc_sum$test$beta[,"L95",]) & (beta_test <= vc_sum$test$beta[,"U95",]), FUN = mean, MARGIN = 2)

beta_int_train[tau_names[tau_ix],] <- 1
beta_int_test[tau_names[tau_ix],] <- 1

ystar_mse_train[tau_names[tau_ix]] <- mean( (Y_train - vc_sum$train$ystar[,"MEAN"])^2 )
ystar_mse_test[tau_names[tau_ix]] <- mean( (Y_test - vc_sum$test$ystar[,"MEAN"])^2 )

ystar_coverage_train[tau_names[tau_ix]] <- mean( (Y_train >= vc_sum$train$ystar[,"L95"]) & (Y_train <= vc_sum$train$ystar[,"U95"]))
ystar_coverage_test[tau_names[tau_ix]] <- mean( (Y_test >= vc_sum$test$ystar[,"L95"]) & (Y_test <= vc_sum$test$ystar[,"U95"]))

ystar_int_train[tau_names[tau_ix]] <- 1
ystar_int_test[tau_names[tau_ix]] <- 1

timing[tau_names[tau_ix]] <- 0.5*(vc_chain1$time + vc_chain2$time) 


for(tau_ix in 2:length(tau_seq)){
  tmp_chain1 <- vc_BART_ind(Y_train, X_train = X_train, Z_train = Z_train,
                            X_test = X_test, Z_test = Z_test, xinfo_list = cutpoints, 
                            sigma_mu_vec = rep(tau_seq[tau_ix]/sqrt(50), times = ncol(X_train)),
                            burn = 250, nd = 1000, verbose = TRUE, print_every = 50)
  tmp_chain2 <-vc_BART_ind(Y_train, X_train = X_train, Z_train = Z_train,
                           X_test = X_test, Z_test = Z_test, xinfo_list = cutpoints, 
                           sigma_mu_vec = rep(tau_seq[tau_ix]/sqrt(50), times = ncol(X_train)),
                           burn = 250, nd = 1000, verbose = TRUE, print_every = 50)
  
  tmp_sum <- get_summary(tmp_chain1, tmp_chain2)
  
  beta_mse_train[tau_names[tau_ix],] <- colMeans( (beta_train - tmp_sum$train$beta[,"MEAN",])^2 , na.rm = TRUE)
  beta_mse_test[tau_names[tau_ix],] <- colMeans( (beta_test - tmp_sum$test$beta[,"MEAN",])^2 , na.rm = TRUE)
  
  beta_coverage_train[tau_names[tau_ix],] <- apply( (beta_train >= tmp_sum$train$beta[,"L95",]) & (beta_train <= tmp_sum$train$beta[,"U95",]), FUN = mean, MARGIN = 2)
  beta_coverage_test[tau_names[tau_ix],] <- apply( (beta_test >= tmp_sum$test$beta[,"L95",]) & (beta_test <= tmp_sum$test$beta[,"U95",]), FUN = mean, MARGIN = 2)
  
  beta_int_train[tau_names[tau_ix],] <- apply( (tmp_sum$train$beta[,"U95",] - tmp_sum$train$beta[,"L95",])/(vc_sum$train$beta[,"U95",] - vc_sum$train$beta[,"L95",]), MARGIN = 2, FUN = mean)
  beta_int_test[tau_names[tau_ix],] <- apply( (tmp_sum$test$beta[,"U95",] - tmp_sum$test$beta[,"L95",])/(vc_sum$test$beta[,"U95",] - vc_sum$test$beta[,"L95",]), MARGIN = 2, FUN = mean)
  
  ystar_mse_train[tau_names[tau_ix]] <- mean( (Y_train - tmp_sum$train$ystar[,"MEAN"])^2 )
  ystar_mse_test[tau_names[tau_ix]] <- mean( (Y_test - tmp_sum$test$ystar[,"MEAN"])^2 )
  
  ystar_coverage_train[tau_names[tau_ix]] <- mean( (Y_train >= tmp_sum$train$ystar[,"L95"]) & (Y_train <= tmp_sum$train$ystar[,"U95"]))
  ystar_coverage_test[tau_names[tau_ix]] <- mean( (Y_test >= tmp_sum$test$ystar[,"L95"]) & (Y_test <= tmp_sum$test$ystar[,"U95"]))
  
  ystar_int_train[tau_names[tau_ix]] <- mean( (tmp_sum$train$ystar[,"U95"] - tmp_sum$train$ystar[,"L95"]) / (vc_sum$train$ystar[,"U95"] - vc_sum$train$ystar[,"L95"]))
  ystar_int_test[tau_names[tau_ix]] <- mean( (tmp_sum$test$ystar[,"U95"] - tmp_sum$test$ystar[,"L95"]) / (vc_sum$test$ystar[,"U95"] - vc_sum$test$ystar[,"L95"]))
  
  timing[tau_names[tau_ix]] <- 0.5 * (tmp_chain1$time + tmp_chain2$time)
}

assign(paste0("beta_mse_train_", rep), beta_mse_train)
assign(paste0("beta_mse_test_", rep), beta_mse_test)
assign(paste0("beta_coverage_train_", rep), beta_coverage_train)
assign(paste0("beta_coverage_test_", rep), beta_coverage_test)
assign(paste0("beta_int_train_", rep), beta_int_train)
assign(paste0("beta_int_test_", rep), beta_int_test)

assign(paste0("ystar_mse_train_", rep), ystar_mse_train)
assign(paste0("ystar_mse_test_", rep), ystar_mse_test)
assign(paste0("ystar_coverage_train_", rep), ystar_coverage_train)
assign(paste0("ystar_coverage_test_", rep), ystar_coverage_test)
assign(paste0("ystar_int_train_", rep), ystar_int_train)
assign(paste0("ystar_int_test_", rep), ystar_int_test)
assign(paste0("timing_", rep), timing)

save_list <- paste0(c("beta_mse_train_", "beta_coverage_train_", "beta_int_train_", 
                      "beta_mse_test_", "beta_coverage_test_", "beta_int_test_",
                      "ystar_mse_train_", "ystar_coverage_train_", "ystar_int_train_", 
                      "ystar_mse_test_", "ystar_coverage_test_", "ystar_int_test_",
                      "timing_"), rep)
save(list = save_list, file = paste0("results/tau_p3R2_", rep, ".RData"))

