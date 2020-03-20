
library(MASS)
library(VCBART)

load("data/p3R2_panel_data.RData")

# Script was designed to run on a high-performance computing cluster.
# Arguments were passed in via command line.
# To run on personal machine or interactive mode, please comment out the following two lines
# and set rep manually (we used rep = 1, 2, ..., 50 in the main paper)

args <- commandArgs(TRUE)
rep <- as.numeric(args[1])

Sigma <- 0.75^abs(outer(1:(ni_train + ni_test), 1:(ni_train + ni_test), FUN = "-"))
rho_true <- 0.75

eps_train <- c()
eps_test <- c()
set.seed(219 + rep)
for(i in 1:n){
  tmp_eps <- mvrnorm(n = 1, mu = rep(0, times = ni_train + ni_test), Sigma = Sigma)
  eps_train <- c(eps_train, tmp_eps[1:ni_train])
  eps_test <- c(eps_test, tmp_eps[(1 + ni_train):(ni_test + ni_train)])
}

Y_train <- mu_train + sigma * eps_train
Y_test <- mu_test + sigma * eps_test

methods <- c("true", "ind", paste0("ar_rho", seq(1, 9, by = 1)), paste0("cs_rho", seq(1, 9, by = 1)))
beta_mse_train <- matrix(nrow = length(methods), ncol = p+1, dimnames = list(methods, c()))
beta_mse_test <- beta_mse_train

beta_coverage_train <- beta_mse_train
beta_coverage_test <- beta_mse_test

beta_int_train <- beta_mse_train
beta_int_test <- beta_mse_train

ystar_mse_train <- rep(NA, times = length(methods))
names(ystar_mse_train) <- methods

ystar_mse_test <- ystar_mse_train

ystar_coverage_train <- ystar_mse_train
ystar_coverage_test <- ystar_mse_train

ystar_int_train <- ystar_mse_train
ystar_int_test <- ystar_int_train


sigma_est <- ystar_mse_train # posterior mean of sigma
timing <- ystar_mse_train


#################
# Run well-specified model
#################
X_train <- cbind(rep(1, times = nrow(X_train)), X_train)
X_test <- cbind(rep(1, times = nrow(X_test)), X_test)

true_chain1 <- vc_BART_ar(Y_train, X_train = X_train, Z_train = Z_train, n_vec_train = n_vec_train, start_index_vec_train = start_index_train,
                          X_test = X_test, Z_test = Z_test, n_vec_test = n_vec_test, start_index_vec_test = start_index_test,
                          cutpoints, nd = 1000, burn = 250, rho = 0.75, verbose = TRUE, print_every = 50)

true_chain2 <- vc_BART_ar(Y_train, X_train = X_train, Z_train = Z_train, n_vec_train = n_vec_train, start_index_vec_train = start_index_train,
                          X_test = X_test, Z_test = Z_test, n_vec_test = n_vec_test, start_index_vec_test = start_index_test,
                          cutpoints, nd = 1000, burn = 250, rho = 0.75, verbose = TRUE, print_every = 50)
true_sum <- get_summary(true_chain1, true_chain2)

beta_mse_train["true",] <- colMeans( (beta_train - true_sum$train$beta[,"MEAN",])^2 )
beta_mse_test["true",] <- colMeans( (beta_test- true_sum$test$beta[,"MEAN",])^2 )

beta_coverage_train["true",] <- colMeans( (beta_train >= true_sum$train$beta[,"L95",]) & (beta_train <= true_sum$train$beta[,"U95",]))
beta_coverage_test["true",] <- colMeans( (beta_test >= true_sum$test$beta[,"L95",]) & (beta_test <= true_sum$test$beta[,"U95",]))

beta_int_train["true",] <- 1
beta_int_test["true",] <- 1

ystar_mse_train["true"] <- mean( (Y_train - true_sum$train$ystar[,"MEAN"])^2 )
ystar_mse_test["true"] <- mean( (Y_test - true_sum$test$ystar[,"MEAN"])^2 )

ystar_coverage_train["true"] <- mean( (Y_train >= true_sum$train$ystar[,"L95"]) & (Y_train <= true_sum$train$ystar[,"U95"]))
ystar_coverage_test["true"] <- mean( (Y_test >= true_sum$test$ystar[,"L95"]) & (Y_test <= true_sum$test$ystar[,"U95"]))

ystar_int_train["true"] <- 1
ystar_int_test["true"] <- 1

sigma_est["true"] <- mean(true_sum$sigma)
timing["true"] <- (true_chain1$time + true_chain2$time)/2

##############
# Run with independent errors
##############

ind_chain1 <- vc_BART_ind(Y_train, X_train = X_train, Z_train = Z_train, n_vec_train = n_vec_train, start_index_vec_train = start_index_train,
                          X_test = X_test, Z_test = Z_test, n_vec_test = n_vec_test, start_index_vec_test = start_index_test,
                          cutpoints, nd = 1000, burn = 250, verbose = TRUE, print_every = 50)

ind_chain2 <- vc_BART_ind(Y_train, X_train = X_train, Z_train = Z_train, n_vec_train = n_vec_train, start_index_vec_train = start_index_train,
                          X_test = X_test, Z_test = Z_test, n_vec_test = n_vec_test, start_index_vec_test = start_index_test,
                          cutpoints, nd = 1000, burn = 250, verbose = TRUE, print_every = 50)

ind_sum <- get_summary(ind_chain1, ind_chain2)


beta_mse_train["ind",] <- colMeans( (beta_train - ind_sum$train$beta[,"MEAN",])^2 )
beta_mse_test["ind",] <- colMeans( (beta_test- ind_sum$test$beta[,"MEAN",])^2 )

beta_coverage_train["ind",] <- colMeans( (beta_train >= ind_sum$train$beta[,"L95",]) & (beta_train <= ind_sum$train$beta[,"U95",]))
beta_coverage_test["ind",] <- colMeans( (beta_test >= ind_sum$test$beta[,"L95",]) & (beta_test <= ind_sum$test$beta[,"U95",]))

beta_int_train["ind",] <- colMeans( (ind_sum$train$beta[,"U95",] - ind_sum$train$beta[,"L95",])/(true_sum$train$beta[,"U95",] - true_sum$train$beta[,"L95",]))
beta_int_test["ind",] <- colMeans( (ind_sum$test$beta[,"U95",] - ind_sum$test$beta[,"L95",])/(true_sum$test$beta[,"U95",] - true_sum$test$beta[,"L95",]))

ystar_mse_train["ind"] <- mean( (Y_train - ind_sum$train$ystar[,"MEAN"])^2 )
ystar_mse_test["ind"] <- mean( (Y_test - ind_sum$test$ystar[,"MEAN"])^2 )

ystar_coverage_train["ind"] <- mean( (Y_train >= ind_sum$train$ystar[,"L95"]) & (Y_train <= ind_sum$train$ystar[,"U95"]))
ystar_coverage_test["ind"] <- mean( (Y_test >= ind_sum$test$ystar[,"L95"]) & (Y_test <= ind_sum$test$ystar[,"U95"]))


ystar_int_train["ind"] <- mean( (ind_sum$train$ystar[,"U95"]- ind_sum$train$ystar[,"L95"])/(true_sum$train$ystar[,"U95"] - true_sum$train$ystar[,"L95"]))
ystar_int_test["ind"] <- mean( (ind_sum$test$ystar[,"U95"]- ind_sum$test$ystar[,"L95"])/(true_sum$test$ystar[,"U95"] - true_sum$test$ystar[,"L95"]))

sigma_est["ind"] <- mean(ind_sum$sigma)

timing["ind"] <- (ind_chain1$time + ind_chain2$time)/2

for(rho_ix in 1:9){
  rho <- rho_ix/10
  
  print(paste("Starting rho = ", rho, "at", Sys.time()))
  
  ar_chain1 <- vc_BART_ar(Y_train, X_train = X_train, Z_train = Z_train, n_vec_train = n_vec_train, start_index_vec_train = start_index_train,
                          X_test = X_test, Z_test = Z_test, n_vec_test = n_vec_test, start_index_vec_test = start_index_test,
                          cutpoints, nd = 1000, burn = 250, rho = rho, verbose = TRUE, print_every = 50)
  ar_chain2 <- vc_BART_ar(Y_train, X_train = X_train, Z_train = Z_train, n_vec_train = n_vec_train, start_index_vec_train = start_index_train,
                          X_test = X_test, Z_test = Z_test, n_vec_test = n_vec_test, start_index_vec_test = start_index_test,
                          cutpoints, nd = 1000, burn = 250, rho = rho, verbose = TRUE, print_every = 50)
  
  ar_sum <- get_summary(ar_chain1, ar_chain2)
  
  beta_mse_train[paste0("ar_rho", rho_ix),] <- colMeans( (beta_train - ar_sum$train$beta[,"MEAN",])^2 )
  beta_mse_test[paste0("ar_rho", rho_ix),] <- colMeans( (beta_test- ar_sum$test$beta[,"MEAN",])^2 )
  
  beta_coverage_train[paste0("ar_rho", rho_ix),] <- colMeans( (beta_train >= ar_sum$train$beta[,"L95",]) & (beta_train <= ar_sum$train$beta[,"U95",]))
  beta_coverage_test[paste0("ar_rho", rho_ix),] <- colMeans( (beta_test >= ar_sum$test$beta[,"L95",]) & (beta_test <= ar_sum$test$beta[,"U95",]))
  
  beta_int_train[paste0("ar_rho", rho_ix),] <- colMeans( (ar_sum$train$beta[,"U95",] - ar_sum$train$beta[,"L95",])/(true_sum$train$beta[,"U95",] - true_sum$train$beta[,"L95",]))
  beta_int_test[paste0("ar_rho", rho_ix),] <- colMeans( (ar_sum$test$beta[,"U95",] - ar_sum$test$beta[,"L95",])/(true_sum$test$beta[,"U95",] - true_sum$test$beta[,"L95",]))
  
  ystar_mse_train[paste0("ar_rho", rho_ix)] <- mean( (Y_train - ar_sum$train$ystar[,"MEAN"])^2 )
  ystar_mse_test[paste0("ar_rho", rho_ix)] <- mean( (Y_test - ar_sum$test$ystar[,"MEAN"])^2 )
  
  ystar_coverage_train[paste0("ar_rho", rho_ix)] <- mean( (Y_train >= ar_sum$train$ystar[,"L95"]) & (Y_train <= ar_sum$train$ystar[,"U95"]))
  ystar_coverage_test[paste0("ar_rho", rho_ix)] <- mean( (Y_test >= ar_sum$test$ystar[,"L95"]) & (Y_test <= ar_sum$test$ystar[,"U95"]))
  
  
  ystar_int_train[paste0("ar_rho", rho_ix)] <- mean( (ar_sum$train$ystar[,"U95"]- ar_sum$train$ystar[,"L95"])/(true_sum$train$ystar[,"U95"] - true_sum$train$ystar[,"L95"]))
  ystar_int_test[paste0("ar_rho", rho_ix)] <- mean( (ar_sum$test$ystar[,"U95"]- ar_sum$test$ystar[,"L95"])/(true_sum$test$ystar[,"U95"] - true_sum$test$ystar[,"L95"]))
  
  sigma_est[paste0("ar_rho", rho_ix)] <- mean(ar_sum$sigma)
  timing[paste0("ar_rho", rho_ix)] <- (ar_chain1$time + ar_chain2$time)/2
  
  # Now do the compound symmetry 
  cs_chain1 <- vc_BART_cs(Y_train, X_train = X_train, Z_train = Z_train, n_vec_train = n_vec_train, start_index_vec_train = start_index_train,
                          X_test = X_test, Z_test = Z_test, n_vec_test = n_vec_test, start_index_vec_test = start_index_test,
                          cutpoints, nd = 1000, burn = 250, rho = rho, verbose = TRUE, print_every = 50)
  cs_chain2 <-  vc_BART_cs(Y_train, X_train = X_train, Z_train = Z_train, n_vec_train = n_vec_train, start_index_vec_train = start_index_train,
                           X_test = X_test, Z_test = Z_test, n_vec_test = n_vec_test, start_index_vec_test = start_index_test,
                           cutpoints, nd = 1000, burn = 250, rho = rho, verbose = TRUE, print_every = 50)
  
  cs_sum <- get_summary(cs_chain1, cs_chain2)
  
  beta_mse_train[paste0("cs_rho", rho_ix),] <- colMeans( (beta_train - cs_sum$train$beta[,"MEAN",])^2 )
  beta_mse_test[paste0("cs_rho", rho_ix),] <- colMeans( (beta_test- cs_sum$test$beta[,"MEAN",])^2 )
  
  beta_coverage_train[paste0("cs_rho", rho_ix),] <- colMeans( (beta_train >= cs_sum$train$beta[,"L95",]) & (beta_train <= cs_sum$train$beta[,"U95",]))
  beta_coverage_test[paste0("cs_rho", rho_ix),] <- colMeans( (beta_test >= cs_sum$test$beta[,"L95",]) & (beta_test <= cs_sum$test$beta[,"U95",]))
  
  beta_int_train[paste0("cs_rho", rho_ix),] <- colMeans( (cs_sum$train$beta[,"U95",] - cs_sum$train$beta[,"L95",])/(true_sum$train$beta[,"U95",] - true_sum$train$beta[,"L95",]))
  beta_int_test[paste0("cs_rho", rho_ix),] <- colMeans( (cs_sum$test$beta[,"U95",] - cs_sum$test$beta[,"L95",])/(true_sum$test$beta[,"U95",] - true_sum$test$beta[,"L95",]))
  
  ystar_mse_train[paste0("cs_rho", rho_ix)] <- mean( (Y_train - cs_sum$train$ystar[,"MEAN"])^2 )
  ystar_mse_test[paste0("cs_rho", rho_ix)] <- mean( (Y_test - cs_sum$test$ystar[,"MEAN"])^2 )
  
  ystar_coverage_train[paste0("cs_rho", rho_ix)] <- mean( (Y_train >= cs_sum$train$ystar[,"L95"]) & (Y_train <= cs_sum$train$ystar[,"U95"]))
  ystar_coverage_test[paste0("cs_rho", rho_ix)] <- mean( (Y_test >= cs_sum$test$ystar[,"L95"]) & (Y_test <= cs_sum$test$ystar[,"U95"]))
  
  
  ystar_int_train[paste0("cs_rho", rho_ix)] <- mean( (cs_sum$train$ystar[,"U95"]- cs_sum$train$ystar[,"L95"])/(true_sum$train$ystar[,"U95"] - true_sum$train$ystar[,"L95"]))
  ystar_int_test[paste0("cs_rho", rho_ix)] <- mean( (cs_sum$test$ystar[,"U95"]- cs_sum$test$ystar[,"L95"])/(true_sum$test$ystar[,"U95"] - true_sum$test$ystar[,"L95"]))
  
  sigma_est[paste0("cs_rho", rho_ix)] <- mean(cs_sum$sigma)
  
  timing[paste0("cs_rho", rho_ix)] <- (cs_chain1$time + cs_chain2$time)/2
}

assign(paste0("beta_mse_train_", sim_type, "_", rep), beta_mse_train)
assign(paste0("beta_coverage_train_", sim_type, "_", rep), beta_coverage_train)
assign(paste0("beta_int_train_", sim_type, "_", rep), beta_int_train)

assign(paste0("beta_mse_test_", sim_type, "_", rep), beta_mse_test)
assign(paste0("beta_coverage_test_", sim_type, "_", rep), beta_coverage_test)
assign(paste0("beta_int_test_", sim_type, "_", rep), beta_int_test)

assign(paste0("ystar_mse_train_", sim_type, "_", rep), ystar_mse_train)
assign(paste0("ystar_coverage_train_", sim_type, "_", rep), ystar_coverage_train)
assign(paste0("ystar_int_train_", sim_type, "_", rep), ystar_int_train)

assign(paste0("ystar_mse_test_", sim_type, "_", rep), ystar_mse_test)
assign(paste0("ystar_coverage_test_", sim_type, "_", rep), ystar_coverage_test)
assign(paste0("ystar_int_test_", sim_type, "_", rep), ystar_int_test)

assign(paste0("timing_", sim_type, "_", rep), timing)
assign(paste0("sigma_est_", sim_type, "_", rep), sigma_est)


save_list <- paste0(c("beta_mse_train_", "beta_mse_test_",
                      "beta_coverage_train_", "beta_coverage_test_",
                      "beta_int_train_", "beta_int_test_",
                      "ystar_mse_train_", "ystar_mse_test_",
                      "ystar_coverage_train_", "ystar_coverage_test_",
                      "ystar_int_train_", "ystar_int_test_",
                      "timing_", "sigma_est_"), sim_type, "_", rep)
save(list = save_list, file = paste0("results/panel_p3R2_", sim_type, "_", rep, ".RData"))




