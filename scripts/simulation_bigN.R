source("scripts/vcbart_wrapper.R")

load("data/bigN_data.RData")

args <- commandArgs(TRUE)
n_train <- as.numeric(args[1])
sim_number <- as.numeric(args[2])
n_test <- floor(n_train/10)

set.seed(129 + sim_number)


# Create the training/testing split
tmp_index <- sample(1:N, size = (n_train + n_test), replace = FALSE)
train_index <- sort(tmp_index[1:n_train])
test_index <- sort(tmp_index[(n_train+1):(n_train + n_test)])

X_train <- X_all[train_index,]
Y_train <- Y_all[train_index]
Z_train <- as.matrix(Z_all[train_index,], ncol = R)

X_test <- X_all[test_index,]
Y_test <- Y_all[test_index]
Z_test <- as.matrix(Z_all[test_index,], ncol = R)


beta_train <- beta_all[train_index,]
beta_test <- beta_all[test_index,]



fit <- vcbart_wrapper(Y_train, X_train, Z_train, n_train, 
               X_test, Z_test, n_test, cutpoints, error_structure = "ind",
               split_probs_type = "adaptive", burn = 500, nd = 1000,
               verbose = TRUE, print_every = 250)
beta_mse_train <- colMeans( (beta_train - fit[["train"]][["beta"]][,"MEAN",])^2 )
beta_mse_test <- colMeans( (beta_test - fit[["test"]][["beta"]][,"MEAN",])^2 )
beta_cov_train <- apply((beta_train >= fit[["train"]][["beta"]][,"L95",]) & (beta_train <= fit[["train"]][["beta"]][,"U95",]),
                        FUN = mean, MARGIN = 2, na.rm = TRUE)
beta_cov_test <- apply((beta_train >= fit[["train"]][["beta"]][,"L95",]) & (beta_train <= fit[["train"]][["beta"]][,"U95",]),
                       FUN = mean, MARGIN = 2, na.rm = TRUE)

ystar_mse_train <- mean( (Y_train - fit[["train"]][["ystar"]][,"MEAN"])^2 )
ystar_mse_test <- mean( (Y_test - fit[["test"]][["ystar"]][,"MEAN"])^2 )
ystar_cov_train <- mean( (Y_train >= fit[["train"]][["ystar"]][,"L95"]) & (Y_train <= fit[["train"]][["ystar"]][,"U95"]))
ystar_cov_test <- mean( (Y_test >= fit[["test"]][["ystar"]][,"L95"]) & (Y_test <= fit[["test"]][["ystar"]][,"U95"]))
time <- fit$time
vcbart_time <- fit$vcbart_time

save_list <- c("beta_mse_train", "beta_mse_test", "beta_cov_train", "beta_cov_test", 
                      "ystar_mse_train", "ystar_mse_test", "ystar_cov_train", "ystar_cov_test",
                      "time", "vcbart_time")
save(list = save_list, file = paste0("results/bigN/results_n", n_train, "_sim", sim_number, ".RData"))
