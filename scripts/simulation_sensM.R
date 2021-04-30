# draw data from original p5R20
# vary M from 1, 10, 50, 100
# track beta_mse_test, ystar_mse_test, time
load("data/p5R20_data.RData")
source("scripts/vcbart_wrapper.R")

args <- commandArgs(TRUE)
M <- as.numeric(args[1])
sim_number <- as.numeric(args[2])

load(paste0("data/sim_p5R20/data_p5R20_", sim_number, "_sigma1.RData"))

fit <- vcbart_wrapper(Y_train, X_train, Z_train, n_train, 
                        X_test, Z_test, n_test, cutpoints, M = M, error_structure = "ind",
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
save(list = save_list, file = paste0("results/sensM/results_sensM", M, "_sim", sim_number, ".RData"))