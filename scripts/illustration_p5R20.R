source("scripts/vcbart_wrapper.R")
load("data/p5R20_data.RData")

vcbart_adapt_subset <- 
  vcbart_wrapper(Y_all, X_all, Z_all, n_all, 
                 X_all, Z_all, n_all, cutpoints, error_structure = "ind",
                 split_probs_type = "adaptive", burn = 250, nd = 1000,
                 verbose = TRUE, print_every = 250)
Z_plot <- Z_all
beta_plot <- beta_all
beta0_hat <- vcbart_adapt_subset$train$beta[,,1]
beta1_hat <- vcbart_adapt_subset$train$beta[,,2]
beta2_hat <- vcbart_adapt_subset$train$beta[,,3]
beta3_hat <- vcbart_adapt_subset$train$beta[,,4]

save(Z_plot, beta_plot, beta0_hat, beta1_hat, beta2_hat, beta3_hat, file = "results/illustration_p5R20.RData")