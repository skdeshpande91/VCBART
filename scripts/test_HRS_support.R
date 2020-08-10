# Run VC-BART with CS errors (rho = 0.5) with adaptive split probabilities
# on the entire HRS dataset
library(VCBARTdev)
load("data/HRS/HRS_data_1.RData")

chain1 <- VCBART(Y_train = Y_train,
                 X_train = X_train, Z_train = Z_train, n_train = n_vec_train,
                 X_test = X_test, Z_test = Z_test, n_test = n_vec_test, 
                 cutpoints = cutpoints, 
                 error_structure = "cs",
                 rho_eps = 0.5,
                 split_probs_type = "adapt",
                 burn = 500, nd = 1000, verbose = TRUE, print_every = 250)
chain2 <- VCBART(Y_train = Y_train,
                 X_train = X_train, Z_train = Z_train, n_train = n_train,
                 X_test = X_test, Z_test = Z_test, n_test = n_test, 
                 cutpoints = cutpoints, 
                 error_structure = "cs",
                 rho_eps = 0.5,
                 split_probs_type = "adapt",
                 burn = 500, nd = 1000, verbose = TRUE, print_every = 250)

beta_summary <- summarize_beta(chain1, chain2, burn = 500)
ystar_summary <- summarize_posterior_predictive(chain1, chain2, burn = 500)
beta_support <- get_beta_support(chain1, chain2, burn = 500, max_cutoff = 10)


save(beta_support, "results/HRS/vcbart_adapt_cs50_support_1.RData")