source("scripts/vcbart_cs_wrapper.R")
args <- commandArgs(TRUE)
sim_number <- as.numeric(args[1])


load(paste0("data/HRS/HRS_data_", sim_number, ".RData"))

assign(paste0("vcbart_adapt_rho_", sim_number),
       vcbart_cs_adapt_rho_wrapper(Y_train, X_train, Z_train, n_vec_train,
                                   X_test, Z_test, n_vec_test, cutpoints,
                                   error_structure = "cs", split_probs_type = "adaptive",
                                   burn = 500, nd = 1000, verbose = TRUE, print_every = 250))

save(list = paste0("vcbart_adapt_rho_", sim_number), 
     file = paste0("results/HRS/vcbart_adapt/vcbart_adapt_rho_", sim_number, ".RData"))

