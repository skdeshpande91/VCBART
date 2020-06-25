
source("scripts/vcbart_wrapper.R")
args <- commandArgs(TRUE)
sim_number <- as.numeric(args[1])


load(paste0("data/HRS/HRS_data_", sim_number, ".RData"))

assign(paste0("vcbart_fixed_", sim_number), 
       vcbart_wrapper(Y_train, X_train, Z_train, n_vec_train, 
                      X_test, Z_test, n_vec_test, cutpoints, error_structure = "ind",
                      split_probs_type = "fixed", burn = 500, nd = 1000,
                      verbose = TRUE, print_every = 100))
save(list = paste0("vcbart_fixed_", sim_number), 
     file = paste0("results/HRS/vcbart_fixed/vcbart_fixed_", sim_number, ".RData"))

