source("scripts/vcbart_cs_wrapper.R")
args <- commandArgs(TRUE)
sim_number <- as.numeric(args[1])
rho_eps <- as.numeric(args[2])/100


load(paste0("data/HRS/HRS_data_", sim_number, ".RData"))

assign(paste0("vcbart_fixed_cs", rho_eps*100, "_", sim_number), 
       vcbart_cs_wrapper(Y_train, X_train, Z_train, n_vec_train, 
                         X_test, Z_test, n_vec_test, cutpoints, error_structure = "cs", rho_eps = rho_eps,
                         split_probs_type = "fixed", burn = 500, nd = 1000,
                         verbose = TRUE, print_every = 100))
save(list = paste0("vcbart_adapt_", sim_number), 
     file = paste0("results/HRS/vcbart_fixed/vcbart_fixed_cs", rho_eps*100,"_", sim_number, ".RData"))

