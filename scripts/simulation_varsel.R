############
# Simulation study to assess ability of VCBART to
# select the modifiers
############

source("scripts/vcbart_wrapper_varsel.R")

args <- commandArgs(TRUE)
sim_number <- as.numeric(args[1])
load("data/p5R20_data.RData")

load(paste0("data/sim_p5R20/data_p5R20_", sim_number, "_sigma2.RData"))

# Run VCBART with fixed, uniform split probabilities
if(!file.exists(paste0("results/variable_selection/vcbart_fixed/vcbart_fixed_", sim_number, ".RData"))){
  assign(paste0("vcbart_fixed_", sim_number),
         vcbart_varsel_wrapper(Y_train, X_train, Z_train, n_train, 
                        X_test, Z_test, n_test, cutpoints, error_structure = "ind",
                        split_probs_type = "fixed", burn = 500, nd = 1000,
                        verbose = TRUE, print_every = 250))
  save(list = paste0("vcbart_fixed_", sim_number), file = paste0("results/variable_selection/vcbart_fixed/vcbart_fixed_", sim_number, ".RData"))
}


# Run VCBART with adaptive split probabilities just like we propose in the paper
if(!file.exists(paste0("results/variable_selection/vcbart_adapt/vcbart_adapt_", sim_number, ".RData"))){
  assign(paste0("vcbart_adapt_", sim_number),
         vcbart_varsel_wrapper(Y_train, X_train, Z_train, n_train, 
                        X_test, Z_test, n_test, cutpoints, error_structure = "ind",
                        split_probs_type = "adaptive", burn = 500, nd = 1000,
                        verbose = TRUE, print_every = 250))
  save(list = paste0("vcbart_adapt_", sim_number), file = paste0("results/variable_selection/vcbart_adapt/vcbart_adapt_", sim_number, ".RData"))
}

# Run VCBART with adaptive split probabilities but fix the dirichlet hyperparameter to alpha_z = 1 (so theta ~ Dirichlet(1/R, ..., 1/R))
if(!file.exists(paste0("results/variable_selection/adapt_dir1/adapt_dir1_", sim_number, ".RData"))){
  assign(paste0("adapt_dir1_", sim_number),
         dirichlet_wrapper(Y_train, X_train, Z_train, n_train, 
                           X_test, Z_test, n_test, cutpoints, init_alpha_z = 1, burn = 500, nd = 1000,
                           verbose = TRUE, print_every = 250))
  save(list = paste0("adapt_dir1_", sim_number), file = paste0("results/variable_selection/adapt_dir1/adapt_dir1_", sim_number, ".RData"))
}


# Run VCBART with adaptive split probabilitie but fix the dirichlet hyperparameter to alpha_z = R (so theta ~ Dirichlet(1, ..., 1))
if(!file.exists(paste0("results/variable_selection/adapt_dirR/adapt_dirR_", sim_number, ".RData"))){
  assign(paste0("adapt_dirR_", sim_number),
         dirichlet_wrapper(Y_train, X_train, Z_train, n_train, 
                           X_test, Z_test, n_test, cutpoints,  init_alpha_z = R, burn = 500, nd = 1000,
                           verbose = TRUE, print_every = 250))
  save(list = paste0("adapt_dirR_", sim_number), file = paste0("results/variable_selection/adapt_dirR/adapt_dirR_", sim_number, ".RData"))
  
}
