source("scripts/vcbart_wrapper.R")

args <- commandArgs(TRUE)
sim_number <- as.numeric(args[1])
load("data/p5R20_data.RData")
rm(X_all, Y_all, Z_all, start_index_all, end_index_all, n_all, N_all)
# We can go ahead and re-load the data used in the simulation from the main text
load(paste0("data/sim_p5R20/data_p5R20_", sim_number, "_sigma1.RData"))

tau_seq <- c(1, 1/4, 1/2, 2/3, 3/2,2, 4)
tau_names <- paste0("tau_", c("1_1", "1_4", "1_2", "2_3", "3_2", "2_1", "4_1"))


for(ix in 1:length(tau_seq)){
  tau <- tau_seq[ix]
  print(paste("Starting tau =", tau, "at", Sys.time()))
  assign(paste0("vcbart_", tau_names[ix], "_", sim_number),
         vcbart_wrapper(Y_train, X_train, Z_train, n_train,
                        X_test, Z_test, n_test, cutpoints,
                        error_structure = "ind", split_probs_type = "adaptive", tau = tau,
                        burn = 500, nd = 1000, verbose = TRUE, print_every = 250))
}

save_list <- paste0("vcbart_", tau_names, "_", sim_number)
save(list = save_list, file = paste0("results/sim_p5R20_tau/results_tau_", sim_number, ".RData"))
