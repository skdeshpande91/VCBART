
source("scripts/bart_wrapper.R")
args <- commandArgs(TRUE)
sim_number <- as.numeric(args[1])


load(paste0("data/HRS/HRS_data_", sim_number, ".RData"))

assign(paste0("bart_", sim_number),
       bart_wrapper(Y_train, X_train, Z_train, X_test, Z_test))
save(list = paste0("bart_", sim_number),
     file = paste0("results/HRS/bart/bart_", sim_number, ".RData"))