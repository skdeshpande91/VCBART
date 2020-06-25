
source("scripts/gbm_wrapper.R")
args <- commandArgs(TRUE)
sim_number <- as.numeric(args[1])


load(paste0("data/HRS/HRS_data_", sim_number, ".RData"))
assign(paste0("gbm_", sim_number),
       gbm_wrapper(Y_train, X_train, Z_train, X_test, Z_test))
save(list = paste0("gbm_", sim_number),
     file = paste0("results/HRS/gbm/gbm_", sim_number, ".RData"))