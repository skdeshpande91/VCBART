
source("scripts/boosted_tvcm_wrapper.R")
args <- commandArgs(TRUE)
sim_number <- as.numeric(args[1])


load(paste0("data/HRS/HRS_data_", sim_number, ".RData"))

assign(paste0("boosted_tvcm_", sim_number),
       boosted_tvcm_wrapper(Y_train, X_train, Z_train, X_test, Z_test, B = 0))
save(list = paste0("boosted_tvcm_", sim_number),
     file = paste0("results/HRS/boosted_tvcm/boosted_tvcm_", sim_number, ".RData"))
print(paste("Finished boosted tvcm at", Sys.time()))