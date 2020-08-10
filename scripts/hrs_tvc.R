source("scripts/tvc_wrapper.R")
args <- commandArgs(TRUE)
sim_number <- as.numeric(args[1])


load(paste0("data/HRS/HRS_data_", sim_number, ".RData"))

assign(paste0("tvc_", sim_number),
       tvc_wrapper(Y_train, X_train, Z_train, X_test, Z_test, B = 0))
save(list = paste0("tvc_", sim_number),
     file = paste0("results/HRS/tvc/tvc_", sim_number, ".RData"))
print(paste("Finished tvc at", Sys.time()))