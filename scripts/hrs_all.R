# Run all of the HRS data
source("scripts/vcbart_cs_wrapper.R")
load("data/HRS/HRS_all.RData")
hrs_vcbart_adapt_rho <- 
  vcbart_cs_adapt_rho_wrapper(Y_all, X_all, Z_all, n_all,
                              X_plot, Z_plot, n_plot, cutpoints,
                              error_structure = "cs", split_probs_type = "adaptive",
                              burn = 500, nd = 1000, verbose = TRUE, print_every = 250)
save(hrs_vcbart_adapt_rho, file = "results/HRS/vcbart_hrs_all.RData")


