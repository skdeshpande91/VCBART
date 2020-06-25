# Do p3 just like in the original paper
#library(Rcpp)
#library(RcppArmadillo)

source("scripts/vcbart_wrapper.R")
source("scripts/lm_wrapper.R")
source("scripts/kernel_smoothing_wrapper.R")
source("scripts/tvc_wrapper.R")
source("scripts/boosted_tvcm_wrapper.R")
source("scripts/bart_wrapper.R")
source('scripts/extraTrees_wrapper.R')
source("scripts/gbm_wrapper.R")

load("data/p3R10_data.RData")


args <- commandArgs(TRUE)
sim_number <- as.numeric(args[1])
set.seed(129 + sim_number)

#train_index <- sort(sample(1:N, size = floor(0.75 * N), replace = FALSE))
#test_index <- (1:N)[which(! (1:N) %in% train_index)]

index <- sample(1:N, size = 500, replace = FALSE)
train_index <- sort(index[1:375])
test_index <- sort(index[376:500])


X_train <- as.matrix(X_all[train_index,], nrow = length(train_index), ncol = p)
Y_train <- Y_all[train_index]
Z_train <- as.matrix(Z_all[train_index,], ncol = R)
n_train <- nrow(X_train)
start_index_train <- 1
beta_train <- beta_all[train_index,]
cov_train <- cov_all[train_index,]
mod_train <- mod_all[train_index,]


X_test <- as.matrix(X_all[test_index,], nrow = length(test_index), ncol = p)
Y_test <- Y_all[test_index]
Z_test <- as.matrix(Z_all[test_index,], ncol = R)
n_test <- nrow(X_test)
start_index_test <- 1
beta_test <- beta_all[test_index,]
cov_test <- cov_all[test_index,]
mod_test <- mod_all[test_index,]


###############
# Run VC-BART
##############

assign(paste0("vcbart_fixed_", sim_number), 
       vcbart_wrapper(Y_train, X_train, Z_train, n_train, 
                      X_test, Z_test, n_test, cutpoints, error_structure = "ind",
                      split_probs_type = "fixed", burn = 500, nd = 1000,
                      verbose = FALSE, print_every = 100))

save(list = paste0("vcbart_fixed_", sim_number),
     file = paste0("results/sim_p3R10/vcbart_fixed/vcbart_fixed_", sim_number, ".RData"))

print(paste("Finished VCBART w/ fixed split probs at", Sys.time()))

assign(paste0("vcbart_adapt_", sim_number),
       vcbart_wrapper(Y_train, X_train, Z_train, n_train, 
                      X_test, Z_test, n_test, cutpoints, error_structure = "ind",
                      split_probs_type = "adaptive", burn = 500, nd = 1000,
                      verbose = FALSE, print_every = 100))
save(list = paste0("vcbart_adapt_", sim_number),
     file = paste0("results/sim_p3R10/vcbart_adapt/vcbart_adapt_", sim_number, ".RData"))
print(paste("Finished VCBART w/ adaptive split probs at", Sys.time()))


####################
# Fit linear model
####################

assign(paste0("lm_", sim_number),
       lm_wrapper(Y_train, X_train, Z_train, X_test, Z_test))
save(list = paste0("lm_", sim_number),
     file = paste0("results/sim_p3R10/lm/lm_", sim_number, ".RData"))
print(paste("Finished lm at", Sys.time()))

####################
# Run kernel smoothing
###################

assign(paste0("kernel_smoothing_",sim_number),
       kernel_smoothing_wrapper(Y_train, cov_train, mod_train, cov_test, mod_test, B = 50))
save(list = paste0("kernel_smoothing_", sim_number),
     file = paste0("results/sim_p3R10/kernel_smoothing/kernel_smoothing_", sim_number, ".RData"))
print(paste("Finished kernel smoothing at", Sys.time()))

#####################
# Run TVC
####################
assign(paste0("tvc_", sim_number),
       tvc_wrapper(Y_train, X_train, Z_train, X_test, Z_test, B = 50))
save(list = paste0("tvc_", sim_number), 
     file = paste0("results/sim_p3R10/tvc/tvc_", sim_number, ".RData"))
print(paste("Finished tvc at", Sys.time()))

##################
# Run boosted tvcm
##################
assign(paste0("boosted_tvcm_", sim_number),
       boosted_tvcm_wrapper(Y_train, X_train, Z_train, X_test, Z_test, B = 50))
save(list = paste0("boosted_tvcm_", sim_number),
     file = paste0("results/sim_p3R10/boosted_tvcm/boosted_tvcm_", sim_number, ".RData"))
print(paste("Finished boosted tvcm at", Sys.time()))

####################
# Run BART
####################
assign(paste0("bart_", sim_number),
       bart_wrapper(Y_train, X_train, Z_train, X_test, Z_test))
save(list = paste0("bart_", sim_number),
     file = paste0("results/sim_p3R10/bart/bart_", sim_number, ".RData"))
print(paste("Finished BART at", Sys.time()))

###################
# Run GBM
###################
assign(paste0("gbm_", sim_number),
       gbm_wrapper(Y_train, X_train, Z_train, X_test, Z_test))
save(list = paste0("gbm_", sim_number),
     file = paste0("results/sim_p3R10/gbm/gbm_", sim_number, ".RData"))
print(paste("Finished GBM at", Sys.time()))

###################
# Run ExtraTrees
###################

assign(paste0("extraTrees_", sim_number),
       extraTrees_wrapper(Y_train, X_train, Z_train, X_test, Z_test))
save(list = paste0("extraTrees_", sim_number),
     file = paste0("results/sim_p3R10/extraTrees/extraTrees_", sim_number, ".RData"))
print(paste("Finished extraTrees at", Sys.time()))


