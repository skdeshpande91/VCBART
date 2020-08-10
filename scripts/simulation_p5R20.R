# Do p5 just like in the original paper
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

load("data/p5R20_data.RData")


args <- commandArgs(TRUE)
sim_number <- as.numeric(args[1])
set.seed(129 + sim_number)

# Randomly select 75 people for training and 25 people for testing
tmp_index <- sample(1:n, size = 100, replace = FALSE)
train_id <- sort(tmp_index[1:75])
test_id <- sort(tmp_index[76:100])

#tmp_index <- sample(1:500, replace = FALSE)
#train_id <- sort(tmp_index[1:375])
#test_id <- sort(tmp_index[376:500])

train_index <- c()
test_index <- c()

for(id in train_id){
  train_index <- c(train_index, start_index_all[id]:end_index_all[id])
}

for(id in test_id){
  test_index <- c(test_index, start_index_all[id]:end_index_all[id])
}

X_train <- X_all[train_index,]
Z_train <- Z_all[train_index,]
Y_train <- Y_all[train_index]
n_train <- n_all[train_id]
cov_train <- cov_all[train_index,]
mod_train <- mod_all[train_index,]
beta_train <- beta_all[train_index,]


X_test <- X_all[test_index,]
Z_test <- Z_all[test_index,]
Y_test <- Y_all[test_index]
n_test <- n_all[test_id]
cov_test <- cov_all[test_index,]
mod_test <- mod_all[test_index,]
beta_test <- beta_all[test_index,]

# Save the data
save(X_train, Z_train, Y_train, n_train, beta_train,
     X_test, Z_test, Y_test, n_test, beta_test,
     file = paste0("data/sim_p5R20/data_p5R20_", sim_number, ".RData"))


###############
# Run VC-BART
##############

assign(paste0("vcbart_fixed_", sim_number), 
       vcbart_wrapper(Y_train, X_train, Z_train, n_train, 
                      X_test, Z_test, n_test, cutpoints, error_structure = "ind",
                      split_probs_type = "fixed", burn = 500, nd = 1000,
                      verbose = TRUE, print_every = 250))

save(list = paste0("vcbart_fixed_", sim_number),
     file = paste0("results/sim_p5R20/vcbart_fixed/vcbart_fixed_", sim_number, ".RData"))

print(paste("Finished VCBART w/ fixed split probs at", Sys.time()))

assign(paste0("vcbart_adapt_", sim_number),
       vcbart_wrapper(Y_train, X_train, Z_train, n_train, 
                      X_test, Z_test, n_test, cutpoints, error_structure = "ind",
                      split_probs_type = "adaptive", burn = 500, nd = 1000,
                      verbose = TRUE, print_every = 250))
save(list = paste0("vcbart_adapt_", sim_number),
     file = paste0("results/sim_p5R20/vcbart_adapt/vcbart_adapt_", sim_number, ".RData"))

print(paste("Finished VCBART w/ adaptive split probs at", Sys.time()))


####################
# Fit linear model
####################

assign(paste0("lm_", sim_number),
       lm_wrapper(Y_train, X_train, Z_train, X_test, Z_test))
save(list = paste0("lm_", sim_number),
     file = paste0("results/sim_p5R20/lm/lm_", sim_number, ".RData"))
print(paste("Finished lm at", Sys.time()))

####################
# Run kernel smoothing
###################

assign(paste0("kernel_smoothing_",sim_number),
       kernel_smoothing_wrapper(Y_train, cov_train, mod_train, cov_test, mod_test, B = 50))
save(list = paste0("kernel_smoothing_", sim_number),
     file = paste0("results/sim_p5R20/kernel_smoothing/kernel_smoothing_", sim_number, ".RData"))
print(paste("Finished kernel smoothing at", Sys.time()))

#####################
# Run TVC
####################
assign(paste0("tvc_", sim_number),
       tvc_wrapper(Y_train, X_train, Z_train, X_test, Z_test, B = 50))
save(list = paste0("tvc_", sim_number), 
     file = paste0("results/sim_p5R20/tvc/tvc_", sim_number, ".RData"))
print(paste("Finished tvc at", Sys.time()))

##################
# Run boosted tvcm
##################
assign(paste0("boosted_tvcm_", sim_number),
       boosted_tvcm_wrapper(Y_train, X_train, Z_train, X_test, Z_test, B = 50))
save(list = paste0("boosted_tvcm_", sim_number),
     file = paste0("results/sim_p5R20/boosted_tvcm/boosted_tvcm_", sim_number, ".RData"))
print(paste("Finished boosted tvcm at", Sys.time()))

####################
# Run BART
####################
assign(paste0("bart_", sim_number),
       bart_wrapper(Y_train, X_train, Z_train, X_test, Z_test))
save(list = paste0("bart_", sim_number),
     file = paste0("results/sim_p5R20/bart/bart_", sim_number, ".RData"))
print(paste("Finished BART at", Sys.time()))

###################
# Run GBM
###################
assign(paste0("gbm_", sim_number),
       gbm_wrapper(Y_train, X_train, Z_train, X_test, Z_test))
save(list = paste0("gbm_", sim_number),
     file = paste0("results/sim_p5R20/gbm/gbm_", sim_number, ".RData"))
print(paste("Finished GBM at", Sys.time()))

###################
# Run ExtraTrees
###################

assign(paste0("extraTrees_", sim_number),
       extraTrees_wrapper(Y_train, X_train, Z_train, X_test, Z_test))
save(list = paste0("extraTrees_", sim_number),
     file = paste0("results/sim_p5R20/extraTrees/extraTrees_", sim_number, ".RData"))
print(paste("Finished extraTrees at", Sys.time()))


