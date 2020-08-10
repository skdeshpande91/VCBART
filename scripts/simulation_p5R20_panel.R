# p5R20 panel simulations

#library(Rcpp)
#library(RcppArmadillo)

source("scripts/vcbart_wrapper.R")
source("scripts/vcbart_cs_wrapper.R")
load("data/p5R20_rho75_data.RData")


args <- commandArgs(TRUE)
sim_number <- as.numeric(args[1])
set.seed(129 + sim_number)


# Randomly select 75 people for training and 25 people for testing
tmp_index <- sample(1:n, size = 100, replace = FALSE)
train_id <- sort(tmp_index[1:75])
test_id <- sort(tmp_index[76:100])

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
     file = paste0("data/sim_p5R20_rho75/data_p5R20_rho75_", sim_number, ".RData"))

rho_list <- c(rho_true, seq(0.0, 0.9, by = 0.1))

for(ix in 1:length(rho_list)){
  rho <- rho_list[ix]
  print(paste("Starting rho = ", rho, "at", Sys.time()))
  
  if(rho == 0){
    assign(paste0("vcbart_adapt_cs", rho*100, "_", sim_number),
           vcbart_wrapper(Y_train, X_train, Z_train, n_train,
                          X_test, Z_test, n_test, cutpoints,
                          error_structure = "ind", split_probs_type = "adaptive",
                          burn = 500, nd = 1000, verbose = TRUE, print_every = 250))
  } else{
    assign(paste0("vcbart_adapt_cs", rho*100, "_", sim_number),
           vcbart_cs_wrapper(Y_train, X_train, Z_train, n_train,
                             X_test, Z_test, n_test, cutpoints,
                             error_structure = "cs", rho_eps = rho, split_probs_type = "adaptive",
                             burn = 500, nd = 1000, verbose = TRUE, print_every = 250))
  }
}

save_list <- paste0("vcbart_adapt_cs", 100*rho_list, "_", sim_number)

save(list = save_list, file = paste0("results/sim_p5R20_rho75/results_sim", sim_number, ".RData"))
