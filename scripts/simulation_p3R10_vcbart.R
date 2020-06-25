# p3R10 simulation for fixed vcbart

source("scripts/vcbart_wrapper.R")
load("data/p3R10_data.RData")


args <- commandArgs(TRUE)
sim_number <- as.numeric(args[1])

set.seed(129 + sim_number)

train_index <- sort(sample(1:N, size = floor(0.75 * N), replace = FALSE))
test_index <- (1:N)[which(! (1:N) %in% train_index)]

X_train <- as.matrix(X_all[train_index,], nrow = length(train_index), ncol = p)
Y_train <- Y_all[train_index]
Z_train <- as.matrix(Z_all[train_index,], ncol = R)
n_train <- nrow(X_train)
start_index_train <- 1

X_test <- as.matrix(X_all[test_index,], nrow = length(test_index), ncol = p)
Y_test <- Y_all[test_index]
Z_test <- as.matrix(Z_all[test_index,], ncol = R)
n_test <- nrow(X_test)
start_index_test <- 1

assign(paste0("vcbart_fixed_", sim_number), 
       vcbart_wrapper(Y_train, X_train, Z_train, n_train, 
                      X_test, Z_test, n_test, cutpoints, error_structure = "ind",
                      split_probs_type = "fixed", burn = 500, nd = 1000,
                      verbose = FALSE, print_every = 100))

assign(paste0("vcbart_adapt_", sim_number),
       vcbart_wrapper(Y_train, X_train, Z_train, n_train, 
                      X_test, Z_test, n_test, cutpoints, error_structure = "ind",
                      split_probs_type = "adaptive", burn = 500, nd = 1000,
                      verbose = FALSE, print_every = 100))

save(list = paste0("vcbart_fixed_", sim_number), 
     file = paste0("results/p3R10/vcbart_fixed/vcbart_fixed_", sim_number, ".RData") )

save(list = paste0("adapt_fixed_", sim_number), 
     file = paste0("results/p3R10/vcbart_adapt/vcbart_adapt_", sim_number, ".RData") )



