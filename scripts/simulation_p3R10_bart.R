# p3R10 simulation for bart
source("scripts/bart_wrapper.R")

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

assign(paste0("bart_", sim_number),
       bart_wrapper(Y_train, X_train, Z_train, X_test, Z_test))
save(list = paste0("bart_", sim_number),
     file = paste0("results/sim_p3R10/bart/bart_", sim_number, ".RData"))

