# Run TVC
source("scripts/tvc_wrapper.R")
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


N_train <- nrow(X_train)
N_test <- nrow(X_test)

assign(paste0("tvc_", sim_number),
       tvc_wrapper(Y_train, X_train, Z_train, X_test, Z_test, B = 50))
save(list = paste0("tvc_", sim_number), 
     file = paste0("results/sim_p3R10/tvc/tvc_", sim_number, ".RData"))
print(paste("Finished tvc at", Sys.time()))

