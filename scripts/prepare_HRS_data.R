##########
# Experiment 2: Train on 75% of subjects, test on 25%
# Generate 20 training/testing splits
##########
library(tidyverse)
load("data/analysis_data.RData")

for(r in 1:25){
  set.seed(518 + r - 1)
  print(paste0("r = ", r))
  train_index <- sample(1:n_all, size = floor(n_all * 0.75), replace = FALSE)
  test_index <- (1:n_all)[!(1:n_all) %in% train_index]
  
  
  Y_train_list <- Y_list[train_index]
  X_train_list <- X_list[train_index]
  Z_train_list <- Z_list[train_index]
  n_vec_train <- n_vec_all[train_index]
  
  Y_train <- unlist(Y_train_list)
  X_train <- as.matrix(bind_rows(X_train_list))
  Z_train <- as.matrix(bind_rows(Z_train_list))
  
  
  Y_test_list <- Y_list[test_index]
  X_test_list <- X_list[test_index]
  Z_test_list <- Z_list[test_index]
  n_vec_test <- n_vec_all[test_index]
  
  Y_test <- unlist(Y_test_list)
  X_test <- as.matrix(bind_rows(X_test_list))
  Z_test <- as.matrix(bind_rows(Z_test_list))
  
  
  n_train <- length(train_index)
  n_test <- length(test_index)
  
  start_index_train <- 1 + c(0, cumsum(n_vec_train)[-n_train])
  start_index_test <- 1 + c(0, cumsum(n_vec_test)[-n_test])
  
  if(sum(is.na(Z_train)) > 0 || sum(is.na(Z_test)) > 0) stop("NA in Z!")
  if(sum(is.na(X_train)) > 0 || sum(is.na(X_test)) > 0) stop("NA in X!")
  
  
  N_train <- nrow(X_train)
  N_test <- nrow(X_test)
  save(X_train, Y_train, Z_train, n_vec_train, start_index_train, n_train, N_train,
       X_test, Y_test, Z_test, n_vec_test, start_index_test, n_test, N_test,
       p, R, cutpoints, file = paste0("data/HRS_data_", r, ".RData"))
}


