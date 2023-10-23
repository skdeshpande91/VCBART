load("philly_crime_data.RData")
load("sim_settings.RData")
source("vcbart_ind_wrapper.R")
source("vcbart_cs_wrapper.R")

#args <- commandArgs(TRUE)
#job_id <- as.numeric(args[1])

#tract <- sim_settings[job_id, "tract"] # heldout tract
#p <- sim_settings[job_id, "p"]
#method <- sim_settings[job_id, "method"]

tract <- 37 # integer b/w 1 and 384, indicating which tract should be heldout
p <- 2 # integer b/w 1 and 4, indicating degree of polynomial approximation
method <- "cs" # either "ind" or "cs" encoding the error structure.

test_index <- which(Z_cat_all[,1] == tract-1) # C++ is 0-index...
train_index <- (1:N)[-test_index]

Y_train <- Y_all[train_index]
Y_test <- Y_all[test_index]

if(p == 1){
  X_train <- matrix(X_all[train_index, 1:p], ncol = p)
  X_test <- matrix(X_all[test_index,1:p], ncol = p)
} else{
  X_train <- X_all[train_index,1:p]
  X_test <- X_all[test_index,1:p]
}

Z_cat_train <- matrix(Z_cat_all[train_index,], ncol = 1)
Z_cat_test <- matrix(Z_cat_all[test_index,], ncol = 1)


ni_train <- ni_all[-tract]
subj_id_train <- rep(1:(n-1), times = ni_train) 

if(method == "ind"){
  fit <- vcbart_ind_wrapper(Y_train = Y_train,
                            subj_id_train = subj_id_train,
                            ni_train = ni_train,
                            X_train = X_train,
                            Z_cat_train = Z_cat_train,
                            X_test = X_test,
                            Z_cat_test = Z_cat_test,
                            unif_cuts = unif_cuts,
                            cutpoints_list = cutpoints_list,
                            cat_levels_list = cat_levels_list,
                            edge_mat_list = edge_mat_list,
                            graph_split = graph_split,
                            sparse = FALSE,
                            M = 50)
} else{
  fit <- vcbart_cs_wrapper(Y_train = Y_train,
                           subj_id_train = subj_id_train,
                           ni_train = ni_train,
                           X_train = X_train,
                           Z_cat_train = Z_cat_train,
                           X_test = X_test,
                           Z_cat_test = Z_cat_test,
                           unif_cuts = unif_cuts,
                           cutpoints_list = cutpoints_list,
                           cat_levels_list = cat_levels_list,
                           edge_mat_list = edge_mat_list,
                           graph_split = graph_split,
                           sparse = FALSE, M = 50)
}

ystar_rmse_train <- sqrt(mean( (Y_train - fit$train$ystar[,"MEAN"])^2 ))
ystar_rmse_test <- sqrt(mean( (Y_test - fit$test$ystar[,"MEAN"])^2 ))

ystar_cov_train <- mean( (Y_train >= fit$train$ystar[,"L95"]) & (Y_train <= fit$train$ystar[,"U95"]))
ystar_cov_test <- mean( (Y_test >= fit$test$ystar[,"L95"]) & (Y_test <= fit$test$ystar[,"U95"]))

results <- 
  list(tract = tract,
       p = p,
       method = method,
       ystar_rmse_train = ystar_rmse_train,
       ystar_rmse_test = ystar_rmse_test,
       ystar_cov_train = ystar_cov_train,
       ystar_cov_test = ystar_cov_test,
       baseline_mse = mean( (Y_test - mean(Y_train))^2 ),
       timing = fit$time)



