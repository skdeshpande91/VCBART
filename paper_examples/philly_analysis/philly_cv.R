load("philly_crime_data.RData")
load("sim_settings.RData")
source("vcbart_ind_wrapper.R")
source("vcbart_cs_wrapper.R")
source("bart_wrapper.R")


# The cross-validation study was run on a high-throughput cluster
# To run locally, you need to manually set job_id, which corresponsd to rows in the 
# sim_settings data frame,
#args <- commandArgs(TRUE)
#job_id <- as.numeric(args[1])
job_id <- 1 # should be between 1 and 450

p <- sim_settings[job_id, "p"]
method <- sim_settings[job_id, "method"]
fold <- sim_settings[job_id, "fold"]

test_index <- test_indices[[fold]]
train_index <- train_indices[[fold]]

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

# How many observations per tract
subj_id_train <- subj_id_all[train_index]
ni_train <- as.vector(table(subj_id_train))



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
                            sparse = FALSE, verbose = FALSE,
                            M = 50)
} else if(method == "cs"){
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
                           sparse = FALSE, verbose = FALSE, M = 50)
} else if(method =="bart"){
  N_train <- length(Y_train)
  N_test <- nrow(X_test)
  Z_train <- data.frame(factor(Z_cat_train, levels = 0:383))
  Z_test <- data.frame(factor(Z_cat_test, levels = 0:383))
  fit <- bart_wrapper(Y_train = Y_train, 
                      X_train = X_train, Z_train = Z_train,
                      X_test = X_test, Z_test = Z_test)
}

ystar_rmse_train <- sqrt(mean( (Y_train - fit$train$ystar[,"MEAN"])^2 ))
ystar_rmse_test <- sqrt(mean( (Y_test - fit$test$ystar[,"MEAN"])^2 ))

ystar_cov_train <- mean( (Y_train >= fit$train$ystar[,"L95"]) & (Y_train <= fit$train$ystar[,"U95"]))
ystar_cov_test <- mean( (Y_test >= fit$test$ystar[,"L95"]) & (Y_test <= fit$test$ystar[,"U95"]))

assign(paste0("philly_cv_", job_id),
  list(tract = fold,
       p = p,
       method = method,
       ystar_rmse_train = ystar_rmse_train,
       ystar_rmse_test = ystar_rmse_test,
       ystar_cov_train = ystar_cov_train,
       ystar_cov_test = ystar_cov_test,
       baseline_mse = mean( (Y_test - mean(Y_train))^2 ),
       timing = fit$time))

save(list = paste0("philly_cv_", job_id),
     file = paste0("philly_cv_", job_id, ".RData"))

