# Simulation for wage1 data
# Make sure to run the script "prepare_wage1_data.R" first!

library(np)
library(gbm)
library(ranger)
library(BART)
library(VCBART)
source("scripts/bart_summary.R")
source("scripts/btvcm_wrapper.R")
source("scripts/kernel_smoothing_wrapper.R")

load("data/wage1.RData")

# Script was designed to run on a high-performance computing cluster.
# Arguments were passed in via command line.
# To run on personal machine or interactive mode, please comment out the following two lines
# and set rep manually (we used rep = 1, 2, ..., 50 in the main paper)

args <- commandArgs(TRUE)
rep <- as.numeric(args[1])


set.seed(129 + rep)

## Generate the test/train split
train_index <- sort(sample(1:N, size = floor(0.75 * N), replace = FALSE))
test_index <- (1:N)[which(! (1:N) %in% train_index)]

X_train <- as.matrix(X_all[train_index,], nrow = length(train_index), ncol = p)
Y_train <- Y_all[train_index]
Z_train <- as.matrix(Z_all[train_index,], ncol = R)


X_test <- as.matrix(X_all[test_index,], nrow = length(test_index), ncol = p)
Y_test <- Y_all[test_index]
Z_test <- as.matrix(Z_all[test_index,], ncol = R)

tmp_data_all <- data.frame("Y" = Y_all, cov_all, mod_all)
tmp_data_train <- tmp_data_all[train_index,]
tmp_data_test <- tmp_data_all[test_index,]

cov_train <- as.data.frame(cov_all[train_index,])
mod_train <- mod_all[train_index,]

cov_test <- as.data.frame(cov_all[test_index,])
mod_test <- mod_all[test_index,]

colnames(cov_train) <- "educ"
colnames(cov_test) <- "educ"

methods <- c("train_mean", "lm", "np", "btvcm", "vc_bart", "bart", "rf", "gbm")

ystar_mse_train <- rep(NA, times = length(methods))
names(ystar_mse_train) <- methods
ystar_mse_test <- ystar_mse_train
ystar_coverage_train <- ystar_mse_train
ystar_coverage_test <- ystar_mse_train
ystar_int_train <- ystar_mse_train
ystar_int_test <- ystar_mse_test

# Timing
timing <- ystar_mse_train

#########################
# Run VC-BART
#########################

# Add column for the intercept
X_train <- cbind(rep(1, times = nrow(X_train)), X_train)
X_test <- cbind(rep(1, times = nrow(X_test)), X_test)

vc_chain1 <- vc_BART_ind(Y_train, X_train = X_train, Z_train = Z_train,
                         X_test = X_test, Z_test = Z_test, xinfo_list = cutpoints, burn = 250, nd = 1000, verbose = TRUE, print_every = 50)

vc_chain2 <- vc_BART_ind(Y_train, X_train = X_train, Z_train = Z_train,
                         X_test = X_test, Z_test = Z_test, xinfo_list = cutpoints, burn = 250, nd = 1000, verbose = TRUE, print_every = 50)

vc_sum <- get_summary(vc_chain1, vc_chain2)

ystar_mse_train["vc_bart"] <- mean( (Y_train - vc_sum$train$ystar[,"MEAN"])^2 )
ystar_mse_test["vc_bart"] <- mean( (Y_test - vc_sum$test$ystar[,"MEAN"])^2 )

ystar_coverage_train["vc_bart"] <- mean( (Y_train >= vc_sum$train$ystar[,"L95"]) & (Y_train <= vc_sum$train$ystar[,"U95"]))
ystar_coverage_test["vc_bart"] <- mean( (Y_test >= vc_sum$test$ystar[,"L95"]) & (Y_test <= vc_sum$test$ystar[,"U95"]))

ystar_int_train["vc_bart"] <- 1
ystar_int_test["vc_bart"] <- 1


timing["vc_bart"] <- (vc_chain1$time + vc_chain2$time)/2

###############
# Run least squares
###############

# we actually want to fit the proper Mincer equation
# log(wages) ~ educ + exper + exper^2
#tmp_data_train[,"exper2"] <- tmp_data_train[,"exper"]^2

lm_data_train <- tmp_data_train
lm_data_test <- tmp_data_test

lm_data_train[,"exper2"] <- lm_data_train[,"exper"] * lm_data_train[,"exper"]
lm_data_test[,"exper2"] <- lm_data_test[,"exper"] * lm_data_test[,"exper"]

lm_time <- system.time(lm_fit <- lm(Y ~ educ + exper + exper2, data = lm_data_train))["elapsed"]

lm_pred_train <- predict(lm_fit, newdata = lm_data_train, interval = "prediction")
lm_pred_test <- predict(lm_fit, newdata = lm_data_test, interval = "prediction")
ystar_mse_train["lm"] <- mean( (Y_train - lm_pred_train[,"fit"])^2 )
ystar_mse_test["lm"] <- mean( (Y_test - lm_pred_test[,"fit"])^2 )

ystar_coverage_train["lm"] <- mean( (Y_train >= lm_pred_train[,"lwr"]) & (Y_train <= lm_pred_train[,"upr"]))
ystar_coverage_test["lm"] <- mean( (Y_test >= lm_pred_test[,"lwr"]) & (Y_test <= lm_pred_test[,"upr"]))

ystar_int_train["lm"] <- mean( (lm_pred_train[,"upr"] - lm_pred_train[,"lwr"])/(vc_sum$train$ystar[,"U95"] - vc_sum$train$ystar[,"L95"]) )
ystar_int_test["lm"] <- mean( (lm_pred_test[,"upr"] - lm_pred_test[,"lwr"])/(vc_sum$test$ystar[,"U95"] - vc_sum$test$ystar[,"L95"]) )

timing["lm"] <- lm_time
########################
# Run kernel smoothing
########################
np_fit <- kernel_smoothing(Y_train, cov_train, mod_train, cov_test, mod_test, B = 50)
if(!is.null(np_fit)){
  
  ystar_mse_train["np"] <- mean( (Y_train - np_fit$train$ystar[,"MEAN"])^2 )
  ystar_mse_test["np"] <- mean( (Y_test - np_fit$test$ystar[,"MEAN"])^2 )
  
  ystar_coverage_train["np"] <- mean( (Y_train >= np_fit$train$ystar[,"L95"]) & (Y_train <= np_fit$train$ystar[,"U95"]))
  ystar_coverage_test["np"] <- mean( (Y_test >= np_fit$test$ystar[,"L95"]) & (Y_test <= np_fit$test$ystar[,"U95"]))
  
  ystar_int_train["np"] <- mean( (np_fit$train$ystar[,"U95"] - np_fit$train$ystar[,"L95"]) / (vc_sum$train$ystar[,"U95"] - vc_sum$train$ystar[,"L95"]))
  ystar_int_test["np"] <- mean( (np_fit$test$ystar[,"U95"] - np_fit$test$ystar[,"L95"]) / (vc_sum$test$ystar[,"U95"] - vc_sum$test$ystar[,"L95"]))
  timing["np"] <- np_fit$time
  
}

###############
# Tree-Boosted VC (bTVCM)
##############
btvcm_fit <- tvcm(Y_train, X_train, Z_train, X_test, Z_test, B = 50)
if(!is.null(tvcm_fit)){
  
  ystar_mse_train["btvcm"] <- mean( (Y_train - btvcm_fit$train$ystar[,"MEAN"])^2 )
  ystar_mse_test["btvcm"] <- mean( (Y_test - btvcm_fit$test$ystar[,"MEAN"])^2 )
  
  ystar_coverage_train["btvcm"] <- mean( (Y_train >= btvcm_fit$train$ystar[,"L95"]) & (Y_train <= btvcm_fit$train$ystar[,"U95"]))
  ystar_coverage_test["btvcm"] <- mean( (Y_test >= btvcm_fit$test$ystar[,"L95"]) & (Y_test <= btvcm_fit$test$ystar[,"U95"]))
  
  ystar_int_train["btvcm"] <- mean( (btvcm_fit$train$ystar[,"U95"] - btvcm_fit$train$ystar[,"L95"]) / (vc_sum$train$ystar[,"U95"] - vc_sum$train$ystar[,"L95"]))
  ystar_int_test["btvcm"] <- mean( (btvcm_fit$test$ystar[,"U95"] - btvcm_fit$test$ystar[,"L95"]) / (vc_sum$test$ystar[,"U95"] - vc_sum$test$ystar[,"L95"]))
  timing["btvcm"] <- btvcm_fit$time
  
}


##############
# Run BART
#############

# Remove the intercept from X_train and X_test
bart_chain1 <- wbart(cbind(X_train[,-1], Z_train), Y_train, cbind(X_test[,-1], Z_test), printevery = 5000)
bart_chain2 <- wbart(cbind(X_train[,-1], Z_train), Y_train, cbind(X_test[,-1], Z_test), printevery = 5000)

bart_sum <- bart_summary(bart_chain1, bart_chain2)

ystar_mse_train["bart"] <- mean( (Y_train - bart_sum$train$fit[,"MEAN"])^2 )
ystar_mse_test["bart"] <- mean( (Y_test - bart_sum$test$fit[,"MEAN"])^2 )

ystar_coverage_train["bart"] <- mean( (Y_train >= bart_sum$train$ystar[,"L95"]) & (Y_train <= bart_sum$train$ystar[,"U95"]))
ystar_coverage_test["bart"] <- mean( (Y_test >= bart_sum$test$ystar[,"L95"]) & (Y_test <= bart_sum$test$ystar[,"U95"]))

ystar_int_train["bart"] <- mean( (bart_sum$train$ystar[,"U95"] - bart_sum$train$ystar[,"L95"]) / (vc_sum$train$ystar[,"U95"] - vc_sum$train$ystar[,"L95"]) )
ystar_int_test["bart"] <- mean( (bart_sum$test$ystar[,"U95"] - bart_sum$test$ystar[,"L95"]) / (vc_sum$test$ystar[,"U95"] - vc_sum$test$ystar[,"L95"]) )

timing["bart"] <- 0.5 * (bart_chain1$proc.time["elapsed"] + bart_chain2$proc.time["elapsed"])

###############
# Extremely Randomized Trees
##############

rf_time <- system.time(rf_fit <- 
                         tryCatch(ranger(Y~ ., data = tmp_data_train, write.forest = TRUE, splitrule = "extratrees"),
                                  error = function(e){return(NULL)}))["elapsed"]
if(!is.null(rf_fit)){
  rf_pred_train <- predict(rf_fit, data = tmp_data_train[,-1])$predictions
  rf_pred_test <- predict(rf_fit, data = tmp_data_test[,-1])$predictions
  
  
  rf_rmse <- sqrt( mean( (Y_train - rf_pred_train)^2 ))
  
  ystar_mse_train["rf"] <- mean( (Y_train - rf_pred_train)^2 )
  ystar_mse_test["rf"] <- mean( (Y_test - rf_pred_test)^2 )
  
  ystar_coverage_train["rf"] <- mean( (Y_train >= rf_pred_train - qnorm(0.975) * rf_rmse) & (Y_train <= rf_pred_train + qnorm(0.975) * rf_rmse) )
  ystar_coverage_test["rf"] <- mean( (Y_test >= rf_pred_test - qnorm(0.975) * rf_rmse) & (Y_test <= rf_pred_test + qnorm(0.975) * rf_rmse) )
  
  ystar_int_train["rf"] <- mean( (2 * qnorm(0.975) * rf_rmse) / (vc_sum$train$ystar[,"U95"] - vc_sum$train$ystar[,"L95"]))
  ystar_int_test["rf"] <- mean( (2 * qnorm(0.975) * rf_rmse) / (vc_sum$test$ystar[,"U95"] - vc_sum$test$ystar[,"L95"]))
  timing["rf"] <- rf_time
} 

#############
# Gradient Boosting Machines
############
gbm_time <- system.time(
  gbm_fit <- tryCatch(gbm(Y ~ ., data = tmp_data_train, distribution = "gaussian", interaction.depth = 6, cv.folds = 5, n.trees = 1000),
                      error = function(e){return(NULL)}))["elapsed"]
if(!is.null(gbm_fit)){
  n_tree_opt <- gbm.perf(gbm_fit, method = "cv", plot.it = FALSE)
  gbm_pred_train <- predict(gbm_fit, newdata = tmp_data_train, n.trees = n_tree_opt)
  gbm_pred_test <- predict(gbm_fit, newdata = tmp_data_test, n.trees = n_tree_opt)
  
  
  gbm_rmse <- sqrt(mean( (Y_train - gbm_pred_train)^2 ))
  
  ystar_mse_train["gbm"] <- mean( (Y_train - gbm_pred_train)^2 )
  ystar_mse_test["gbm"] <- mean( (Y_test - gbm_pred_test)^2 )
  
  ystar_coverage_train["gbm"] <- mean( (Y_train >= gbm_pred_train - qnorm(0.975) * gbm_rmse) & (Y_train <= gbm_pred_train + qnorm(0.975) * gbm_rmse) )
  ystar_coverage_test["gbm"] <- mean( (Y_test >= gbm_pred_test - qnorm(0.975) * gbm_rmse) & (Y_test <= gbm_pred_test + qnorm(0.975) * gbm_rmse) )
  
  ystar_int_train["gbm"] <- mean( (2 * qnorm(0.975) * gbm_rmse) / (vc_sum$train$ystar[,"U95"] - vc_sum$train$ystar[,"L95"]))
  ystar_int_test["gbm"] <- mean( (2 * qnorm(0.975) * gbm_rmse) / (vc_sum$test$ystar[,"U95"] - vc_sum$test$ystar[,"L95"]))
  timing["gbm"] <- gbm_time
}

########
# One last stupid method: predicting out-of-sample with the mean
rmse_mean <- sqrt(mean((Y_train - mean(Y_train))^2))
ystar_mse_train["train_mean"] <- mean( (Y_train - mean(Y_train))^2 )
ystar_mse_test["train_mean"] <- mean( (Y_test - mean(Y_train))^2 )

ystar_coverage_train["train_mean"] <- mean( (Y_train >= mean(Y_train) - qnorm(0.975) * rmse_mean) & (Y_train <= mean(Y_train) + qnorm(0.975) * rmse_mean) )
ystar_coverage_test["train_mean"] <- mean( (Y_test >= mean(Y_train) - qnorm(0.975) * rmse_mean) & (Y_test <= mean(Y_train) + qnorm(0.975) * rmse_mean) )

ystar_int_train["train_mean"] <- mean( 2 * qnorm(0.975) * rmse_mean / (vc_sum$train$ystar[,"U95"] - vc_sum$train$ystar[,"L95"]) )
ystar_int_test["train_mean"] <- mean( 2 * qnorm(0.975) * rmse_mean / (vc_sum$test$ystar[,"U95"] - vc_sum$test$ystar[,"L95"]) )

##########

assign(paste0("ystar_mse_train_", rep), ystar_mse_train)
assign(paste0("ystar_mse_test_", rep), ystar_mse_test)
assign(paste0("ystar_coverage_train_", rep), ystar_coverage_train)
assign(paste0("ystar_coverage_test_", rep), ystar_coverage_test)
assign(paste0("ystar_int_train_", rep), ystar_int_train)
assign(paste0("ystar_int_test_", rep), ystar_int_test)
assign(paste0("timing_", rep), timing)


save_list <- paste0(c("ystar_mse_train_", "ystar_coverage_train_", "ystar_int_train_", 
                      "ystar_mse_test_", "ystar_coverage_test_", "ystar_int_test_",
                      "timing_"), rep)

save(list = save_list, file = paste0("results/wage1_", rep, ".RData"))
