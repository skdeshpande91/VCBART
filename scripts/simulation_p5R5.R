# p5R5 simulations
# Be sure to run the script "generate_p5R5_data.R" before running this script
library(np)
library(gbm)
library(ranger)
library(BART)
library(VCBART)
source("scripts/bart_summary.R")
source("scripts/btvcm_wrapper.R")
source("scripts/kernel_smoothing_wrapper.R")

load("data/p5R5_data.RData")

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

#######
# Compile the data into the data frames necessary for some competiting methods
#######
tmp_data_all <- data.frame("Y" = Y_all, X_all, Z_all)
colnames(tmp_data_all) <- c("Y", "X1", "X2", "X3", "X4", "X5", "Z1", "Z2", "Z3", "Z4", "Z5")

cov_all <- tmp_data_all[, c("X1", "X2", "X3", "X4", "X5")]
mod_all <- tmp_data_all[,c("Z1", "Z2", "Z3", "Z4", "Z5")]

tmp_data_train <- tmp_data_all[train_index,]
tmp_data_test <- tmp_data_all[test_index,]

cov_train <- cov_all[train_index,]
cov_test <- cov_all[test_index,]

mod_train <- mod_all[train_index,]
mod_test <- mod_all[test_index,]

#########
# Collect the true betas together in a matrix
beta_all <- cbind(beta0_all, beta1_all, beta2_all, beta3_all, beta4_all, beta5_all)
beta_train <- beta_all[train_index,]
beta_test <- beta_all[test_index,]
#########

######
# Prepare the containers for saving output
######

methods <- c("train_mean", "lm", "np", "btvcm", "vc_bart", "bart", "rf", "gbm")
lin_methods <- c("lm", "np", "btvcm", "vc_bart")

#lin_methods <- c("vc_bart", "np", "btvcm", "lm")

# For beta recovery
beta_mse_train <- matrix(nrow = length(lin_methods), ncol = p + 1, dimnames = list(lin_methods, c()))
beta_mse_test <- matrix(nrow = length(lin_methods), ncol = p + 1, dimnames = list(lin_methods, c()))
beta_coverage_train <- matrix(nrow = length(lin_methods), ncol = p + 1, dimnames = list(lin_methods, c()))
beta_coverage_test <- matrix(nrow = length(lin_methods), ncol = p + 1, dimnames = list(lin_methods, c()))
# ratio of pointwise interval lengths compared to VC-BART
beta_int_train <- matrix(nrow = length(lin_methods), ncol = p + 1, dimnames = list(lin_methods, c()))
beta_int_test <- matrix(nrow = length(lin_methods), ncol = p + 1, dimnames = list(lin_methods, c()))
# For predictions
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

beta_mse_train["vc_bart",] <- colMeans( (beta_train - vc_sum$train$beta[,"MEAN",])^2 , na.rm = TRUE)
beta_mse_test["vc_bart",] <- colMeans( (beta_test - vc_sum$test$beta[,"MEAN",])^2 , na.rm = TRUE)

beta_coverage_train["vc_bart",] <- apply( (beta_train >= vc_sum$train$beta[,"L95",]) & (beta_train <= vc_sum$train$beta[,"U95",]), FUN = mean, MARGIN = 2)
beta_coverage_test["vc_bart",] <- apply( (beta_test >= vc_sum$test$beta[,"L95",]) & (beta_test <= vc_sum$test$beta[,"U95",]), FUN = mean, MARGIN = 2)

beta_int_train["vc_bart",] <- 1
beta_int_test["vc_bart",] <- 1

ystar_mse_train["vc_bart"] <- mean( (Y_train - vc_sum$train$ystar[,"MEAN"])^2 )
ystar_mse_test["vc_bart"] <- mean( (Y_test - vc_sum$test$ystar[,"MEAN"])^2 )

ystar_coverage_train["vc_bart"] <- mean( (Y_train >= vc_sum$train$ystar[,"L95"]) & (Y_train <= vc_sum$train$ystar[,"U95"]))
ystar_coverage_test["vc_bart"] <- mean( (Y_test >= vc_sum$test$ystar[,"L95"]) & (Y_test <= vc_sum$test$ystar[,"U95"]))

ystar_int_train["vc_bart"] <- 1
ystar_int_test["vc_bart"] <- 1

######################
# Run least squares
######################

lm_time <- system.time(lm_fit <- lm(Y ~ X1 + X2 + X3  + X4 + X5, data = tmp_data_train))["elapsed"]

lm_pred_train <- predict(lm_fit, newdata = tmp_data_train, interval = "prediction")
lm_pred_test <- predict(lm_fit, newdata = tmp_data_test, interval = "prediction")

lm_coef <- as.data.frame(summary(lm_fit)$coef)
lm_coef[,"L95"] <- lm_coef[,"Estimate"] - qnorm(0.975) * lm_coef[,"Std. Error"]
lm_coef[,"U95"] <- lm_coef[,"Estimate"] + qnorm(0.975) * lm_coef[, "Std. Error"]
for(j in 1:(p+1)){
  beta_mse_train["lm",j] <- mean( (beta_train[,j] - lm_coef[j, "Estimate"])^2 )
  beta_mse_test["lm",j] <- mean( (beta_test[,j] - lm_coef[j, "Estimate"])^2 )
  
  beta_coverage_train["lm",j] <- mean( (beta_train[,j] >= lm_coef[j,"L95"]) & (beta_train[,j] <= lm_coef[j, "U95"]))
  beta_coverage_test["lm",j] <- mean( (beta_test[,j] >= lm_coef[j, "L95"]) & (beta_test[,j] <= lm_coef[j, "U95"]))
  
  beta_int_train["lm",j] <- mean( (lm_coef[j, "U95"] - lm_coef[j, "L95"])/(vc_sum$train$beta[,"U95",j] - vc_sum$train$beta[,"L95",j]) )
  beta_int_test["lm",j] <- mean( (lm_coef[j, "U95"] - lm_coef[j, "L95"])/(vc_sum$test$beta[,"U95",j] - vc_sum$test$beta[,"L95",j]) )
}

ystar_mse_train["lm"] <- mean( (Y_train - lm_pred_train[,"fit"])^2 )
ystar_mse_test["lm"] <- mean( (Y_test - lm_pred_test[,"fit"])^2 )

ystar_coverage_train["lm"] <- mean( (Y_train >= lm_pred_train[,"lwr"]) & (Y_train <= lm_pred_train[,"upr"]))
ystar_coverage_test["lm"] <- mean( (Y_test >= lm_pred_test[,"lwr"]) & (Y_test <= lm_pred_test[,"upr"]))

ystar_int_train["lm"] <- mean( (lm_pred_train[,"upr"] - lm_pred_train[,"lwr"])/(vc_sum$train$ystar[,"U95"] - vc_sum$train$ystar[,"L95"]) )
ystar_int_test["lm"] <- mean( (lm_pred_test[,"upr"] - lm_pred_test[,"lwr"])/(vc_sum$test$ystar[,"U95"] - vc_sum$test$ystar[,"L95"]) )

timing["lm"] <- lm_time

########################
# Run Kernel Smoothing
########################
np_fit <- kernel_smoothing(Y_train, cov_train, mod_train, cov_test, mod_test, B = 50)
if(!is.null(np_fit)){
  
  beta_mse_train["np",] <- colMeans( (beta_train - np_fit$train$beta[,"MEAN",])^2)
  beta_mse_test["np",] <- colMeans( (beta_test - np_fit$test$beta[,"MEAN",])^2)
  
  beta_coverage_train["np",] <- apply( (beta_train >= np_fit$train$beta[,"L95",]) & (beta_train <= np_fit$train$beta[,"U95",]), FUN = mean, MARGIN = 2)
  beta_coverage_test["np",] <- apply( (beta_test >= np_fit$test$beta[,"L95",]) & (beta_test <= np_fit$test$beta[,"U95",]), FUN = mean, MARGIN = 2)
  
  beta_int_train["np",] <- apply( (np_fit$train$beta[,"U95",] - np_fit$train$beta[,"L95",])/(vc_sum$train$beta[,"U95",] - vc_sum$train$beta[,"L95",]), MARGIN = 2, FUN = mean)
  beta_int_test["np",] <- apply( (np_fit$test$beta[,"U95",] - np_fit$test$beta[,"L95",])/(vc_sum$test$beta[,"U95",] - vc_sum$test$beta[,"L95",]), MARGIN = 2, FUN = mean)
  
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

btvcm_fit <- btvcm(Y_train, X_train, Z_train, X_test, Z_test, B = 50)
if(!is.null(btvcm_fit)){
  
  beta_mse_train["btvcm",] <- colMeans( (beta_train - btvcm_fit$train$beta[,"MEAN",])^2)
  beta_mse_test["btvcm",] <- colMeans( (beta_test - btvcm_fit$test$beta[,"MEAN",])^2)
  
  beta_coverage_train["btvcm",] <- apply( (beta_train >= btvcm_fit$train$beta[,"L95",]) & (beta_train <= btvcm_fit$train$beta[,"U95",]), FUN = mean, MARGIN = 2)
  beta_coverage_test["btvcm",] <- apply( (beta_test >= btvcm_fit$test$beta[,"L95",]) & (beta_test <= btvcm_fit$test$beta[,"U95",]), FUN = mean, MARGIN = 2)
  
  beta_int_train["btvcm",] <- apply( (btvcm_fit$train$beta[,"U95",] - btvcm_fit$train$beta[,"L95",])/(vc_sum$train$beta[,"U95",] - vc_sum$train$beta[,"L95",]), MARGIN = 2, FUN = mean)
  beta_int_test["btvcm",] <- apply( (btvcm_fit$test$beta[,"U95",] - btvcm_fit$test$beta[,"L95",])/(vc_sum$test$beta[,"U95",] - vc_sum$test$beta[,"L95",]), MARGIN = 2, FUN = mean)
  
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
# One last method: predicting out-of-sample with the mean
rmse_mean <- sqrt(mean((Y_train - mean(Y_train))^2))
ystar_mse_train["train_mean"] <- mean( (Y_train - mean(Y_train))^2 )
ystar_mse_test["train_mean"] <- mean( (Y_test - mean(Y_train))^2 )

ystar_coverage_train["train_mean"] <- mean( (Y_train >= mean(Y_train) - qnorm(0.975) * rmse_mean) & (Y_train <= mean(Y_train) + qnorm(0.975) * rmse_mean) )
ystar_coverage_test["train_mean"] <- mean( (Y_test >= mean(Y_train) - qnorm(0.975) * rmse_mean) & (Y_test <= mean(Y_train) + qnorm(0.975) * rmse_mean) )

ystar_int_train["train_mean"] <- mean( 2 * qnorm(0.975) * rmse_mean / (vc_sum$train$ystar[,"U95"] - vc_sum$train$ystar[,"L95"]) )
ystar_int_test["train_mean"] <- mean( 2 * qnorm(0.975) * rmse_mean / (vc_sum$test$ystar[,"U95"] - vc_sum$test$ystar[,"L95"]) )

########

assign(paste0("beta_mse_train_", rep), beta_mse_train)
assign(paste0("beta_mse_test_", rep), beta_mse_test)
assign(paste0("beta_coverage_train_", rep), beta_coverage_train)
assign(paste0("beta_coverage_test_", rep), beta_coverage_test)
assign(paste0("beta_int_train_", rep), beta_int_train)
assign(paste0("beta_int_test_", rep), beta_int_test)

assign(paste0("ystar_mse_train_", rep), ystar_mse_train)
assign(paste0("ystar_mse_test_", rep), ystar_mse_test)
assign(paste0("ystar_coverage_train_", rep), ystar_coverage_train)
assign(paste0("ystar_coverage_test_", rep), ystar_coverage_test)
assign(paste0("ystar_int_train_", rep), ystar_int_train)
assign(paste0("ystar_int_test_", rep), ystar_int_test)
assign(paste0("timing_", rep), timing)

save_list <- paste0(c("beta_mse_train_", "beta_coverage_train_", "beta_int_train_", 
                      "beta_mse_test_", "beta_coverage_test_", "beta_int_test_",
                      "ystar_mse_train_", "ystar_coverage_train_", "ystar_int_train_", 
                      "ystar_mse_test_", "ystar_coverage_test_", "ystar_int_test_",
                      "timing_"), rep)
save(list = save_list, file = paste0("results/results_p5R5_", rep, ".RData"))






