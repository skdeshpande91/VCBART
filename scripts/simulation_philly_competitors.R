# Run competing methods (i.e. lm, BART, GBM, ERT, bTVCM) on the Philadelphia data 

library(MASS)
library(gbm)
library(ranger)
library(BART)
source("scripts/btvcm_wrapper.R")
source("scripts/get_summary.R")
load("data/philly_monthly_cc.RData")

p <- ncol(X_train)
R <- ncol(Z_train)

###############
# Set up data containers
###############

tmp_data_all <- data.frame("Y" = c(Y_train, Y_test), rbind(X_train, X_test), rbind(Z_train, Z_test))
cov_names <- colnames(tmp_data_all)[2:(p+1)]
mod_names <- colnames(tmp_data_all)[(p+2):(p+4)]

cov_all <- tmp_data_all[, cov_names]
mod_all <- tmp_data_all[, mod_names]

tmp_data_train <- data.frame("Y" = Y_train, X_train, Z_train)
tmp_data_test <- data.frame("Y" = Y_test, X_test, Z_test)

cov_train <- tmp_data_train[, cov_names]
cov_test <- tmp_data_test[, cov_names]

mod_train <- tmp_data_train[, mod_names]
mod_test <- tmp_data_test[, mod_names]

#####################
# Output containers
#####################
save_list <- c()

############
# Run BART
############

if(!file.exists("results/philly_bart.RData")){
  bart_chain1 <- wbart(cbind(X_train, Z_train), Y_train, cbind(X_test, Z_test), printevery = 5000)
  bart_chain2 <- wbart(cbind(X_train, Z_train), Y_train, cbind(X_test, Z_test), printevery = 5000)
  
  bart_sum <- bart_summary(bart_chain1, bart_chain2)
  bart_time <-  0.5 * (bart_chain1$proc.time["elapsed"] + bart_chain2$proc.time["elapsed"])
  save(bart_sum, bart_time, file = "results/philly_bart.RData")
}



###################################
# Run extremely randomized trees
###################################

if(!file.exists("results/philly_rf.RData")){
  rf_time <- system.time(rf_fit <- 
                           tryCatch(ranger(Y~ ., data = tmp_data_train, write.forest = TRUE, splitrule = "extratrees"),
                                    error = function(e){return(NULL)}))["elapsed"]
  if(!is.null(rf_fit)){
    rf_pred_train <- predict(rf_fit, data = tmp_data_train[,-1])$predictions
    rf_pred_test <- predict(rf_fit, data = tmp_data_test[,-1])$predictions
    rf_rmse <- sqrt( mean( (Y_train - rf_pred_train)^2 ))
    save(rf_pred_train, rf_pred_test, rf_rmse, rf_time, file = "results/philly_rf.RData")
  } 
}



#########
# Run GBM
#########
if(!file.exists("results/philly_gbm.RData")){
  gbm_time <- system.time(
    gbm_fit <- tryCatch(gbm(Y ~ ., data = tmp_data_train, distribution = "gaussian", interaction.depth = 6, cv.folds = 5, n.trees = 1000),
                        error = function(e){return(NULL)}))["elapsed"]
  if(!is.null(gbm_fit)){
    n_tree_opt <- gbm.perf(gbm_fit, method = "cv", plot.it = FALSE)
    gbm_pred_train <- predict(gbm_fit, newdata = tmp_data_train, n.trees = n_tree_opt)
    gbm_pred_test <- predict(gbm_fit, newdata = tmp_data_test, n.trees = n_tree_opt)
    
    gbm_rmse <- sqrt( mean( (Y_train - gbm_pred_train)^2 ))
    save(gbm_pred_train, gbm_pred_test, gbm_rmse, gbm_time, file = "results/philly_gbm.RData")
  }
}



#########
# Run least squares
##########
if(!file.exists("results/philly_lm.RData")){
  lm_time <- system.time(lm_fit <- lm(Y ~ log.income13 + sqrt.poverty13 + sqrt.vacantprop + 
                                        sqrt.comresprop + pop_total + segregation, data = tmp_data_train))["elapsed"]
  lm_pred_train <- predict(lm_fit, newdata = tmp_data_train, interval = "prediction")
  lm_pred_test <- predict(lm_fit, newdata = tmp_data_test, interval = "prediction")
  
  save(lm_pred_train, lm_pred_test, lm_time, file = "results/philly_lm.RData")
  
}


############
# Run bTVCM
if(!file.exists("results/philly_btvcm.RData")){
  btvcm_fit <- btvcm(Y_train, X_train, Z_train, X_test, Z_test, B = 0)
  save(tvcm_fit, file = "results/philly_btvcm.RData")
}


# Run np
if(!file.exists("results/philly_np.RData")){
  np_fit <- kernel_smoothing(Y_train, cov_train, mod_train, cov_test, mod_test, B = 0)
  save(np_fit, file = "results/philly_np.RData")
}
