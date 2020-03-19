# Assess the out-of-sample performance of all methods on the Philadelphia crime data

methods <- c(paste0("cs_rho", 0:9), "bart", "gbm", "rf", "lm", "tvcm")
load("data/philly_monthly_cc.RData")
n_train <- n_vec_train[1]
n_test <- n_vec_test[1]
n <- length(n_vec_train) # number of census tracts

for(t in 1:n_test){
  assign(paste0("test_index",t), (0:254)*n_test + t) # indices of test set observations
}

ystar_mse_train <- rep(NA, times = length(methods))
ystar_coverage_train <- rep(NA, times = length(methods))
ystar_int_train <- rep(NA, times = length(methods))

names(ystar_mse_train) <- methods
names(ystar_coverage_train) <- methods
names(ystar_int_train) <- methods

ystar_mse_test <- matrix(nrow = length(methods), ncol = n_test, dimnames = list(methods, c()))
ystar_coverage_test <- matrix(nrow = length(methods), ncol = n_test, dimnames = list(methods, c()))
ystar_int_test <- matrix(nrow = length(methods), ncol = n_test, dimnames = list(methods, c()))

load("results/philly_cs_0.RData") # load independent VC-BART

for(rho in 0:9){
  if(file.exists(paste0("results/philly_cs_", rho, ".RData"))){
    load(paste0("results/philly_cs_", rho, ".RData"))
    fit_sum <- get(paste0("cs_rho", rho, "_sum"))
    
    ystar_mse_train[paste0("cs_rho",rho)] <- mean( (Y_train - fit_sum$train$ystar[,"MEAN"])^2)
    ystar_coverage_train[paste0("cs_rho",rho)] <- mean( (Y_train >= fit_sum$train$ystar[,"L95"]) & (Y_train <= fit_sum$train$ystar[,"U95"]))
    ystar_int_train[paste0("cs_rho",rho)] <- mean( (fit_sum$train$ystar[,"U95"] - fit_sum$train$ystar[,"L95"])/ (cs_rho0_sum$train$ystar[,"U95"] - cs_rho0_sum$train$ystar[,"L95"]))
    
    for(t in 1:n_test){
      index <- get(paste0("test_index",t))
      ystar_mse_test[paste0("cs_rho",rho), t] <- mean( (Y_test[index] - fit_sum$test$ystar[index,"MEAN"])^2 )
      ystar_coverage_test[paste0("cs_rho",rho),t] <- mean( (Y_test[index] >= fit_sum$test$ystar[index,"L95"]) & (Y_test[index] <= fit_sum$test$ystar[index,"U95"]))
      ystar_int_test[paste0("cs_rho", rho),t] <- mean( (fit_sum$test$ystar[index,"U95"] - fit_sum$test$ystar[index,"L95"])/(cs_rho0_sum$test$ystar[index,"U95"] - cs_rho0_sum$test$ystar[index,"L95"]))
    }
    if(rho != 0) rm(list = paste0("cs_rho",rho, "_sum"))
  }
}
rm(fit_sum)
# Get the BART fits
if(file.exists(paste0("results/philly_bart.RData"))){
  load("results/philly_bart.RData")
  
  ystar_mse_train["bart"] <- mean( (Y_train - bart_sum$train$ystar[,"MEAN"])^2)
  ystar_coverage_train["bart"] <- mean( (Y_train >= bart_sum$train$ystar[,"L95"]) & (Y_train <= bart_sum$train$ystar[,"U95"]))
  ystar_int_train["bart"] <- mean( (bart_sum$train$ystar[,"U95"] - bart_sum$train$ystar[,"L95"])/ (cs_rho0_sum$train$ystar[,"U95"] - cs_rho0_sum$train$ystar[,"L95"]))
  
  for(t in 1:n_test){
    index <- get(paste0("test_index",t))
    ystar_mse_test["bart", t] <- mean( (Y_test[index] - bart_sum$test$ystar[index,"MEAN"])^2 )
    ystar_coverage_test["bart",t] <- mean( (Y_test[index] >= bart_sum$test$ystar[index,"L95"]) & (Y_test[index] <= bart_sum$test$ystar[index,"U95"]))
    ystar_int_test["bart",t] <- mean( (bart_sum$test$ystar[index,"U95"] - bart_sum$test$ystar[index,"L95"])/(cs_rho0_sum$test$ystar[index,"U95"] - cs_rho0_sum$test$ystar[index,"L95"]))
  }
  rm(bart_sum)
}

# Get RF
if(file.exists(paste0("results/philly_rf.RData"))){
  load("results/philly_rf.RData")
  ystar_mse_train["rf"] <- mean( (Y_train - rf_pred_train)^2)
  ystar_coverage_train["rf"] <- mean( (Y_train >= rf_pred_train - qnorm(0.975) * rf_rmse) & (Y_train <= rf_pred_train + qnorm(0.975) *rf_rmse))
  ystar_int_train["rf"] <- mean( 2* qnorm(0.975) * rf_rmse/(cs_rho0_sum$train$ystar[,"U95"] - cs_rho0_sum$train$ystar[,"L95"]))
  for(t in 1:n_test){
    index <- get(paste0("test_index",t))
    ystar_mse_test["rf",t] <- mean( (Y_test[index] - rf_pred_test[index])^2 )
    ystar_coverage_test["rf",t] <- mean( (Y_test[index] >= rf_pred_test - qnorm(0.975) * rf_rmse) & (Y_test[index] <= rf_pred_test + qnorm(0.975) *rf_rmse))
    ystar_int_test["rf",t] <- mean( 2 * qnorm(0.975) * rf_rmse/(cs_rho0_sum$test$ystar[index,"U95"] - cs_rho0_sum$test$ystar[index,"L95"]))
  }
}


if(file.exists(paste0("results/philly_gbm.RData"))){
  load("results/philly_gbm.RData")
  
  ystar_mse_train["gbm"] <- mean( (Y_train - gbm_pred_train)^2)
  ystar_coverage_train["gbm"] <- mean( (Y_train >= gbm_pred_train - qnorm(0.975) * gbm_rmse) & (Y_train <= gbm_pred_train + qnorm(0.975) *gbm_rmse))
  ystar_int_train["gbm"] <- mean( 2* qnorm(0.975) * gbm_rmse/(cs_rho0_sum$train$ystar[,"U95"] - cs_rho0_sum$train$ystar[,"L95"]))
  for(t in 1:n_test){
    index <- get(paste0("test_index",t))
    ystar_mse_test["gbm",t] <- mean( (Y_test[index] - gbm_pred_test[index])^2 )
    ystar_coverage_test["gbm",t] <- mean( (Y_test[index] >= gbm_pred_test - qnorm(0.975) * gbm_rmse) & (Y_test[index] <= gbm_pred_test + qnorm(0.975) *gbm_rmse))
    ystar_int_test["gbm",t] <- mean( 2 * qnorm(0.975) * gbm_rmse/(cs_rho0_sum$test$ystar[index,"U95"] - cs_rho0_sum$test$ystar[index,"L95"]))
  }
}


if(file.exists(paste0("results/philly_tvcm.RData"))){
  load("results/philly_tvcm.RData")
  ystar_mse_train["tvcm"] <- mean( (Y_train - tvcm_fit$train$ystar[,"MEAN"])^2)
  ystar_coverage_train["tvcm"] <- mean( (Y_train >= tvcm_fit$train$ystar[,"L95"]) & (Y_train <= tvcm_fit$train$ystar[,"U95"]))
  ystar_int_train["tvcm"] <- mean((tvcm_fit$train$ystar[,"U95"] - tvcm_fit$train$ystar[,"L95"])/(cs_rho0_sum$train$ystar[,"U95"] - cs_rho0_sum$train$ystar[,"L95"]))
  for(t in 1:n_test){
    index <- get(paste0("test_index",t))
    ystar_mse_test["tvcm",t] <- mean( (Y_test[index] - tvcm_fit$test$ystar[index, "MEAN"])^2 )
    ystar_coverage_test["tvcm",t] <- mean( (Y_test[index] >= tvcm_fit$test$ystar[index,"L95"]) & (Y_test[index] <= tvcm_fit$test$ystar[index, "U95"]))
    ystar_int_test["tvcm",t] <- mean( (tvcm_fit$test$ystar[index, "U95"] - tvcm_fit$test$ystar[index,"L95"])/(cs_rho0_sum$test$ystar[index,"U95"] - cs_rho0_sum$test$ystar[index,"L95"]))
    
  }
}

if(file.exists("results/philly_lm.RData")){
  load("results/philly_lm.RData")
  
  ystar_mse_train["lm"] <- mean( (Y_train - lm_pred_train[,"fit"])^2 )
  ystar_coverage_train["lm"] <- mean( (Y_train >= lm_pred_train[,"lwr"]) & (Y_train <= lm_pred_train[,"upr"]))
  ystar_int_train["lm"] <- mean( (lm_pred_train[,"upr"] - lm_pred_train[,"lwr"])/(cs_rho0_sum$train$ystar[,"U95"] - cs_rho0_sum$train$ystar[,"L95"]))
  for(t in 1:n_test){
    index <- get(paste0("test_index",t))
    ystar_mse_test["lm",t] <- mean( (Y_test[index] - lm_pred_test[index,"fit"])^2 )
    ystar_coverage_test["lm",t] <- mean( (Y_test[index] >= lm_pred_test[index,"lwr"]) & (Y_test[index] <= lm_pred_test[index, "upr"]))
    ystar_int_test["lm",t] <- mean( (lm_pred_test[,"upr"] - lm_pred_test[,"lwr"])/(cs_rho0_sum$test$ystar[index,"U95"] - cs_rho0_sum$test$ystar[index,"L95"]))
    
  }
  
}

new_order <- c("lm", "tvcm", "bart", "rf", "gbm", "cs_rho0", "cs_rho1", "cs_rho8")

round(ystar_mse_test[new_order,], digits = 2)
round(rowMeans(ystar_mse_test[new_order,]), digits = 2)