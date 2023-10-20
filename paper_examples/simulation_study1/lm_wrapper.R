# quick wrapper for lm
lm_wrapper <- function(Y_train, X_train, Z_train, X_test, Z_test){
  
  n_obs_train <- nrow(X_train)
  n_obs_test <- nrow(X_test)
  
  p <- ncol(X_train)
  R <- ncol(Z_train)
  
  # Define necessary data frames
  colnames(X_train) <- paste0("X", 1:p)
  colnames(X_test) <- paste0("X",1:p)
  
  colnames(Z_train) <- paste0("Z", 1:R)
  colnames(Z_test) <- paste0("Z", 1:R)
  
  tmp_data_train <- data.frame(Y_train,X_train, Z_train)
  tmp_data_test <- data.frame(Y_test,X_test, Z_test)
  colnames(tmp_data_train)[1] <- "Y"
  colnames(tmp_data_test)[1] <- "Y"
  
  
  form <- formula(paste0("Y ~", paste0(colnames(X_train), collapse = " + ")))
  train_time <- system.time(
    lm_fit <- lm(form, data = tmp_data_train)
  )["elapsed"]
  
  lm_coef <- as.data.frame(summary(lm_fit)$coef)
  
  beta_summary_train <- array(dim = c(n_obs_train, 4, 1 + p), dimnames = list(c(), c("MEAN", "SD", "L95", "U95"), c()))
  beta_summary_test <- array(dim = c(n_obs_test, 4, 1 + p), dimnames = list(c(), c("MEAN", "SD", "L95", "U95"), c()))
  
  for(k in 1:(p+1)){
    beta_summary_train[,"MEAN",k] <- lm_coef[k, "Estimate"]
    beta_summary_train[,"SD",k] <- lm_coef[k,"Std. Error"]
    beta_summary_test[,"MEAN",k] <- lm_coef[k, "Estimate"]
    beta_summary_test[,"SD",k] <- lm_coef[k,"Std. Error"]
  }
  
  beta_summary_train[,"L95",] <- beta_summary_train[,"MEAN",] - qnorm(0.975) * beta_summary_train[,"SD",]
  beta_summary_train[,"U95",] <- beta_summary_train[,"MEAN",] + qnorm(0.975) * beta_summary_train[,"SD",]
  
  beta_summary_test[,"L95",] <- beta_summary_test[,"MEAN",] - qnorm(0.975) * beta_summary_test[,"SD",]
  beta_summary_test[,"U95",] <- beta_summary_test[,"MEAN",] + qnorm(0.975) * beta_summary_test[,"SD",]
  
  ystar_summary_train <- predict(lm_fit, newdata = tmp_data_train, interval = "prediction")
  ystar_summary_test <- predict(lm_fit, newdata = tmp_data_test, interval = "prediction")
  
  colnames(ystar_summary_train) <- c("MEAN", "L95", "U95")
  colnames(ystar_summary_test) <- c("MEAN", "L95", "U95")
  
  results <- list(train = list(beta = beta_summary_train, ystar = ystar_summary_train),
                  test = list(beta = beta_summary_test, ystar = ystar_summary_test),
                  time = train_time, train_time = train_time)
  return(results)
}
                  