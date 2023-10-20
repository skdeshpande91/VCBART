# wrapper for gbm
library(gbm)

gbm_wrapper <- function(Y_train, X_train, Z_train, X_test, Z_test){
  
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
  
  
  train_time <- system.time(
    gbm_fit <- tryCatch(
      gbm(Y ~ ., data = tmp_data_train, distribution = "gaussian",
          interaction.depth = 6, cv.folds = 5, n.trees = 1000)
    )
  )["elapsed"]
  

  if(!is.null(gbm_fit)){
    n_tree_opt <- gbm.perf(gbm_fit, method = "cv", plot.it = FALSE)
    ystar_train <- predict(gbm_fit, newdata = tmp_data_train, n.trees = n_tree_opt)
    ystar_test <- predict(gbm_fit, newdata = tmp_data_test, n.trees = n_tree_opt)
    
    rmse_train <- sqrt(mean( (Y_train - ystar_train)^2 ))
    
    ystar_summary_train <- matrix(nrow = n_obs_train, ncol = 3, dimnames = list(c(), c("MEAN", "L95", "U95")))
    ystar_summary_test <- matrix(nrow = n_obs_test, ncol = 3, dimnames = list(c(), c("MEAN", "L95", "U95")))
    ystar_summary_train[,"MEAN"] <- ystar_train
    ystar_summary_test[,"MEAN"] <- ystar_test
    
    ystar_summary_train[,"L95"] <- ystar_summary_train[,"MEAN"] - qnorm(0.975) * rmse_train
    ystar_summary_train[,"U95"] <- ystar_summary_train[,"MEAN"] + qnorm(0.97) * rmse_train
    ystar_summary_test[,"L95"] <- ystar_summary_test[,"MEAN"] - qnorm(0.975) * rmse_train
    ystar_summary_test[,"U95"] <- ystar_summary_test[,"MEAN"] + qnorm(0.975 ) * rmse_train
  
    results <- list(train = list(ystar = ystar_summary_train),
                    test = list(ystar = ystar_summary_test),
                    time = train_time, train_time = train_time) 
  } else{
    results <- NULL
  }

  return(results)
}
