# Wrapper function for the tvcglm 
library(vcrpart)
tvc_wrapper <- function(Y_train, X_train, Z_train, X_test, Z_test,B = 50){
  
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
  
  # Define the formula
  form_list <- rep(NA, times = p+1)
  
  form_list[1] <- paste0("vc(", paste0(colnames(Z_train), collapse = ", "), ")")
  for(j in 1:p){
    form_list[j+1] <- paste0("vc(", paste0(colnames(Z_train), collapse = ", "), ", by = X", j, ")")
  }
  
  tvc_form <- formula(paste0("Y ~ -1 + ", paste0(form_list, collapse = "+")))
  
  # Get the initial fit
  train_time <- system.time(
    tvc_fit <- tryCatch(tvcglm(formula = tvc_form,
                               family = gaussian(),
                               data = tmp_data_train,
                               control = tvcglm_control(cv = TRUE)),
                        error = function(e){return(NULL)}))["elapsed"]
  
  if(!is.null(tvc_fit)){
    # call cvloss in order to get the cross-validated cp parameters
    tvc_cp <- cvloss(tvc_fit)$cp.hat
    
    # Get the initial estimates of beta
    beta_hat_train <- predict(tvc_fit, newdata = tmp_data_train, type = "coef")
    beta_hat_test <- predict(tvc_fit, newdata = tmp_data_test, type = "coef")
    
    # Get the initial point estimates of ystar
    ystar_train <- predict(tvc_fit, newdata = tmp_data_train, type = "response")
    ystar_test <- predict(tvc_fit, newdata = tmp_data_test, type = "response")
    
    if(B == 0){
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
                      time = train_time,
                      train_time = train_time)     
    } else{
      # Bootstrap time
      boot_beta <- array(dim = c(n_obs_train + n_obs_test, p + 1, B))
      
      boot_time <- system.time({
        for(b in 1:B){
          if(b %% 5 == 0) print(paste("tvc-boot: b = ", b, "at", Sys.time()))
          
          boot_index <- sample(1:n_obs_train, size = n_obs_train,replace = TRUE)
          tryCatch(
            {
              boot_fit <- tvcglm(formula = tvc_form,
                                 family = gaussian(),
                                 data = tmp_data_train[boot_index,],
                                 control = tvcglm_control(cv = FALSE))
              boot_pruned_fit <- prune(boot_fit, cp = tvc_cp)
              tmp_beta <- predict(boot_pruned_fit, newdata = rbind(tmp_data_train, tmp_data_test), type = "coef")
            }, 
            error = function(e){return(NULL)})
          if(!is.null(tmp_beta)) boot_beta[,,b] <- tmp_beta
        }
      })["elapsed"]
      
        boot_beta_se <- apply(boot_beta, FUN = sd, MAR = c(1,2), na.rm = TRUE)
      
        beta_summary_train <- array(dim = c(n_obs_train, 4, 1 + p), dimnames = list(c(), c("MEAN", "SD", "L95", "U95"), c()))
        beta_summary_test <- array(dim = c(n_obs_test, 4, 1 + p), dimnames = list(c(), c("MEAN", "SD", "L95", "U95"), c()))
        
        ystar_summary_train <- matrix(nrow = n_obs_train, ncol = 3, dimnames = list(c(), c("MEAN", "L95", "U95")))
        ystar_summary_test <- matrix(nrow = n_obs_test, ncol = 3, dimnames = list(c(), c("MEAN", "L95", "U95")))
        
        beta_summary_train[,"MEAN",] <- beta_hat_train
        beta_summary_test[,"MEAN",] <- beta_hat_test
        
        beta_summary_train[,"SD",] <- boot_beta_se[1:n_obs_train,]
        beta_summary_test[,"SD",] <- boot_beta_se[(1 + n_obs_train):(n_obs_train + n_obs_test),]
        
        beta_summary_train[,"L95",] <- beta_summary_train[,"MEAN",] - qnorm(0.975) * beta_summary_train[,"SD",]
        beta_summary_train[,"U95",] <- beta_summary_train[,"MEAN",] + qnorm(0.975) * beta_summary_train[,"SD",]
        
        beta_summary_test[,"L95",] <- beta_summary_test[,"MEAN",] - qnorm(0.975) * beta_summary_test[,"SD",]
        beta_summary_test[,"U95",] <- beta_summary_test[,"MEAN",] + qnorm(0.975) * beta_summary_test[,"SD",]
        
        rmse_train <- sqrt(mean( (Y_train - ystar_train)^2))
        ystar_summary_train[,"MEAN"] <- ystar_train
        ystar_summary_test[,"MEAN"] <- ystar_test
        
        ystar_summary_train[,"L95"] <- ystar_summary_train[,"MEAN"] - qnorm(0.975) * rmse_train
        ystar_summary_train[,"U95"] <- ystar_summary_train[,"MEAN"] + qnorm(0.97) * rmse_train
      
        ystar_summary_test[,"L95"] <- ystar_summary_test[,"MEAN"] - qnorm(0.975) * rmse_train
        ystar_summary_test[,"U95"] <- ystar_summary_test[,"MEAN"] + qnorm(0.97) * rmse_train
        
        results <- list(train = list(beta = beta_summary_train, ystar = ystar_summary_train),
                        test = list(beta = beta_summary_test, ystar = ystar_summary_test),
                        time = train_time + boot_time,
                        train_time = train_time,
                        boot_time = boot_time)
      
    } # close else checking that we need to do bootstrap
    
  }else{
    # if the initial fit failed, return null
    results <- NULL
  }
  return(results)
 
}