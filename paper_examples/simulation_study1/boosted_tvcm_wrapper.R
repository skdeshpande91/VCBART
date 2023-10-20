# Wrapper function for bTVCM
source("honestRpart.R") # this is to compare with Zhou & Hooker (2019)
source("lgd.R") # this is to compare with Zhou & Hooker (2019)

boosted_tvcm_wrapper <- function(Y_train, X_train, Z_train, X_test, Z_test,
                                 intercept = TRUE, lambda = 0.05, ntree = 200, B = 50)
{
  n_obs_train <- nrow(X_train)
  n_obs_test <- nrow(X_test)
  p <- ncol(X_train)
  
  model <- list()
  model$dummy <- intercept
  model$xscale <- TRUE
  model$n <- n_obs_train
  model$p <- p + model$dummy
  model$q <- ncol(Z_train)
  model$diff <- ols.diff
  model$diff <- ols.diff
  model$init <- ols.init.beta
  model$pred <- ols.pred
  model$loss <- ols.loss
  model$lambda <- lambda
  model$subsample <- 0.5
  model$ntree <- ntree
  model$control <- rpart.control(maxdepth = 6, cp=0.0001)
  model$woods <- list()
  model$savetree <- TRUE
  
  time_train <- system.time(  
    tvcm_fit <- tryCatch(lgd(y=Y_train, x=X_train, z=Z_train, model=model),
                         error = function(e){return(NULL)}))["elapsed"]
  
  if(!is.null(tvcm_fit)){
    yhat_train <- tvcm_fit$yhat
    time_test <- system.time(yhat_test <- lgd.predict(X_test, Z_test, tvcm_fit))["elapsed"]
    beta_train <- tvcm_fit$beta
    
    beta_test <- matrix(nrow = nrow(X_test), ncol = 1 + p)
    beta_time <- 0
    tmp_X <- matrix(0, nrow = nrow(X_test), ncol = p)
    beta_time <- system.time(beta_test[,1] <- lgd.predict(tmp_X, Z_test, tvcm_fit))["elapsed"]
    for(j in 1:ncol(X_test)){
      tmp_X <- matrix(0, nrow = nrow(X_test), ncol = p)
      tmp_X[,j] <- 1
      tmp_time <- system.time(beta_test[,j+1] <- lgd.predict(tmp_X, Z_test, tvcm_fit) - beta_test[,1])["elapsed"] # remember to remove the intercept!
      beta_time <- beta_time + tmp_time
    }
    if(B == 0){
      ystar_summary_train <- matrix(nrow = n_obs_train, ncol = 3, dimnames = list(c(), c("MEAN", "L95", "U95")))
      ystar_summary_test <- matrix(nrow = n_obs_test, ncol = 3, dimnames = list(c(), c("MEAN", "L95", "U95")))
      rmse_train <- sqrt(mean( (Y_train - yhat_train)^2))
      ystar_summary_train[,"MEAN"] <- yhat_train
      ystar_summary_test[,"MEAN"] <- yhat_test
      
      ystar_summary_train[,"L95"] <- ystar_summary_train[,"MEAN"] - qnorm(0.975) * rmse_train
      ystar_summary_train[,"U95"] <- ystar_summary_train[,"MEAN"] + qnorm(0.97) * rmse_train
      
      ystar_summary_test[,"L95"] <- ystar_summary_test[,"MEAN"] - qnorm(0.975) * rmse_train
      ystar_summary_test[,"U95"] <- ystar_summary_test[,"MEAN"] + qnorm(0.975 ) * rmse_train
      results <- list(train = list(ystar = ystar_summary_train),
                      test = list(ystar = ystar_summary_test),
                      time = time_train + beta_time,
                      train_time = time_train)      
    } else{
      ######
      # Do the bootstrap now
      ######
      boot_beta <- array(dim = c(n_obs_train + n_obs_test, p + 1, B))
      
      boot_time <- tryCatch(
        {
          system.time(
            {
              for(b in 1:B){
                if(b %% 5 == 0) print(paste("tvcm-boot: b = ", b, Sys.time()))
                boot_train <- sample(1:n_obs_train, size = n_obs_train, replace = TRUE)
                boot_fit <- lgd(y=Y_train[boot_train], x=X_train[boot_train,], z=Z_train[boot_train,], model=model)
                
                boot_beta[1:n_obs_train,,b] <- boot_fit$beta
                
                tmp_X <- matrix(0, nrow = n_obs_train + n_obs_test, ncol = ncol(X_train))
                tmp_int <- lgd.predict(tmp_X, rbind(Z_train, Z_test), boot_fit)
                
                boot_beta[,1,b] <- tmp_int
                for(j in 1:p){
                  tmp_X <- matrix(0, nrow = n_obs_train + n_obs_test, ncol = p)
                  tmp_X[,j] <- 1
                  boot_beta[,j+1,b] <- lgd.predict(tmp_X, rbind(Z_train, Z_test), boot_fit) - tmp_int
                }
              }
            }
          )["elapsed"] 
        }, error = function(e){return(NULL)})
      if(!is.null(boot_time)){
        boot_beta_se <- apply(boot_beta, FUN = sd, MAR = c(1,2))
        
        beta_summary_train <- array(dim = c(n_obs_train, 4, 1 + p), dimnames = list(c(), c("MEAN", "SD", "L95", "U95"), c()))
        beta_summary_test <- array(dim = c(n_obs_test, 4, 1 + p), dimnames = list(c(), c("MEAN", "SD", "L95", "U95"), c()))
        
        ystar_summary_train <- matrix(nrow = n_obs_train, ncol = 3, dimnames = list(c(), c("MEAN", "L95", "U95")))
        ystar_summary_test <- matrix(nrow = n_obs_test, ncol = 3, dimnames = list(c(), c("MEAN", "L95", "U95")))
        
        beta_summary_train[,"MEAN",] <- beta_train
        beta_summary_test[,"MEAN",] <- beta_test
        
        beta_summary_train[,"SD",] <- boot_beta_se[1:n_obs_train,]
        beta_summary_test[,"SD",] <- boot_beta_se[(1 + n_obs_train):(n_obs_train + n_obs_test),]
        
        beta_summary_train[,"L95",] <- beta_summary_train[,"MEAN",] - qnorm(0.975) * beta_summary_train[,"SD",]
        beta_summary_train[,"U95",] <- beta_summary_train[,"MEAN",] + qnorm(0.975) * beta_summary_train[,"SD",]
        
        beta_summary_test[,"L95",] <- beta_summary_test[,"MEAN",] - qnorm(0.975) * beta_summary_test[,"SD",]
        beta_summary_test[,"U95",] <- beta_summary_test[,"MEAN",] + qnorm(0.975) * beta_summary_test[,"SD",]
        
        rmse_train <- sqrt(mean( (Y_train - yhat_train)^2))
        ystar_summary_train[,"MEAN"] <- yhat_train
        ystar_summary_test[,"MEAN"] <- yhat_test
        
        ystar_summary_train[,"L95"] <- ystar_summary_train[,"MEAN"] - qnorm(0.975) * rmse_train
        ystar_summary_train[,"U95"] <- ystar_summary_train[,"MEAN"] + qnorm(0.975) * rmse_train
        
        
        ystar_summary_test[,"L95"] <- ystar_summary_test[,"MEAN"] - qnorm(0.975) * rmse_train
        ystar_summary_test[,"U95"] <- ystar_summary_test[,"MEAN"] + qnorm(0.975 ) * rmse_train
        
        
        results <- list(train = list(beta = beta_summary_train, ystar = ystar_summary_train),
                        test = list(beta = beta_summary_test, ystar = ystar_summary_test),
                        time = time_train + beta_time + boot_time,
                        train_time = time_train,
                        boot_time = boot_time)
      } else results <- NULL  
    }
  } else results <- NULL
  return(results)

}