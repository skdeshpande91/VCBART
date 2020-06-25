# wrapper for kernel smoothing

# Prepare the outcomes
# beta_summary_ gives the point estimate, bootstrap standard error, and bootstrap CI for each beta
# fit_summary_ gives the point estimate, SE returned by npscoef, and approximate CI for E[y | x,z]
# ypred_summary_ gives the point estimate and upper/lower 95 prediction interval. Here we just point estimate +/- 2 * RMSE.

library(np)
kernel_smoothing_wrapper <- function(Y_train, cov_train, mod_train, cov_test, mod_test, B = 50)
{
  n_train <- nrow(cov_train)
  n_test <- nrow(cov_test)
  
  bw_time <- system.time(
    bw <- tryCatch(npscoefbw(xdat = cov_train, ydat = Y_train, zdat = mod_train),
                   error = function(e){return(NULL)}))["elapsed"]
  if(!is.null(bw)){
    ks_fit <- npscoef(bw, betas = TRUE, residuals = TRUE, errors = TRUE,
                      txdat = cov_train, tydat = Y_train, tzdat = mod_train,
                      exdat = rbind(cov_train, cov_test), eydat = c(Y_train, Y_test), ezdat = rbind(mod_train, mod_test))
    yhat_train <- ks_fit$mean[1:length(Y_train)]
    yhat_test <- ks_fit$mean[(1 + length(Y_train)):(length(Y_test) + length(Y_train))]
    
    if(B == 0){
      
      ystar_summary_train <- matrix(nrow = n_train, ncol = 3, dimnames = list(c(), c("MEAN", "L95", "U95")))
      ystar_summary_test <- matrix(nrow = n_test, ncol = 3, dimnames = list(c(), c("MEAN", "L95", "U95")))
      rmse_train <- sqrt(mean( (Y_train - yhat_train)^2))
      ystar_summary_train[,"MEAN"] <- yhat_train
      ystar_summary_test[,"MEAN"] <- yhat_test
      
      ystar_summary_train[,"L95"] <- ystar_summary_train[,"MEAN"] - qnorm(0.975) * rmse_train
      ystar_summary_train[,"U95"] <- ystar_summary_train[,"MEAN"] + qnorm(0.97) * rmse_train
      
      ystar_summary_test[,"L95"] <- ystar_summary_test[,"MEAN"] - qnorm(0.975) * rmse_train
      ystar_summary_test[,"U95"] <- ystar_summary_test[,"MEAN"] + qnorm(0.975 ) * rmse_train
      results <- list(train = list(ystar = ystar_summary_train),
                      test = list(ystar = ystar_summary_test),
                      time = bw_time)      
    } else{
      # Let's try a bootstrap
      boot_beta <- array(dim = c(n_train + n_test, ncol(cov_train) + 1, B))
      boot_time <- tryCatch(
        {
          system.time(
            {
              for(b in 1:B){
                if(b %% 5 == 0) print(paste("np-boot: b = ", b, Sys.time()))
                boot_train <- sample(1:n_train, size = n_train, replace = TRUE)
                boot_fit <- npscoef(bw, betas = TRUE,
                                    txdat = cov_train[boot_train,], tydat = Y_train[boot_train], tzdat = mod_train[boot_train,],
                                    exdat = rbind(cov_train, cov_test), eydat = c(Y_train, Y_test), ezdat = rbind(mod_train, mod_test))
                boot_beta[,,b] <- boot_fit$beta
              }
            }
          )["elapsed"]
        }, error = function(e){return(NULL)})
      if(!is.null(boot_time)){
        boot_beta_se <- apply(boot_beta, FUN = sd, MAR = c(1,2))
        
        beta_summary_train <- array(dim = c(n_train, 4, 1 + ncol(cov_train)), dimnames = list(c(), c("MEAN", "SD", "L95", "U95"), c()))
        beta_summary_test <- array(dim = c(n_test, 4, 1 + ncol(cov_test)), dimnames = list(c(), c("MEAN", "SD", "L95", "U95"), c()))
        
        ystar_summary_train <- matrix(nrow = n_train, ncol = 4, dimnames = list(c(), c("MEAN", "SD", "L95", "U95")))
        ystar_summary_test <- matrix(nrow = n_test, ncol = 4, dimnames = list(c(), c("MEAN", "SD", "L95", "U95")))
        
        
        # We can avoid the for loops entirely
        beta_summary_train[,"MEAN",] <- ks_fit$beta[1:n_train,]
        beta_summary_test[,"MEAN",] <- ks_fit$beta[(1 + n_train):(n_test + n_train),]
        
        beta_summary_train[,"SD",] <- boot_beta_se[1:n_train,]
        beta_summary_test[,"SD",] <- boot_beta_se[(1 + n_train):(n_test + n_train),]
        
        beta_summary_train[,"L95",] <- beta_summary_train[,"MEAN",] - qnorm(0.975) * beta_summary_train[,"SD",]
        beta_summary_train[,"U95",] <- beta_summary_train[,"MEAN",] + qnorm(0.975) * beta_summary_train[,"SD",]
        
        beta_summary_test[,"L95",] <- beta_summary_test[,"MEAN",] - qnorm(0.975) * beta_summary_test[,"SD",]
        beta_summary_test[,"U95",] <- beta_summary_test[,"MEAN",] + qnorm(0.975) * beta_summary_test[,"SD",]
        
        ystar_summary_train[,"MEAN"] <- ks_fit$mean[1:n_train]
        rmse_train <- sqrt(mean( (Y_train - ystar_summary_train[,"MEAN"])^2 ))
        
        ystar_summary_train[,"L95"] <- ystar_summary_train[,"MEAN"] - qnorm(0.975) * rmse_train
        ystar_summary_train[,"U95"] <- ystar_summary_train[,"MEAN"] + qnorm(0.975) * rmse_train
        
        ystar_summary_test[,"MEAN"] <- ks_fit$mean[(1 + n_train):(n_test + n_train)]
        ystar_summary_test[,"L95"] <- ystar_summary_test[,"MEAN"] - qnorm(0.975) * rmse_train
        ystar_summary_test[,"U95"] <- ystar_summary_test[,"MEAN"] + qnorm(0.975 ) * rmse_train
        
        results <- list(train = list(beta = beta_summary_train, ystar = ystar_summary_train),
                        test = list(beta = beta_summary_test, ystar = ystar_summary_test),
                        time = bw_time + boot_time)
      } else results <- NULL # closes if/else checking that bootstrapping went fine
    }
    
    
  } else results <- NULL # closes if/else checking that we computed bandwidth fine
  
  return(results)
  
}
