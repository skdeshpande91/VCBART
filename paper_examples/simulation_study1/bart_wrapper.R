# Bart wrapper

bart_wrapper <- function(Y_train, X_train, Z_train, X_test, Z_test, nd = 1000, burn = 1000, thin = 1){
  
  bart_time <- system.time(
    bart_fit <- BART::wbart(x.train = cbind(X_train, Z_train), 
                            y.train = Y_train, 
                            x.test = cbind(X_test, Z_test),
                            ndpost = nd, nskip = burn, keepevery = thin,
                            printevery=nd*thin+burn+1))["elapsed"]
  sigma_samples <- bart_fit$sigma[-(1:burn)]
  f_samples_train <- bart_fit$yhat.train
  f_samples_test <- bart_fit$yhat.test
  
  # Summarize posterior over regression function
  fit_summary_train <- matrix(nrow = N_train, ncol = 3, dimnames = list(c(), c("MEAN","L95", "U95")))
  fit_summary_test <- matrix(nrow = N_test, ncol = 3, dimnames = list(c(), c("MEAN", "L95", "U95")))
  
  fit_summary_train[,"MEAN"] <- apply(f_samples_train, FUN = mean, MARGIN = 2, na.rm = TRUE)
  fit_summary_train[,"L95"] <- apply(f_samples_train, FUN = quantile, MARGIN = 2, probs = 0.025)
  fit_summary_train[,"U95"] <- apply(f_samples_train, FUN = quantile, MARGIN = 2, probs = 0.975)
  
  fit_summary_test[,"MEAN"] <- apply(f_samples_test, FUN = mean, MARGIN = 2, na.rm = TRUE)
  fit_summary_test[,"L95"] <- apply(f_samples_test, FUN = quantile, MARGIN = 2, probs = 0.025)
  fit_summary_test[,"U95"] <- apply(f_samples_test, FUN = quantile, MARGIN = 2, probs = 0.975)
  
  # Summarize posterior predictive
  ystar_summary_train <- matrix(nrow = N_train, ncol = 3, dimnames = list(c(), c("MEAN", "L95", "U95")))
  ystar_summary_test <- matrix(nrow = N_test, ncol = 3, dimnames = list(c(), c("MEAN", "L95", "U95")))
  
  N_train <- ncol(f_samples_train)
  N_test <- ncol(f_samples_test)
  
  for(i in 1:N_train){
    tmp_ystar_train <- f_samples_train[,i] + sigma_samples * rnorm(nrow(f_samples_train), 0, 1)
    ystar_summary_train[i,"MEAN"] <- mean(tmp_ystar_train)
    ystar_summary_train[i,"L95"] <- quantile(tmp_ystar_train, probs = 0.025)
    ystar_summary_train[i,"U95"] <- quantile(tmp_ystar_train, probs = 0.975)
  }
  
  for(i in 1:N_test){
    tmp_ystar_test <- f_samples_test[,i] + sigma_samples * rnorm(nrow(f_samples_test), 0, 1)
    ystar_summary_test[i,"MEAN"] <- mean(tmp_ystar_test)
    ystar_summary_test[i,"L95"] <- quantile(tmp_ystar_test, probs = 0.025)
    ystar_summary_test[i,"U95"] <- quantile(tmp_ystar_test, probs = 0.975)
  }
  
  return(list(time = bart_time, train_time = bart_time,
              train = list(fit = fit_summary_train, ystar = ystar_summary_train),
              test = list(fit = fit_summary_test, ystar = ystar_summary_test)))
}