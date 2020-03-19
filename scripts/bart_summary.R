# Summarize the output of two BART chains
bart_summary <- function(chain1, chain2)
{
  sigma_samples <- c(chain1$sigma[-(1:100)], chain2$sigma[-(1:100)]) # remove the first 100 burn-in samples
  f_samples_train <- cbind(t(chain1$yhat.train), t(chain2$yhat.train))
  f_samples_test <- cbind(t(chain1$yhat.test), t(chain2$yhat.test))
  
  N_train <- nrow(f_samples_train)
  N_test <- nrow(f_samples_test)
  
  
  fit_summary_train <- matrix(nrow = N_train, ncol = 4, dimnames = list(c(), c("MEAN", "SD", "L95", "U95")))
  fit_summary_test <- matrix(nrow = N_test, ncol = 4, dimnames = list(c(), c("MEAN", "SD", "L95", "U95")))
  
  fit_summary_train[,"MEAN"] <- apply(f_samples_train, FUN = mean, MARGIN = 1, na.rm = TRUE)
  fit_summary_train[,"SD"] <- apply(f_samples_train, FUN = sd, MARGIN = 1, na.rm = TRUE)
  fit_summary_train[,"L95"] <- apply(f_samples_train, FUN = quantile, MARGIN = 1, probs = 0.025)
  fit_summary_train[,"U95"] <- apply(f_samples_train, FUN = quantile, MARGIN = 1, probs = 0.975)
  
  fit_summary_test[,"MEAN"] <- apply(f_samples_test, FUN = mean, MARGIN = 1, na.rm = TRUE)
  fit_summary_test[,"SD"] <- apply(f_samples_test, FUN = sd, MARGIN = 1, na.rm = TRUE)
  fit_summary_test[,"L95"] <- apply(f_samples_test, FUN = quantile, MARGIN = 1, probs = 0.025)
  fit_summary_test[,"U95"] <- apply(f_samples_test, FUN = quantile, MARGIN = 1, probs = 0.975)
  
  
  
  ypred_summary_train <- matrix(nrow = N_train, ncol = 4, dimnames = list(c(), c("MEAN", "SD", "L95", "U95")))
  ypred_summary_test <- matrix(nrow = N_test, ncol = 4, dimnames = list(c(), c("MEAN", "SD", "L95", "U95")))
  
  
  for(i in 1:N_train){
    tmp_ypred_train <- f_samples_train[i,] + sigma_samples * rnorm(ncol(f_samples_train), 0, 1)
    ypred_summary_train[i,"MEAN"] <- mean(tmp_ypred_train)
    ypred_summary_train[i,"SD"] <- sd(tmp_ypred_train)
    ypred_summary_train[i,"L95"] <- quantile(tmp_ypred_train, probs = 0.025)
    ypred_summary_train[i,"U95"] <- quantile(tmp_ypred_train, probs = 0.975)
  }
  
  for(i in 1:N_test){
    tmp_ypred_test <- f_samples_test[i,] + sigma_samples * rnorm(ncol(f_samples_test), 0, 1)
    ypred_summary_test[i,"MEAN"] <- mean(tmp_ypred_test)
    ypred_summary_test[i,"SD"] <- sd(tmp_ypred_test)
    ypred_summary_test[i,"L95"] <- quantile(tmp_ypred_test, probs = 0.025)
    ypred_summary_test[i,"U95"] <- quantile(tmp_ypred_test, probs = 0.975)
  }
  
  return(list(train = list(fit = fit_summary_train, ystar = ypred_summary_train),
              test = list(fit = fit_summary_test, ystar = ypred_summary_test)))
  
}
