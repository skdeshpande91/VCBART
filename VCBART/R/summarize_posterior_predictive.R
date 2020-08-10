summarize_posterior_predictive <- function(chain1, chain2, burn){
  N_train <- dim(chain1$beta_train_samples)[1]
  N_test <- dim(chain1$beta_test_samples)[1]
  p <- dim(chain1$beta_train_samples)[2]
  
  f_samples_train <- cbind(chain1$f_train_samples, chain2$f_train_samples) #N_train x iter draws 
  f_samples_test <- cbind(chain1$f_test_samples, chain2$f_test_samples)
  sigma_samples <- c(chain1$sigma_samples[-(1:burn)], chain2$sigma_samples[-(1:burn)])
  
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
  return(list(train = ypred_summary_train, test = ypred_summary_test))
}