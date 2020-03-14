get_summary <- function(chain1, chain2, burn = 250){
  N_train <- dim(chain1$beta_train_samples)[1]
  N_test <- dim(chain1$beta_test_samples)[1]
  
  p <- dim(chain1$beta_train_samples)[2]
  
  beta_summary_train <- array(dim = c(N_train, 4, p), dimnames = list(c(), c("MEAN", "SD", "L95", "U95"), c()))
  beta_summary_test <- array(dim = c(N_test, 4, p), dimnames = list(c(), c("MEAN", "SD", "L95", "U95"), c()))
  
  sigma_summary <- cbind("Chain1" = chain1$sigma_samples, "Chain2" = chain2$sigma_samples)

  for(k in 1:p){
    beta_summary_train[,"MEAN",k] <- apply(cbind(chain1$beta_train_samples[,k,], chain2$beta_train_samples[,k,]), FUN = mean, MARGIN = 1, na.rm = TRUE)
    beta_summary_train[,"SD",k] <- apply(cbind(chain1$beta_train_samples[,k,], chain2$beta_train_samples[,k,]), FUN = sd, MARGIN = 1, na.rm = TRUE)
    beta_summary_train[,"L95", k] <- apply(cbind(chain1$beta_train_samples[,k,], chain2$beta_train_samples[,k,]), FUN = quantile, MARGIN = 1, probs = 0.025)
    beta_summary_train[,"U95", k] <- apply(cbind(chain1$beta_train_samples[,k,], chain2$beta_train_samples[,k,]), FUN = quantile, MARGIN = 1, probs = 0.975)
    
    beta_summary_test[,"MEAN",k] <- apply(cbind(chain1$beta_test_samples[,k,], chain2$beta_test_samples[,k,]), FUN = mean, MARGIN = 1, na.rm = TRUE)
    beta_summary_test[,"SD",k] <- apply(cbind(chain1$beta_test_samples[,k,], chain2$beta_test_samples[,k,]), FUN = sd, MARGIN = 1, na.rm = TRUE)
    beta_summary_test[,"L95", k] <- apply(cbind(chain1$beta_test_samples[,k,], chain2$beta_test_samples[,k,]), FUN = quantile, MARGIN = 1, probs = 0.025)
    beta_summary_test[,"U95", k] <- apply(cbind(chain1$beta_test_samples[,k,], chain2$beta_test_samples[,k,]), FUN = quantile, MARGIN = 1, probs = 0.975)
  }
  
  # Now get the posterior summaries of the conditional mean function E[y | x,z]

  fit_summary_train <- matrix(nrow = N_train, ncol = 4, dimnames = list(c(), c("MEAN", "SD", "L95", "U95")))
  fit_summary_test <- matrix(nrow = N_test, ncol = 4, dimnames = list(c(), c("MEAN", "SD", "L95", "U95")))
  
  fit_summary_train[,"MEAN"] <- apply(cbind(chain1$f_train_samples, chain2$f_train_samples), FUN = mean, MARGIN = 1, na.rm = TRUE)
  fit_summary_train[,"SD"] <- apply(cbind(chain1$f_train_samples, chain2$f_train_samples), FUN = sd, MARGIN = 1, na.rm = TRUE)
  fit_summary_train[,"L95"] <- apply(cbind(chain1$f_train_samples, chain2$f_train_samples), FUN = quantile, MARGIN = 1, probs = 0.025)
  fit_summary_train[,"U95"] <- apply(cbind(chain1$f_train_samples, chain2$f_train_samples), FUN = quantile, MARGIN = 1, probs = 0.975)
  
  fit_summary_test[,"MEAN"] <- apply(cbind(chain1$f_test_samples, chain2$f_test_samples), FUN = mean, MARGIN = 1, na.rm = TRUE)
  fit_summary_test[,"SD"] <- apply(cbind(chain1$f_test_samples, chain2$f_test_samples), FUN = sd, MARGIN = 1, na.rm = TRUE)
  fit_summary_test[,"L95"] <- apply(cbind(chain1$f_test_samples, chain2$f_test_samples), FUN = quantile, MARGIN = 1, probs = 0.025)
  fit_summary_test[,"U95"] <- apply(cbind(chain1$f_test_samples, chain2$f_test_samples), FUN = quantile, MARGIN = 1, probs = 0.975)
  
  
  
  # Now we need to worry about the posterior predictive distribution for y*
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
  

  results <- list(train = list(fit = fit_summary_train, beta = beta_summary_train, ystar = ypred_summary_train),
                  test = list(fit = fit_summary_test, beta = beta_summary_test, ystar = ypred_summary_test),
                  sigma = sigma_summary)
  return(results)
 
}
