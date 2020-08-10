# Summarize the beta fits
summarize_beta <- function(chain1, chain2, burn){
  N_train <- dim(chain1$beta_train_samples)[1]
  N_test <- dim(chain1$beta_test_samples)[1]
  
  p <- dim(chain1$beta_train_samples)[2]
  
  beta_summary_train <- array(dim = c(N_train, 4, p), dimnames = list(c(), c("MEAN", "SD", "L95", "U95"), c()))
  beta_summary_test <- array(dim = c(N_test, 4, p), dimnames = list(c(), c("MEAN", "SD", "L95", "U95"), c()))
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
  return(list(train = beta_summary_train, test = beta_summary_test))
}


