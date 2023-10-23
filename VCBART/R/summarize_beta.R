summarize_beta <- function(beta_samples){
  nd <- dim(beta_samples)[1]
  n <- dim(beta_samples)[2]
  p <- dim(beta_samples)[3]
  
  beta_sum <- array(dim = c(n, 3, p), dimnames = list(c(), c("MEAN", "L95", "U95"), c()))
  for(j in 1:p){
    tmp_betas <- beta_samples[,,j]
    # tmp_betas is nd x n. we want column-wise means & quantiles
    beta_sum[,"MEAN",j] <- apply(tmp_betas, MARGIN = 2, FUN = mean)
    beta_sum[,"L95",j] <- apply(tmp_betas, MARGIN = 2, FUN = stats::quantile, probs = 0.025)
    beta_sum[,"U95",j] <- apply(tmp_betas, MARGIN = 2, FUN = stats::quantile, probs = 0.975)
  }
  return(beta_sum)
}