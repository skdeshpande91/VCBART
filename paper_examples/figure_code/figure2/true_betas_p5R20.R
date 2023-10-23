
beta0_true <- function(Z){
  tmp_Z <- (Z+1)/2
  return( 3 * tmp_Z[,1] + (2 - 5 * (tmp_Z[,2] > 0.5)) * sin(pi * tmp_Z[,1]) - 2 * (tmp_Z[,2] > 0.5))
}
beta1_true <- function(Z){
  tmp_Z <- (Z+1)/2
  return(sin(2*tmp_Z[,1] + 0.5)/(4*tmp_Z[,1] + 1) + (2*tmp_Z[,1] - 0.5)^3)
}
beta2_true <- function(Z){
  tmp_Z <- (Z+1)/2
  return( (3 - 3*cos(6*pi*tmp_Z[,1]) * tmp_Z[,1]^2) * (tmp_Z[,1] > 0.6) - (10 * sqrt(tmp_Z[,1])) * (tmp_Z[,1] < 0.25) )
}
beta3_true <- function(Z){
  #return(rep(0, times = nrow(Z)))
  return(rep(1, times = nrow(Z)))
}
beta4_true <- function(Z){
  tmp_Z <- (Z+1)/2
  return(10 * sin(pi * tmp_Z[,1] * tmp_Z[,2]) + 20 * (tmp_Z[,3] - 0.5)^2 + 10 * tmp_Z[,4] + 5 * tmp_Z[,5])
}
beta5_true <- function(Z){
  tmp_Z <- (Z+1)/2
  return(exp(sin((0.9 * (tmp_Z[,1] + 0.48))^10)) + tmp_Z[,2] * tmp_Z[,3] + tmp_Z[,4])
}

true_support <- list()
true_support[[1]] <- as.integer(c(1,2))
true_support[[2]] <- as.integer(c(1))
true_support[[3]] <- as.integer(c(1))
true_support[[4]] <- integer(0)
true_support[[5]] <- as.integer(c(1,2,3,4,5))
true_support[[6]] <- as.integer(c(1,2,3,4))