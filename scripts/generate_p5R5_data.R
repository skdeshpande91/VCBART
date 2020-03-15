# Generate the p5R5 example for the bakeoff in the paper
library(MASS)
n <- 1 # A single subject
N <- 2500 
p <- 5
R <- 5
sigma <- 1

beta0_true <- function(Z){
  return( 3 * Z[,1] + (2 - 5 * (Z[,2] > 0.5)) * sin(pi * Z[,1]) - 2 * (Z[,2] > 0.5))
  
}
beta1_true <- function(Z){
  k_se <- 2 * exp(-1/(2 *0.05 * 0.05) * outer(Z[,1], Z[,1], FUN = "-") * outer(Z[,1], Z[,1], FUN = "-"))
  k_per <- exp(-2/(0.1*0.1) * sin(pi*outer(Z[,1], Z[,1], FUN = "-")/2) * sin(pi*outer(Z[,1], Z[,1], FUN = "-")/2))
  k <- k_se * k_per
  return(mvrnorm(n = 1, mu = rep(0, times = length(Z[,1])), Sigma = k))
}
beta2_true <- function(Z){
  return( (3 - 3*cos(6*pi*Z[,1]) * Z[,1]^2) * (Z[,1] > 0.6) - (10 * sqrt(Z[,1])) * (Z[,1] < 0.25) )
}
beta3_true <- function(Z){
  (1 - Z[,3])^(1/3) * sin(3 *pi*(1 - Z[,4])) - sqrt(1 - Z[,1])
}
beta4_true <- function(Z){
  return(10 * sin(pi * Z[,1] * Z[,2]) + 20 * (Z[,3] - 0.5)^2 + 10 * Z[,4] + 5 * Z[,5])
}
beta5_true <- function(Z){
  return(exp(sin((0.9 * (Z[,1] + 0.48))^10)) + Z[,2] * Z[,3] + Z[,4])
}

set.seed(121719)
Sigma_X <- (0.75)^(abs(outer(1:p, 1:p, FUN = "-")))
X_all <- mvrnorm(N, mu = rep(0, times = p), Sigma = Sigma_X)
Z_all <- matrix(runif(N * R), nrow = N, ncol = R)

X_all <- mvrnorm(N, mu = rep(0, times = p), Sigma = Sigma_X)
Z_all <- matrix(runif(N * R), nrow = N, ncol = R)

beta0_all <- beta0_true(Z_all)
beta1_all <- beta1_true(Z_all)
beta2_all <- beta2_true(Z_all)
beta3_all <- beta3_true(Z_all)
beta4_all <- beta4_true(Z_all)
beta5_all <- beta5_true(Z_all)

mu_all <- beta0_all + X_all[,1] * beta1_all + X_all[,2] * beta2_all + 
  X_all[,3] * beta3_all + X_all[,4] * beta4_all + X_all[,5] * beta5_all

eps_all <- rnorm(N, 0, 1)
Y_all <- mu_all + sigma * eps_all

cutpoints <- list()
for(r in 1:R) cutpoints[[r]] <- seq(0, 1, length = 10000)

save(N, p, R, X_all, Z_all, Y_all,
     beta0_all, beta1_all, beta2_all, beta3_all, beta4_all, beta5_all,
     mu_all, eps_all, cutpoints, file = "data/p5R5_data.RData")
