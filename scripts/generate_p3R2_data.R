#####
# Generate the data used for the p3R2 cross-sectional simulation setting.
# Saves an .RData object to the data/ sub-directory

library(MASS)
n <- 1 # A single subject
N <- 2500 # total number of subjects
p <- 3
R <- 2
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
  return(rep(0, times = nrow(Z)))
}

set.seed(121719)
Sigma_X <- (0.75)^(abs(outer(1:p, 1:p, FUN = "-")))
X_all <- mvrnorm(N, mu = rep(0, times = p), Sigma = Sigma_X)
Z_all <- matrix(runif(N * R), nrow = N, ncol = R)
Z_all[,2] <- 1 * (Z_all[,2] > 0.5)

beta0_all <- beta0_true(Z_all)
beta1_all <- beta1_true(Z_all)
beta2_all <- beta2_true(Z_all)
beta3_all <- beta3_true(Z_all)

### Now generate outcomes
mu_all <- beta0_all + X_all[,1] * beta1_all + X_all[,2] * beta2_all + X_all[,3] * beta3_all
eps_all <- rnorm(N, 0, 1)
Y_all <- mu_all + sigma * eps_all

cutpoints <- list(seq(0,1,length = 10000),  c(0,1))

save(N, p, R, X_all, Z_all, Y_all,
     beta0_all, beta1_all, beta2_all, beta3_all,
     mu_all, eps_all, cutpoints, file = "data/p3R2_data.RData")

