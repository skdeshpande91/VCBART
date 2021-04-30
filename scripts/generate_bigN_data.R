# Generate data for illustration

library(MASS)
n <- 1 # total number of subjects
N <- 100000 * 1.1
p <- 5
R <- 20

beta0_true <- function(Z){
  return( 3 * Z[,1] + (2 - 5 * (Z[,2] > 0.5)) * sin(pi * Z[,1]) - 2 * (Z[,2] > 0.5))
}
beta1_true <- function(Z){
  tmp_z <- 2 * Z[,1] + 0.5
  return( sin(10 * pi * tmp_z)/(2 * tmp_z) + (tmp_z - 1)^4)
}
beta2_true <- function(Z){
  return( (3 - 3*cos(6*pi*Z[,1]) * Z[,1]^2) * (Z[,1] > 0.6) - (10 * sqrt(Z[,1])) * (Z[,1] < 0.25) )
}
beta3_true <- function(Z){
  return(rep(0, times = nrow(Z)))
}
beta4_true <- function(Z){
  return(10 * sin(pi * Z[,1] * Z[,2]) + 20 * (Z[,3] - 0.5)^2 + 10 * Z[,4] + 5 * Z[,5])
}
beta5_true <- function(Z){
  return(exp(sin((0.9 * (Z[,1] + 0.48))^10)) + Z[,2] * Z[,3] + Z[,4])
}


set.seed(6820)
Sigma_X <-  (0.75)^(abs(outer(1:p, 1:p, FUN = "-")))

X_all <- mvrnorm(N, mu = rep(0, times = p), Sigma = Sigma_X)

Z_all <- matrix(runif(N * R), nrow = N, ncol = R)
Z_all[,2] <- 1*(Z_all[,2] > 0.5)
Z_all[,16] <- 1*(Z_all[,16] > 0.5)
Z_all[,17] <- 1*(Z_all[,17] > 0.5)
Z_all[,18] <- 1*(Z_all[,18] > 0.5)
Z_all[,19] <- 1*(Z_all[,19] > 0.5)
Z_all[,20] <- 1*(Z_all[,20] > 0.5)

beta0_all <- beta0_true(Z_all)
beta1_all <- beta1_true(Z_all)
beta2_all <- beta2_true(Z_all)
beta3_all <- beta3_true(Z_all)
beta4_all <- beta4_true(Z_all)
beta5_all <- beta5_true(Z_all)
##############

beta_all <- cbind(beta0_all, beta1_all, beta2_all, beta3_all, beta4_all, beta5_all)
colnames(beta_all) <- c("beta0", "beta1", "beta2", "beta3", "beta4", "beta5")

mu_all <- beta0_all + X_all[,1] * beta1_all + X_all[,2] * beta2_all +
  X_all[,3] * beta3_all + X_all[,4] * beta4_all + X_all[,5] * beta5_all

Y_all <- mu_all + rnorm(N, 0, 1)

cutpoints <- list()
for(r in 1:R) cutpoints[[r]] <- seq(0, 1, length = 1000)
cutpoints[[2]] <- c(0,1)
cutpoints[[16]] <- c(0,1)
cutpoints[[17]] <- c(0,1)
cutpoints[[18]] <- c(0,1)
cutpoints[[19]] <- c(0,1)
cutpoints[[20]] <- c(0,1)

cov_all <- data.frame(X_all)
colnames(cov_all) <- paste0("X",1:p)
mod_all <- data.frame(Z_all)
mod_all[,2] <- factor(Z_all[,2], levels = c(0,1), ordered = FALSE)
mod_all[,16] <- factor(Z_all[,16], levels = c(0,1), ordered = FALSE)
mod_all[,17] <- factor(Z_all[,17], levels = c(0,1), ordered = FALSE)
mod_all[,18] <- factor(Z_all[,18], levels = c(0,1), ordered = FALSE)
mod_all[,19] <- factor(Z_all[,19], levels = c(0,1), ordered = FALSE)
mod_all[,20] <- factor(Z_all[,20], levels = c(0,1), ordered = FALSE)

colnames(mod_all) <- paste0("Z", 1:R)

save(N, n,  p, R, X_all, Z_all, Y_all, mu_all, beta_all,
     cutpoints, cov_all, mod_all, file = "data/bigN_data.RData")



