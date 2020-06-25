# Generate data for illustration

library(MASS)
N <- 2500 # total number of subjects
p <- 3
R <- 10
sigma <- 1

set.seed(21820)
Sigma_X <-  (0.75)^(abs(outer(1:p, 1:p, FUN = "-")))

##############
# Define the true covariate effect functions
##############
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

# Make a list letting us which variables
# each beta function depends on
true_support <- list()
true_support[[1]] <- as.integer(c(1,2))  # intercept beta0 depends on Z1, Z2
true_support[[2]] <- as.integer(c(1))
true_support[[3]] <- as.integer(c(1))
true_support[[4]] <- integer(0)


#############
# Generate X, Z, etc.
X_all <- mvrnorm(N, mu = rep(0, times = p), Sigma = Sigma_X)
Z_all <- matrix(runif(N* R), nrow = N, ncol = R)
Z_all[,2] <- 1*(Z_all[,2] > 0.5)


beta0_all <- beta0_true(Z_all)
beta1_all <- beta1_true(Z_all)
beta2_all <- beta2_true(Z_all)
beta3_all <- beta3_true(Z_all)

##############

beta_all <- cbind(beta0_all, beta1_all, beta2_all, beta3_all)

mu_all <- beta0_all + X_all[,1] * beta1_all + X_all[,2] * beta2_all + X_all[,3] * beta3_all

Y_all <- mu_all + sigma * rnorm(N, 0, 1)

cutpoints <- list()
for(r in 1:R) cutpoints[[r]] <- seq(0, 1, length = 1000)
cutpoints[[2]] <- c(0,1)

##########
# In order to run kernel smoothing
# with the appropriate kernels
# we need to convert Z[,2] to a factor
###########
cov_all <- data.frame(X_all)
colnames(cov_all) <- paste0("X",1:p)
mod_all <- data.frame(Z_all)
mod_all[,2] <- factor(Z_all[,2], levels = c(0,1), ordered = FALSE)
colnames(mod_all) <- paste0("Z", 1:R)

save(N, p, R, X_all, Z_all, Y_all, cov_all, mod_all, mu_all, beta_all, cutpoints, true_support, file = "data/p3R10_data.RData")



