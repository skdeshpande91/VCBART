# Generate panel data for p3R2

library(MASS)
n <- 25 # 25 subjects

ni_train <- 20 # 20 training observations per subject
ni_test <- 5 # 5 testing observations per subject


p <- 3
R <- 2
sigma <- 1

n_vec_train <- rep(ni_train, times = n) # ni_train training observations per subject
n_vec_test <- rep(ni_test, times = n) # ni_test testing observations per subject

N_train <- sum(n_vec_train) # total number of training observations
N_test <- sum(n_vec_test) # total number of testing observations

tmp_train <- cumsum(n_vec_train)
start_index_train <- seq(1, N_train, by = ni_train)

tmp_test <- cumsum(n_vec_test)
start_index_test <- seq(1, N_test, by = ni_test)

############
# Generate X_train, X_test
###########

set.seed(21820)
Sigma_X <-  (0.75)^(abs(outer(1:p, 1:p, FUN = "-")))
X_train <- mvrnorm(N_train, mu = rep(0, times = p), Sigma = Sigma_X)
X_test <- mvrnorm(N_test, mu = rep(0, times = p), Sigma = Sigma_X)

Z_train <- matrix(runif(N_train * R), nrow = N_train, ncol = R)
Z_test <- matrix(runif(N_test * R), nrow = N_test, ncol = R)

Z_train[,2] <- 1 * (Z_train[,2] > 0.5)
Z_test[,2] <- 1 * (Z_test[,2] > 0.5)

X_all <- rbind(X_train, X_test)
Z_all <- rbind(Z_train, Z_test)

cutpoints <- list(seq(0, 1, length = 10000), c(0,1))

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

beta0_all <- beta0_true(Z_all)
beta1_all <- beta1_true(Z_all)
beta2_all <- beta2_true(Z_all)
beta3_all <- beta3_true(Z_all)

beta0_train <- beta0_all[1:N_train]
beta1_train <- beta1_all[1:N_train]
beta2_train <- beta2_all[1:N_train]
beta3_train <- beta3_all[1:N_train]

beta0_test <- beta0_all[(1 + N_train):(N_test + N_train)]
beta1_test <- beta1_all[(1 + N_train):(N_test + N_train)]
beta2_test <- beta2_all[(1 + N_train):(N_test + N_train)]
beta3_test <- beta3_all[(1 + N_train):(N_test + N_train)]

#################

beta_train <- cbind(beta0_train, beta1_train, beta2_train, beta3_train)
beta_test <- cbind(beta0_test, beta1_test, beta2_test, beta3_test)


mu_train <- beta0_train + X_train[,1] * beta1_train + X_train[,2] * beta2_train + X_train[,3] * beta3_train 
mu_test <- beta0_test + X_test[,1] * beta1_test + X_test[,2] * beta2_test + X_test[,3] * beta3_test

save(X_train, X_test, Z_train, Z_test, beta_train, beta_test, mu_train, mu_test,
     n, p, R, ni_train, ni_test, N_train, N_test, n_vec_train, n_vec_test, start_index_train, start_index_test,
     sigma, cutpoints,
     file = "data/p3R2_panel_data.RData")
