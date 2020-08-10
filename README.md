# VC-BART: Bayesian trees for varying coefficients

An R package for fitting a linear varying coefficient model using Bayesian Additive Regression Trees.


For more details about the VC-BART procedure, see [our paper](https://arxiv.org/abs/2003.06416).


This is the development branch.
If you just want to use VC-BART, you should download the package source that is available on the main branch.
A copy of that package source is maintained in the directory VCBART.

The directory VCBARTdev is a prototype package in which updates and other development occurs.
Do not install directly from that source. 




### Details

#### Installation

The package source files are contained in the sub-directory VCBART/ .
To install, you can either download that directory and then build and install the package from the command line (e.g. `R CMD BUILD ...` and `R CMD INSTALL ...`).
You can also install the package using `devtools::install_github` as follows.

```r
library(devtools)
devtools::install_github(repo = "skdeshpande91/VCBART/VCBART")
```

#### Examples


Scripts to reproduce some of the examples in the paper are avilable in the scripts/ directory.

The following code chunk shows how to run VCBART  with independent errors and adaptive split probabilities.

```r
library(VCBART)
library(MASS)

n <- 1 # A single subject
N <- 500 # total number of observations
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

##################################
# Generate X and Z and cutpoints
##################################
set.seed(121719)
Sigma_X <- (0.75)^(abs(outer(1:p, 1:p, FUN = "-")))
X_all <- mvrnorm(N, mu = rep(0, times = p), Sigma = Sigma_X)
Z_all <- matrix(runif(N * R), nrow = N, ncol = R)
Z_all[,2] <- 1 * (Z_all[,2] > 0.5)

beta0_all <- beta0_true(Z_all)
beta1_all <- beta1_true(Z_all)
beta2_all <- beta2_true(Z_all)
beta3_all <- beta3_true(Z_all)

#######################
# Generate outcomes
#######################
mu_all <- beta0_all + X_all[,1] * beta1_all + X_all[,2] * beta2_all + X_all[,3] * beta3_all
eps_all <- rnorm(N, 0, 1)
Y_all <- mu_all + sigma * eps_all

#######################
# Make a training/testing split
# Run VC-BART with independent errors
######################
set.seed(130)
train_index <- sort(sample(1:N, size = floor(0.75 * N), replace = FALSE))
test_index <- (1:N)[which(! (1:N) %in% train_index)]

X_train <- X_all[train_index,]
Y_train <- Y_all[train_index]
Z_train <- as.matrix(Z_all[train_index,], ncol = R)

X_test <- X_all[test_index,]
Y_test <- Y_all[test_index]
Z_test <- as.matrix(Z_all[test_index,], ncol = R)


beta0_train <- beta0_all[train_index]
beta1_train <- beta1_all[train_index]
beta2_train <- beta2_all[train_index]
beta3_train <- beta3_all[train_index]


beta0_test <- beta0_all[test_index]
beta1_test <- beta1_all[test_index]
beta2_test <- beta2_all[test_index]
beta3_test <- beta3_all[test_index]



n_train <- nrow(X_train)
n_test <- nrow(X_test)

cutpoints <- list(seq(0,1, length = 10000), c(0,1))


##########################
# Run 2 chains of VCBART
##########################

chain1 <- VCBART(Y_train, X_train, Z_train, n_train,
                 X_test, Z_test, n_test, cutpoints,
                 intercept = TRUE, M = 50, error_structure = "ind", 
                 split_probs_type = "adaptive", ht_sigma_y = TRUE,
                 nd = 1000, burn = 500, verbose = TRUE, print_every = 50)

chain2 <- VCBART(Y_train, X_train, Z_train, n_train,
                 X_test, Z_test, n_test, cutpoints,
                 intercept = TRUE, M = 50, error_structure = "ind", 
                 split_probs_type = "adaptive", ht_sigma_y = TRUE,
                 nd = 1000, burn = 500, verbose = TRUE, print_every = 50)

################
# Compute posterior mean and credible intervals for beta_j
# Compute posterior predictive mean and 95% predictive intervals
# Perfor modifier selection
################

beta_summary <- summarize_beta(chain1, chain2, burn = 500)
ystar_summary <- summarize_posterior_predictive(chain1, chain2, burn = 500)
beta_support <- get_beta_support(chain1, chain2, burn = burn, max_cutoff = 10)

# Plot the actual out-of-sample observations against the posterior predictive means
plot(Y_test, ystar_summary$test[,"MEAN"])

# Plot the values of beta1 against the posterior means
plot(beta1_test, beta_summary$test[,"MEAN",2])

# Check the median probability models for cut-off of 1
# We select Z_1 and Z_2 for beta0, Z_1 for beta_1, Z_1, for beta_2, and Z_1 & Z_2 for beta_3

beta_support$support[[1]]
```
