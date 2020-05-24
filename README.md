# VC-BART: Bayesian trees for varying coefficients

An R package for fitting a linear varying coefficient model using Bayesian Additive Regression Trees.

## To-do (dev)

1. Add sparse version (independent error and compound symmetry error)
2. Change export names (e.g. .sparse_vcbart_ind) and add wrapper function
3. Add comparison to vcrpart package.
4. Adjust rho-sensitivity simulation (only do CS)
5. Re-run p3R2, p5R5 simulation with more modifiers (to show off the sparsity-inducing behavior)


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

The following code chunk sets up a really basic example of `vc_BART_ind`, which fits a VC-BART model with independent errors.

```r
library(VCBART)
library(MASS)

n <- 1 # A single subject
N <- 2500 # total number of observations
p <- 3
R <- 2
sigma <- 1

####################################
# True covariate effect functions
####################################

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


# Add a column for the intercept
X_train <- cbind(rep(1, times = nrow(X_train)), X_train)
X_test <- cbind(rep(1, times = nrow(X_test)), X_test)

cutpoints <- list(seq(0,1,length = 10000),  c(0,1))

chain1 <- vc_BART_ind(Y = Y_train, X_train = X_train, Z_train = Z_train,
                      X_test = X_test, Z_test = Z_test, xinfo_list = cutpoints, verbose = TRUE, print_every = 50)
chain2 <- vc_BART_ind(Y = Y_train, X_train = X_train, Z_train = Z_train,
                      X_test = X_test, Z_test = Z_test, xinfo_list = cutpoints, verbose = TRUE, print_every = 50)

fit_sum <- get_summary(chain1, chain2)

#######
# Make a plot of actual testing responses against posterior predictive mean
plot(Y_test, fit_sum$test$ystar[,"MEAN"])

# Make a plot of actual beta1 values plotted against posterior means (on testing data)
# Note that the indexing is off-set by 1.
plot(beta1_test, fit_sum$test$beta[,"MEAN",2])
```




