\name{VCBART_ind}
\alias{VCBART_ind}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit a VCBART model with independent error structure
}
\description{
Fit a varying coefficient model to panel data. Assumes residual errors are independent within and between subjects.
See Deshpande et al. (2025) for details about the model and MCMC sampler.
}
\usage{
VCBART_ind(Y_train,subj_id_train, ni_train,X_train,
           Z_cont_train = matrix(0, nrow = 1, ncol = 1),
           Z_cat_train = matrix(0L, nrow = 1, ncol = 1),
           X_test = matrix(0, nrow = 1, ncol = 1),
           Z_cont_test = matrix(0, nrow = 1, ncol = 1),
           Z_cat_test = matrix(0, nrow = 1, ncol = 1),
           unif_cuts = rep(TRUE, times = ncol(Z_cont_train)),
           cutpoints_list = NULL,
           cat_levels_list = NULL,
           edge_mat_list = NULL,
           graph_split = rep(FALSE, times = ncol(Z_cat_train)),
           sparse = TRUE,
           M = 50,
           mu0 = NULL, tau = NULL, nu = NULL, lambda = NULL,
           nd = 1000, burn = 1000, thin = 1,
           save_samples = TRUE, save_trees = TRUE,
           verbose = TRUE, print_every = floor( (nd*thin + burn)/10))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y_train}{Vector of continous responses for training data}
  \item{ni_train}{Vector containing the number of observations per subject in the training data.}
  \item{subj_id_train}{Vector of length \code{length(Y_train)} that records which subject contributed each observation. Subjects should be numbered sequentially from 1 to \code{length(ni_train)}.}
  \item{X_train}{Matrix of covariates for training observations. Do not include intercept as the first column.}
  \item{Z_cont_train}{Matrix of continuous modifiers for training data. Note, modifiers must be rescaled to lie in the interval [-1,1]. Default is a 1x1 matrix, which signals that there are no continuous modifiers in the training data.}
  \item{Z_cat_train}{Integer matrix of categorical modifiers for training data. Note categorical levels should be 0-indexed. That is, if a categorical modifier has 10 levels, the values should run from 0 to 9. Default is a 1x1 matrix, which signals that there are no categorical modifiers in the training data.}
  \item{X_test}{Matrix of covariate for testing observations. Default is a 1x1 matrix, which signals that testing data is not provided.}
  \item{Z_cont_test}{Matrix of continuous modifiers for testing data. Default is a 1x1 matrix, which signals that testing data is not provided.}
  \item{Z_cat_test}{Integer matrix of categorical modifiers for testing data. Default is a 1x1 matrix, which signals that testing data is not provided.}
  \item{unif_cuts}{Vector of logical values indicating whether cutpoints for each continuous modifier should be drawn from a continuous uniform distribution (\code{TRUE}) or a discrete set (\code{FALSE}) specified in \code{cutpoints_list}. Default is \code{TRUE} for each variable in \code{Z_cont_train}}
  \item{cutpoints_list}{List of length \code{ncol(Z_cont_train)} containing a vector of cutpoints for each continuous modifier. By default, this is set to \code{NULL} so that cutpoints are drawn uniformly from a continuous distribution.}
  \item{cat_levels_list}{List of length \code{ncol(Z_cat_train)} containing a vector of levels for each categorical modifier. If the j-th categorical modifier contains L levels, \code{cat_levels_list[[j]]} should be the vector \code{0:(L-1)}. Default is \code{NULL}, which corresponds to the case that no categorical modifiers are available.}
  \item{edge_mat_list}{List of adjacency matrices if any of the categorical modifiers are network-structured. Default is \code{NULL}, which corresponds to the case that there are no network-structured categorical modifiers.}
  \item{graph_split}{Vector of logicals indicating whether each categorical modifier is network-structured. Default is \code{rep(FALSE, times = ncol(Z_cat_train))}.}
  \item{sparse}{Logical, indicating whether or not to perform variable selection in each tree ensemble based on a sparse Dirichlet prior rather than uniform prior; see Linero 2018. Default is \code{TRUE}}
  \item{M}{Number of trees in each ensemble. Default is 50.}
  \item{mu0}{Prior mean for jumps/leaf parameters. Default is 0 for each beta function. If supplied, must be a vector of length \code{1 + ncol(X_train)}.}
  \item{tau}{Prior standard deviation for jumps/leaf parameters. Default is \code{1/sqrt(M)} for each beta function. If supplied, must be a vector of length \code{1 + ncol(X_train)}.}
  \item{nu}{Degrees of freedom for scaled-inverse chi-square prior on sigma^2. Default is 3.}
  \item{lambda}{Scale hyperparameter for scaled-inverse chi-square prior on sigma^2. Default places 90\% prior probability that sigma is less than \code{sd(Y_train)}.}
  \item{nd}{Number of posterior draws to return. Default is 1000.}
  \item{burn}{Number of MCMC iterations to be treated as "warmup" or "burn-in". Default is 1000.}
  \item{thin}{Number of post-warmup MCMC iteration by which to thin. Default is 1.}
  \item{save_samples}{Logical, indicating whether to return all posterior samples. Default is \code{TRUE}. If \code{FALSE}, only posterior mean is returned.}
  \item{save_trees}{Logical, indicating whether or not to save a text-based representation of the tree samples. This representation can be passed to \code{predict_flexBART} to make predictions at a later time. Default is \code{FALSE}.}
  \item{verbose}{Logical, inciating whether to print progress to R console. Default is \code{TRUE}.}
  \item{print_every}{As the MCMC runs, a message is printed every \code{print_every} iterations. Default is \code{floor( (nd*thin + burn)/10)} so that only 10 messages are printed.}
}
\details{
Given \eqn{p}{p} covariates \eqn{X_{1}, \ldots, X_{p}}{X1, ..., Xp} and \eqn{r}{r} effect modifiers \eqn{Z_{1}, \ldots, Z_{r}}{Z1, ..., Zr}, the varying coefficient model asserts that

\eqn{E[Y \vert X = x, Z = ] = \beta_0(z) + \beta_1(z) * x_1 + ... \beta_p(z) * X_p.}{E[Y | X = x, Z = ] = beta_0(z) + beta_1(z) * x1 + ... beta_p(z) * xp.}

That is, for any r-vector \eqn{Z}{Z}, the relationships between \eqn{X}{X} and \eqn{Y}{Y} is linear. 
However, the specific relationship is allowed to vary with respect tp \eqn{Z}{Z}.
VCBART approximates the covariate effect functions \eqn{\beta_0(Z), \ldots, \beta_p(Z)}{beta_0, ..., beta_p} using ensembles of regression trees.
This function assumes that the within-subject errors are independent.
}
\value{A list containing
\item{y_mean}{Mean of the training observations (needed by \code{predict_VCBART})}
\item{y_sd}{Standard deviation of the training observations (needed by \code{predict_VCBART})}
\item{x_mean}{Vector of means of columns of \code{X_train}, including the intercept (needed by \code{predict_VCBART}).}
\item{x_sd}{Vector of standard deviations of \code{X_trian}, including the intercept (needed by \code{predict_VCBART}).}
\item{yhat.train.mean}{Vector containing posterior mean of evaluations of regression function E[y|x,z] on training data.}
\item{betahat.train.mean}{Matrix with \code{length(Y_train)} rows and \code{ncol(X_train)+1} columns containing the posterior mean of evaluations of each coefficient function evaluated on the training data. Each row corresponds to a training set observation and each colunn corresponds to a coefficient function. Note the first column is for the intercept function.}
\item{yhat.train}{Matrix with \code{nd} rows and \code{length(Y_train)} columns. Each row corresponds to a posterior sample of the regression function E[y|x,z] and each column corresponds to a training set observation. Only returned if \code{save_samples == TRUE}.}
\item{betahat.train}{Array of dimension with \code{nd} x \code{length(Y_train)} x \code{ncol(X_train)+1} containing posterior samples of evaluations of the coefficient functions. The first dimension corresponds to posterior samples/MCMC iterations, the second dimension corresponds to individual training set observations, and the third dimension corresponds to coefficient functions. Only returned if \code{save_samples == TRUE}.}
\item{yhat.test.mean}{Vector containing posterior mean of evaluations of regression function E[y|x,z] on testing data.}
\item{betahat.test.mean}{Matrix with \code{nrow(X_test)} rows and \code{ncol(X_testn)+1} columns containing the posterior mean of evaluations of each coefficient function evaluated on the training data. Each row corresponds to a training set observation and each colunn corresponds to a coefficient function. Note the first column is for the intercept function.}
\item{yhat.test}{Matrix with \code{nd} rows and \code{nrow(X_test)} columns. Each row corresponds to a posterior sample of the regression function E[y|x,z] and each column corresponds to a testing set observation. Only returned if \code{save_samples == TRUE}.}
\item{betahat.test}{Array of size \code{nd} x \code{nrow(X_test)} x \code{ncol(X_test)+1} containing posterior samples of evaluations of the coefficient functions. The first dimension corresponds to posterior samples/MCMC iterations, the second dimension corresponds to individual training set observations, and the third dimension corresponds to coefficient functions. Only returned if \code{save_samples == TRUE}.}
\item{sigma}{Vector containing ALL samples of the residual standard deviation, including warmup.}
\item{varcounts}{Array of size \code{nd} x R x \code{ncol(X)+1} that counts the number of times a variable was used in a decision rule in each posterior sample of each ensemble. Here R is the total number of potential modifiers (i.e. \code{R = ncol(Z_cont_train) + ncol(Z_cat_train)}).}
\item{theta}{If \code{sparse=TRUE}, an array of size \code{nd} x R \code{ncol(X)+1} containing samples of the variable splitting probabilities.}
\item{trees}{A list (of length \code{nd}) of lists (of length \code{ncol(X_train)+1}) of character vectors (of length \code{M}) containing textual representations of the regression trees. The string for the s-th sample of the m-th tree in the j-th ensemble is contaiend in \code{trees[[s]][[j]][m]}. These strings are parsed by \code{predict_VCBART} to reconstruct the C++ representations of the sampled trees.}
}
\examples{
\dontrun{
############
# True beta functions
beta0_true <- function(Z){
  tmp_Z <- (Z+1)/2
  return( 3 * tmp_Z[,1] + 
  (2 - 5 * (tmp_Z[,2] > 0.5)) * sin(pi * tmp_Z[,1]) - 
  2 * (tmp_Z[,2] > 0.5))
}
beta1_true <- function(Z){
  tmp_Z <- (Z+1)/2
  return(sin(2*tmp_Z[,1] + 0.5)/(4*tmp_Z[,1] + 1) + (2*tmp_Z[,1] - 0.5)^3)
}
beta2_true <- function(Z){
  tmp_Z <- (Z+1)/2
  return( (3 - 3*cos(6*pi*tmp_Z[,1]) * tmp_Z[,1]^2) * (tmp_Z[,1] > 0.6) - 
  (10 * sqrt(tmp_Z[,1])) * (tmp_Z[,1] < 0.25) )
}
beta3_true <- function(Z){
  #return(rep(0, times = nrow(Z)))
  return(rep(1, times = nrow(Z)))
}
beta4_true <- function(Z){
  tmp_Z <- (Z+1)/2
  return(10 * sin(pi * tmp_Z[,1] * tmp_Z[,2]) + 
  20 * (tmp_Z[,3] - 0.5)^2 + 
  10 * tmp_Z[,4] + 5 * tmp_Z[,5])
}
beta5_true <- function(Z){
  tmp_Z <- (Z+1)/2
  return(exp(sin((0.9 * (tmp_Z[,1] + 0.48))^10)) + 
  tmp_Z[,2] * tmp_Z[,3] + tmp_Z[,4])
}

################
# Set up problem dimension
###############
n_train <- 250 
n_test <- 25
sigma <- 1

set.seed(417)
n_all <- n_train + n_test
ni_all <- rep(4, times = n_all) # 4 observations per subject
subj_id_all <- rep(1:n_all, each = 4) # give every subject an id number
N_all <- sum(ni_all) # total number of observations

p <- 5 # number of covariates
R_cont <- 20 # number of continuous modifiers
R_cat <- 0 # number of categorical modifiers
R <- R_cont + R_cat

Sigma_X <-  (0.5)^(abs(outer(1:p, 1:p, FUN = "-"))) # covariates are all correlated

X_all <- MASS::mvrnorm(N_all, mu = rep(0, times = p), Sigma = Sigma_X)
Z_cont_all <- matrix(runif(N_all * R_cont, min = -1, max = 1), nrow = N_all, ncol = R_cont)
beta0_all <- beta0_true(Z_cont_all)
beta1_all <- beta1_true(Z_cont_all)
beta2_all <- beta2_true(Z_cont_all)
beta3_all <- beta3_true(Z_cont_all)
beta4_all <- beta4_true(Z_cont_all)
beta5_all <- beta5_true(Z_cont_all)

beta_all <- cbind(beta0_all, beta1_all, beta2_all, 
beta3_all, beta4_all, beta5_all)

mu_all <- beta0_all + X_all[,1] * beta1_all + X_all[,2] * beta2_all +
  X_all[,3] * beta3_all + X_all[,4] * beta4_all + X_all[,5] * beta5_all
Y_all <- mu_all + sigma * rnorm(n = N_all, mean = 0, sd = 1)


#############
# Create a training/testing split
# since we randomly generated X and Z, 
# we can just use first n_train subjects as train
############
train_subjects <- 1:n_train
test_subjects <- (n_train+1):(n_train + n_test)

test_index <- which(subj_id_all \%in\% test_subjects)
train_index <- which(subj_id_all \%in\% train_subjects)

ni_train <- ni_all[train_subjects]
ni_test <- ni_all[test_subjects]

N_train <- sum(ni_train)
N_test <- sum(ni_test)

X_train <- X_all[train_index,]
X_test <- X_all[test_index,]

Z_cont_train <- Z_cont_all[train_index,]
Z_cont_test <- Z_cont_all[test_index,]
subj_id_train <- rep(1:n_train, times = ni_train)

Y_train <- Y_all[train_index]
Y_test <- Y_all[test_index]

beta_train <- beta_all[train_index,]
beta_test <- beta_all[test_index,]

fit <-
  VCBART::VCBART_ind(Y_train = Y_train,
                     subj_id_train = subj_id_train,
                     ni_train = ni_train,
                     X_train = X_train,
                     Z_cont_train = Z_cont_train,
                     X_test = X_test,
                     Z_cont_test = Z_cont_test,
                     verbose = FALSE)
plot(c(beta_train, beta_test),
     c(fit$betahat.train.mean, fit$betahat.test.mean),
     xlab = "Actual", ylab = "Fitted")
abline(a = 0, b = 1, col = 'blue')
}
}
\references{
Deshpande, S.K, Bai, R., Balocchi, C., Starling, J., and Weiss, J. (2025). VCBART: Bayesian trees for varying coefficients.
\emph{Bayesian Analysis}. \doi{doi:10.1214/24-BA1470}
}