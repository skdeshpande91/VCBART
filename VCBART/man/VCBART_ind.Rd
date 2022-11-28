\name{VCBART_ind}
\alias{VCBART_ind}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit a VCBART model with independent error structure
}
\description{
Fit a varying coefficient model to panel data. Assumes residual errors are independent within and between subjects.
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
           sparse = TRUE,
           M = 200,
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
  \item{sparse}{Logical, indicating whether or not to perform variable selection in each tree ensemble based on a sparse Dirichlet prior rather than uniform prior; see Linero 2018. Default is \code{TRUE}}
  \item{M}{Number of trees in each ensemble. Default is 200.}
  \item{nd}{Number of posterior draws to return. Default is 1000.}
  \item{burn}{Number of MCMC iterations to be treated as "warmup" or "burn-in". Default is 1000.}
  \item{thin}{Number of post-warmup MCMC iteration by which to thin. Default is 1.}
  \item{save_samples}{Logical, indicating whether to return all posterior samples. Default is \code{TRUE}. If \code{FALSE}, only posterior mean is returned.}
  \item{save_trees}{Logical, indicating whether or not to save a text-based representation of the tree samples. This representation can be passed to \code{predict_flexBART} to make predictions at a later time. Default is \code{FALSE}.}
  \item{verbose}{Logical, inciating whether to print progress to R console. Default is \code{TRUE}.}
  \item{print_every}{As the MCMC runs, a message is printed every \code{print_every} iterations. Default is \code{floor( (nd*thin + burn)/10)} so that only 10 messages are printed.}
}
\details{
Given p covariates X1, ..., Xp and r modifiers Z1, ..., Zr, the varying coefficient model asserts that

E[Y | X = x, Z = ] = beta_0(z) + beta_1(z) * x1 + ... beta_p(z) * xp.

That is, for any r-vector z, the relationships between X and Y is linear. However, the specific relationship is allowed to vary across z's.
VCBART attempts to learn the covariate effect functions beta_0(*) ... beta_p(*) using ensembles of regression trees.
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
\item{trees}{A list (of length \code{nd}) of lists (of length \code{ncol(X_train)+1}) of character vectors (of length \code{M}) containing textual representations of the regression trees. The string for the s-th sample of the m-th tree in the j-th ensemble is contaiend in \code{trees[[s]][[j]][m]}. These strings are parsed by \code{predict_VCBART} to reconstruct the C++ representations of the sampled trees.}
}
\examples{
}