\name{predict_betas}
\alias{predict_betas}
\title{
Compute posterior predictive evaluates of covariate effect functions.
}
\description{
Given an object returned by \code{VCBART_ind} or \code{VCBART_cs} and matrices of continuous and categorical modifiers, returns MCMC samples of the coefficient functions evaluated the provided points.
}
\usage{
predict_betas(fit,
              Z_cont = matrix(0, nrow = 1, ncol = 1),
              Z_cat =  matrix(0, nrow = 1, ncol = 1),
              verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{A list returned by \code{\link{VCBART_ind}} or \code{\link{VCBART_cs}}}
  \item{Z_cont}{Matrix of continuous modifiers at which you wish to evaluate the covariate effect functions. Default is a 1x1 matrix, which signals that no continuous modifiers are required for these evaluations.}
  \item{Z_cat}{Integer matrix of categorical modifiers at which you wish to evaluate the covariate effect functions. Default is a 1x1 matrix, which signals that no continuous modifiers are required for these evaluations.}
  \item{verbose}{Boolean indicating whether the code should print its progress (\code{TRUE}). Default is \code{TRUE}.}
}

\value{
An array of size nd x N x (p+1) where nd is the total number of MCMC draws, N is the total number of points at which you are evaluating the covariate effect functions (i.e. \code{nrow(Z_cont)} or \code{nrow(Z_cat)}), and p is the number of covariates.
Note that the intercept function is included as the first slice in the third dimension.
}