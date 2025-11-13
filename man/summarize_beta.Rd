\name{summarize_beta}
\alias{summarize_beta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute posterior mean and 95\% credible interval for evaluations of each coefficient function.
}
\description{
Given an array of posterior samples of coefficient function evaluations, returns the posterior mean and 95\% credible interval for each evaluation.
}
\usage{
summarize_beta(beta_samples)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta_samples}{An array, returned by \code{VCBART_ind}, \code{VCBART_cs}, or \code{predict_betas} of posterior samples of coefficient function evaluations}
}

\value{
An array of size N x 3 x p where N is the number of inputs at which the coefficient functions are evaluated (i.e. N = \code{dim(beta_samples)[2]}) and p is the total number of coefficient functions including the intercept (i.e. p = \code{dim(beta_samples)[3]}). The j-th slice is an N x 3 matrix where the columns correspond to the posterior mean, 2.5\% quantile, and 97.5\% quantile of each evaluation of the (j-1)-th coefficient function. Note the effect of predictor X_j is the (j+1)-st coefficient function.
}