\name{summarize_posterior_predictive}
\alias{summarize_posterior_predictive}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute posterior preditive mean 95\% posterior predictive interval 
}
\description{
iven two chains (output by two runs of \code{VCBART}), computes the posterior predictive mean and 95\% interval for each new observation. Note these intervals include the observation noise. 
}
\usage{
summarize_posterior_predictive(chain1, chain2, burn)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chain1}{List returned by \code{VCBART} for first MCMC chain}
  \item{chain2}{List returned by \code{VCBART} for second MCMC chain}
  \item{burn}{Number of iterations to be discarded as burn-in}
%%     ~~Describe \code{x} here~~
}

\value{
A list with two named elements, \code{train} and \code{test}, for in-sample and out-of-sample predictions. Each of these elements is a dataframe with 4 columns, one each for the posterior predictive mean, standard deviation, lower bound of 95\% interval, and upper bound of 95\% interval. The rows of \code{train} (resp. \code{test}) corresponds to the rows of the \code{Z_train} (resp. \code{Z_test}) matrices
}
