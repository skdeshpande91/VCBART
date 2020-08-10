\name{summarize_beta}
\alias{summarize_beta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute posterior mean and 95\% credible interval for covariate effects
}
\description{
Given two chains (output by two runs of \code{VCBART}), computes the posterior mean and 95\% credible intervals for each \eqn{\beta_{j}(z)}. 
}
\usage{
summarize_beta(chain1, chain2, burn)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chain1}{List returned by \code{VCBART} for first MCMC chain}
  \item{chain2}{List returned by \code{VCBART} for second MCMC chain}
  \item{burn}{Number of iterations to be discarded as burn-in}
%%     ~~Describe \code{x} here~~
}


\value{
A list with the following elements
  \item{exceed_cutoff_probs}{An R x p x max_cutoff dimensional array. The (r,j) entry in slice k reports the estimated posterior probability that Z_r is selected at least k times in the ensemble for beta_j}
  \item{support}{A list of length K, whose k-th element is a list of length p that lists all indices r such that the (r,j,k) entry of \code{exceed_cutoff_probs} is greater than 0.5. These are the analogs of the median probability model}
}