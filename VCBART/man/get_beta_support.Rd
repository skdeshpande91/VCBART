\name{get_beta_support}
\alias{get_beta_support}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to perform modifier selection.
}
\description{
Given two chains (output by two runs of \code{VCBART}), computes the posterior probability that each modifier is selected at least K times in each tree ensemble.
Then it reports the modifiers for which this selection probability exceeds 0.5.
When K = 1, the returned set of modifiers can be viewed as an analog of a median probability model. 

}
\usage{
get_beta_support(chain1, chain2, burn = burn, max_cutoff = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chain1}{List returned by \code{VCBART} for first MCMC chain}
  \item{chain2}{List returned by \code{VCBART} for second MCMC chain}
  \item{burn}{Number of iterations to be discarded as burn-in}
  \item{max_cutoff}{Largest K to consider. Default is 10.}
}