\name{get_summary}
\alias{get_summary}
\title{
Function to combine and summarize two chains of MCMC output for the VC-BART model.
}
\description{
Take the output from multiple runs of each of \code{\link{vc_bart_ind}}, \code{\link{vc_bart_ar}}, and \code{\link{vc_bart_cs}} and summarizes the results.
Specifically, it computes the posterior mean, standard deviation, and 95\% credible interval for the covariate effect functions and the regression function.
It further computes the posterior predictive mean, standard deviation, and 95\% predictive interval for new testing observations.
It also returns the sigma samples.
}
\usage{
get_summary(chain1, chain2, burn = 250)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chain1}{Output from the first run of \code{\link{vc_bart_ind}}, \code{\link{vc_bart_ar}}, or \code{\link{vc_bart_cs}}}
  \item{chain2}{Output from the second run of \code{\link{vc_bart_ind}}, \code{\link{vc_bart_ar}}, or \code{\link{vc_bart_cs}}}
  \item{burn}{Number of burn-in iterations. Default is 250.}
}
\details{
}
\value{
A list with three named elements, ``train'', ``test'', and ``sigma.''. The last, ``sigma'', contains the post-burn-in samples of sigma. Each of ``train'' and ``test'' are themselves lists with the following elements:
\item{fit}{A matrix containing the mean, standard deviation, lower 2.5\% quantile, and upper 97.5\% quantile of the regression function.}
\item{beta}{An array containing the mean, standard deviation, lower 2.5\% quantile, and upper 96.5\% quantile of each covariate effect function. The second dimension of the array corresponds to the different functions}
\item{ystar}{A matrix containing the mean, standard deviation, lower 2.5\% quantile, and upper 97.5\% quantile of the posterior predictive distribution.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Sameer K. Deshpande
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x) 
{
  }
}
