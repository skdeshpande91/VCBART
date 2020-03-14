\name{vc_bart_cs}
\alias{vc_bart_cs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to fit a panel VC-BART model with compoundy symmetry-structured errors
}
\description{
Runs a conditional Gibbs sampler to draw the trees and residual standard deviation conditionally on a fixed autocorrelation parameter rho.
Assumes a compound symmetry error correlation structure in which all off-diagonal entries in the subject-specific correlation matrix are equal to rho.
}
\usage{
vc_BART_cs(Y, X_train, Z_train, 
           n_vec_train = c(nrow(X_train)), start_index_vec_train = c(1), 
           X_test, Z_test, n_vec_test = c(nrow(X_train)), start_index_vec_test = c(1), 
           xinfo_list, m = 50, sigma_mu_vec = rep(1/sqrt(m), times = nrow(X_train)), 
           alpha_vec = rep(0.95, times = ncol(X_train)), beta_vec = rep(2, times = ncol(X_train)), 
           burn = 250L, nd = 1000L, 
           rho = 0.5, 
           sigma_hat = 1.0, ht_sigma_y = TRUE, ht_sigma_mu = FALSE, 
           nu_sigma = 7, nu_sigma_mu = 7, variance_prob = 0.9, 
           verbose = FALSE, print_every = 500L)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{Vector containing the training observations.}
  \item{X_train}{Matrix of covariates for training observations. Assumes that the first column contains all 1's, for the intercept.}
  \item{Z_train}{Matrix of effect modifiers for training observations.}
  \item{n_vec_train}{Vector containing number of training obesrvations for each subject. Default is c(nrow(X_train)), corresponding to scenario in which all observations come from a single subject.}
  \item{start_index_train}{Vector of indices for first observation of each subject. Default is c(1), corresponding to scenario in which all observatiosn come from a single subject.}
  \item{X_test}{Matrix of covariates for testing observations. Assumes that the first column contains all 1's, for the intercept.}
  \item{Z_test}{Matrix of effect modifiers for testing observations.}
  \item{n_vec_test}{Vector containing number of testing observations for each subject. Default is c(nrow(X_test)).}
  \item{start_index_test}{Vector of indices for first observation of each subject. Default is c(1).}
  \item{xinfo_list}{List of length R containing all allowable cutpoints for decision tree}
  \item{m}{Number of trees used to approximate each beta_j Default is 50.}
  \item{sigma_mu_vec}{Vector of length p taining scale parameter for prior on jumps. Default is 1/sqrt(m), where m is number of trees}
  \item{alpha_vec}{Vector of length p containing the parameter alpha in the tree splitting probabilties. Default is a vector with constant entries 0.95}
  \item{beta_vec}{Vector of length p containing the parameter beta in the tree splitting probabilties. Default is a vector with constant entries 2}
  \item{burn}{Number of MCMC iteratiosn to discard as burn-in. Default is 1000.}
  \item{nd}{Number of MCMC samples to output. Default is 1000.}
  \item{rho}{Value of autocorrelation parameter rho. Default is 0.5.}
  \item{sigma_hat}{Initial over-estimate of sigma. Default is one. Used only if ht_sigma_y = FALSE.}
  \item{ht_sigma_y}{Logical value indicating whether half-t prior should be used for sigma. Default is TRUE. Otherwise scaled inverse chi-squared prior is used.}
  \item{ht_sigma_mu}{Logical value indicating whether half-t prior should be used for sigma_mu. Default is FALSE.}
  \item{nu_sigma}{Degrees of freedom for sigma prior. Default is 7}
  \item{nu_sigma_mu}{Degrees of freedom for sigma_mu prior. Default is 7.}
  \item{variance_prop}{Amount of prior probability that scaled-inverse chi-square prior on sigma (if used) places on the event sigma < sigma_hat. Default is 0.9. Used only if ht_sigma_y = FALSE}
  \item{verbose}{Logical value indicating whether to print progress. Default is FALSE.}
  \item{print_every}{Number of iterations after which to print MCMC progress.}
%%     ~~Describe \code{x} here~~
}
\details{

}
\value{
A list containing
\item{f_train_samples}{An nrow(X_train) x nd matrix containing posterior samples of the regression function evaluated at training observations.}
\item{f_test_samples}{An nrow(X_train) x nd matrix containing posterior samples of the regression function evaluated at testing observations.}
\item{beta_train_samples}{An nrow(X_train) x ncol(X_train) x nd array containing posterior samples of each covariate effect function evaluated at training observations.}
\item{beta_test_samples}{An nrow(X_test) x ncol(X_test) x nd array containing posterior samples of each covariate effect function evaluated at testing observations.}
\item{sigma_samples}{A vector of length nd containing posterior samples of sigma.}
\item{time}{Run-time (in seconds) of MCMC chain.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Sameer K. Deshpande 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\seealso{\code{\link{vc_bart_ar}} and \code{\link{vc_bart_ind}}}

\examples{

}

