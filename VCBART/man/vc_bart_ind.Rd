\name{vc_bart_ind}
\alias{vc_bart_ind}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to fit a VC-BART model with independent errors
}
\description{
The same function is used to fit a VC-BART model with independent errors in both the cross-sectional and panel settings.
}
\usage{
vc_BART_ind(Y_train, X_train, Z_train, n_vec_train, start_index_train,
X_test, Z_test,  n_vec_test, start_index_test ,
cutpoints, sigma_mu_vec, 
burn, nd, alpha_vec , beta_vec ,nu_sigma = 7, ht_sigma_y = TRUE, ht_sigma_mu = FALSE, verbose = TRUE, print_every = 50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y_train}{A vector of length n containing the training observations}
  \item{X_train}{An n x p matrix of covariates for training observations. Assumes that the first column contains all 1's, for the intercept}
  \item{Z_train}{An n x R matrix of effect modifiers for training observations.}
  \item{n_vec_train}{Vector containing number of obesrvations for each subject. For cross-sectional data, this should be a vector of all 1's.}
  \item{start_index_train}{Vector of indices for first observation of each subject. For cross-sectional models this should be equal to 1:n_train, where n_train = length(Y_train).}
  \item{cutpoints}{List of length R containing all allowable cutpoints for decision tree}
  \item{sigma_mu_vec}{Vector of length p taining scale parameter for prior on jumps. Default is 1/sqrt(m), where m is number of trees}
  \item{burn}{Number of MCMC iteratiosn to discard as burn-in}
  \item{nd}{Number of MCMC samples to output}
  \item{alpha_vec}{Vector of length p containing the parameter alpha in the tree splitting probabilties. Default is a vector with constant entries 0.95}
  \item{beta_vec}{Vector of length p containing the parameter beta in the tree splitting probabilties. Default is a vector with constant entries 2}
  \item{nu_sigma}{Scale parameter for sigma prior}
%%     ~~Describe \code{x} here~~
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}

A list!
\item{f_train_samples}{MCMC samples for regression function evaluated at training observations}
\item{f_test_samples}{MCMC samples for regression function evaluated at testing observations}
\item{beta_train_samples}{MCMC samples of each covariate effect function evaluated at training observations}
\item{beta_test_samples}{MCMC samples of each covaritae effect function evaluated at testing observations}
\item{sigma_samples}{MCMC samples of sigma}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}

