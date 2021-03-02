setwd("~/Google Drive/My Drive/STAT/Research/2VC_Bart/R package/VCBART/")
tmp = load("data/p5R20_rho75_data.RData")

## to try the method with a larger number of per person observations 
## (20 instead of 4, with a total population of 100 instead of 500, so the number of total Y's is the same)
## I ran generate_p5R20_panel_data.R changing "n <- 100" and "n_all <- rep(20, times = n)" together with
## "Sigma_eps <- matrix(rho_true, nrow = 20, ncol = 20) + diag(1 - rho_true, nrow = 20)"
# tmp = load("data/p5R20_rho75_data2.RData")

## I modified functions in VCBARTdev/src (now saved in the repo)
## then compiled the package following instructions here: https://kbroman.org/pkg_primer/pages/build.html
## Specifically: in a terminal, navigate to VCBART, R CMD build VCBARTdev, then R CMD INSTALL VCBARTdev_1.0.tar.gz
## I did the same with the (non modified) VCBART, to compare results on the same data.

## To look at the distribution of rho I did not use the wrapper, but the VCBART function directly
## However, we could modify the wrapper function to output the rho chains too, and call that in a simulation
# source("scripts/vcbart_cs_wrapper.R")
# vcbart_adapt_subset <- 
#   vcbart_cs_wrapper(Y_all, X_all, Z_all, n_all, 
#                  X_all, Z_all, n_all, cutpoints, 
#                  error_structure = "cs",rho_eps = 0.5,
#                  split_probs_type = "fixed", burn = 250, nd = 1000,
#                  verbose = TRUE, print_every = 250)

library(VCBARTdev)

Y_train = Y_all;
X_train = X_all; Z_train = Z_all; n_train = n_all;
X_test = X_all; Z_test = Z_all; n_test = n_all;
error_structure = "cs"
rho_eps = 0.9
split_probs_type = "adaptive"
burn = 250; nd = 1000
verbose = TRUE; print_every = 250

chain1 <- VCBART(Y_train = Y_train,
                 X_train = X_train, Z_train = Z_train, n_train = n_train,
                 X_test = X_test, Z_test = Z_test, n_test = n_test, 
                 cutpoints = cutpoints, 
                 error_structure = error_structure,
                 rho_eps = rho_eps,
                 split_probs_type = split_probs_type,
                 burn = burn, nd = nd, verbose = verbose, print_every = print_every)
plot(chain1$rho_samples, type = "l"); abline(h = rho_true, lty = 2)

mean((rowMeans(chain1$f_test_samples)-Y_train)^2)
mean((rowMeans(chain1$f_test_samples)-mu_all)^2)

colMeans((apply(chain1$beta_train_samples, MARGIN = c(1,2), mean)- beta_all)^2)

## let's compare with the standard VCBART in which we fix rho_eps to 0.9
library(VCBART)
Y_train = Y_all;
X_train = X_all; Z_train = Z_all; n_train = n_all;
X_test = X_all; Z_test = Z_all; n_test = n_all;
error_structure = "cs"
rho_eps = 0.9
split_probs_type = "adaptive"
burn = 250; nd = 1000
verbose = TRUE; print_every = 250

chain1 <- VCBART(Y_train = Y_train,
                 X_train = X_train, Z_train = Z_train, n_train = n_train,
                 X_test = X_test, Z_test = Z_test, n_test = n_test, 
                 cutpoints = cutpoints, 
                 error_structure = error_structure,
                 rho_eps = rho_eps,
                 split_probs_type = split_probs_type,
                 burn = burn, nd = nd, verbose = verbose, print_every = print_every)
plot(chain1$rho_samples, type = "l"); abline(h = rho_true, lty = 2)

mean((rowMeans(chain1$f_test_samples)-Y_train)^2)
mean((rowMeans(chain1$f_test_samples)-mu_all)^2)

colMeans((apply(chain1$beta_train_samples, MARGIN = c(1,2), mean)- beta_all)^2)
