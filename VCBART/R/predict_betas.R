predict_betas <- function(fit, 
                          Z_cont = matrix(0, nrow = 1, ncol = 1),
                          Z_cat =  matrix(0, nrow = 1, ncol = 1),
                          verbose = TRUE)
{
  p <- length(fit$x_mean)
  M <- length(fit$trees[[1]][[1]])
  
  tmp <- .predict_vcbart(tree_draws = fit$trees,
                         p = p, M = M,
                         tZ_cont = t(Z_cont), tZ_cat = t(Z_cat), verbose = verbose)
  out <- rescale_beta(tmp, fit$y_mean, fit$y_sd, fit$x_mean, fit$x_sd)
  return(out)
}