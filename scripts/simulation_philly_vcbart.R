# Fit Philadelphia crime data using VC-BART CS



library(VCBART)


load("data/philly_monthly_cc.RData")

p <- ncol(X_train)
R <- ncol(Z_train)

# Script designed to be run on a high-performance computing cluster
# in R batch mode
# If running on personal machine or in interactive mode
# ignore the following two lines and set rho manually.

args <- commandArgs(TRUE)
rho <- as.numeric(args[1])/10

# Add a column of 1's for the intercept
X_train <- cbind(rep(1, times = nrow(X_train)), X_train)
X_test <- cbind(rep(1, times = nrow(X_test)), X_test)

if(rho == 0){
  tmp_chain1 <- vc_BART_ind(Y = Y_train, X_train = X_train, Z_train = Z_train,
                            n_vec_train = n_vec_train, start_index_vec_train = start_index_train,
                            X_test = X_test, Z_test = Z_test, n_vec_test = n_vec_test, 
                            start_index_vec_test = start_index_test,
                            xinfo_list = cutpoints, nd = 1000, burn = 250, verbose = TRUE, print_every = 50)
  
  
  tmp_chain2 <- vc_BART_ind(Y = Y_train, X_train = X_train, Z_train = Z_train,
                            n_vec_train = n_vec_train, start_index_vec_train = start_index_train,
                            X_test = X_test, Z_test = Z_test, n_vec_test = n_vec_test, 
                            start_index_vec_test = start_index_test,
                            xinfo_list = cutpoints, nd = 1000, burn = 250, verbose = TRUE, print_every = 50)
  
  tmp_sum <- get_summary(tmp_chain1, tmp_chain2)
  
} else{
  tmp_chain1 <- vc_BART_cs(Y = Y_train, X_train = X_train, Z_train = Z_train,
                           n_vec_train = n_vec_train, start_index_vec_train = start_index_train,
                           X_test = X_test, Z_test = Z_test, n_vec_test = n_vec_test, 
                           start_index_vec_test = start_index_test,
                           xinfo_list = cutpoints, nd = 1000, burn = 250, rho = rho, verbose = TRUE, print_every = 50)
  tmp_chain2 <- vc_BART_cs(Y = Y_train, X_train = X_train, Z_train = Z_train,
                           n_vec_train = n_vec_train, start_index_vec_train = start_index_train,
                           X_test = X_test, Z_test = Z_test, n_vec_test = n_vec_test, 
                           start_index_vec_test = start_index_test,
                           xinfo_list = cutpoints, nd = 1000, burn = 250, rho = rho, verbose = TRUE, print_every = 50)
  tmp_sum <- get_summary(tmp_chain1, tmp_chain2)
}


assign(paste0("cs_rho", rho*10, "_sum"), tmp_sum)
save(list = paste0("cs_rho", rho*10, "_sum"), file = paste0("results/philly_cs_", rho*10, ".RData"))
