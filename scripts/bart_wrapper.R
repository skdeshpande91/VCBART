# Bart wrapper
library(BART)
source("scripts/bart_summary.R")


bart_wrapper <- function(Y_train, X_train, Z_train, X_test, Z_test){
  
  bart_time <- system.time(
    {
      chain1 <- wbart(cbind(X_train, Z_train), Y_train, cbind(X_test, Z_test), printevery = 5000)
      chain2 <- wbart(cbind(X_train, Z_train), Y_train, cbind(X_test, Z_test), printevery = 5000)
    })["elapsed"]
  
  results <- bart_summary(chain1, chain2)
  results["time"] <- bart_time
  return(results)
}