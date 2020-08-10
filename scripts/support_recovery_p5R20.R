load("data/p5R20_data.RData")
source("scripts/assess_support_recovery.R")

support_recovery <- array(dim = c(4, 10, 25), dimnames = list(c("tp", "fp", "fn", "tn"), paste0("cutoff_", 1:10), c()))

get_precision <- function(x){
  prec <- x["tp"]/(x["tp"] + x["fp"])
  return(prec)
}
get_accuracy <- function(x){
  acc <- (x["tp"] + x["tn"])/sum(x)
}
get_sensitivity <- function(x){
  sens <- x["tp"]/(x["tp"] + x["fn"])
}
get_specificity <- function(x){
  spec <- x["tn"]/(x["tn"] + x["fp"])
}
get_f1 <- function(x){
  f1 <- 2 * x["tp"]/(2 * x["tp"] + x["fp"] + x["fn"])
}

precision <- matrix(nrow = 25, ncol = 10, dimnames = list(c(), paste0("cutoff", 1:10)))
accuracy <- matrix(nrow = 25, ncol = 10, dimnames = list(c(), paste0("cutoff", 1:10)))
sensitivty <- matrix(nrow = 25, ncol = 10, dimnames = list(c(), paste0("cutoff", 1:10)))
specifity <- matrix(nrow = 25, ncol = 10, dimnames = list(c(), paste0("cutoff", 1:10)))
f1 <- matrix(nrow = 25, ncol = 10, dimnames = list(c(), paste0("cutoff", 1:10)))

for(sim_number in 1:25){
  if(file.exists(paste0("results/sim_p5R20/vcbart_adapt/vcbart_adapt_", sim_number, ".RData"))){
    load(paste0("results/sim_p5R20/vcbart_adapt/vcbart_adapt_", sim_number, ".RData"))
    tmp_fit <- get(paste0("vcbart_adapt_", sim_number))
    rm(list = paste0("vcbart_adapt_", sim_number))
    
    for(cutoff in 1:10){
      support_recovery[,cutoff,sim_number] <- assess_support_recovery(true_support, tmp_fit$beta_support[["support"]][[cutoff]], R)
    }
    
    precision[sim_number,] <- apply(support_recovery[,,sim_number], MARGIN = 2, FUN = get_precision)
    accuracy[sim_number,] <- apply(support_recovery[,,sim_number], MARGIN = 2, FUN = get_accuracy)
    sensitivty[sim_number,] <- apply(support_recovery[,,sim_number], MARGIN = 2, FUN = get_sensitivity)
    specifity[sim_number,] <- apply(support_recovery[,,sim_number], MARGIN = 2, FUN = get_specificity)
    f1[sim_number,] <- apply(support_recovery[,,sim_number], MARGIN = 2, FUN = get_f1)
  }
}

save(precision, accuracy, sensitivty, specifity, f1, support_recovery, file = "results/sim_p5R20/support_recovery_p5R20.RData")
