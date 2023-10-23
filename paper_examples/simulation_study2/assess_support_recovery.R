get_confusion_matrix <- function(true_support, est_support, R){
  if(length(true_support) != length(est_support)) stop("True support and estimated support must be of same length")
  else{
    
    
    support_perf <- matrix(nrow = length(true_support), ncol = 4,
                           dimnames = list(c(), c("tp", "fp", "fn", "tn")))
    
    for(k in 1:length(true_support)){
      support_perf[k,"tp"] <- sum(est_support[[k]] %in% true_support[[k]])
      support_perf[k, "fp"] <- sum(!est_support[[k]] %in% true_support[[k]])
      support_perf[k, "fn"] <- sum(!true_support[[k]] %in% est_support[[k]])
      support_perf[k, "tn"] <- R - support_perf[k, "fn"] - length(est_support[[k]])
    }
    total_perf <- colSums(support_perf)
    total_perf["sen"] <- total_perf["tp"]/(total_perf["tp"] + total_perf["fn"])
    total_perf["spec"] <- total_perf["tn"]/(total_perf["tn"] + total_perf["fp"])
    total_perf["prec"] <- total_perf["tp"]/(total_perf["tp"] + total_perf["fp"])
    total_perf["acc"] <- (total_perf["tp"] + total_perf["tn"])/(total_perf["tp"] + total_perf["fp"] + total_perf["fn"] + total_perf["tn"])
    
    mcc_num <- (total_perf["tp"] * total_perf["tn"] - total_perf["fp"] * total_perf["fn"])
    mcc_denom <- sqrt( (total_perf["tp"] + total_perf["fp"]) * (total_perf["tp"] + total_perf["fn"]) * (total_perf["tn"] + total_perf["fp"]) * (total_perf["tn"] + total_perf["fn"]))
    
    total_perf["mcc"] <- mcc_num/mcc_denom
    total_perf["f1"] <- (2 * total_perf["tp"])/(2 * total_perf["tp"]+ total_perf["fp"] + total_perf["fn"])
    
    return(total_perf) # return totals
  }
}

assess_support_recovery <- function(varcount_mean, varcount_prob, true_support, R){
  
  # Method 1 (varcount_mean): r is in active set j iff E[# times Z_j split in ensemble j | ...] > critical threshold
  # Method 2 (varcount_probs): r is in active set j iff P(Z_j split on at least once in ensemble j | ...) > critical threshold
  
  metrics <- c("tp", "fp", "fn", "tn", "sen", "spec", "prec", "acc", "mcc", "f1")
  
  recovery_varcount_mean <- matrix(nrow = 101, ncol = length(metrics), dimnames = list(c(), metrics))
  recovery_varcount_prob <- matrix(nrow = 101, ncol = length(metrics), dimnames = list(c(), metrics))
  
  varcount_mean_seq <- seq(0, R, length = 101)
  varcount_prob_seq <- seq(0, 1, by = 0.01)
  
  
  for(i in 1:101){
    varcount_mean_support <- list()
    varcount_prob_support <- list()
    
    for(j in 1:6){
      varcount_mean_support[[j]] <- which(varcount_mean[,j] > varcount_mean_seq[i])
      varcount_prob_support[[j]] <- which(varcount_prob[,j] > varcount_prob_seq[i])
    }
    recovery_varcount_mean[i,] <- get_confusion_matrix(true_support, varcount_mean_support, R)
    recovery_varcount_prob[i,] <- get_confusion_matrix(true_support, varcount_prob_support, R)
  }
  
  return(list("varcount_mean" = recovery_varcount_mean,
              "varcount_prob" = recovery_varcount_prob))
  
  
}