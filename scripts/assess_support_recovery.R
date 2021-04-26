


assess_support_recovery <- function(true_support, est_support, R){
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
    return(colSums(support_perf)) # return totals
  }
  
  
}