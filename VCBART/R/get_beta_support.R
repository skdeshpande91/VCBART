# Function to select modifiers

get_beta_support <- function(chain1, chain2, burn, max_cutoff = 10){
  if(("var_counts_samples" %in% names(chain1)) && ("var_counts_samples" %in% names(chain2))){

    R <- dim(chain1$var_counts_samples)[1]
    p <- dim(chain1$var_counts_samples)[2]
    exceed_cutoff_probs <- array(dim = c(R,p,max_cutoff))
    support <- list()
    # exceed_cutoff_probs[r,k,c] reports posterior probability
    # modifier r is used in ensemble k at least c times
    for(cutoff in 1:max_cutoff){
      tmp_support <- list()
      exceed_cutoff_probs[,,cutoff] <- 
        0.5 * apply(chain1$var_counts_samples[,,-(1:burn)] >= cutoff, MARGIN = c(1,2), FUN = mean) +
        0.5 * apply(chain2$var_counts_samples[,,-(1:burn)] >= cutoff, MARGIN = c(1,2), FUN = mean)
      
      for(k in 1:p) tmp_support[[k]] <- which(exceed_cutoff_probs[,k,cutoff] >= 0.5)
      support[[cutoff]] <- tmp_support
    }
    
    return(list(exceed_cutoff_probs = exceed_cutoff_probs, support = support))
 
    # now look for the entries > 0.5 and write the results back as a list
    
  } else{
    stop("chain1 and chain2 must have named an element named var_counts_samples")
  }
}