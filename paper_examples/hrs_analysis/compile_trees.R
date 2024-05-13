load("../hrs_analysis/hrs_data1.RData")

nd <- 1000
n_chains <- 4
R_cont <- ncol(Z_cont_all)
R_cat <- ncol(Z_cat_all)
R <- R_cont + R_cat

Z_names <- c(colnames(Z_cont_all), colnames(Z_cat_all))
p <- ncol(X1_all)


varcounts <- array(dim = c(n_chains * nd, R, p+1),
                   dimnames = list(c(), Z_names, 
                                   c("Intercept",colnames(X1_all))))
all_trees <- list()

for(chain_num in 1:n_chains){
  file_name <- paste0("hrs_chain", chain_num, ".RData")
  if(file.exists(file_name)){
    load(file_name)
    start_index <- (chain_num-1)*nd + 1
    end_index <- chain_num*nd
    varcounts[start_index:end_index,,] <- get(paste0("chain", chain_num))$varcounts
    
    all_trees <- c(all_trees, get(paste0("chain", chain_num))$trees)
    
    rm(list = paste0("chain", chain_num))
  }
}

selection_prob <- apply(varcounts >= 1, FUN = mean, MARGIN = c(2,3))
save(all_trees, file = "hrs_all_trees.RData")
save(selection_prob, file = "hrs_selection_probs.RData")
