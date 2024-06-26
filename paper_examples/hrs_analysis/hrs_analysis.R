#args <- commandArgs(TRUE)
#chain_num <- as.numeric(args[1])

# In paper, we ran 4 chains in parallel in our HTC environment
# The chain number was passed as a command line argument
# To run locally, just increment the chain number manually or run in a loop
chain_num <- 1
set.seed(424+chain_num)

load("hrs_data1.RData")

p <- ncol(X1_all)
R_cont <- ncol(Z_cont_all)

R_cat <- ncol(Z_cat_all)

R <- R_cont + R_cat

M <- 50
mu0 <- rep(0, times = p+1)
tau <- rep(0.5/sqrt(M), times = p+1)

subj_id_all <- c()
for(i in 1:n_all){
  subj_id_all <- c(subj_id_all, rep(i, times = ni_all[i]))
}

nd <- 1000
burn <- 1000
thin <- 1

fit <- VCBART::VCBART_cs(Y_train = Y_all,
                         subj_id_train = subj_id_all,
                         ni_train = ni_all,
                         X_train = X1_all,
                         Z_cont_train = Z_cont_all,
                         Z_cat_train = Z_cat_all,
                         unif_cuts = unif_cuts,
                         cutpoints_list = cutpoints_list,
                         cat_levels_list = cat_levels_list,
                         sparse = TRUE,
                         M = 50, mu0 = mu0, tau = tau,
                         nd = nd, burn = burn, thin = thin,
                         save_trees = TRUE, save_samples = FALSE,
                         verbose = FALSE)

assign(paste0("chain", chain_num),
       list(trees = fit$trees,
            varcounts = fit$varcounts[-(1:burn),,],
            sigma = fit$sigma))

save(list = paste0("chain", chain_num), file = paste0("hrs_chain", chain_num, ".RData"))
