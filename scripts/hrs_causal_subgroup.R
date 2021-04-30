library(VCBARTdev)

load("data/HRS_causal_vcbart_data.RData")


age_start <- c(60, 65, 70, 75, 80) * 12
age_end <- c(65 * 12 - 1, 70 * 12 - 1, 75 * 12 - 1, 80 * 12 - 1, 85 * 12)
age_index <- list()
for(aix in 1:5) age_index[[aix]] <- which(age_start[aix] <= Z_all[,"AGE"] & Z_all[,"AGE"] < age_end[aix])

educ_index <- list()
educ_index[[1]] <- which(Z_all[,"EDUC"] < 12) # LT HS
educ_index[[2]] <- which(Z_all[,"EDUC"] == 12) # HS degree but no college
educ_index[[3]] <- which(Z_all[,"EDUC"] > 12 & Z_all[,"EDUC"] < 16) # some college
educ_index[[4]] <- which(Z_all[,"EDUC"] >= 16) # college degree


X_test <- X_all[1:10,]
Z_test <- Z_all[1:10,]
n_test <- 10

chain1 <- VCBART(Y_all, X_all, Z_all, n_all,
                 X_test, Z_test, n_test, cutpoints,
                 intercept = TRUE, split_probs_type = "adaptive",
                 error_structure = "ind", burn = 1000, nd = 1000)
chain2 <- VCBART(Y_all, X_all, Z_all, n_all,
                 X_test, Z_test, n_test, cutpoints,
                 intercept = TRUE, split_probs_type = "adaptive",
                 error_structure = "ind", burn = 1000, nd = 1000)

beta1 <- chain1$beta_train_samples
beta2 <- chain2$beta_train_samples

save(beta1, file = "results/hrs_causal_beta_chain1.RData")
save(beta2, file = "results/hrs_causal_beta_chain2.RData")


age_samples <- list()
for(aix in 1:5){
  # beta_train_samples is n x (p+1) x nd
  # we want to average each covariate effect over the age groupings not over the samples
  # apply(, MARGIN = c(2,3), FUN = mean) returns a (p+1) x nd matrix containing draws of each effect, averaged over the individuals in the subgroup
  age_samples[[aix]] <- cbind(apply(beta1[age_index[[aix]],,], MARGIN = c(2,3), FUN = mean), 
                              apply(beta2[age_index[[aix]],,], MARGIN = c(2,3), FUN = mean))
  
}
educ_samples <- list()
for(eix in 1:4){
  educ_samples[[eix]] <- cbind(apply(beta1[educ_index[[eix]],,], MARGIN = c(2,3), FUN = mean),
                               apply(beta2[educ_index[[eix]],,], MARGIN = c(2,3), FUN = mean))
}

ate_samples <- cbind(apply(beta1, MARGIN = c(2,3), FUN = mean),
                     apply(beta2, MARGIN = c(2,3), FUN = mean))

save(age_samples, educ_samples, ate_samples, file = "results/hrs_causal_subgroup_samples.RData")


