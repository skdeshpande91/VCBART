source("scripts/vcbart_wrapper.R")
load("data/p5R20_data.RData")

#sim_number <- 1
#set.seed(129 + sim_number)

#tmp_index <- sample(1:n, replace = FALSE)
#train_id <- sort(tmp_index[1:floor(n*0.75)])
#test_id <- sort(tmp_index[(1 + floor(n*0.75)):n])

#train_index <- c()
#test_index <- c()

#for(id in train_id){
#  train_index <- c(train_index, start_index_all[id]:end_index_all[id])
#}
#
#for(id in test_id){
##  test_index <- c(test_index, start_index_all[id]:end_index_all[id])
#}

#X_train <- X_all[train_index,]
#Z_train <- Z_all[train_index,]
#Y_train <- Y_all[train_index]
#n_train <- n_all[train_id]
#beta_train <- beta_all[train_index,]


#X_test <- X_all[test_index,]
#Z_test <- Z_all[test_index,]
#Y_test <- Y_all[test_index]
#n_test <- n_all[test_id]
#beta_test <- beta_all[test_index,]



vcbart_adapt_subset <- 
  vcbart_wrapper(Y_all, X_all, Z_all, n_all, 
                 X_all, Z_all, n_all, cutpoints, error_structure = "ind",
                 split_probs_type = "adaptive", burn = 250, nd = 1000,
                 verbose = TRUE, print_every = 250)



# Make a plot with the subsetted data!

#Z_plot <- rbind(Z_train, Z_test)
#beta_plot <- rbind(beta_train, beta_test)
#beta0_hat <- rbind(vcbart_adapt_subset$train$beta[,,1], vcbart_adapt_subset$test$beta[,,1])
#beta1_hat <- rbind(vcbart_adapt_subset$train$beta[,,2], vcbart_adapt_subset$test$beta[,,2])
#beta2_hat <- rbind(vcbart_adapt_subset$train$beta[,,3], vcbart_adapt_subset$test$beta[,,3])
#beta3_hat <- rbind(vcbart_adapt_subset$train$beta[,,4], vcbart_adapt_subset$test$beta[,,4])

Z_plot <- Z_all
beta_plot <- beta_all
beta0_hat <- vcbart_adapt_subset$train$beta[,,1]
beta1_hat <- vcbart_adapt_subset$train$beta[,,2]
beta2_hat <- vcbart_adapt_subset$train$beta[,,3]
beta3_hat <- vcbart_adapt_subset$train$beta[,,4]

save(Z_plot, beta_plot, beta0_hat, beta1_hat, beta2_hat, beta3_hat, file = "results/illustration_p5R20.RData")


###########
# Plot beta0
############
ix0 <- which(Z_plot[,2] < 0.5)
ix1 <- which(Z_plot[,2] >= 0.5)
plot0_ix0 <- ix0[order(Z_plot[ix0,1])]
plot0_ix1 <- ix1[order(Z_plot[ix1,1])]
ylim0 <- max(abs(c(beta_plot[,1], beta0_plot[,c("L95", "U95")])))

plot(1,type= "n", xlim = c(0,1), ylim = c(-1.05, 1.05) * ylim0,
     main = expression("Intercept"), ylab = expression(beta[0]), xlab = expression(Z[1]))
lines(Z_plot[plot0_ix0,1], beta_plot[plot0_ix0,1])
lines(Z_plot[plot0_ix1,1], beta_plot[plot0_ix1,1], col = 'black', lty = 2)

polygon(c(Z_plot[plot0_ix0,1], rev(Z_plot[plot0_ix0,1])),
        c(beta0_hat[plot0_ix0,"L95"], rev(beta0_hat[plot0_ix0, "U95"])),
        col = rgb(0,0,1,1/3), border = NA)
polygon(c(Z_plot[plot0_ix1,1], rev(Z_plot[plot0_ix1,1])),
        c(beta0_hat[plot0_ix1,"L95"], rev(beta0_hat[plot0_ix1, "U95"])),
        col = rgb(0,0,1,1/3), border = NA)


legend("bottomright", cex = 0.8, legend = expression(Z[2] == 0), lty = 2, bty = "n")
legend("topleft", cex = 0.8, legend = expression(Z[2] == 1), lty = 1, bty = "n")
lines(Z_plot[plot0_ix0,1], beta0_hat[plot0_ix0,"MEAN"], col = 'blue')
lines(Z_plot[plot0_ix1,1], beta0_hat[plot0_ix1,"MEAN"], col = 'blue', lty = 2)

#############
# Plot beta1
#############
plot1_ix <- order(Z_plot[,1])
ylim1 <- max(abs(c(beta_plot[,2], beta1_hat[,c("L95", "U95")])))

plot(1, type = "n", xlim = c(0,1), ylim = c(-1.05, 1.05)*ylim1,
     main = expression("Effect of"~X[1]), xlab = expression(Z[1]), ylab = expression(beta[1]))
polygon(c(Z_plot[plot1_ix,1], rev(Z_plot[plot1_ix,1])),
        c(beta1_hat[plot1_ix,"L95"], rev(beta1_hat[plot1_ix, "U95"])),
        col = rgb(0,0,1,1/3), border = NA)
lines(Z_plot[plot1_ix,1], beta1_hat[plot1_ix,"MEAN"], col = 'blue')
lines(Z_plot[plot1_ix,1], beta_plot[plot1_ix,2])



##################
vcbart_adapt_all <- 
       vcbart_wrapper(Y_all, X_all, Z_all, n_all, 
                      X_all, Z_all, n_all, cutpoints, error_structure = "ind",
                      split_probs_type = "adaptive", burn = 250, nd = 1000,
                      verbose = TRUE, print_every = 250)



plot(Z_all[,1], vcbart_adapt$train$beta[,"MEAN",1])

plot(Z_all[,1], vcbart_adapt$train$beta[,"MEAN",2])

plot(Z_all[,1], vcbart_adapt$train$beta[,"MEAN",3])

plot(Z_all[,1], vcbart_adapt$train$beta[,"MEAN",4])

tmp_Z <- rbind(Z_train, Z_test)
tmp_beta <- rbind(beta_train, beta_test)




plot_ix0 <- order(tmp_Z[tmp_Z[,2] > 0.5,1])
plot_ix1 <- order(tmp_Z[tmp_Z[,2] < 0.5, 1])
plot_ix <- order(tmp_Z[,1])


plot(tmp_Z[plot_ix,1],tmp_beta[plot_ix,2], type = "l")
points(tmp_Z[plot_ix,1], beta_hat[plot_ix,2], pch = 16, cex = 0.5, col = 'blue')

beta_hat <- rbind(vcbart_adapt$train$beta[,"MEAN",], vcbart_adapt$test$beta[,"MEAN",])




