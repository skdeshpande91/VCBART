###### 
# Code to produce Figure S4 in the Supplementary Materials
######

load("results/results_p5R20_sigma1.RData")
load("results/results_p5R20_sigma2.RData")
load("results/results_p5R20_sigma3.RData")
load("results/results_p5R20_sigma4.RData")

tmp_beta_mse1 <- apply(beta_mse_test_sigma1, MARGIN = c(1,3), FUN = mean, na.rm = TRUE)
tmp_beta_mse2 <- apply(beta_mse_test_sigma2, MARGIN = c(1,3), FUN = mean, na.rm = TRUE)
tmp_beta_mse3 <- apply(beta_mse_test_sigma3, MARGIN = c(1,3), FUN = mean, na.rm = TRUE)
tmp_beta_mse4 <- apply(beta_mse_test_sigma4, MARGIN = c(1,3), FUN = mean, na.rm = TRUE)

beta_mse_mean <- matrix(nrow = 4, ncol = dim(beta_mse_test_sigma1)[1], dimnames = list(c(), dimnames(beta_mse_test_sigma1)[[1]]))
beta_mse_mean[1,] <- apply(tmp_beta_mse1, MARGIN = 1, FUN = mean, na.rm = TRUE)
beta_mse_mean[2,] <- apply(tmp_beta_mse2, MARGIN = 1, FUN = mean, na.rm = TRUE)
beta_mse_mean[3,] <- apply(tmp_beta_mse3, MARGIN = 1, FUN = mean, na.rm = TRUE)
beta_mse_mean[4,] <- apply(tmp_beta_mse4, MARGIN = 1, FUN = mean, na.rm = TRUE)

beta_mse_sd <- matrix(nrow = 4, ncol = dim(beta_mse_test_sigma1)[1], dimnames = list(c(), dimnames(beta_mse_test_sigma1)[[1]]))
beta_mse_sd[1,] <- apply(tmp_beta_mse1, MARGIN = 1, FUN = sd, na.rm = TRUE)
beta_mse_sd[2,] <- apply(tmp_beta_mse2, MARGIN = 1, FUN = sd, na.rm = TRUE)
beta_mse_sd[3,] <- apply(tmp_beta_mse3, MARGIN = 1, FUN = sd, na.rm = TRUE)
beta_mse_sd[4,] <- apply(tmp_beta_mse4, MARGIN = 1, FUN = sd, na.rm = TRUE)




sigma_list <- c(0.5,1, 2, 4)

pch_list <- c(16, 17, 18, 15, 8)
col_list <- c("black", "green", "red", "blue", "orange")
method_list <- c("vcbart_adapt", "boosted_tvcm", "tvc", "lm", "kernel_smoothing")

png("figures/beta_mse_sigma.png", width = 4.5, height = 4.5, units = "in", res = 300)
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0))
plot(1, type = "n", xlim = log(c(0.4, 4.1)), ylim = c(0, 20), xlab = expression(sigma), ylab = "MSE", xaxt = "n",
     main = "Covariate effect recovery")
axis(side = 1, at = log(sigma_list), labels = sigma_list)

for(i in 1:5){
  lines(log(sigma_list), beta_mse_mean[,method_list[i]], col = col_list[i])
  
  for(j in 1:4){
    lines(x = log(c(sigma_list[j], sigma_list[j])), y = beta_mse_mean[j, method_list[i]] + beta_mse_sd[j, method_list[i]] * c(-1,1), lty = 2, col = col_list[i])
  }
  points(log(sigma_list), beta_mse_mean[,method_list[i]], col = col_list[i], pch = 16, cex = 1)
  
  
}
#legend("topleft", legend = c("VCBART", "BTVCM", "TVCM", "lm", "KS"), pch = pch_list, col = col_list)
legend("topleft", legend = c("VCBART", "BTVCM", "TVCM", "lm", "KS"), pch = 16, col = col_list)

dev.off()


ystar_smse_mean <- matrix(nrow = 4, ncol = dim(ystar_smse_test_sigma1)[1], dimnames = list(c(), dimnames(ystar_smse_test_sigma1)[[1]]))
ystar_smse_mean[1,] <- apply(ystar_smse_test_sigma1, MARGIN = 1, FUN = mean, na.rm = TRUE)
ystar_smse_mean[2,] <- apply(ystar_smse_test_sigma2, MARGIN = 1, FUN = mean, na.rm = TRUE)
ystar_smse_mean[3,] <- apply(ystar_smse_test_sigma3, MARGIN = 1, FUN = mean, na.rm = TRUE)
ystar_smse_mean[4,] <- apply(ystar_smse_test_sigma4, MARGIN = 1, FUN = mean, na.rm = TRUE)

ystar_smse_sd <- ystar_smse_mean
ystar_smse_sd[1,] <- apply(ystar_smse_test_sigma1, MARGIN = 1, FUN = sd, na.rm = TRUE)
ystar_smse_sd[2,] <- apply(ystar_smse_test_sigma2, MARGIN = 1, FUN = sd, na.rm = TRUE)
ystar_smse_sd[3,] <- apply(ystar_smse_test_sigma3, MARGIN = 1, FUN = sd, na.rm = TRUE)
ystar_smse_sd[4,] <- apply(ystar_smse_test_sigma4, MARGIN = 1, FUN = sd, na.rm = TRUE)


pch_list2 <- c(16, 17, 18, 15, 8, 3, 4, 10)
method_list2 <- c("vcbart_adapt", "boosted_tvcm", "tvc", "lm", "kernel_smoothing", "bart", "extraTrees", "gbm")
col_list2 <- c("black", "green", "red", "blue", "orange", "cyan", "purple", "pink")
png("figures/ystar_smse_sigma.png", width = 4.5, height = 4.5, units = "in", res = 300)
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0))
plot(1, type = "n", xlim = log(c(0.4, 4.1)), ylim = c(0, 0.5), xlab = expression(sigma), ylab = "SMSE", xaxt = "n",
     main = "Predictive performance")
axis(side = 1, at = log(sigma_list), labels = sigma_list)

for(i in 1:8){
  lines(log(sigma_list), ystar_smse_mean[,method_list2[i]], col = col_list2[i])
  
  for(j in 1:4){
    lines(x = log(c(sigma_list[j], sigma_list[j])), y = ystar_smse_mean[j, method_list2[i]] + ystar_smse_sd[j, method_list2[i]] * c(-1,1), lty = 2, col = col_list2[i])
  }
  #points(log(sigma_list), ystar_smse_mean[,method_list2[i]], pch = pch_list2[i], cex = 0.8, col = col_list2[i])
  points(log(sigma_list), ystar_smse_mean[,method_list2[i]], pch = 16, cex = 0.8, col = col_list2[i])
  
  
}
legend("topleft", legend = c("VCBART", "BTVCM", "TVCM", "lm"), pch = 16, col = col_list2[1:4])
legend("topright", legend = c("KS", "BART", "ERT", "GBM"), pch = 16, col = col_list2[5:8])
dev.off()


png("~/Documents/Research/vc_bart/figures/sim_p5R20_sigma.png", width = 8, height = 4, units = "in", res = 300)
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0), mfrow = c(1,2))

pch_list <- c(16, 17, 18, 15, 8)
col_list <- c("black", "green", "red", "blue", "orange")
method_list <- c("vcbart_adapt", "boosted_tvcm", "tvc", "lm", "kernel_smoothing")

plot(1, type = "n", xlim = log(c(0.4, 4.1)), ylim = c(0, 20), xlab = expression(sigma), ylab = "MSE", xaxt = "n",
     main = "Covariate effect recovery")
axis(side = 1, at = log(sigma_list), labels = sigma_list)

for(i in 1:5){
  lines(log(sigma_list), beta_mse_mean[,method_list[i]], col = col_list[i])
  
  #for(j in 1:4){
  #  lines(x = log(c(sigma_list[j], sigma_list[j])), y = beta_mse_mean[j, method_list[i]] + beta_mse_sd[j, method_list[i]] * c(-1,1), lty = 2, col = col_list[i])
  #}
  points(log(sigma_list), beta_mse_mean[,method_list[i]], col = col_list[i], pch = 16, cex = 1)
  
  
}
legend("topleft", legend = c("VCBART", "BTVCM", "TVCM", "lm", "KS"), pch = 16, col = col_list)

pch_list2 <- c(16, 17, 18, 15, 8, 3, 4, 10)
method_list2 <- c("vcbart_adapt", "boosted_tvcm", "tvc", "lm", "kernel_smoothing", "bart", "extraTrees", "gbm")
col_list2 <- c("black", "green", "red", "blue", "orange", "cyan", "purple", "pink")
plot(1, type = "n", xlim = log(c(0.4, 4.1)), ylim = c(0, 0.5), xlab = expression(sigma), ylab = "SMSE", xaxt = "n",
     main = "Predictive performance")
axis(side = 1, at = log(sigma_list), labels = sigma_list)

for(i in 1:8){
  lines(log(sigma_list), ystar_smse_mean[,method_list2[i]], col = col_list2[i])
  
  #for(j in 1:4){
  #  lines(x = log(c(sigma_list[j], sigma_list[j])), y = ystar_smse_mean[j, method_list2[i]] + ystar_smse_sd[j, method_list2[i]] * c(-1,1), lty = 2, col = col_list2[i])
  #}
  #points(log(sigma_list), ystar_smse_mean[,method_list2[i]], pch = pch_list2[i], cex = 0.8, col = col_list2[i])
  points(log(sigma_list), ystar_smse_mean[,method_list2[i]], pch = 16, cex = 0.8, col = col_list2[i])
  
  
}
legend("topleft", legend = c("VCBART", "BTVCM", "TVCM", "lm"), pch = 16, col = col_list2[1:4])
legend("topright", legend = c("KS", "BART", "ERT", "GBM"), pch = 16, col = col_list2[5:8])
dev.off()