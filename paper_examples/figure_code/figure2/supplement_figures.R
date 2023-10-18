#################
# Figure for the Supplementary materials
# True betas plotted against estimated betas

my_colors <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


load("beta_hats.RData")

supp_ylim0 <- range(c(beta0_all, beta0_hat[,1]))
supp_ylim1 <- range(c(beta1_all, beta1_hat[,1]))
supp_ylim2 <- range(c(beta2_all, beta2_hat[,1]))
supp_ylim3 <- range(c(beta3_all, beta3_hat[,1]))
supp_ylim4 <- range(c(beta4_all, beta4_hat[,1]))
supp_ylim5 <- range(c(beta5_all, beta5_hat[,1]))

#########
pdf("../../figures/p5R20_demo_beta0.pdf", width = 4.5, height = 4.5)
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0))
plot(1, type = "n",
     xlim = supp_ylim0 , ylim = supp_ylim0 ,
     main = expression("Intercept"), 
     xlab = expression("Actual"),
     ylab = expression("Estimated"))
points(beta0_all, beta0_hat[,1], pch = 16, cex = 0.75)
abline(a = 0, b = 1, col = my_colors[7], lty = 2)
dev.off()

pdf("../../figures/p5R20_demo_beta1.pdf", width = 4.5, height = 4.5)
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0))
plot(1, type = "n",
     xlim = supp_ylim1, ylim = supp_ylim1,
     main = expression(beta[1]), 
     xlab = expression("Actual"),
     ylab = expression("Estimated"))
points(beta1_all, beta1_hat[,1], pch = 16, cex = 0.75)
abline(a = 0, b = 1, col = my_colors[7], lty = 2)
dev.off()

pdf("../../figures/p5R20_demo_beta2.pdf", width = 4.5, height = 4.5)
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0))
plot(1, type = "n",
     xlim = supp_ylim2, ylim = supp_ylim2,
     main = expression(beta[2]), 
     xlab = expression("Actual"),
     ylab = expression("Estimated"))
points(beta2_all, beta2_hat[,1], pch = 16, cex = 0.75)
abline(a = 0, b = 1, col = my_colors[7], lty = 2)
dev.off()

pdf("../../figures/p5R20_demo_beta3.pdf", width = 4.5, height = 4.5)
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0))
plot(1, type = "n",
     xlim = supp_ylim3, ylim = supp_ylim3,
     main = expression(beta[3]), 
     xlab = expression("Actual"),
     ylab = expression("Estimated"))
points(beta3_all, beta3_hat[,1], pch = 16, cex = 0.75)
abline(a = 0, b = 1, col = my_colors[7], lty = 2)
dev.off()

pdf("../../figures/p5R20_demo_beta4.pdf", width = 4.5, height = 4.5)
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0))
plot(1, type = "n",
     xlim = supp_ylim4, ylim = supp_ylim4,
     main = expression(beta[4]), 
     xlab = expression("Actual"),
     ylab = expression("Estimated"))
points(beta4_all, beta4_hat[,1], pch = 16, cex = 0.75)
abline(a = 0, b = 1, col = my_colors[7], lty = 2)
dev.off()

pdf("../../figures/p5R20_demo_beta5.pdf", width = 4.5, height = 4.5)
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0))
plot(1, type = "n",
     xlim = supp_ylim5, ylim = supp_ylim5,
     main = expression(beta[5]), 
     xlab = expression("Actual"),
     ylab = expression("Estimated"))
points(beta5_all, beta5_hat[,1], pch = 16, cex = 0.75)
abline(a = 0, b = 1, col = my_colors[7], lty = 2)
dev.off()

##########

par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0), mfrow = c(2,3))

plot(1, type = "n", 
     xlim = supp_ylim0 * c(-1,1), ylim = supp_ylim0 * c(-1,1),
     main = expression("Intercept"), xlab = "Actual", ylab = "Estimated")
points(beta_all[train_index,1], fit$train$beta[,"MEAN",1],
       pch = 16, cex = 0.75, col = my_colors[3])
points(beta_all[test_index,1], fit$test$beta[,"MEAN",1],
       pch = 15, cex = 0.75, col = my_colors[2])
abline(a = 0, b = 1)
legend("bottomright", legend = c("Train", "Test"),
       pch = c(16, 15), col = my_colors[c(3,2)], cex = 0.75,
       bty = "n")

# beta1
plot(1, type = "n", 
     xlim = supp_ylim1 * c(-1,1), ylim = supp_ylim1 * c(-1,1),
     main = expression("Effect of"~X[1]), 
     xlab = "Actual", ylab = "Estimated")
points(beta_all[train_index,2], fit$train$beta[,"MEAN",2],
       pch = 16, cex = 0.75, col = my_colors[3])
points(beta_all[test_index,2], fit$test$beta[,"MEAN",2],
       pch = 15, cex = 0.75, col = my_colors[2])
abline(a = 0, b = 1)
legend("bottomleft", legend = c("Train", "Test"),
       pch = c(16, 15), col = my_colors[c(3,2)], cex = 0.75,
       bty = "n")

# beta2
plot(1, type = "n", 
     xlim = supp_ylim2 * c(-1,1), ylim = supp_ylim2 * c(-1,1),
     main = expression("Effect of"~X[2]), 
     xlab = "Actual", ylab = "Estimated")
points(beta_all[train_index,3], fit$train$beta[,"MEAN",3],
       pch = 16, cex = 0.75, col = my_colors[3])
points(beta_all[test_index,3], fit$test$beta[,"MEAN",3],
       pch = 15, cex = 0.75, col = my_colors[2])
abline(a = 0, b = 1)
legend("bottomright", legend = c("Train", "Test"),
       pch = c(16, 15), col = my_colors[c(3,2)], cex = 0.75,
       bty = "n")

# beta3
plot(1, type = "n", 
     xlim = supp_ylim3 * c(-1,1), ylim = supp_ylim3 * c(-1,1),
     main = expression("Effect of"~X[2]), 
     xlab = "Actual", ylab = "Estimated")
points(beta_all[train_index,4], fit$train$beta[,"MEAN",4],
       pch = 16, cex = 0.75, col = my_colors[3])
points(beta_all[test_index,4], fit$test$beta[,"MEAN",4],
       pch = 15, cex = 0.75, col = my_colors[2])
abline(a = 0, b = 1)
legend("bottomright", legend = c("Train", "Test"),
       pch = c(16, 15), col = my_colors[c(3,2)], cex = 0.75,
       bty = "n")






