source("true_betas_p5R20.R")

n_train <- 250
n_test <- 25
sigma <- 1
sim_number <- 1
set.seed(417)
source("figure2_generate_data.R")

my_colors <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


######
source("vcbart_wrapper.R")

M <- 50
tau <- rep(0.5/sqrt(M), times = p+1)

fit <- vcbart_wrapper(Y_train = Y_train,
                      subj_id_train = subj_id_train,
                      ni_train = ni_train,
                      X_train = X_train,
                      Z_cont_train = Z_cont_train,
                      X_test = X_test,
                      Z_cont_test = Z_cont_test,
                      unif_cuts = unif_cuts,
                      sparse = TRUE,
                      M = M,tau = tau,
                      verbose = TRUE,
                      n_chains = 4)



###############

Z_plot <- (1+Z_cont_all)/2
beta_plot <- beta_all

beta0_hat <- rbind(fit$train$beta[,,1], fit$test$beta[,,1])
beta1_hat <- rbind(fit$train$beta[,,2], fit$test$beta[,,2])
beta2_hat <- rbind(fit$train$beta[,,3], fit$test$beta[,,3])
beta3_hat <- rbind(fit$train$beta[,,4], fit$test$beta[,,4])
beta4_hat <- rbind(fit$train$beta[,,5], fit$test$beta[,,5])
beta5_hat <- rbind(fit$train$beta[,,6], fit$test$beta[,,6])

save(beta0_hat, beta1_hat, beta2_hat, beta3_hat, beta4_hat, beta5_hat,
     beta0_all, beta1_all, beta2_all, beta3_all, beta4_all, beta5_all,
     file = "beta_hats.RData")


ix0 <- which(Z_plot[,2] < 0.5)
ix1 <- which(Z_plot[,2] >= 0.5)
plot0_ix0 <- ix0[order(Z_plot[ix0,1])]
plot0_ix1 <- ix1[order(Z_plot[ix1,1])]
ylim0 <- max(abs(c(beta_plot[,1], beta0_hat[,c("L95", "U95")])))


pdf("p5R20_beta.pdf", width = 6, height = 3)
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0), mfrow = c(2,4))

plot(1,type= "n", xlim = c(0,1), ylim = c(-1.05, 1.05) * ylim0,
     main = expression("Intercept"), ylab = expression(beta[0]), xlab = expression(Z[1]))
lines(Z_plot[plot0_ix0,1], beta_plot[plot0_ix0,1], lty = 1, lwd = 0.75)
lines(Z_plot[plot0_ix1,1], beta_plot[plot0_ix1,1], col = 'black', lty = 2, lwd = 0.75)
legend("topleft", cex = 0.8, legend = expression(Z[2] >= 0.5), lty = 1, bty = "n")
legend("bottomright", cex = 0.8, legend = expression(Z[2] < 0.5), lty = 2, bty = "n")

############
# Plot 6: beta1 plus VCBART results
plot1_ix <- order(Z_plot[,1])
ylim1 <- max(abs(c(beta_plot[,2], beta1_hat[,c("L95", "U95")])))

plot(1, type = "n", xlim = c(0,1), ylim = c(-1.05, 1.05)*ylim1,
     main = expression("Effect of"~X[1]), xlab = expression(Z[1]), ylab = expression(beta[1]))
lines(Z_plot[plot1_ix,1], beta_plot[plot1_ix,2], lwd = 0.75)

# Plot 7: beta2 & VCBART results
ix0 <- which(Z_plot[,1] > 0.6)
ix1 <- which(Z_plot[,1] <= 0.25)
ix2 <- which(Z_plot[,1] > 0.25 & Z_plot[,1] < 0.6)

plot2_ix0 <- ix0[order(Z_plot[ix0,1])]
plot2_ix1 <- ix1[order(Z_plot[ix1,1])]
plot2_ix2 <- ix2[order(Z_plot[ix2,1])]

ylim2 <- max(abs(c(beta_plot[,3], beta2_hat[,c("L95", "U95")])))

plot(1, type = "n", xlim = c(0,1), ylim = c(-1.05, 1.05)*ylim2,
     main = expression("Effect of"~X[2]), xlab = expression(Z[1]), ylab = expression(beta[2]))

lines(Z_plot[plot2_ix0,1], beta_plot[plot2_ix0, 3], lwd = 0.75)
lines(Z_plot[plot2_ix1,1], beta_plot[plot2_ix1, 3], lwd = 0.75)
lines(Z_plot[plot2_ix2,1], beta_plot[plot2_ix2, 3], lwd = 0.75)

# Plot 8: beta3 & VCBART results
plot3_ix <- order(Z_plot[,1])
ylim3 <- max(abs(c(beta_plot[,4], beta3_hat[,c("L95", "U95")])))

plot(1, type = "n", xlim = c(0,1), ylim = c(-1.05, 1.05)*ylim3,
     main = expression("Effect of"~X[3]), xlab = expression(Z[1]), ylab = expression(beta[3]))
lines(Z_plot[plot3_ix,1], beta_plot[plot3_ix,4], lwd = 0.75)

plot(1,type= "n", xlim = c(0,1), ylim = c(-1.05, 1.05) * ylim0,
     main = expression("Intercept"), ylab = expression(beta[0]), xlab = expression(Z[1]))
lines(Z_plot[plot0_ix0,1], beta_plot[plot0_ix0,1], lty = 1, lwd = 0.75)
lines(Z_plot[plot0_ix1,1], beta_plot[plot0_ix1,1], col = 'black', lty = 2, lwd = 0.75)
legend("topleft", cex = 0.8, legend = expression(Z[2] >= 0.5), lty = 1, bty = "n")
legend("bottomright", cex = 0.8, legend = expression(Z[2] < 0.5), lty = 2, bty = "n")
polygon(c(Z_plot[plot0_ix0,1], rev(Z_plot[plot0_ix0,1])),
        c(beta0_hat[plot0_ix0,"L95"], rev(beta0_hat[plot0_ix0, "U95"])),
        col = rgb(0.5,0.5,0.5,1/3), border = NA)
polygon(c(Z_plot[plot0_ix1,1], rev(Z_plot[plot0_ix1,1])),
        c(beta0_hat[plot0_ix1,"L95"], rev(beta0_hat[plot0_ix1, "U95"])),
        col = rgb(0.5,0.5,0.5,1/3), border = NA)
lines(Z_plot[plot0_ix0,1], beta0_hat[plot0_ix0,"MEAN"], lwd = 1)
lines(Z_plot[plot0_ix1,1], beta0_hat[plot0_ix1,"MEAN"], lwd = 1, lty = 2)

############
# Plot 6: beta1 plus VCBART results
plot1_ix <- order(Z_plot[,1])
ylim1 <- max(abs(c(beta_plot[,2], beta1_hat[,c("L95", "U95")])))

plot(1, type = "n", xlim = c(0,1), ylim = c(-1.05, 1.05)*ylim1,
     main = expression("Effect of"~X[1]), xlab = expression(Z[1]), ylab = expression(beta[1]))
polygon(c(Z_plot[plot1_ix,1], rev(Z_plot[plot1_ix,1])),
        c(beta1_hat[plot1_ix,"L95"], rev(beta1_hat[plot1_ix, "U95"])),
        col = rgb(0.5,0.5,0.5,1/3), border = NA)
lines(Z_plot[plot1_ix,1], beta1_hat[plot1_ix,"MEAN"], lwd = 1)
lines(Z_plot[plot1_ix,1], beta_plot[plot1_ix,2], lwd = 0.75)

# Plot 7: beta2 & VCBART results
ix0 <- which(Z_plot[,1] > 0.6)
ix1 <- which(Z_plot[,1] <= 0.25)
ix2 <- which(Z_plot[,1] > 0.25 & Z_plot[,1] < 0.6)

plot2_ix0 <- ix0[order(Z_plot[ix0,1])]
plot2_ix1 <- ix1[order(Z_plot[ix1,1])]
plot2_ix2 <- ix2[order(Z_plot[ix2,1])]

ylim2 <- max(abs(c(beta_plot[,3], beta2_hat[,c("L95", "U95")])))

plot(1, type = "n", xlim = c(0,1), ylim = c(-1.05, 1.05)*ylim2,
     main = expression("Effect of"~X[2]), xlab = expression(Z[1]), ylab = expression(beta[2]))

polygon(c(Z_plot[plot2_ix0,1], rev(Z_plot[plot2_ix0,1])),
        c(beta2_hat[plot2_ix0,"L95"], rev(beta2_hat[plot2_ix0, "U95"])),
        col = rgb(0.5,0.5,0.5,1/3), border = NA)
polygon(c(Z_plot[plot2_ix1,1], rev(Z_plot[plot2_ix1,1])),
        c(beta2_hat[plot2_ix1,"L95"], rev(beta2_hat[plot2_ix1, "U95"])),
        col = rgb(0.5,0.5,0.5,1/3), border = NA)
polygon(c(Z_plot[plot2_ix2,1], rev(Z_plot[plot2_ix2,1])),
        c(beta2_hat[plot2_ix2,"L95"], rev(beta2_hat[plot2_ix2, "U95"])),
        col = rgb(0.5,0.5,0.5,1/3), border = NA)
lines(Z_plot[plot2_ix0,1], beta2_hat[plot2_ix0, "MEAN"], lwd = 1)
lines(Z_plot[plot2_ix1,1], beta2_hat[plot2_ix1, "MEAN"], lwd = 1)
lines(Z_plot[plot2_ix2,1], beta2_hat[plot2_ix2, "MEAN"], lwd = 1)

lines(Z_plot[plot2_ix0,1], beta_plot[plot2_ix0, 3], lwd = 0.75)
lines(Z_plot[plot2_ix1,1], beta_plot[plot2_ix1, 3], lwd = 0.75)
lines(Z_plot[plot2_ix2,1], beta_plot[plot2_ix2, 3], lwd = 0.75)

# Plot 8: beta3 & VCBART results
plot3_ix <- order(Z_plot[,1])
ylim3 <- max(abs(c(beta_plot[,4], beta3_hat[,c("L95", "U95")])))

plot(1, type = "n", xlim = c(0,1), ylim = c(-1.05, 1.05)*ylim3,
     main = expression("Effect of"~X[3]), xlab = expression(Z[1]), ylab = expression(beta[3]))
polygon(c(Z_plot[plot3_ix,1], rev(Z_plot[plot3_ix,1])),
        c(beta3_hat[plot3_ix,"L95"], rev(beta3_hat[plot3_ix, "U95"])),
        col = rgb(0.5,0.5,0.5,1/3), border = NA)
lines(Z_plot[plot3_ix,1], beta3_hat[plot3_ix,"MEAN"], lwd = 1)
lines(Z_plot[plot3_ix,1], beta_plot[plot3_ix,4], lwd = 0.75)

dev.off()

