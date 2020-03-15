# Script to re-create figure 1 of the main text

library(VCBART)
load("data/p3R2_data.RData")

# run on the full dataset
# add a column of 1's to X_all for the intercept
X_all <- cbind(rep(1, times = N), X_all)

chain1 <- vc_BART_ind(Y = Y_all, X_train = X_all, Z_train = Z_all,
                      X_test = X_all, Z_test = Z_all, xinfo_list = cutpoints, 
                      nd = 1000, burn = 250, verbose = TRUE, print_every = 50)

chain2 <- vc_BART_ind(Y = Y_all, X_train = X_all, Z_train = Z_all,
                      X_test = X_all, Z_test = Z_all, xinfo_list = cutpoints, 
                      nd = 1000, burn = 250, verbose = TRUE, print_every = 50)

fit_sum <- get_summary(chain1, chain2)

beta0_hat <- fit_sum$train$beta[,,1]
beta1_hat <- fit_sum$train$beta[,,2]
beta2_hat <- fit_sum$train$beta[,,3]
beta3_hat <- fit_sum$train$beta[,,4]

ylim0 <- max(abs(c(beta0_all, beta0_hat)))
ylim1 <- max(abs(c(beta1_all, beta1_hat)))
ylim2 <- max(abs(c(beta2_all, beta2_hat)))
ylim3 <- max(abs(c(beta3_all, beta3_hat)))

ylim_all <- max(ylim0, ylim1, ylim2, ylim3)

png("figures/beta_p3R2.png", width = 6.5, height = 3.25, units = "in", res = 600)
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0), mfrow = c(2,4), cex.main = 0.95, cex.lab = 0.95, cex.axis = 0.9)

# The true function beta_0
plot(1, type = "n", xlim = c(0,1), ylim = c(-1.01, 1.01) * ylim_all, main = expression("True intercept"), ylab = expression(beta[0]), xlab = expression(Z[1]))

ix0 <- which(Z_all[,2] < 0.5)
ix1 <- which(Z_all[,2] >= 0.5)

plot0_ix0 <- ix0[order(Z_all[ix0,1])]
plot0_ix1 <- ix1[order(Z_all[ix1,1])]

lines(Z_all[plot0_ix0,1], beta0_all[plot0_ix0], col = 'black')
lines(Z_all[plot0_ix1,1], beta0_all[plot0_ix1], col = 'black', lty = 2)

legend("bottomright", cex = 0.8, legend = expression(Z[2] >= 0.5), lty = 2, bty = "n")
legend("topleft", cex = 0.8, legend = expression(Z[2] < 0.5), lty = 1, bty = "n")

# The true function beta_1
plot(1, type = "n", xlim = c(0,1), ylim = c(-1.01, 1.01) * ylim_all, main = expression("True effect of"~X[1]), 
     ylab = expression(beta[1]), xlab = expression(Z[1]))
plot1_ix <- order(Z_all[,1])
lines(Z_all[plot1_ix,1], beta1_all[plot1_ix], col = 'black')


# The true function beta_2
plot(1, type = "n", xlim = c(0,1), ylim = c(-1.01, 1.01) * ylim_all, main = expression("True effect of"~X[2]), 
     ylab = expression(beta[2]), xlab = expression(Z[1]))

ix0 <- which(Z_all[,1] > 0.6)
ix1 <- which(Z_all[,1] < 0.25)
ix2 <- which(Z_all[,1] >= 0.25 & Z_all[,1] <= 0.6)


plot2_ix0 <- ix0[order(Z_all[ix0,1])]
plot2_ix1 <- ix1[order(Z_all[ix1,1])]
plot2_ix2 <- ix2[order(Z_all[ix2,1])]


lines(Z_all[plot2_ix0,1], beta2_all[plot2_ix0], col = 'black')
lines(Z_all[plot2_ix1,1], beta2_all[plot2_ix1], col = 'black')
lines(Z_all[plot2_ix2,1], beta2_all[plot2_ix2], col = 'black')


# The true function beta_3
plot(1, type = "n", xlim = c(0,1), ylim = c(-1.01, 1.01) * ylim_all, main = expression("True effect of"~X[3]), 
     ylab = expression(beta[3]), xlab = expression(Z[1]))
plot1_ix <- order(Z_all[,1])

lines(Z_all[plot1_ix,1], beta3_all[plot1_ix], col = 'black')

# Estimate of beta_0

plot(1, type = "n", xlim = c(0,1), ylim = c(-1.01, 1.01) * ylim_all, main = expression("Estimated intercept"), ylab = expression(beta[0]), xlab = expression(Z[1]))

ix0 <- which(Z_all[,2] < 0.5)
ix1 <- which(Z_all[,2] >= 0.5)
plot0_ix0 <- ix0[order(Z_all[ix0,1])]
plot0_ix1 <- ix1[order(Z_all[ix1,1])]
polygon(c(Z_all[plot0_ix0,1], rev(Z_all[plot0_ix0,1])),
        c(beta0_hat[plot0_ix0,"L95"], rev(beta0_hat[plot0_ix0,"U95"])),
        col = rgb(0,0,1,1/3), border = NA)
polygon(c(Z_all[plot0_ix1,1], rev(Z_all[plot0_ix1,1])),
        c(beta0_hat[plot0_ix1,"L95"], rev(beta0_hat[plot0_ix1,"U95"])),
        col = rgb(0,0,1,1/3), border = NA)

lines(Z_all[plot0_ix0,1], beta0_all[plot0_ix0], col = 'black')
lines(Z_all[plot0_ix0,1], beta0_hat[plot0_ix0], col = 'blue')

lines(Z_all[plot0_ix1,1], beta0_all[plot0_ix1], col = 'black', lty = 2)
lines(Z_all[plot0_ix1,1], beta0_hat[plot0_ix1], col = 'blue', lty = 2)
legend("bottomright", cex = 0.8, legend = expression(Z[2] >= 0.5), lty = 2, bty = "n")
legend("topleft", cex = 0.8, legend = expression(Z[2] < 0.5), lty = 1, bty = "n")

# Estimate of beta1
plot(1, type = "n", xlim = c(0,1), ylim = c(-1.01, 1.01) * ylim_all, main = expression("Estimated effect of"~X[1]), 
     ylab = expression(beta[1]), xlab = expression(Z[1]))
plot1_ix <- order(Z_all[,1])

polygon(c(Z_all[plot1_ix,1], rev(Z_all[plot1_ix,1])),
        c(beta1_hat[plot1_ix,"L95"], rev(beta1_hat[plot1_ix,"U95"])),
        col = rgb(0,0,1,1/3), border = NA)

lines(Z_all[plot1_ix,1], beta1_all[plot1_ix], col = 'black')
lines(Z_all[plot1_ix,1], beta1_hat[plot1_ix], col = 'blue')

# Estimate of beta2
plot(1, type = "n", xlim = c(0,1), ylim = c(-1.01, 1.01) * ylim_all, main = expression("Estimated effect of"~X[2]), 
     ylab = expression(beta[2]), xlab = expression(Z[1]))

ix0 <- which(Z_all[,1] > 0.6)
ix1 <- which(Z_all[,1] < 0.25)
ix2 <- which(Z_all[,1] >= 0.25 & Z_all[,1] <= 0.6)


plot2_ix0 <- ix0[order(Z_all[ix0,1])]
plot2_ix1 <- ix1[order(Z_all[ix1,1])]
plot2_ix2 <- ix2[order(Z_all[ix2,1])]


polygon(c(Z_all[plot2_ix0,1], rev(Z_all[plot2_ix0,1])),
        c(beta2_hat[plot2_ix0,"L95"], rev(beta2_hat[plot2_ix0,"U95"])),
        col = rgb(0,0,1,1/3), border = NA)
polygon(c(Z_all[plot2_ix1,1], rev(Z_all[plot2_ix1,1])),
        c(beta2_hat[plot2_ix1,"L95"], rev(beta2_hat[plot2_ix1,"U95"])),
        col = rgb(0,0,1,1/3), border = NA)
polygon(c(Z_all[plot2_ix2,1], rev(Z_all[plot2_ix2,1])),
        c(beta2_hat[plot2_ix2,"L95"], rev(beta2_hat[plot2_ix2,"U95"])),
        col = rgb(0,0,1,1/3), border = NA)

lines(Z_all[plot2_ix0,1], beta2_all[plot2_ix0], col = 'black')
lines(Z_all[plot2_ix1,1], beta2_all[plot2_ix1], col = 'black')
lines(Z_all[plot2_ix2,1], beta2_all[plot2_ix2], col = 'black')

lines(Z_all[plot2_ix0,1], beta2_hat[plot2_ix0], col = 'blue')
lines(Z_all[plot2_ix1,1], beta2_hat[plot2_ix1], col = 'blue')
lines(Z_all[plot2_ix2,1], beta2_hat[plot2_ix2], col = 'blue')

# Estimate of beta3
plot(1, type = "n", xlim = c(0,1), ylim = c(-1.01,1.01) * ylim_all, main = expression("Estimated effect of"~X[3]), 
     ylab = expression(beta[3]), xlab = expression(Z[1]))
plot1_ix <- order(Z_all[,1])

polygon(c(Z_all[plot1_ix,1], rev(Z_all[plot1_ix,1])),
        c(beta3_hat[plot1_ix,"L95"], rev(beta3_hat[plot1_ix,"U95"])),
        col = rgb(0,0,1,1/3), border = NA)

lines(Z_all[plot1_ix,1], beta3_all[plot1_ix], col = 'black')
lines(Z_all[plot1_ix,1], beta3_hat[plot1_ix], col = 'blue')

dev.off()