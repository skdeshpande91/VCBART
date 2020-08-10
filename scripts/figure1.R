
load("results/illustration_p5R20.RData")

############
# Figure 1 is 2 x 4
############

png("figures/p5R20_beta.png", width = 8, height = 4, units = "in", res = 400)
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0), cex.main = 0.9, cex.lab = 0.9, cex.axis = 0.9, mfrow = c(2,4))

# Plot 1: beta0
ix0 <- which(Z_plot[,2] < 0.5)
ix1 <- which(Z_plot[,2] >= 0.5)
plot0_ix0 <- ix0[order(Z_plot[ix0,1])]
plot0_ix1 <- ix1[order(Z_plot[ix1,1])]
ylim0 <- max(abs(c(beta_plot[,1], beta0_hat[,c("L95", "U95")])))

plot(1,type= "n", xlim = c(0,1), ylim = c(-1.05, 1.05) * ylim0,
     main = expression("Intercept"), ylab = expression(beta[0]), xlab = expression(Z[1]))
lines(Z_plot[plot0_ix0,1], beta_plot[plot0_ix0,1])
lines(Z_plot[plot0_ix1,1], beta_plot[plot0_ix1,1], col = 'black', lty = 2)
legend("bottomright", cex = 0.8, legend = expression(Z[2] == 0), lty = 2, bty = "n")
legend("topleft", cex = 0.8, legend = expression(Z[2] == 1), lty = 1, bty = "n")

# Plot 2: beta1
plot1_ix <- order(Z_plot[,1])
ylim1 <- max(abs(c(beta_plot[,2], beta1_hat[,c("L95", "U95")])))
plot(1, type = "n", xlim = c(0,1), ylim = c(-1.05, 1.05)*ylim1,
     main = expression("Effect of"~X[1]), xlab = expression(Z[1]), ylab = expression(beta[1]))
lines(Z_plot[plot1_ix,1], beta_plot[plot1_ix,2])

# Plot 3 beta2
ix0 <- which(Z_plot[,1] > 0.6)
ix1 <- which(Z_plot[,1] <= 0.25)
ix2 <- which(Z_plot[,1] > 0.25 & Z_plot[,1] < 0.6)

plot2_ix0 <- ix0[order(Z_plot[ix0,1])]
plot2_ix1 <- ix1[order(Z_plot[ix1,1])]
plot2_ix2 <- ix2[order(Z_plot[ix2,1])]

ylim2 <- max(abs(c(beta_plot[,3], beta2_hat[,c("L95", "U95")])))

plot(1, type = "n", xlim = c(0,1), ylim = c(-1.05, 1.05)*ylim2,
     main = expression("Effect of"~X[2]), xlab = expression(Z[1]), ylab = expression(beta[2]))

lines(Z_plot[plot2_ix0,1], beta_plot[plot2_ix0, 3])
lines(Z_plot[plot2_ix1,1], beta_plot[plot2_ix1, 3])
lines(Z_plot[plot2_ix2,1], beta_plot[plot2_ix2, 3])

# Plot 4: beta3
plot3_ix <- order(Z_plot[,1])
ylim3 <- max(abs(c(beta_plot[,4], beta3_hat[,c("L95", "U95")])))

plot(1, type = "n", xlim = c(0,1), ylim = c(-1.05, 1.05)*ylim3,
     main = expression("Effect of"~X[3]), xlab = expression(Z[1]), ylab = expression(beta[3]))
lines(Z_plot[plot3_ix,1], beta_plot[plot3_ix,4])

# Plot 5: beta0 & VC-BART estimates
# Plot 5: beta0 + VC-BART estimates
ix0 <- which(Z_plot[,2] < 0.5)
ix1 <- which(Z_plot[,2] >= 0.5)
plot0_ix0 <- ix0[order(Z_plot[ix0,1])]
plot0_ix1 <- ix1[order(Z_plot[ix1,1])]
ylim0 <- max(abs(c(beta_plot[,1], beta0_hat[,c("L95", "U95")])))

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
lines(Z_plot[plot0_ix0,1], beta0_hat[plot0_ix0,"MEAN"], col = 'blue')
lines(Z_plot[plot0_ix1,1], beta0_hat[plot0_ix1,"MEAN"], col = 'blue', lty = 2)


legend("bottomright", cex = 0.8, legend = expression(Z[2] == 0), lty = 2, bty = "n")
legend("topleft", cex = 0.8, legend = expression(Z[2] == 1), lty = 1, bty = "n")

# Plot 6: beta1 plus VC-BART results
plot1_ix <- order(Z_plot[,1])
ylim1 <- max(abs(c(beta_plot[,2], beta1_hat[,c("L95", "U95")])))

plot(1, type = "n", xlim = c(0,1), ylim = c(-1.05, 1.05)*ylim1,
     main = expression("Effect of"~X[1]), xlab = expression(Z[1]), ylab = expression(beta[1]))
polygon(c(Z_plot[plot1_ix,1], rev(Z_plot[plot1_ix,1])),
        c(beta1_hat[plot1_ix,"L95"], rev(beta1_hat[plot1_ix, "U95"])),
        col = rgb(0,0,1,1/3), border = NA)
lines(Z_plot[plot1_ix,1], beta1_hat[plot1_ix,"MEAN"], col = 'blue')
lines(Z_plot[plot1_ix,1], beta_plot[plot1_ix,2])

# Plot 7: beta2 & VC-BART results
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
        col = rgb(0,0,1,1/3), border = NA)
polygon(c(Z_plot[plot2_ix1,1], rev(Z_plot[plot2_ix1,1])),
        c(beta2_hat[plot2_ix1,"L95"], rev(beta2_hat[plot2_ix1, "U95"])),
        col = rgb(0,0,1,1/3), border = NA)
polygon(c(Z_plot[plot2_ix2,1], rev(Z_plot[plot2_ix2,1])),
        c(beta2_hat[plot2_ix2,"L95"], rev(beta2_hat[plot2_ix2, "U95"])),
        col = rgb(0,0,1,1/3), border = NA)
lines(Z_plot[plot2_ix0,1], beta2_hat[plot2_ix0, "MEAN"], col = 'blue')
lines(Z_plot[plot2_ix1,1], beta2_hat[plot2_ix1, "MEAN"], col = 'blue')
lines(Z_plot[plot2_ix2,1], beta2_hat[plot2_ix2, "MEAN"], col = 'blue')

lines(Z_plot[plot2_ix0,1], beta_plot[plot2_ix0, 3])
lines(Z_plot[plot2_ix1,1], beta_plot[plot2_ix1, 3])
lines(Z_plot[plot2_ix2,1], beta_plot[plot2_ix2, 3])

# Plot 8: beta3 & VC-BART results
plot3_ix <- order(Z_plot[,1])
ylim3 <- max(abs(c(beta_plot[,4], beta3_hat[,c("L95", "U95")])))

plot(1, type = "n", xlim = c(0,1), ylim = c(-1.05, 1.05)*ylim3,
     main = expression("Effect of"~X[3]), xlab = expression(Z[1]), ylab = expression(beta[3]))
polygon(c(Z_plot[plot3_ix,1], rev(Z_plot[plot3_ix,1])),
        c(beta3_hat[plot3_ix,"L95"], rev(beta3_hat[plot3_ix, "U95"])),
        col = rgb(0,0,1,1/3), border = NA)
lines(Z_plot[plot3_ix,1], beta3_hat[plot3_ix,"MEAN"], col = 'blue')
lines(Z_plot[plot3_ix,1], beta_plot[plot3_ix,4])

dev.off()
