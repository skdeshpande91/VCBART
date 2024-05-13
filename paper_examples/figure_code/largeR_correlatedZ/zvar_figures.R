load("p5_z_variants_results.RData")
beta_mse_mean <- apply(beta_mse_test, FUN = mean, MARGIN = c(1,3))
beta_cov_mean <- apply(beta_cov_test, FUN = mean, MARGIN = c(1,3))


pdf("../../figures/zvar_beta_mse.pdf", width = 4, height = 4)
par(mar = c(3, 1, 2, 1), mgp = c(1.8, 0.5, 0), cex.lab = 1, cex.axis = 1, cex.main = 1)
boxplot(beta_mse_mean, use.cols = TRUE, horizontal = TRUE,
        main = expression("Covariate effect recovery"), xlab = "", 
        yaxt = "n",
        pch = 16, cex = 0.5, medlwd = 0.5, lwd = 0.75)
abline(h = 0.5 + c(4, 8), col = 'gray', lwd = 2, lty = 2)

mtext(text = "MSE", side = 1, line = 2, cex = par('cex'))
text(x = 0.2, y = 1, labels = "rho = 0", cex = 0.85 * par('usr')["cex"])
text(x = 0.2, y = 2, labels = "rho = 0.5", cex = 0.85 * par('usr')["cex"])
text(x = 0.2, y = 3, labels = "rho = 0.75", cex = 0.85 * par('usr')["cex"])
text(x = 0.2, y = 4, labels = "rho = 0.9", cex = 0.85 * par('usr')["cex"])
text(x = 0.45, y = 1.3, labels = "R = 20", cex = par('usr')["cex"])


text(x = 0.2, y = 5, labels = "rho = 0", cex = 0.85 * par('usr')["cex"])
text(x = 0.2, y = 6, labels = "rho = 0.5", cex = 0.85 * par('usr')["cex"])
text(x = 0.2, y = 7, labels = "rho = 0.75", cex = 0.85 * par('usr')["cex"])
text(x = 0.2, y = 8, labels = "rho = 0.9", cex = 0.85 * par('usr')["cex"])
text(x = 0.45, y = 5.3, labels = "R = 50", cex = par("usr")["cex"])

text(x = 0.32, y = 9, labels = "rho = 0", cex = 0.85 * par('usr')["cex"])
text(x = 0.32, y = 10, labels = "rho = 0.5", cex = 0.85 * par('usr')["cex"])
text(x = 0.4, y = 11, labels = "rho = 0.75", cex = 0.85 * par('usr')["cex"])
text(x = 0.32, y = 12, labels = "rho = 0.9", cex = 0.85 * par('usr')["cex"])
text(x = 0.45, y = 9.3, labels = "R = 100", cex = par("usr")["cex"])

dev.off()


###################
# beta coverage
pdf("../../figures/zvar_beta_cov.pdf", width = 4, height = 4)
par(mar = c(3, 1, 2, 1), mgp = c(1.8, 0.5, 0), cex.lab = 1, cex.axis = 1, cex.main = 1)
boxplot(beta_cov_mean, use.cols = TRUE, horizontal = TRUE,
        main = expression("Uncertainty interval coverage"), xlab = "", 
        yaxt = "n", ylim = c(0.90, 1),
        pch = 16, cex = 0.5, medlwd = 0.5, lwd = 0.75)
abline(h = 0.5 + c(4, 8), col = 'gray', lwd = 2, lty = 2)
abline(v = 0.95, col = 'red', lty = 2)

mtext(text = "Coverage", side = 1, line = 2, cex = par('cex'))
text(x = 0.96, y = 1, labels = "rho = 0", cex = 0.85 * par('usr')["cex"])
text(x = 0.96, y = 2, labels = "rho = 0.5", cex = 0.85 * par('usr')["cex"])
text(x = 0.96, y = 3, labels = "rho = 0.75", cex = 0.85 * par('usr')["cex"])
text(x = 0.96, y = 4, labels = "rho = 0.9", cex = 0.85 * par('usr')["cex"])
text(x = 0.92, y = 1.3, labels = "R = 20", cex = par('usr')["cex"])


text(x = 0.96, y = 5, labels = "rho = 0", cex = 0.85 * par('usr')["cex"])
text(x = 0.96, y = 6, labels = "rho = 0.5", cex = 0.85 * par('usr')["cex"])
text(x = 0.96, y = 7, labels = "rho = 0.75", cex = 0.85 * par('usr')["cex"])
text(x = 0.96, y = 8, labels = "rho = 0.9", cex = 0.85 * par('usr')["cex"])
text(x = 0.92, y = 5.3, labels = "R = 50", cex = par("usr")["cex"])

text(x = 0.96, y = 9, labels = "rho = 0", cex = 0.85 * par('usr')["cex"])
text(x = 0.96, y = 10, labels = "rho = 0.5", cex = 0.85 * par('usr')["cex"])
text(x = 0.96, y = 11, labels = "rho = 0.75", cex = 0.85 * par('usr')["cex"])
text(x = 0.96, y = 12, labels = "rho = 0.9", cex = 0.85 * par('usr')["cex"])
text(x = 0.92, y = 9.3, labels = "R = 100", cex = par("usr")["cex"])
dev.off()


########
# Ystar rmse
pdf("../../figures/zvar_ystar_rmse.pdf", width = 4, height = 4)
par(mar = c(3, 1, 2, 1), mgp = c(1.8, 0.5, 0), cex.lab = 1, cex.axis = 1, cex.main = 1)
boxplot(ystar_rmse_test, use.cols = TRUE, horizontal = TRUE,
        main = expression("Predictive performance"), xlab = "", 
        yaxt = "n",
        pch = 16, cex = 0.5, medlwd = 0.5, lwd = 0.75)
abline(h = 0.5 + c(4, 8), col = 'gray', lwd = 2, lty = 2)

mtext(text = "RMSE", side = 1, line = 2, cex = par('cex'))
text(x = 1.8, y = 1, labels = "rho = 0", cex = 0.85 * par('usr')["cex"])
text(x = 1.8, y = 2, labels = "rho = 0.5", cex = 0.85 * par('usr')["cex"])
text(x = 1.8, y = 3, labels = "rho = 0.75", cex = 0.85 * par('usr')["cex"])
text(x = 1.8, y = 4, labels = "rho = 0.9", cex = 0.85 * par('usr')["cex"])
text(x = 2.0, y = 1.3, labels = "R = 20", cex = par('usr')["cex"])


text(x = 1.8, y = 5, labels = "rho = 0", cex = 0.85 * par('usr')["cex"])
text(x = 1.8, y = 6, labels = "rho = 0.5", cex = 0.85 * par('usr')["cex"])
text(x = 1.8, y = 7, labels = "rho = 0.75", cex = 0.85 * par('usr')["cex"])
text(x = 1.8, y = 8, labels = "rho = 0.9", cex = 0.85 * par('usr')["cex"])
text(x = 2, y = 5.3, labels = "R = 50", cex = par("usr")["cex"])

text(x = 1.9, y = 9, labels = "rho = 0", cex = 0.85 * par('usr')["cex"])
text(x = 1.9, y = 10, labels = "rho = 0.5", cex = 0.85 * par('usr')["cex"])
text(x = 1.9, y = 11, labels = "rho = 0.75", cex = 0.85 * par('usr')["cex"])
text(x = 1.9, y = 12, labels = "rho = 0.9", cex = 0.85 * par('usr')["cex"])
text(x = 2, y = 9.3, labels = "R = 100", cex = par("usr")["cex"])
dev.off()
########
# Ystar coverage
pdf("../../figures/zvar_ystar_cov.pdf", width = 4, height = 4)
par(mar = c(3, 1, 2, 1), mgp = c(1.8, 0.5, 0), cex.lab = 1, cex.axis = 1, cex.main = 1)
boxplot(ystar_cov_test, use.cols = TRUE, horizontal = TRUE,
        main = expression("Predictive interval coverage"), xlab = "", 
        yaxt = "n", ylim = c(0.90, 1),
        pch = 16, cex = 0.5, medlwd = 0.5, lwd = 0.75)
abline(h = 0.5 + c(4, 8), col = 'gray', lwd = 2, lty = 2)
abline(v = 0.95, col = 'red', lty = 2)

mtext(text = "Coverage", side = 1, line = 2, cex = par('cex'))
text(x = 0.94, y = 1, labels = "rho = 0", cex = 0.85 * par('usr')["cex"])
text(x = 0.94, y = 2, labels = "rho = 0.5", cex = 0.85 * par('usr')["cex"])
text(x = 0.94, y = 3, labels = "rho = 0.75", cex = 0.85 * par('usr')["cex"])
text(x = 0.94, y = 4, labels = "rho = 0.9", cex = 0.85 * par('usr')["cex"])
text(x = 0.92, y = 1.3, labels = "R = 20", cex = par('usr')["cex"])


text(x = 0.94, y = 5, labels = "rho = 0", cex = 0.85 * par('usr')["cex"])
text(x = 0.94, y = 6, labels = "rho = 0.5", cex = 0.85 * par('usr')["cex"])
text(x = 0.94, y = 7, labels = "rho = 0.75", cex = 0.85 * par('usr')["cex"])
text(x = 0.94, y = 8, labels = "rho = 0.9", cex = 0.85 * par('usr')["cex"])
text(x = 0.92, y = 5.3, labels = "R = 50", cex = par("usr")["cex"])

text(x = 0.94, y = 9, labels = "rho = 0", cex = 0.85 * par('usr')["cex"])
text(x = 0.94, y = 10, labels = "rho = 0.5", cex = 0.85 * par('usr')["cex"])
text(x = 0.94, y = 11, labels = "rho = 0.75", cex = 0.85 * par('usr')["cex"])
text(x = 0.94, y = 12, labels = "rho = 0.9", cex = 0.85 * par('usr')["cex"])
text(x = 0.92, y = 9.3, labels = "R = 100", cex = par("usr")["cex"])
dev.off()

# F1
pdf("../../figures/zvar_f1.pdf", width = 4, height = 4)
par(mar = c(3, 1, 2, 1), mgp = c(1.8, 0.5, 0), cex.lab = 1, cex.axis = 1, cex.main = 1)
boxplot(varsel[,"f1",], use.cols = TRUE, horizontal = TRUE,
        main = expression("Modifier Selection"), xlab = "", 
        yaxt = "n",
        pch = 16, cex = 0.5, medlwd = 0.5, lwd = 0.75)
abline(h = 0.5 + c(4, 8), col = 'gray', lwd = 2, lty = 2)

mtext(text = "F1", side = 1, line = 2, cex = par('cex'))
text(x = 0.6, y = 1, labels = "rho = 0", cex = 0.85 * par('usr')["cex"])
text(x = 0.6, y = 2, labels = "rho = 0.5", cex = 0.85 * par('usr')["cex"])
text(x = 0.6, y = 3, labels = "rho = 0.75", cex = 0.85 * par('usr')["cex"])
text(x = 0.6, y = 4, labels = "rho = 0.9", cex = 0.85 * par('usr')["cex"])
text(x = 0.5, y = 1.3, labels = "R = 20", cex = par('usr')["cex"])


text(x = 0.6, y = 5, labels = "rho = 0", cex = 0.85 * par('usr')["cex"])
text(x = 0.6, y = 6, labels = "rho = 0.5", cex = 0.85 * par('usr')["cex"])
text(x = 0.6, y = 7, labels = "rho = 0.75", cex = 0.85 * par('usr')["cex"])
text(x = 0.6, y = 8, labels = "rho = 0.9", cex = 0.85 * par('usr')["cex"])
text(x = 0.5, y = 5.3, labels = "R = 50", cex = par("usr")["cex"])

text(x = 0.6, y = 9, labels = "rho = 0", cex = 0.85 * par('usr')["cex"])
text(x = 0.6, y = 10, labels = "rho = 0.5", cex = 0.85 * par('usr')["cex"])
text(x = 0.6, y = 11, labels = "rho = 0.75", cex = 0.85 * par('usr')["cex"])
text(x = 0.6, y = 12.3, labels = "rho = 0.9", cex = 0.85 * par('usr')["cex"])
text(x = 0.5, y = 9.3, labels = "R = 100", cex = par("usr")["cex"])
dev.off()

# F1
pdf("../../figures/zvar_timing.pdf", width = 4, height = 4)
par(mar = c(3, 1, 2, 1), mgp = c(1.8, 0.5, 0), cex.lab = 1, cex.axis = 1, cex.main = 1)
boxplot(timing, use.cols = TRUE, horizontal = TRUE,
        main = expression("Run time"), xlab = "", 
        yaxt = "n", ylim = c(50, 200),
        pch = 16, cex = 0.5, medlwd = 0.5, lwd = 0.75)
abline(h = 0.5 + c(4, 8), col = 'gray', lwd = 2, lty = 2)

mtext(text = "Seconds", side = 1, line = 2, cex = par('cex'))
text(x = 55, y = 1, labels = "rho = 0", cex = 0.85 * par('usr')["cex"])
text(x = 55, y = 2, labels = "rho = 0.5", cex = 0.85 * par('usr')["cex"])
text(x = 55, y = 3, labels = "rho = 0.75", cex = 0.85 * par('usr')["cex"])
text(x = 55, y = 4, labels = "rho = 0.9", cex = 0.85 * par('usr')["cex"])
text(x = 190, y = 1.3, labels = "R = 20", cex = par('usr')["cex"])

text(x = 55, y = 5, labels = "rho = 0", cex = 0.85 * par('usr')["cex"])
text(x = 55, y = 6, labels = "rho = 0.5", cex = 0.85 * par('usr')["cex"])
text(x = 55, y = 7, labels = "rho = 0.75", cex = 0.85 * par('usr')["cex"])
text(x = 55, y = 8, labels = "rho = 0.9", cex = 0.85 * par('usr')["cex"])
text(x = 190, y = 5.3, labels = "R = 50", cex = par("usr")["cex"])

text(x = 55, y = 9, labels = "rho = 0", cex = 0.85 * par('usr')["cex"])
text(x = 55, y = 10, labels = "rho = 0.5", cex = 0.85 * par('usr')["cex"])
text(x = 55, y = 11, labels = "rho = 0.75", cex = 0.85 * par('usr')["cex"])
text(x = 55, y = 12.3, labels = "rho = 0.9", cex = 0.85 * par('usr')["cex"])
text(x = 190, y = 9.3, labels = "R = 100", cex = par("usr")["cex"])
dev.off()
