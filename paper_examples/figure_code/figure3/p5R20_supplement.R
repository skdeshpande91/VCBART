load("p5R20_results.RData")
vc_methods <- c("vcbart", "btvcm", "tvc", "ks", "lm")

########

cex_lab <- 1
cex_axis <- 1
cex_main <- 1.5

pdf("beta0_mse_supp.pdf", width = 4.5, height = 4.5)
par(mar = c(3, 1, 2, 1), mgp = c(1.8, 0.5, 0), 
    cex.lab = cex_lab, cex.axis = cex_axis, cex.main = cex_main)
text_cex = par('cex')

boxplot(beta_mse_test[,1,vc_methods], use.cols = TRUE, horizontal = TRUE,
        main = expression("Intercept"), xlab = "", 
        yaxt = "n",
        pch = 16, cex = 0.5, medlwd = 0.5, lwd = 0.75)
mtext(text = "MSE", side = 1, line = 2, cex = text_cex)
text(x = 4, y = 1, labels = "VCBART", cex = text_cex)
text(x = 4, y = 2, labels = "BTVCM", cex = text_cex)
text(x = 4, y = 3, labels = "TVCM",cex = text_cex)
text(x = 4, y = 4, labels = "KS",cex = text_cex)
text(x = 4, y = 5, labels = "lm", cex = text_cex)
dev.off()


pdf("beta1_mse_supp.pdf", width = 4.5, height = 4.5)
par(mar = c(3, 1, 2, 1), mgp = c(1.8, 0.5, 0), 
    cex.lab = cex_lab, cex.axis = cex_axis, cex.main = cex_main)
text_cex <- par("cex")
boxplot(beta_mse_test[,2,vc_methods], use.cols = TRUE, horizontal = TRUE,
        main = expression(beta[1]), xlab = "", 
        yaxt = "n",
        pch = 16, cex = 0.5, medlwd = 0.5, lwd = 0.75)
mtext(text = "MSE", side = 1, line = 2, cex = text_cex)
text(x = 2, y = 1, labels = "VCBART", cex = text_cex)
text(x = 2, y = 2, labels = "BTVCM", cex = text_cex)
text(x = 2, y = 3, labels = "TVCM",cex = text_cex)
text(x = 2.5, y = 4.2, labels = "KS",cex = text_cex)
text(x = 2, y = 5, labels = "lm", cex = text_cex)
dev.off()


pdf("beta2_mse_supp.pdf", width = 4.5, height = 4.5)
par(mar = c(3, 1, 2, 1), mgp = c(1.8, 0.5, 0), 
    cex.lab = cex_lab, cex.axis = cex_axis, cex.main = cex_main)
text_cex <- par("cex")

boxplot(beta_mse_test[,3,vc_methods], use.cols = TRUE, horizontal = TRUE,
        main = expression(beta[2]), xlab = "", 
        yaxt = "n",
        pch = 16, cex = 0.5, medlwd = 0.5, lwd = 0.75)
mtext(text = "MSE", side = 1, line = 2, cex = text_cex)
text(x = 4.5, y = 1, labels = "VCBART", cex = text_cex)
text(x = 4.5, y = 2, labels = "BTVCM", cex = text_cex)
text(x = 4.5, y = 3, labels = "TVCM",cex = text_cex)
text(x = 5.75, y = 4.2, labels = "KS",cex = text_cex)
text(x = 4.5, y = 5, labels = "lm", cex = text_cex)
dev.off()

pdf("beta3_mse_supp.pdf", width = 4.5, height = 4.5)
par(mar = c(3, 1, 2, 1), mgp = c(1.8, 0.5, 0), 
    cex.lab = cex_lab, cex.axis = cex_axis, cex.main = cex_main)
text_cex <- par("cex")
boxplot(beta_mse_test[,4,vc_methods], use.cols = TRUE, horizontal = TRUE,
        main = expression(beta[3]), xlab = "", 
        yaxt = "n",
        pch = 16, cex = 0.5, medlwd = 0.5, lwd = 0.75)
mtext(text = "MSE", side = 1, line = 2, cex = text_cex)
text(x = 3, y = 1, labels = "VCBART", cex = text_cex)
text(x = 3, y = 2, labels = "BTVCM", cex = text_cex)
text(x = 3, y = 3, labels = "TVCM",cex = text_cex)
text(x = 3.5, y = 4.2, labels = "KS",cex = text_cex)
text(x = 3, y = 5, labels = "lm", cex = text_cex)
dev.off()


pdf("beta4_mse_supp.pdf", width = 4.5, height = 4.5)
par(mar = c(3, 1, 2, 1), mgp = c(1.8, 0.5, 0), 
    cex.lab = cex_lab, cex.axis = cex_axis, cex.main = cex_main)
text_cex <- par("cex")
boxplot(beta_mse_test[,5,vc_methods], use.cols = TRUE, horizontal = TRUE,
        main = expression(beta[4]), xlab = "", 
        yaxt = "n",
        pch = 16, cex = 0.5, medlwd = 0.5, lwd = 0.75)
mtext(text = "MSE", side = 1, line = 2, cex = text_cex)
text(x = 5, y = 1, labels = "VCBART", cex = text_cex)
text(x = 3, y = 2, labels = "BTVCM", cex = text_cex)
text(x = 3, y = 3, labels = "TVCM",cex = text_cex)
text(x = 2, y = 4, labels = "KS",cex = text_cex)
text(x = 3, y = 5, labels = "lm", cex = text_cex)
dev.off()

pdf("beta5_mse_supp.pdf", width = 4.5, height = 4.5)
par(mar = c(3, 1, 2, 1), mgp = c(1.8, 0.5, 0), 
    cex.lab = cex_lab, cex.axis = cex_axis, cex.main = cex_main)
text_cex <- par("cex")
boxplot(beta_mse_test[,6,vc_methods], use.cols = TRUE, horizontal = TRUE,
        main = expression(beta[5]), xlab = "", 
        yaxt = "n",
        pch = 16, cex = 0.5, medlwd = 0.5, lwd = 0.75)
mtext(text = "MSE", side = 1, line = 2, cex = text_cex)
text(x = 1.5, y = 1, labels = "VCBART", cex = text_cex)
text(x = 1.5, y = 2, labels = "BTVCM", cex = text_cex)
text(x = 2, y = 3, labels = "TVCM",cex = text_cex)
text(x = 2.75, y = 4.2, labels = "KS",cex = text_cex)
text(x = 1.5, y = 5, labels = "lm", cex = text_cex)
dev.off()


###############
# Coverages
cex_lab <- 1
cex_axis <- 1
cex_main <- 1.5

pdf("beta0_cov_supp.pdf", width = 4.5, height = 4.5)
par(mar = c(3, 1, 2, 1), mgp = c(1.8, 0.5, 0), 
    cex.lab = cex_lab, cex.axis = cex_axis, cex.main = cex_main)
text_cex = par('cex')

boxplot(beta_cov_test[,1,vc_methods], use.cols = TRUE, horizontal = TRUE,
        main = expression("Intercept"), xlab = "", 
        yaxt = "n", ylim = c(0,1),
        pch = 16, cex = 0.5, medlwd = 0.5, lwd = 0.75)
mtext(text = "Coverage", side = 1, line = 2, cex = text_cex)
text(x = 0.8, y = 1, labels = "VCBART", cex = text_cex)
text(x = 0.8, y = 2, labels = "BTVCM", cex = text_cex)
text(x = 0.8, y = 3, labels = "TVCM",cex = text_cex)
text(x = 0.5, y = 4, labels = "KS",cex = text_cex)
text(x = 0.5, y = 5, labels = "lm", cex = text_cex)
abline(v = 0.95, col = 'red', lty = 2)
dev.off()


pdf("beta1_cov_supp.pdf", width = 4.5, height = 4.5)
par(mar = c(3, 1, 2, 1), mgp = c(1.8, 0.5, 0), 
    cex.lab = cex_lab, cex.axis = cex_axis, cex.main = cex_main)
text_cex <- par("cex")
boxplot(beta_cov_test[,2,vc_methods], use.cols = TRUE, horizontal = TRUE,
        main = expression(beta[1]), xlab = "", 
        yaxt = "n", ylim = c(0,1),
        pch = 16, cex = 0.5, medlwd = 0.5, lwd = 0.75)
mtext(text = "Coverage", side = 1, line = 2, cex = text_cex)
text(x = 0.8, y = 1, labels = "VCBART", cex = text_cex)
text(x = 0.6, y = 2, labels = "BTVCM", cex = text_cex)
text(x = 0.8, y = 3, labels = "TVCM",cex = text_cex)
text(x = 0.8, y = 4, labels = "KS",cex = text_cex)
text(x = 0.85, y = 5, labels = "lm", cex = text_cex)
abline(v = 0.95, col = 'red', lty = 2)

dev.off()


pdf("beta2_cov_supp.pdf", width = 4.5, height = 4.5)
par(mar = c(3, 1, 2, 1), mgp = c(1.8, 0.5, 0), 
    cex.lab = cex_lab, cex.axis = cex_axis, cex.main = cex_main)
text_cex <- par("cex")

boxplot(beta_cov_test[,3,vc_methods], use.cols = TRUE, horizontal = TRUE,
        main = expression(beta[2]), xlab = "", 
        yaxt = "n", ylim = c(0,1),
        pch = 16, cex = 0.5, medlwd = 0.5, lwd = 0.75)
mtext(text = "Coverage", side = 1, line = 2, cex = text_cex)
text(x = 0.8, y = 1, labels = "VCBART", cex = text_cex)
text(x = 0.85, y = 2, labels = "BTVCM", cex = text_cex)
text(x = 0.8, y = 3, labels = "TVCM",cex = text_cex)
text(x = 0.5, y = 4, labels = "KS",cex = text_cex)
text(x = 0.55, y = 5, labels = "lm", cex = text_cex)
abline(v = 0.95, col = 'red', lty = 2)

dev.off()

pdf("beta3_cov_supp.pdf", width = 4.5, height = 4.5)
par(mar = c(3, 1, 2, 1), mgp = c(1.8, 0.5, 0), 
    cex.lab = cex_lab, cex.axis = cex_axis, cex.main = cex_main)
text_cex <- par("cex")
boxplot(beta_cov_test[,4,vc_methods], use.cols = TRUE, horizontal = TRUE,
        main = expression(beta[3]), xlab = "", 
        yaxt = "n", ylim = c(0,1),
        pch = 16, cex = 0.5, medlwd = 0.5, lwd = 0.75)
mtext(text = "Coverage", side = 1, line = 2, cex = text_cex)
text(x = 0.8, y = 1, labels = "VCBART", cex = text_cex)
text(x = 0.5, y = 2, labels = "BTVCM", cex = text_cex)
text(x = 0.7, y = 3, labels = "TVCM",cex = text_cex)
text(x = 0.7, y = 4.2, labels = "KS",cex = text_cex)
text(x = 0.8, y = 5, labels = "lm", cex = text_cex)
dev.off()


pdf("beta4_cov_supp.pdf", width = 4.5, height = 4.5)
par(mar = c(3, 1, 2, 1), mgp = c(1.8, 0.5, 0), 
    cex.lab = cex_lab, cex.axis = cex_axis, cex.main = cex_main)
text_cex <- par("cex")
boxplot(beta_cov_test[,5,vc_methods], use.cols = TRUE, horizontal = TRUE,
        main = expression(beta[4]), xlab = "", 
        yaxt = "n", ylim = c(0,1),
        pch = 16, cex = 0.5, medlwd = 0.5, lwd = 0.75)
mtext(text = "Coverage", side = 1, line = 2, cex = text_cex)
text(x = 0.8, y = 1, labels = "VCBART", cex = text_cex)
text(x = 0.5, y = 2, labels = "BTVCM", cex = text_cex)
text(x = 0.6, y = 3, labels = "TVCM",cex = text_cex)
text(x = 0.85, y = 4, labels = "KS",cex = text_cex)
text(x = 0.2, y = 5, labels = "lm", cex = text_cex)
abline(v = 0.95, col = 'red', lty = 2)
dev.off()

pdf("beta5_cov_supp.pdf", width = 4.5, height = 4.5)
par(mar = c(3, 1, 2, 1), mgp = c(1.8, 0.5, 0), 
    cex.lab = cex_lab, cex.axis = cex_axis, cex.main = cex_main)
text_cex <- par("cex")
boxplot(beta_cov_test[,6,vc_methods], use.cols = TRUE, horizontal = TRUE,
        main = expression(beta[5]), xlab = "", 
        yaxt = "n", ylim = c(0,1),
        pch = 16, cex = 0.5, medlwd = 0.5, lwd = 0.75)
mtext(text = "Coverage", side = 1, line = 2, cex = text_cex)
text(x = 0.8, y = 1, labels = "VCBART", cex = text_cex)
text(x = 0.45, y = 2, labels = "BTVCM", cex = text_cex)
text(x = 0.8, y = 3, labels = "TVCM",cex = text_cex)
text(x = 0.65, y = 4.2, labels = "KS",cex = text_cex)
text(x = 0.7, y = 5, labels = "lm", cex = text_cex)
abline(v = 0.95, col = 'red', lty = 2)
dev.off()

