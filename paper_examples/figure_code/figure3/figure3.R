load("p5R20_results.RData")
vc_methods <- c("vcbart", "btvcm", "tvc", "ks", "lm")
ml_methods <- c("bart", "ert")
all_methods <- c(vc_methods, ml_methods)

beta_mse_mean <- apply(beta_mse_test, FUN = mean, MARGIN = c(1,3))
beta_cov_mean <- apply(beta_cov_test, FUN = mean, MARGIN = c(1,3))

pdf("p5R20_beta_mse.pdf", width = 3, height = 3)
par(mar = c(3, 1, 2, 1), mgp = c(1.8, 0.5, 0), cex.lab = 1, cex.axis = 1, cex.main = 1)
boxplot(beta_mse_mean[,vc_methods], use.cols = TRUE, horizontal = TRUE,
        main = expression("Covariate effect recovery"), xlab = "", 
        yaxt = "n",
        pch = 16, cex = 0.5, medlwd = 0.5, lwd = 0.75)
mtext(text = "MSE", side = 1, line = 2, cex = par('cex'))
text(x = 2, y = 1, labels = "VCBART", cex = 0.85*par('cex'))
text(x = 5, y = 2, labels = "BTVCM", cex = 0.85*par('cex'))
text(x = 4, y = 3, labels = "TVCM",cex = 0.85*par('cex'))
text(x = 0.5, y = 4, labels = "KS",cex = 0.85*par('cex'))
text(x = 5, y = 5, labels = "lm", cex = 0.85*par('cex'))
dev.off()

pdf("p5R20_beta_cov.pdf", width = 3, height = 3)
par(mar = c(3, 1, 2, 1), mgp = c(1.8, 0.5, 0), cex.lab = 1, cex.axis = 1, cex.main = 1)
boxplot(beta_cov_mean[,vc_methods], use.cols = TRUE, horizontal = TRUE,
        main = expression("Uncertainty interval coverage"), xlab = "", 
        yaxt = "n", 
        ylim = c(0,1),
        pch = 16, cex = 0.5, medlwd = 0.5, lwd = 0.75)
mtext(text = "Coverage", side = 1, line = 2, cex =  par('cex'))
text(x = 0.8, y = 1, labels = "VCBART", cex = 0.85*par('cex'))
text(x = 0.8, y = 2, labels = "BTVCM", cex = 0.85*par('cex'))
text(x = 0.8, y = 3, labels = "TVCM",cex = 0.85*par('cex'))
text(x = 0.65, y = 4, labels = "KS",cex = 0.85*par('cex'))
text(x = 0.65, y = 5, labels = "lm", cex = 0.85*par('cex'))
abline(v = 0.95, col = 'red', lty = 2)
dev.off()

pdf("p5R20_ystar_rmse.pdf", width = 3, height = 3)
par(mar = c(3, 1, 2, 1), mgp = c(1.8, 0.5, 0), cex.lab = 1, cex.axis = 1, cex.main = 1)
boxplot.matrix(ystar_rmse_test[,all_methods], use.cols = TRUE, horizontal = TRUE,
               main = expression("Predictive performance"), xlab = "", 
               yaxt = "n",
               pch = 16, cex = 0.5, medlwd = 0.5)
mtext(text = "RMSE", side = 1, line = 2, cex = par('cex'))
text(x = 3, y = 1, labels = "VCBART", cex = 0.85*par('cex'))
text(x = 6, y = 2, labels = "BTVCM", cex = 0.85*par('cex'))
text(x = 6, y = 3, labels = "TVCM",cex = 0.85*par('cex'))
text(x = 6, y = 4, labels = "KS",cex = 0.85*par('cex'))
text(x = 5, y = 5, labels = "lm", cex = 0.85*par('cex'))
text(x = 6, y = 6, labels = "BART", cex = 0.85*par('cex'))
text(x = 6, y = 7, labels = "ERT", cex = 0.85*par('cex'))
dev.off()


pdf("p5R20_ystar_cov.pdf", width = 3, height = 3)
par(mar = c(3, 1, 2, 1), mgp = c(1.8, 0.5, 0), cex.lab = 1, cex.axis = 1, cex.main = 1)
boxplot.matrix(ystar_cov_test[,all_methods], use.cols = TRUE, horizontal = TRUE,
               main = expression("Predictive interval coverage"), xlab = "", 
               yaxt = "n",
               pch = 16, cex = 0.5, medlwd = 0.5)
mtext(text = "Coverage", side = 1, line = 2, cex = 1 * par('cex'))
text(x = 0.8, y = 1, labels = "VCBART", cex = 0.85*par('cex'))
text(x = 0.7, y = 2, labels = "BTVCM", cex = 0.85*par('cex'))
text(x = 0.6, y = 3, labels = "TVCM",cex = 0.85*par('cex'))
text(x = 0.9, y = 4, labels = "KS",cex = 0.85*par('cex'))
text(x = 0.8, y = 5, labels = "lm", cex = 0.85*par('cex'))
text(x = 0.65, y = 6, labels = "BART", cex = 0.85*par('cex'))
text(x = 0.6, y = 7, labels = "ERT", cex = 0.85*par('cex'))
abline(v = 0.95, col = 'red', lty = 2)
dev.off()