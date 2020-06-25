load("results/sim_p5R20/results_p5R20.RData")

lin_method_names <- c("vcbart_adapt", "lm", "kernel_smoothing", "tvc", "boosted_tvcm")
method_names <- c("vcbart_adapt", "lm", "kernel_smoothing", "tvc", "boosted_tvcm", "bart", "extraTrees", "gbm")


beta_mse_mean <- apply(beta_mse_test, FUN = mean, MARGIN = c(1,3))
beta_cov_mean <- apply(beta_cov_test, FUN = mean, MARGIN = c(1,3))

png("~/Documents/Research/vc_bart/figures/p5R20_performance.png", 
    width = 6.5, height = 6.5, units = "in", res = 300)
par(mar = c(4,1,2,1), mgp = c(1.8, 0.5, 0), mfrow = c(2,2), cex.main = 1.1, cex.axis = 0.95, cex.lab = 1)
boxplot.matrix(beta_mse_mean[lin_method_names,], use.cols = FALSE, horizontal = TRUE,
               main = expression("Covariate effect recovery"), xlab = "", yaxt = "n",
               pch = 16, cex = 0.5, medlwd = 0.5, at = 1:5)
mtext(text = "MSE\n(a)", side = 1, line = 3, cex = 0.95)
text(x = 7.5, y = 1, labels = expression("VC-BART"~(rho == 0)), cex = 0.95)
text(x = 12.5, y = 2, labels = expression("lm"), cex = 0.95)
text(x = 2.5, y = 3, labels = expression("KS"), cex = 0.95)
text(x = 10, y = 4, labels = expression("TVCM"), cex = 0.95)
text(x = 10, y = 5, labels = expression("BTVCM"), cex = 0.95)

boxplot.matrix(beta_cov_mean[lin_method_names,], use.cols = FALSE, horizontal = TRUE,
               main = expression("Uncertainty interval coverage"~(beta)), xlab = "", yaxt = "n",
               pch = 16, cex = 0.5, medlwd = 0.5, at = 1:5)
mtext(text = "Coverage\n(a)", side = 1, line = 3, cex = 0.95)
text(x = 0.6, y = 1, labels = expression("VC-BART"~(rho == 0)), cex = 0.95)
text(x = 0.65, y = 2, labels = expression("lm"), cex = 0.95)
text(x = 0.65, y = 3, labels = expression("KS"), cex = 0.95)
text(x = 0.65, y = 4, labels = expression("TVCM"), cex = 0.95)
text(x = 0.8, y = 5, labels = expression("BTVCM"), cex = 0.95)
abline(v = 0.95, lty = 2, col = 'red')

boxplot.matrix(ystar_smse_test[method_names,], use.cols = FALSE, horizontal = TRUE,
               main = expression("Predictive performance"), xlab = "", yaxt = "n",
               pch = 16, cex = 0.7, medlwd = 0.5, at = 1:8)
mtext(text = "SMSE\n(c)", side = 1, line = 3, cex = 0.95)
text(x = 0.125, y = 1, labels = expression("VC-BART"~(rho == 0)), cex = 0.95)
text(x = 0.15, y = 2, labels = expression("lm"), cex = 0.95)
text(x = 0.05, y = 3, labels = expression("KS"), cex = 0.95)
text(x = 0.175, y = 4, labels = expression("TVCM"), cex = 0.95)
text(x = 0.175, y = 5, labels = expression("BTVCM"), cex = 0.95)
text(x = 0.025, y = 6, labels = expression("BART"), cex = 0.95)
text(x = 0.175, y = 7, labels = expression("ERT"), cex = 0.95)
text(x = 0.05, y = 8, labels = expression("GBM"), cex = 0.95)


boxplot.matrix(ystar_cov_test[method_names,], use.cols = FALSE, horizontal = TRUE,
               main = expression("Prediction interval coverage"), xlab = "", yaxt = "n",
               pch = 16, cex = 0.7, medlwd = 0.5, at = 1:8)
mtext(text = "Coverage\n(d)", side = 1, line = 3, cex = 0.95)
abline(v= 0.95, lty = 2, col = 'red')
text(x = 0.65, y = 1, labels = expression("VC-BART"~(rho == 0)), cex = 0.95)
text(x = 0.8, y = 2, labels = expression("lm"), cex = 0.95)
text(x = 0.2, y = 3, labels = expression("KS"), cex = 0.95)
text(x = 0.5, y = 4, labels = expression("TVCM"), cex = 0.95)
text(x = 0.5, y = 5, labels = expression("BTVCM"), cex = 0.95)
text(x = 0.7, y = 6, labels = expression("BART"), cex = 0.95)
text(x = 0.5, y = 7, labels = expression("ERT"), cex = 0.95)
text(x = 0.2, y = 8.5, labels = expression("GBM"), cex = 0.95)


dev.off()

png("~/Documents/Research/vc_bart/figures/p5R20_ystar_performance.png", 
    width = 6.5, height = 3.25, units = "in", res = 300)
par(mar = c(4,1,2,1), mgp = c(1.8, 0.5, 0), mfrow = c(1,2), cex.main = 1.1, cex.axis = 0.95, cex.lab = 1)
boxplot.matrix(ystar_smse_test[method_names,], use.cols = FALSE, horizontal = TRUE,
               main = expression("Predictive performance"), xlab = "", yaxt = "n",
               pch = 16, cex = 0.7, medlwd = 0.5, at = 1:8)
mtext(text = "SMSE\n(a)", side = 1, line = 3, cex = 0.95)
text(x = 0.15, y = 1, labels = expression("VC-BART"~(rho == 0)), cex = 0.95)
text(x = 0.175, y = 2, labels = expression("lm"), cex = 0.95)
text(x = 0.05, y = 3, labels = expression("KS"), cex = 0.95)
text(x = 0.175, y = 4, labels = expression("TVCM"), cex = 0.95)
text(x = 0.175, y = 5, labels = expression("BTVCM"), cex = 0.95)

boxplot.matrix(ystar_cov_test[method_names,], use.cols = FALSE, horizontal = TRUE,
               main = expression("Prediction interval coverage"), xlab = "", yaxt = "n",
               pch = 16, cex = 0.7, medlwd = 0.5, at = 1:8)
mtext(text = "Coverage\n(b)", side = 1, line = 3, cex = 0.95)

abline(v = 0.95, lty = 2, col = 'red')

dev.off()


boxplot.matrix(beta_int_mean, use.cols = FALSE, horizontal = TRUE,
               main = expression("Uncertainty interval length ("~beta~")"), xlab = "Coverage", yaxt = "n",
               pch = 16, cex = 0.5, medlwd = 0.5, ylim = c(0,max(beta_int_mean)))

boxplot.matrix(ystar_rmse_test^2, use.cols = FALSE, horizontal = TRUE,
               main = expression("Predictive performance"), xlab = "MSE",
               pch = 16, cex = 0.5, medlwd = 0.5)
boxplot.matrix(ystar_cov_test, use.cols = FALSE, horizontal = TRUE,
               main = expression("Uncertainly interval coverage ("~y^{"*"}~")"), xlab = "Coverage",
               pch = 16, cex = 0.5, medlwd = 0.5)
boxplot.matrix(ystar_int_test, use.cols = FALSE, horizontal = TRUE,
               main = expression("Uncertainly interval length ("~y^{"*"}~")"), xlab = "Coverage",
               pch = 16, cex = 0.5, medlwd = 0.5)
