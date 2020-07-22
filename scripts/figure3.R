# Make some figures about the predictive accuracy of VC-BART
load("results/results_HRS.RData")

method_names <- c("vcbart_adapt", "vcbart_adapt_cs25", "vcbart_adapt_cs50", "vcbart_adapt_cs75",
                  "lm", "boosted_tvcm", "bart", "extraTrees", "gbm")

#png("figures/hrs_performance.png", width = 6, height = 3, units = "in", res = 300)
#png("~/Documents/Research/vc_bart/figures/hrs_performance.png", width = 8, height = 8*1/3, units = "in", res =400)
png("figures/hrs_performance.png", width = 8, height = 8*1/3, units = "in", res =400)

par(mar = c(4,1,2,1), mgp = c(1.8, 0.5, 0), mfrow = c(1,3), 
    cex.main = 1.5, cex.axis = 1, cex.lab = 1)
boxplot.matrix(ystar_rmse_test[method_names,], use.cols = FALSE, yaxt ="n",
               horizontal = TRUE, ylim = c(3,4.5), pch = 16, cex = 0.5, medlwd = 0.75,
               main = "Predictive performance", xlab = "", at = 1:9)
mtext(text = "RMSE\n(a)", side = 1, line = 3, cex = 1.25 * par("cex"))
# add custom labels
text(x = 3.3, y = 1, labels = expression("VC-BART"~(rho == 0)), cex = 0.95)
text(x = 3.3, y = 2, labels = expression("VC-BART"~(rho == 0.25)), cex = 0.95)
text(x = 3.3, y = 3, labels = expression("VC-BART"~(rho == 0.5)), cex = 0.95)
text(x = 3.3, y = 4, labels = expression("VC-BART"~(rho == 0.75)), cex = 0.95)
text(x = 3.5, y = 5, labels = expression("lm"), cex = 0.95)
text(x = 3.5, y = 6, labels = expression("BTVCM"), cex = 0.95)
text(x = 3.5, y = 7, labels = expression("BART"), cex = 0.95)
text(x = 3.5, y = 8, labels = expression("ERT"), cex = 0.95)
text(x = 3.5, y = 9, labels = expression("GBM"), cex = 0.95)

boxplot.matrix(ystar_cov_test[method_names,], use.cols = FALSE, yaxt ="n",
               horizontal = TRUE, ylim = c(0.7,1), pch = 16, cex = 0.5, medlwd = 0.75,
               main = "Predictive interval coverage", xlab = "")
mtext(text = "Coverage\n(b)", side = 1, line = 3, cex = 1.25 * par("cex"))
text(x = 0.775, y = 1, labels = expression("VC-BART"~(rho == 0)), cex = 0.95)
text(x = 0.775, y = 2, labels = expression("VC-BART"~(rho == 0.25)), cex = 0.95)
text(x = 0.775, y = 3, labels = expression("VC-BART"~(rho == 0.5)), cex = 0.95)
text(x = 0.775, y = 4, labels = expression("VC-BART"~(rho == 0.75)), cex = 0.95)
text(x = 0.75, y = 5, labels = expression("lm"), cex = 0.95)
text(x = 0.75, y = 6, labels = expression("BTVCM"), cex = 0.95)
text(x = 0.75, y = 7, labels = expression("BART"), cex = 0.95)
text(x = 0.75, y = 8, labels = expression("ERT"), cex = 0.95)
text(x = 0.75, y = 9, labels = expression("GBM"), cex = 0.95)
abline(v = 0.95, col = 'red', lty = 2)

boxplot.matrix(ystar_int_test[method_names,], use.cols = FALSE, yaxt ="n",
               horizontal = TRUE, pch = 16, cex = 0.5, medlwd = 0.75,
               main = "Prediction interval length", xlab = "", ylim = c(0.5, 1.5), at = 1:9)
mtext(text = "Relative length \n(c)", side = 1, line = 3, cex = 1.25 * par("cex"))
text(x = 1.3, y = 1, labels = expression("VC-BART"~(rho == 0)), cex = 0.95)
text(x = 1.3, y = 2, labels = expression("VC-BART"~(rho == 0.25)), cex = 0.95)
text(x = 1.3, y = 3, labels = expression("VC-BART"~(rho == 0.5)), cex = 0.95)
text(x = 1.1, y = 4, labels = expression("VC-BART"~(rho == 0.75)), cex = 0.95)
text(x = 1.3, y = 5, labels = expression("lm"), cex = 0.95)
text(x = 0.9, y = 6, labels = expression("BTVCM"), cex = 0.95)
text(x = 1.1, y = 7, labels = expression("BART"), cex = 0.95)
text(x = 0.85, y = 8, labels = expression("ERT"), cex = 0.95)
text(x = 1, y = 9, labels = expression("GBM"), cex = 0.95)

dev.off()

png("~/Documents/Research/vc_bart/figures/hrs_performance_alt.png", 
    width = 6.5, height = 6.5*1/2, units = "in", res =400)
par(mar = c(4,1,2,1), mgp = c(1.8, 0.5, 0), mfrow = c(1,2), 
    cex.main = 1, cex.axis = 1, cex.lab = 1)
boxplot.matrix(ystar_rmse_test[method_names,], use.cols = FALSE, yaxt ="n",
               horizontal = TRUE, ylim = c(3,4.5), pch = 16, cex = 0.5, medlwd = 0.75,
               main = "Predictive performance", xlab = "", at = 1:9)
mtext(text = "RMSE\n(a)", side = 1, line = 3, cex = 1 * par("cex"))
# add custom labels
text(x = 3.3, y = 1, labels = expression("VC-BART"~(rho == 0)), cex = 0.8)
text(x = 3.3, y = 2, labels = expression("VC-BART"~(rho == 0.25)), cex = 0.8)
text(x = 3.3, y = 3, labels = expression("VC-BART"~(rho == 0.5)), cex = 0.8)
text(x = 3.3, y = 4, labels = expression("VC-BART"~(rho == 0.75)), cex = 0.8)
text(x = 3.5, y = 5, labels = expression("lm"), cex = 0.8)
text(x = 3.5, y = 6, labels = expression("BTVCM"), cex = 0.8)
text(x = 3.5, y = 7, labels = expression("BART"), cex = 0.8)
text(x = 3.5, y = 8, labels = expression("ERT"), cex = 0.8)
text(x = 3.5, y = 9, labels = expression("GBM"), cex = 0.8)

boxplot.matrix(ystar_cov_test[method_names,], use.cols = FALSE, yaxt ="n",
               horizontal = TRUE, ylim = c(0.7,1), pch = 16, cex = 0.5, medlwd = 0.75,
               main = "Predictive interval coverage", xlab = "")
mtext(text = "Coverage\n(b)", side = 1, line = 3, cex = 1 * par("cex"))
text(x = 0.775, y = 1, labels = expression("VC-BART"~(rho == 0)), cex = 0.8)
text(x = 0.775, y = 2, labels = expression("VC-BART"~(rho == 0.25)), cex = 0.8)
text(x = 0.775, y = 3, labels = expression("VC-BART"~(rho == 0.5)), cex = 0.8)
text(x = 0.775, y = 4, labels = expression("VC-BART"~(rho == 0.75)), cex = 0.8)
text(x = 0.75, y = 5, labels = expression("lm"), cex = 0.8)
text(x = 0.75, y = 6, labels = expression("BTVCM"), cex = 0.8)
text(x = 0.75, y = 7, labels = expression("BART"), cex = 0.8)
text(x = 0.75, y = 8, labels = expression("ERT"), cex = 0.8)
text(x = 0.75, y = 9, labels = expression("GBM"), cex = 0.8)
abline(v = 0.95, col = 'red', lty = 2)

dev.off()
