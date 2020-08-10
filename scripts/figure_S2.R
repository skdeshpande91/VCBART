# Produces figure s2 of the suppelementary materials
load("results/sim_p5R20/results_p5R20.RData")

lin_method_names <- c("vcbart_adapt", "lm", "kernel_smoothing", "tvc", "boosted_tvcm")
method_names <- c("vcbart_adapt", "lm", "kernel_smoothing", "tvc", "boosted_tvcm", "bart", "extraTrees", "gbm")

#png("~/Documents/Research/vc_bart/figures/beta_mse_p5R20.png", width = 8, height = 8*2/3, units = "in", res = 400)
png("figures/beta_mse_p5R20.png", width = 8, height = 8*2/3, units = "in", res = 400)

par(mar = c(4.2,1,2,1), mgp  = c(1.8, 0.5, 0), mfrow = c(2,3), cex.main = 1.5, cex.axis = 1.1, cex.lab = 1.1)
boxplot(beta_mse_test[lin_method_names,1,], use.cols = FALSE, horizontal = TRUE,
        main = expression("Intercept"), xlab = "", yaxt = "n",
        pch = 16, cex = 0.5, medlwd = 0.5, at = 1:5)
mtext(text = "MSE\n(a)", side = 1, line = 3, cex = 1.25 * par('cex'))
text(x = 3, y = 1, labels = expression("VCBART"~(rho == 0)), cex = 0.95)
text(x = 6.5, y = 2, labels = expression("lm"), cex = 0.95)
text(x = 7.5, y = 3, labels = expression("KS"), cex = 0.95)
text(x = 4, y = 4, labels = expression("TVCM"), cex = 0.95)
text(x = 4, y = 5, labels = expression("BTVCM"), cex = 0.95)

boxplot(beta_mse_test[lin_method_names,2,], use.cols = FALSE, horizontal = TRUE,
        main = expression("Effect of"~X[1]), xlab = "", yaxt = "n",
        pch = 16, cex = 0.5, medlwd = 0.5, at = 1:5)
mtext(text = "MSE\n(b)", side = 1, line = 3, cex = 1.25 * par('cex'))
text(x = 6.5, y = 1, labels = expression("VCBART"~(rho == 0)), cex = 0.95)
text(x = 6.5, y = 2, labels = expression("lm"), cex = 0.95)
text(x = 2, y = 3, labels = expression("KS"), cex = 0.95)
text(x = 6.5, y = 4, labels = expression("TVCM"), cex = 0.95)
text(x = 6.5, y = 5, labels = expression("BTVCM"), cex = 0.95)

boxplot(beta_mse_test[lin_method_names,3,], use.cols = FALSE, horizontal = TRUE,
        main = expression("Effect of"~X[2]), xlab = "", yaxt = "n",
        pch = 16, cex = 0.5, medlwd = 0.5, at = 1:5)
mtext(text = "MSE\n(c)", side = 1, line = 3, cex = 1.25 * par('cex'))
text(x = 10, y = 1, labels = expression("VCBART"~(rho == 0)), cex = 0.95)
text(x = 2.5, y = 2, labels = expression("lm"), cex = 0.95)
text(x = 2, y = 3, labels = expression("KS"), cex = 0.95)
text(x = 10, y = 4, labels = expression("TVCM"), cex = 0.95)
text(x = 10, y = 5, labels = expression("BTVCM"), cex = 0.95)

boxplot(beta_mse_test[lin_method_names,4,], use.cols = FALSE, horizontal = TRUE,
        main = expression("Effect of"~X[3]), xlab = "", yaxt = "n",
        pch = 16, cex = 0.5, medlwd = 0.5, at = 1:5)
mtext(text = "MSE\n(d)", side = 1, line = 3, cex = 1.25 * par('cex'))
text(x = 6.5, y = 1, labels = expression("VC-BART"~(rho == 0)), cex = 0.95)
text(x = 6.5, y = 2, labels = expression("lm"), cex = 0.95)
text(x = 1, y = 3, labels = expression("KS"), cex = 0.95)
text(x = 10, y = 4, labels = expression("TVCM"), cex = 0.95)
text(x = 7.5, y = 5, labels = expression("BTVCM"), cex = 0.95)

boxplot(beta_mse_test[lin_method_names,5,], use.cols = FALSE, horizontal = TRUE,
        main = expression("Effect of"~X[4]), xlab = "", yaxt = "n",
        pch = 16, cex = 0.5, medlwd = 0.5, at = 1:5)
mtext(text = "MSE\n(e)", side = 1, line = 3, cex = 1.25 * par('cex'))
text(x = 15, y = 1, labels = expression("VCBART"~(rho == 0)), cex = 0.95)
text(x = 20, y = 2, labels = expression("lm"), cex = 0.95)
text(x = 3, y = 3, labels = expression("KS"), cex = 0.95)
text(x = 3, y = 4, labels = expression("TVCM"), cex = 0.95)
text(x = 5, y = 5, labels = expression("BTVCM"), cex = 0.95)

boxplot(beta_mse_test[lin_method_names,6,], use.cols = FALSE, horizontal = TRUE,
        main = expression("Effect of"~X[5]), xlab = "", yaxt = "n",
        pch = 16, cex = 0.5, medlwd = 0.5, at = 1:5)
mtext(text = "MSE\n(f)", side = 1, line = 3, cex = 1.25 * par('cex'))
text(x = 7.5, y = 1, labels = expression("VCBART"~(rho == 0)), cex = 0.95)
text(x = 5, y = 2, labels = expression("lm"), cex = 0.95)
text(x = 1, y = 3, labels = expression("KS"), cex = 0.95)
text(x = 12, y = 4, labels = expression("TVCM"), cex = 0.95)
text(x = 5, y = 5, labels = expression("BTVCM"), cex = 0.95)

dev.off()
