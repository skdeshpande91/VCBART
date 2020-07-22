# Produces Figure S4 of the supplementary materials
load("results/sim_p5R20/results_p5R20.RData")

lin_method_names <- c("vcbart_adapt", "lm", "kernel_smoothing", "tvc", "boosted_tvcm")
method_names <- c("vcbart_adapt", "lm", "kernel_smoothing", "tvc", "boosted_tvcm", "bart", "extraTrees", "gbm")

png("figures/beta_int_p5R20.png", width = 8, height = 8*2/3, units = "in", res = 400)
par(mar = c(4.2,1,2,1), mgp  = c(1.8, 0.5, 0), mfrow = c(2,3), cex.main = 1.5, cex.axis = 1.1, cex.lab = 1.1)
boxplot(beta_int_test[lin_method_names,1,], use.cols = FALSE, horizontal = TRUE,
        main = expression("Intercept"), xlab = "", yaxt = "n",
        pch = 16, cex = 0.5, medlwd = 0.5, at = 1:5, ylim = c(0, 7))
mtext(text = "Relative lenght\n(a)", side = 1, line = 3, cex = 1.25 * par('cex'))
text(x = 3, y = 1, labels = expression("VC-BART"~(rho == 0)), cex = 0.95)
text(x = 3, y = 2, labels = expression("lm"), cex = 0.95)
text(x = 0.5, y = 3, labels = expression("KS"), cex = 0.95)
text(x = 0.5, y = 4, labels = expression("TVCM"), cex = 0.95)
text(x = 4, y = 5, labels = expression("BTVCM"), cex = 0.95)

boxplot(beta_int_test[lin_method_names,2,], use.cols = FALSE, horizontal = TRUE,
        main = expression("Effect of"~X[1]), xlab = "", yaxt = "n",
        pch = 16, cex = 0.5, medlwd = 0.5, at = 1:5, ylim = c(0,8))
mtext(text = "Relative length\n(b)", side = 1, line = 3, cex = 1.25 * par('cex'))
text(x = 3, y = 1, labels = expression("VC-BART"~(rho == 0)), cex = 0.95)
text(x = 3, y = 2, labels = expression("lm"), cex = 0.95)
text(x = 0.5, y = 3, labels = expression("KS"), cex = 0.95)
text(x = 4, y = 4, labels = expression("TVCM"), cex = 0.95)
text(x = 3, y = 5, labels = expression("BTVCM"), cex = 0.95)

boxplot(beta_int_test[lin_method_names,3,], use.cols = FALSE, horizontal = TRUE,
        main = expression("Effect of"~X[2]), xlab = "", yaxt = "n",
        pch = 16, cex = 0.5, medlwd = 0.5, at = 1:5, ylim = c(0,8))
mtext(text = "Relative length\n(c)", side = 1, line = 3, cex = 1.25 * par('cex'))
text(x =3, y = 1, labels = expression("VC-BART"~(rho == 0)), cex = 0.95)
text(x = 3, y = 2, labels = expression("lm"), cex = 0.95)
text(x = 1, y = 3, labels = expression("KS"), cex = 0.95)
text(x = 0.6, y = 4, labels = expression("TVCM"), cex = 0.95)
text(x = 3, y = 5, labels = expression("BTVCM"), cex = 0.95)

boxplot(beta_int_test[lin_method_names,4,], use.cols = FALSE, horizontal = TRUE,
        main = expression("Effect of"~X[3]), xlab = "", yaxt = "n",
        pch = 16, cex = 0.5, medlwd = 0.5, at = 1:5, ylim = c(0,12))
mtext(text = "Relative length\n(d)", side = 1, line = 3, cex = 1.25 * par('cex'))
text(x = 4, y = 1, labels = expression("VC-BART"~(rho == 0)), cex = 0.95)
text(x = 4, y = 2, labels = expression("lm"), cex = 0.95)
text(x = 1, y = 3, labels = expression("KS"), cex = 0.95)
text(x = 7, y = 4, labels = expression("TVCM"), cex = 0.95)
text(x = 4, y = 5, labels = expression("BTVCM"), cex = 0.95)

boxplot(beta_int_test[lin_method_names,5,], use.cols = FALSE, horizontal = TRUE,
        main = expression("Effect of"~X[4]), xlab = "", yaxt = "n",
        pch = 16, cex = 0.5, medlwd = 0.5, at = 1:5, ylim = c(0,6))
mtext(text = "Relative length\n(e)", side = 1, line = 3, cex = 1.25 * par('cex'))
text(x = 3, y = 1, labels = expression("VC-BART"~(rho == 0)), cex = 0.95)
text(x = 2, y = 2, labels = expression("lm"), cex = 0.95)
text(x = 0.3, y = 3, labels = expression("KS"), cex = 0.95)
text(x = 1, y = 4, labels = expression("TVCM"), cex = 0.95)
text(x = 2, y = 5, labels = expression("BTVCM"), cex = 0.95)

boxplot(beta_int_test[lin_method_names,6,], use.cols = FALSE, horizontal = TRUE,
        main = expression("Effect of"~X[5]), xlab = "", yaxt = "n",
        pch = 16, cex = 0.5, medlwd = 0.5, at = 1:5, ylim = c(0,10))
mtext(text = "Relative length\n(f)", side = 1, line = 3, cex = 1.25 * par('cex'))
text(x = 3, y = 1, labels = expression("VC-BART"~(rho == 0)), cex = 0.95)
text(x = 3, y = 2, labels = expression("lm"), cex = 0.95)
text(x = 0.6, y = 3, labels = expression("KS"), cex = 0.95)
text(x = 7, y = 4, labels = expression("TVCM"), cex = 0.95)
text(x = 3, y = 5, labels = expression("BTVCM"), cex = 0.95)

dev.off()
