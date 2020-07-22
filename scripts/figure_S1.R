# Script to reproduce Figure S1 in the Supplementary Materials.
load("results/sim_p5R20_tau/results_p5R20_tau.RData")

beta_mse_mean <- apply(beta_mse_test, FUN = mean, MARGIN = c(1,3))
beta_cov_mean <- apply(beta_cov_test, FUN = mean, MARGIN = c(1,3))
beta_int_mean <- apply(beta_int_test, FUN = mean, MARGIN = c(1,3))

new_names <- paste0("tau_", c("1_4", "1_2", "2_3", "1_1", "3_2", "2_1", "4_1"))

png("figures/tau_sensitivity_p5R20.png", width = 8, height = 8 * 2/3, units = "in", res = 400)

par(mar = c(4.2,1,2,1), mgp  = c(1.8, 0.5, 0), mfrow = c(2,3), cex.main = 1.5, cex.axis = 1.1, cex.lab = 1.1)
boxplot.matrix(beta_mse_mean[new_names,], use.cols = FALSE, horizontal = TRUE,
               main = expression("Covariate effect recovery"), xlab = "", yaxt = "n",
               pch = 16, cex = 0.5, medlwd = 0.5, at = 1:7, ylim = c(0, 2.5))
mtext(text = "MSE\n(a)", side = 1, line = 3, cex = 1.25 * par('cex'))

text(x = 0.2, y = 1, labels = expression(tau==0.25), cex = 0.95)
text(x = 0.2, y = 2, labels = expression(tau==0.5), cex = 0.95)
text(x = 0.2, y = 3, labels = expression(tau==0.67), cex = 0.95)
text(x = 0.2, y = 4, labels = expression(tau==1), cex = 0.95)
text(x = 0.2, y = 5, labels = expression(tau==1.5), cex = 0.95)
text(x = 0.2, y = 6, labels = expression(tau==2), cex = 0.95)
text(x = 0.2, y = 7, labels = expression(tau==4), cex = 0.95)

boxplot.matrix(beta_cov_mean[new_names,], use.cols = FALSE, horizontal = TRUE,
               main = expression("Uncertainty interval coverage"), xlab = "", yaxt = "n",
               pch = 16, cex = 0.5, medlwd = 0.5, at = 1:7, ylim = c(0.5, 1))
mtext(text = "Coverage\n(b)", side = 1, line = 3, cex = 1.25 * par('cex'))

text(x = 0.55, y = 1, labels = expression(tau==0.25), cex = 0.95)
text(x = 0.55, y = 2, labels = expression(tau==0.5), cex = 0.95)
text(x = 0.55, y = 3, labels = expression(tau==0.67), cex = 0.95)
text(x = 0.55, y = 4, labels = expression(tau==1), cex = 0.95)
text(x = 0.55, y = 5, labels = expression(tau==1.5), cex = 0.95)
text(x = 0.55, y = 6, labels = expression(tau==2), cex = 0.95)
text(x = 0.55, y = 7, labels = expression(tau==4), cex = 0.95)
abline(v = 0.95, lty = 2, col = 'red')

boxplot.matrix(beta_int_mean[new_names,], use.cols = FALSE, horizontal = TRUE,
               main = expression("Uncertainty interval length"), xlab = "", yaxt = "n",
               pch = 16, cex = 0.5, medlwd = 0.5, at = 1:7, ylim = c(0.5, 2))
mtext(text = "Relative Length \n(c)", side = 1, line = 3, cex = 1.25 * par('cex'))

text(x = 0.65, y = 1, labels = expression(tau==0.25), cex = 0.95)
text(x = 0.65, y = 2, labels = expression(tau==0.5), cex = 0.95)
text(x = 0.65, y = 3, labels = expression(tau==0.67), cex = 0.95)
text(x = 0.65, y = 4, labels = expression(tau==1), cex = 0.95)
text(x = 0.65, y = 5, labels = expression(tau==1.5), cex = 0.95)
text(x = 0.65, y = 6, labels = expression(tau==2), cex = 0.95)
text(x = 0.65, y = 7, labels = expression(tau==4), cex = 0.95)

boxplot.matrix(ystar_rmse_test[new_names,], use.cols = FALSE, horizontal = TRUE,
               main = expression("Predictive performance"), xlab = "", yaxt = "n",
               pch = 16, cex = 0.7, medlwd = 0.5, at = 1:7, ylim = c(1, 3.5))
mtext(text = "RMSE \n(d)", side = 1, line = 3, cex = 1.25 * par('cex'))
text(x = 1.2, y = 1, labels = expression(tau==0.25), cex = 0.95)
text(x = 1.2, y = 2, labels = expression(tau==0.5), cex = 0.95)
text(x = 1.2, y = 3, labels = expression(tau==0.67), cex = 0.95)
text(x = 1.2, y = 4, labels = expression(tau==1), cex = 0.95)
text(x = 1.2, y = 5, labels = expression(tau==1.5), cex = 0.95)
text(x = 1.2, y = 6, labels = expression(tau==2), cex = 0.95)
text(x = 1.2, y = 7, labels = expression(tau==4), cex = 0.95)

boxplot.matrix(ystar_cov_test[new_names,], use.cols = FALSE, horizontal = TRUE,
               main = expression("Predictive interval coverage"), xlab = "", yaxt = "n",
               pch = 16, cex = 0.7, medlwd = 0.5, at = 1:7, ylim = c(0.8, 1))
mtext(text = "Coverage \n(e)", side = 1, line = 3, cex = 1.25 * par('cex'))
text(x = 0.825, y = 1, labels = expression(tau==0.25), cex = 0.95)
text(x = 0.825, y = 2, labels = expression(tau==0.5), cex = 0.95)
text(x = 0.825, y = 3, labels = expression(tau==0.67), cex = 0.95)
text(x = 0.825, y = 4, labels = expression(tau==1), cex = 0.95)
text(x = 0.825, y = 5, labels = expression(tau==1.5), cex = 0.95)
text(x = 0.825, y = 6, labels = expression(tau==2), cex = 0.95)
text(x = 0.825, y = 7, labels = expression(tau==4), cex = 0.95)
abline(v = 0.95, lty = 2, col = 'red')


boxplot.matrix(ystar_int_test[new_names,], use.cols = FALSE, horizontal = TRUE,
               main = expression("Predictive interval length"), xlab = "", yaxt = "n",
               pch = 16, cex = 0.7, medlwd = 0.5, at = 1:7, ylim = c(0.5, 1.5))
mtext(text = "Relative Length \n(f)", side = 1, line = 3, cex = 1.25 * par('cex'))
text(x = 0.65, y = 1, labels = expression(tau==0.25), cex = 0.95)
text(x = 0.65, y = 2, labels = expression(tau==0.5), cex = 0.95)
text(x = 0.65, y = 3, labels = expression(tau==0.67), cex = 0.95)
text(x = 0.65, y = 4, labels = expression(tau==1), cex = 0.95)
text(x = 0.65, y = 5, labels = expression(tau==1.5), cex = 0.95)
text(x = 0.65, y = 6, labels = expression(tau==2), cex = 0.95)
text(x = 0.65, y = 7, labels = expression(tau==4), cex = 0.95)

dev.off()