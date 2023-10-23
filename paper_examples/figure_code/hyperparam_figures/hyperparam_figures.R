load("hyperparam_results.RData")

# we will use sigma2 results since
# that corresponds to sigma = 1

###################
beta_mse_mean <- apply(results_sigma2$beta_mse_test, FUN = mean, MARGIN = c(1,3))
beta_cov_mean <- apply(results_sigma2$beta_cov_test, FUN = mean, MARGIN = c(1,3))


pdf("../../figures/hyperparam_beta_mse.pdf",
    width = 6.5, height = 6.5)
par(mar = c(3,2,2,1), mgp = c(1.8, 0.5, 0))
boxplot(beta_mse_mean,
        use.cols = TRUE, horizontal = TRUE,
        main = expression("Covariate effect recovery"), xlab = "", 
        yaxt = "n",
        pch = 16, cex = 0.5, medlwd = 0.5, lwd = 0.75)
#axis(side = 2, at = c(3, 8, 13, 18), 
#     labels = paste("M =", M_list))
abline(h = c(5.5, 10.5, 15.5), col = 'gray', lty = 2)
text(x = 1.9, y = c(3,8,13,18),
     labels = paste("M =",M_list))
text(x = 2.75, y = 1:5, 
     labels = paste("tau =", tau_list))
text(x = 2.75, y = 1:5 + 5, 
     labels = paste("tau =", tau_list))
text(x = 2.75, y = 1:5 + 10, 
     labels = paste("tau =", tau_list))
text(x = 2.75, y = 1:5 + 15, 
     labels = paste("tau =", tau_list))
dev.off()

pdf("../../figures/hyperparam_ystar_rmse.pdf",
    width = 6.5, height = 6.5)
par(mar = c(3,2,2,1), mgp = c(1.8, 0.5, 0))
boxplot(results_sigma2$ystar_rmse_test,
        use.cols = TRUE, horizontal = TRUE,
        main = expression("Predictive performance"), xlab = "", 
        yaxt = "n",
        pch = 16, cex = 0.5, medlwd = 0.5, lwd = 0.75)
#axis(side = 2, at = c(3, 8, 13, 18), 
#     labels = paste("M =", M_list))
abline(h = c(5.5, 10.5, 15.5), col = 'gray', lty = 2)
text(x = c(2.7, 2.7, 2.7, 3.7), y = c(3,8,13, 18),
     labels = paste("M =",M_list))
text(x = 3.2, y = 1:5, 
     labels = paste("tau =", tau_list))
text(x = 3.2, y = 1:5 + 5, 
     labels = paste("tau =", tau_list))
text(x = 3.2, y = 1:5 + 10, 
     labels = paste("tau =", tau_list))
text(x = 4.2, y = 1:5 + 15, 
     labels = paste("tau =", tau_list))
dev.off()

pdf("../../figures/hyperparam_varsel_f1.pdf",
    width = 6.5, height = 6.5)
par(mar = c(3,2,2,1), mgp = c(1.8, 0.5, 0))
boxplot(results_sigma2$varsel[,"f1",],
        use.cols = TRUE, horizontal = TRUE,
        main = expression("F1 score"), xlab = "", 
        yaxt = "n", ylim = c(0, 1),
        pch = 16, cex = 0.5, medlwd = 0.5, lwd = 0.75)
#axis(side = 2, at = c(3, 8, 13, 18), 
#     labels = paste("M =", M_list))
abline(h = c(5.5, 10.5, 15.5), col = 'gray', lty = 2)

text(x = 0.025, y = c(3,8,13,18),
     labels = paste("M =",M_list))
text(x = 0.15, y = 1:5, 
     labels = paste("tau =", tau_list))
text(x = 0.15, y = 1:5 + 5, 
     labels = paste("tau =", tau_list))
text(x = 0.15, y = 1:5 + 10, 
     labels = paste("tau =", tau_list))
text(x = 0.15, y = 1:5 + 15, 
     labels = paste("tau =", tau_list))
dev.off()


