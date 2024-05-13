####################################################################
# Define a (mostly) color-blind friendly palette
####################################################################

my_colors <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
               "#8F2727")
my_rgb <- col2rgb(my_colors, alpha = FALSE)/255

####################################################################
# Load the compiled results from running the MCMC for longer
####################################################################
load("longer_run_results.RData")

nd_list <- c(500, 1000, 2500, 5000, 10000, 25000, 50000)

beta_mse_mean <- apply(beta_mse_test, FUN = mean, MARGIN = c(2,3))
beta_cov_mean <- apply(beta_cov_test, FUN = mean, MARGIN = c(2,3))


pdf(file = "../../figures/mixing_rhat.pdf", width = 6, height = 6*9/16)
par(mar = c(3,1,2,1), mgp = c(1.8, 0.5, 0))
boxplot(sigma_rhat, ylim = c(0.875, 1.5), main = "Sigma R-hat",
        pch = 16, cex = 0.75, medlwd = 0.5, horizontal = TRUE,
        names = NA,
        yaxt = "n",
        xlab = "R-hat")
text(x = 0.925, y = 1:7, labels = paste(2*nd_list, "iterations"), cex = 0.8)
abline(v = 1.1, col = my_colors[4])
dev.off()

pdf(file = "../../figures/mixing_beta_mse.pdf", width = 6, height = 6*9/16)
par(mar = c(3,1,2,1), mgp = c(1.8, 0.5, 0))
boxplot(beta_mse_mean, ylim = c(0, 0.5), main = "Beta MSE",
        pch = 16, cex = 0.75, medlwd = 0.5, horizontal = TRUE,
        names = NA,
        yaxt = "n",
        xlab = "Mean square error")
text(x = 0.25, y = 1:7, labels = paste(2*nd_list, "iterations"), cex = 0.8)
dev.off()

pdf(file = "../../figures/mixing_ystar_rmse.pdf", width = 6, height = 6*9/16)
par(mar = c(3,1,2,1), mgp = c(1.8, 0.5, 0))
boxplot(ystar_rmse_test, ylim = c(1, 1.75), main = "Y* RMSE",
        pch = 16, cex = 0.75, medlwd = 0.5, horizontal = TRUE,
        names = NA,
        yaxt = "n",
        xlab = "Root mean square error")
text(x = 1.05, y = 1:7, labels = paste(2*nd_list, "iterations"), cex = 0.8)
dev.off()


pdf(file = "../../figures/mixing_beta_cov.pdf", width = 6, height = 6*9/16)
par(mar = c(3,1,2,1), mgp = c(1.8, 0.5, 0))
boxplot(beta_cov_mean, ylim = c(0.95, 1), main = "Beta Coverage",
        pch = 16, cex = 0.75, medlwd = 0.5, horizontal = TRUE,
        names = NA,
        yaxt = "n",
        xlab = "Coverage")
text(x = 0.96, y = 1:7, labels = paste(2*nd_list, "iterations"), cex = 0.8)
dev.off()

pdf(file = "../../figures/mixing_ystar_cov.pdf", width = 6, height = 6*9/16)
par(mar = c(3,1,2,1), mgp = c(1.8, 0.5, 0))
boxplot(ystar_cov_test, ylim = c(0.92, 1), main = "Y* Coverage",
        pch = 16, cex = 0.75, medlwd = 0.5, horizontal = TRUE,
        names = NA,
        yaxt = "n",
        xlab = "Coverage")
text(x = 0.94, y = 1:7, labels = paste(2*nd_list, "iterations"), cex = 0.8)
dev.off()

pdf(file = "../../figures/mixing_f1.pdf", width = 6, height = 6*9/16)
par(mar = c(3,1,2,1), mgp = c(1.8, 0.5, 0))
boxplot(varsel[,"f1",], ylim = c(0.75, 1), main = "Modifier selection (F1)",
        pch = 16, cex = 0.75, medlwd = 0.5, horizontal = TRUE,
        names = NA,
        yaxt = "n",
        xlab = "F1 Score")
text(x = 0.8, y = 1:7, labels = paste(2*nd_list, "iterations"), cex = 0.8)
dev.off()
