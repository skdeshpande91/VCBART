########
# Script to produce Figure S7 of the Supplemental materials

M_list <- c(1, 25, 50, 100)
load("results/sensM_results.RData")

beta_mse <- data.frame("M1" = rowMeans(M1_beta_mse),
                       "M25" = rowMeans(M25_beta_mse),
                       "M50" = rowMeans(M50_beta_mse),
                       "M100" = rowMeans(M100_beta_mse))
ystar_rmse <- data.frame("M1" = sqrt(M1_ystar_mse),
                         "M25" = sqrt(M25_ystar_mse),
                         "M50" = sqrt(M50_ystar_mse),
                         "M100" = sqrt(M100_ystar_mse))
time <- data.frame("M1" = M1_time, 
                   "M25" = M25_time,
                   "M50" = M50_time,
                   "M100" = M100_time)

f1 <- data.frame("M1" = M1_varsel[,"f1"],
                 "M25" = M25_varsel[,"f1"],
                 "M50" = M50_varsel[,"f1"],
                 "M100" = M100_varsel[,"f1"])

#png("~/Documents/Research/vc_bart/figures/sensM_boxplot.png", width = 8, height = 8/3, units = "in", res = 400)
#png("figures/sensM_boxplot.png", width = 8, height = 8/3, units = "in", res = 400)
#par(mar = c(4.2,1,2,1), mgp  = c(1.8, 0.5, 0), mfrow = c(1,3), cex.main = 1, cex.axis = 1, cex.lab = 1)
png("figures/sensM_boxplot.png", width = 8, height = 8/3, units = "in", res = 400)
par(mar = c(4.2,1,2,1), mgp  = c(1.8, 0.5, 0), mfrow = c(1,3), cex.main = 1, cex.axis = 1, cex.lab = 1)
boxplot(beta_mse, horizontal = TRUE, main = expression("Covariate effect recovery"), 
        yaxt = "n", xlab = "", at = 1:4, medlwd = 0.5, pch = 16, cex = 0.8)
text(x = 1, y = 1, labels = "M = 1")
text(x = 3, y = 2, labels = "M = 25")
text(x = 3, y = 3, labels = "M = 50")
text(x = 3, y = 4, labels = "M = 100")
mtext(text = "MSE\n(a)", side = 1, line = 3, cex = 1 * par('cex'))

boxplot(ystar_rmse, horizontal = TRUE, main = expression("Predictive performance"), 
        yaxt = "n", xlab = "", medlwd = 0.5, pch = 16, cex = 0.8)
text(x = 2, y = 1, labels = "M = 1")
text(x = 3.5, y = 2, labels = "M = 25")
text(x = 3.5, y = 3, labels = "M = 50")
text(x = 4, y = 4, labels = "M = 100")
mtext(text = "RMSE\n(b)", side = 1, line = 3, cex = 1 * par('cex'))

boxplot(f1, horizontal = TRUE, main = expression("F1 score"),
        yaxt = "n", xlab = "", at = 1:4, medlwd = 0.5, pch = 16, cex = 0.8, ylim = c(0,1))
text(x = 0.3, y = 1, labels = "M = 1")
text(x = 0.5, y = 2, labels = "M = 25")
text(x = 0.5, y = 3, labels = "M = 50")
text(x = 0.45, y = 4, labels = "M = 100")
mtext(text = "F1\n(c)", side = 1, line = 3, cex = 1 * par('cex'))

dev.off()


