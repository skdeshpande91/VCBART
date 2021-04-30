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

#png("~/Documents/Research/vc_bart/figures/sensM_boxplot.png", width = 8, height = 8/3, units = "in", res = 400)
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
text(x = 3.5, y = 4, labels = "M = 100")
mtext(text = "RMSE\n(b)", side = 1, line = 3, cex = 1 * par('cex'))

boxplot(time/60, horizontal = TRUE, main = expression("Run time"), 
        yaxt = "n", xlab = "", at = 1:4, medlwd = 0.5, pch = 16, cex = 0.8)
text(x = 15, y = 1, labels = "M = 1")
text(x = 20, y = 2, labels = "M = 25")
text(x = 45, y = 3, labels = "M = 50")
text(x = 50, y = 4, labels = "M = 100")
mtext(text = "Minutes\n(c)", side = 1, line = 3, cex = 1 * par('cex'))
dev.off()

##############
beta_mse_mean <- rep(NA, times = length(M_list))
ystar_mse_mean <- rep(NA, times = length(M_list))
time_mean <- rep(NA, times = length(M_list))

for(i in 1:length(M_list)){
  M <- M_list[i]
  
  tmp_beta <- get(paste0("M", M, "_beta_mse"))
  tmp_beta_mean <- rowMeans(tmp_beta)
  
  beta_mse_mean[i] <- mean(tmp_beta_mean)
  
  ystar_mse_mean[i] <- mean(get(paste0("M", M, "_ystar_mse")))
  
  time_mean[i] <- mean(get(paste0("M", M, "_time")))
}

boxplot(beta_mse_mean[lin_method_names,], use.cols = FALSE, horizontal = TRUE,
        main = expression("Covariate effect recovery"), xlab = "", yaxt = "n",
        pch = 16, cex = 0.5, medlwd = 0.5, at = 1:5)
log_M_list <- log(M_list)

png("~/Documents/Research/vc_bart/figures/sensM_performance.png", width = 9, height = 3, units = "in", res = 600)
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0), mfrow = c(1,3), cex.main = 1.1, cex.lab = 1.1, cex.axis = 1.1)
plot(1, type = "n", xlim = range(log_M_list), ylim = c(0, max(beta_mse_mean)),
     main = "Covariate effect recovery", ylab = "MSE", xlab = "M", xaxt = "n")
axis(side = 1, at = log_M_list, labels = M_list)
lines(log_M_list, beta_mse_mean)
points(log_M_list, beta_mse_mean, pch = 16)

plot(1, type = "n", xlim = range(log_M_list), ylim = c(0, max(ystar_mse_mean)),
     main = "Predictive performance", ylab = "MSE", xlab = "M", xaxt = "n")
axis(side = 1, at = log_M_list, labels = M_list)
lines(log_M_list, ystar_mse_mean)
points(log_M_list, ystar_mse_mean, pch = 16)

log_time_list <- log(c(1, 5, 15, 30, 60, 90)) # time in minutes

plot(1, type = "n", xlim = range(log_M_list), ylim = range(log_time_list),
     main = "Run time (minutes)", ylab = "Time", xlab = "N", xaxt = "n", yaxt = "n")

points(log_M_list, log(time_mean/60), pch = 16)

#tmp_fit <- lm(log(time_mean/60) ~ log_M_list)
#abline(a = tmp_fit$coefficients[1], b = tmp_fit$coefficients[2], lty = 2)

axis(side = 1, at = log_M_list, labels = M_list)
axis(side = 2, at = log_time_list, labels = c(1, 5, 15, 30, 60, 90))
#legend("topleft", legend = paste0("slope = ", round(tmp_fit$coefficients[2], digits = 2)), lty = 2, cex = 1.1)

dev.off()

