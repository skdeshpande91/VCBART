n_list <- c(100, 500, 1000, 5000, 10000, 50000)
load("results/bigN_results.RData")

beta_mse_mean <- rep(NA, times = length(n_list))
beta_mse_sd <- rep(NA, times = length(n_list))
ystar_rmse_mean <- rep(NA, times = length(n_list))
ystar_rmse_sd <- rep(NA, times = length(n_list))
time_mean <- rep(NA, times = length(n_list))
time_sd <- rep(NA, times = length(n_list))

for(i in 1:length(n_list)){
  n <- n_list[i]
  
  tmp_beta <- get(paste0("bigN_n", n, "_beta_mse"))
  tmp_beta_mean <- rowMeans(tmp_beta)
  
  beta_mse_mean[i] <- mean(tmp_beta_mean)
  beta_mse_sd[i] <- sd(tmp_beta_mean)
  
  ystar_rmse_mean[i] <- mean(sqrt(get(paste0("bigN_n", n, "_ystar_mse"))))
  ystar_rmse_sd[i] <- sd(sqrt(get(paste0("bigN_n", n, "_ystar_mse"))))
  
  time_mean[i] <- mean(get(paste0("bigN_n", n, "_time")))
  time_sd[i] <- sd(get(paste0("bigN_n", n, "_time")))
  
}

log_n_list <- log(n_list)

png("figures/bigN_performance.png", width = 8, height = 8*1/3, units = "in", res = 400)
par(mar = c(4.2,3,2,1), mgp = c(1.8, 0.5, 0), mfrow = c(1,3), cex.main = 1.2, cex.lab = 1.1, cex.axis = 1.1)
plot(1, type = "n", xlim = range(log_n_list), ylim = range(beta_mse_mean, na.rm = TRUE),
     main = "Covariate effect recovery", ylab = "MSE", xlab = "", xaxt = "n")
axis(side = 1, at = log_n_list, labels = n_list)
mtext(text = "N (log-scale)\n(a)", side = 1, line = 3, cex = 1.2 * par('cex'))
lines(log_n_list, beta_mse_mean)
points(log_n_list, beta_mse_mean, pch = 16)

plot(1, type = "n", xlim = range(log_n_list), ylim = range(ystar_rmse_mean, na.rm = TRUE),
     main = "Predictive performance", ylab = "RMSE", xlab = "", xaxt = "n")
axis(side = 1, at = log_n_list, labels = n_list)
mtext(text = "N (log-scale)\n(b)", side = 1, line = 3, cex = 1.2 * par('cex'))
lines(log_n_list, ystar_rmse_mean)
points(log_n_list, ystar_rmse_mean, pch = 16)

log_time_list <- log(c(1, 5, 15, 30, 60, 240, 720)) # time in minutes

plot(1, type = "n", xlim = range(log_n_list), ylim = range(log_time_list),
     main = "Run time (minutes)", ylab = "Minutes (log-scale)", xlab = "", xaxt = "n", yaxt = "n")

points(log_n_list, log(time_mean/60), pch = 16)
tmp_fit <- lm(log(time_mean/60) ~ log_n_list)
abline(a = tmp_fit$coefficients[1], b = tmp_fit$coefficients[2], lty = 2)

axis(side = 1, at = log_n_list, labels = n_list)
axis(side = 2, at = log_time_list, labels = c(1, 5, 15, 30, 60, 240, 720))
mtext(text = "N (log-scale)\n(c)", side = 1, line = 3, cex = 1.2 * par('cex'))

legend("topleft", legend = paste0("slope = ", round(tmp_fit$coefficients[2], digits = 2)), lty = 2, cex = 1.1)

dev.off()

