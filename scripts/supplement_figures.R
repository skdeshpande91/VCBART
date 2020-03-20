# Figures S1 and S2
load("results/tau_p3R2.RData")
load("results/tau_p5R5.RData")

new_order <- c("tau_1_4", "tau_1_2", "tau_1_1", "tau_2_1", "tau_4_1")

beta_avg_mse_p3R2 <- apply(beta_mse_test_tau_p3R2, FUN = mean, MARGIN = c(1,3), na.rm = TRUE)[new_order,] # computes average MSE across betas for p3R2
beta_avg_mse_p3R2[is.nan(beta_avg_mse_p3R2)] <- NA # needed in case one of the simulation runs errored out

beta_avg_mse_p5R5 <- apply(beta_mse_test_tau_p5R5, FUN = mean, MARGIN = c(1,3), na.rm = TRUE)[new_order,] # computes average MSE across betas for p5R5
beta_avg_mse_p5R5[is.nan(beta_avg_mse_p5R5)] <- NA

beta_avg_coverage_p3R2 <- apply(beta_coverage_test_tau_p3R2, FUN = mean, MARGIN = c(1,3), na.rm = TRUE)[new_order,] # computes the average coverage across betas for p3R2
beta_avg_coverage_p3R2[is.nan(beta_avg_coverage_p3R2)] <- NA

beta_avg_coverage_p5R5 <- apply(beta_coverage_test_tau_p5R5, FUN = mean, MARGIN = c(1,3), na.rm = TRUE)[new_order,] # computes the average coverage across betas for p5R5
beta_avg_coverage_p5R5[is.nan(beta_avg_coverage_p5R5)] <- NA

beta_avg_int_p3R2 <- apply(beta_int_test_tau_p3R2, FUN = mean, MARGIN = c(1,3), na.rm = TRUE)[new_order,] # computes the relative interval length across betas for p3R2
beta_avg_int_p3R2[is.nan(beta_avg_int_p3R2)] <- NA

beta_avg_int_p5R5 <- apply(beta_int_test_tau_p5R5, FUN = mean, MARGIN = c(1,3), na.rm = TRUE)[new_order,] # computes the relative interval length across betas for p5R5
beta_avg_int_p5R5[is.nan(beta_avg_int_p5R5)] <- NA


ystar_mse_p3R2 <- ystar_mse_test_tau_p3R2[,new_order]
ystar_mse_p5R5 <- ystar_mse_test_tau_p5R5[,new_order]

ystar_coverage_p3R2 <- ystar_coverage_test_tau_p3R2[,new_order]
ystar_coverage_p5R5 <- ystar_coverage_test_tau_p5R5[,new_order]

ystar_int_p3R2 <- ystar_int_test_tau_p3R2[,new_order]
ystar_int_p5R5 <- ystar_int_test_tau_p5R5[,new_order]

tau_seq <- c(1/4, 1/2, 1, 2, 4)
tau_labs <- c("1/4", "1/2", "1", "2", "4")

png("figures/hyperparameter_sensitivity_p3R2.png", width = 6, height = 4, units = "in", res = 600)
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0), mfrow = c(2,3), cex.main = 0.9, cex.lab = 0.9, cex.axis = 0.9)
#tau_seq <- c(1/4, 1/2, 2/3, 1, 3/2, 2, 4)
#tau_labs <- c("1/4", "1/2", "2/3", "1", "3/2", "2", "4")


# beta mse for p3R2
plot(1, type = "n", xlim = c(0, 6), ylim = range(beta_avg_mse_p3R2), xlab = expression(tau), xaxt = "n", 
     ylab = expression("MSE"), main = expression("MSE"~(beta)))
for(tau_ix in 1:length(tau_seq)){
  boxplot(beta_avg_mse_p3R2[tau_ix,], at = tau_ix, add = TRUE, pch = 16, cex = 0.7, medlwd = 0.5)
}
axis(side = 1, at = 1:length(tau_seq), labels = tau_labs)

# beta coverage for p3R2
plot(1, type = "n", xlim = c(0, 6), ylim = c(0.5, 1), xlab = expression(tau), xaxt = "n", 
     ylab = expression("Coverage"), main = expression("Coverage"~(beta)))
for(tau_ix in 1:length(tau_seq)){
  boxplot(beta_avg_coverage_p3R2[tau_ix,], at = tau_ix, add = TRUE, pch = 16, cex = 0.7, medlwd = 0.5)
}
#abline(h = 0.95, lty = 2)
axis(side = 1, at = 1:length(tau_seq), labels = tau_labs)

# beta interval length for p3R2
plot(1, type = "n", xlim = c(0, 6), ylim = c(0.8, 1.2), xlab = expression(tau), xaxt = "n", 
     ylab = expression("Relative Interval Length"), main = expression("Relative Length"~(beta)))
for(tau_ix in 1:length(tau_seq)){
  boxplot(beta_avg_int_p3R2[tau_ix,], at = tau_ix, add = TRUE, pch = 16, cex = 0.7, medlwd = 0.5)
}
abline(h = 1, lty = 2)
axis(side = 1, at = 1:length(tau_seq), labels = tau_labs)

# ystar mse for p3R2
plot(1, type = "n", xlim = c(0, 6), ylim = range(ystar_mse_p3R2), xlab = expression(tau), xaxt = "n", 
     ylab = expression("MSE"), main = expression("MSE"~(y^{"*"})))
for(tau_ix in 1:length(tau_seq)){
  boxplot(ystar_mse_p3R2[,tau_ix], at = tau_ix, add = TRUE, pch = 16, cex = 0.7, medlwd = 0.5)
}
axis(side = 1, at = 1:length(tau_seq), labels = tau_labs)

# ystar coverage for p3R2
plot(1, type = "n", xlim = c(0, 6), ylim = c(0.9, 1), xlab = expression(tau), xaxt = "n", 
     ylab = expression("Coverage"), main = expression("Coverage"~(y^{"*"})))
for(tau_ix in 1:length(tau_seq)){
  boxplot(ystar_coverage_p3R2[,tau_ix], at = tau_ix, add = TRUE, pch = 16, cex = 0.7, medlwd = 0.5)
}
#abline(h = 0.95, lty = 2)
axis(side = 1, at = 1:length(tau_seq), labels = tau_labs)

# ystar interval length for p3R2
plot(1, type = "n", xlim = c(0, 6), ylim = c(0.9, 1.1), xlab = expression(tau), xaxt = "n", 
     ylab = expression("Relative Interval Length"), main = expression("Relative Length"~(y^{"*"})))
for(tau_ix in 1:length(tau_seq)){
  boxplot(ystar_int_p3R2[,tau_ix], at = tau_ix, add = TRUE, pch = 16, cex = 0.7, medlwd = 0.5)
}
#abline(h = 1, lty = 2)
axis(side = 1, at = 1:length(tau_seq), labels = tau_labs)

dev.off()

png("figures/hyperparameter_sensitivity_p5R5.png", width = 6, height = 4, units = "in", res = 600)
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0), mfrow = c(2,3), cex.main = 0.9, cex.lab = 0.9, cex.axis = 0.9)
#tau_seq <- c(1/4, 1/2, 2/3, 1, 3/2, 2, 4)
#tau_labs <- c("1/4", "1/2", "2/3", "1", "3/2", "2", "4")


# beta mse for p5R5
plot(1, type = "n", xlim = c(0, 6), ylim = range(beta_avg_mse_p5R5), xlab = expression(tau), xaxt = "n", 
     ylab = expression("MSE"), main = expression("MSE"~(beta)))
for(tau_ix in 1:length(tau_seq)){
  boxplot(beta_avg_mse_p5R5[tau_ix,], at = tau_ix, add = TRUE, pch = 16, cex = 0.7, medlwd = 0.5)
}
axis(side = 1, at = 1:length(tau_seq), labels = tau_labs)

# beta coverage for p5R5
plot(1, type = "n", xlim = c(0, 6), ylim = c(0.75, 1), xlab = expression(tau), xaxt = "n", 
     ylab = expression("Coverage"), main = expression("Coverage"~(beta)))
for(tau_ix in 1:length(tau_seq)){
  boxplot(beta_avg_coverage_p5R5[tau_ix,], at = tau_ix, add = TRUE, pch = 16, cex = 0.7, medlwd = 0.5)
}
#abline(h = 0.95, lty = 2)
axis(side = 1, at = 1:length(tau_seq), labels = tau_labs)

# beta interval length for p5R5
plot(1, type = "n", xlim = c(0, 6), ylim = c(0.5, 1.5), xlab = expression(tau), xaxt = "n", 
     ylab = expression("Relative Interval Length"), main = expression("Relative Length"~(beta)))
for(tau_ix in 1:length(tau_seq)){
  boxplot(beta_avg_int_p5R5[tau_ix,], at = tau_ix, add = TRUE, pch = 16, cex = 0.7, medlwd = 0.5)
}
#abline(h = 1, lty = 2)
axis(side = 1, at = 1:length(tau_seq), labels = tau_labs)

# ystar mse for p5R5
plot(1, type = "n", xlim = c(0, 6), ylim = range(ystar_mse_p5R5), xlab = expression(tau), xaxt = "n", 
     ylab = expression("MSE"), main = expression("MSE"~(y^{"*"})))
for(tau_ix in 1:length(tau_seq)){
  boxplot(ystar_mse_p5R5[,tau_ix], at = tau_ix, add = TRUE, pch = 16, cex = 0.7, medlwd = 0.5)
}
axis(side = 1, at = 1:length(tau_seq), labels = tau_labs)

# ystar coverage for p5R5
plot(1, type = "n", xlim = c(0, 6), ylim = c(0.9, 1), xlab = expression(tau), xaxt = "n", 
     ylab = expression("Coverage"), main = expression("Coverage"~(y^{"*"})))
for(tau_ix in 1:length(tau_seq)){
  boxplot(ystar_coverage_p5R5[,tau_ix], at = tau_ix, add = TRUE, pch = 16, cex = 0.7, medlwd = 0.5)
}
#abline(h = 0.95, lty = 2)
axis(side = 1, at = 1:length(tau_seq), labels = tau_labs)

# ystar interval length for p5R5
plot(1, type = "n", xlim = c(0, 6), ylim = c(0.9, 1.1), xlab = expression(tau), xaxt = "n", 
     ylab = expression("Relative Interval Length"), main = expression("Relative Length"~(y^{"*"})))
for(tau_ix in 1:length(tau_seq)){
  boxplot(ystar_int_p5R5[,tau_ix], at = tau_ix, add = TRUE, pch = 16, cex = 0.7, medlwd = 0.5)
}
#abline(h = 1, lty = 2)
axis(side = 1, at = 1:length(tau_seq), labels = tau_labs)

dev.off()