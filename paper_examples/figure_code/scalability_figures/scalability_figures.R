n_train_list <- c(25, 50, 100, 
                  125,250, 500, 1000,
                  1250, 2500, 5000, 10000, 12500)

ystar_rmse_test <- rep(NA, times = length(n_train_list))
ystar_cov_test <- rep(NA, times = length(n_train_list))

beta_mse_test <- rep(NA, times = length(n_train_list))
beta_cov_test <- rep(NA, times = length(n_train_list))

timing <- rep(NA, times = length(n_train_list))

sen <- rep(NA, times = length(n_train_list))
spec <- rep(NA, times = length(n_train_list))
prec <- rep(NA, times = length(n_train_list))
f1 <- rep(NA, times = length(n_train_list))
fp <- rep(NA, times = length(n_train_list))
fn <- rep(NA, times = length(n_train_list))

###########################
# Collate the results
###########################

for(nix in 1:length(n_train_list)){
  n_train <- n_train_list[nix]
  if(file.exists(paste0("p5R20_demo_n", n_train, ".RData"))){
    load(paste0("p5R20_demo_n", n_train, ".RData"))
    tmp <- get(paste0("results_n", n_train))
    
    ystar_rmse_test[nix] <- mean(tmp$ystar_rmse_test)
    beta_mse_test[nix] <- mean(rowMeans(tmp$beta_mse_test))
    
    #ystar_cov_test[nix] <- mean(tmp$ystar_cov_test)
    
    #beta_cov_test[nix] <- mean(rowMeans(tmp$beta_cov_test))
    
    timing[nix] <- mean(tmp$timing)
    
    sen[nix] <- mean(tmp$varsel[,"sen"])
    spec[nix] <- mean(tmp$varsel[,"spec"])
    prec[nix] <- mean(tmp$varsel[,"prec"])
    f1[nix] <- mean(tmp$varsel[,"f1"])
    fp[nix] <- mean(tmp$varsel[,"fp"])
    fn[nix] <- mean(tmp$varsel[,"fn"])
    
    rm(list = paste0("results_n", n_train))
  }
}

###################################
# Total observations on log-scale
###################################

N_list <- 4 * n_train_list
log_N_list <- log(N_list, base = 10)
log_N_lim <- c(1, 5)

xlabels <- c(50, 500, 5000, 50000)
ticks <- log(xlabels, base = 10)


###################################
# Covariate effect recovery
###################################

pdf("../../figures/p5R20_scale_beta.pdf", width = 3, height = 3)
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0))
plot(1, type = "n", 
     xlim = log_N_lim, ylim = c(0,5),
     xaxt = "n", xlab = expression("Total observations (N)"),
     ylab = expression("MSE"), main = expression("Covariate effect recovery"))
axis(side = 1, at = ticks, labels = xlabels)
lines(log_N_list, beta_mse_test)
points(log_N_list, beta_mse_test, pch = 16, cex = 0.7)
dev.off()

#################################
# Predictive RMSE
#################################
pdf("../../figures/p5R20_scale_ystar.pdf", width = 3, height = 3)
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0))
plot(1, type = "n", 
     xlim = log_N_lim, ylim = c(0,5),
     xaxt = "n", xlab = expression("Total observations (N)"),
     ylab = expression("RMSE"), main = expression("Predictive performance"))
axis(side = 1, at = ticks, labels = xlabels)
lines(log_N_list, ystar_rmse_test)
points(log_N_list, ystar_rmse_test, pch = 16, cex = 0.7)
dev.off()

###############################
# Variable selection
###############################
my_colors <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


pdf("../../figures/p5R20_scale_varsel.pdf", width = 3, height = 3)
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0))
plot(1, type = "n", 
     xlim = log_N_lim, ylim = c(0,1),
     xaxt = "n", xlab = expression("Total observations (N)"),
     ylab = "", main = expression("Variable selection performance"))
axis(side = 1, at = ticks, labels = xlabels)

lines(log_N_list, sen, col = my_colors[2])
points(log_N_list, sen, col = my_colors[2], pch = 15)

lines(log_N_list, spec, col = my_colors[3])
points(log_N_list, spec, col = my_colors[3], pch = 16)

lines(log_N_list, prec, col = my_colors[4])
points(log_N_list, prec, col = my_colors[4], pch = 17)

lines(log_N_list, f1, col = my_colors[8])
points(log_N_list, f1, col = my_colors[8], pch = 18)

legend("bottomright", bty = "n", pch = 15:18, 
       col = my_colors[c(2,3,4,8)],
       legend = c("Sensitivity", "Specificity", "Precision", "F1"))
dev.off()


###################
# Runtime
###################

pdf("../../figures/p5R20_scale_time.pdf", width = 3, height = 3)
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0))
plot(1, type = "n", 
     xlim = log_N_lim, ylim = c(0, 4800),
     xaxt = "n", xlab = expression("Total observations (N)"),
     ylab = expression("Seconds"), main = expression("Runtime"))
axis(side = 1, at = ticks, labels = xlabels)

lines(log_N_list, timing)
points(log_N_list, timing, pch = 16)

dev.off()
