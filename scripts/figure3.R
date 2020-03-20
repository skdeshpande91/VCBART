# Make Figure 3
# Must run example_wage1.R to get the two chains run on the entire wage 1 dataset
# Must run simulation_wage1.R and compile_wage1.R to do the cross-validation predictive simulations

load("results/wage1_chains.RData")
load("results/results_wage1.RData")

load("data/wage1.RData")

train_index <- 1:N
X_train <- as.matrix(X_all[train_index,], nrow = length(train_index), ncol = p)
Y_train <- Y_all[train_index]
Z_train <- as.matrix(Z_all[train_index,], ncol = R)
n_vec_train <- length(train_index)
start_index_train <- 1

fem_grid <- c(0,1)
mar_grid <- c(0,1)
nonwhite_grid <- c(0,1)
tenure_grid <- 0:10
exper_grid <- 1:40
Z_test <- as.matrix(expand.grid("female"= fem_grid, "nonwhite" = nonwhite_grid, "married" = mar_grid, "tenure" = tenure_grid, "exper" = exper_grid))

N_test <- nrow(Z_test)
X_test <- matrix(sample(X_train[,1], size = N_test, replace = TRUE), nrow = N_test, ncol = 1) # just randomly sample some X's
n_vec_test <- N_test
start_index_test <- 1

methods <- c("lm", "np", "tvcm", "bart", "rf", "gbm", "vc_bart")
method_labs <- c("lm", "np", "bTVCM", "BART", "ERT", "GBM", "VC-BART")
png("figures/wage1_beta_mse.png", width = 8, height = 4, units = "in", res = 600)
par(mar = c(4,3,2,1), mgp = c(1.8, 0.5, 0), cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8, mfrow = c(1,2))

plot(1, type = "n", xlim = c(0,8), ylim = c(0, max(ystar_mse_test_wage1[,methods])), 
     ylab = "MSE", xlab = "Method", xaxt = "n", main = "Predictive Performance")
boxplot(ystar_mse_test_wage1[,"lm"], at = 1, add = TRUE, pch = 16, cex = 0.5, medlwd = 0.5)
boxplot(ystar_mse_test_wage1[,"np"], at = 2, add = TRUE, pch = 16, cex = 0.5, medlwd = 0.5)
boxplot(ystar_mse_test_wage1[,"tvcm"], at = 3, add = TRUE, pch = 16, cex = 0.5, medlwd = 0.5)
boxplot(ystar_mse_test_wage1[,"bart"], at = 4, add = TRUE, pch = 16, cex = 0.5, medlwd = 0.5)
boxplot(ystar_mse_test_wage1[,"rf"], at = 5, add = TRUE, pch = 16, cex = 0.5, medlwd = 0.5)
boxplot(ystar_mse_test_wage1[,"gbm"], at = 6, add = TRUE, pch = 16, cex = 0.5, medlwd = 0.5)
boxplot(ystar_mse_test_wage1[,"vc_bart"], at  = 7, add = TRUE, pch = 16, cex = 0.5, medlwd = 0.5)
axis(side = 1, at = 1:7, labels = rep("", times = 7), tick = TRUE)
#axis(side = 1, at = 1:7, method_labs, cex = 0.775)

# female, white, married
index1_0 <- which(Z_test[,"female"] == 0 & Z_test[,"nonwhite"] == 1 & Z_test[,"married"] == 0 & Z_test[,"tenure"] == 0)
index1_1 <- which(Z_test[,"female"] == 0 & Z_test[,"nonwhite"] == 1 & Z_test[,"married"] == 0 & Z_test[,"tenure"] == 5)
index1_2 <- which(Z_test[,"female"] == 0 & Z_test[,"nonwhite"] == 1 & Z_test[,"married"] == 0 & Z_test[,"tenure"] == 10)

# male, nonwhite, single
index2_0 <- which(Z_test[,"female"] == 1 & Z_test[,"nonwhite"] == 0 & Z_test[,"married"] == 1 & Z_test[,"tenure"] == 0)
index2_1 <- which(Z_test[,"female"] == 1 & Z_test[,"nonwhite"] == 0 & Z_test[,"married"] == 1 & Z_test[,"tenure"] == 5)
index2_2 <- which(Z_test[,"female"] == 1 & Z_test[,"nonwhite"] == 0 & Z_test[,"married"] == 1 & Z_test[,"tenure"] == 10)

plot(1, type = "n", xlim = c(5, 38), ylim = c(0, 16), xlab = "Experience", ylab = "Return (%)", main = "Return to Education")
legend("topleft", legend = c( "Unmarried nonwhite male, 0 years in job", "Married white female, 5 years in job","Mincer Estimate"), col = c( "red","blue", "black"),
       pch = c(17, 16, NA), lty = c(NA, NA, 2), cex = 0.75, bty = "n")

abline(h = 100 * mincer_fit$coefficients["educ"], lty = 2)
polygon(c(Z_test[index1_1, "exper"], rev(Z_test[index1_1, "exper"])),
        100 * c(vc_sum$test$beta[index1_1, "L95",2], rev(vc_sum$test$beta[index1_1, "U95",2])),
        col = rgb(0,0,1,1/5), border = NA)
polygon(c(Z_test[index2_0, "exper"], rev(Z_test[index2_0, "exper"])),
        100 * c(vc_sum$test$beta[index2_1, "L95",2], rev(vc_sum$test$beta[index2_1, "U95",2])),
        col = rgb(1,0,0,1/5), border = NA)

lines(Z_test[index1_1, "exper"], 100 * vc_sum$test$beta[index1_1, "MEAN",2], col = 'blue') # white married female, 5 years on the job
lines(Z_test[index2_0, "exper"], 100 * vc_sum$test$beta[index2_0, "MEAN", 2], col = 'red') # nonwhite single male, 0 years of experience

points(Z_test[index1_1, "exper"], 100 * vc_sum$test$beta[index1_1, "MEAN",2], pch = 16, col = 'blue', cex = 0.5)
points(Z_test[index2_0, "exper"], 100 * vc_sum$test$beta[index2_0, "MEAN",2], pch = 17, col = 'red', cex = 0.5)


dev.off()


