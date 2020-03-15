library(VCBART)
load("data/wage1.RData")


N_train <- nrow(X_all)
X_train <- cbind(rep(1, times = N_train), X_all)
Y_train <- Y_all
Z_train <- Z_all
######
# Create a grid of testing points
# Really interested in the modifiers
#######

fem_grid <- c(0,1)
mar_grid <- c(0,1)
nonwhite_grid <- c(0,1)
tenure_grid <- 0:10
exper_grid <- 1:40
Z_test <- as.matrix(expand.grid("female"= fem_grid, "nonwhite" = nonwhite_grid, "married" = mar_grid, "tenure" = tenure_grid, "exper" = exper_grid))

N_test <- nrow(Z_test)
X_test <- X_train[sample(1:N_train, size = N_test,replace = TRUE), ]


chain1 <- vc_BART_ind(Y = Y_train, X_train = X_train, Z_train = Z_train,
                      X_test = X_test, Z_test = Z_test, xinfo_list = cutpoints, 
                      nd = 1000, burn = 250, verbose = TRUE, print_every = 50)

chain2 <- vc_BART_ind(Y = Y_train, X_train = X_train, Z_train = Z_train,
                      X_test = X_test, Z_test = Z_test, xinfo_list = cutpoints, 
                      nd = 1000, burn = 250, verbose = TRUE, print_every = 50)

fit_sum <- get_summary(chain1, chain2)

######
# Fit to Mincer's original specification
#####

lm_data <- data.frame("Y" = Y_all, cov_all, mod_all)
lm_data[,"exper2"] <- lm_data[,"exper"] * lm_data[,"exper"]

mincer_fit <- lm(Y ~ educ + exper + exper2, data = lm_data)
marginal_fit <- lm(Y ~ educ, data = lm_data)


# female, white, married
index1 <- which(Z_test[,"female"] == 0 & Z_test[,"nonwhite"] == 1 & Z_test[,"married"] == 0 & Z_test[,"tenure"] == 5)

# male, nonwhite, single
index2 <- which(Z_test[,"female"] == 1 & Z_test[,"nonwhite"] == 0 & Z_test[,"married"] == 1 & Z_test[,"tenure"] == 0)


png("figures/wage1_beta.png", width = 6, height = 4, units = "in", res = 600)
par(mar = c(3,3,2,1), mgp = c(1.5, 0.8, 0), cex.main = 0.9, cex.lab = 0.85, cex.axis = 0.85)
plot(1, type = "n", xlim = c(5, 38), ylim = c(0, 16), xlab = "Potential Years of Work Experience", ylab = "Return (%)", main = "Return to Education")
legend("topleft", legend = c("Married white female, 5 years with current employer", "Unmarried nonwhite male, 0 years with employer", "Marginal return"), col = c("blue", "red", "black"),
       pch = c(16, 17, NA), lty = c(NA, NA, 2), cex = 0.7, bty = "n")

abline(h = 100 * marginal_fit$coefficients["educ"], lty = 2)
polygon(c(Z_test[index1, "exper"], rev(Z_test[index1, "exper"])),
        100 * c(fit_sum$test$beta[index1, "L95",2], rev(fit_sum$test$beta[index1, "U95",2])),
        col = rgb(0,0,1,1/5), border = NA)
polygon(c(Z_test[index2, "exper"], rev(Z_test[index2, "exper"])),
        100 * c(fit_sum$test$beta[index2, "L95",2], rev(fit_sum$test$beta[index2, "U95",2])),
        col = rgb(1,0,0,1/5), border = NA)

lines(Z_test[index1, "exper"], 100 * fit_sum$test$beta[index1, "MEAN",2], col = 'blue') # white married female, 5 years on the job
lines(Z_test[index2, "exper"], 100 * fit_sum$test$beta[index2, "MEAN", 2], col = 'red') # nonwhite single male, 0 years of experience

points(Z_test[index1, "exper"], 100 * fit_sum$test$beta[index1, "MEAN",2], pch = 16, col = 'blue', cex = 0.7)
points(Z_test[index2, "exper"], 100 * fit_sum$test$beta[index2, "MEAN",2], pch = 17, col = 'red', cex = 0.7)


dev.off()
