# Tabulate the results from the p3R2 simulation and make Table 1
load("results/results_p3R2.RData")

# Compute the average of the metrics across the beta's
beta_avg_mse <- apply(beta_mse_test_p3R2, FUN = mean, MARGIN = c(1,3), na.rm = TRUE) # computes average MSE across beta
beta_avg_mse[is.nan(beta_avg_mse)] <- NA

beta_avg_mse_mean <- round(apply(beta_avg_mse, FUN = mean, MARGIN = 1, na.rm = TRUE), digits = 2)
beta_avg_mse_sd <- round(apply(beta_avg_mse, FUN = sd, MARGIN = 1, na.rm = TRUE), digits = 2)


beta_avg_coverage <- apply(beta_coverage_test_p3R2, FUN = mean, MARGIN = c(1,3), na.rm = TRUE) # computes the avg MSE 
beta_avg_coverage[is.nan(beta_avg_coverage)] <- NA

beta_avg_coverage_mean <- round(apply(beta_avg_coverage, FUN = mean, MARGIN = 1, na.rm = TRUE), digits = 2)
beta_avg_coverage_sd <- round(apply(beta_avg_coverage, FUN = sd, MARGIN = 1, na.rm = TRUE), digits = 2)

beta_avg_int <- apply(beta_int_test_p3R2, FUN = mean, MARGIN = c(1,3), na.rm = TRUE) # computes the avg MSE 
beta_avg_int[is.nan(beta_avg_int)] <- NA

beta_avg_int_mean <- round(apply(beta_avg_int, FUN = mean, MARGIN = 1, na.rm = TRUE), digits = 2)
beta_avg_int_sd <- round(apply(beta_avg_int, FUN = sd, MARGIN = 1, na.rm = TRUE), digits = 2)


beta_avg_table <- data.frame("MSE" = paste0(beta_avg_mse_mean, " (", beta_avg_mse_sd, ")"),
                             "COV" = paste0(beta_avg_coverage_mean, " (", beta_avg_coverage_sd, ")"),
                             "INT" = paste0(beta_avg_int_mean, " (", beta_avg_int_sd, ")"))
rownames(beta_avg_table) <- dimnames(beta_mse_test_p3R2)[[1]]


# Average predictive performance
new_order <- c("lm", "np", "tvcm","bart", "rf", "gbm", "vc_bart")
ystar_mse_mean <- round(apply(ystar_mse_test_p3R2, FUN = mean, MARGIN = 2, na.rm = TRUE), digits = 2)[new_order]
ystar_mse_sd <- round(apply(ystar_mse_test_p3R2, FUN = sd, MARGIN = 2, na.rm = TRUE), digits = 2)[new_order]

ystar_coverage_mean <- round(apply(ystar_coverage_test_p3R2, FUN = mean, MARGIN = 2, na.rm = TRUE), digits = 2)[new_order]
ystar_coverage_sd <- round(apply(ystar_coverage_test_p3R2, FUN = sd, MARGIN = 2, na.rm = TRUE), digits = 2)[new_order]

ystar_int_mean <- round(apply(ystar_int_test_p3R2, FUN = mean, MARGIN = 2, na.rm = TRUE), digits = 2)[new_order]
ystar_int_sd <- round(apply(ystar_int_test_p3R2, FUN = sd, MARGIN = 2, na.rm = TRUE), digits = 2)[new_order]

timing_mean <- round(apply(timing_p3R2, FUN = mean, MARGIN = 2, na.rm = TRUE), digits = 2)[new_order]
timing_sd <- round(apply(timing_p3R2, FUN = sd, MARGIN = 2, na.rm = TRUE), digits = 2)[new_order]

ystar_table <- data.frame("MSE" = paste0(ystar_mse_mean, " (", ystar_mse_sd, ")"),
                          "COV" = paste0(ystar_coverage_mean, " (", ystar_coverage_sd, ")"),
                          "INT" = paste0(ystar_int_mean, " (", ystar_int_sd, ")"),
                          "TIME" = paste0(timing_mean, " (", timing_sd,")"))
rownames(ystar_table) <- new_order


# Assess performance function-by-function

beta_mse_mean <- round(apply(beta_mse_test_p3R2, MARGIN = c(1,2), FUN = mean, na.rm = TRUE), digits = 2)
beta_mse_sd <- round(apply(beta_mse_test_p3R2, MARGIN = c(1,2), FUN = sd, na.rm = TRUE), digits = 2)

beta_mse_table <- data.frame("beta0" = paste0(beta_mse_mean[,1], " (", beta_mse_sd[,1], ")"),
                             "beta1" = paste0(beta_mse_mean[,2], " (", beta_mse_sd[,2], ")"),
                             "beta2" = paste0(beta_mse_mean[,3], " (", beta_mse_sd[,3], ")"),
                             "beta3" = paste0(beta_mse_mean[,4], " (", beta_mse_sd[,4], ")"))
rownames(beta_mse_table) <- dimnames(beta_mse_test_p3R2)[[1]]

beta_coverage_mean <- round(apply(beta_coverage_test_p3R2, MARGIN = c(1,2), FUN = mean, na.rm = TRUE), digits = 2)
beta_coverage_sd <- round(apply(beta_coverage_test_p3R2, MARGIN = c(1,2), FUN = sd, na.rm = TRUE), digits = 2)

beta_coverage_table <- data.frame("beta0" = paste0(beta_coverage_mean[,1], " (", beta_coverage_sd[,1], ")"),
                                  "beta1" = paste0(beta_coverage_mean[,2], " (", beta_coverage_sd[,2], ")"),
                                  "beta2" = paste0(beta_coverage_mean[,3], " (", beta_coverage_sd[,3], ")"),
                                  "beta3" = paste0(beta_coverage_mean[,4], " (", beta_coverage_sd[,4], ")"))
rownames(beta_coverage_table) <- dimnames(beta_mse_test_p3R2)[[1]]

beta_int_mean <- round(apply(beta_int_test_p3R2, MARGIN = c(1,2), FUN = mean, na.rm = TRUE), digits = 2)
beta_int_sd <- round(apply(beta_int_test_p3R2, MARGIN = c(1,2), FUN = sd, na.rm = TRUE), digits = 2)

beta_int_table <- data.frame("beta0" = paste0(beta_int_mean[,1], " (", beta_int_sd[,1], ")"),
                             "beta1" = paste0(beta_int_mean[,2], " (", beta_int_sd[,2], ")"),
                             "beta2" = paste0(beta_int_mean[,3], " (", beta_int_sd[,3], ")"),
                             "beta3" = paste0(beta_int_mean[,4], " (", beta_int_sd[,4], ")"))
rownames(beta_int_table) <- dimnames(beta_mse_test_p3R2)[[1]]


beta_mse_table

beta_coverage_table


beta_int_table
