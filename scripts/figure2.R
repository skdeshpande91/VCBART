# Make Figure 2
sim_type <- "ar75"
rho_true <- 0.75

if(file.exists(paste0("results/panel_p3R2_", sim_type, ".RData"))){
  load(paste0("results/panel_p3R2_", sim_type, ".RData"))
  
  beta_mse_all <- apply(get(paste0("beta_mse_test_p3R2_", sim_type)), FUN = mean, MARGIN = c(1,3), na.rm = TRUE)
  
  beta_mse_avg <- apply(beta_mse_all, FUN = mean, MARGIN = 1, na.rm = TRUE) # average MSE across all functions and simulation settings
  beta_mse_sd <- apply(beta_mse_all, FUN = sd, MARGIN = 1, na.rm = TRUE) # standard deviation across simulations
  beta_mse_upper <- beta_mse_avg + beta_mse_sd
  beta_mse_lower <- beta_mse_avg - beta_mse_sd
  
  
  
  beta_coverage_all <- apply(get(paste0("beta_coverage_test_p3R2_", sim_type)), FUN = mean, MARGIN = c(1,3), na.rm = TRUE)
  
  beta_coverage_avg <- apply(beta_coverage_all, FUN = mean, MARGIN = 1, na.rm = TRUE)
  beta_coverage_sd <- apply(beta_coverage_all, FUN = sd, MARGIN = 1, na.rm = TRUE)
  beta_coverage_upper <- beta_coverage_avg + beta_coverage_sd
  beta_coverage_lower <- beta_coverage_avg - beta_coverage_sd
  
  
  beta_int_all <- apply(get(paste0("beta_int_test_p3R2_", sim_type)), FUN = mean, MARGIN = c(1,3), na.rm = TRUE)
  
  beta_int_avg <- apply(beta_int_all, FUN = mean, MARGIN = 1, na.rm = TRUE)
  beta_int_sd <- apply(beta_int_all, FUN = sd, MARGIN = 1, na.rm = TRUE)
  beta_int_upper <- beta_int_avg + beta_int_sd
  beta_int_lower <- beta_int_avg - beta_int_sd
  
  ystar_mse_avg <- apply(get(paste0("ystar_mse_test_p3R2_", sim_type)), MARGIN = 2, FUN = mean, na.rm = TRUE)
  ystar_mse_sd <- apply(get(paste0("ystar_mse_test_p3R2_", sim_type)), MARGIN = 2, FUN = sd, na.rm = TRUE)
  ystar_mse_upper <- ystar_mse_avg + ystar_mse_sd
  ystar_mse_lower <- ystar_mse_avg - ystar_mse_sd
  
  
  
  ystar_coverage_avg <- apply(get(paste0("ystar_coverage_test_p3R2_", sim_type)), MARGIN = 2, FUN = mean, na.rm = TRUE)
  ystar_coverage_sd <- apply(get(paste0("ystar_coverage_test_p3R2_", sim_type)), MARGIN = 2, FUN = sd, na.rm = TRUE)
  ystar_coverage_upper <- ystar_coverage_avg + ystar_coverage_sd
  ystar_coverage_lower <- ystar_coverage_avg - ystar_coverage_sd
  
  ystar_int_avg <- apply(get(paste0("ystar_int_test_p3R2_", sim_type)), MARGIN = 2, FUN = mean, na.rm = TRUE)
  ystar_int_sd <- apply(get(paste0("ystar_int_test_p3R2_", sim_type)), MARGIN = 2, FUN = sd, na.rm = TRUE)
  ystar_int_upper <- ystar_int_avg + ystar_int_sd
  ystar_int_lower <- ystar_int_avg - ystar_int_sd
  
  
  
  ######
  # Get rid of some files
  ######
  rm_list <- paste0(c("beta_mse_train_", "beta_coverage_train_", "beta_int_train_",
                      "ystar_mse_train_", "ystar_coverage_train_", "ystar_int_train_",
                      "beta_mse_test_", "beta_coverage_test_", "beta_int_test_",
                      "ystar_mse_test_", "ystar_coverage_test_", "ystar_int_test_",
                      "sigma_est_", "timing_"), "p3R2_", sim_type)
  rm(list = rm_list)
  
  #######################
  # Now make the plot
  #######################
  png(file= paste0("figures/panel_p3R2_", sim_type, "_performance.png"), width = 6.5, height = 6.5/3 * 2, units = "in", res = 600)
  par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0.0), mfrow = c(2,3), cex.main = 0.95, cex.lab = 0.95, cex.axis = 0.9)
  
  # BETA MSE
  beta_mse_ylim <- c(0, max(beta_mse_upper) * 1.05)
  plot(1, type = "n", xlim = c(0,0.9), ylim = beta_mse_ylim, xlab = expression(rho), ylab = "MSE", main = expression("Average MSE"~(beta)))
  polygon( c(seq(0, 0.9, by = 0.1), rev(seq(0, 0.9, by = 0.1))), 
           c(beta_mse_upper[c("ind", paste0("ar_rho",1:9))], rev(beta_mse_lower[c("ind", paste0("ar_rho", 1:9))])),
           col = rgb(0,0,1,1/5), border = NA)
  lines(seq(0, 0.9, by = 0.1), beta_mse_avg[c("ind", paste0("ar_rho",1:9))], col = 'blue')        
  points(seq(0, 0.9, by = 0.1), beta_mse_avg[c("ind", paste0("ar_rho",1:9))], col = 'blue', pch = 16, cex = 0.8)
  
  polygon(c(seq(0.0, 0.9, by = 0.1), rev(seq(0., 0.9, by = 0.1))),
          c(beta_mse_upper[c("ind", paste0("cs_rho",1:9))], rev(beta_mse_lower[c("ind", paste0("cs_rho", 1:9))])),
          col = rgb(1,0,0,1/5), border = NA)
  lines(seq(0.0, 0.9, by = 0.1), beta_mse_avg[c("ind", paste0("cs_rho",1:9))], col = 'red')        
  points(seq(0.0, 0.9, by = 0.1), beta_mse_avg[c("ind",paste0("cs_rho",1:9))], col = 'red', pch = 17, cex = 0.8)
  
  legend("bottomleft", horiz = TRUE, bty = "n", 
         legend = c("AR", "CS", "Truth"), pch = c(16, 17, 8), col = c("blue", "red", "black"), cex = 0.9)
  #lines(x = c(rho_true, rho_true), y = c(beta_mse_upper["true"], beta_mse_lower["true"]))
  points(x = rho_true, y = beta_mse_avg["true"], pch = 8)
  
  
  # BETA COVERAGE
  beta_coverage_ylim <- c(0.95 * min(beta_coverage_lower), 1)
  plot(1, type = "n", xlim = c(0,0.9), ylim = beta_coverage_ylim, xlab = expression(rho), ylab = "Coverage", main = expression("Interval Coverage"~(beta)))
  polygon( c(seq(0, 0.9, by = 0.1), rev(seq(0, 0.9, by = 0.1))), 
           c(beta_coverage_upper[c("ind", paste0("ar_rho",1:9))], rev(beta_coverage_lower[c("ind", paste0("ar_rho", 1:9))])),
           col = rgb(0,0,1,1/5), border = NA)
  lines(seq(0, 0.9, by = 0.1), beta_coverage_avg[c("ind", paste0("ar_rho",1:9))], col = 'blue')        
  points(seq(0, 0.9, by = 0.1), beta_coverage_avg[c("ind", paste0("ar_rho",1:9))], col = 'blue', pch = 16, cex = 0.8)
  
  polygon(c(seq(0.0, 0.9, by = 0.1), rev(seq(0., 0.9, by = 0.1))),
          c(beta_coverage_upper[c("ind", paste0("cs_rho",1:9))], rev(beta_coverage_lower[c("ind", paste0("cs_rho", 1:9))])),
          col = rgb(1,0,0,1/5), border = NA)
  lines(seq(0.0, 0.9, by = 0.1), beta_coverage_avg[c("ind", paste0("cs_rho",1:9))], col = 'red')        
  points(seq(0.0, 0.9, by = 0.1), beta_coverage_avg[c("ind",paste0("cs_rho",1:9))], col = 'red', pch = 17, cex = 0.8)
  
  legend("bottomleft", horiz = TRUE, bty = "n", 
         legend = c("AR", "CS", "Truth"), pch = c(16, 17, 8), col = c("blue", "red", "black"), cex = 0.9)
  #lines(x = c(rho_true, rho_true), y = c(beta_coverage_upper["true"], beta_coverage_lower["true"]))
  points(x = rho_true, y = beta_coverage_avg["true"], pch = 8)
  
  # BETA INTERVAL
  beta_int_ylim <- c(0.95 * min(beta_int_lower), 1.05 * max(beta_int_upper))
  plot(1, type = "n", xlim = c(0,0.9), ylim = beta_int_ylim, xlab = expression(rho), ylab = "Relative Length", main = expression("Relative Interval Length"~(beta)))
  polygon( c(seq(0, 0.9, by = 0.1), rev(seq(0, 0.9, by = 0.1))), 
           c(beta_int_upper[c("ind", paste0("ar_rho",1:9))], rev(beta_int_lower[c("ind", paste0("ar_rho", 1:9))])),
           col = rgb(0,0,1,1/5), border = NA)
  lines(seq(0, 0.9, by = 0.1), beta_int_avg[c("ind", paste0("ar_rho",1:9))], col = 'blue')        
  points(seq(0, 0.9, by = 0.1), beta_int_avg[c("ind", paste0("ar_rho",1:9))], col = 'blue', pch = 16, cex = 0.8)
  
  polygon(c(seq(0.0, 0.9, by = 0.1), rev(seq(0., 0.9, by = 0.1))),
          c(beta_int_upper[c("ind", paste0("cs_rho",1:9))], rev(beta_int_lower[c("ind", paste0("cs_rho", 1:9))])),
          col = rgb(1,0,0,1/5), border = NA)
  lines(seq(0.0, 0.9, by = 0.1), beta_int_avg[c("ind", paste0("cs_rho",1:9))], col = 'red')        
  points(seq(0.0, 0.9, by = 0.1), beta_int_avg[c("ind",paste0("cs_rho",1:9))], col = 'red', pch = 17, cex = 0.8)
  
  legend("topleft", horiz = TRUE, bty = "n", 
         legend = c("AR", "CS", "Truth"), pch = c(16, 17, 8), col = c("blue", "red", "black"), cex = 0.9)
  #lines(x = c(rho_true, rho_true), y = c(beta_int_upper["true"], beta_int_lower["true"]))
  points(x = rho_true, y = beta_int_avg["true"], pch = 8)
  
  
  # ystar mse
  ystar_mse_ylim <- c(0, max(ystar_mse_upper) * 1.05)
  plot(1, type = "n", xlim = c(0,0.9), ylim = ystar_mse_ylim, xlab = expression(rho), ylab = "MSE", main = expression("MSE"~(y~"*")))
  polygon( c(seq(0, 0.9, by = 0.1), rev(seq(0, 0.9, by = 0.1))), 
           c(ystar_mse_upper[c("ind", paste0("ar_rho",1:9))], rev(ystar_mse_lower[c("ind", paste0("ar_rho", 1:9))])),
           col = rgb(0,0,1,1/5), border = NA)
  lines(seq(0, 0.9, by = 0.1), ystar_mse_avg[c("ind", paste0("ar_rho",1:9))], col = 'blue')        
  points(seq(0, 0.9, by = 0.1), ystar_mse_avg[c("ind", paste0("ar_rho",1:9))], col = 'blue', pch = 16, cex = 0.8)
  
  polygon(c(seq(0.0, 0.9, by = 0.1), rev(seq(0., 0.9, by = 0.1))),
          c(ystar_mse_upper[c("ind", paste0("cs_rho",1:9))], rev(ystar_mse_lower[c("ind", paste0("cs_rho", 1:9))])),
          col = rgb(1,0,0,1/5), border = NA)
  lines(seq(0.0, 0.9, by = 0.1), ystar_mse_avg[c("ind", paste0("cs_rho",1:9))], col = 'red')        
  points(seq(0.0, 0.9, by = 0.1), ystar_mse_avg[c("ind",paste0("cs_rho",1:9))], col = 'red', pch = 17, cex = 0.8)
  
  legend("bottomleft", horiz = TRUE, bty = "n", 
         legend = c("AR", "CS", "Truth"), pch = c(16, 17, 8), col = c("blue", "red", "black"), cex = 0.9)
  #lines(x = c(rho_true, rho_true), y = c(ystar_mse_upper["true"], ystar_mse_lower["true"]))
  points(x = rho_true, y = ystar_mse_avg["true"], pch = 8)
  
  
  # ystar COVERAGE
  ystar_coverage_ylim <- c(0.95 * min(ystar_coverage_lower), 1.0)
  plot(1, type = "n", xlim = c(0,0.9), ylim = ystar_coverage_ylim, xlab = expression(rho), ylab = "Coverage", main = expression("Interval Coverage"~(y~"*")))
  polygon( c(seq(0, 0.9, by = 0.1), rev(seq(0, 0.9, by = 0.1))), 
           c(ystar_coverage_upper[c("ind", paste0("ar_rho",1:9))], rev(ystar_coverage_lower[c("ind", paste0("ar_rho", 1:9))])),
           col = rgb(0,0,1,1/5), border = NA)
  lines(seq(0, 0.9, by = 0.1), ystar_coverage_avg[c("ind", paste0("ar_rho",1:9))], col = 'blue')        
  points(seq(0, 0.9, by = 0.1), ystar_coverage_avg[c("ind", paste0("ar_rho",1:9))], col = 'blue', pch = 16, cex = 0.8)
  
  polygon(c(seq(0.0, 0.9, by = 0.1), rev(seq(0., 0.9, by = 0.1))),
          c(ystar_coverage_upper[c("ind", paste0("cs_rho",1:9))], rev(ystar_coverage_lower[c("ind", paste0("cs_rho", 1:9))])),
          col = rgb(1,0,0,1/5), border = NA)
  lines(seq(0.0, 0.9, by = 0.1), ystar_coverage_avg[c("ind", paste0("cs_rho",1:9))], col = 'red')        
  points(seq(0.0, 0.9, by = 0.1), ystar_coverage_avg[c("ind",paste0("cs_rho",1:9))], col = 'red', pch = 17, cex = 0.8)
  
  legend("bottomleft", horiz = TRUE, bty = "n", 
         legend = c("AR", "CS", "Truth"), pch = c(16, 17, 8), col = c("blue", "red", "black"), cex = 0.9)
  #lines(x = c(rho_true, rho_true), y = c(ystar_coverage_upper["true"], ystar_coverage_lower["true"]))
  points(x = rho_true, y = ystar_coverage_avg["true"], pch = 8)
  
  # ystar INTERVAL
  ystar_int_ylim <- c(0.95 * min(ystar_int_lower), 1.05 * max(ystar_int_upper))
  plot(1, type = "n", xlim = c(0,0.9), ylim = ystar_int_ylim, xlab = expression(rho), ylab = "Relative Length", main = expression("Relative Interval Length"~(y~"*")))
  polygon( c(seq(0, 0.9, by = 0.1), rev(seq(0, 0.9, by = 0.1))), 
           c(ystar_int_upper[c("ind", paste0("ar_rho",1:9))], rev(ystar_int_lower[c("ind", paste0("ar_rho", 1:9))])),
           col = rgb(0,0,1,1/5), border = NA)
  lines(seq(0, 0.9, by = 0.1), ystar_int_avg[c("ind", paste0("ar_rho",1:9))], col = 'blue')        
  points(seq(0, 0.9, by = 0.1), ystar_int_avg[c("ind", paste0("ar_rho",1:9))], col = 'blue', pch = 16, cex = 0.8)
  
  polygon(c(seq(0.0, 0.9, by = 0.1), rev(seq(0., 0.9, by = 0.1))),
          c(ystar_int_upper[c("ind", paste0("cs_rho",1:9))], rev(ystar_int_lower[c("ind", paste0("cs_rho", 1:9))])),
          col = rgb(1,0,0,1/5), border = NA)
  lines(seq(0.0, 0.9, by = 0.1), ystar_int_avg[c("ind", paste0("cs_rho",1:9))], col = 'red')        
  points(seq(0.0, 0.9, by = 0.1), ystar_int_avg[c("ind",paste0("cs_rho",1:9))], col = 'red', pch = 17, cex = 0.8)
  
  legend("topleft", horiz = TRUE, bty = "n", 
         legend = c("AR", "CS", "Truth"), pch = c(16, 17, 8), col = c("blue", "red", "black"), cex = 0.9)
  #lines(x = c(rho_true, rho_true), y = c(ystar_int_upper["true"], ystar_int_lower["true"]))
  points(x = rho_true, y = ystar_int_avg["true"], pch = 8)
  dev.off()
}