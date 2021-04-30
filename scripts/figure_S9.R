# Script to produce Figure S9
load("results/hrs_causal_subgroup_samples.RData")


##############
# Table of ATE means and credible intervals

comps <- c("1-0", "2-0", "3-0", "2-1", "3-1", "3-2")
ate_table <- data.frame("MEAN" = rep(NA, times = length(comps)),
                        "L95" = rep(NA, times = length(comps)),
                        "U95" = rep(NA, times = length(comps)),
                        "Pneg" = rep(NA, times = length(comps)))
rownames(ate_table) <- comps

for(t in 1:3){
  comp_name <- paste0(t, "-", 0)
  tmp_samples <- ate_samples[t+1,]
  ate_table[comp_name,"MEAN"] <- round(mean(tmp_samples), digits = 2)
  ate_table[comp_name, c("L95", "U95")] <- round(quantile(tmp_samples, probs = c(0.025, 0.975)), digits = 2)
  ate_table[comp_name, "Pneg"] <- round(mean(tmp_samples < 0) * 100, digits = 1)
}


for(t in 2:3){
  for(tt in 1:(t-1)){
    comp_name <- paste0(t, "-", tt)
    tmp_samples <- ate_samples[t+1,] - ate_samples[tt+1,]
    ate_table[comp_name,"MEAN"] <- round(mean(tmp_samples), digits = 2)
    ate_table[comp_name, c("L95", "U95")] <- round(quantile(tmp_samples, probs = c(0.025, 0.975)), digits = 2)
    ate_table[comp_name, "Pneg"] <- round(100 * mean(tmp_samples < 0), digits = 1)
    
    
  }
}




# Effect of inactivity and depression (relative to healthy controls)


summary(educ_samples[[1]][4,]) # LT HS
summary(educ_samples[[4]][4,]) # HS only


beta3_educ1_dens <- density(educ_samples[[1]][4,])
beta3_educ2_dens <- density(educ_samples[[2]][4,])
beta3_educ3_dens <- density(educ_samples[[3]][4,])
beta3_educ4_dens <- density(educ_samples[[4]][4,])

beta3_age2_dens <- density(age_samples[[2]][4,])
beta3_age3_dens <- density(age_samples[[3]][4,])
beta3_age4_dens <- density(age_samples[[4]][4,])
beta3_age5_dens <- density(age_samples[[5]][4,])


y_max_educ <- max(c(beta3_educ1_dens$y, beta3_educ2_dens$y, beta3_educ3_dens$y, beta3_educ4_dens$y))
x_range_educ <- range(c(beta3_educ1_dens$x, beta3_educ2_dens$x, beta3_educ3_dens$x, beta3_educ4_dens$x))


y_max_age <- max(c(beta3_age2_dens$y, beta3_age3_dens$y, beta3_age4_dens$y, beta3_age5_dens$y))
x_range_age <- range(c(beta3_age2_dens$x, beta3_age3_dens$x, beta3_age4_dens$x, beta3_age5_dens$x))


png("figures/hrs_causal_beta3_dens.png", width = 8, height = 4, units = "in", res = 400)
par(mar = c(4.2,3,2,1), mgp = c(1.8, 0.5, 0), mfrow = c(1,2))
plot(1, type = "n", main = "Effect of inactivity & depression", xlab = "", ylab = "Density",
     xlim = x_range_educ, ylim = c(1e-3, y_max_educ*1.10))
mtext(text = "Effect\n(a)",side = 1, line = 3)
legend("topleft", legend = c("Less than HS", "HS degree & no college", "College degree"), 
       lty = 1, col = c("red", "orange", "blue"), bty = "n", cex = 0.8)
lines(beta3_educ1_dens$x, beta3_educ1_dens$y, col = 'red')
lines(beta3_educ2_dens$x, beta3_educ2_dens$y, col = 'orange')
#lines(beta3_educ3_dens$x, beta3_educ3_dens$y, col = 'black')
lines(beta3_educ4_dens$x, beta3_educ4_dens$y, col = 'blue')


plot(1, type = "n", main = "Effect of inactivity & depression", xlab = "", ylab = "Density",
     xlim = x_range_age, ylim = c(1e-3, y_max_age*1.10))
mtext(text = "Effect\n(b)",side = 1, line = 3)

legend("topleft", legend = c("65 - 69 years", "70 - 74 years", "80 - 85 years"), 
       lty = 1, col = c("black", "purple", "cyan"), bty = "n", cex = 0.8)
lines(beta3_age2_dens$x, beta3_age2_dens$y, col = 'black')
lines(beta3_age3_dens$x, beta3_age3_dens$y, col = 'purple')
#lines(beta3_educ3_dens$x, beta3_educ3_dens$y, col = 'black')
lines(beta3_age5_dens$x, beta3_age5_dens$y, col = 'cyan')
dev.off()

