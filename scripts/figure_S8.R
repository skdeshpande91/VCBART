# Make a 1 1 x 3 panels love plot of HRS causal imbalances

load("results/HRS_causal_balance.RData")

vars <- c("Age", "Female", "NH White", "NH Black", "Hispanic", "cSEP", "Educ", 
          "Poor chld health", "Excellent chld health", "Smoked")




L <- length(vars)

png("~/Documents/Research/vc_bart/figures/hrs_love_plot.png", width = 7.5, height = 2.5, units = "in", res= 400)
par(mar = c(3,1,2,1), mgp = c(1.8, 0.5, 0), mfrow = c(1,3))

plot(1, xaxt = "n", yaxt = "n", xlim = c(-1, 0.6), ylim = c(0, L+1), ylab = "",
     main = "T = 1 vs T = 0", xlab = "Standardized Difference")
abline(v = c(-0.1, 0, 0.1), col = 'gray')
for(l in 1:L){
  lines(c(std_diff_before[vars[l], "1-0"], std_diff_after[vars[l], "1-0"]), c(l,l),
        col = 'lightgray', lty = 2)
  points(c(std_diff_before[vars[l],"1-0"], std_diff_after[vars[l],"1-0"]),
         c(l,l), pch = c(4,16), col = c(rgb(1,0,0,1/2), rgb(0,0,1,1/2)))
  text(-0.7, l, labels = vars[l], cex = 0.9)
}
axis(side = 1, at = seq(-0.6, 0.6, by = 0.2), labels = round(seq(-0.6, 0.6, by = 0.2), digits = 1))


plot(1, xaxt = "n", yaxt = "n", xlim = c(-1, 0.6), ylim = c(0, L+1), ylab = "",
     main = "T = 2 vs T = 0", xlab = "Standardized Difference")
abline(v = c(-0.1, 0, 0.1), col = 'gray')
for(l in 1:L){
  lines(c(std_diff_before[vars[l], "2-0"], std_diff_after[vars[l], "2-0"]), c(l,l),
        col = 'lightgray', lty = 2)
  points(c(std_diff_before[vars[l],"2-0"], std_diff_after[vars[l],"2-0"]),
         c(l,l), pch = c(4,16), col = c(rgb(1,0,0,1/2), rgb(0,0,1,1/2)))
  text(-0.7, l, labels = vars[l], cex = 0.9)
}
axis(side = 1, at = seq(-0.6, 0.6, by = 0.2), labels = round(seq(-0.6, 0.6, by = 0.2), digits = 1))

plot(1, xaxt = "n", yaxt = "n", xlim = c(-1, 0.6), ylim = c(0, L+1), ylab = "",
     main = "T = 3 vs T = 0", xlab = "Standardized Difference")
abline(v = c(-0.1, 0, 0.1), col = 'gray')
for(l in 1:L){
  lines(c(std_diff_before[vars[l], "3-0"], std_diff_after[vars[l], "3-0"]), c(l,l),
        col = 'lightgray', lty = 2)
  points(c(std_diff_before[vars[l],"3-0"], std_diff_after[vars[l],"3-0"]),
         c(l,l), pch = c(4,16), col = c(rgb(1,0,0,1/2), rgb(0,0,1,1/2)))
  text(-0.7, l, labels = vars[l], cex = 0.9)
}
axis(side = 1, at = seq(-0.6, 0.6, by = 0.2), labels = round(seq(-0.6, 0.6, by = 0.2), digits = 1))

dev.off()




# 19 people 