load("hrs_pred_beta.RData")
my_colors <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


##################
# intercept plot
##################
ylim0 <- range(c(beta_sum_id1[,,1], beta_sum_id2[,,1]))

par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0))
plot(1, type = "n", 
     xlim = range(age_seq), ylim = ylim0,
     xlab = expression("Age (months)"), 
     ylab = expression ("Cognitive score"),
     main = expression("Intercept"))
polygon(x = c(age_seq, rev(age_seq)),
        y = c(beta_sum_id1[,"L95",1], rev(beta_sum_id1[,"U95",1])),
        border = NA, col = adjustcolor(my_colors[1], alpha.f = 0.35))

polygon(x = c(age_seq, rev(age_seq)),
        y = c(beta_sum_id2[,"L95",1], rev(beta_sum_id2[,"U95",1])),
        border = NA, col = adjustcolor(my_colors[3], alpha.f = 0.35))

lines(x = age_seq, y = beta_sum_id1[,"MEAN",1], col = my_colors[1], lwd = 2.5)
lines(x = age_seq, y = beta_sum_id2[,"MEAN",1], col = my_colors[3], lwd = 2.5)

abline(v = age_quants, col = my_colors[8], lty = 2, lwd = 1)

legend("bottomleft", legend = c("Black", "White"), col = my_colors[c(3,1)],
       lwd = 2.5, lty = 1, bty = "n")

par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0))
plot(1, type = "n",
     xlim = range(age_seq), ylim = c(-1.01, 1.01) * max(abs(diff_sum[,,1])),
     ylab = expression("Difference"), xlab = expression("Age (months)"),
     main = expression("Difference in Intercepts"))
polygon(x = c(age_seq, rev(age_seq)),
        y = c(diff_sum[,"L95",1], rev(diff_sum[,"U95",1])),
        border = NA, col = adjustcolor(my_colors[2], alpha.f = 0.25))
lines(x = age_seq, y = diff_sum[,"MEAN",1], col = my_colors[2], lwd = 2.5)
abline(h = 0, col = my_colors[1], lty = 2)
#########
# childhood sep
#########

par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0))
ylim1 <- range(c(beta_sum_id1[,,2], beta_sum_id2[,,2]))
plot(1, type = "n", 
     xlim = range(age_seq), ylim = ylim1,
     xlab = expression("Age (months)"), 
     ylab = expression ("Cognitive score"),
     main = expression("Childhood SEP"))
polygon(x = c(age_seq, rev(age_seq)),
        y = c(beta_sum_id1[,"L95",2], rev(beta_sum_id1[,"U95",2])),
        border = NA, col = adjustcolor(my_colors[1], alpha.f = 0.35))

polygon(x = c(age_seq, rev(age_seq)),
        y = c(beta_sum_id2[,"L95",2], rev(beta_sum_id2[,"U95",2])),
        border = NA, col = adjustcolor(my_colors[3], alpha.f = 0.35))

lines(x = age_seq, y = beta_sum_id1[,"MEAN",2], col = my_colors[1], lwd = 2.5)
lines(x = age_seq, y = beta_sum_id2[,"MEAN",2], col = my_colors[3], lwd = 2.5)

abline(v = age_quants, col = my_colors[8], lty = 2, lwd = 1)

legend("bottomleft", legend = c("Black", "White"), col = my_colors[c(3,1)],
       lwd = 2.5, lty = 1, bty = "n")

par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0))
plot(1, type = "n",
     xlim = range(age_seq), ylim = c(-1.01, 1.01) * max(abs(diff_sum[,,2])),
     ylab = expression("Difference"), xlab = expression("Age (months)"),
     main = expression("Difference in cSEP effect"))
polygon(x = c(age_seq, rev(age_seq)),
        y = c(diff_sum[,"L95",2], rev(diff_sum[,"U95",2])),
        border = NA, col = adjustcolor(my_colors[2], alpha.f = 0.25))
lines(x = age_seq, y = diff_sum[,"MEAN",2], col = my_colors[2], lwd = 2.5)
abline(h = 0, col = my_colors[1], lty = 2)

# HS completition
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0))
ylim2 <- range(c(beta_sum_id1[,,3], beta_sum_id2[,,3]))
plot(1, type = "n", 
     xlim = range(age_seq), ylim = ylim2,
     xlab = expression("Age (months)"), 
     ylab = expression ("Cognitive score"),
     main = expression("High school completition"))
polygon(x = c(age_seq, rev(age_seq)),
        y = c(beta_sum_id1[,"L95",3], rev(beta_sum_id1[,"U95",3])),
        border = NA, col = adjustcolor(my_colors[1], alpha.f = 0.35))

polygon(x = c(age_seq, rev(age_seq)),
        y = c(beta_sum_id2[,"L95",3], rev(beta_sum_id2[,"U95",3])),
        border = NA, col = adjustcolor(my_colors[3], alpha.f = 0.35))

lines(x = age_seq, y = beta_sum_id1[,"MEAN",3], col = my_colors[1], lwd = 2.5)
lines(x = age_seq, y = beta_sum_id2[,"MEAN",3], col = my_colors[3], lwd = 2.5)

abline(v = age_quants, col = my_colors[8], lty = 2, lwd = 1)

legend("bottomleft", legend = c("Black", "White"), col = my_colors[c(3,1)],
       lwd = 2.5, lty = 1, bty = "n")

par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0))
plot(1, type = "n",
     xlim = range(age_seq), ylim = c(-1.01, 1.01) * max(abs(diff_sum[,,2])),
     ylab = expression("Difference"), xlab = expression("Age (months)"),
     main = expression("Difference in HS effect"))
polygon(x = c(age_seq, rev(age_seq)),
        y = c(diff_sum[,"L95",3], rev(diff_sum[,"U95",3])),
        border = NA, col = adjustcolor(my_colors[2], alpha.f = 0.25))
lines(x = age_seq, y = diff_sum[,"MEAN",3], col = my_colors[2], lwd = 2.5)
abline(h = 0, col = my_colors[1], lty = 2)



