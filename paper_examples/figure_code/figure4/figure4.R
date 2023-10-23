load("hrs_pred_beta.RData")
my_colors <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


##################
# intercept plot
##################


ylim0 <- range(c(beta_sum_id1[,,1], beta_sum_id2[,,1]))

pdf("hrs_intercept.pdf", width = 6, height = 6)
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
dev.off()

#########
# childhood sep
#########

pdf("hrs_cSEP.pdf", width = 6, height = 6)
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
dev.off()

# HS completition
pdf("hrs_hs.pdf", width = 6, height = 6)
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
dev.off()