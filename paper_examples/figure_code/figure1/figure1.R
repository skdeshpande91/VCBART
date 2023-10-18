my_colors <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


pdf("step_new.pdf", width = 3, height = 3)
par(mar = c(3,1,2,1), mgp = c(1.8, 0.5,0), mfrow = c(2,1),
    cex.lab = 0.8, cex.axis = 0.7, cex.main = 1)
plot(1, type = "n",
     main = expression(Z[2] %in% "{"~c[2]~","~c[5]~"}"),
     xlab = expression(Z[1]),
     yaxt = "n",
     xlim = c(0,1), ylim = c(-0.75, 0.75))
lines(x = c(0, 0.7), y = c(-0.3, -0.3), col = my_colors[1], lwd = 2)
lines(x = c(0.7, 1), y = c(0.4, 0.4), col = my_colors[8], lwd = 2)
lines(x = c(0.7, 0.7), y = c(-0.3, 0.4), lty = 2, lwd = 0.5)

text(x = 0.35, y = 0, labels = expression(mu[3]))
text(x = 0.85, y = 0.1, labels = expression(mu[4]))


plot(1, type = "n",
     main = expression(Z[2] %in% "{"~c[1]~","~c[3]~","~c[4]~"}"),
     xlab = expression(Z[1]),
     yaxt = "n",
     xlim = c(0,1), ylim = c(-0.75, 0.75))
lines(x = c(0.7, 1), y = c(0.4, 0.4), col = my_colors[8], lwd = 2)
text(x = 0.85, y = 0.1, labels = expression(mu[4]))

lines(x = c(0, 0.3), y = c(0.3, 0.3), col = my_colors[2], lwd = 2)
text(x = 0.15, y = 0.05, labels = expression(mu[1]))

lines(x = c(0.3, 0.7), y = c(0, 0), col = my_colors[3], lwd = 2)
text(x = 0.5, y = -0.25, labels = expression(mu[2]))

lines(x = c(0.3, 0.3), y= c(0.3, 0), lty = 2, lwd = 0.5)
lines(x = c(0.7, 0.7), y = c(0,0.4), lty = 2, lwd = 0.5)

dev.off()