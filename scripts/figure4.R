# Figure 4
# Effects of individual covariates for HRS data
load("data/HRS/HRS_all.RData")
load("results/vcbart_hrs_all.RData")

# vector of distinct ages for hypothetical individuals
age_plot <- Z_plot[1:544,"age"]

# Get the betas for each individual
for(id in 1:8){
  tmp_beta <- hrs_vcbart_adapt_cs50$test$beta[544*(id-1) + 1:544, ,]
  dimnames(tmp_beta)[[3]] <- c("int", colnames(X_all))
  
  assign(paste0("ind",id, "_beta"),
         tmp_beta)
}

dimnames(ind1_beta)

# Compare individual 1 in blue (white, unmarried) to individual 7 in red (black, married)

png("figures/hrs_beta.png", width = 8, height = 2, units = "in", res = 400)
par(mar = c(4,3,2,1), mgp = c(1.8, 0.5, 0), mfrow = c(1,4))

# Intercept
plot(1, type = "n", 
     xlim = range(age_plot/12), 
     ylim = range(c(ind1_beta[,c("L95", "U95"), "int"], ind7_beta[,c("L95", "U95"), "int"])),
     xlab = "", ylab = expression(beta[0]), main = "Intercept")
mtext(text = "Age\n(a)", side = 1, line = 3, cex = par("cex"))
polygon(c(age_plot/12, rev(age_plot/12)),
        c(ind1_beta[,"L95","int"], rev(ind1_beta[,"U95", "int"])),
        col = rgb(0,0,1,1/5), border = NA)
polygon(c(age_plot/12, rev(age_plot/12)),
        c(ind7_beta[,"L95","int"], rev(ind7_beta[,"U95", "int"])),
        col = rgb(1,0,0,1/5), border = NA)

lines(age_plot/12, ind1_beta[,"MEAN", "int"], col = rgb(0,0,1,1/3))
lines(age_plot/12, ind7_beta[,"MEAN", "int"], col = rgb(1,0,0,1/3))
abline(h = 0, lty = 2)
axis(side = 1, at = quantile(Z_all[,1]/12, probs = c(0.05, 0.95)), labels = NA)
abline(v = quantile(Z_all[,1]/12, probs = c(0.05, 0.95)), lty = 2, col = 'gray')


# cSEP
plot(1, type = "n", 
     xlim = range(age_plot/12), 
     ylim = range(c(ind1_beta[,c("L95", "U95"), "cSES"], ind7_beta[,c("L95", "U95"), "cSES"])),
     xlab = "", ylab = expression(beta["cSEP"]), main = "Effect of cSEP")
mtext(text = "Age\n(b)", side = 1, line = 3, cex = par("cex"))
polygon(c(age_plot/12, rev(age_plot/12)),
        c(ind1_beta[,"L95","cSES"], rev(ind1_beta[,"U95", "cSES"])),
        col = rgb(0,0,1,1/5), border = NA)
polygon(c(age_plot/12, rev(age_plot/12)),
        c(ind7_beta[,"L95","cSES"], rev(ind7_beta[,"U95", "cSES"])),
        col = rgb(1,0,0,1/5), border = NA)

lines(age_plot/12, ind1_beta[,"MEAN", "cSES"], col = rgb(0,0,1,1/3))
lines(age_plot/12, ind7_beta[,"MEAN", "cSES"], col = rgb(1,0,0,1/3))
abline(h = 0, lty = 2)
axis(side = 1, at = quantile(Z_all[,1]/12, probs = c(0.05, 0.95)), labels = NA)
abline(v = quantile(Z_all[,1]/12, probs = c(0.025, 0.975)), lty = 2, col = 'gray')

# Education
plot(1, type = "n", 
     xlim = range(age_plot/12), 
     ylim = range(c(ind1_beta[,c("L95", "U95"), "education"], ind7_beta[,c("L95", "U95"), "education"])),
     xlab = "", ylab = expression(beta["educ"]), main = "Effect of Education")
mtext(text = "Age\n(c)", side = 1, line = 3, cex = par("cex"))

polygon(c(age_plot/12, rev(age_plot/12)),
        c(ind1_beta[,"L95","education"], rev(ind1_beta[,"U95", "education"])),
        col = rgb(0,0,1,1/5), border = NA)
polygon(c(age_plot/12, rev(age_plot/12)),
        c(ind7_beta[,"L95","education"], rev(ind7_beta[,"U95", "education"])),
        col = rgb(1,0,0,1/5), border = NA)

lines(age_plot/12, ind1_beta[,"MEAN", "education"], col = rgb(0,0,1,1/3))
lines(age_plot/12, ind7_beta[,"MEAN", "education"], col = rgb(1,0,0,1/3))
abline(h = 0, lty = 2)
axis(side = 1, at = quantile(Z_all[,1]/12, probs = c(0.05, 0.95)), labels = NA)
abline(v = quantile(Z_all[,1]/12, probs = c(0.025, 0.975)), lty = 2, col = 'gray')


# Diabetes
plot(1, type = "n", 
     xlim = range(age_plot/12), 
     ylim = range(c(ind1_beta[,c("L95", "U95"), "diabetes"], ind7_beta[,c("L95", "U95"), "diabetes"])),
     xlab = "", ylab = expression(beta["Diabetes"]), main = "Effect of Diabetes")
mtext(text = "Age\n(d)", side = 1, line = 3, cex = par("cex"))

polygon(c(age_plot/12, rev(age_plot/12)),
        c(ind1_beta[,"L95","diabetes"], rev(ind1_beta[,"U95", "diabetes"])),
        col = rgb(0,0,1,1/5), border = NA)
polygon(c(age_plot/12, rev(age_plot/12)),
        c(ind7_beta[,"L95","diabetes"], rev(ind7_beta[,"U95", "diabetes"])),
        col = rgb(1,0,0,1/5), border = NA)

lines(age_plot/12, ind1_beta[,"MEAN", "diabetes"], col = rgb(0,0,1,1/3))
lines(age_plot/12, ind7_beta[,"MEAN", "diabetes"], col = rgb(1,0,0,1/3))
abline(h = 0, lty = 2)
axis(side = 1, at = quantile(Z_all[,1]/12, probs = c(0.05, 0.95)), labels = NA)
abline(v = quantile(Z_all[,1]/12, probs = c(0.025, 0.975)), lty = 2, col = 'gray')

dev.off()

