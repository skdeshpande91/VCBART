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

# Compare individual 4 in blue (white, married, foodstamps) to individual 5 in red (black, unmarried, no foodstamps)

# Plot effects of intercept, cSEP, education, diabetes, BMI, depression



png("~/Documents/Research/vc_bart/figures/hrs_beta.png", width = 6, height = 4, units = "in", res = 300)
par(mar = c(4,3,2,1), mgp = c(1.8, 0.5, 0), mfrow = c(2,3))

# Intercept
plot(1, type = "n", 
     xlim = range(age_plot/12), 
     ylim = range(c(ind4_beta[,c("L95", "U95"), "int"], ind5_beta[,c("L95", "U95"), "int"])),
     xlab = "", ylab = expression(beta[0]), main = "Intercept")
mtext(text = "Age\n(a)", side = 1, line = 3, cex = par("cex"))
polygon(c(age_plot/12, rev(age_plot/12)),
        c(ind4_beta[,"L95","int"], rev(ind4_beta[,"U95", "int"])),
        col = rgb(0,0,1,1/5), border = NA)
polygon(c(age_plot/12, rev(age_plot/12)),
        c(ind5_beta[,"L95","int"], rev(ind5_beta[,"U95", "int"])),
        col = rgb(1,0,0,1/5), border = NA)

lines(age_plot/12, ind4_beta[,"MEAN", "int"], col = rgb(0,0,1,1/3))
lines(age_plot/12, ind5_beta[,"MEAN", "int"], col = rgb(1,0,0,1/3))
abline(h = 0, lty = 2)
axis(side = 1, at = quantile(Z_all[,1]/12, probs = c(0.05, 0.95)), labels = NA)
#abline(v = quantile(Z_all[,1]/12, probs = c(0.025, 0.975)), lty = 2)


# cSEP
plot(1, type = "n", 
     xlim = range(age_plot/12), 
     ylim = range(c(ind4_beta[,c("L95", "U95"), "cSES"], ind5_beta[,c("L95", "U95"), "cSES"])),
     xlab = "", ylab = expression(beta["cSEP"]), main = "Effect of cSEP")
mtext(text = "Age\n(b)", side = 1, line = 3, cex = par("cex"))
polygon(c(age_plot/12, rev(age_plot/12)),
        c(ind4_beta[,"L95","cSES"], rev(ind4_beta[,"U95", "cSES"])),
        col = rgb(0,0,1,1/5), border = NA)
polygon(c(age_plot/12, rev(age_plot/12)),
        c(ind5_beta[,"L95","cSES"], rev(ind5_beta[,"U95", "cSES"])),
        col = rgb(1,0,0,1/5), border = NA)

lines(age_plot/12, ind4_beta[,"MEAN", "cSES"], col = rgb(0,0,1,1/3))
lines(age_plot/12, ind5_beta[,"MEAN", "cSES"], col = rgb(1,0,0,1/3))
abline(h = 0, lty = 2)
axis(side = 1, at = quantile(Z_all[,1]/12, probs = c(0.05, 0.95)), labels = NA)

# Education
plot(1, type = "n", 
     xlim = range(age_plot/12), 
     ylim = range(c(ind4_beta[,c("L95", "U95"), "education"], ind5_beta[,c("L95", "U95"), "education"])),
     xlab = "", ylab = expression(beta["educ"]), main = "Effect of Education")
mtext(text = "Age\n(c)", side = 1, line = 3, cex = par("cex"))

polygon(c(age_plot/12, rev(age_plot/12)),
        c(ind4_beta[,"L95","education"], rev(ind4_beta[,"U95", "education"])),
        col = rgb(0,0,1,1/5), border = NA)
polygon(c(age_plot/12, rev(age_plot/12)),
        c(ind5_beta[,"L95","education"], rev(ind5_beta[,"U95", "education"])),
        col = rgb(1,0,0,1/5), border = NA)

lines(age_plot/12, ind4_beta[,"MEAN", "education"], col = rgb(0,0,1,1/3))
lines(age_plot/12, ind5_beta[,"MEAN", "education"], col = rgb(1,0,0,1/3))
abline(h = 0, lty = 2)
axis(side = 1, at = quantile(Z_all[,1]/12, probs = c(0.05, 0.95)), labels = NA)

# BMI
plot(1, type = "n", 
     xlim = range(age_plot/12), 
     ylim = range(c(ind4_beta[,c("L95", "U95"), "bmi"], ind5_beta[,c("L95", "U95"), "bmi"])),
     xlab = "", ylab = expression(beta["BMI"]), main = "Effect of BMI")
mtext(text = "Age\n(d)", side = 1, line = 3, cex = par("cex"))

polygon(c(age_plot/12, rev(age_plot/12)),
        c(ind4_beta[,"L95","bmi"], rev(ind4_beta[,"U95", "bmi"])),
        col = rgb(0,0,1,1/5), border = NA)
polygon(c(age_plot/12, rev(age_plot/12)),
        c(ind5_beta[,"L95","bmi"], rev(ind5_beta[,"U95", "bmi"])),
        col = rgb(1,0,0,1/5), border = NA)

lines(age_plot/12, ind4_beta[,"MEAN", "bmi"], col = rgb(0,0,1,1/3))
lines(age_plot/12, ind5_beta[,"MEAN", "bmi"], col = rgb(1,0,0,1/3))
abline(h = 0, lty = 2)
axis(side = 1, at = quantile(Z_all[,1]/12, probs = c(0.05, 0.95)), labels = NA)


# CESD
plot(1, type = "n", 
     xlim = range(age_plot/12), 
     ylim = range(c(ind4_beta[,c("L95", "U95"), "bmi"], ind5_beta[,c("L95", "U95"), "cesd"])),
     xlab = "", ylab = expression(beta["CESD"]), main = "Effect of CESD")
mtext(text = "Age\n(d)", side = 1, line = 3, cex = par("cex"))

polygon(c(age_plot/12, rev(age_plot/12)),
        c(ind4_beta[,"L95","cesd"], rev(ind4_beta[,"U95", "cesd"])),
        col = rgb(0,0,1,1/5), border = NA)
polygon(c(age_plot/12, rev(age_plot/12)),
        c(ind5_beta[,"L95","cesd"], rev(ind5_beta[,"U95", "cesd"])),
        col = rgb(1,0,0,1/5), border = NA)

lines(age_plot/12, ind4_beta[,"MEAN", "cesd"], col = rgb(0,0,1,1/3))
lines(age_plot/12, ind5_beta[,"MEAN", "cesd"], col = rgb(1,0,0,1/3))
abline(h = 0, lty = 2)
axis(side = 1, at = quantile(Z_all[,1]/12, probs = c(0.05, 0.95)), labels = NA)

# Diabetes
plot(1, type = "n", 
     xlim = range(age_plot/12), 
     ylim = range(c(ind4_beta[,c("L95", "U95"), "diabetes"], ind5_beta[,c("L95", "U95"), "diabetes"])),
     xlab = "", ylab = expression(beta["Diabetes"]), main = "Effect of Diabetes")
mtext(text = "Age\n(d)", side = 1, line = 3, cex = par("cex"))

polygon(c(age_plot/12, rev(age_plot/12)),
        c(ind4_beta[,"L95","diabetes"], rev(ind4_beta[,"U95", "diabetes"])),
        col = rgb(0,0,1,1/5), border = NA)
polygon(c(age_plot/12, rev(age_plot/12)),
        c(ind5_beta[,"L95","diabetes"], rev(ind5_beta[,"U95", "diabetes"])),
        col = rgb(1,0,0,1/5), border = NA)

lines(age_plot/12, ind4_beta[,"MEAN", "diabetes"], col = rgb(0,0,1,1/3))
lines(age_plot/12, ind5_beta[,"MEAN", "diabetes"], col = rgb(1,0,0,1/3))
abline(h = 0, lty = 2)
axis(side = 1, at = quantile(Z_all[,1]/12, probs = c(0.05, 0.95)), labels = NA)

dev.off()



# Individual plots
# intercept
png("figures/hrs_intercept.png", width = 6, height = 6, units = "in", res = 300)
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0))
plot(1, type = "n", 
     xlim = range(age_plot), 
     ylim = range(c(ind1_beta[,c("L95", "U95"), "int"], ind7_beta[,c("L95", "U95"), "int"])),
     xlab = "Age", ylab = expression(beta[0]), main = "Intercept")
polygon(c(age_plot, rev(age_plot)),
        c(ind1_beta[,"L95","int"], rev(ind1_beta[,"U95", "int"])),
        col = rgb(0,0,1,1/5), border = NA)
polygon(c(age_plot, rev(age_plot)),
        c(ind7_beta[,"L95","int"], rev(ind7_beta[,"U95", "int"])),
        col = rgb(1,0,0,1/5), border = NA)
lines(age_plot, ind7_beta[,"MEAN", "int"], col = rgb(1,0,0,1/3))
lines(age_plot, ind1_beta[,"MEAN", "int"], col = rgb(0,0,1,1/3))
abline(h = 0, lty = 2)
dev.off()

# cses
png("figures/hrs_cses.png", width = 6, height = 6, units = "in", res = 300)
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0))

plot(1, type = "n", 
     xlim = range(age_plot), 
     ylim = range(c(ind1_beta[,c("L95", "U95"), "cSES"], ind7_beta[,c("L95", "U95"), "cSES"])),
     xlab = "Age", ylab = expression(beta["cSEP"]), main = "Effect of childhood SEP")
polygon(c(age_plot, rev(age_plot)),
        c(ind1_beta[,"L95","cSES"], rev(ind1_beta[,"U95", "cSES"])),
        col = rgb(0,0,1,1/5), border = NA)
polygon(c(age_plot, rev(age_plot)),
        c(ind7_beta[,"L95","cSES"], rev(ind7_beta[,"U95", "cSES"])),
        col = rgb(1,0,0,1/5), border = NA)
lines(age_plot, ind7_beta[,"MEAN", "cSES"], col = rgb(1,0,0,1/3))
lines(age_plot, ind1_beta[,"MEAN", "cSES"], col = rgb(0,0,1,1/3))
dev.off()

# education
png("figures/hrs_educ.png", width = 6, height = 6, units = "in", res = 300)
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0))

plot(1, type = "n", 
     xlim = range(age_plot), 
     ylim = range(c(0,ind1_beta[,c("L95", "U95"), "education"], ind7_beta[,c("L95", "U95"), "education"])),
     xlab = "Age", ylab = expression(beta["education"]), main = "Effect of education")
polygon(c(age_plot, rev(age_plot)),
        c(ind1_beta[,"L95","education"], rev(ind1_beta[,"U95", "education"])),
        col = rgb(0,0,1,1/3), border = NA)
polygon(c(age_plot, rev(age_plot)),
        c(ind7_beta[,"L95","education"], rev(ind7_beta[,"U95", "education"])),
        col = rgb(1,0,0,1/3), border = NA)
lines(age_plot, ind7_beta[,"MEAN", "education"], col = 'red')
lines(age_plot, ind1_beta[,"MEAN", "education"], col = 'blue')
dev.off()

# waelth
png("figures/hrs_wealth.png", width = 6, height = 6, units = "in", res = 300)
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0))

plot(1, type = "n",
     xlim = range(age_plot),
     ylim = 1e5 * range(c(ind1_beta[,c("L95", "U95"), "wealth"], ind7_beta[,c("L95", "U95"), "wealth"])),
     xlab = "Age", ylab = expression(beta["wealth"]), main = "Effect of wealth")
polygon(x = c(age_plot, rev(age_plot)),
        y = c(ind1_beta[,"L95", "wealth"], rev(ind1_beta[,"U95", "wealth"])) * 1e5,
        col = rgb(0,0,1,1/5), border = NA)
polygon(x = c(age_plot, rev(age_plot)),
        y = c(ind7_beta[,"L95", "wealth"], rev(ind7_beta[,"U95", "wealth"])) * 1e5,
        col = rgb(1,0,0,1/5), border = NA)
lines(age_plot, 1e5*ind1_beta[,"MEAN", "wealth"], col = rgb(0,0,1,1/2))
lines(age_plot, 1e5*ind7_beta[,"MEAN", "wealth"], col = rgb(1,0,0,1/2))

dev.off()


# diabetes
png("figures/hrs_diab.png", width = 6, height = 6, units = "in", res = 300)
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0))

plot(1, type = "n", 
     xlim = range(age_plot), 
     ylim = range(c(ind1_beta[,c("L95","U95"), "diabetes"], ind7_beta[,c("L95", "U95"), "diabetes"])), 
     xlab = "Age (months)", ylab = expression(beta["diabetes"]),
     main = "Effect of diabetes")
abline(h = 0, lty = 2)
polygon(c(age_plot, rev(age_plot)), c(ind1_beta[,"L95","diabetes"], rev(ind1_beta[,"U95", "diabetes"])),
        col = rgb(0,0,1,1/5), border = NA)
lines(age_plot, ind1_beta[,"MEAN", "diabetes"], col = 'blue')
polygon(c(age_plot, rev(age_plot)), c(ind7_beta[,"L95", "diabetes"], rev(ind7_beta[,"U95", "diabetes"])),
        col = rgb(1,0,0,1/5), border = NA)
lines(age_plot, ind7_beta[,"MEAN", "diabetes"], col = 'red')
dev.off()


# BMI
png("figures/hrs_bmi.png", width = 6, height = 6, units = "in", res = 300)
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0))

plot(1, type = "n", 
     xlim = range(age_plot), 
     ylim = range(c(ind1_beta[,c("L95","U95"), "bmi"], ind7_beta[,c("L95", "U95"), "bmi"])), 
     xlab = "Age (months)", ylab = expression(beta["bmi"]),
     main = "Effect of BMI")
abline(h = 0, lty = 2)
polygon(c(age_plot, rev(age_plot)), c(ind1_beta[,"L95","bmi"], rev(ind1_beta[,"U95", "bmi"])),
        col = rgb(0,0,1,1/5), border = NA)
lines(age_plot, ind1_beta[,"MEAN", "bmi"], col = 'blue')
polygon(c(age_plot, rev(age_plot)), c(ind7_beta[,"L95", "bmi"], rev(ind7_beta[,"U95", "bmi"])),
        col = rgb(1,0,0,1/5), border = NA)
lines(age_plot, ind7_beta[,"MEAN", "bmi"], col = 'red')
dev.off()

# physical activity
png("figures/hrs_phys_act.png", width = 6, height = 6, units = "in", res = 300)
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0))

plot(1, type = "n", 
     xlim = range(age_plot), 
     ylim = range(c(ind1_beta[,c("L95", "U95"), "phys_activity"], ind7_beta[,c("L95", "U95"), "phys_activity"])),
     xlab = "Age", ylab = expression(beta["phys_act"]), main = "Effect of Physical Activity")
polygon(c(age_plot, rev(age_plot)),
        c(ind1_beta[,"L95","phys_activity"], rev(ind1_beta[,"U95", "phys_activity"])),
        col = rgb(0,0,1,1/5), border = NA)
polygon(c(age_plot, rev(age_plot)),
        c(ind7_beta[,"L95","phys_activity"], rev(ind7_beta[,"U95", "phys_activity"])),
        col = rgb(1,0,0,1/5), border = NA)
lines(age_plot, ind7_beta[,"MEAN", "phys_activity"], col = rgb(1,0,0,1/3))
lines(age_plot, ind1_beta[,"MEAN", "phys_activity"], col = rgb(0,0,1,1/3))
dev.off()

# Hi BP
png("figures/hrs_hibp.png", width = 6, height = 6, units = "in", res = 300)
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0))

plot(1, type = "n", 
     xlim = range(age_plot), 
     ylim = range(c(ind1_beta[,c("L95", "U95"), "hi_bp"], ind7_beta[,c("L95", "U95"), "hi_bp"])),
     xlab = "Age", ylab = expression(beta["hi_bp"]), main = "Effect of High Blood Pressure")
polygon(c(age_plot, rev(age_plot)),
        c(ind1_beta[,"L95","hi_bp"], rev(ind1_beta[,"U95", "hi_bp"])),
        col = rgb(0,0,1,1/5), border = NA)
polygon(c(age_plot, rev(age_plot)),
        c(ind7_beta[,"L95","hi_bp"], rev(ind7_beta[,"U95", "hi_bp"])),
        col = rgb(1,0,0,1/5), border = NA)
lines(age_plot, ind7_beta[,"MEAN", "hi_bp"], col = rgb(1,0,0,1/3))
lines(age_plot, ind1_beta[,"MEAN", "hi_bp"], col = rgb(0,0,1,1/3))

dev.off()
# loneliness
png("figures/hrs_loneliness.png", width = 6, height = 6, units = "in", res = 300)
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0))
plot(1, type = "n", 
     xlim = range(age_plot), 
     ylim = range(c(ind1_beta[,c("L95", "U95"), "loneliness"], ind7_beta[,c("L95", "U95"), "loneliness"])),
     xlab = "Age", ylab = expression(beta["loneliness"]), main = "Effect of Loneliness")
polygon(c(age_plot, rev(age_plot)),
        c(ind1_beta[,"L95","loneliness"], rev(ind1_beta[,"U95", "loneliness"])),
        col = rgb(0,0,1,1/5), border = NA)
polygon(c(age_plot, rev(age_plot)),
        c(ind7_beta[,"L95","loneliness"], rev(ind7_beta[,"U95", "loneliness"])),
        col = rgb(1,0,0,1/5), border = NA)
lines(age_plot, ind7_beta[,"MEAN", "loneliness"], col = rgb(1,0,0,1/3))
lines(age_plot, ind1_beta[,"MEAN", "loneliness"], col = rgb(0,0,1,1/3))
dev.off()

png("figures/hrs_stroke.png", width = 6, height = 6, units = "in", res = 300)
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0))
plot(1, type = "n", 
     xlim = range(age_plot), 
     ylim = range(c(ind1_beta[,c("L95", "U95"), "stroke"], ind7_beta[,c("L95", "U95"), "stroke"])),
     xlab = "Age", ylab = expression(beta["stroke"]), main = "Effect of Stroke")
polygon(c(age_plot, rev(age_plot)),
        c(ind1_beta[,"L95","stroke"], rev(ind1_beta[,"U95", "stroke"])),
        col = rgb(0,0,1,1/5), border = NA)
polygon(c(age_plot, rev(age_plot)),
        c(ind7_beta[,"L95","stroke"], rev(ind7_beta[,"U95", "stroke"])),
        col = rgb(1,0,0,1/5), border = NA)
lines(age_plot, ind7_beta[,"MEAN", "stroke"], col = rgb(1,0,0,1/3))
lines(age_plot, ind1_beta[,"MEAN", "stroke"], col = rgb(0,0,1,1/3))
dev.off()

png("figures/hrs_cesd.png", width = 6, height = 6, units = "in", res = 300)
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0))
plot(1, type = "n", 
     xlim = range(age_plot), 
     ylim = range(c(ind1_beta[,c("L95", "U95"), "cesd"], ind7_beta[,c("L95", "U95"), "cesd"])),
     xlab = "Age", ylab = expression(beta["cesd"]), main = "Effect of CESD")
polygon(c(age_plot, rev(age_plot)),
        c(ind1_beta[,"L95","cesd"], rev(ind1_beta[,"U95", "cesd"])),
        col = rgb(0,0,1,1/5), border = NA)
polygon(c(age_plot, rev(age_plot)),
        c(ind7_beta[,"L95","cesd"], rev(ind7_beta[,"U95", "cesd"])),
        col = rgb(1,0,0,1/5), border = NA)
lines(age_plot, ind7_beta[,"MEAN", "cesd"], col = rgb(1,0,0,1/3))
lines(age_plot, ind1_beta[,"MEAN", "cesd"], col = rgb(0,0,1,1/3))
abline(h = 0, lty = 2)
dev.off()


png("figures/hrs_heart.png", width = 6, height = 6, units = "in", res = 300)
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0))
plot(1, type = "n", 
     xlim = range(age_plot), 
     ylim = range(c(ind1_beta[,c("L95", "U95"), "heart_problems"], ind7_beta[,c("L95", "U95"), "heart_problems"])),
     xlab = "Age", ylab = expression(beta["heart_problems"]), main = "Effect of Heart problems")
polygon(c(age_plot, rev(age_plot)),
        c(ind1_beta[,"L95","heart_problems"], rev(ind1_beta[,"U95", "heart_problems"])),
        col = rgb(0,0,1,1/5), border = NA)
polygon(c(age_plot, rev(age_plot)),
        c(ind7_beta[,"L95","heart_problems"], rev(ind7_beta[,"U95", "heart_problems"])),
        col = rgb(1,0,0,1/5), border = NA)
lines(age_plot, ind7_beta[,"MEAN", "heart_problems"], col = rgb(1,0,0,1/3))
lines(age_plot, ind1_beta[,"MEAN", "heart_problems"], col = rgb(0,0,1,1/3))
abline(h = 0, lty = 2)
dev.off()


# cSES
png("~/Documents/Research/vc_bart/figures/hrs_betas.png",
    width = 6, height = 6, units = "in", res = 300)

par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0), mfrow = c(2,2))

plot(1, type = "n", 
     xlim = range(age_plot), 
     ylim = range(c(ind1_beta[,c("L95", "U95"), "cSES"], ind7_beta[,c("L95", "U95"), "cSES"])),
     xlab = "Age", ylab = expression(beta["cSEP"]), main = "Effect of childhood SEP")
polygon(c(age_plot, rev(age_plot)),
        c(ind1_beta[,"L95","cSES"], rev(ind1_beta[,"U95", "cSES"])),
        col = rgb(0,0,1,1/5), border = NA)
polygon(c(age_plot, rev(age_plot)),
        c(ind7_beta[,"L95","cSES"], rev(ind7_beta[,"U95", "cSES"])),
        col = rgb(1,0,0,1/5), border = NA)
lines(age_plot, ind7_beta[,"MEAN", "cSES"], col = rgb(1,0,0,1/3))
lines(age_plot, ind1_beta[,"MEAN", "cSES"], col = rgb(0,0,1,1/3))

# education
plot(1, type = "n", 
     xlim = range(age_plot), 
     ylim = range(c(0,ind1_beta[,c("L95", "U95"), "education"], ind7_beta[,c("L95", "U95"), "education"])),
     xlab = "Age", ylab = expression(beta["education"]), main = "Effect of education")
polygon(c(age_plot, rev(age_plot)),
        c(ind1_beta[,"L95","education"], rev(ind1_beta[,"U95", "education"])),
        col = rgb(0,0,1,1/3), border = NA)
polygon(c(age_plot, rev(age_plot)),
        c(ind7_beta[,"L95","education"], rev(ind7_beta[,"U95", "education"])),
        col = rgb(1,0,0,1/3), border = NA)
lines(age_plot, ind7_beta[,"MEAN", "education"], col = 'red')
lines(age_plot, ind1_beta[,"MEAN", "education"], col = 'blue')

# wealth: show effect of an additional 100k dollars
plot(1, type = "n",
     xlim = range(age_plot),
     ylim = 1e5 * range(c(ind1_beta[,c("L95", "U95"), "wealth"], ind7_beta[,c("L95", "U95"), "wealth"])),
     xlab = "Age", ylab = expression(beta["wealth"]), main = "Effect of wealth")
polygon(x = c(age_plot, rev(age_plot)),
        y = c(ind1_beta[,"L95", "wealth"], rev(ind1_beta[,"U95", "wealth"])) * 1e5,
        col = rgb(0,0,1,1/5), border = NA)
polygon(x = c(age_plot, rev(age_plot)),
        y = c(ind7_beta[,"L95", "wealth"], rev(ind7_beta[,"U95", "wealth"])) * 1e5,
        col = rgb(1,0,0,1/5), border = NA)
lines(age_plot, 1e5*ind1_beta[,"MEAN", "wealth"], col = rgb(0,0,1,1/2))
lines(age_plot, 1e5*ind7_beta[,"MEAN", "wealth"], col = rgb(1,0,0,1/2))

# diabetes
plot(1, type = "n", 
     xlim = range(age_plot), 
     ylim = range(c(ind1_beta[,c("L95","U95"), "diabetes"], ind7_beta[,c("L95", "U95"), "diabetes"])), 
     xlab = "Age (months)", ylab = expression(beta["diabetes"]),
     main = "Effect of diabetes")
abline(h = 0, lty = 2)
polygon(c(age_plot, rev(age_plot)), c(ind1_beta[,"L95","diabetes"], rev(ind1_beta[,"U95", "diabetes"])),
        col = rgb(0,0,1,1/5), border = NA)
lines(age_plot, ind1_beta[,"MEAN", "diabetes"], col = 'blue')
polygon(c(age_plot, rev(age_plot)), c(ind7_beta[,"L95", "diabetes"], rev(ind7_beta[,"U95", "diabetes"])),
        col = rgb(1,0,0,1/5), border = NA)
lines(age_plot, ind7_beta[,"MEAN", "diabetes"], col = 'red')

dev.off()







# Look at cSES -- very small effects across the board
cses_ylim <- range(c(ind1_beta[, c("MEAN"),"cSES"], ind2_beta[,c("MEAN"), "cSES"]))
plot(1, type = "n", xlim = range(age_plot), ylim = cses_ylim,
     xlab = "Age", ylab = expression(beta["cses"]), main = "Effect of childhood SEP")
#polygon(c(age_plot, rev(age_plot)),
#        c(ind1_beta[,"L95", "cSES"], rev(ind1_beta[,"U95", "cSES"])),
#        col = rgb(0,0,1,1/5), border = NA)
#polygon(c(age_plot, rev(age_plot)),
#        c(ind2_beta[,"L95", "cSES"], rev(ind2_beta[,"U95", "cSES"])),
#        col = rgb(1,0,0,1/5), border = NA)
lines(age_plot, ind1_beta[,"MEAN", "cSES"], col = 'blue')
lines(age_plot, ind2_beta[,"MEAN", "cSES"], col = 'red')

# wealth things look like a constant shift -- so how the effect 


# education
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0), mfrow = c(2,2))

plot(1, type = "n", 
     xlim = range(age_plot), 
     ylim = range(c(ind1_beta[,c("L95", "U95"), "bmi"], ind7_beta[,c("L95", "U95"), "bmi"])),
     xlab = "Age", ylab = expression(beta["BMI"]), main = "Effect of BMI")
polygon(c(age_plot, rev(age_plot)),
        c(ind1_beta[,"L95","bmi"], rev(ind1_beta[,"U95", "bmi"])),
        col = rgb(0,0,1,1/5), border = NA)
polygon(c(age_plot, rev(age_plot)),
        c(ind7_beta[,"L95","bmi"], rev(ind7_beta[,"U95", "bmi"])),
        col = rgb(1,0,0,1/5), border = NA)
lines(age_plot, ind7_beta[,"MEAN", "bmi"], col = rgb(1,0,0,1/3))
lines(age_plot, ind1_beta[,"MEAN", "bmi"], col = rgb(0,0,1,1/3))

plot(1, type = "n", 
     xlim = range(age_plot), 
     ylim = range(c(ind1_beta[,c("L95", "U95"), "phys_activity"], ind7_beta[,c("L95", "U95"), "phys_activity"])),
     xlab = "Age", ylab = expression(beta["phys_act"]), main = "Effect of Physical Activity")
polygon(c(age_plot, rev(age_plot)),
        c(ind1_beta[,"L95","phys_activity"], rev(ind1_beta[,"U95", "phys_activity"])),
        col = rgb(0,0,1,1/5), border = NA)
polygon(c(age_plot, rev(age_plot)),
        c(ind7_beta[,"L95","phys_activity"], rev(ind7_beta[,"U95", "phys_activity"])),
        col = rgb(1,0,0,1/5), border = NA)
lines(age_plot, ind7_beta[,"MEAN", "phys_activity"], col = rgb(1,0,0,1/3))
lines(age_plot, ind1_beta[,"MEAN", "phys_activity"], col = rgb(0,0,1,1/3))


plot(1, type = "n", 
     xlim = range(age_plot), 
     ylim = range(c(ind1_beta[,c("L95", "U95"), "hi_bp"], ind7_beta[,c("L95", "U95"), "hi_bp"])),
     xlab = "Age", ylab = expression(beta["hi_bp"]), main = "Effect of High Blood Pressure")
polygon(c(age_plot, rev(age_plot)),
        c(ind1_beta[,"L95","hi_bp"], rev(ind1_beta[,"U95", "hi_bp"])),
        col = rgb(0,0,1,1/5), border = NA)
polygon(c(age_plot, rev(age_plot)),
        c(ind7_beta[,"L95","hi_bp"], rev(ind7_beta[,"U95", "hi_bp"])),
        col = rgb(1,0,0,1/5), border = NA)
lines(age_plot, ind7_beta[,"MEAN", "hi_bp"], col = rgb(1,0,0,1/3))
lines(age_plot, ind1_beta[,"MEAN", "hi_bp"], col = rgb(0,0,1,1/3))


plot(1, type = "n", 
     xlim = range(age_plot), 
     ylim = range(c(ind1_beta[,c("L95", "U95"), "heart_problems"], ind7_beta[,c("L95", "U95"), "heart_problems"])),
     xlab = "Age", ylab = expression(beta["heart_prob"]), main = "Effect of Heart Problems")
polygon(c(age_plot, rev(age_plot)),
        c(ind1_beta[,"L95","heart_problems"], rev(ind1_beta[,"U95", "heart_problems"])),
        col = rgb(0,0,1,1/5), border = NA)
polygon(c(age_plot, rev(age_plot)),
        c(ind7_beta[,"L95","heart_problems"], rev(ind7_beta[,"U95", "heart_problems"])),
        col = rgb(1,0,0,1/5), border = NA)
lines(age_plot, ind7_beta[,"MEAN", "heart_problems"], col = rgb(1,0,0,1/3))
lines(age_plot, ind1_beta[,"MEAN", "heart_problems"], col = rgb(0,0,1,1/3))


plot(1, type = "n", 
     xlim = range(age_plot), 
     ylim = range(c(ind1_beta[,c("L95", "U95"), "loneliness"], ind7_beta[,c("L95", "U95"), "loneliness"])),
     xlab = "Age", ylab = expression(beta["loneliness"]), main = "Effect of Loneliness")
polygon(c(age_plot, rev(age_plot)),
        c(ind1_beta[,"L95","loneliness"], rev(ind1_beta[,"U95", "loneliness"])),
        col = rgb(0,0,1,1/5), border = NA)
polygon(c(age_plot, rev(age_plot)),
        c(ind7_beta[,"L95","loneliness"], rev(ind7_beta[,"U95", "loneliness"])),
        col = rgb(1,0,0,1/5), border = NA)
lines(age_plot, ind7_beta[,"MEAN", "loneliness"], col = rgb(1,0,0,1/3))
lines(age_plot, ind1_beta[,"MEAN", "loneliness"], col = rgb(0,0,1,1/3))


plot(1, type = "n", 
     xlim = range(age_plot), 
     ylim = range(c(ind1_beta[,c("L95", "U95"), "stroke"], ind7_beta[,c("L95", "U95"), "stroke"])),
     xlab = "Age", ylab = expression(beta["stroke"]), main = "Effect of Stroke")
polygon(c(age_plot, rev(age_plot)),
        c(ind1_beta[,"L95","stroke"], rev(ind1_beta[,"U95", "stroke"])),
        col = rgb(0,0,1,1/5), border = NA)
polygon(c(age_plot, rev(age_plot)),
        c(ind7_beta[,"L95","stroke"], rev(ind7_beta[,"U95", "stroke"])),
        col = rgb(1,0,0,1/5), border = NA)
lines(age_plot, ind7_beta[,"MEAN", "stroke"], col = rgb(1,0,0,1/3))
lines(age_plot, ind1_beta[,"MEAN", "stroke"], col = rgb(0,0,1,1/3))

plot(1, type = "n", 
     xlim = range(age_plot), 
     ylim = range(c(ind1_beta[,c("L95", "U95"), "cesd"], ind7_beta[,c("L95", "U95"), "cesd"])),
     xlab = "Age", ylab = expression(beta["cesd"]), main = "Effect of CESD")
polygon(c(age_plot, rev(age_plot)),
        c(ind1_beta[,"L95","cesd"], rev(ind1_beta[,"U95", "cesd"])),
        col = rgb(0,0,1,1/5), border = NA)
polygon(c(age_plot, rev(age_plot)),
        c(ind7_beta[,"L95","cesd"], rev(ind7_beta[,"U95", "cesd"])),
        col = rgb(1,0,0,1/5), border = NA)
lines(age_plot, ind7_beta[,"MEAN", "cesd"], col = rgb(1,0,0,1/3))
lines(age_plot, ind1_beta[,"MEAN", "cesd"], col = rgb(0,0,1,1/3))
abline(h = 0, lty = 2)


plot(1, type = "n", 
     xlim = range(age_plot), 
     ylim = range(c(ind1_beta[,c("L95", "U95"), "int"], ind7_beta[,c("L95", "U95"), "int"])),
     xlab = "Age", ylab = expression(beta[0]), main = "Intercept")
polygon(c(age_plot, rev(age_plot)),
        c(ind1_beta[,"L95","int"], rev(ind1_beta[,"U95", "int"])),
        col = rgb(0,0,1,1/5), border = NA)
polygon(c(age_plot, rev(age_plot)),
        c(ind7_beta[,"L95","int"], rev(ind7_beta[,"U95", "int"])),
        col = rgb(1,0,0,1/5), border = NA)
lines(age_plot, ind7_beta[,"MEAN", "int"], col = rgb(1,0,0,1/3))
lines(age_plot, ind1_beta[,"MEAN", "int"], col = rgb(0,0,1,1/3))
abline(h = 0, lty = 2)


start_index_all <- 1 + c(0, cumsum(n_all)[-length(n_all)])
end_index_all <- cumsum(n_all)
n_all[1]
plot(1, type = "n", xlim = range(age_plot), ylim = range(hrs_vcbart_adapt_cs50$train$ystar[,c("L95", "U95")]),
     main = "Predicted Trajectories", xlab = "Age", ylab = "Cognitive Score")

i <- 3



points(Z_all[start_index_all[i]:end_index_all[i],1], Y_all[start_index_all[i]:end_index_all[i]], pch = 16)
lines(Z_all[start_index_all[i]:end_index_all[i],1], hrs_vcbart_adapt_cs50$train$ystar[start_index_all[i]:end_index_all[i],"MEAN"], col = 'red')
polygon(c(Z_all[start_index_all[i]:end_index_all[i]], rev(Z_all[start_index_all[i]:end_index_all[i]])),
        c(hrs_vcbart_adapt_cs50$train$ystar[start_index_all[i]:end_index_all[i],"L95"], 
          rev(hrs_vcbart_adapt_cs50$train$ystar[start_index_all[i]:end_index_all[i],"U95"])),
        col = rgb(1, 0, 0, 1/4), border = NA)

fit <- hrs_vcbart_adapt_cs50

training_loss <- rep(NA, times = length(n_all))
for(i in 1:length(n_all)){
  training_loss[i] <- 
    sqrt(mean( (Y_all[start_index_all[i]:end_index_all[i]] - fit$train$ystar[start_index_all[i]:end_index_all[i]])^2 ))
  
}


