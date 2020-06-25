load("data/HRS/HRS_all.RData")
load("results/vcbart_hrs_all.RData")
age_plot <- Z_plot[1:544,"age"]

for(id in 1:8){
  
  tmp_beta <- hrs_vcbart_adapt_cs50$test$beta[544*(id-1) + 1:544, ,]
  dimnames(tmp_beta)[[3]] <- c("int", colnames(X_all))
  
  assign(paste0("ind",id, "_beta"),
         tmp_beta)
}







# plot the effect of diabetes for individual 1 and 8
# diabetes effect modified by 
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0))
plot(1, type = "n", xlim = range(age_plot), ylim = c(-6,6), 
     xlab = "Age (months)", ylab = expression(beta["diab"]),
     main = "Effect of diabetes")
abline(h = 0, lty = 2)
polygon(c(age_plot, rev(age_plot)), c(ind1_beta[,"L95","diabetes"], rev(ind1_beta[,"U95", "diabetes"])),
        col = rgb(0,0,1,1/3), border = NA)
lines(age_plot, ind1_beta[,"MEAN", "diabetes"], col = 'blue')
polygon(c(age_plot, rev(age_plot)), c(ind3_beta[,"L95", "diabetes"], rev(ind3_beta[,"U95", "diabetes"])),
        col = rgb(1,0,0,1/3), border = NA)
lines(age_plot, ind3_beta[,"MEAN", "diabetes"], col = 'red')

###### Effect of education
# varies wrt race and food stamp

plot(1, type = "n", xlim = range(age_plot), ylim = c(0.25,1.25), 
     xlab = "Age", ylab = expression(beta["educ"]), main = "Effect of education")
polygon(c(age_plot, rev(age_plot)),
        c(ind1_beta[,"L95","education"], rev(ind1_beta[,"U95", "education"])),
        col = rgb(0,0,1,1/3), border = NA)
polygon(c(age_plot, rev(age_plot)),
        c(ind2_beta[,"L95","education"], rev(ind2_beta[,"U95", "education"])),
        col = rgb(1,0,0,1/3), border = NA)
lines(age_plot, ind2_beta[,"MEAN", "education"], col = 'red')
lines(age_plot, ind1_beta[,"MEAN", "education"], col = 'blue')



##########
# Look at cSES -- there is some variation but it's all's non-significant and very close to 0

