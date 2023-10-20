load("hrs_data1.RData")
load("chains/hrs_all_trees.RData")


#########

y_mean <- mean(Y_all)
y_sd <- sd(Y_all)
x_mean <- c(0, apply(X1_all, MARGIN = 2, FUN = mean))
x_sd <- c(1, apply(X1_all, MARGIN = 2, FUN = sd))

fit <- list(trees = all_trees,
            y_mean = y_mean, y_sd = y_sd, 
            x_mean = x_mean, x_sd = x_sd)
#############

age_seq <- cutpoints_list[["AGE"]]
R_cont <- ncol(Z_cont_all)
R_cat <- ncol(Z_cat_all)
cont_names <- colnames(Z_cont_all)
cat_names <- colnames(Z_cat_all)

########################
# Individual 1: Gender = 0 (female), Hispanic = 0, Race = 2 (white), Birthplace = 0 (E/N Central)
Z_cont1 <- matrix(nrow = length(age_seq), ncol = R_cont, dimnames = list(c(), cont_names))
Z_cont1[,"AGE"] <- age_seq
Z_cont1[,"GENDER_num"] <- 0
Z_cont1[,"HISPANIC_num"] <- 0

Z_cat1 <- matrix(nrow = length(age_seq), ncol = R_cat, dimnames = list(c(), cat_names))
Z_cat1[,"RACE_num"] <- 2
Z_cat1[,"BIRTH_PLACE_num"] <- 0
########
# Individual 2: Gender = 0 (female), Hispanic = 0, Race = 0 (Black), Birthplace = 0 (EN Central)
Z_cont2 <- matrix(nrow = length(age_seq), ncol = R_cont, dimnames = list(c(), cont_names))
Z_cont2[,"AGE"] <- age_seq
Z_cont2[,"GENDER_num"] <- 0
Z_cont2[,"HISPANIC_num"] <- 0

Z_cat2 <- matrix(nrow = length(age_seq), ncol = R_cat, dimnames = list(c(), cat_names))
Z_cat2[,"RACE_num"] <- 0
Z_cat2[,"BIRTH_PLACE_num"] <- 0

##############
beta_id1 <- VCBART::predict_betas(fit, Z_cont = Z_cont1, Z_cat = Z_cat1)
beta_sum_id1 <- VCBART::summarize_beta(beta_id1)

beta_id2 <- VCBART::predict_betas(fit, Z_cont = Z_cont2, Z_cat = Z_cat2)
beta_sum_id2 <- VCBART::summarize_beta(beta_id2)
age_quants <- quantile(Z_cont_all[,"AGE"], probs = c(0.025, 0.975))

save(beta_sum_id1, Z_cont1, Z_cat1,
     beta_sum_id2, Z_cont2, Z_cat2, 
     age_seq, age_quants,
     file = "~/Documents/Research/vc_bart/writing/figure_scripts/figure7/hrs_pred_beta.RData")

#########################
# Intercept (a baseline cognitive score)
#########################
xlim <- quantile(Z_cont_all[,"AGE"], probs = c(0.025, 0.975))
ylim0 <- range(c(beta_sum_id1[,,1], beta_sum_id2[,,1]))

pdf("~/Documents/Research/vc_bart/")
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

abline(v = xlim, col = my_colors[8], lty = 2, lwd = 1)

legend("bottomleft", legend = c("Black", "White"), col = my_colors[c(3,1)],
       lwd = 2.5, lty = 1, bty = "n")
########################

# Effect of cSEP

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

legend("bottomleft", legend = c("Black", "White"), col = my_colors[c(3,1)],
       lwd = 2.5, lty = 1, bty = "n")

##############
# HS completition
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

legend("bottomleft", legend = c("Black", "White"), col = my_colors[c(3,1)],
       lwd = 2.5, lty = 1, bty = "n")
############


##############
# Labor force participation

ylim5 <- range(c(beta_sum_id1[,,6], beta_sum_id2[,,6]))

plot(1, type = "n", 
     xlim = xlim, ylim = ylim5,
     xlab = expression("Age (months)"), 
     ylab = expression ("Cognitive score"),
     main = expression("Labor force participation"))
polygon(x = c(age_seq, rev(age_seq)),
        y = c(beta_sum_id1[,"L95",6], rev(beta_sum_id1[,"U95",6])),
        border = NA, col = adjustcolor(my_colors[1], alpha.f = 0.35))

polygon(x = c(age_seq, rev(age_seq)),
        y = c(beta_sum_id2[,"L95",6], rev(beta_sum_id2[,"U95",6])),
        border = NA, col = adjustcolor(my_colors[3], alpha.f = 0.35))

lines(x = age_seq, y = beta_sum_id1[,"MEAN",6], col = my_colors[1], lwd = 2.5)
lines(x = age_seq, y = beta_sum_id2[,"MEAN",6], col = my_colors[3], lwd = 2.5)

legend("bottomleft", legend = c("Black", "White"), col = my_colors[c(3,1)],
       lwd = 2.5, lty = 1, bty = "n")

