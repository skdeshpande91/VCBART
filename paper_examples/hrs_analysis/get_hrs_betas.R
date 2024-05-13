load("../hrs_analysis/hrs_data1.RData")
load("hrs_all_trees.RData")
#########
y_mean <- mean(Y_all)
y_sd <- sd(Y_all)
x_mean <- c(0, apply(X1_all, MARGIN = 2, FUN = mean))
x_sd <- c(1, apply(X1_all, MARGIN = 2, FUN = sd))

fit <- list(trees = all_trees,
            y_mean = y_mean, y_sd = y_sd, 
            x_mean = x_mean, x_sd = x_sd)
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
beta_id2 <- VCBART::predict_betas(fit, Z_cont = Z_cont2, Z_cat = Z_cat2)

beta_sum_id1 <- VCBART::summarize_beta(beta_id1)
beta_sum_id2 <- VCBART::summarize_beta(beta_id2)

diff <- beta_id1 - beta_id2
diff_sum <- VCBART::summarize_beta(diff)

age_quants <- quantile(Z_cont_all[,"AGE"], probs = c(0.025, 0.975))
save(beta_sum_id1, Z_cont1, Z_cat1,
     beta_sum_id2, Z_cont2, Z_cat2, 
     diff_sum,
     age_seq, age_quants,
     file = "~/Documents/Research/vc_bart/writing/figure_scripts/figure4/hrs_pred_beta.RData")
