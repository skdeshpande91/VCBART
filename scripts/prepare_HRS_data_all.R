load("data/HRS_data_1.RData")
X_all <- rbind(X_train, X_test)
Z_all <- rbind(Z_train, Z_test)
Y_all <- c(Y_train, Y_test)
n_all <- c(n_vec_train, n_vec_test)


#############
# Make up a few individuals
# Male, southern born, not foreign born, working full time at 50
# childhood health = 5
##############
X_plot <- matrix(nrow = 8 * length(cutpoints[[1]]), ncol = p, dimnames = list(c(), colnames(X_train)))
Z_plot <- matrix(nrow = 8 * length(cutpoints[[1]]), ncol = R, dimnames = list(c(), colnames(Z_train)))
n_plot <- rep(length(cutpoints[[1]]), times = 8)
start_index_plot <- 1 + c(0, cumsum(n_plot)[-8])
end_index_plot <- cumsum(n_plot)

binary_covs <- c("phys_activity", "diabetes", "hi_bp", "heart_problems", "stroke", "loneliness")
cont_covs <- c("education", "cesd", "bmi", "cSES", "wealth")

X_plot[,cont_covs] <- colMeans(X_all[,cont_covs])
for(k in binary_covs){
  if(mean(X_all[,k]) > 0.5) X_plot[,k] <- 1
  else X_plot[,k] <- 0
}

# Now do Z_plot
Z_plot[,"race2"] <- 0
Z_plot[,"race3"] <- 0
Z_plot[,"childhood_health"] <- 5
Z_plot[,"gender"] <- 1
Z_plot[,"souther_born"] <- 1
Z_plot[,"foreign_born"] <- 0
Z_plot[,"veteran"] <- 0
Z_plot[,"work_full_time"] <- 1
Z_plot[,"work_part_time"] <- 0
Z_plot[,"unemployed"] <- 0
Z_plot[,"food_insecurity"] <- 0

# Individual 1:
# race0 = 1, race1 = 0, married = 0, food stamp = 0
Z_plot[start_index_plot[1]:end_index_plot[1],"age"] <- cutpoints[[1]]
Z_plot[start_index_plot[1]:end_index_plot[1],"race0"] <- 1
Z_plot[start_index_plot[1]:end_index_plot[1],"race1"] <- 0
Z_plot[start_index_plot[1]:end_index_plot[1],"married_partnered"] <- 0
Z_plot[start_index_plot[1]:end_index_plot[1], "food_stamp"] <- 0
Z_plot[start_index_plot[1]:end_index_plot[1], "food_insecurity"] <- 0

# Individiual 2
# race0 = 1, race1 = 0, married = 0, food stamp = 1
Z_plot[start_index_plot[2]:end_index_plot[2],"age"] <- cutpoints[[1]]
Z_plot[start_index_plot[2]:end_index_plot[2],"race0"] <- 1
Z_plot[start_index_plot[2]:end_index_plot[2],"race1"] <- 0
Z_plot[start_index_plot[2]:end_index_plot[2],"married_partnered"] <- 0
Z_plot[start_index_plot[2]:end_index_plot[2], "food_stamp"] <- 1

# Individiual 3
# race0 = 1, race1 = 0, married = 1, food stamp = 0
Z_plot[start_index_plot[3]:end_index_plot[3],"age"] <- cutpoints[[1]]
Z_plot[start_index_plot[3]:end_index_plot[3],"race0"] <- 1
Z_plot[start_index_plot[3]:end_index_plot[3],"race1"] <- 0
Z_plot[start_index_plot[3]:end_index_plot[3],"married_partnered"] <- 1
Z_plot[start_index_plot[3]:end_index_plot[3], "food_stamp"] <- 0

# Individiual 4
# race0 = 1, race1 = 0, married = 1, food stamp = 1
Z_plot[start_index_plot[4]:end_index_plot[4],"age"] <- cutpoints[[1]]
Z_plot[start_index_plot[4]:end_index_plot[4],"race0"] <- 1
Z_plot[start_index_plot[4]:end_index_plot[4],"race1"] <- 0
Z_plot[start_index_plot[4]:end_index_plot[4],"married_partnered"] <- 1
Z_plot[start_index_plot[4]:end_index_plot[4], "food_stamp"] <- 1



# Individual 5:
# race0 = 1, race1 = 0, married = 0, food stamp = 0
Z_plot[start_index_plot[5]:end_index_plot[5],"age"] <- cutpoints[[1]]
Z_plot[start_index_plot[5]:end_index_plot[5],"race0"] <- 0
Z_plot[start_index_plot[5]:end_index_plot[5],"race1"] <- 1
Z_plot[start_index_plot[5]:end_index_plot[5],"married_partnered"] <- 0
Z_plot[start_index_plot[5]:end_index_plot[5], "food_stamp"] <- 0
Z_plot[start_index_plot[5]:end_index_plot[5], "food_insecurity"] <- 0

# Individiual 6
# race0 = 1, race1 = 0, married = 0, food stamp = 1
Z_plot[start_index_plot[6]:end_index_plot[6],"age"] <- cutpoints[[1]]
Z_plot[start_index_plot[6]:end_index_plot[6],"race0"] <- 0
Z_plot[start_index_plot[6]:end_index_plot[6],"race1"] <- 1
Z_plot[start_index_plot[6]:end_index_plot[6],"married_partnered"] <- 0
Z_plot[start_index_plot[6]:end_index_plot[6], "food_stamp"] <- 1

# Individiual 7
# race0 = 0, race1 = 1, married = 1, food stamp = 0
Z_plot[start_index_plot[7]:end_index_plot[7],"age"] <- cutpoints[[1]]
Z_plot[start_index_plot[7]:end_index_plot[7],"race0"] <- 0
Z_plot[start_index_plot[7]:end_index_plot[7],"race1"] <- 1
Z_plot[start_index_plot[7]:end_index_plot[7],"married_partnered"] <- 1
Z_plot[start_index_plot[7]:end_index_plot[7], "food_stamp"] <- 0

# Individiual 8
# race0 = 0, race1 = 1, married = 1, food stamp = 1
Z_plot[start_index_plot[8]:end_index_plot[8],"age"] <- cutpoints[[1]]
Z_plot[start_index_plot[8]:end_index_plot[8],"race0"] <- 0
Z_plot[start_index_plot[8]:end_index_plot[8],"race1"] <- 1
Z_plot[start_index_plot[8]:end_index_plot[8],"married_partnered"] <- 1
Z_plot[start_index_plot[8]:end_index_plot[8], "food_stamp"] <- 1


save(X_all, Z_all, Y_all, n_all, X_plot, Z_plot, n_plot, cutpoints, 
     file = "~/Documents/Research/VCBART/data/HRS_all.RData")
