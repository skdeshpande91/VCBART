library(dbarts)
library(CBPS)
load("data/HRS_causal_df.RData")

n <- nrow(data_df)
data_df[,"tmpTreat"] <- rep(NA, times = n)
data_df[which(data_df[,"PHYS_INACT"] == 0 & data_df[,"DEPRESS"] == 0), "tmpTreat"] <- 0 # control -- not inactive, not depressed
data_df[which(data_df[,"PHYS_INACT"] == 1 & data_df[,"DEPRESS"] == 0), "tmpTreat"] <- 1 # inactive & not depressed
data_df[which(data_df[,"PHYS_INACT"] == 0 & data_df[,"DEPRESS"] == 1), "tmpTreat"] <- 2 # active & depressed
data_df[which(data_df[,"PHYS_INACT"] == 1 & data_df[,"DEPRESS"] == 1), "tmpTreat"] <- 3 # inactive & depressed


ps_df <- data_df[,!colnames(data_df) %in% c("Y", "PHYS_INACT", "DEPRESS")]
ps_df[,"tmpTreat"] <- factor(data_df[,"tmpTreat"], levels = 0:3, labels = 0:3)


cbps_fit <- CBPS(tmpTreat ~ ., data = ps_df, ATT = 0, standardize = FALSE, iterations = 5000)
save(cbps_fit, file = "results/HRS_causal_cbps_fit.RData")
prop_score <- cbps_fit$fitted.values
colnames(prop_score) <- paste0("prop", 0:3)
weights <- cbps_fit$weights

####################
# Set up what we need to look at the balance table
tmp_df <- data_df[,!colnames(data_df) %in% c("Y", "tmpTreat", "PHYS_INACT", "DEPRESS")]
tmp_df[,"CHLD_HLTH"] <- as.numeric(tmp_df[,"CHLD_HLTH"])
Z_all <- makeModelMatrixFromDataFrame(tmp_df)

Z_all <- cbind(Z_all, prop_score[,paste0("prop", 1:3)])
ix0 <- which(data_df[,"tmpTreat"] == 0)
ix1 <- which(data_df[,"tmpTreat"] == 1)
ix2 <- which(data_df[,"tmpTreat"] == 2)
ix3 <- which(data_df[,"tmpTreat"] == 3)

Y_all <- data_df[,"Y"]
X_all <- matrix(0, nrow = n, ncol = 3)
X_all[ix1, 1] <- 1
X_all[ix2, 2] <- 1
X_all[ix3, 3] <- 1

R <- ncol(Z_all)
cutpoints <- list()
for(r in 1:R){
  if(colnames(Z_all)[r] == "AGE") cutpoints[[r]] <- 720:1020
  else if (!colnames(Z_all)[r] %in% c("cSEP", "prop1", "prop2", "prop3")) cutpoints[[r]] <- sort(unique(Z_all[,r]))
  else cutpoints[[r]] <- seq(floor(min(Z_all[,r])), ceiling(max(Z_all[,r])), length = 1000)
}
n_all <- length(Y_all)
save(X_all, Y_all, Z_all, R, n_all, cutpoints, file = "data/HRS_causal_vcbart_data.RData")
############################
# Compute pairwise differences in unweighted and weighted means
# now we want chld health treated as a factor!

tmp_df <- data_df[,!colnames(data_df) %in% c("Y", "PHYS_INACT", "DEPRESS", "tmpTreat")]

tmp_df <- data_df[,c("AGE", "GENDER", "RACE", "cSEP", "EDUC", "CHLD_HLTH", "SOUTHERN", "FOREIGN", "SMOKE")]

tmp_Z <- makeModelMatrixFromDataFrame(tmp_df)
colnames(tmp_Z) <- c("Age", "Female","NH White", "NH Black", "Hispanic", "NH Other", "cSEP", "Educ", 
                     "Poor chld health", "Fair chld health", "Good chld health", "Very good chld health", "Excellent chld health",
                     "Born in Southern US", "Born outside US", "Smoked")
p <- ncol(tmp_Z)
std_diff_before <- matrix(nrow = p, ncol = 6, dimnames = list(colnames(tmp_Z), c("1-0", "2-0", "3-0", "2-1", "3-1", "3-2")))
std_diff_after <- matrix(nrow = p, ncol = 6, dimnames = list(colnames(tmp_Z), c("1-0", "2-0", "3-0", "2-1", "3-1", "3-2")))

n0 <- length(ix0)
n1 <- length(ix1)
n2 <- length(ix2)
n3 <- length(ix3)

for(j in 1:p){
  z_name <- colnames(tmp_Z)[j]
  z <- tmp_Z[,j]
  
  var0 <- var(z[ix0])
  var1 <- var(z[ix1])
  var2 <- var(z[ix2])
  var3 <- var(z[ix3])
  
  pooled_sd <- sqrt( ((n0 - 1)*var0 + (n1 - 1) * var1 + (n2 - 1) * var2 + (n3 - 1) * var3) / (n0 + n1 + n2 + n3 - 4))
  
  before_mean0 <- mean(z[ix0])
  before_mean1 <- mean(z[ix1])
  before_mean2 <- mean(z[ix2])
  before_mean3 <- mean(z[ix3])
  
  after_mean0 <- weighted.mean(z[ix0], w = weights[ix0])
  after_mean1 <- weighted.mean(z[ix1], w = weights[ix1])
  after_mean2 <- weighted.mean(z[ix2], w = weights[ix2])
  after_mean3 <- weighted.mean(z[ix3], w = weights[ix3])
  
  
  
  std_diff_before[z_name, "1-0"] <- (before_mean1 - before_mean0)/pooled_sd
  std_diff_before[z_name, "2-0"] <- (before_mean2 - before_mean0)/pooled_sd
  std_diff_before[z_name, "3-0"] <- (before_mean3 - before_mean0)/pooled_sd
  
  std_diff_before[z_name, "2-1"] <- (before_mean2 - before_mean1)/pooled_sd
  std_diff_before[z_name, "3-1"] <- (before_mean3 - before_mean1)/pooled_sd
  std_diff_before[z_name, "3-2"] <- (before_mean3 - before_mean2)/pooled_sd
  
  std_diff_after[z_name, "1-0"] <- (after_mean1 - after_mean0)/pooled_sd
  std_diff_after[z_name, "2-0"] <- (after_mean2 - after_mean0)/pooled_sd
  std_diff_after[z_name, "3-0"] <- (after_mean3 - after_mean0)/pooled_sd
  
  std_diff_after[z_name, "2-1"] <- (after_mean2 - after_mean1)/pooled_sd
  std_diff_after[z_name, "3-1"] <- (after_mean3 - after_mean1)/pooled_sd
  std_diff_after[z_name, "3-2"] <- (after_mean3 - after_mean2)/pooled_sd
}
save(tmp_Z, std_diff_before, std_diff_after, file = "results/HRS_causal_balance.RData")