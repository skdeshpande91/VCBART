n_all <- n_train + n_test
ni_all <- rep(4, times = n_all) # 4 observations per subject
subj_id_all <- rep(1:n_all, each = 4) # give every subject an id number
N_all <- sum(ni_all) # total number of observations

p <- 5 # number of covariates
R_cont <- 20 # number of continuous modifiers
R_cat <- 0 # number of categorical modifiers
R <- R_cont + R_cat

Sigma_X <-  (0.5)^(abs(outer(1:p, 1:p, FUN = "-"))) # covariates are all correlated

#############
# Generate X, Z, etc.
sigma <- 1
X_all <- MASS::mvrnorm(N_all, mu = rep(0, times = p), Sigma = Sigma_X)
Z_cont_all <- matrix(runif(N_all * R_cont, min = -1, max = 1), nrow = N_all, ncol = R_cont)
beta0_all <- beta0_true(Z_cont_all)
beta1_all <- beta1_true(Z_cont_all)
beta2_all <- beta2_true(Z_cont_all)
beta3_all <- beta3_true(Z_cont_all)
beta4_all <- beta4_true(Z_cont_all)
beta5_all <- beta5_true(Z_cont_all)

beta_all <- cbind(beta0_all, beta1_all, beta2_all, beta3_all, beta4_all, beta5_all)

mu_all <- beta0_all + X_all[,1] * beta1_all + X_all[,2] * beta2_all +
  X_all[,3] * beta3_all + X_all[,4] * beta4_all + X_all[,5] * beta5_all
Y_all <- mu_all + sigma * rnorm(n = N_all, mean = 0, sd = 1)

# For the purposes of testing other methods
# we need a data frame containing all of the covariates and modifiers
cov_all <- data.frame(X_all)
colnames(cov_all) <- paste0("X",1:p)
mod_all <- data.frame(Z_cont_all)
colnames(mod_all) <- paste0("Z",1:R_cont)



# Create a training/testing split
# since we randomly generated X and Z, we can just use first n_train subjects as train
train_subjects <- 1:n_train
test_subjects <- (n_train+1):(n_train + n_test)

test_index <- which(subj_id_all %in% test_subjects)
train_index <- which(subj_id_all %in% train_subjects)

ni_train <- ni_all[train_subjects]
ni_test <- ni_all[test_subjects]

N_train <- sum(ni_train)
N_test <- sum(ni_test)

X_train <- X_all[train_index,]
X_test <- X_all[test_index,]

Z_cont_train <- Z_cont_all[train_index,]
Z_cont_test <- Z_cont_all[test_index,]
subj_id_train <- rep(1:n_train, times = ni_train)

Y_train <- Y_all[train_index]
Y_test <- Y_all[test_index]

beta_train <- beta_all[train_index,]
beta_test <- beta_all[test_index,]

cov_train <- cov_all[train_index,]
cov_test <- cov_all[test_index,]

mod_train <- mod_all[train_index,]
mod_test <- mod_all[test_index,]



###########
# Set a few arguments for VCBART
cutpoints_list <- NULL
unif_cuts <- rep(TRUE, times = R_cont)
cat_levels_list <- NULL

#M <- 200
#tau <- rep(0.5/sqrt(M), times = p+1)