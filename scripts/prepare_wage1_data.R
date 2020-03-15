# Prepare the wage1 dataset
library(np)
data(wage1)

N <- nrow(wage1)

cov_names <- c("educ")
mod_names <- c("female", "nonwhite", "married", "tenure", "exper")

p <- length(cov_names)
R <- length(mod_names)
cov_all <- as.data.frame( wage1[,cov_names])
colnames(cov_all) <- cov_names
mod_all <- wage1[,mod_names]


X_all <- as.matrix(cov_all, nrow = N, ncol = p)
Z_all <- mod_all
Z_all[,"female"] <- as.numeric(Z_all[,"female"]) - 1 # 0 for female, 1 for male
Z_all[,"nonwhite"] <- as.numeric(Z_all[,"nonwhite"]) - 1 # 0 for nonwhite, 1 for white
Z_all[,"married"] <- as.numeric(Z_all[,"married"]) - 1 # 0 for married, 1 for not married
Z_all <- as.matrix(Z_all, nrow  = N, ncol = R)
Y_all <- wage1[,"lwage"]

cutpoints <- list()
cutpoints[[1]] <- c(0,1)
cutpoints[[2]] <- c(0,1)
cutpoints[[3]] <- c(0,1)
cutpoints[[4]] <- seq(min(Z_all[,"tenure"]), max(Z_all[,"tenure"]), length = 10000)
cutpoints[[5]] <- seq(min(Z_all[,"exper"]), max(Z_all[,"exper"]), length = 10000)

save(N, p, R, Y_all, X_all, Z_all, mod_all, cov_all, cutpoints, file = "data/wage1.RData")
