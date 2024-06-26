load("philly_crime_data.RData")
source("vcbart_ind_wrapper.R")
p <- 4


X_train <- X_all[,1:p]

ni_train <- ni_all
subj_id_train <- rep(1:n, times = ni_train) 

fit <- 
  vcbart_ind_wrapper(Y_train = Y_all,
                     subj_id_train = subj_id_train,
                     ni_train = ni_train,
                     X_train = X_train,
                     Z_cat_train = Z_cat_all,
                     unif_cuts = unif_cuts,
                     cutpoints_list = cutpoints_list,
                     cat_levels_list = cat_levels_list,
                     edge_mat_list = edge_mat_list,
                     graph_split = graph_split,
                     sparse = FALSE, M = 50)

beta0_hat <- unique(fit$train$beta[,,1])
beta1_hat <- unique(fit$train$beta[,,2])
beta2_hat <- unique(fit$train$beta[,,3])
beta3_hat <- unique(fit$train$beta[,,4])
beta4_hat <- unique(fit$train$beta[,,5])

sigma_chain <- 
  coda::mcmc.list(
    coda::mcmc(fit[[4]][,1]),
    coda::mcmc(fit[[4]][,2]),
    coda::mcmc(fit[[4]][,3]),
    coda::mcmc(fit[[4]][,4]))

coda::gelman.diag(sigma_chain)

save(beta0_hat, beta1_hat, beta2_hat, 
     beta3_hat, beta4_hat, sigma_chain,
     file = "philly_betas.RData")



###############
library(ggplot2)
#library(ggpattern)
load("../figure5/tracts.RData")
tract_fort <- fortify(tracts)
id <- as.character(unique(tract_fort$id))

tmp_data <- data.frame(id = id, 
                       tract = rownames(A_tract),
                       beta0_mean = beta0_hat[,1],
                       beta1_mean = beta1_hat[,1],
                       beta2_mean = beta2_hat[,1],
                       beta3_mean = beta3_hat[,1],
                       beta4_mean = beta4_hat[,1],
                       beta0_sig = 1*(beta0_hat[,2] * beta0_hat[,3] >= 0),
                       beta1_sig = 1*(beta1_hat[,2] * beta1_hat[,3] >= 0),
                       beta2_sig = 1*(beta2_hat[,2] * beta2_hat[,3] >= 0),
                       beta3_sig = 1*(beta3_hat[,2] * beta3_hat[,3] >= 0),
                       beta4_sig = 1*(beta4_hat[,2] * beta4_hat[,3] >= 0))

plot_data <- dplyr::inner_join(tract_fort, tmp_data, by = "id")
######################
# Limits
beta_range <- range(c(beta1_hat, beta2_hat, beta3_hat, beta4_hat))

####### Intercept plot
beta0_plot <- 
  ggplot(data = plot_data) + 
  geom_polygon(aes(x = long, y = lat, group = group, fill = beta0_mean),
               color = alpha("black", 0.25), alpha = 1, show.legend = TRUE) + 
  scale_fill_distiller(type = "div", palette = "PRGn") + 
  coord_map() + 
  theme(panel.background = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(),
        axis.title = element_blank(), 
        legend.title = element_blank())

# First-order term plot
beta1_plot <- 
  ggplot(data = plot_data) + 
  geom_polygon(aes(x = long, y = lat, group = group, fill = beta1_mean),
               color = alpha("black", 0.25), alpha = 1, show.legend = TRUE) + 
  scale_fill_distiller(type = "div", palette = "RdBu", 
                       limits = beta_range) + 
  coord_map() + 
  theme(panel.background = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(),
        axis.title = element_blank(), 
        legend.title = element_blank())

beta2_plot <- 
  ggplot(data = plot_data) + 
  geom_polygon(aes(x = long, y = lat, group = group, fill = beta2_mean),
               color = alpha("black", 0.25), alpha = 1, show.legend = TRUE) + 
  coord_map() + 
  scale_fill_distiller(type = "div", palette = "RdBu", 
                       limits = beta_range) +
  theme(panel.background = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(),
        axis.title = element_blank(), 
        legend.title = element_blank())


beta3_plot <- 
  ggplot(data = plot_data) + 
  geom_polygon(aes(x = long, y = lat, group = group, fill = beta3_mean),
               color = alpha("black", 0.25), alpha = 1, show.legend = TRUE) + 
  coord_map() + 
  scale_fill_distiller(type = "div", palette = "RdBu", 
                       limits = beta_range) +
  theme(panel.background = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(),
        axis.title = element_blank(), 
        legend.title = element_blank())

beta4_plot <- 
  ggplot(data = plot_data) + 
  geom_polygon(aes(x = long, y = lat, group = group, fill = beta4_mean),
               color = alpha("black", 0.25), alpha = 1, show.legend = TRUE) + 
  coord_map() + 
  scale_fill_distiller(type = "div", palette = "RdBu", 
                       limits = beta_range) +
  theme(panel.background = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(),
        axis.title = element_blank(), 
        legend.title = element_blank())
