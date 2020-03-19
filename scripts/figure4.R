load("data/philly_cc_shp.Rdata")
source("scripts/maps_function.R")

library(ggplot2)
library(dplyr)
library(gridExtra)
library(scales)
library(RColorBrewer)
philly_cc_fort <- fortify(philly_cc_poly)

load("results/philly_cs_0.RData")
load("data/philly_monthly_cc.RData")
n <- 255
p <- ncol(X_train) + 1
# Collect the means 
beta_mean_train <- array(dim = c(n, 42, p), dimnames = list(c(), c(), c("intercept", colnames(X_train))))
beta_mean_test <- array(dim = c(n, 6, p), dimnames = list(c(), c(), c("intercept", colnames(X_train))))

beta_L95_train <- array(dim = c(n, 42, p), dimnames = list(c(), c(), c("intercept", colnames(X_train))))
beta_L95_test <- array(dim = c(n, 6, p), dimnames = list(c(), c(), c("intercept", colnames(X_train))))

beta_U95_train <- array(dim = c(n, 42, p), dimnames = list(c(), c(), c("intercept", colnames(X_train))))
beta_U95_test <- array(dim = c(n, 6, p), dimnames = list(c(), c(), c("intercept", colnames(X_train))))

for(i in 1:n){
  beta_mean_train[i,,] <- cs_rho0_sum$train$beta[start_index_train[i]:(start_index_train[i] + n_vec_train[i] - 1),"MEAN",]
  beta_mean_test[i,,] <- cs_rho0_sum$test$beta[start_index_test[i]:(start_index_test[i] + n_vec_test[i] - 1),"MEAN",]
  
  beta_L95_train[i,,] <- cs_rho0_sum$train$beta[start_index_train[i]:(start_index_train[i] + n_vec_train[i] - 1),"L95",]
  beta_L95_test[i,,] <- cs_rho0_sum$test$beta[start_index_test[i]:(start_index_test[i] + n_vec_test[i] - 1),"L95",]
  
  beta_U95_train[i,,] <- cs_rho0_sum$train$beta[start_index_train[i]:(start_index_train[i] + n_vec_train[i] - 1),"U95",]
  beta_U95_test[i,,] <- cs_rho0_sum$test$beta[start_index_test[i]:(start_index_test[i] + n_vec_test[i] - 1),"U95",]
  
}

#######
# Plot effect of log income at time t = 1 (January 2015) and t = 19 (July 2016)
beta1_raw <- beta_mean_train[,,"log.income13"]
group_ind <- as.numeric(unique(philly_cc_fort$id))

beta1_t1_t19_t37_raw <- beta1_raw[,c(1,19,37)]
beta1_t1_t19_t37_rescaled <- beta1_t1_t19_t37_raw

beta1_t1_t19_t37_rescaled[beta1_t1_t19_t37_raw <= 0] <- rescale(beta1_t1_t19_t37_raw[beta1_t1_t19_t37_raw <= 0], to = c(-1,0))
beta1_t1_t19_t37_rescaled[beta1_t1_t19_t37_raw > 0] <- rescale(beta1_t1_t19_t37_raw[beta1_t1_t19_t37_raw > 0], to = c(0,1))

tmp_data <- data.frame(id = as.character(group_ind),
                       beta1_t1 = beta1_t1_t19_t37_rescaled[,1],
                       beta1_t19 = beta1_t1_t19_t37_rescaled[,2],
                       beta1_t37 = beta1_t1_t19_t37_rescaled[,3])
plotData <- inner_join(philly_cc_fort, tmp_data, by = "id")



beta1_t1_plot <- ggplot() +
  geom_polygon(data = plotData, aes(x = long, y = lat, group = group, fill = beta1_t1),
               color = alpha("black", 0.15), alpha = 1, show.legend = FALSE) +
  scale_fill_distiller(type = "div", palette = "RdBu", limits = c(-1,1)) + 
  theme(panel.background = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), 
        axis.title = element_blank(), plot.title = element_text(hjust = 0.5))
beta1_t19_plot <- ggplot() +
  geom_polygon(data = plotData, aes(x = long, y = lat, group = group, fill = beta1_t19),
               color = alpha("black", 0.15), alpha = 1, show.legend = FALSE) +
  scale_fill_distiller(type = "div", palette = "RdBu", limits = c(-1,1)) + 
  theme(panel.background = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), 
        axis.title = element_blank(), plot.title = element_text(hjust = 0.5))
beta1_t37_plot <- ggplot() +
  geom_polygon(data = plotData, aes(x = long, y = lat, group = group, fill = beta1_t37),
               color = alpha("black", 0.15), alpha = 1, show.legend = FALSE) +
  scale_fill_distiller(type = "div", palette = "RdBu", limits = c(-1,1)) + 
  theme(panel.background = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), 
        axis.title = element_blank(), plot.title = element_text(hjust = 0.5))


tmp_data <- data.frame(id = as.character(group_ind),
                       beta1_t1 = beta1_t1_rescaled)
plotData <- inner_join(philly_cc_fort, tmp_data, by = "id")

ggsave(plot = beta1_t1_plot, file = "figures/philly_beta1_t1.png", 
       width = 4, height = 4, units = "in", scale = 3)
ggsave(plot = beta1_t19_plot, file = "figures/philly_beta1_t19.png", 
       width = 4, height = 4, units = "in", scale = 3)
ggsave(plot = beta1_t37_plot, file = "figures/philly_beta1_t37.png", 
       width = 4, height = 4, units = "in", scale = 3)


ggsave("figures/philly_beta1_t1_t19_t37.png", width = 6, height = 2, units = "in", scale = 3,
       plot = arrangeGrob(grobs = list(beta1_t1_plot, beta1_t19_plot, beta1_t37_plot),nrow = 1, ncol = 3))


# Make the legend
col_list <- rev(brewer.pal(n = 5, name = "RdBu"))
lower_range <- c(min(beta1_t1_t19_t37_raw[beta1_t1_t19_t37_raw <= 0]), 0)
upper_range <- c(0,max(beta1_t1_t19_t37_raw[beta1_t1_t19_t37_raw > 0]))

lower_lim <- min(beta1_t1_t19_t37_raw)
upper_lim <- max(beta1_t1_t19_t37_raw)

png("figures/philly_legend.png", width = 4, height = 0.33, units = "in", res = 300)
par(mar = c(0,0,0,0))
plot(1, type = "n", xlim = c(-1,1), ylim = c(0, 1), xaxt = "n", yaxt = "n", xlab = "", ylab = "")

legend_seq <- seq(par("usr")[1], par("usr")[2], length = 500)
for(leg_ix in 1:499){
  rect(legend_seq[leg_ix], par("usr")[3], legend_seq[leg_ix+1], par("usr")[4],
       border = NA, col = rgb(colorRamp(col_list, bias = 1)((leg_ix-1)/500)/255))
}

dev.off()