library(igraph)
library(dplyr)
library(ggplot2)

###################
# Map of census tracts colored by average density
##################

load("tracts.RData") 
load("yearly_crimes.RData")
load("philly_crime_data.RData")
tract_data <- 
  as_tibble(tracts@data) %>%
  mutate("TRACT" = paste0("tr.", GEOID10))

tract_fort <- ggplot2::fortify(tracts)
id <- as.character(unique(tract_fort$id))
tract_names <- paste0("tr.", tracts@data[,"GEOID10"])

#######
# Get yearly crime counts & density

ihst <- function(x){
  return(log(x + (x^2 +1)^2) - log(2))
}

n_tracts <- nrow(tracts@data)
n_times <- 16
Y_mat <- matrix(nrow = n_tracts, ncol = n_times, dimnames = list(tract_names, c()))


for(tr in rownames(Y_mat)){
  tmp_y <- as.vector(unlist(crimes[which(crimes$TRACT == tr),"transformed_density"]))
  Y_mat[tr,] <- tmp_y
}



mean_levels <- rowMeans(Y_mat)
tmp_data <- data.frame(id = id, 
                       tract = rownames(Y_mat),
                       mean_levels = mean_levels)
plot_data <- dplyr::inner_join(tract_fort, tmp_data, by = "id")


mean_levels_plot <- 
  ggplot() + 
  geom_polygon(data = plot_data, aes(x = long, y = lat, group = group, fill = mean_levels),
               color = alpha("black", 0.15), alpha = 1, show.legend = TRUE) + 
  scale_fill_distiller(type = "div", palette = "PRGn") + 
  theme(panel.background = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(),
        axis.title = element_blank(), 
        aspect.ratio = 0.8,
        legend.title = element_blank())

ggsave(filename = "philly_map.pdf", width = 6, height = 6, units = "in",
       plot = mean_levels_plot)
##########################
# Selected time series

index <- which(tracts@data[,"NAME10"] == "9800")
my_tract <- paste0("tr.", tracts@data[index, "GEOID10"])
adj_tracts <- rownames(A_tract)[A_tract[,my_tract] == 1]
my_tracts <- c(my_tract, adj_tracts)


col_list <- rev(RColorBrewer::brewer.pal(n = 11, name = "PRGn")) # define a color palette

rescaled_means <- scales::rescale(mean_levels, to = c(0,1))

pdf("philly_time_series.pdf",
    width = 4.5, height = 3)
par(mar = c(3,3,1,1), mgp = c(1.8, 0.5, 0))
plot(1, type = "n",xlim = c(1, n_times), xlab = "Time", xaxt = "n", 
     ylab = "Density",ylim = range(Y_mat), main = "")
for(tr in my_tracts){
  if(mean_levels[tr] >= 25){
    lines(1:n_times, Y_mat[tr,], lwd = 1,
          col = rgb(colorRamp(col_list,bias=1)(rescaled_means[tr])/255))
    points(1:n_times, Y_mat[tr,], 
           col = rgb(colorRamp(col_list,bias=1)(rescaled_means[tr])/255),
           pch = 16, cex = 0.75)
  } else if(tr == my_tract){
    lines(1:n_times, Y_mat[tr,], lwd = 1,
          col = rgb(colorRamp(col_list,bias=1)(rescaled_means[tr])/255))
    points(1:n_times, Y_mat[tr,],
           col = rgb(colorRamp(col_list,bias=1)(rescaled_means[tr])/255),
           pch = 16)
  }
}
dev.off()


# Network representation
pdf("philly_network.pdf", width = 6, height = 6)
par(mar = c(1,1,1,1), mgp = c(1.8, 0.5, 0))
plot(g, layout = tract_layout,
     vertex.size = 4, vertex.label = NA, vertex.color = "gray")
dev.off()
