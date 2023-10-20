load("philly_crime_data.RData")
method <- c("ind", "cs")
tract <- 1:n
p <- 1:p_max

sim_settings <- expand.grid(tract = tract, p = p, method = method)
save(sim_settings, file = "sim_settings.RData")
