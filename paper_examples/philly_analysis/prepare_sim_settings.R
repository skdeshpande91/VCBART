load("philly_crime_data.RData")
n_folds <- 50
p_max <- 4
# randomly hold out 10% of all data
# and try to predict that
set.seed(220)
test_indices <- list()
train_indices <- list()

for(fold_id in 1:n_folds){
  test_indices[[fold_id]] <- sort(sample(1:N, size = floor(N/10), replace = FALSE))
  train_indices[[fold_id]] <- sort( (1:N)[-test_indices[[fold_id]]])
}

method <- c("ind", "cs", "bart")
fold <- 1:n_folds
p <- 1:p_max
sim_settings <- expand.grid(method = method, p = p, fold = fold)
drop_index <- which(sim_settings[,"method"] == "bart" & sim_settings[,"p"] > 1)
sim_settings <- sim_settings[-drop_index,]

save(sim_settings, test_indices, train_indices, file = "sim_settings.RData")
