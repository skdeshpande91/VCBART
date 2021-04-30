source("scripts/assess_support_recovery.R")
load("data/p5R20_data.RData")

n_sim <- 25
method_list <- c("vcbart_fixed", "vcbart_adapt", "adapt_dir1", "adapt_dirR")

metrics <- c("tp", "fp", "fn", "tn", "sen", "spec", "prec", "acc", "mcc", "f1")
tmp_results <- matrix(0, nrow = 101, ncol = length(metrics), dimnames = list(c(), metrics))


for(m in method_list){
  
  recovery_theta_mean <- tmp_results
  recovery_varcount_mean <- tmp_results
  recovery_varcount_prob <- tmp_results

  counter <- 0
  
  for(sim_number in 1:n_sim){
    file_name <- paste0("results/variable_selection/", m,"/", m, "_", sim_number, ".RData")
    if(file.exists(file_name)){
      counter <- counter + 1
      load(paste0("results/variable_selection/", m,"/", m, "_", sim_number, ".RData"))
      fit <- get(paste0(m, "_", sim_number))
      tmp_recovery <- assess_support_recovery(fit$varcount_mean, fit$varcount_prob, fit$theta_mean, true_support, R)
      recovery_theta_mean <- recovery_theta_mean + tmp_recovery[["theta_mean"]]
      recovery_varcount_mean <- recovery_varcount_mean + tmp_recovery[["varcount_mean"]]
      recovery_varcount_prob <- recovery_varcount_prob + tmp_recovery[["varcount_prob"]]
      
      rm(list = paste0(m, "_", sim_number))
    }
    

  }
  assign(paste0(m, "_support_recovery"),
         list("theta_mean" = recovery_theta_mean/counter,
              "varcount_mean" = recovery_varcount_mean/counter,
              "varcount_prob" = recovery_varcount_prob/counter))
  rm(recovery_theta_mean, recovery_varcount_mean, recovery_varcount_prob)
}
save(list = paste0(method_list, "_support_recovery"), file = "results/variable_selection/support_recovery_performance.RData")


