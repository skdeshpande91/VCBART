# Need to fix some filenames


for(sim_number in 1:25){
  fixed_cs75_filename <- paste0("results/HRS/vcbart_fixed/vcbart_fixed_cs75_", sim_number, ".RData")
  if(file.exists(fixed_cs75_filename)){
    load(fixed_cs75_filename)
    assign(paste0("vcbart_fixed_cs75_", sim_number), get(paste0("vcbart_adapt_", sim_number)))
    rm(list = paste0("vcbart_adapt_", sim_number))
    save(list = paste0("vcbart_fixed_cs75_", sim_number), file = fixed_cs75_filename)
  }
}