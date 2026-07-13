library(dplyr)

setwd("C:/Users/01459189/OneDrive/phd/hpc/RAW26/app")

datasets  <- c("nids", "ghs", "tb5", "ecd")
sample    <- c(1000, 2000, 500, 3000)
num_strat <- c(4, 5, 3, 6)

methods   <- c("med_ney_", "kmeans_", "mclust_")
method_nm <- c("med_ney", "kmeans", "mclust")

res_list <- list()
counter  <- 1

base_in  <- "C:/Users/01459189/OneDrive/phd/hpc/RAW26/app/"
base_out <- "C:/Users/01459189/OneDrive/phd/hpc/WIP26/app/"

for (d in seq_along(datasets)) {
  for (m in seq_along(methods)) {
    
    filename <- paste0(base_in, methods[m], datasets[d], "_results.RData")
    
    # skip if file doesn't exist (important for HPC runs)
    if (!file.exists(filename)) {
      message("Skipping missing file: ", filename)
      next
    }
    
    load(filename)
    
    # ---- extract ----
    n      <- store[[3]]
    strata <- store[[4]]
    cv     <- store[[5]]
    trace  <- store[[6]]
    alloc  <- store[[10]]
    time   <- store[[11]]
    
    # ---- compute ----
    cv_tot <- rowSums(cv)
    n_tot  <- rowSums(n)
    
    # ---- tidy ----
    res_tmp <- tibble(
      method  = method_nm[m],
      dataset = datasets[d],
      strata  = num_strat[d],
      N       = sample[d],
      alloc   = unname(as.vector(alloc)),
      trace   = unname(as.vector(trace)),
      cv_tot  = unname(as.vector(cv_tot)),
      n_tot   = unname(as.vector(n_tot)),
      time    = unname(as.vector(time))
    )
    
    res_list[[counter]] <- res_tmp
    counter <- counter + 1
    
    rm(store)
  }
}

# combine
res <- bind_rows(res_list)

res$alloc[res$method=="med_ney"] <- "neyman"
res <- res[res$alloc!="prop",]

save(res, file = paste0(base_out, "app_all_res.RData"))