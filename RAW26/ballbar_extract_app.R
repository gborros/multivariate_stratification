## ballbar data extraction (application)

library(dplyr)
library(tidyr)

setwd("C:/Users/01459189/OneDrive/phd/hpc/RAW26/ballbar_app")

datasets  <- c("nids", "ghs", "tb5", "ecd")
sample    <- c(1000, 2000, 500, 3000)
num_strat <- c(4, 5, 3, 6)


res_list <- list()
counter  <- 1

base_in  <- "C:/Users/01459189/OneDrive/phd/hpc/RAW26/ballbar_app/"
base_out <- "C:/Users/01459189/OneDrive/phd/hpc/WIP26/ballbar_app/"
methods   <- c("kmeans", "kmedoids", "mclust")

res_list <- list()
counter  <- 1

for (d in seq_along(datasets)) {
  for (m in seq_along(methods)) {
    
    filename <- paste0(base_in,
                       "ballbar_", datasets[d],
                       "_results_", methods[m], ".RData")
    
    if (!file.exists(filename)) {
      message("Skipping missing file: ", filename)
      next
    }
    
    load(filename)   # loads `store`
    
    # ---- extract ----
    n      <- data.frame(store[[3]])
    strata <- store[[4]]
    cv     <- data.frame(store[[5]])
    trace  <- store[[6]]
    time   <- unlist(store[[7]])
    
    # ---- compute ----
    cv_tot <- rowSums(cv)
    n_tot  <- rowSums(n)
    
    # ---- tidy ----
    res_tmp <- tibble(
      method  = methods[m],   # <-- now varies
      dataset = datasets[d],
      strata  = num_strat[d],
      N       = sample[d],
      alloc   = "Bethel-Chromy",
      trace   = unname(as.vector(trace)),
      cv_tot  = unname(as.vector(cv_tot)),
      n_tot   = unname(as.vector(n_tot))
    )
    
    res_list[[counter]] <- res_tmp
    counter <- counter + 1
    
    rm(store)
  }
}

res <- bind_rows(res_list)

save(res, file = paste0(base_out, "ballbar_app_all_res.RData"))

