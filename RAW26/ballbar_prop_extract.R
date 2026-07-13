## ballbar prop data extraction

library(dplyr)
library(tidyr)
setwd("C:/Users/01459189/OneDrive/phd/hpc/RAW26/ballbar_prop")

files <- list.files(
  "C:/Users/01459189/OneDrive/phd/hpc/RAW26/ballbar_prop",
  pattern = "ballbar_.*\\.Rdata$",
  full.names = TRUE
)

parse_settings <- function(name) {
  list(
    dataset = sub("^ballbar_(.*?)_results_prop.*", "\\1", name),
    strata = ifelse(grepl("4strata", name), 4,
                    ifelse(grepl("6strata", name), 6, 5)),
    N = ifelse(grepl("_750", name), 750,
               ifelse(grepl("_250", name), 250, 500)),
    ref_method = sub(".*_(kmedoids|mclust|kmeans)$", "\\1", name)
  )
}

res_list <- list()

base_out <- "C:/Users/01459189/OneDrive/phd/hpc/WIP26/ballbar_prop/"

for (f in files) {
  
  load(f)  # loads object called `store`
  
  fname <- tools::file_path_sans_ext(basename(f))
  settings <- parse_settings(fname)
  
  # ---- extract cleanly ----
  n      <- store[[3]]
  strata <- store[[4]]
  cv     <- store[[5]]
  trace  <- unlist(store[[6]])
  alloc  <- "beth-chrom"
  time   <- unlist(store[[7]])
  
  cv_tot <- rowSums(as.data.frame(cv))
  n_tot <- rowSums(as.data.frame(n))
  
  cv <- data.frame(cv)          # convert matrix to data.frame
  cv$alloc <- alloc             # add character column
  
  # ---- save side objects ----
  save(n, file = paste0(base_out, fname, "_n.RData"))
  save(cv, file = paste0(base_out, fname, "_cv.RData"))
  save(strata, file = paste0(base_out, fname, "_strata.RData"))
  
  # ---- clean tibble row ----
  res_tmp <- tibble(
    method  = "ballbar",
    dataset = settings$dataset,
    strata  = settings$strata,
    nstart  = NA,
    N       = settings$N,
    alloc   = unname(as.vector(alloc)),
    trace   = unname(as.vector(trace)),
    cv_tot  = unname(as.vector(cv_tot)),
    n_tot = unname(as.vector(n_tot)),
    time    = unname(as.vector(time)),
    ref_method = settings$ref_method, 
    run_id  = fname   # VERY useful
  )
  
  res_list[[fname]] <- res_tmp
  
  rm(store)  # avoid accidental reuse
}

# combine everything
res <- bind_rows(res_list)

save(res, file = paste0(base_out, "ballbar_prop_all_res.RData"))