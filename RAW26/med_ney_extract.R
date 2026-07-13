## med ney extract

library(dplyr)
setwd("C:/Users/01459189/OneDrive/phd/hpc/RAW26/med_ney")

files <- list.files(
  "C:/Users/01459189/OneDrive/phd/hpc/RAW26/med_ney",
  pattern = "med_ney.*\\.Rdata$",
  full.names = TRUE
)

parse_settings <- function(name) {
  list(
    dataset = sub("^med_ney(.*?)_results.*", "\\1", name),
    strata = ifelse(grepl("4strata", name), 4,
                    ifelse(grepl("6strata", name), 6, 5)),
    nstart = if (grepl("10start", name)) {
      10
    } else if (grepl("25start", name)) {
      25
    } else if (grepl("35start", name)) {
      35
    } else {
      25  # default if nothing found
    },
    N = ifelse(grepl("_750", name), 750,
               ifelse(grepl("_250", name), 250, 500))
  )
}

res_list <- list()

base_out <- "C:/Users/01459189/OneDrive/phd/hpc/WIP26/med_ney/"

for (f in files) {
  
  load(f)  # loads object called `store`
  
  fname <- tools::file_path_sans_ext(basename(f))
  settings <- parse_settings(fname)
  
  # ---- extract cleanly ----
  n      <- store[[3]]
  strata <- store[[4]]
  cv     <- store[[5]]
  trace  <- store[[6]]
  alloc  <- "neyman"
  time   <- store[[11]]
  distance <- store[[10]]
  
  cv_tot <- rowSums(cv)
  
  cv <- data.frame(cv)          # convert matrix to data.frame
  cv$alloc <- alloc             # add character column
  cv$distance <- distance
  # ---- save side objects ----
  save(n, file = paste0(base_out, fname, "_n.RData"))
  save(cv, file = paste0(base_out, fname, "_cv.RData"))
  save(strata, file = paste0(base_out, fname, "_strata.RData"))
  
  # ---- clean tibble row ----
  res_tmp <- tibble(
    method  = "med_ney",
    dataset = settings$dataset,
    strata  = settings$strata,
    nstart  = settings$nstart,
    N       = settings$N,
    alloc   = unname(as.vector(alloc)),
    trace   = unname(as.vector(trace)),
    cv_tot  = unname(as.vector(cv_tot)),
    time    = unname(as.vector(time)),
    distance    = unname(as.vector(distance)),
    run_id  = fname   # VERY useful
  )
  
  res_list[[fname]] <- res_tmp
  
  rm(store)  # avoid accidental reuse
}

# combine everything
res <- bind_rows(res_list)

save(res, file = paste0(base_out, "med_ney_all_res.RData"))