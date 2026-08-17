## Run File - loops over methods automatically

library(dplyr)
library(haven)
library(doParallel)
library(SamplingStrata)
source("calculate_cv_nofpc.R")
source("function_calc_variance.R")
source("function_sample_size.R")

ncores <- 40

## ---- methods to run, one after another ----
methods <- c("kmeans_", "med_ney", "mclust_")
#methods <- c("kmeans_", "mclust_")


for (method in methods) {

  cl <- makeCluster(ncores)

  clusterEvalQ(cl, {
    library(dplyr)
    library(haven)
    library(SamplingStrata)
    source("calculate_cv_nofpc.R")
    source("function_calc_variance.R")
  })
  registerDoParallel(cl)

  ## comb5
  load("multivar_datasets/d_comb5.RData")
  d_comb5 <- as.data.frame(d_comb5)
  load(paste0("CV/", method, "d_comb5_results_cv.RData"))

  cv <- cv[cv$alloc == "neyman", ]
  df_cv <- cv[, 1:4]
  row_sums <- rowSums(df_cv)
  min_row_index <- which.min(row_sums)
  df_cv <- as.matrix(df_cv[min_row_index, ])

  datasets <- c("d_comb5")
  dfs <- list(d_comb5)
  sample <- c(500)
  cv_list <- list(df_cv)

  for (dta in 1:length(datasets)) {
    ssize <- sample[1]
    vars <- character()
    cfs <- character()
    sf <- dfs[[dta]]
    error <- as.matrix(cv_list[[dta]], nrow = 1)

    for (z in 1:ncol(sf)) {
      var <- paste0("X", z)
      vars <- c(vars, var)
      cf <- paste0("CV", z)
      cfs <- c(cfs, cf)
    }

    colnames(sf) <- vars
    colnames(error) <- cfs
    #------------------------ Section 2: Dataset set up --------------------------

    set_up <- function(df, sort_var) {
      df <- as.data.frame(na.omit(df))
      df <- df[order(df[[sort_var]]), ]
      rownames(df) <- 1:nrow(df)
      return(df)
    }

    df <- set_up(df = sf, sort_var = "X1")

    set_up_ws <- function(df, sort_var) {
      df <- as.data.frame(scale(na.omit(df)))
      df <- df[order(df[[sort_var]]), ]
      rownames(df) <- 1:nrow(df)
      return(df)
    }

    df_ws <- set_up_ws(df = sf, sort_var = "X1")

    num_strata <- 5

    #------------------------ Section 3: Run GGA ----------------------------------
    store_n_all      <- list()
    store_cv_all     <- list()
    store_var_all    <- list()
    store_trace_all  <- list()
    store_strata_all <- list()
    store_time_all   <- list()

    df <- as.data.frame(df)
    error <- as.data.frame(error)

    df$id <- 1:nrow(df)
    df$dom <- 1

    frame <- buildFrameDF(df = df,
                           id = "id",
                           X = vars,
                           Y = vars,
                           domainvalue = "dom")

    ndom <- length(unique(df$dom))
    error <- as.data.frame(list(DOM = rep("DOM1", ndom),
                                 error,
                                 domainvalue = c(1)))

    # ------------------------ Parallelised runs ------------------------
    results <- foreach(iter = 1:40, .packages = c("dplyr", "haven", "SamplingStrata")) %dopar% {

      set.seed(iter)
      start.time <- Sys.time()

      ## Try the KmeansSolution2 warm-start first. If it errors, skip the
      ## warm-start and let optimStrata run without prescribed suggestions,
      ## instead of crashing the whole iteration.
      init_attempt <- tryCatch({
        init_sol3 <- KmeansSolution2(frame = frame,
                                      errors = error,
                                      maxclusters = num_strata,
                                      minnumstrat = num_strata)

        nstrata3 <- tapply(init_sol3$suggestions,
                            init_sol3$domainvalue,
                            FUN = function(x) length(unique(x)))

        initial_solution3 <- prepareSuggestion(init_sol3, frame, nstrata3)

        list(ok = TRUE, nstrata3 = nstrata3, initial_solution3 = initial_solution3)
      }, error = function(e) {
        list(ok = FALSE, error_msg = conditionMessage(e))
      })

      set.seed(iter)
      if (init_attempt$ok && init_attempt$nstrata3 == num_strata) {
        solution <- optimStrata(method = "continuous",
                                 errors = error,
                                 framesamp = frame,
                                 iter = 2000,
                                 pops = 50,
                                 nStrata = num_strata,
                                 mut_chance = 0.3,
                                 suggestions = init_attempt$initial_solution3)
      } else {
        solution <- optimStrata(method = "continuous",
                                 errors = error,
                                 framesamp = frame,
                                 iter = 2000,
                                 pops = 50,
                                 nStrata = num_strata,
                                 mut_chance = 0.3)
      }
      end.time <- Sys.time()
      time <- difftime(end.time, start.time, units = "mins")

      strataStructure <- summaryStrata(solution$framenew,
                                        solution$aggr_strata,
                                        progress = FALSE)

      n <- as.vector(strataStructure$Allocation)
      n_strata_realized <- length(n)

      df_sol <- solution$framenew
      df_sol$strata <- df_sol$STRATO

      cv_val <- calculate_cv(df = df_sol, strata = df_sol$strata, vars = vars, n = n)

      trace <- calculate_variance(df = df_ws,
                                   strata = df_sol$strata,
                                   vars = vars,
                                   n = n,
                                   objective = 'trace')

      list(
        n = n,
        cv = cv_val,
        trace = trace,
        objective = solution$objective,
        strata = df_sol$strata,
        time = time,
        n_strata_realized = n_strata_realized,
        warmstart_ok = init_attempt$ok,          # TRUE/FALSE flag per iteration
        warmstart_error = if (!init_attempt$ok) init_attempt$error_msg else NA
      )
    } # end foreach

    # ------------------------ Store results ------------------------
    store_n_all                       <- lapply(results, `[[`, "n")
    store_cv_all[[datasets[dta]]]     <- do.call(rbind, lapply(results, `[[`, "cv"))
    store_trace_all[[datasets[dta]]]  <- do.call(rbind, lapply(results, `[[`, "trace"))
    store_strata_all[[datasets[dta]]] <- lapply(results, `[[`, "strata")
    store_time_all[[datasets[dta]]]   <- do.call(rbind, lapply(results, `[[`, "time"))
    store_n_strata_all                <- sapply(results, `[[`, "n_strata_realized")
    store_warmstart_ok                <- sapply(results, `[[`, "warmstart_ok")
  }

  stopCluster(cl)

  # ------------------------ Save results ------------------------
  store <- list(
    df = df,
    df_ws = df_ws,
    store_n_all = store_n_all,
    store_strata_all = store_strata_all,
    store_cv_all = store_cv_all,
    store_trace_all = store_trace_all,
    store_time_all = store_time_all,
    store_n_strata_all = store_n_strata_all,
    store_warmstart_ok = store_warmstart_ok
  )

  # dynamic output filename per method
  filename <- paste0("OUTPUT/ballbar_", datasets[1], "_results_", method, ".Rdata")
  save(store, file = filename)

  cat("Finished method:", method, "-> saved to", filename, "\n")

} # end methods loop
