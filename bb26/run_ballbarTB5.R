## ============================================================
## Ballin-Barcaroli GGA benchmark (optimStrata "continuous")
## ============================================================

library(dplyr)
library(haven)
library(doParallel)
library(foreach)
library(SamplingStrata)

source("calculate_cv_nofpc.R")
source("function_calc_variance.R")
source("function_sample_size.R")


ncores <- 40

methods <- c("kmeans", "med_ney", "mclust")

n_iter <- 40

dataset_specs <- list(
  tb5 = list(
    path        = "multivar_datasets/tb5.RData",
    obj         = "df",
    ssize       = 500,
    num_strata  = 3,
    subset_cols = 2:7
  )
  
)

# Variable combinations per dataset. 
var_combos <- list(
  tb5 = list(
    c("age", "domain_1", "domain_2", "domain_3"),
    c("age", "domain_1", "domain_2", "domain_3", "domain_4"),
    c("age", "domain_1", "domain_2", "domain_3", "domain_4", "domain_5")
  )
  
)

## ------------------------------------------------------------------
## Load consolidated CV targets
## ------------------------------------------------------------------

load("cv_all_res.RData")   # provides `best_by_group`

## ------------------------------------------------------------------
## Data setup helpers (generalised from your two scripts)
## ------------------------------------------------------------------

convert_labelled <- function(df) {
  df <- haven::zap_labels(df)
  df <- as.data.frame(lapply(df, function(x) {
    if (is.factor(x)) as.numeric(as.character(x)) else as.numeric(x)
  }))
  df
}

set_up <- function(df, sort_var) {
  df <- as.data.frame(na.omit(df))
  df <- df[order(df[[sort_var]]), ]
  rownames(df) <- 1:nrow(df)
  df
}

set_up_ws <- function(df, sort_var) {
  df <- as.data.frame(scale(na.omit(df)))
  df <- df[order(df[[sort_var]]), ]
  rownames(df) <- 1:nrow(df)
  df
}

## ------------------------------------------------------------------
## Load all datasets up front into a single list, keyed by name.
## ------------------------------------------------------------------

dfs_all <- list()
for (dsname in names(dataset_specs)) {
  spec <- dataset_specs[[dsname]]
  local_env <- new.env()
  load(spec$path, envir = local_env)
  raw <- local_env[[spec$obj]]
  if (!is.null(spec$subset_cols)) {
    raw <- raw[, spec$subset_cols]
  }
  dfs_all[[dsname]] <- convert_labelled(raw)
}

## ------------------------------------------------------------------
## Error-vector lookup from best_by_group.
## ------------------------------------------------------------------

get_error_vector <- function(best_by_group, this_method, this_dataset,
                             num_vars, this_num_strata, this_ssize) {
  
  row <- best_by_group %>%
    filter(
      method  == this_method,
      dataset == this_dataset,
      strata  == this_num_strata,
      N       == this_ssize,
      vars    == num_vars
    )
  
  if (nrow(row) != 1) {
    stop(sprintf(
      "get_error_vector: expected exactly 1 matching row for method=%s dataset=%s strata=%s N=%s num_vars=%s, got %d",
      this_method, this_dataset, this_num_strata, this_ssize, num_vars, nrow(row)
    ))
  }
  
  # Take exactly one cv column per variable, in order (cv_vec_1 ... cv_vec_k)
  cv_cols <- paste0("cv_vec_", 1:num_vars)
  as.numeric(row[1, cv_cols])
}

## ------------------------------------------------------------------
## Build the task grid: method x dataset x combo_id x iter
## ------------------------------------------------------------------

task_list <- list()
for (m in methods) {
  for (dsname in names(dataset_specs)) {
    combos <- var_combos[[dsname]]
    if (is.null(combos)) next
    for (combo_id in seq_along(combos)) {
      for (iter in 1:n_iter) {
        task_list[[length(task_list) + 1]] <- data.frame(
          method   = m,
          dataset  = dsname,
          combo_id = combo_id,
          iter     = iter,
          stringsAsFactors = FALSE
        )
      }
    }
  }
}
task_grid <- do.call(rbind, task_list)

cat("Total tasks:", nrow(task_grid), "\n")

## ------------------------------------------------------------------
## Parallel cluster - single cluster
## ------------------------------------------------------------------

cl <- makeCluster(ncores)
clusterEvalQ(cl, {
  library(dplyr)
  library(haven)
  library(SamplingStrata)
  source("calculate_cv_nofpc.R")
  source("function_calc_variance.R")
})
registerDoParallel(cl)

results <- foreach(
  t = 1:nrow(task_grid),
  .packages = c("dplyr", "haven", "SamplingStrata"),
  .export = c(
    "dataset_specs", "var_combos", "dfs_all", "best_by_group",
    "set_up", "set_up_ws", "get_error_vector",
    "calculate_cv", "calculate_variance"
  )
) %dopar% {
  
  method   <- task_grid$method[t]
  dsname   <- task_grid$dataset[t]
  combo_id <- task_grid$combo_id[t]
  iter     <- task_grid$iter[t]
  
  spec       <- dataset_specs[[dsname]]
  ssize      <- spec$ssize
  num_strata <- spec$num_strata
  vars_orig  <- var_combos[[dsname]][[combo_id]]
  num_vars   <- length(vars_orig)
  
  iter_file <- paste0(
    "OUTPUT/ballbar_", method, dsname,
    "_vars", num_vars, "_combo", combo_id,
    "_iter", iter, "_results.RData"
  )
  
  sf <- dfs_all[[dsname]][, vars_orig, drop = FALSE]
  vars <- paste0("X", seq_along(vars_orig))
  colnames(sf) <- vars
  
  df    <- set_up(sf, "X1")
  df_ws <- set_up_ws(sf, "X1")
  
  err_lookup <- tryCatch(
    get_error_vector(best_by_group, method, dsname, num_vars, num_strata, ssize),
    error = function(e) e
  )
  
  if (inherits(err_lookup, "error")) {
    iter_result <- list(
      method = method, dataset = dsname, combo_id = combo_id,
      num_vars = num_vars, iter = iter,
      error = paste("error_vector_lookup_failed:", conditionMessage(err_lookup))
    )
    save(iter_result, file = iter_file)
    return(iter_result)
  }
  
  error_row <- as.data.frame(t(err_lookup))
  cfs <- paste0("CV", seq_along(vars))
  colnames(error_row) <- cfs
  
  df$id  <- 1:nrow(df)
  df$dom <- 1
  
  frame <- buildFrameDF(df = df, id = "id", X = vars, Y = vars, domainvalue = "dom")
  
  error <- as.data.frame(list(
    DOM = "DOM1",
    error_row,
    domainvalue = 1
  ))
  
  set.seed(iter)
  init_attempt <- tryCatch({
    init_sol3 <- KmeansSolution2(
      frame = frame, errors = error,
      maxclusters = num_strata, nstrata = num_strata
    )
    nstrata3 <- tapply(
      init_sol3$suggestions, init_sol3$domainvalue,
      FUN = function(x) length(unique(x))
    )
    initial_solution3 <- prepareSuggestion(init_sol3, frame, nstrata3)
    list(ok = TRUE, nstrata3 = nstrata3, initial_solution3 = initial_solution3)
  }, error = function(e) {
    list(ok = FALSE, error_msg = conditionMessage(e))
  })
  
  set.seed(iter)
  if (init_attempt$ok && init_attempt$nstrata3 == num_strata) {
    solution <- optimStrata(
      method = "continuous", errors = error, framesamp = frame,
      iter = 2000, pops = 50, nStrata = num_strata,
      mut_chance = 0.3, suggestions = init_attempt$initial_solution3
    )
  } else {
    solution <- optimStrata(
      method = "continuous", errors = error, framesamp = frame,
      iter = 2000, pops = 50, nStrata = num_strata, mut_chance = 0.3
    )
  }
  
  strataStructure <- summaryStrata(solution$framenew, solution$aggr_strata, progress = FALSE)
  n <- as.vector(strataStructure$Allocation)
  n_strata_realized <- length(n)
  
  df_sol <- solution$framenew
  df_sol$strata <- df_sol$STRATO
  
  cv_val <- calculate_cv(df = df_sol, strata = df_sol$strata, vars = vars, n = n)
  trace  <- calculate_variance(df = df_ws, strata = df_sol$strata, vars = vars, n = n, objective = "trace")
  det    <- calculate_variance(df = df_ws, strata = df_sol$strata, vars = vars, n = n, objective = "determinant")
  varcov <- calculate_variance(df = df_ws, strata = df_sol$strata, vars = vars, n = n, objective = "varcov")
  pca    <- calculate_variance(df = df_ws, strata = df_sol$strata, vars = vars, n = n, objective = "pca")
  
  iter_result <- list(
    method = method, dataset = dsname, combo_id = combo_id, num_vars = num_vars,
    iter = iter, n = n, cv = cv_val, trace = trace, det = det, varcov = varcov,
    pca = pca, objective = solution$objective, strata = df_sol$strata,
    n_strata_realized = n_strata_realized,
    warmstart_ok = init_attempt$ok,
    warmstart_error = if (!init_attempt$ok) init_attempt$error_msg else NA
  )
  
  # Save immediately from inside the worker
  save(iter_result, file = iter_file)
  
  iter_result
}

stopCluster(cl)