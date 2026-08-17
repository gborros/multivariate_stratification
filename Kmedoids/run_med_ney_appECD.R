library(dplyr)
library(cluster)
library(haven)
library(doParallel)
library(foreach)

############# source functions #############

source("function_calc_variance.R")
source("calculate_cv_nofpc.R")
source("function_sample_size.R")


#------------------------ Load Data ---------------------------------#
load("multivar_datasets/ecd_census_multivariate.RData")
ecd <- df

#------------------------ Variable combinations ---------------------#
var_combos <- list(
  c("fees_amount", "count_register_all", "age"),
  c("fees_amount", "count_register_all", "age", "count_staff_all"),
  c("fees_amount", "count_register_all", "age", "count_staff_all", "hours")
)

datasets <- c("ecd")
dfs <- list(ecd)

sample <- 3000
num_strata <- 6
#------------------------ Setup functions --------------------------#

convert_labelled <- function(df){
  
  df <- haven::zap_labels(df)
  
  df <- as.data.frame(lapply(df, function(x){
    
    if(is.factor(x)){
      as.numeric(as.character(x))
    } else {
      as.numeric(x)
    }
    
  }))
  
  return(df)
}



set_up <- function(df, sort_var){
  
  df <- convert_labelled(df)
  
  df <- as.data.frame(scale(na.omit(df)))
  
  df <- df[order(df[[sort_var]]),]
  
  rownames(df) <- NULL
  
  df
  
}



set_up_ns <- function(df, sort_var){
  
  df <- convert_labelled(df)
  
  df <- as.data.frame(na.omit(df))
  
  df <- df[order(df[[sort_var]]),]
  
  rownames(df) <- NULL
  
  df
  
}

#------------------------ Parallel setup ---------------------------#

cl <- makeCluster(12)

registerDoParallel(cl)



#------------------------ Run -------------------------------------#

for(dta in seq_along(datasets)){
  
  
  sf_full <- dfs[[dta]]
  
  
  task_grid <- expand.grid(
    combo_id = seq_along(var_combos),
    iter = 1:40
  )
  
  
  results <- foreach(
    t = 1:nrow(task_grid),
    .packages = c("cluster","dplyr"),
    .combine = "list",
    .multicombine = TRUE,
    .export = c(
      "convert_labelled",
      "set_up",
      "set_up_ns",
      "calculate_sample_size",
      "calculate_cv",
      "calculate_variance"
    )
  ) %dopar% {
    
    
    combo_id <- task_grid$combo_id[t]
    
    iter <- task_grid$iter[t]
    
    
    vars_orig <- var_combos[[combo_id]]
    
    num_vars <- length(vars_orig)
    
    
    sf <- sf_full[,vars_orig]
    
    
    vars <- paste0("X",1:ncol(sf))
    
    colnames(sf) <- vars
    
    
    df <- set_up(sf,"X1")
    
    df_ns <- set_up_ns(sf,"X1")
    
    
    set.seed(iter)
    
    
    start <- Sys.time()
    
    
    med <- tryCatch(
      pam(
        x=df,
        k=num_strata,
        metric="euclidean",
        nstart=25
      ),
      error=function(e) NULL
    )
    
    
    time <- as.numeric(
      difftime(Sys.time(),start,units="mins")
    )
    
    
    # Build per-(combo_id, iter) output filename - unique per task,

    iter_file <- paste0(
      "OUTPUT/med_ney_",
      datasets[dta],
      "_vars", num_vars,
      "_combo", combo_id,
      "_iter", iter,
      "_results.RData"
    )
    
    
    if(is.null(med)){
      
      iter_result <- list(
        dataset=datasets[dta],
        method="med_ney",
        combo_id=combo_id,
        num_vars=num_vars,
        iter=iter,
        error="med_ney_failed"
      )
      
      # Save immediately, even on failure - so a missing file on disk
      # unambiguously means the worker died before getting here,
      # not that it failed cleanly.
      save(iter_result, file=iter_file)
      
      return(iter_result)
      
    }
    
    
    strata <- med$clustering
    
    
    df_clust <- data.frame(
      df,
      strata=strata
    )
    
    df_ns$strata <- strata
    
    
    n <- calculate_sample_size(
      df=df_clust,
      strata=strata,
      vars=vars,
      method="neyman",
      ssize=sample
    )
    
    
    cv <- calculate_cv(
      df=df_ns,
      strata=strata,
      vars=vars,
      n=n
    )
    
    
    trace <- calculate_variance(df_clust,strata,vars,n,"trace")
    
    det <- calculate_variance(df_clust,strata,vars,n,"determinant")
    
    varcov <- calculate_variance(df_clust,strata,vars,n,"varcov")
    
    pca <- calculate_variance(df_clust,strata,vars,n,"pca")
    
    
    iter_result <- list(
      dataset=datasets[dta],
      method="med_ney",
      combo_id=combo_id,
      num_vars=num_vars,
      iter=iter,
      n=list(n),
      cv=list(cv),
      trace=trace,
      det=det,
      varcov=varcov,
      pca=pca,
      time=time,
      strata=list(strata)
    )
    
    # Save this iteration's result immediately, from inside the worker.
    # This is the key change: each (combo_id, iter) pair gets its own
    # file, written as soon as it's computed - so if this worker is
    # killed (e.g. SLURM OOM) before the foreach() call as a whole
    # finishes, everything it already completed is still on disk.
    save(iter_result, file=iter_file)
    
    iter_result
    
  }
  
  

  }


stopCluster(cl)