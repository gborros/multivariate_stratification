## Run File

library(dplyr)
library(haven)
library(doParallel)
library(SamplingStrata)
source("calculate_cv.R")
source("function_calc_variance.R")
source("function_sample_size.R")

ncores <- 40
cl <- makeCluster(ncores)

clusterEvalQ(cl, {
  library(dplyr)
  library(haven)
  library(SamplingStrata)
  source("calculate_cv.R")
  source("function_calc_variance.R")
})
registerDoParallel(cl)

## COMB1
load("multivar_datasets/d_comb1.RData")
d_comb1 <- as.data.frame(d_comb1)
load("CV/kmeans_d_comb1_results_6strata_cv.RData")
cv <- cv[cv$alloc=="prop",]
df_cv <- cv[,1:4]
row_sums <- rowSums(df_cv)
min_row_index <- which.min(row_sums)
df_cv <- as.matrix(df_cv[min_row_index, ])

datasets <- c("d_comb1")
dfs <- list(d_comb1)
sample <- c(500)
cv <- list(df_cv)

for (dta in 1:length(datasets)) {
  ssize = sample[1]
  vars <- character()
  cfs <- character()
  sf <- dfs[[dta]]
  error <- as.matrix(cv[[dta]], nrow=1)
  
  for (z in 1:ncol(sf)) {
    
    var <- paste0("X", z)
    vars <- c(vars, var)
    cf <- paste0("CV", z)
    cfs <- c(cfs, cf)
  }
  
  colnames(sf) <- vars
  colnames(error) <- cfs
  #------------------------ Section 2: Dataset set up --------------------------
  
  set_up <- function(df, sort_var) { ## set up function to ensure scaled and no na data
    df <- as.data.frame(na.omit(df))
    df <- df[order(df[[sort_var]]), ]
    rownames(df) <- 1:nrow(df)
    return(df)
  }
  
  df <- set_up(df = sf, sort_var = "X1")
  
  set_up_ws <- function(df, sort_var) { ## set up function to ensure scaled and no na data
    df <- as.data.frame(scale(na.omit(df)))
    df <- df[order(df[[sort_var]]), ]
    rownames(df) <- 1:nrow(df)
    return(df)
  }
  
  df_ws <- set_up_ws(df = sf, sort_var = "X1")
  
  num_strata <- 6
  
  #------------------------ Section 3: Run GGA ----------------------------------
  # Initialize storage lists
  store_n_all      <- list()
  store_cv_all     <- list()
  store_var_all    <- list()
  store_trace_all  <- list()
  store_strata_all <- list()
  store_time_all <- list()
  
  df <- as.data.frame(df)
  error <- as.data.frame(error)
  
  df$id <- 1:nrow(df) ## adding identifier
  df$dom <- 1 ## adding domain
  
  frame <- buildFrameDF(df=df, ## now taking full dataset
                        id = "id", ## id variable
                        X = vars,
                        Y = vars, ## 
                        domainvalue = "dom")
  
  ## specify error thresholds:
  ndom <- length(unique(df$dom))
  error <- as.data.frame(list(DOM=rep("DOM1",ndom),
                              error,
                              domainvalue=c(1) ))
  
  
  # ------------------------ Parallelised runs ------------------------
  results <- foreach(iter = 1:40, .packages = c("dplyr", "haven", "SamplingStrata")) %dopar% {
    
    set.seed(iter)
    start.time <- Sys.time()
    ## kmeans
    # doesn't want to work with prescribed strata #
       
    ## GGA
    set.seed(iter)

    solution <- optimStrata(method = "continuous",
                              errors = error, 
                              framesamp = frame,
                              iter = 2000,
                              pops = 50,
                              nStrata = num_strata,
                              mut_chance = 0.3)
    
    end.time <- Sys.time()
    time <- difftime(end.time, start.time, units="mins")
    
    strataStructure <- summaryStrata(solution$framenew,
                                     solution$aggr_strata,
                                     progress = FALSE)
    
    n <- as.vector(strataStructure$Allocation)
    n_strata_realized <- length(n)   # <-- add this: actual number of strata this run produced
    
    df_sol <- solution$framenew
    df_sol$strata <- df_sol$STRATO
    
    ## CV and other objectives
    cv <- calculate_cv(df = df_sol, strata = df_sol$strata, vars = vars, n = n)
    
    trace <- calculate_variance(df = df_ws,
                                strata = df_sol$strata,
                                vars = vars,
                                n = n,
                                objective = 'trace')
    
    
    # Return all results as a list
    list(
      n = n,
      cv = cv,
      trace = trace,
      objective = solution$objective,
      strata = df_sol$strata,
      time=time,
      n_strata_realized = n_strata_realized
    )
  } # end foreach
  
  # ------------------------ Store results ------------------------
  store_n_all           <- lapply(results, `[[`, "n")  
  store_cv_all[[datasets[dta]]]     <- do.call(rbind, lapply(results, `[[`, "cv"))
  store_trace_all[[datasets[dta]]]  <- do.call(rbind, lapply(results, `[[`, "trace"))
  store_strata_all[[datasets[dta]]] <- lapply(results, `[[`, "strata")
  store_time_all[[datasets[dta]]]  <- do.call(rbind, lapply(results, `[[`, "time"))
  store_n_strata_all    <- sapply(results, `[[`, "n_strata_realized") 
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
  store_n_strata_all = store_n_strata_all
)

filename <- paste0("OUTPUT/ballbar_", datasets[dta], "_results_6strata_prop.Rdata")
save(store, file = filename)