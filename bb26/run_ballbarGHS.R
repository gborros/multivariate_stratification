## Run File

library(dplyr)
library(haven)
library(doParallel)
library(SamplingStrata)
source("calculate_cv.R")
source("function_calc_variance.R")
source("function_sample_size.R")

#ncores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", 1))
ncores <- 20
cl <- makeCluster(ncores)
#wd <- getwd()   # your current working directory
#clusterExport(cl, "wd")
#clusterEvalQ(cl, setwd(wd))
clusterEvalQ(cl, {
  library(dplyr)
  library(haven)
  library(SamplingStrata)
  source("calculate_cv.R")
  source("function_calc_variance.R")
})
registerDoParallel(cl)

## ECD Census
# load("multivar_datasets/ecd_census_multivariate.RData")
# ecd <- df
# load("OUTPUT/GGA_ney_ecd_results.Rdata")
# ecd_cv<- store[[5]]
# ecd_cv <- ecd_cv[[1]]
# row_sums <- rowSums(ecd_cv)
# min_row_index <- which.min(row_sums)
# ecd_cv <- t(as.matrix(ecd_cv[min_row_index, ]))

# datasets <- c("ecd")
# dfs <- list(ecd)
# sample <- c(3000)
# cv <- list(ecd_cv)


# ## TB5
# load("multivar_datasets/tb5.RData")
# tb5 <- df[2:7]
# load("OUTPUT/GGA_ney_tb5_results.Rdata")
# tb5_cv<- store[[5]]
# tb5_cv <- tb5_cv[[1]]
# row_sums <- rowSums(tb5_cv)       # sum across each row
# min_row_index <- which.min(row_sums)  # index of row with smallest sum
# tb5_cv <- t(as.matrix(tb5_cv[min_row_index, ]))
# 
# 
# datasets <- c("tb5")
# dfs <- list(tb5)
# sample <- c(500)
# cv <- list(tb5_cv)


# ## GHS
load("multivar_datasets/ghs_2024.RData")
df$hhinc_pc <- df$totmhinc/df$hholdsz
df <- df[df$fin_reqinc!=9999999,]
ghs <- df[,-4]
load("CV/kmeans_ghs_results.Rdata")
ghs_cv<- store[[5]]
row_sums <- rowSums(ghs_cv)
min_row_index <- which.min(row_sums)
ghs_cv <- t(as.matrix(ghs_cv[min_row_index, ]))

datasets <- c("ghs")
dfs <- list(ghs)
sample <- c(2000)
cv <- list(ghs_cv)

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
  
  num_strata <- 5
  
  #------------------------ Section 3: Run GGA ----------------------------------
  # Initialize storage lists
  store_n_all      <- list()
  store_cv_all     <- list()
  store_var_all    <- list()
  store_trace_all  <- list()
  store_det_all    <- list()
  store_varcov_all <- list()
  store_pca_all    <- list()
  store_strata_all <- list()
  
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
  results <- foreach(iter = 1:10, .packages = c("dplyr", "haven", "SamplingStrata")) %dopar% {
    
    set.seed(iter)
    
    ## kmeans
    init_sol3 <- KmeansSolution2(frame=frame,
                                 errors=error,
                                 maxclusters = num_strata,
                                 nstrata=num_strata)
    
    nstrata3 <- tapply(init_sol3$suggestions,
                       init_sol3$domainvalue,
                       FUN=function(x) length(unique(x)))
    
    initial_solution3 <- prepareSuggestion(init_sol3, frame, nstrata3)
    
    ## GGA
    set.seed(iter)
    if(nstrata3 != num_strata) {
      solution <- optimStrata(method = "continuous",
                              errors = error, 
                              framesamp = frame,
                              iter = 2000,
                              pops = 50,
                              nStrata = num_strata,
                              mut_chance = 0.3)
    } else {
      solution <- optimStrata(method = "continuous",
                              errors = error, 
                              framesamp = frame,
                              iter = 2000,
                              pops = 50,
                              nStrata = num_strata,
                              mut_chance = 0.3,
                              suggestions = initial_solution3)
    }
    
    strataStructure <- summaryStrata(solution$framenew,
                                     solution$aggr_strata,
                                     progress = FALSE)
    
    n <- as.vector(strataStructure$Allocation)
    df_sol <- solution$framenew
    df_sol$strata <- df_sol$STRATO
    
    ## CV and other objectives
    cv <- calculate_cv(df = df_sol, strata = df_sol$strata, vars = vars, n = n)
    
    trace <- calculate_variance(df = df_ws,
                                strata = df_sol$strata,
                                vars = vars,
                                n = n,
                                objective = 'trace')
    
    det <- calculate_variance(df = df_ws,
                              strata = df_sol$strata,
                              vars = vars,
                              n = n,
                              objective = 'determinant')
    
    varcov <- calculate_variance(df = df_ws,
                                 strata = df_sol$strata,
                                 vars = vars,
                                 n = n,
                                 objective = 'varcov')
    
    pca <- calculate_variance(df = df_ws,
                              strata = df_sol$strata,
                              vars = vars,
                              n = n,
                              objective = 'pca')
    
    # Return all results as a list
    list(
      n = n,
      cv = cv,
      trace = trace,
      det = det,
      varcov = varcov,
      pca = pca,
      objective = solution$objective,
      strata = df_sol$strata
    )
  } # end foreach
  
  # ------------------------ Store results ------------------------
  store_n_all[[datasets[dta]]]      <- do.call(rbind, lapply(results, `[[`, "n"))
  store_cv_all[[datasets[dta]]]     <- do.call(rbind, lapply(results, `[[`, "cv"))
  store_trace_all[[datasets[dta]]]  <- do.call(rbind, lapply(results, `[[`, "trace"))
  store_det_all[[datasets[dta]]]    <- do.call(rbind, lapply(results, `[[`, "det"))
  store_varcov_all[[datasets[dta]]] <- do.call(rbind, lapply(results, `[[`, "varcov"))
  store_pca_all[[datasets[dta]]]    <- do.call(rbind, lapply(results, `[[`, "pca"))
  store_strata_all[[datasets[dta]]] <- lapply(results, `[[`, "strata")
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
  store_det_all = store_det_all,
  store_varcov_all = store_varcov_all,
  store_pca_all = store_pca_all
)

filename <- paste0("OUTPUT/ballbar_", datasets[dta], "_results.Rdata")
save(store, file = filename)