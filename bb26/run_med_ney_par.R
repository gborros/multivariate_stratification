## Run MClust

library(dplyr)
library(haven)
library(cluster)
library(doParallel)
############# source functions #############
source("calculate_cv.R")
source("function_calc_variance.R")
source("function_sample_size.R")

ncores <- 40
cl <- makeCluster(ncores)
clusterEvalQ(cl, {
  library(dplyr)
  library(haven)
  library(cluster)
  source("function_sample_size.R")
  source("calculate_cv.R")
  source("function_calc_variance.R")
})
registerDoParallel(cl)

#------------------------ Section 1: Datasets ---------------------------------

## NORM
load("multivar_datasets/d_norm.RData")
d_norm <- as.data.frame(d_norm)

## NORM LOW
load("multivar_datasets/d_norm_low.RData")
d_norm_low <- as.data.frame(d_norm_low)

## LOGNORM
load("multivar_datasets/d_lognorm.RData")
d_lognorm <- as.data.frame(d_lognorm)

## LOGNORM LOW
load("multivar_datasets/d_lognorm_low.RData")
d_lognorm_low <- as.data.frame(d_lognorm_low)

## CHISQ
load("multivar_datasets/d_chisq.RData")
d_chisq <- as.data.frame(d_chisq)

## COMB1
load("multivar_datasets/d_comb1.RData")
d_comb1 <- as.data.frame(d_comb1)

## COMB2
load("multivar_datasets/d_comb2.RData")
d_comb2 <- as.data.frame(d_comb2)

## COMB3
load("multivar_datasets/d_comb3.RData")
d_comb3 <- as.data.frame(d_comb3)

## COMB4
load("multivar_datasets/d_comb4.RData")
d_comb4 <- as.data.frame(d_comb4)

## COMB5
load("multivar_datasets/d_comb5.RData")
d_comb5 <- as.data.frame(d_comb5)

datasets <- c("d_norm", "d_norm_low", "d_lognorm", "d_lognorm_low",
              "d_chisq", "d_comb1", "d_comb2", "d_comb3", "d_comb4",
              "d_comb5")
dfs <- list(d_norm, d_norm_low, d_lognorm, d_lognorm_low, d_chisq,
            d_comb1, d_comb2, d_comb3, d_comb4, d_comb5)

sample <- rep(500,10)
num_strat <- rep(5, 10)

for (dta in 1:length(datasets)) {
  ssize = sample[1]
  vars <- character()
  sf <- dfs[[dta]]
  
  for (z in 1:ncol(sf)) {
    
    var <- paste0("X", z)
    vars <- c(vars, var)
    
  }
  
  colnames(sf) <- vars
#------------------------ Section 2: Dataset set up --------------------------

set_up <- function(df, sort_var) { ## set up function to ensure scaled and no na data
  df <- as.data.frame(scale(na.omit(df)))
  df <- df[order(df[[sort_var]]), ]
  rownames(df) <- 1:nrow(df)
  return(df)
}

df <- set_up(df = sf, sort_var = "X1")

set_up_ns <- function(df, sort_var) { ## set up function with no scaling
  df <- as.data.frame(na.omit(df))
  df <- df[order(df[[sort_var]]), ]
  rownames(df) <- 1:nrow(df)
  return(df)
}

df_ns <- set_up_ns(df = sf, sort_var = "X1")


#------------------------ Section 3: Run Spec ----------------------------------
data <- df
num_strata <- num_strat[dta]

store_n_all <- list()
store_cv_all <- list()
store_trace_all <- list()
store_det_all <- list()
store_varcov_all <- list()
store_pca_all <- list()
store_obj_all <- list()
store_strata_all <- list()
store_time_all <- list()

type <- c("manhattan", "euclidean")

for (i in 1:length(type)) {
  print(type[i])
  results <- foreach(iter = 1:40, .packages = c("dplyr", "haven", "cluster")) %dopar% {
  set.seed(iter)
  result <- pam(data, k=num_strata, metric=type[i], nstart=25)

  df <- data.frame(data, result$clustering)
  names(df) <- c(vars, "strata")
  strata <- df$strata
  df_ns$strata <- strata 

  n <- calculate_sample_size(df=df, strata=df$strata, vars=vars, method='neyman', ssize=ssize)
  print(sum(n))
  cv <- calculate_cv(df=df_ns, strata=df_ns$strata, vars=vars, n=n)
  strata <- as.factor(df$strata)
  
  ########## Calculate other objectives for best solution #################
  
  trace <- calculate_variance(df=df,
                              strata=df$strata,
                              vars=vars,
                              n=n,
                              objective='trace')
  
  list(
    n=n, cv=cv, trace=trace, strata=strata, objective=type[i]
  )
  
  }
  
  store_n_all[[ type[i] ]]      <- do.call(rbind, lapply(results, `[[`, "n"))
  store_cv_all[[ type[i] ]]     <- do.call(rbind, lapply(results, `[[`, "cv"))
  store_trace_all[[ type[i] ]]  <- do.call(rbind, lapply(results, `[[`, "trace"))
  store_obj_all[[ type[i] ]]    <- sapply(results, `[[`, "objective")
  store_strata_all[[ type[i] ]] <- lapply(results, `[[`, "strata")
}

# ----------------------- Section 5: Store Results --------------------------
stopCluster(cl)

store <- list(
  df_ns, df,
  store_n_all,
  store_strata_all,
  store_cv_all,
  store_trace_all,
  store_obj_all
)

filename = paste0("OUTPUT/med_ney", datasets[dta], "_results.Rdata")

save(store, file = filename)

}


