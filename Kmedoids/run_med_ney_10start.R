## Run K-medoids

library(dplyr)
library(haven)
library(cluster)
############# source functions #############
source("calculate_cv.R")
source("function_calc_variance.R")
source("function_sample_size.R")

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
  num_strata <- num_strat[dta]
  
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
store_n <- numeric()
store_cv <- numeric()
store_trace <- numeric()
store_det <- numeric()
store_varcov <- numeric()
store_pca <- numeric()
store_obj <- numeric()
store_strata <- list()
store_time <- numeric()
type <- c("manhattan", "euclidean")

for (i in 1:length(type)) {
  print(type[i])
  for (iter in 1:40) {
    cat("Simulation Number:", print(iter), "\n")
    set.seed(iter)
  
  start.time <- Sys.time()
  result <- pam(data, k=num_strata, metric=type[i], nstart=10)
  end.time <- Sys.time()
  time <- difftime(end.time, start.time, units="mins")
  print("just finished pam")
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
  
  det <- calculate_variance(df=df,
                            strata=df$strata,
                            vars=vars,
                            n=n,
                            objective='determinant')
  
  varcov <- calculate_variance(df=df,
                               strata=df$strata,
                               vars=vars,
                               n=n,
                               objective='varcov')
  
  pca <- calculate_variance(df=df,
                            strata=df$strata,
                            vars=vars,
                            n=n,
                            objective='pca')
  
  store_n <- rbind(store_n, n)
  store_cv <- rbind(store_cv, cv)
  store_trace <- rbind(store_trace, trace)
  store_det <- rbind(store_det, det)
  store_varcov <- rbind(store_varcov, varcov)
  store_pca <- rbind(store_pca, pca)
  store_obj <- rbind(store_obj, type[i])
  store_strata[[iter]] <- strata
  store_time <- rbind(store_time, time)
}
store <- list(df_ns, df, store_n, store_strata, store_cv, store_trace, store_det, store_varcov, store_pca, store_obj, store_time)

filename = paste0("OUTPUT/med_ney", datasets[dta], "_results_10start.Rdata")

save(store, file = filename)

}

}

