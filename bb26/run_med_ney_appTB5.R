## Run MClust

library(dplyr)
library(haven)
library(cluster)
############# source functions #############
source("calculate_cv.R")
source("function_calc_variance.R")
source("function_sample_size.R")

#------------------------ Section 1: Datasets ---------------------------------

## TB5
load("multivar_datasets/tb5.RData")
tb5 <- df[2:7]


datasets <- c("tb5")
dfs <- list(tb5)
sample <- c(500)
num_strat <- c(3) 

for (dta in 1:length(datasets)) {
  ssize = sample[dta]
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
  type <- c("euclidean")
  
  for (i in 1:length(type)) {
    print(type[i])
    for (iter in 1:40) {
      cat("Simulation Number:", print(iter), "\n")
      set.seed(iter)
      
      start.time <- Sys.time()
      result <- pam(data, k=num_strata, metric=type[i], nstart=25)
      end.time <- Sys.time()
      time <- difftime(end.time, start.time, units="mins")
      
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
    
    filename = paste0("OUTPUT/med_ney_app", datasets[dta], "_results.Rdata")
    
    save(store, file = filename)
    
  }
  
}
