## Run MClust

library(dplyr)
library(cluster)
library(haven)
par(family="serif", bty="L", mfrow = c(1, 1))

############# source functions #############

source("function_calc_variance.R")
source("calculate_cv.R")
source("function_sample_size.R")

#------------------------ Section 1: Datasets ---------------------------------


## GHS
load("multivar_datasets/ghs_2024.RData")
df$hhinc_pc <- df$totmhinc/df$hholdsz
df <- df[df$fin_reqinc!=9999999,]
ghs <- df[,-4]


datasets <- c("ghs")
dfs <- list(ghs)
sample <- c(2000)
num_strata <- c(5) 


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


#------------------------ Section 3: Run mclust ----------------------------------
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
type <- c('neyman', 'prop')

for (i in 1:length(type)) {
print(type[i])
  for (iter in 1:10) {
    cat("Simulation Number:", print(iter), "\n")
set.seed(iter)

start.time <- Sys.time()
m <- kmeans(x=data, centers=num_strata, iter.max = 200, nstart = 25)
end.time <- Sys.time()
time <- difftime(end.time, start.time, units="mins")

df <- data.frame(data, m$cluster)
names(df) <- c(vars, "strata")
strata <- df$strata

# ----------------------- Section 4: Report Results --------------------------
##########

n <- calculate_sample_size(df=df, strata=df$strata, vars=vars, method=type[i], ssize=ssize)
print(sum(n))
cv <- calculate_cv(df=df_ns, strata=df$strata, vars=vars, n=n)

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
}
store <- list(df_ns, df, store_n, store_strata, store_cv, store_trace, store_det, store_varcov, store_pca, store_obj, store_time)

filename = paste0("OUTPUT/kmeans_", datasets[dta], "_results.Rdata")

save(store, file = filename)

}
