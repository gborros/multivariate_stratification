## ============================================================
## Local testing script: K-medoids stratification, single run
## by: dataset, num_strata, sample_size, seed, metric
## ============================================================

library(dplyr)
library(haven)
library(cluster)

setwd("C:/Users/01459189/OneDrive/phd/multivariate_stratification")
source("calculate_cv_nofpc.R")
source("function_calc_variance.R")
source("function_sample_size.R")

#------------------------ Section 1: Parameters --------------------------

dataset_name    <- "d_norm"
num_strata <- 5
ssize <- 500
seed       <- 1
metric     <- "euclidean"   # "euclidean" or "manhattan"

cat("Running local test with:\n")
cat("  dataset     :", dataset_name, "\n")
cat("  num_strata  :", num_strata, "\n")
cat("  sample_size :", ssize, "\n")
cat("  seed        :", seed, "\n")
cat("  metric      :", metric, "\n\n")

#------------------------ Section 2: Load the single dataset -------------

data_path <- file.path("multivar_datasets", paste0(dataset_name, ".RData"))
if (!file.exists(data_path)) {
  stop("Could not find dataset file: ", data_path)
}
load(data_path)
sf <- as.data.frame(get(dataset_name))

vars <- paste0("X", seq_len(ncol(sf)))
colnames(sf) <- vars

#------------------------ Section 3: Dataset set up -----------------------

set_up <- function(df, sort_var) {
  df <- as.data.frame(scale(na.omit(df)))
  df <- df[order(df[[sort_var]]), ]
  rownames(df) <- 1:nrow(df)
  return(df)
}

set_up_ns <- function(df, sort_var) {
  df <- as.data.frame(na.omit(df))
  df <- df[order(df[[sort_var]]), ]
  rownames(df) <- 1:nrow(df)
  return(df)
}

df    <- set_up(df = sf, sort_var = "X1")
df_ns <- set_up_ns(df = sf, sort_var = "X1")

#------------------------ Section 4: Single run ----------------------------

data <- df

set.seed(seed)

start.time <- Sys.time()
result <- pam(data, k = num_strata, metric = metric, nstart = 25)
end.time <- Sys.time()
run_time <- difftime(end.time, start.time, units = "mins")

cat("pam() finished in", round(as.numeric(run_time), 3), "minutes\n\n")

df_clustered <- data.frame(data, result$clustering)
names(df_clustered) <- c(vars, "strata")

df_ns$strata <- df_clustered$strata

n <- calculate_sample_size(
  df = df_clustered,
  strata = df_clustered$strata,
  vars = vars,
  method = "neyman",
  ssize = ssize
)
cat("Total allocated sample size:", sum(n), "\n")

cv <- calculate_cv(
  df = df_ns,
  strata = df_ns$strata,
  vars = vars,
  n = n
)

strata_factor <- as.factor(df_clustered$strata)

trace <- calculate_variance(df = df_clustered, strata = df_clustered$strata, vars = vars, n = n, objective = "trace")
det   <- calculate_variance(df = df_clustered, strata = df_clustered$strata, vars = vars, n = n, objective = "determinant")
varcov <- calculate_variance(df = df_clustered, strata = df_clustered$strata, vars = vars, n = n, objective = "varcov")
pca   <- calculate_variance(df = df_clustered, strata = df_clustered$strata, vars = vars, n = n, objective = "pca")

#------------------------ Section 5: Console summary -----------------------

cat("\n---- Summary ----\n")
cat("n per stratum:\n"); print(n)
cat("\ncv:\n"); print(cv)
cat("\ntrace:", trace, "\n")
cat("determinant:", det, "\n")
cat("varcov:", varcov, "\n")
cat("pca:", pca, "\n")

