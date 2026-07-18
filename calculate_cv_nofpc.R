calculate_cv <- function(df,   ## data frame with y variables
                         strata, ## data frame with strata variable
                         vars,   ## vector of column names for y variables
                         n)      ## sample allocation across strata 
{
  df$strata <- strata
  df <- df %>%
    group_by(strata) %>%
    summarise(across(all_of(vars), list(mean = mean, sd = sd)),
              N_h = n(), .groups = "drop")
  
  N <- df$N_h 
  W <- df$N_h / sum(N) ## population proportion of stratum
  L <- nrow(df) ## number of strata
  k <- length(vars) ## number of variables
  
  S <- df %>% dplyr::select(ends_with("_sd")) %>% as.matrix()
  M <- df %>% dplyr::select(ends_with("_mean")) %>% as.matrix()
  
  fitness <- matrix(0, nrow = L, ncol = k)
  results <- numeric(k)
  Mk <- matrix(0, nrow = L, ncol = k)
  meanst <- numeric(k)
  for (i in 1:k) {
    Vk <- S[, i]^2 ## Variance for variable i
    #Mk <- M[, i] ## Mean for variable i
    for (h in 1:L) {
      
      if (n[h]!= N[h])
      fitness[h, i] <- (W[h]^2 * (Vk[h] / n[h])) ## Cochran's variance formula
      print(fitness[h, i])
      
      if (n[h]==N[h]) {
        fitness[h, i] = 0 ## take-all stratum
      }
      
      Mk[h,i] <- W[h]*M[h,i]
    }
    

    
    variance_i <- sum(fitness[, i])
    results[i] <- sqrt(variance_i) / sum(Mk[,i]) ## Compute CV
    meanst[i] <- sum(Mk[,i]) 
    print(results[i])
  }
  print("Summary of df for cv calculation (should be unscaled):")
  print(df)
  total_cv <- sum(results) ## SUM CV ACROSS K FOR MULTIVARIATE
  # return(list(
  #   results=results,
  #   mst = meanst))
  
  return(results)
}

