## Calculate Variance Function

calculate_variance <- function(df,   ## data frame with y variables
                               strata, ## data frame with strata variable
                               vars,   ## vector of column names for y variables
                               n,  ## sample allocation across strata 
                               objective)      ## objective function (determinant of var-cov, total variance, sum of squares)
{
  #print(n)
  df$strata <- strata
  dfg <- df %>%
    group_by(strata) %>%
    summarise(across(all_of(vars), list(mean = mean, sd = sd)),
              N_h = n(), .groups = "drop")
  
  N <- dfg$N_h 
  W <- dfg$N_h / sum(N) ## population proportion of stratum
  L <- nrow(dfg) ## number of strata
  k <- length(vars) ## number of variables
  
  S <- dfg %>% dplyr::select(dplyr::ends_with("_sd")) %>% as.matrix()
  
  fitness <- matrix(0, nrow = L, ncol = k)
  results <- numeric(k)
  
  if (objective == 'trace') {
  
  for (i in 1:k) {
    Vk <- S[, i]^2 ## Variance for variable i
    
    for (h in 1:L) {
      fitness[h, i] <- (W[h]^2 * (1 - (n[h] / N[h])) * (Vk[h] / n[h])) ## Kozak (5)
    }
    
    results[i] <- sum(fitness[, i])
  }
  
  #print(dfg)
  total_variance <- sum(results) ## SUM VARIANCE ACROSS K FOR MULTIVARIATE
  #print(results)
  
  return(total_variance)
  
  }
  
  if (objective == 'determinant') {
    
    varcov = list(0)
    for(h in 1:L){
      varcov[[h]]=cov(df[df$strata==h,vars])
    }
    varcov
    
    weighted1 <- Map(`*`, varcov, W^2/n) 
    weighted1.sum <- Reduce("+", weighted1)
    
    weighted2 <- Map(`*`,varcov, W)
    weighted2.sum <- Reduce("+", weighted2)/sum(N)
    matrix.optim = weighted1.sum-weighted2.sum
    
    determinant= det(matrix.optim)
    
    return(determinant)
  }
  
  if (objective == 'varcov') {
    
    varcov = list(0)
    for(h in 1:L){
      varcov[[h]]=cov(df[df$strata==h,vars])
    }
    varcov
    
    weighted1 <- Map(`*`, varcov, W^2/n) 
    weighted1.sum <- Reduce("+", weighted1)
    
    weighted2 <- Map(`*`,varcov, W)
    weighted2.sum <- Reduce("+", weighted2)/sum(N)
    matrix.optim = weighted1.sum-weighted2.sum
    
    sum <- Reduce("+", matrix.optim)
    
    return(sum)
  }
  
  ## weighted sum / trace
  if (objective == 'pca') {
    
    for (i in 1:k) {
      Vk <- S[, i]^2 ## Standard deviation for variable i
      
      for (h in 1:L) {
        fitness[h, i] <- (W[h]^2 * (1 - (n[h] / N[h])) * (Vk[h] / n[h])) ## Kozak (5)
      }
      
      results[i] <- sum(fitness[, i])
    }
    
    
    x <- tibble(df[,vars])
    
    ### PCA ###############################################################
    pca <- princomp(x, cor=TRUE)
    load <- pca$loadings ## loadings from comp 1
    w <- abs(load[,1]/sum(abs(load[,1]))) ## weights from loadings, not direction specific
    
    total_variance <- w %*% results ## weighted fitness ## SUM VARIANCE ACROSS K FOR MULTIVARIATE
    #print(results)
    
    return(total_variance)
    
  }
  
  ## between and within SS
  if (objective == 'ss') {

    grand_mean <- colMeans(df[,vars])
    
    within_ss <- 0
    between_ss <- 0
    
    for (h in L) {
      mean <- colMeans(df[df$strata==h,vars])
      within_ss <- within_ss + colSums((df[df$strata==h,vars] - mean)^2)
      between_ss <- between_ss + n[h] * (mean - grand_mean)^2
    }
    
    ratio <- sum(between_ss/within_ss)
    
    return(ratio)
  }
}
