## Sample size

calculate_sample_size <- function(df,   ## data frame with y variables
                               strata, ## data frame with strata variable
                               vars,   ## vector of column names for y variables
                               method, ## allocation method
                               ssize)  ## sample size
{
  df$strata <- strata
  dfg <- df %>%
    group_by(strata) %>%
    summarise(across(all_of(vars), list(mean = mean, sd = sd)),
              N_h = n(), .groups = "drop")
  
  N <- dfg$N_h 
  W <- dfg$N_h / sum(N) ## population proportion of stratum
  L <- nrow(dfg) ## number of strata
  k <- length(vars) ## number of variables
  
  S <- dfg %>% dplyr::select(ends_with("_sd")) %>% as.matrix()

  if (method == 'prop') {
  
    n <- numeric()
    for (h in 1:L){
      n[h] <- as.integer(ssize*(N[h]/sum(N)))
    }
    
    for (h in 1:L){
      if (n[h] %in% c(0,1)) { ## ensure smallest stratum is bigger than 1 
        n[h] = 2
        n[which.max(n)] <- ssize-sum(n)+n[which.max(n)]
      }
    }
    
    z <- which.max(N)
    n[z] <- n[z] + (ssize-sum(n)) ## allocate any remaining sample to largest stratum
  }
  
  if (method == 'neyman') {
    
    V <- S^2
    n <- numeric()
    
    for (h in 1:L){ 
      
      Vh <- V[h, ] ## Variance for stratum h 
      
      n[h] <- as.integer(ssize*((N[h]*sqrt(sum(Vh))))/(sum(N*sqrt(rowSums(V))))) 
    }
    
    z <- which.max(rowSums(V))
    r <- which.max(N)
    
    ## Constraint 1 check
    
    for (h in 1:L){
      if (n[h] %in% c(0,1)) { ## ensure smallest stratum is bigger than 1 
        n[h] = 2
        n[which.max(n)] <- ssize-sum(n)+n[which.max(n)]
      }
    }
    
    ## Constraint 2 check
    
    for (h in 1:L){
      if (n[h] >= N[h]) { ## ensure no stratum sample is bigger than stratum size
        n[h] = N[h]
        n[which.max(N)] <- ssize-sum(n)+n[which.max(N)]
      }
    }
    
    if (n[z] + (ssize-sum(n)) <= N[z]) {
      n[z] <- n[z] + (ssize-sum(n))
      } ## allocating remaining sample to strata with most variance
    
    if (n[z] + (ssize-sum(n)) > N[z]) {
      n[r] <- n[r] + (ssize-sum(n))
      } ## allocating remaining sample to largest stratum
    
  }
  
  if (method == 'neyman outlier') {
    
    V <- S^2
    n <- numeric()
    n[L] <- N[L] ## NOTE: In most of our cases, this will be too big to satisfy sample size requirements.
    
    for (h in 1:L){ 
      
      Vh <- V[h, ] ## Variance for stratum h 
      
      n[h] <- as.integer(ssize*((N[h]*sqrt(sum(Vh))))/(sum(N*sqrt(rowSums(V))))) 
    }
    
    z <- which.max(rowSums(V))
    r <- which.max(N)
    
    ## Constraint 1 check
    
    for (h in 1:L){
      if (n[h] %in% c(0,1)) { ## ensure smallest stratum is bigger than 1 
        n[h] = 2
        n[which.max(n)] <- ssize-sum(n)+n[which.max(n)]
      }
    }
    
    ## Constraint 2 check
    
    for (h in 1:L){
      if (n[h] >= N[h]) { ## ensure no stratum sample is bigger than stratum size
        n[h] = N[h]
        n[which.max(N)] <- ssize-sum(n)+n[which.max(N)]
      }
    }
    
    if (n[z] + (ssize-sum(n)) <= N[z]) {
      n[z] <- n[z] + (ssize-sum(n))
    } ## allocating remaining sample to strata with most variance
    
    if (n[z] + (ssize-sum(n)) > N[z]) {
      n[r] <- n[r] + (ssize-sum(n))
    } ## allocating remaining sample to largest stratum
    
  }
  
  return(n)
}

#calculate_sample_size(df, eg, vars, 'prop', 100)
#calculate_sample_size(df, eg, vars, 'neyman', 100)