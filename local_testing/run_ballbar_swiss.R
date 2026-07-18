## Run File - Local Test

library(dplyr)
library(haven)
library(SamplingStrata)
library(cluster)

setwd("C:/Users/01459189/OneDrive/phd/multivariate_stratification")
source("calculate_cv_nofpc.R")
source("function_calc_variance.R")
source("function_sample_size.R")


library(SamplingStrata)
data(swissmunicipalities)
swissmun <- swissmunicipalities[swissmunicipalities$REG ==1,
                                c("REG","COM","Nom","HApoly",
                                  "Surfacesbois","Surfacescult",
                                  "Airbat","POPTOT")]
head(swissmun)

frame3 <- buildFrameDF(df = swissmun,
                       id = "COM",
                       X = c("POPTOT","HApoly"),
                       Y = c("POPTOT","HApoly"),
                       domainvalue = "REG")

ndom <- length(unique(swissmun$REG))

cv <- as.data.frame(list(DOM=rep("DOM1",ndom),
                         CV1=rep(0.05340251,ndom),
                         CV2=rep(0.09218627,ndom),
                         domainvalue=c(1:ndom) ))

set.seed(1234)
init_sol3 <- KmeansSolution2(frame=frame3,
                             errors=cv,
                             maxclusters = 10)  

nstrata3 <- tapply(init_sol3$suggestions,
                   init_sol3$domainvalue,
                   FUN=function(x) length(unique(x)))
nstrata3

initial_solution3 <- prepareSuggestion(init_sol3,frame3,nstrata3)

set.seed(1234)
solution3 <- optimStrata(method = "continuous",
                         errors = cv, 
                         framesamp = frame3,
                         iter = 50,
                         pops = 10,
                         nStrata = nstrata3,
                         suggestions = initial_solution3)

## Domain 1 = 46.93

strataStructure <- summaryStrata(solution3$framenew,
                                 solution3$aggr_strata,
                                 progress=FALSE)
head(strataStructure)

n <- as.vector(strataStructure$Allocation)
n

df_sol <- solution3$framenew
df_sol$strata <- df_sol$STRATO

vars = c("Y1", "Y2")

cv_val <- calculate_cv(df = df_sol, strata = df_sol$strata, vars = vars, n = n)
cv_val
############################### PROPOSED METHOD ##
ssize = 34

data = df_sol[,2:3]
result <- pam(data, k=10, metric="euclidean", nstart=1)

df <- data.frame(data, result$clustering)
vars = c("X1", "X2")

names(df) <- c(vars, "strata")
strata <- df$strata

n <- calculate_sample_size(df=df, strata=df$strata, vars=vars, method='neyman', ssize=ssize)
n

df_sol$strata = result$clustering
vars = c("Y1", "Y2")
cv <- calculate_cv(df=df_sol, strata=df_sol$strata, vars=vars, n=n)
cv