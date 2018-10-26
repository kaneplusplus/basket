
rm(list=ls())

library(foreach)
library(doMC)
registerDoMC()

### Seq Vermu Sim ###
source("../Utils/parSim-functions.R")

########################################
## 1) Define Data and hyperparameters ##
########################################
Data <- dget("Vemu-data.txt") ## Vectors of Observed number of Responses (X) and Patients (N)
p0 <- 0.15				## Null response rate for Posterior Probability Calculation

pr.Inclus <- matrix(0.5, length(Data$X), length(Data$X)); diag(pr.Inclus) <- 1 
alp=0.05; N=Data$N; p0=p0


prob=c(0.45,0.45,0.45,0.45,0.45,0.45)

 for(H in 3619:5000){
 # no_cores <- detectCores() - 1;  cl <- makeCluster(no_cores); registerDoParallel(cl)
 #  foreach( ITER = ( (H-1)*no_cores ):((H-1)*no_cores + no_cores), .combine = list, 
 #           .multicombine = TRUE)  %dopar% 
 #   runSim(ITER, N, prob, pr.Inclus, p0, alp)
    runSim(H, N, prob, pr.Inclus, p0, alp)
  #stopCluster(cl)
 }

