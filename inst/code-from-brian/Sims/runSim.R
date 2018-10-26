library(devtools)
library(foreach)
library(doMC)
registerDoMC()

document()
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
  runSim(H, N, prob, pr.Inclus, p0, alp)
  print(H)
}

