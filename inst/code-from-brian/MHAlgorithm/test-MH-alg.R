
library(tibble)
library(knitr)
library(basket)
library(readr)
library(ggplot2)
library(foreach)
foreach::registerDoSEQ()

###################################
### Vemu Study Data ###############
###################################

Data <- dget("Vemu-data.txt") ## Vectors of Observed number of Responses (X) and Patients (N)

### Exact Computation ##
fit <- mem_full_bayes(responses=Data$X, size=Data$N, name=c("NSCLC ","CRC.v ","CRC.vc","  BD  ","ED.LH "," ATC  "),
                      p0=0.15 )

### MCMC ###
source("mem-updates.R")

fit.MCMC.1 <- mem_full_bayes_mcmc(responses=Data$X, size=Data$N, name=c("NSCLC ","CRC.v ","CRC.vc","  BD  ","ED.LH "," ATC  "),
                                null = 0.15, Initial = NA)

init <- fit$maximizer
fit.MCMC.2 <- mem_full_bayes_mcmc(responses=Data$X, size=Data$N, name=c("NSCLC ","CRC.v ","CRC.vc","  BD  ","ED.LH "," ATC  "),
                                null = 0.15, Initial = init)

print(fit$CDF)
print(fit.MCMC.1$PostProb)
print(fit.MCMC.2$PostProb)

print(fit$PEP)
print(fit.MCMC.1$PEP)
print(fit.MCMC.2$PEP)

print(fit$HPD)
print(fit.MCMC.1$HPD)
print(fit.MCMC.2$HPD)

print(fit$ESS)
print(fit.MCMC.1$ESS)
print(fit.MCMC.1$ESS2)

print(fit.MCMC.2$ESS)
print(fit.MCMC.2$ESS2)

