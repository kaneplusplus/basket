library(testthat)
library(basket)
library(igraph)
library(foreach)


## Vectors of Observed number of Responses (X) and Patients (N)
Data <- dget("../inst/code-from-brian/MHAlgorithm/Vemu-data.txt") 

MHResult1 <- mem_full_bayes_mcmc(responses=Data$X, size=Data$N, 
                                 name=c("NSCLC ","CRC.v ","CRC.vc","  BD  ","ED.LH "," ATC  "),
                                 p0 = 0.15, Initial = NA)


print(MHResult1$PEP)
print(MHResult1$HPD)
print(MHResult1$ESS)
print(MHResult1$ESS2)

source("pep.r")
result <- cluster_PEP(responses=Data$X, size=Data$N, 
            name=c("NSCLC ","CRC.v ","CRC.vc","  BD  ","ED.LH "," ATC  "),
            p0 = 0.15, models=MHResult1$models, pweights=MHResult1$pweights,
            PEP=MHResult1$PEP)

print(result$clusters)
print(result$HPD)
print(result$mean_est)
print(result$median_est)
print(result$ESS)
print(result$ESS2)
plot(density(result$samples[,1]))
plot(density(result$samples[,2]))







