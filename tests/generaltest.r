#library(testthat)
library(basket)
library(igraph)
#library(foreach)


data(vemu_wide)

baskets <- 1:5

vemu_wide1 <- vemu_wide[baskets,]

# Full Bayes
exactRes <- mem_full_bayes_exact(responses = vemu_wide1$responders, 
                     size = vemu_wide1$evaluable,
                     name = vemu_wide1$baskets, p0=rep(0.15, length(basket)))



print(exactRes$call)

print(exactRes$basketwise$PEP)
print(exactRes$basketwise$MAP)
print(exactRes$basketwise$CDF)
print(exactRes$basketwise$HPD)
print(exactRes$basketwise$ESS)
print(exactRes$basketwise$mean_est)
print(exactRes$basketwise$median_est)


print(exactRes$clusterwise$cluster)
print(exactRes$clusterwise$CDF)
print(exactRes$clusterwise$HPD)
print(exactRes$clusterwise$ESS)
print(exactRes$clusterwise$mean_est)
print(exactRes$clusterwise$median_est)

OCTable(exactRes$basketwise)
OCTable(exactRes$clusterwise)


exactResNew <- updateResult(exactRes, 0.25)

OCTable(exactRes$basketwise)
OCTable(exactResNew$basketwise)

OCTable(exactRes$clusterwise)
OCTable(exactResNew$clusterwise)

plot_density(exactRes$basketwise)
plot_posterior_exchangeability(exactRes$basketwise)
plot_all_exchangeability(exactRes$basketwise, c("PRIOR", "MAP", "PEP"))
plot_all_exchangeability(exactRes$basketwise, c("MAP", "PEP"))

#rm(list=ls())
## Vectors of Observed number of Responses (X) and Patients (N)
##Data <- dget("../inst/code-from-brian/MHAlgorithm/Vemu-data.txt") 

baskets <- 1:6
vemu_wide2 <- vemu_wide[baskets,]

MHResult1 <- mem_full_bayes_mcmc(responses = vemu_wide2$responders, 
                                 size = vemu_wide2$evaluable,
                                 name=c("NSCLC ","CRC.v ","CRC.vc","  BD  ","ED.LH "," ATC  "),
                                 p0 = c(0.15, 0.15, 0.15, 0.2, 0.15, 0.15), Initial = NA)


print(MHResult1$call)

print(MHResult1$basketwise$PEP)
print(MHResult1$basketwise$MAP)
print(MHResult1$basketwise$CDF)
print(MHResult1$basketwise$HPD)
print(MHResult1$basketwise$ESS)
print(MHResult1$basketwise$ESS2)
print(MHResult1$basketwise$mean_est)
print(MHResult1$basketwise$median_est)


print(MHResult1$clusterwise$cluster)
print(MHResult1$clusterwise$CDF)
print(MHResult1$clusterwise$HPD)
print(MHResult1$clusterwise$ESS)
print(MHResult1$clusterwise$ESS2)
print(MHResult1$clusterwise$mean_est)
print(MHResult1$clusterwise$median_est)

OCTable(MHResult1$basketwise)
OCTable(MHResult1$clusterwise)

MHResult1New <- updateResult(MHResult1, 0.25)

OCTable(MHResult1$basketwise)
OCTable(MHResult1New$basketwise)

OCTable(MHResult1$clusterwise)
OCTable(MHResult1New$clusterwise)

