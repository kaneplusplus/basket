
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

### MCMC ###
source("mem-updates.R")

fit.MCMC.1 <- mem_full_bayes_mcmc(responses=Data$X, size=Data$N, name=c("NSCLC ","CRC.v ","CRC.vc","  BD  ","ED.LH "," ATC  "),
                                null = 0.15, Initial = NA)

##################################################
###### Nan: Updated 1/6 for MEM clustering #######
##################################################
fit.MCMC.1$PEP
library(igraph)
g <- igraph::graph_from_adjacency_matrix(fit.MCMC.1$PEP, mode = "undirected", weighted = TRUE, diag = FALSE)
A <- factor( cluster_louvain(g, weights=E(g)$weight)$membership )

for(i in 1:length(levels(A))){
 print(paste0("#### Cluster Assignment: ",levels(A)[i]))
 print( fit.MCMC.1$name[which(A==levels(A)[i])] )
}





