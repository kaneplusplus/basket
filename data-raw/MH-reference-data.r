library(tibble)
library(knitr)
library(basket)
library(readr)
library(ggplot2)
library(foreach)
foreach::registerDoSEQ()

source("..\\inst\\code-from-brian\\MHAlgorithm\\mem-updates.R")

Data <- dget("..\\inst\\code-from-brian\\MHAlgorithm\\Vemu-data.txt")


fit.MCMC.1 <- mem_full_bayes_mcmc(
  responses = Data$X, size = Data$N, name = c("NSCLC ", "CRC.v ", "CRC.vc", "  BD  ", "ED.LH ", " ATC  "),
  null = 0.15, Initial = NA
)
save(fit.MCMC.1, file = "../data/fitMCMC1_reference.rda")


init <- dget("..\\inst\\code-from-brian\\MHAlgorithm\\MCMCtestInitial.txt")
fit.MCMC.2 <- mem_full_bayes_mcmc(
  responses = Data$X, size = Data$N, name = c("NSCLC ", "CRC.v ", "CRC.vc", "  BD  ", "ED.LH ", " ATC  "),
  null = 0.15, Initial = init
)
save(fit.MCMC.2, file = "../data/fitMCMC2_reference.rda")
