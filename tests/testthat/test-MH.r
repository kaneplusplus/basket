library(testthat)
library(basket)
library(foreach)
foreach::registerDoSEQ()

context("Metropolis-Hasting test")

# Load the reference outputs.
load(file = "../../Data/fitMCMC1_reference.rda")
load(file = "../../Data/fitMCMC2_reference.rda")


## Vectors of Observed number of Responses (X) and Patients (N)
Data <- dget("../../inst/code-from-brian/MHAlgorithm/Vemu-data.txt") 

MHResult1 <- mem_full_bayes_mcmc(responses=Data$X, size=Data$N, 
                                 name=c("NSCLC ","CRC.v ","CRC.vc","  BD  ","ED.LH "," ATC  "),
                                 p0 = 0.15, Initial = NA)


expect_equal(MHResult1$PEP, fit.MCMC.1$PEP)
expect_equal(MHResult1$HPD, fit.MCMC.1$HPD, tolerance = 5e-2)
expect_equal(MHResult1$CDF, fit.MCMC.1$CDF)
expect_equal(MHResult1$ESS, fit.MCMC.1$ESS, tolerance = 5e-2)
expect_equal(MHResult1$ESS2, fit.MCMC.1$ESS2, tolerance = 5e-2)



init <- dget("../../inst/code-from-brian/MHAlgorithm/MCMCtestInitial.txt")

MHResult2 <- mem_full_bayes_mcmc(responses=Data$X, size=Data$N, name=c("NSCLC ","CRC.v ","CRC.vc","  BD  ","ED.LH "," ATC  "),
                                 p0 = 0.15, Initial = init)


expect_equal(MHResult2$PEP, fit.MCMC.2$PEP)
expect_equal(MHResult2$HPD, fit.MCMC.2$HPD, tolerance = 5e-2)
expect_equal(MHResult2$CDF, fit.MCMC.2$CDF)
expect_equal(MHResult2$ESS, fit.MCMC.2$ESS, tolerance = 5e-2)
expect_equal(MHResult2$ESS2, fit.MCMC.2$ESS2, tolerance = 5e-2)



