

library(basket)
#library(foreach)
#foreach::registerDoSEQ()


trial_sizes <- rep(5, 5)

# The response rates for the baskets.
resp_rate <- 0.15

trials <-
  data.frame(responses = rbinom(trial_sizes, trial_sizes, resp_rate),
             size = trial_sizes)
mem_full_bayes_mcmc(trials$responses,
                    trials$size,
                    name = c("  T1  ", "  T2  ", "  T3  ", "  T4  ", "  T5  "))



# Trial test
Data <- dget("Vemu-data.txt") ## Vectors of Observed number of Responses (X) and Patients (N)

MHResult1 <- mem_full_bayes_mcmc(responses=Data$X, size=Data$N, name=c("NSCLC ","CRC.v ","CRC.vc","  BD  ","ED.LH "," ATC  "),
                                p0 = 0.15, Initial = NA)
print(MHResult1$ESS)
print(MHResult1$ESS2)


init <- dget("MCMCtestInitial.txt")
MHResult2 <- mem_full_bayes_mcmc(responses=Data$X, size=Data$N, name=c("NSCLC ","CRC.v ","CRC.vc","  BD  ","ED.LH "," ATC  "),
                                p0 = 0.15, Initial = init)

print(MHResult1$PostProb)
print(MHResult2$PostProb)


print(MHResult1$PEP)
print(MHResult2$PEP)


print(MHResult1$HPD)
print(MHResult2$HPD)


print(MHResult1$ESS)
print(MHResult1$ESS2)

print(MHResult2$ESS)
print(MHResult2$ESS2)

