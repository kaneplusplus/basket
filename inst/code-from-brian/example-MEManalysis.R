
source("MEM-singleArm-Functions.R")
source("MEM-FullBayes-function.R")

seed <- round(runif(1,1,1000000000))
set.seed(seed)

########################################
## Define Data and EB hyperparameters ##
########################################
N <- rep(25, 8)

## Generate responses for N baskets with prob vector prob
prob <- rep(0.15, 8) 

## Vectors of Observed number of Responses (X) and Patients (N)
Data <- list(X = rbinom(N, N, prob), N = N ) 
pars <- list(UB=1, 						## MMLE upper boundary
             LB=0) 						## MMLE lower boundary
p0 <- 0.25 # Null response rate for Posterior Probability Calculation

## Emp Bayes optimal MMLE MEM ##

marg.M <- MEM_marginalLoop(Data$X, Data$N, rep(0.5, length(Data$X)), 
                           rep(0.5, length(Data$X)), 
                           INITIAL = diag(1, length(Data$X)))
marg.M$maximizer
round(MEM.cdf(Data,pars,p0,marg.M),3)
round(post.ESS(Data,pars,marg.M),3)
t(round(post.HPD(Data,pars,marg.M,0.05),3))

## constained Emp Bayes optimal ##
pars$UB <- 0.5
marg.M <- MEM_marginalLoop(Data$X, Data$N, rep(0.5,length(Data$X)), 
                           rep(0.5, length(Data$X)), 
                           INITIAL = diag(1,length(Data$X)))
marg.M$maximizer

round(MEM.cdf(Data,pars,p0,marg.M),3)
round(post.ESS(Data,pars,marg.M),3)
t(round(post.HPD(Data,pars,marg.M,0.05),3))

## Full Bayesian with Reference Prior ##
pr.Inclus <- matrix(0.5, length(Data$X), length(Data$X))
diag(pr.Inclus) <- 1
FullBayes.anal <- MEM_FullBayes(Data$X, Data$N, rep(0.5, length(Data$X)), 
                                rep(0.5, length(Data$X)), pr.Inclus, 0.05)
