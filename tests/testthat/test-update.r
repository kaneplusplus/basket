library(testthat)
library(basket)

context("Update models")

data(vemu_wide)


mh1 <- mem_mcmc(responses = vemu_wide$responders, 
                size = vemu_wide$evaluable,
                name= vemu_wide$baskets,
                p0 = c(0.15, 0.15, 0.15, 0.2, 0.15, 0.15),
                mcmc_iter = 100)

mh2 <- update_p0(mh1, 0.25)

