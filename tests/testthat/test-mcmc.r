library(testthat)
library(basket)

context("Fit MCMC models")

data(vemu_wide)

baskets <- c(1, 5, 6)

vemu_wide1 <- vemu_wide[baskets,]

mcmc_ref1 <- readRDS("reference_data/mem-mcmc1.rds")

mcmc_res1 <- mem_mcmc(
  responses = vemu_wide1$responders, 
  size = vemu_wide1$evaluable,
  name = vemu_wide1$baskets,
  p0 = 0.15)

summary(mcmc_res1$basket)
summary(mcmc_res1$cluster)
summary(mcmc_res1)

b1 <- mcmc_ref1$basket

expect_equal(mcmc_res1$basketwise$PEP, mcmc_ref1$PEP)
expect_equal(mcmc_res1$basketwise$HPD, mcmc_ref1$HPD, tolerance = 5e-2)
expect_equal(mcmc_res1$basketwise$CDF, mcmc_ref1$CDF)
expect_equal(mcmc_res1$basketwise$ESS, mcmc_ref1$ESS, tolerance = 5e-2)

