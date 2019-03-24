library(testthat)
library(basket)

context("Fit MCMC models")

data(vemu_wide)

baskets <- c(1, 2, 3)

vemu_wide1 <- vemu_wide[baskets, ]

mcmc_ref1 <- readRDS("reference_data/mem-mcmc1.rds")

mcmc_ref2 <- mcmc_ref1$basketwise

mcmc_res1 <- mem_mcmc(
  responses = vemu_wide1$responders,
  size = vemu_wide1$evaluable,
  name = vemu_wide1$baskets,
  p0 = 0.15
)

summary(mcmc_res1$basket)
summary(mcmc_res1$cluster)
summary(mcmc_res1)


expect_equal(mcmc_res1$basket$PEP, mcmc_ref2$PEP)
expect_equal(mcmc_res1$basket$HPD, mcmc_ref2$HPD, tolerance = 5e-2)
expect_equal(mcmc_res1$basket$post.prob, mcmc_ref2$post.prob)
expect_equal(mcmc_res1$basket$ESS, mcmc_ref2$ESS, tolerance = 5e-2)
