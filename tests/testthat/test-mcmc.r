library(testthat)
library(basket)

context("Fit MCMC models")

data(vemu_wide)

baskets <- c(1, 2, 3)
#baskets <- 1:6

vemu_wide1 <- vemu_wide[baskets, ]

mcmc_ref1 <- readRDS("reference_data/mem-mcmc1.rds")

mcmc_ref2 <- mcmc_ref1$basketwise
start.time <- Sys.time()
mcmc_res1 <- mem_mcmc(
  responses = vemu_wide1$responders,
  size = vemu_wide1$evaluable,
  name = vemu_wide1$baskets,
  p0 = 0.15
)
end.time <- Sys.time()
time.taken <- end.time - start.time
cat ("Time used: ", time.taken, "\n")
t<-sample_posterior(mcmc_res1$basket)

summary(mcmc_res1$basket)
summary(mcmc_res1$cluster)
summary(mcmc_res1)
#print(mcmc_res1)
plot_density(mcmc_res1$basket)
plot_mem(mcmc_res1$basket,  type=c("pep", "map"))
plot_pep(mcmc_res1$basket)
plot_map(mcmc_res1$basket)
expect_equal(mcmc_res1$basket$PEP, mcmc_ref2$PEP)
expect_equal(mcmc_res1$basket$HPD, mcmc_ref2$HPD, tolerance = 5e-2)
expect_equal(mcmc_res1$basket$post.prob, mcmc_ref2$post.prob)
expect_equal(mcmc_res1$basket$ESS, mcmc_ref2$ESS, tolerance = 5e-2)


mcmc_lower <- mem_mcmc(
  responses = vemu_wide1$responders,
  size = vemu_wide1$evaluable,
  name = vemu_wide1$baskets,
  alternative = "less",
  p0 = 0.15
)
summary(mcmc_lower)
mcmc_res2 <- update_p0(mcmc_res1, p0=0.18)
#print(mcmc_res2$basket)
cluster_baskets(mcmc_res2)
cluster_map(mcmc_res2)
cluster_pep(mcmc_res2)
basket_map(mcmc_res2)
basket_pep(mcmc_res2)
