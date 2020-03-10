
context("Fit MCMC models")

data(vemu_wide)

baskets <- c(1, 5, 6)

vemu_wide1 <- vemu_wide[baskets, ]

mcmc_ref1 <- readRDS("reference_data/mem-mcmc2.rds")
mcmc_ref2 <- mcmc_ref1$basket

time_taken <- system.time({
  mcmc_res1 <- mem_mcmc(
    responses = vemu_wide1$responders,
    size = vemu_wide1$evaluable,
    name = vemu_wide1$baskets,
    cluster_analysis = TRUE,
    p0 = 0.15,
    mcmc_iter = 100
  )
})

test_that("Exact corner case models", {

  skip_on_cran()

  mcmc_res1_basket <- 
    basket(responses = vemu_wide1$responders,
      size = vemu_wide1$evaluable,
      name = vemu_wide1$baskets,
      cluster_analysis = TRUE,
      p0 = 0.15,
      mcmc_iter = 100
    )

  mcmc_res2 <- mem_mcmc(
    responses = c(4, 3, 0),
    size = c(10, 3, 0),
    name = letters[1:3],
    cluster_analysis = TRUE,
    p0 = 0.25)

  mcmc_res3 <- mem_mcmc(
    responses = c(4, 3),
    size = c(10, 3),
    name = letters[1:2],
    cluster_analysis = TRUE,
    p0 = 0.25)

  mcmc_res3 <- mem_mcmc(
    responses = c(4, 3),
    size = c(10, 3),
    name = letters[1:2],
    cluster_analysis = FALSE,
    p0 = 0.25)

  expect_equal(mcmc_res1$basket[-10], mcmc_res1_basket$basket[-10], 
               tolerance = 0.5)

})

expect_true(is.matrix(cluster_map(mcmc_res1)))

expect_equal(class(summary(mcmc_res1$basket)), "mem_basket_summary")
expect_equal(class(summary(mcmc_res1$cluster)), "mem_cluster_summary")
expect_equal(class(summary(mcmc_res1)), "mem_summary")

expect_true(inherits(plot_density(mcmc_res1$basket), "ggplot"))
expect_true(inherits(
  plot_mem(mcmc_res1$basket, type = c("pep", "map")),
  "gtable"
))
expect_true(inherits(plot_pep(mcmc_res1$basket), "ggplot"))
expect_true(inherits(plot_map(mcmc_res1$basket), "ggplot"))

expect_equivalent(mcmc_res1$basket$pep, mcmc_ref2$PEP, tolerance = 5e-2)
expect_equivalent(mcmc_res1$basket$hpd, mcmc_ref2$HPD, tolerance = 5e-2)
expect_equivalent(mcmc_res1$basket$post_prob, mcmc_ref2$post.prob,
  tolerance = 5e-2
)
expect_equivalent(mcmc_res1$basket$ess, mcmc_ref2$ESS, tolerance = 5e-2)

mcmc_lower <- mem_mcmc(
  responses = vemu_wide1$responders,
  size = vemu_wide1$evaluable,
  name = vemu_wide1$baskets,
  alternative = "less",
  p0 = 0.15, mcmc_iter = 10000
)

expect_equal(class(summary(mcmc_lower)), "mem_summary")

mcmc_res2 <- update_p0(mcmc_res1, p0 = 0.18)

expect_true(is.list(cluster_baskets(mcmc_res2)))
expect_equivalent(cluster_map(mcmc_res2), matrix(1, nrow = 3, ncol = 3))
expect_true(is.matrix(cluster_pep(mcmc_res2)))
expect_equivalent(basket_map(mcmc_res2), matrix(1, nrow = 3, ncol = 3))
expect_true(is.matrix(basket_pep(mcmc_res2)))

# Single basket test
time_taken <- system.time({
  mcmc_res3 <- mem_mcmc(
    responses = c(0),
    size = c(1),
    name = "Basket 1",
    p0 = 0.15,
    mcmc_iter = 100,
    mcmc_burnin = 100
  )
})
