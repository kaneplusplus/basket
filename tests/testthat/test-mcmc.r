
data(vemu_wide)

baskets <- c(1, 5, 6)

vemu_wide1 <- vemu_wide[baskets, ]


test_that("MCMC can be calculated correctly.",
  {
    set.seed(123)
    mcmc_res1 <- mem_mcmc(
      responses = vemu_wide1$responders,
      size = vemu_wide1$evaluable,
      name = vemu_wide1$baskets,
      cluster_analysis = TRUE,
      p0 = 0.15,
      mcmc_burnin = 100,
      mcmc_iter = 100
    )
    expect_snapshot(unclass(summary(mcmc_res1)))

    expect_true(is.matrix(cluster_map(mcmc_res1)))
    expect_equal(class(summary(mcmc_res1$basket)), "mem_basket_summary")
    expect_equal(class(summary(mcmc_res1$cluster)), "mem_cluster_summary")
    expect_equal(class(summary(mcmc_res1)), "mem_summary")

    pd <- plot_density(mcmc_res1$basket)
    expect_true(inherits(pd, "ggplot"))
#    pm <- plot_mem(mcmc_res1$basket, type = c("pep", "map"))
#    expect_true(inherits(pm, "gtable"))
    expect_true(inherits(plot_pep(mcmc_res1$basket), "ggplot"))
    expect_true(inherits(plot_map(mcmc_res1$basket), "ggplot"))
  }
)


test_that("Exact corner case models", {

  skip_on_cran()

  set.seed(123)
  mcmc_res1 <- mem_mcmc(
    responses = vemu_wide1$responders,
    size = vemu_wide1$evaluable,
    name = vemu_wide1$baskets,
    cluster_analysis = TRUE,
    p0 = 0.15,
    mcmc_burnin = 100,
    mcmc_iter = 100
  )

  set.seed(123)
  mcmc_res1_basket <- 
    basket(responses = vemu_wide1$responders,
      size = vemu_wide1$evaluable,
      name = vemu_wide1$baskets,
      cluster_analysis = TRUE,
      p0 = 0.15,
      mcmc_iter = 100,
      mcmc_burnin = 100
    )

  set.seed(123)
  expect_snapshot(
    summary(mem_mcmc(
      responses = c(4, 3, 0),
      size = c(10, 3, 0),
      name = letters[1:3],
      cluster_analysis = TRUE,
      mcmc_iter = 100,
      mcmc_burnin = 100,
      p0 = 0.25
    ))
  )

  set.seed(123)
  expect_snapshot(
    summary(mem_mcmc(
      responses = c(4, 3),
      size = c(10, 3),
      name = letters[1:2],
      cluster_analysis = TRUE,
      mcmc_iter = 100,
      mcmc_burnin = 100,
      p0 = 0.25)
    )
  )

  set.seed(123)
  expect_snapshot(
    summary(mem_mcmc(
      responses = c(4, 3),
      size = c(10, 3),
      name = letters[1:2],
      cluster_analysis = FALSE,
      mcmc_iter = 100,
      mcmc_burnin = 100,
      p0 = 0.25))
  )

  expect_equal(mcmc_res1$basket[-10], mcmc_res1_basket$basket[-10], 
               tolerance = 0.5)

  set.seed(123)
  mcmc_lower <- mem_mcmc(
    responses = vemu_wide1$responders,
    size = vemu_wide1$evaluable,
    name = vemu_wide1$baskets,
    alternative = "less",
    p0 = 0.15, 
    mcmc_iter = 100,
    mcmc_burnin = 100
  )

  expect_equal(class(summary(mcmc_lower)), "mem_summary")

  mcmc_res2 <- update_p0(mcmc_res1, p0 = 0.18)

  expect_true(is.list(cluster_baskets(mcmc_res2)))
  expect_equal(cluster_map(mcmc_res2), matrix(1, nrow = 3, ncol = 3),
               ignore_attr = TRUE)
  expect_true(is.matrix(cluster_pep(mcmc_res2)))
  expect_equal(basket_map(mcmc_res2), matrix(1, nrow = 3, ncol = 3),
               ignore_attr = TRUE)
  expect_true(is.matrix(basket_pep(mcmc_res2)))

})

