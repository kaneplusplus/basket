
context("Plotting")

data(vemu_wide)

mh1 <- mem_mcmc(
  responses = vemu_wide$responders,
  size = vemu_wide$evaluable,
  name = vemu_wide$baskets,
  cluster_analysis = TRUE,
  p0 = c(0.15, 0.15, 0.15, 0.2, 0.15, 0.15),
  mcmc_iter = 500,
  mcmc_burnin = 100
)

expect_true(inherits(plot_density(mh1), "ggplot"))
expect_true(inherits(plot_density(mh1$basket), "ggplot"))
expect_true(inherits(plot_mem(mh1, type = "prior"), "gtable"))
expect_true(inherits(plot_pep_graph(mh1), "ggplot"))
