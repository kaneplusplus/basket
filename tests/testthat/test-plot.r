
test_that("Plotting works.",
  {
    data(vemu_wide)

    set.seed(123)
    mh1 <- mem_mcmc(
      responses = vemu_wide$responders,
      size = vemu_wide$evaluable,
      name = vemu_wide$baskets,
      cluster_analysis = TRUE,
      p0 = c(0.15, 0.15, 0.15, 0.2, 0.15, 0.15),
      mcmc_iter = 1000,
      mcmc_burnin = 10
    )

    expect_snapshot(summary(mh1))
    pd <- plot_density(mh1)
    expect_true(inherits(pd, "ggplot"))
    pdb <- plot_density(mh1$basket)
    expect_true(inherits(pdb, "ggplot"))
#    pm <- plot_mem(mh1, type = "prior")
#    expect_true(inherits(pm, "gtable"))
    ppg <- plot_pep_graph(mh1)
    expect_true(inherits(ppg, "ggplot"))
  }
)
