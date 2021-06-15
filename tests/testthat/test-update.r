
test_that("Update works.",
  {

    data(vemu_wide)

    set.seed(123)
    mh1 <- mem_mcmc(
      responses = vemu_wide$responders,
      size = vemu_wide$evaluable,
      name = vemu_wide$baskets,
      p0 = c(0.15, 0.15, 0.15, 0.2, 0.15, 0.15),
      mcmc_iter = 10,
      mcmc_burnin = 10
    )


    set.seed(123)
    expect_snapshot(summary(update_p0(mh1, 0.25)))
  }
)
