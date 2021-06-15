
test_that(
  "cluster_basket() error when given NULL",
  {
    expect_error(cluster_baskets(NULL))
  }
)

test_that(
  "baket_pep() error when given NULL",
  {
    expect_error(basket_pep(NULL))
  }
)

test_that(
  "baket_map() error when given NULL",
  {
    expect_error(basket_map(NULL))
  }
)

test_that(
  "cluster_pep() error when given NULL",
  {
    expect_error(cluster_pep(NULL))
  }
)

test_that(
  "cluster_map() error when given NULL",
  {
    expect_error(cluster_map(NULL))
  }
)

test_that(
  "sample_posterior() error when given NULL",
  {
    expect_error(sample_posterior(NULL))
  }
)

test_that(
  "basket_name() error when given NULL",
  {
    expect_error(basket_name(NULL))
  }
)

test_that(
  "update_p0() error when given NULL",
  {
    expect_error(update_p0(NULL))
  }
)

test_that(
  "plot_density() error when given NULL",
  {
    expect_error(plot_density(NULL))
  }
)

test_that(
  "plot_pep() error when given NULL",
  {
    expect_error(plot_pep(NULL))
  }
)

test_that(
  "exchangeogram() error when given bunk",
  {
    expect_error(exchangeogram(matrix(c(1, 2, 3, 4), ncol = 2)))
  }
)

test_that(
  "plot_mem() error when given NULL",
  {
    expect_error(plot_mem(NULL))
  }
)

test_that(
  "plot_map() error when given NULL",
  {
    expect_error(plot_map(NULL))
  }
)

test_that(
  "mem_exact() error when given bunk",
  {
    expect_error(mem_exact(c(1, 2, 3), 1))
  }
)

test_that(
  "mem_mcmc() error when given bunk",
  {
    expect_error(mem_mcmc(c(1, 2, 3), 1))
  }
)
