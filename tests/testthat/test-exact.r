
data(vemu_wide)

baskets <- c(1, 2, 3)

vemu_wide1 <- vemu_wide[baskets, ]


# Full Bayes
test_that("Full Bayes can be calculated correctly.",
  {

    set.seed(123)
    exact_res <- mem_exact(
      responses = vemu_wide1$responders,
      size = vemu_wide1$evaluable,
      name = vemu_wide1$baskets,
      cluster_analysis = TRUE,
      p0 = 0.25)

    expect_snapshot(unclass(summary(exact_res)))

    expect_true(class(summary(exact_res)) == "mem_summary")
    expect_output(print(summary(exact_res)), "The Null Response Rates")

    expect_true(class(summary(exact_res$basket)) == "mem_basket_summary")
    expect_output(print(summary(exact_res$basket)), "The Null Response Rates")
    expect_true(class(summary(exact_res$cluster)) == "mem_cluster_summary")

    expect_output(
      print(summary(exact_res$cluster)), "The Null Response Rates")

    res1 <- update(exact_res, p0 = 0.15)
    expect_true(inherits(res1, "mem_exact"))

    expect_equal(basket_name(res1), c("NSCLC", "CRC (vemu)", "CRC (vemu+cetu)"))

    expect_true(inherits(
      exact_res_up <- update_p0(exact_res, alternative = "less"),
      "exchangeability_model"
    ))
  }
)

test_that("Exact corner case models", {

  skip_on_cran()

  set.seed(123)
  exact_res <- mem_exact(
    responses = vemu_wide1$responders,
    size = vemu_wide1$evaluable,
    name = vemu_wide1$baskets,
    cluster_analysis = TRUE,
    p0 = 0.25)

  set.seed(123)
  exact_res_basket <- basket(
    responses = vemu_wide1$responders,
    size = vemu_wide1$evaluable,
    name = vemu_wide1$baskets,
    cluster_analysis = TRUE,
    p0 = 0.25,
    method = "exact")

  set.seed(123)
  expect_snapshot(
    basket(
      responses = vemu_wide1$responders,
      size = vemu_wide1$evaluable,
      name = vemu_wide1$baskets,
      cluster_analysis = TRUE,
      p0 = 0.25,
      method = "exact"
    )
  )

  
  set.seed(123)
  expect_snapshot(
    mem_exact(
      responses = c(4, 3, 0),
      size = c(10, 3, 0),
      name = letters[1:3],
      cluster_analysis = TRUE,
      p0 = 0.25
    )
  )

  set.seed(123)
  expect_snapshot(
    mem_exact(
      responses = c(4, 3),
      size = c(10, 3),
      name = letters[1:2],
      cluster_analysis = TRUE,
      p0 = 0.25
    )
  )

  # Remove the call and check for equality.
  expect_equal(exact_res$basket[-18], 
               exact_res_basket$basket[-18], 
               tolerance = 0.5)

})

