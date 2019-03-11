library(testthat)
library(basket)
context("Fit exact models")

data(vemu_wide)

fb_reference <- readRDS("reference_data/mem-exact.rds")

baskets <- c(1, 5, 6)

vemu_wide1 <- vemu_wide[baskets, ]

# Full Bayes
exact_res <- mem_exact(
  responses = vemu_wide1$responders,
  size = vemu_wide1$evaluable,
  name = vemu_wide1$baskets,
  p0 = 0.25
)

fb <- exact_res$basketwise

expect_equal(exact_res$PEP, fb_reference$PEP)
expect_equal(exact_res$HPD, fb_reference$HPD, tolerance = 5e-2)
expect_equal(exact_res$post.prob, fb_reference$CDF)
expect_equal(exact_res$ESS, fb_reference$ESS)

