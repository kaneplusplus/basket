library(testthat)
library(basket)
context("fit full Bayes exact models")


data(vemu_wide)
data(fb_reference)

baskets <- 1:3

vemu_wide1 <- vemu_wide[baskets, ]


# Full Bayes
exact_res <- mem_full_bayes_exact(
  responses = vemu_wide1$responders,
  size = vemu_wide1$evaluable,
  name = c("8", "0", "1"), p0 = rep(0.25, length(baskets))
)


fb <- exact_res$basketwise

expect_equal(fb$PEP, fb_reference$PEP)
expect_equal(fb$HPD, fb_reference$HPD, tolerance = 5e-2)
expect_equal(fb$post.prob, fb_reference$CDF)
expect_equal(fb$ESS, fb_reference$ESS)

