library(testthat)
library(basket)
context("Fit exact models")

data(vemu_wide)

fb_reference <- readRDS("reference_data/mem-exact.rds")

baskets <- c(1, 2, 3)

vemu_wide1 <- vemu_wide[baskets, ]

# Full Bayes
exact_res <- mem_exact(
  responses = vemu_wide1$responders,
  size = vemu_wide1$evaluable,
  name = c("NSCLC", "CRC (vemu)", "CRC (vemu+cetu)"),
  p0 = 0.25, seed=1000
)

fb <- exact_res$basket

expect_equal(fb$PEP, fb_reference$basketwise$PEP)
expect_equal(fb$HPD, fb_reference$basketwise$HPD, tolerance = 5e-2)
expect_equal(fb$post.prob, fb_reference$basketwise$post.prob)
expect_equal(fb$ESS, fb_reference$basketwise$ESS)

res1 <- update(exact_res, p0=0.15)

basket_name(res1)
