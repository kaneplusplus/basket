context("fit full Bayes exact models")

# Load the reference outputs.
#load("../../inst/test-data/eb_reference.rda")
#load("../../inst/test-data/eb_reference_ub.rda")
load("../../inst/test-data/fb_reference.rda")

data(vemu_wide)

baskets <- 1:3

vemu_wide1 <- vemu_wide[baskets, ]


# Full Bayes
exact_res <- mem_full_bayes_exact(
  responses = vemu_wide1$responders,
  size = vemu_wide1$evaluable,
  name = c("8", "0", "1"), p0 = rep(0.25, length(basket))
)

fb <- exact_res$basketwise

expect_equal(fb$PEP, fb_reference$PEP)
expect_equal(fb$HPD, fb_reference$HPD, tolerance = 5e-2)
expect_equal(fb$post.prob, fb_reference$CDF)
expect_equal(fb$ESS, fb_reference$ESS)

