context("fit MEM models")

# Load the reference outputs.
load("../../inst/test-data/eb_reference.rda")
load("../../inst/test-data/eb_reference_ub.rda")
load("../../inst/test-data/fb_reference.rda")

data(vemu_wide)

baskets <- 1:3

vemu_wide <- vemu_wide[baskets, ]

# Full Bayes
#fb <- mem_full_bayes(
#  responses = vemu_wide$responders,
#  size = vemu_wide$evaluable,
#  name = vemu_wide$baskets
#)

# Visualization
# plot_posterior_exchangeability(fb,
#   basket_name_hoffset = c(0, -0.1, -0.4, 0, -0.1, 0))


#expect_equal(fb$maximizer, fb_reference$maximizer)
#expect_equal(fb$PEP, fb_reference$PEP)
#expect_equal(fb$HPD, fb_reference$HPD, tolerance = 5e-2)
#expect_equal(fb$CDF, fb_reference$CDF)
#expect_equal(fb$ESS, fb_reference$ESS)

# Empirical Bayes
eb <- mem_empirical_bayes(
  responses = vemu_wide$responders,
  size = vemu_wide$evaluable,
  name = vemu_wide$baskets
)

expect_equal(eb$maximizer, eb_reference$maximizer)
expect_equal(eb$PEP, eb_reference$PEP)
expect_equal(eb$HPD, eb_reference$hpd, tolerance = 5e-2)
expect_equal(eb$CDF, eb_reference$cdf)
expect_equal(eb$ESS, eb_reference$ess)

# Constrained Empirical Bayes
ebc <- mem_empirical_bayes(
  responses = vemu_wide$responders,
  size = vemu_wide$evaluable,
  name = vemu_wide$baskets,
  upper_bound = 0.2
)

expect_equal(ebc$maximizer, eb_reference_ub$maximizer)
expect_equal(ebc$PEP, eb_reference_ub$PEP)
expect_equal(ebc$HPD, eb_reference_ub$hpd, tolerance = 5e-2)
expect_equal(ebc$CDF, eb_reference_ub$cdf)
expect_equal(ebc$ESS, eb_reference_ub$ess)
