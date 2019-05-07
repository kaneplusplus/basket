library(testthat)
library(basket)

context("Plotting")

data(vemu_wide)

mh1 <- mem_mcmc(responses = vemu_wide$responders, 
                size = vemu_wide$evaluable,
                name= vemu_wide$baskets,
                p0 = c(0.15, 0.15, 0.15, 0.2, 0.15, 0.15))

expect_true(inherits(plot_density(mh1), "ggplot"))
expect_true(inherits(plot_density(mh1$basket), "ggplot"))
expect_true(inherits(plot_density(mh1$cluster), "ggplot"))
expect_true(inherits(plot_density(mh1$cluster), "ggplot"))
expect_true(inherits(plot_mem(mh1, type = "prior"), "gtable"))

