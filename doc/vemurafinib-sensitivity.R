## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- warning = FALSE, message = FALSE-----------------------------------
# Load the basket library.
library(basket)
library(tibble)
library(dplyr)

# Set to sequential execution to get rid of warning.
if (require(foreach)) {
  foreach::registerDoSEQ()
}

# Take a look at the Vemurafinib data.
data(vemu_wide)
vemu_wide %>% glimpse

## ------------------------------------------------------------------------
fit <- mem_full_bayes(responses = vemu_wide$responders,
                      size = vemu_wide$evaluable,
                      name = vemu_wide$baskets,
                      p0 = 0.15)
fit

## ------------------------------------------------------------------------
plot_density(fit)

## ------------------------------------------------------------------------
plot_posterior_exchangeability(fit, size = 1.5, label = TRUE)

## ------------------------------------------------------------------------
plot_exchangeability(fit, size = 1.5)

