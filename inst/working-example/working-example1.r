library(basket)
library(tibble)
library(dplyr)

# Load the vemurafinib data for analysis.
data(vemu_wide)
vemu_wide

# Fit the MEM model using Empirical Bayes.
eb <- mem_empirical_bayes(responses = vemu_wide$responders, 
                          size = vemu_wide$evaluable,
                          name = vemu_wide$baskets,
                          p0 = 0.15)
eb

# Plot the posteriors
plot_density(eb)

plot_exchangeability(eb, size = 1.5, label = TRUE)

# Fit the MEM model using full Bayes.
fb <- mem_full_bayes(responses = vemu_wide$responders, 
                     size = vemu_wide$evaluable,
                     name = vemu_wide$baskets
                     p0 = 0.15)


# Plot the posteriors
plot_density(fb)

plot_posterior_exchangeability(fb, size = 1.5, label = TRUE)

plot_exchangeability(fb, size = 1.5, label = TRUE)




