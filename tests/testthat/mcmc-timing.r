
data(vemu_wide)

library(profvis)

profvis({
mcmc_res1 <- mem_mcmc(
  responses = vemu_wide$responders,
  size = vemu_wide$evaluable,
  name = vemu_wide$baskets,
  p0 = 0.15
)
})
