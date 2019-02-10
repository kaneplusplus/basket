library(devtools)

document()
data(vemu_wide)

source("../inst/code-from-brian/Utils/R-functions-v11.R")
source("../inst/code-from-brian/MEM-FullBayes-function.R")

prob <- rep(0.15, nrow(vemu_wide))

p0 <- 0.25

# Only test the first 3 baskets to save time.
baskets <- 1:3
vemu_wide <- vemu_wide[baskets, ]

# EB unconstrained.
marg_m <- MEM_marginal(
  vemu_wide$responders, vemu_wide$evaluable,
  rep(0.5, nrow(vemu_wide)),
  rep(0.5, nrow(vemu_wide))
)

Data <- list(X = vemu_wide$responders, N = vemu_wide$evaluable)
pars <- list(UB = 1, LB = 0)
eb_reference <- list(
  marg_m = marg_m, maximizer = marg_m$maximizer,
  cdf = MEM.cdf(Data, pars, p0, marg_m),
  ess = post.ESS(Data, pars, marg_m),
  hpd = t(post.HPD(Data, pars, marg_m, 0.05))
)

save(eb_reference, file = "../data/eb_reference.rda")

# EB constrained with upper bound of 0.2.
pars <- list(UB = 0.2, LB = 0)
eb_reference_ub <- list(
  marg_m = marg_m, maximizer = marg_m$maximizer,
  cdf = MEM.cdf(Data, pars, p0, marg_m),
  ess = post.ESS(Data, pars, marg_m),
  hpd = t(post.HPD(Data, pars, marg_m, 0.05))
)

save(eb_reference_ub, file = "../data/eb_reference_ub.rda")

# Full Bayes
inclusion_prob <- matrix(0.5, length(Data$X), length(Data$X))
diag(inclusion_prob) <- 1
fb_reference <- MEM_FullBayes(
  Data$X, Data$N, rep(0.5, length(Data$X)),
  rep(0.5, length(Data$X)), inclusion_prob, 0.05
)

save(fb_reference, file = "../data/fb_reference.rda")
