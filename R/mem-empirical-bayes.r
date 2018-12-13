
#' @title MEM Empirical Bayes
#' 
#' @description Fit the MEM model using empirical Bayesian inference.
#' @param responses the number of responses in each basket.
#' @param size the size of each basket.
#' @param name the name of each basket (default: NULL - no basket names).
#' @param p0 the null response rate for the poster probability calculation
#' (default 0.25).
#' @param shape1 the first shape parameter(s) for the prior of each basket
#' (default 0.5).
#' @param shape2 the second shape parameter(s) for the prior of each basket
#' (default 0.5).
#' @param alpha the trial significance.
#' @param upper_bound for constrained empirical Bayes, the upper bound on the
#' inclusion probability (default 1).
#' @param lower_bound for constrained empirical Bayes, the lower bound on the
#' inclusion probability (default 0).
#' @param call the call of the function (default NULL).
#' @examples
#' # 5 baskets, each with enrollement size 5
#' trial_sizes <- rep(5, 5) 
#' 
#' # The response rates for the baskets.
#' resp_rate <- 0.15
#' 
#' # The trials: a column of the number of responses and a column of the
#' # the size of each trial.
#' trials <- data.frame(responses=rbinom(trial_sizes, trial_sizes, resp_rate), 
#'                      size = trial_sizes)
#' 
#' mem_empirical_bayes(trials$responses, trials$size)
#' @importFrom stats quantile rbinom
#' @export
mem_empirical_bayes <- function(
  responses, 
  size, 
  name = NULL,
  p0 = 0.25, # Null response rate for Posterior Probability Calc.
  shape1 = 0.5, 
  shape2 = 0.5, 
  alpha = 0.05, 
  upper_bound = 1, 
  lower_bound = 0, 
  call = NULL) {

  pars <- list(UB = upper_bound, LB = lower_bound)

  if (length(responses) != length(size)) {
    stop(paste("The length of the responses and size parameters",
               "must be equal."))
  }
  if (length(shape1) == 1) {
    shape1 <- rep(shape1, length(responses))  
  } 
  if (length(shape2) == 1) {
    shape2 <- rep(shape2, length(responses))  
  } 
  
  Data <- list(X = responses, N = size, Basket = name)

  marg.M <- MEM_marginal(Data$X, Data$N, shape1, shape2)
  CDF <- MEM.cdf(Data, pars, p0, marg.M)
  ESS <- post.ESS(Data, pars, marg.M)
  HPD <- t(post.HPD(Data, pars, marg.M, alpha))

  if (is.null(call)) {
    call <- match.call()
  }

  browser()
  ret <- list(marg.M = marg.M, maximizer = marg.M$maximizer, CDF = CDF,
              ESS = ESS, HPD = HPD, responses = responses, size = size,
              name = name, p0 = p0, shape1 = shape2, call = call,
              U = MEM.w(Data, pars, marg.M))

  class(ret) <- c("empirical_bayes", "exchangeability_model")
  ret$samples <- sample_posterior(ret)
  ret$mean_est <- colMeans(ret$samples)
  ret$median_est <- apply(ret$samples, 2, median)
  ret
}

#' @importFrom foreach foreach %do%
#' @export
sample_posterior.empirical_bayes <- function(model, num_samples = 10000) {
  ret <- replicate(num_samples, samp.Post(model$responses, model$size,
                   model$U$models, model$U$weights[[1]]))
  K <- length(model$responses)
  ret <- rbind(ret, 
    foreach(j = 2:(K-1), .combine = rbind) %do% {
      Ii <- c(j,1:(j-1),(j+1):K)
      replicate(num_samples, 
                samp.Post(model$responses[Ii], model$size[Ii],
                          model$U$models, model$U$weights[[j]]))
    })
  j <- K
  Ii <- c(j,1:(j-1))
  ret <- rbind(ret,
               replicate(num_samples, 
                         samp.Post(model$responses[Ii], model$size[Ii], 
                                   model$U$models, model$U$weights[[j]]) ))
#  ret <- foreach(w = model$U$weights, .combine = cbind) %do% {
#    replicate(num_samples, samp.Post(model$responses, model$size, 
#      model$U$models, w))
#  }
  ret <- t(ret)
  dimnames(ret) <- list(NULL, model$name)
  ret
}
