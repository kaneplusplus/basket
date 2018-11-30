
#' @title MEM Empirical Bayes
#' 
#' @description Fit the MEM model using empirical Bayesian inference.
#' @param responses the number of responses in each basket.
#' @param basket_size the size of each basket.
#' @param p0 the null response rate for the poster probability calculation.
#' @param shape1 the first shape parameter(s) for the prior of each basket
#' (default 0.5).
#' @param shape2 the second shape parameter(s) for the prior of each basket
#' (default 0.5).
#' @param prior_inclusion the matrix giving the prior inclusion probability
#' for each pair of baskets. The default is on on the main diagonal and 0.5
#' elsewhere.
#' @param alpha the trial significance.
#' @param upper_bound for constrained empirical Bayes, the upper bound on the
#' inclusion probability (default 1).
#' @param lower_bound for constrained empirical Bayes, the lower bound on the
#' inclusion probability (default 0).
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
  basket_size, 
  p0 = 0.25, # Null response rate for Posterior Probability Calc.
  shape1 = 0.5, 
  shape2 = 0.5, 
  prior_inclusion = diag(length(responses))/2 + 
                      matrix(0.5, nrow = length(responses), 
                             ncol = length(responses)),
  alpha = 0.05, 
  upper_bound = 1, 
  lower_bound = 0) {

  pars <- list(UB = upper_bound, LB = lower_bound)

  if (length(responses) != length(basket_size)) {
    stop(paste("The length of the responses and basket_size parameters",
               "must be equal."))
  }
  if (length(shape1) == 1) {
    avec <- rep(shape1, length(responses))  
  } else {
    avec <- shape1
  }
  if (length(shape2) == 1) {
    bvec <- rep(shape2, length(responses))  
  } else {
    bvec <- shape2
  }
  INITIAL <- prior_inclusion

  xvec <- responses
  nvec <- basket_size

  mod.mat <- list()
  k <- length(xvec)-1
  j <- 1

  while(k > 0) {
    mod.mat[[j]] <- as.matrix(expand.grid( rep(list(c(0,1)), k) ))[
        order(rowSums(as.matrix(expand.grid( rep(list(c(0,1)), k) )))),]

    k <- k - 1
    j <- j + 1
  }

  if( is.na(INITIAL[1]) ) {
    M.init <- diag(1,length(xvec))
  } else {
    M.init <- INITIAL
  }

  MOld <- M.init
  d <- length(xvec)

  for(i in 1:(d-1)) {
    for(k in (i+1):d) {
      if (MOld[i,k]!=MOld[k,i]) {
        MOld[i,k]<-MOld[k,i]<-rbinom(1,1,0.5)
      }
    }
  }

  logDenOld <- logMarg.DensSA(MOld, mod.mat, xvec, nvec, avec, bvec)
  for(iter in 1:10000) {
    changed <- FALSE
    for( i in 1:(d-1)) {
      for(k in (i+1):d) {
        loc <- c(k,i)
        MNew <- MOld
        MNew[loc[1],loc[2]] <- MNew[loc[2],loc[1]] <- 1 - MOld[loc[1],loc[2]]
        logDenNew <- logMarg.DensSA(MNew, mod.mat, xvec, nvec, avec, bvec)
        if (logDenNew > logDenOld) {
          MOld<-MNew
          changed<-TRUE
          logDenOld<-logDenNew
        }
      }
    }
    if (!changed) {
      break
    }
  }

  colnames(M.init) <- rownames(M.init) <- xvec
  colnames(MOld) <- rownames(MOld) <- xvec
  ret <- list(mod.mat = mod.mat, init=M.init, maximizer=MOld)
  ret$CDF <- MEM.cdf(list(X=xvec, N=nvec), pars, p0, ret)
  ret$ESS <- post.ESS(list(X=xvec, N=nvec), pars, ret)
  ret$HPD <- t(post.HPD(list(X=xvec, N=nvec), pars, ret, alpha))
  U <- MEM.w(list(X=xvec, N=nvec), pars, ret)
  # TODO: Again, make sure with Brian that this is correct.
  ret$mean_est <- vapply(U$weights, 
    function(w) {
      mean(replicate(10000, samp.Post(xvec, nvec, U$models, w)))
    }, as.numeric(NA))
  ret$median_est <- vapply(U$weights,
    function(w) {
      quantile(replicate(10000, samp.Post(xvec, nvec, U$models, w)), 
               0.5)
    }, as.numeric(NA))

  class(ret) <- c("emprical_bayes", "mem")
  ret
}

