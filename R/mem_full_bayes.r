
#' @title MEM Full Bayes
#' 
#' @description Fit the MEM model using full Bayesian inference.
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
#' mem_full_bayes(trials$responses, trials$size)
#' @importFrom foreach foreach %dopar%
#' @export
mem_full_bayes <- function(
  responses, 
  basket_size, 
  p0 = 0.25, # Null response rate for Posterior Probability Calc.
  shape1 = 0.5, 
  shape2 = 0.5, 
  prior_inclusion = diag(length(responses))/2 + 
                      matrix(0.5, nrow = length(responses), 
                             ncol = length(responses)),
  alpha = 0.05) {

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
  pr.Inclus <- prior_inclusion
  alp <- alpha

  xvec <- responses
  nvec <- basket_size
  mod.mat <- list()
  k <- length(xvec)-1
  j <- 1

  # TODO: turn this into a foreach loop
  while (k > 0) { 
    m <- as.matrix(expand.grid( rep(list(c(0,1)), k) ))

    mod.mat[[j]] <- m[order(rowSums(m)),]
#( as.matrix(expand.grid( rep(list(c(0,1)), k) )) )[order(rowSums(as.matrix(expand.grid( rep(list(c(0,1)), k) )))),]
    k <- k - 1
    j <- j + 1 
  }

  H <- length(mod.mat)
  temp <- list()
  
  for (h in 1:(H-1)) {
    temp[[h]] <- 1:dim(mod.mat[[h]])[1] 
  }
  temp[[H]] <- 1:2
  Mod.I <- as.matrix(expand.grid( temp ))

  ## identify maximizer of marginal density ##
  log.Marg <- apply(Mod.I, MARGIN = 1, FUN = logMarg.Dens, mod.mat, xvec, 
                    nvec, avec, bvec)

  max.i <- order(log.Marg, decreasing=TRUE)[1]

  MAX <- MEM.mat(Mod.I[max.i,], mod.mat, length(xvec))
  colnames(MAX) <- rownames(MAX) <- xvec

  ## compute prior + posterior MEM probabilities ##
  PRIOR <- apply(Mod.I, MARGIN = 1, FUN = mem.Prior, mod.mat, pr.Inclus)
  POST <- ( exp(log.Marg) * PRIOR ) / sum( exp(log.Marg) * PRIOR )
  map.i <- order(POST, decreasing = TRUE)[1]

  MAP <- MEM.mat(Mod.I[map.i,], mod.mat, length(xvec))
  colnames(MAP) <- rownames(MAP) <- xvec

  ## Posterior Exchangeability Prob ##
  PEP. <- matrix(NA, length(xvec), length(xvec))
  colnames(PEP.) <- rownames(PEP.) <- xvec; diag(PEP.) <- 1
  i <- NA
  pep <- foreach(i = 1:(nrow(PEP.)-1), combine = list, 
                 .multicombine = TRUE) %dopar% {
    compPEP(i, length(xvec), Mod.I, mod.mat, pr.Inclus, log.Marg, PRIOR)
  }
  l <- lower.tri(PEP.)
  PEP.[l] <- unlist(pep)
  PEP <- t(PEP.)

  ## Marginal CDF and ESS ##
  pweights <- foreach(j = 1:length(xvec), .combine = list, 
                      .multicombine = TRUE) %dopar% {
    post.Weights(j, length(xvec), Mod.I, mod.mat, pr.Inclus, log.Marg, PRIOR)
  }

  pESS <- CDF <- rep(NA, length(xvec))
  names(pESS) <- names(CDF) <- xvec
  HPD <- matrix(NA,2,length(xvec))
  rownames(HPD)=c("lower","upper")
  colnames(HPD) <- xvec
  models <- cbind(rep(1, dim(mod.mat[[1]])[1]), mod.mat[[1]])

  CDF[1] <- eval.Post(p0, xvec, nvec, models, pweights[[1]])
  pESS[1] <- pweights[[1]] %*% ESS(xvec, nvec, models)
  HPD[,1] <- boa.hpd(
    replicate(10000, samp.Post(xvec, nvec, models, pweights[[1]]) ), 
    alp)

  K <- length(xvec)
  for(j in 2:(K-1)) { 
    Ii <- c(j,1:(j-1),(j+1):K)
    CDF[j] <- eval.Post(p0, xvec[Ii], nvec[Ii], models, pweights[[j]])
    pESS[j] <- pweights[[j]]%*%ESS(xvec[Ii], nvec[Ii], models)
    HPD[,j] <- boa.hpd(
      replicate(10000, samp.Post(xvec[Ii], nvec[Ii], models, pweights[[j]]) ), 
      alp)
  }
  j <- j + 1
  Ii <- c(j,1:(j-1))
  CDF[j] <- eval.Post(p0, xvec[Ii], nvec[Ii], models, pweights[[j]])
  pESS[j] <- pweights[[j]] %*% ESS(xvec[Ii], nvec[Ii], models)
  HPD[,j] <- boa.hpd( 
    replicate(10000, samp.Post(xvec[Ii], nvec[Ii], models, pweights[[j]]) ), 
    alp)

  # TODO: check with Brian to make sure this is correct.
  median_est <- vapply(pweights, 
    function(w) {
      quantile(replicate(10000, samp.Post(xvec, nvec, models, w)), 0.5)
    }, as.numeric(NA))

  mean_est <- vapply(pweights, 
    function(w) {
      mean(replicate(10000, samp.Post(xvec, nvec, models, w)))
    }, as.numeric(NA))

  ret <- list(mod.mat = mod.mat, maximizer = MAX, MAP = MAP, PEP = PEP, 
             CDF = CDF, ESS = pESS, HPD = HPD, mean_est = mean_est,
             median_est = median_est, responses = responses,
             basket_size = basket_size, p0 = p0, shape1 = shape1,
             shape2 = shape2)

  class(ret) <- c("full_bayes", "mem")
  ret 
}
