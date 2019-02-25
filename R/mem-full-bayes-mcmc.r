


#' @title MEM Full Bayes MCMC method
#'
#' @description Fit the MEM model using full Bayesian Metropolis-Hasting MCMC inference.
#' @param responses the number of responses in each basket.
#' @param size the size of each basket.
#' @param name the name of each basket.
#' @param p0 the null response rate for the poster probability calculation
#' (default 0.15).
#' @param shape1 the first shape parameter(s) for the prior of each basket
#' (default 0.5).
#' @param shape2 the second shape parameter(s) for the prior of each basket
#' (default 0.5).
#' @param Prior the matrix giving the prior inclusion probability
#' for each pair of baskets. The default is on on the main diagonal and 0.5
#' elsewhere.
#' @param HPD.alpha the highest posterior density trial significance.
#' @param alternative the alternative case defination (default greater)
#' @param niter.MCMC the number of MCMC iterations.
#' @param Initial the initial MEM matrix.
#' @param seed the random number seed.
#' @param call the call of the function.
#' @examples
#' # 5 baskets, each with enrollement size 5
#' trial_sizes <- rep(5, 5)
#' 
#' # The response rates for the baskets.
#' resp_rate <- 0.15
#' 
#' # The trials: a column of the number of responses and a column of the
#' # the size of each trial.
#' trials <- data.frame(
#'   responses = rbinom(trial_sizes, trial_sizes, resp_rate),
#'   size = trial_sizes
#' )
#' mem_full_bayes_mcmc(trials$responses, trials$size)
#' @importFrom stats median
#' @importFrom igraph graph_from_adjacency_matrix cluster_louvain E
#' @export
mem_full_bayes_mcmc <- function(responses, size, name, p0 = 0.15, shape1 = 0.5,
  shape2 = 0.5, Prior = diag(length(responses)) / 2 + 
                        matrix(0.5, nrow = length(responses), 
  ncol = length(responses)), HPD.alpha = 0.05, alternative = "greater",
  niter.MCMC = 10000, Initial = NA, seed = 1000, call = NULL) {

  set.seed(seed)
  if (is.null(getDoParName())) {
    registerDoSEQ()
  }
  if (length(responses) != length(size)) {
    stop(paste(
      "The length of the responses and size parameters",
      "must be equal."
    ))
  }
  if (length(shape1) == 1) {
    shape1 <- rep(shape1, length(responses))
  }
  if (length(shape2) == 1) {
    shape2 <- rep(shape2, length(responses))
  }

  isP0Vector <- TRUE
  if (length(p0) == 1) {
    isP0Vector <- FALSE
    p0 <- rep(p0, length(responses))
  }
  ### Produce sample space of MEM ###
  mod.mat <- list()
  k <- length(responses) - 1
  j <- 1
  while (k > 0) {
    mod.mat[[j]] <-
      (as.matrix(expand.grid(rep(list(
        c(0, 1)
      ), k))))[order(rowSums(as.matrix(expand.grid(rep(
        list(c(0, 1)), k
      ))))), ]
    k <- k - 1
    j <- j + 1
  }

  ### Set initial MEM for MCMC ###
  if (is.na(Initial)[1]) {
    M.init <- diag(1, length(responses))
  } else {
    M.init <- Initial
  }

  MOld <- M.init

  ### Check symmetry of Initial MEM ###
  d <- length(responses)
  update.Init <- FALSE
  for (i in 1:(d - 1))
    for (k in (i + 1):d)
    {
      if (MOld[i, k] != MOld[k, i]) {
        MOld[i, k] <- MOld[k, i] <- rbinom(1, 1, 0.5)
        update.Init <- TRUE
      }
    }
  if (update.Init) {
    M.init <- MOld
  }

  ### Create Map for Proposal Distribution ###
  M <- diag(NA, nrow(MOld))
  K <- 1
  for (ii in 1:(nrow(MOld) - 1)) {
    for (jj in (ii + 1):ncol(MOld)) {
      M[ii, jj] <- M[jj, ii] <- K
      K <- K + 1
    }
  }

  ### Implement Metropolis-Hastings Alg ###
  n.chg <- 0
  models <- cbind(rep(1, dim(mod.mat[[1]])[1]), mod.mat[[1]])
  mweights <-
    matrix(0, nrow(models), length(responses))
  colnames(mweights) <- name
  mem.Samp <-
    list(MOld)
  mweights <-
    mweights + models.Count(Samp = mem.Samp[[1]], models = models)
  MAP.list <- list(mem.Samp[[1]])
  MAP.count <- c(1)
  mem.Samp[[2]] <-
    update.MH(MOld, M, responses, size, shape1, shape2, mod.mat, Prior)
  mweights <-
    mweights + models.Count(Samp = mem.Samp[[2]], models = models)
  Samp.Sum <- mem.Samp[[1]] + mem.Samp[[2]]
  if (sum(mem.Samp[[2]] == mem.Samp[[1]]) < length(mem.Samp[[2]])) {
    n.chg <- n.chg + 1
  }
  i.Map <-
    which(unlist(lapply(
      MAP.list,
      FUN = function(old, new) {
        return(sum(new == old) == length(new))
      },
      new = mem.Samp[[2]]
    )))
  if (length(i.Map) > 0) {
    MAP.count[i.Map] <- MAP.count[i.Map] + 1
  } else {
    MAP.list[[length(MAP.list) + 1]] <-
      mem.Samp[[2]]
    MAP.count <- c(MAP.count, 1)
  }
  for (KK in 3:niter.MCMC) {
    # print(KK)
    mem.Samp[[KK]] <- update.MH(mem.Samp[[KK - 1]], M, responses, size, 
      shape1, shape2, mod.mat, Prior)

    # print(mem.Samp[[KK]])
    mweights <- mweights + models.Count(Samp = mem.Samp[[KK]], models = models)
    Samp.Sum <- Samp.Sum + mem.Samp[[KK]]
    if (sum(mem.Samp[[KK]] == mem.Samp[[KK - 1]]) < 
        length(mem.Samp[[KK - 1]])) {
      n.chg <- n.chg + 1
    }
    i.Map <-
      which(unlist(lapply(
        MAP.list,
        FUN = function(old, new) {
          return(sum(new == old) == length(new))
        },
        new = mem.Samp[[KK]]
      )))
    if (length(i.Map) > 0) {
      MAP.count[i.Map] <- MAP.count[i.Map] + 1
    } else {
      MAP.list[[length(MAP.list) + 1]] <-
        mem.Samp[[KK]]
      MAP.count <- c(MAP.count, 1)
    }
  }

  ### Compute Posterior Model Weights ###
  pweights <-
    list()
  for (KK in 1:ncol(mweights)) {
    pweights[[KK]] <- mweights[, KK] / niter.MCMC
  }

  ### List for post-processing ###
  MODEL <-
    list(
      responses = responses,
      size = size,
      name = name,
      shape1 = shape1,
      shape2 = shape2,
      models = models,
      pweights = pweights,
      p0 = p0,
      alpha = HPD.alpha,
      alternative = alternative
    )

  ### Compute and output results ###
  PEP <-
    Samp.Sum / niter.MCMC
  rownames(PEP) <- colnames(PEP) <- MODEL$name
  MAP <-
    MAP.list[[order(MAP.count, decreasing = TRUE)[1]]]
  rownames(MAP) <- colnames(MAP) <- MODEL$name

  if (is.null(call)) {
    call <- match.call()
  }

  ret <-
    list(
      responses = responses,
      size = size,
      name = name,
      p0 = p0,
      alpha = HPD.alpha,
      alternative = alternative,
      shape1 = shape1,
      shape2 = shape2,
      Prior = Prior,
      call = call
    )
  ret$mod.mat <- mod.mat
  ret$Inital <-
    M.init
  rownames(ret$Inital) <- colnames(ret$Inital) <- MODEL$name
  ret$models <- models
  ret$pweights <- pweights
  ret$samples <- sample_posterior.full_bayes(MODEL)
  ret$accept.rate <- (n.chg) / niter.MCMC
  ret$mean_est <- colMeans(ret$samples)
  ret$median_est <- apply(ret$samples, 2, median)
  ret$PEP <- PEP
  ret$MAP <- MAP
  ret$HPD <-
    apply(
      ret$samples,
      MARGIN = 2,
      FUN = boa.hpd,
      alpha = MODEL$alpha
    )
  ret$post.prob <- mem.PostProb(MODEL, fit = ret)
  ret$ESS <-
    calc.ESS.from.HPD(fit = ret, alpha = MODEL$alpha)
  names(ret$ESS) <- MODEL$name
  # ret$ESS2 <-
  #  calc.ESS.from.HPDwid(fit = ret, alpha = MODEL$alpha)
  # names(ret$ESS2) <- MODEL$name
  # ret$PostProb <-
  #  mem.PostProb(MODEL, method = "samples", fit = ret)
  class(ret) <- c("full_bayes", "exchangeability_model")

  clusterRet <- clusterComp(ret)
  result <- list(call = call, basketwise = ret, clusterwise = clusterRet)
  return(result)
}
