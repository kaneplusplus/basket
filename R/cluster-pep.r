

#' @title Cluster results based on PEP
#'
#' @description Conductthe the basket clustering based on the PEP information and provide the information on the cluster level. .
#' @param responses the number of responses in each basket.
#' @param size the size of each basket.
#' @param name the name of each basket.
#' @param models the models calculated from the MEM computation.
#' @param pweights the weights calculated from the MEM computation.
#' @param PEP the PEP matrix calculated from the MEM computation.
#' @param p0 the null response rate for the poster probability calculation
#' (default 0.15).
#' @param shape1 the first shape parameter(s) for the prior of each basket
#' (default 0.5).
#' @param shape2 the second shape parameter(s) for the prior of each basket
#' (default 0.5).
#' @param HPD.alpha the highest posterior density trial significance.
#' @param alternative the alternative case definition (default "greater"). 
#' TODO: WHAT DOES THIS MEAN? WHAT ARE THE ALTERNATIVES?
#' @param seed the random number seed.
#' @param num_samples the number of sample iterations for each cluster.
#' @param call the call of the function.
#' @examples
#' # 5 baskets, each with enrollement size 5
#' trial_sizes <- rep(15, 6)
#' # The response rates for the baskets.
#' resp_rate <- 0.15
#' # The trials: a column of the number of responses, a column of the
#' # the size of each trial.
#' trials <- data.frame(
#'   responses = c(1, 4, 5, 0, 1, 6),
#'   size = trial_sizes
#' )
#' MHResult2 <- mem_full_bayes_mcmc(trials$responses, trials$size,
#'   name = c(" D1 ", " D2 ", " D3 ", " D4 ", " D5 ", " D6 ")
#' )
#' 
#' #result <- cluster_PEP(
#' #  responses = trials$responses, size = trials$size,
#' #  name = c(" D1 ", " D2 ", " D3 ", " D4 ", " D5 ", " D6 "),
#' #  models = MHResult2$models, pweights = MHResult2$pweights, p0 = 0.15,
#' #  PEP = MHResult2$PEP
#' #)
#' print(trials)
#' #print(result$clusters)
#' #print(result$HPD)
#' #print(result$mean_est)
#' @importFrom stats median
#' @importFrom foreach foreach %dopar% getDoParName getDoSeqName registerDoSEQ
#' %do%
#' @importFrom igraph graph_from_adjacency_matrix cluster_louvain E
#' @export
cluster_PEP <- function(responses,
                        size,
                        name,
                        models,
                        pweights,
                        PEP,
                        p0 = 0.15,
                        shape1 = 0.5,
                        shape2 = 0.5,
                        HPD.alpha = 0.05,
                        alternative = "greater",
                        seed = 1000,
                        num_samples = 20000,
                        call = NULL) {
  j <- NULL
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
  if (length(p0) == 1) {
    p0 <- rep(p0, length(responses))
  }

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

  # browser()
  if (is.null(call)) {
    call <- match.call()
  }

  graph <- graph_from_adjacency_matrix(PEP, mode = "undirected", 
      weighted = TRUE, diag = FALSE)

  result <- factor(cluster_louvain(graph, weights = E(graph)$weight)$membership)
  # samples <- sample_posterior.full_bayes(MODEL)
  # sample2 <- foreach(j = 1:6, .combine = cbind) %do%
  # {
  #  samplePostOneBasket(j, MODEL)
  # }


  numClusters <- length(levels(result))
  sampleC <- matrix(0, num_samples, numClusters)
  cName <- c()
  clusterElement <- list()
  for (k in 1:numClusters) {
    rank <- which(result == levels(result)[k])

    cBasket <- name[rank]
    clusterElement[[k]] <- cBasket
#    print(paste0("#### Cluster Assignment: ", k))
#    print(rank)
#    print(name[rank])
    cName <- c(cName, paste0("Cluster ", k))
    numSamp <- num_samples / length(rank) + 1
    samples <- foreach(j = rank, .combine = cbind) %do% {
      samplePostOneBasket(j, MODEL, num_samples = numSamp)
    }
    sampV <- as.vector(samples)
    sampleC[, k] <- sampV[sample(1:num_samples)]
  }
  colnames(sampleC) <- cName

  ret <-
    list(
      responses = responses,
      size = size,
      basketname = name,
      name = cName,
      p0 = p0,
      alternative = alternative,
      shape1 = shape1,
      shape2 = shape2,
      clusters = clusterElement,
      # Prior = Prior,
      call = call
    )
  ret$models <- models
  ret$pweights <- pweights
  ret$samples <- sampleC
  # ret$accept.rate <- (n.chg) / niter.MCMC
  ret$mean_est <- colMeans(ret$samples)
  ret$median_est <- apply(ret$samples, 2, median)
  ret$PEP <- PEP
  ret$HPD <-
    apply(
      ret$samples,
      MARGIN = 2,
      FUN = boa.hpd,
      alpha = MODEL$alpha
    )
  ret$ESS <-
    calc.ESS.from.HPD(fit = ret, alpha = MODEL$alpha)
  names(ret$ESS) <- cName
  ret$ESS2 <-
    calc.ESS.from.HPDwid(fit = ret, alpha = MODEL$alpha)
  names(ret$ESS2) <- cName
  # ret$PostProb <- mem.PostProb(MODEL, method = "samples", fit = ret)
  class(ret) <- c("full_bayes", "exchangeability_model")
  return(ret)
}
