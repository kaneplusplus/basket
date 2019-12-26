
#' @title Fit the Exact MEM Model
#'
#' @description Fit the MEM model using full Bayesian inference.
#' @param responses the number of responses in each basket.
#' @param size the size of each basket.
#' @param name the name of each basket.
#' @param p0 the null response rate for the poster probability calculation
#' (default 0.15).
#' @param shape1 the first shape parameter(s) for the prior of each basket
#' (default 0.5).
#' @param shape2 the second shape parameter(s) for the prior of each basket
#' (default 0.5).
#' @param prior the matrix giving the prior inclusion probability
#' for each pair of baskets. The default is on on the main diagonal and 0.5
#' elsewhere.
#' @param hpd_alpha the highest posterior density trial significance.
#' @param alternative the alternative case definition (default greater)
#' @param seed the random number seed.
#' @param cluster_analysis if the cluster analysis is conducted.
#' @param call the call of the function (default NULL).
#' @param cluster_function a function to cluster baskets 
#' @seealso cluster_membership
#' @importFrom stats rbinom
#' @examples
#' \donttest{
#' # 3 baskets, each with enrollement size 5
#' trial_sizes <- rep(5, 3)
#' 
#' # The response rates for the baskets.
#' resp_rate <- 0.15
#' 
#' # The trials: a column of the number of responses and a column of the
#' # the size of each trial.
#' trials <- data.frame(
#'   responses = rbinom(trial_sizes, trial_sizes, resp_rate),
#'   size = trial_sizes,
#'   name = letters[1:3]
#' )
#' 
#' summary(mem_exact(trials$responses, trials$size, trials$name))
#' }
#' @importFrom foreach foreach %dopar% getDoParName getDoSeqName registerDoSEQ
#' %do%
#' @importFrom stats median
#' @importFrom crayon red
#' @importFrom itertools isplitRows
#' @export
mem_exact <- function(responses,
                      size,
                      name,
                      p0 = 0.15,
                      shape1 = 0.5,
                      shape2 = 0.5,
                      prior = diag(length(responses)) / 2 +
                        matrix(0.5,
                          nrow = length(responses),
                          ncol = length(responses)
                        ),
                      hpd_alpha = 0.05,
                      alternative = "greater",
                      seed = 1000,
                      cluster_analysis = FALSE,
                      call = NULL,
                      cluster_function = cluster_membership) {
  set.seed(seed)
  h <- mod_i <- NULL
  if (is.null(getDoParName())) {
    registerDoSEQ()
  }

  if (length(responses) != length(size)) {
    stop(red(
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
  
  
  size1 <- size[size != 0]
  alp <- hpd_alpha
  if(length(size1) < 1){
    stop(red(
      "The length of the responses must be equal or greater than 1"
    ))
  }
  if(length(size1) == 1)
  {
    ind <- which(size != 0)
    nVec <- length(size)
    prior_inclusion <- prior
    tM <- matrix(1, nVec, nVec) #diag(nVec)
    colnames(tM) <- rownames(tM) <- name
    if(nVec > 1)
    {
      tM[, ind] <- 0
      tM[ind, ] <- 0
      tM[ind, ind] <- 1
    }
    MAX <- MAP <- PEP <- tM 
    pweights <- rep(0, nVec)
    pweights[ind] <- 1
    HPD <- matrix(NA, 2, nVec)
    pESS <- post.prob <- rep(NA, nVec)
    names(pESS) <- names(post.prob) <- name
    rownames(HPD) <- c("lower", "upper")
    colnames(HPD) <- name    
    a <- shape1
    b <- shape2
    samp <- samp_one_group(responses[1], size[1], a[1], b[1])
    HPD[, 1] <- boa.hpd(samp,alp)
    pESS[1] <- a[1] + b[1] + size[1]
    t<-eval_post_one_group(p0[1], responses[1], size[1], a[1], b[1], 
                           alternative)
    post.prob[1] <- t 


    if (nVec > 1)
    {
      for(i in 2:nVec){
        sampG <- samp_one_group(responses[i], size[i], a[i], b[i])
        samp <-cbind(samp, sampG)
        HPD[, i] <- boa.hpd(sampG,alp)
        pESS[i] <- a[i] + b[i] + size[i]
        post.prob[i] <- eval_post_one_group(p0[i], responses[i], size[i], 
                                            a[i], b[i], alternative)
      }
      colnames(samp)<-name
    }
    if (is.null(call)) {
      call <- match.call()
    }
    ret <-
      list(maximizer = MAX,
           PRIOR = prior_inclusion,
           MAP = MAP,
           PEP = PEP,
           post.prob = post.prob,
           #ESS = pESS,
           HPD = HPD,
           responses = responses,
           size = size,
           name = name,
           p0 = p0,
           alpha = hpd_alpha,
           alternative = alternative,
           pweights = pweights,
           shape1 = shape1,
           shape2 = shape2,
           call = call
      )
    
    ret$samples <- samp
    if(nVec > 1)
    {
      ret$mean_est <- colMeans(ret$samples)
      ret$median_est <- apply(ret$samples, 2, median)
    }else{
      ret$mean_est <- mean(ret$samples)
      ret$median_est <- median(ret$samples)
    }

    
    ret$ESS <- pESS
    
    class(ret) <- c("mem_basket", "mem")
    
    if(cluster_analysis)
    {
      clusterRet <- clusterComp(ret, cluster_function)
      class(clusterRet) <- c("mem_cluster", "mem")
      result <-
        list(
          call = call,
          basket = ret,
          cluster = clusterRet
        )
    }else{
      result <-
        list(
          call = call,
          basket = ret,
          cluster = NA
        )
    }
    class(result) <- c("mem_exact", "exchangeability_model")
    return(result)
  }
    
  #   clusterRet <- clusterComp(ret, cluster_function)
  #   class(clusterRet) <- c("mem_cluster", "mem")
  #   result <-
  #     list(
  #       call = call,
  #       basket = ret,
  #       cluster = clusterRet
  #     )
  #   class(result) <- c("mem_exact", "exchangeability_model")
  #   return(result)
  # }
  
  prior_inclusion <- prior


  xvec <- responses
  nvec <- size
  mod.mat <- list()
  k <- length(xvec) - 1
  j <- 1

  mod.mat <- foreach(k = rev(seq_len(length(xvec) - 1))) %do% {
    m <- as.matrix(expand.grid(rep(list(c(
      0, 1
    )), k)))
    m[order(rowSums(m)), ]
  }

  H <- length(mod.mat)
  temp <- foreach(h = seq_len(H - 1)) %do% {
    seq_len(dim(mod.mat[[h]])[1])
  }
  temp[[H]] <- 1:2

  Mod.I <- as.matrix(expand.grid(temp))

  ## identify maximizer of marginal density ##
  log.Marg <-
    foreach(
      mod_i = isplitRows(Mod.I, chunks = num_workers()),
      .combine = c
    ) %dopar% {
      apply(mod_i,
        MARGIN = 1,
        FUN = logMarg.Dens,
        mod.mat,
        xvec,
        nvec,
        shape1,
        shape2
      )
    }

  max.i <- order(log.Marg, decreasing = TRUE)[1]

  MAX <- MEM.mat(Mod.I[max.i, ], mod.mat, length(xvec))
  colnames(MAX) <- rownames(MAX) <- name

  ## compute prior + posterior MEM probabilities ##
  PRIOR <-
    apply(Mod.I, MARGIN = 1, FUN = mem.Prior, mod.mat, prior_inclusion)
  POST <- (exp(log.Marg) * PRIOR) / sum(exp(log.Marg) * PRIOR)
  map.i <- order(POST, decreasing = TRUE)[1]

  MAP <- MEM.mat(Mod.I[map.i, ], mod.mat, length(xvec))
  colnames(MAP) <- rownames(MAP) <- name

  ## Posterior Exchangeability Prob ##
  PEP. <- matrix(NA, length(xvec), length(xvec))
  colnames(PEP.) <- rownames(PEP.) <- name
  diag(PEP.) <- 1
  i <- NA
  pep <- foreach(
    i = 1:(nrow(PEP.) - 1),
    combine = list,
    .multicombine = TRUE
  ) %dopar% {
    compPEP(i, length(xvec), Mod.I, mod.mat, prior_inclusion, log.Marg, PRIOR)
  }
  l <- lower.tri(PEP.)
  PEP.[l] <- unlist(pep)
  PEP <- t(PEP.)
  colnames(prior_inclusion) <- rownames(prior_inclusion) <- name

  ## Marginal CDF and ESS ##
  # Here's where all of the compute goes.

  pweights <- foreach(
    j = seq_along(xvec),
    .combine = list,
    .multicombine = TRUE
  ) %dopar% {
    post.Weights(j, length(xvec), Mod.I, mod.mat, prior_inclusion,
                 log.Marg, PRIOR)
  }

  pESS <- post.prob <- rep(NA, length(xvec))
  names(pESS) <- names(post.prob) <- name
  HPD <- matrix(NA, 2, length(xvec))
  rownames(HPD) <- c("lower", "upper")
  colnames(HPD) <- name
 
  if (length(xvec) == 2)
  {
    models <- cbind(rep(1, 2), mod.mat[[1]])
  }else{
    models <- cbind(rep(1, dim(mod.mat[[1]])[1]), mod.mat[[1]])    
  }
  

  post.prob[1] <-
    eval.Post(
      p0[1],
      xvec,
      nvec,
      models,
      pweights[[1]],
      shape1[1],
      shape2[1],
      alternative
    )
  # pESS[1] <-
  #   pweights[[1]] %*% ESS(xvec, nvec, models, shape1[1], shape2[1])
  HPD[, 1] <- boa.hpd(
    replicate(
      10000,
      samp.Post(xvec, nvec, models, pweights[[1]], shape1[1], shape2[1])
    ),
    alp
  )

  K <- length(xvec)
  if (K > 2)
  {
    for (j in 2:(K - 1)) {
      Ii <- c(j, 1:(j - 1), (j + 1):K)
      post.prob[j] <-
        eval.Post(p0[j],
                  xvec[Ii],
                  nvec[Ii],
                  models,
                  pweights[[j]],
                  shape1[j],
                  shape2[j],
                  alternative)
      # pESS[j] <-
      # pweights[[j]] %*% ESSi(xvec[Ii], nvec[Ii], models, shape1[j], shape2[j])
      HPD[, j] <- boa.hpd(replicate(
        10000,
        samp.Post(xvec[Ii], nvec[Ii], models, pweights[[j]],
                  shape1[j], shape2[j])
      ),
      alp)
    }
    
    j <- j + 1
    Ii <- c(j, 1:(j - 1))
    post.prob[j] <-
      eval.Post(p0[j],
                xvec[Ii],
                nvec[Ii],
                models,
                pweights[[j]],
                shape1[j],
                shape2[j],
                alternative)
  } else{
    j <- 2
    Ii <- c(2, 1)
    post.prob[j] <-
      eval.Post(p0[j],
                xvec[Ii],
                nvec[Ii],
                models,
                pweights[[j]],
                shape1[j],
                shape2[j],
                alternative)
    
  }
  # pESS[j] <- pweights[[j]] %*%
  #   ESS(xvec[Ii], nvec[Ii], models, shape1[j], shape2[j])
  
  HPD[, j] <- boa.hpd(replicate(
    10000,
    samp.Post(xvec[Ii], nvec[Ii], models, pweights[[j]],
              shape1[j], shape2[j])
  ),
  alp)
  
  if (missing(name)) {
    name <- paste("basket", seq_along(size))
  } else {
    if (is.factor(name)) {
      name <- as.character(name)
    }
    if (!is.character(name) ||
        length(name) != length(size)) {
      stop(red(
        "The basket name argument must be a character vector with",
        "one name per basket."
      ))
    }
  }
  
  if (is.null(call)) {
    call <- match.call()
  }
  ret <-
    list(mod.mat = mod.mat,
      maximizer = MAX,
      PRIOR = prior_inclusion,
      MAP = MAP,
      PEP = PEP,
      post.prob = post.prob,
      #ESS = pESS,
      HPD = HPD,
      responses = responses,
      size = size,
      name = name,
      p0 = p0,
      alpha = hpd_alpha,
      alternative = alternative,
      pweights = pweights,
      shape1 = shape1,
      shape2 = shape2,
      models = models,
      call = call
    )

  ret$samples <- sample_posterior_model(ret)
  ret$mean_est <- colMeans(ret$samples)
  ret$median_est <- apply(ret$samples, 2, median)
  ret$ESS <- calc.ESS.from.HPD(fit = ret, alpha = hpd_alpha)
  class(ret) <- c("mem_basket", "mem")
  
  if(cluster_analysis)
  {
    clusterRet <- clusterComp(ret, cluster_function)
    class(clusterRet) <- c("mem_cluster", "mem")
    result <-
      list(
        call = call,
        basket = ret,
        cluster = clusterRet
      )
  }else{
    result <-
      list(
        call = call,
        basket = ret,
        cluster = NA
      )
  }
  class(result) <- c("mem_exact", "exchangeability_model")
  return(result)
}
