
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

  if (length(p0) == 1) {
    p0 <- rep(p0, length(responses))
  }

  size1 <- size[size != 0]
  alp <- hpd_alpha

  if (length(size1) < 1) {
    stop(red(
      "The length of the responses must be equal or greater than 1"
    ))
  }

  if (length(size1) == 1) {
    ind <- which(size != 0)
    vec_len <- length(size)
    prior_inclusion <- prior
    pep <- matrix(1, vec_len, vec_len)
    colnames(pep) <- rownames(pep) <- name
    if (vec_len > 1) {
      pep[, ind] <- 0
      pep[ind, ] <- 0
      pep[ind, ind] <- 1
    }
    maximizer <- map <- pep
    pweights <- rep(0, vec_len)
    pweights[ind] <- 1
    hpd <- matrix(NA, 2, vec_len)
    p_ess <- post_prob <- rep(NA, vec_len)
    names(p_ess) <- names(post_prob) <- name
    rownames(hpd) <- c("lower", "upper")
    colnames(hpd) <- name
    a <- shape1
    b <- shape2
    samp <- samp_one_group(responses[1], size[1], a[1], b[1])
    hpd[, 1] <- boa_hpd(samp, alp)
    p_ess[1] <- a[1] + b[1] + size[1]
    t <- eval_post_one_group(
      p0[1], responses[1], size[1], a[1], b[1],
      alternative
    )
    post_prob[1] <- t


    if (vec_len > 1) {
      for (i in 2:vec_len) {
        group_sample <- samp_one_group(responses[i], size[i], a[i], b[i])
        samp <- cbind(samp, group_sample)
        hpd[, i] <- boa_hpd(group_sample, alp)
        p_ess[i] <- a[i] + b[i] + size[i]
        post_prob[i] <- eval_post_one_group(
          p0[i], responses[i], size[i],
          a[i], b[i], alternative
        )
      }
      colnames(samp) <- name
    }
    if (is.null(call)) {
      call <- match.call()
    }
    ret <-
      list(
        maximizer = maximizer,
        prior = prior_inclusion,
        map = map,
        pep = pep,
        post_prob = post_prob,
        hpd = hpd,
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
    if (vec_len > 1) {
      ret$mean_est <- colMeans(ret$samples)
      ret$median_est <- apply(ret$samples, 2, median)
    } else {
      ret$mean_est <- mean(ret$samples)
      ret$median_est <- median(ret$samples)
    }


    ret$ess <- p_ess
    names(ret$ess) <- ret$name

    class(ret) <- c("mem_basket", "mem")

    if (cluster_analysis) {
      cluster_ret <- cluster_comp(ret, cluster_function)
      class(cluster_ret) <- c("mem_cluster", "mem")
      result <-
        list(
          call = call,
          basket = ret,
          cluster = cluster_ret
        )
    } else {
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

  prior_inclusion <- prior

  xvec <- responses
  nvec <- size
  mod_mat <- list()
  k <- length(xvec) - 1
  j <- 1

  mod_mat <- foreach(k = rev(seq_len(length(xvec) - 1))) %do% {
    m <- as.matrix(expand.grid(rep(list(c(
      0, 1
    )), k)))
    m[order(rowSums(m)), ]
  }

  temp <- foreach(h = seq_len(length(mod_mat) - 1)) %do% {
    seq_len(dim(mod_mat[[h]])[1])
  }
  temp[[length(mod_mat)]] <- 1:2

  mod_i_mat <- as.matrix(expand.grid(temp))

  ## identify maximizer of marginal density ##
  log_marg <-
    foreach(
      mod_i = isplitRows(mod_i_mat, chunks = num_workers()),
      .combine = c
    ) %dopar% {
      apply(mod_i,
        MARGIN = 1,
        FUN = log_marg_dens,
        mod_mat,
        xvec,
        nvec,
        shape1,
        shape2
      )
    }

  max_i <- order(log_marg, decreasing = TRUE)[1]

  maximizer <- mem_mat(mod_i_mat[max_i, ], mod_mat, length(xvec))
  colnames(maximizer) <- rownames(maximizer) <- name

  ## compute prior + posterior MEM probabilities ##
  prior <- apply(mod_i_mat, MARGIN = 1, FUN = mem_prior, mod_mat, 
                 prior_inclusion)
  post <- (exp(log_marg) * prior) / sum(exp(log_marg) * prior)
  map_i <- order(post, decreasing = TRUE)[1]

  map <- mem_mat(mod_i_mat[map_i, ], mod_mat, length(xvec))
  colnames(map) <- rownames(map) <- name

  ## Posterior Exchangeability Prob ##
  tpep <- matrix(NA, length(xvec), length(xvec))
  colnames(tpep) <- rownames(tpep) <- name
  diag(tpep) <- 1
  i <- NA
  pep <- foreach(
    i = 1:(nrow(tpep) - 1),
    combine = list,
    .multicombine = TRUE
  ) %dopar% {
    comp_pep(i, length(xvec), mod_i_mat, mod_mat, prior_inclusion, log_marg,
             prior)
  }
  l <- lower.tri(tpep)
  tpep[l] <- unlist(pep)
  pep <- t(tpep)
  colnames(prior_inclusion) <- rownames(prior_inclusion) <- name

  ## Marginal CDF and ESS ##
  # Here's where all of the compute goes.

  pweights <- foreach(
    j = seq_along(xvec),
    .combine = list,
    .multicombine = TRUE
  ) %dopar% {
    post_weights(
      j, length(xvec), mod_i_mat, mod_mat, prior_inclusion,
      log_marg, prior
    )
  }

  p_ess <- post_prob <- rep(NA, length(xvec))
  names(p_ess) <- names(post_prob) <- name
  hpd <- matrix(NA, 2, length(xvec))
  rownames(hpd) <- c("lower", "upper")
  colnames(hpd) <- name

  if (length(xvec) == 2) {
    models <- cbind(rep(1, 2), mod_mat[[1]])
  } else {
    models <- cbind(rep(1, dim(mod_mat[[1]])[1]), mod_mat[[1]])
  }

  post_prob[1] <-
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

  hpd[, 1] <- boa_hpd(
    replicate(
      10000,
      samp_post(xvec, nvec, models, pweights[[1]], shape1[1], shape2[1])
    ),
    alp
  )

  k <- length(xvec)

  if (k > 2) {

    for (j in 2:(k - 1)) {

      ii <- c(j, 1:(j - 1), (j + 1):k)

      post_prob[j] <-
        eval.Post(
          p0[j],
          xvec[ii],
          nvec[ii],
          models,
          pweights[[j]],
          shape1[j],
          shape2[j],
          alternative
        )

      hpd[, j] <- boa_hpd(
        replicate(
          10000,
          samp_post(
            xvec[ii], nvec[ii], models, pweights[[j]],
            shape1[j], shape2[j]
          )
        ),
        alp
      )
    }

    j <- j + 1
    ii <- c(j, 1:(j - 1))
    post_prob[j] <-
      eval.Post(
        p0[j],
        xvec[ii],
        nvec[ii],
        models,
        pweights[[j]],
        shape1[j],
        shape2[j],
        alternative
      )
  } else {

    j <- 2
    ii <- c(2, 1)

    post_prob[j] <-
      eval.Post(
        p0[j],
        xvec[ii],
        nvec[ii],
        models,
        pweights[[j]],
        shape1[j],
        shape2[j],
        alternative
      )
  }

  hpd[, j] <- boa_hpd(
    replicate(
      10000,
      samp_post(
        xvec[ii], nvec[ii], models, pweights[[j]],
        shape1[j], shape2[j]
      )
    ),
    alp
  )

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
    list(
      mod_mat = mod_mat,
      maximizer = maximizer,
      prior = prior_inclusion,
      map = map,
      pep = pep,
      post_prob = post_prob,
      hpd = hpd,
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
  ret$ess <- ess_from_hpd(fit = ret, alpha = hpd_alpha)
  names(ret$ess) <- colnames(ret$hpd)
  class(ret) <- c("mem_basket", "mem")

  if (cluster_analysis) {
    cluster_ret <- cluster_comp(ret, cluster_function)
    class(cluster_ret) <- c("mem_cluster", "mem")
    result <-
      list(
        call = call,
        basket = ret,
        cluster = cluster_ret
      )
  } else {
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
