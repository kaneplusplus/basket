
#' @title Fit the MEM Model using MCMC
#'
#' @description Fit the MEM model using Bayesian Metropolis-Hasting MCMC
#' inference.
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
#' @param mcmc_iter the number of MCMC iterations.
#' @param mcmc_burnin the number of MCMC Burn_in iterations.
#' @param initial_mem the initial MEM matrix.
#' @param seed the random number seed.
#' @param cluster_analysis if the cluster analysis is conducted.
#' @param call the call of the function.
#' @param cluster_function a function to cluster baskets
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
#' res <- mem_mcmc(trials$responses, trials$size)
#' }
#' @importFrom stats median
#' @importFrom crayon red
#' @importFrom itertools isplitVector
#' @export
mem_mcmc <- function(responses,
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
                     mcmc_iter = 200000,
                     mcmc_burnin = 50000,
                     initial_mem = round(prior - 0.001),
                     seed = 1000,
                     cluster_analysis = FALSE,
                     call = NULL,
                     cluster_function = cluster_membership) {
  set.seed(seed)
  k <- NULL
  if (is.null(getDoParName())) {
    registerDoSEQ()
  }
  if (length(responses) != length(size)) {
    stop(red(
      "The length of the responses and size parameters",
      "must be equal."
    ))
  }

  # If the shape and p0 argument is a single value, make it a vector of
  # appropriate length.
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
    n_vec <- length(size)
    prior_inclusion <- prior
    pep <- matrix(1, n_vec, n_vec)
    colnames(pep) <- rownames(pep) <- name
    if (n_vec > 1) {
      pep[, ind] <- 0
      pep[ind, ] <- 0
      pep[ind, ind] <- 1
    }
    maximizer <- map <- pep
    pweights <- rep(0, n_vec)
    pweights[ind] <- 1
    hpd <- matrix(NA, 2, n_vec)
    p_ess <- post_prob <- rep(NA, n_vec)
    names(p_ess) <- names(post_prob) <- name
    rownames(hpd) <- c("lower", "upper")
    colnames(hpd) <- name
    a <- shape1
    b <- shape2
    samp <- samp_one_group(responses[1], size[1], a[1], b[1])
    hpd[, 1] <- boa_hpd(samp, alp)
    p_ess[1] <- a[1] + b[1] + size[1]
    t <-
      eval_post_one_group(p0[1], responses[1], size[1], a[1], b[1], alternative)
    post_prob[1] <- t


    if (n_vec > 1) {
      for (i in 2:n_vec) {
        group_sample <- samp_one_group(responses[i], size[i], a[i], b[i])
        samp <- cbind(samp, group_sample)
        hpd[, i] <- boa_hpd(group_sample, alp)
        p_ess[i] <- a[i] + b[i] + size[i]
        post_prob[i] <- eval_post_one_group(p0[i], responses[i], size[i],
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
        # ESS = p_ess,
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
    if (n_vec > 1) {
      ret$mean_est <- colMeans(ret$samples)
      ret$median_est <- apply(ret$samples, 2, median)
    } else {
      ret$mean_est <- mean(ret$samples)
      ret$median_est <- median(ret$samples)
    }


    ret$ess <- p_ess

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

  if (!isTRUE(all.equal(diag(prior), rep(1, ncol(prior))))) {
    stop(red("Elements on the main diagonal of `prior` must be 1."))
  }
  ### Produce sample space of MEM ###
  mod_mat <- foreach(k = rev(seq_len(length(responses) - 1))) %do% {
    mem_sample_space <- as.matrix(expand.grid(rep(list(c(
      0, 1
    )), k)))
    mem_sample_space[order(rowSums(mem_sample_space)), ]
  }

  # Error if the inital_mem isn't symmetric.
  if (!isTRUE(all.equal(initial_mem, t(initial_mem)))) {
    stop(red("The `initial_mem` matrix must be symmetric."))
  }
  if (!isTRUE(all(diag(initial_mem) == 1))) {
    stop(red("The main diagonal of the `initial_mem` matrix must be 1's."))
  }

  m_init <- initial_mem
  m_old <- m_init

  ### Create Map for Proposal Distribution ###
  m <- diag(NA, nrow(m_old))
  k <- 1
  for (ii in seq_len(nrow(m_old) - 1)) {
    for (jj in (ii + 1):ncol(m_old)) {
      m[ii, jj] <- m[jj, ii] <- k
      k <- k + 1
    }
  }
  ### Implement Metropolis-Hastings Alg ###
  n_chg <- 0
  mod_mat[[1]] <- as.matrix(mod_mat[[1]])
  models <- cbind(rep(1, dim(mod_mat[[1]])[1]), mod_mat[[1]])
  mweights <- matrix(0, nrow(models), length(responses))
  if (missing(name)) {
    name <- paste("basket", seq_along(size))
  }
  if (is.factor(name)) {
    name <- as.character(name)
  }
  colnames(mweights) <- name
  mem_samp <- list(m_old)
  mweights <- mweights + models_count(samp = mem_samp[[1]], models = models)
  map_list <- list(mem_samp[[1]])
  map_count <- c(1)


  map_hash <- new.env()
  map_hash[[toString(m_old)]] <- 1

  old_dens <- NA
  xvec <- responses
  nvec <- size
  beta_vec <- beta(shape1, shape2)
  prod.vec <- beta(xvec + shape1, nvec + shape2 - xvec) / beta(shape1, shape2)


  for (i in 1:mcmc_burnin) {
    t <- update_mh(
      m_old, m, responses, size,
      shape1, shape2, mod_mat, prior, beta_vec, old_dens, prod.vec
    )
    m_old <- t[[1]]
    old_dens <- t[[2]]
  }

  t <- update_mh(
    m_old, m, responses, size,
    shape1, shape2, mod_mat, prior, beta_vec, old_dens, prod.vec
  )
  mem_samp[[2]] <- t[[1]]
  old_dens <- t[[2]]
  mweights <- mweights + models_count(samp = mem_samp[[2]], models = models)
  samp_sum <- mem_samp[[1]] + mem_samp[[2]]

  if (sum(mem_samp[[2]] == mem_samp[[1]]) < length(mem_samp[[2]])) {
    n_chg <- n_chg + 1
  }

  new <- mem_samp[[2]]
  key <- toString(new)
  if (!is.null(map_hash[[key]])) {
    index <- map_hash[[key]]
    map_count[index] <- map_count[index] + 1
  } else {
    map_list[[length(map_list) + 1]] <-
      mem_samp[[2]]
    map_count <- c(map_count, 1)
    map_hash[[key]] <- length(map_list)
  }


  for (kk in seq_len(mcmc_iter)[-(1:2)]) {
    t <- update_mh(
      mem_samp[[kk - 1]], m, responses, size,
      shape1, shape2, mod_mat, prior, beta_vec, old_dens, prod.vec
    )
    mem_samp[[kk]] <- t[[1]]
    old_dens <- t[[2]]
  }

  it <- NULL
  models_count <- foreach(
    it = isplitVector(seq_len(mcmc_iter)[-(1:2)],
      chunks = num_workers()
    ),
    .combine = c
  ) %dopar% {
    foreach(k = it) %do% {
      models_count(samp = mem_samp[[k]], models = models)
    }
  }

  for (kk in seq_len(mcmc_iter)[-(1:2)]) {
    mweights <- mweights + models_count[[kk - 2]]
    samp_sum <- samp_sum + mem_samp[[kk]]
    if (sum(mem_samp[[kk]] == mem_samp[[kk - 1]]) <
      length(mem_samp[[kk - 1]])) {
      n_chg <- n_chg + 1
    }

    new <- mem_samp[[kk]]
    key <- toString(new)
    if (!is.null(map_hash[[key]])) {
      index <- map_hash[[key]]
      map_count[index] <- map_count[index] + 1
    } else {
      map_list[[length(map_list) + 1]] <-
        mem_samp[[kk]]
      map_count <- c(map_count, 1)
      map_hash[[key]] <- length(map_list)
    }
  }

  # Compute posterior model weights
  pweights <- list()
  for (kk in seq_len(ncol(mweights))) {
    pweights[[kk]] <- mweights[, kk] / mcmc_iter
  }

  # List for post-processing
  model <-
    list(
      responses = responses,
      size = size,
      name = name,
      shape1 = shape1,
      shape2 = shape2,
      models = models,
      pweights = pweights,
      p0 = p0,
      alpha = hpd_alpha,
      alternative = alternative
    )

  ### Compute and output results ###
  pep <- samp_sum / mcmc_iter
  rownames(pep) <- colnames(pep) <- model$name
  map <- map_list[[order(map_count, decreasing = TRUE)[1]]]
  rownames(map) <- colnames(map) <- model$name

  if (is.null(call)) {
    call <- match.call()
  }

  ret <-
    list(
      responses = responses,
      size = size,
      name = name,
      p0 = p0,
      alpha = hpd_alpha,
      alternative = alternative,
      shape1 = shape1,
      shape2 = shape2,
      prior = prior,
      call = call
    )
  ret$mod_mat <- mod_mat
  ret$initial <- m_init

  rownames(ret$initial) <- colnames(ret$initial) <- model$name
  ret$models <- models
  ret$pweights <- pweights

  # Ret doesn't have class information yet so, we'll call
  # sample_posterior.exchangeability_model directly.
  ret$samples <- sample_posterior_model(model)
  ret$accept_rate <- (n_chg) / mcmc_iter
  ret$mean_est <- colMeans(ret$samples)
  ret$median_est <- apply(ret$samples, 2, median)
  ret$pep <- pep
  ret$map <- map
  ret$hpd <- apply(ret$samples, MARGIN = 2, FUN = boa_hpd, alpha = model$alpha)
  ret$post_prob <- mem_post_prob(model, fit = ret)
  #ret$ess <- ess_from_hpd(fit = ret, alpha = model$alpha)
  ret$ess <- ess_from_qnt(fit = ret)
  names(ret$ess) <- model$name
  class(ret) <- c("mem_basket", "mem")

  if (cluster_analysis) {
    cluster_ret <- cluster_comp(ret, cluster_function)
    class(cluster_ret) <- c("mem_cluster", "mem")
    result <-
      list(
        call = call,
        basket = ret,
        cluster = cluster_ret,
        seed = seed
      )
  } else {
    result <-
      list(
        call = call,
        basket = ret,
        cluster = NA,
        seed = seed
      )
  }

  class(result) <- c("mem_mcmc", "exchangeability_model")
  result
}
