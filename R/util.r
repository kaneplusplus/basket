
#' @importFrom foreach getDoParWorkers
num_workers <- function(factor = 1) {
  getDoParWorkers() * factor
}

boa_hpd <- function(x, alpha) {
  n <- length(x)
  m <- max(1, ceiling(alpha * n))
  y <- sort(x)
  a <- y[seq_len(m)]
  b <- y[(n - m + 1):n]
  i <- order(b - a)[1]
  structure(c(a[i], b[i]), names = c("Lower Bound", "Upper Bound"))
}

mem_mat <- function(indices, mod_mat, h) {
  m <- matrix(NA, h, h)
  diag(m) <- rep(1, dim(m)[1])
  if (h == 2) {
    m[1, 2] <- m[2, 1] <- mod_mat[[1]][indices[1]]
    return(m)
  }
  for (i in seq_len(length(indices) - 1)) {
    m[(i + 1):dim(m)[2], i] <- m[i, (i + 1):dim(m)[2]] <-
      mod_mat[[i]][indices[i], ]
  }
  m[dim(m)[2], i + 1] <- m[i + 1, dim(m)[2]] <-
    c(0, 1)[indices[length(indices)]]
  m
}

log_marg_dens <- function(i_indices, mod_mat, xvec, nvec, avec, bvec) {
  m <- mem_mat(i_indices, mod_mat, length(xvec))
  marg_vec <- rep(NA, dim(m)[1])
  # calculate the product portion of integrated marginal likelihood
  prod_vec <- beta(xvec + avec, nvec + bvec - xvec) / beta(avec, bvec)
  for (i in seq_len(dim(m)[1])) {
    p_vec <- prod(prod_vec^(1 - m[i, ]))
    marg_vec[i] <- (beta(avec[i] + m[i, ] %*% xvec, bvec[i] + m[i, ] %*%
      (nvec - xvec)) / beta(avec[i], bvec[i])) *
      p_vec
  }
  sum(log(marg_vec))
}

ess <- function(x, n, omega, a, b) {
  alpha <- a + omega %*% x
  beta <- b + (omega %*% n - omega %*% x)
  alpha + beta
}

mem_prior <- function(i_indices, mod_mat, prior_inclusion) {
  m <- mem_mat(i_indices, mod_mat, nrow(prior_inclusion))
  mem <- m[upper.tri(m)]
  source_vec <- prior_inclusion[upper.tri(prior_inclusion)]
  s_in <- source_vec # source inclusion probability
  s_ex <- 1 - source_vec # source exclusion probability
  prod(s_in^(mem) * s_ex^(1 - mem))
}

mem_ij <- function(i_indices, mod_mat, prior_inclusion, i, j) {
  mem_mat(i_indices, mod_mat, nrow(prior_inclusion))[i, j]
}

ij_pep <- function(j, mod_i, mod_mat, prior_inclusion, i, log_marg, prior) {
  i_1 <- apply(mod_i, MARGIN = 1, FUN = mem_ij, mod_mat, prior_inclusion, i, j)
  ((exp(log_marg) * prior) %*% i_1) / sum(exp(log_marg) * prior)
}

comp_pep <- function(i, j_indices, mod_i, mod_mat, prior_inclusion, 
                    log_marg, prior) {
  vapply((i + 1):j_indices,
    FUN = ij_pep, FUN.VALUE = NA_real_,
    mod_i, mod_mat, prior_inclusion, i, log_marg, prior
  )
}

mem_j <- function(i_indices, mod_mat, prior_inclusion, j, m) {
  mm <- mem_mat(i_indices, mod_mat, nrow(prior_inclusion))
  as.numeric(sum(as.numeric(mm[j, ] == m)) == length(m))
}

j_weight_1 <- function(i, mod_i, mod_mat, prior_inclusion, j, log_marg, prior) {
  if (length(log_marg) == 2) {
    i_m <- apply(mod_i,
      MARGIN = 1, FUN = mem_j, mod_mat, prior_inclusion, j,
      c(1, mod_mat[[1]][i])
    )
  } else {
    i_m <- apply(mod_i,
      MARGIN = 1, FUN = mem_j, mod_mat, prior_inclusion, j,
      c(1, mod_mat[[1]][i, ])
    )
  }
  ((exp(log_marg) * prior) %*% i_m) / sum(exp(log_marg) * prior)
}

j_weight_j <- function(i, mod_i, mod_mat, prior_inclusion, j, log_marg, prior) {
  u <- mod_mat[[1]][i, ]
  i_m <- apply(mod_i,
    MARGIN = 1, FUN = mem_j, mod_mat, prior_inclusion, j,
    c(u[1:(j - 1)], 1, u[j:length(u)])
  )
  ((exp(log_marg) * prior) %*% i_m) / sum(exp(log_marg) * prior)
}

j_weight_j_indices  <- function(i, mod_i, mod_mat, prior_inclusion, j, 
                                log_marg, prior) {
  if (length(log_marg) == 2) {
    i_m <- apply(mod_i,
      MARGIN = 1, FUN = mem_j, mod_mat, prior_inclusion, j,
      c(mod_mat[[1]][i], 1)
    )
  } else {
    i_m <- apply(mod_i,
      MARGIN = 1, FUN = mem_j, mod_mat, prior_inclusion, j,
      c(mod_mat[[1]][i, ], 1)
    )
  }
  ((exp(log_marg) * prior) %*% i_m) / sum(exp(log_marg) * prior)
}

post_weights <- function(j, j_indices, mod_i, mod_mat, prior_inclusion,
                         log_marg, prior) {

  if (length(log_marg) == 2) {
    if (j == 1) {
      out <- vapply(1:2,
        FUN = j_weight_1, FUN.VALUE = NA_real_,
        mod_i, mod_mat, prior_inclusion, j, log_marg, prior
      )
    } else if (j == j_indices) {
      out <- vapply(1:2,
        FUN = j_weight_j_indices, FUN.VALUE = NA_real_,
        mod_i, mod_mat, prior_inclusion, j, log_marg, prior
      )
    }
    return(out)
  }

  if (j == 1) {
    out <- vapply(seq_len(nrow(mod_mat[[1]])),
      FUN = j_weight_1, FUN.VALUE = NA_real_,
      mod_i, mod_mat, prior_inclusion, j, log_marg, prior
    )
  } else if (j == j_indices) {
    out <- vapply(seq_len(nrow(mod_mat[[1]])),
      FUN = j_weight_j_indices, FUN.VALUE = NA_real_,
      mod_i, mod_mat, prior_inclusion, j, log_marg, prior
    )
  } else {
    out <- vapply(seq_len(nrow(mod_mat[[1]])),
      FUN = j_weight_j, FUN.VALUE = NA_real_,
      mod_i, mod_mat, prior_inclusion, j, log_marg, prior
    )
  }
  out
}

euc_dist <- function(x1, x2, w = c(1, 1)) {
  if (sum(is.na(x1)) > 1) {
    Inf
  } else {
    sqrt(sum(w * ((x1 - x2)^2)))
  }
}

#' @importFrom stats qbeta
dist_beta_hpd <- function(ess, fit, alpha, jj) {
  al <- fit$mean_est[jj] * ess
  al <- max(1e-2, al)
  be <- ess - al
  be <- max(1e-2, be)
  euc_dist(fit$hpd[, jj], qbeta(
    c(alpha / 2, 1 - alpha / 2),
    al, be
  ))
}

#' @importFrom GenSA GenSA
ess_from_hpd_i <- function(jj, fit, alpha) {
  opt <-
    GenSA(
      par = 1,
      fn = dist_beta_hpd,
      lower = 0,
      upper = 10000000,
      fit = fit,
      alpha = alpha,
      jj = jj
    )
  opt$par
}

#' @importFrom foreach %dopar%
ess_from_hpd <- function(fit, alpha) {
  ## fit is list with median vec and HPD vec ##
  i <- NULL
  foreach(i = seq_along(fit$mean_est), .combine = c) %dopar% {
    ess_from_hpd_i(i, fit, alpha)
  }
}

####################################################################
########## MCMC for Bayes with Metropolis-Hastings #################
####################################################################

flip_mem <- function(v, m_old, m) {
  out <- m_old
  if (m_old[which(m == v)[1]] == 1) {
    out[which(m == v)] <- 0
  } else {
    out[which(m == v)] <- 1
  }
  out
}

mem_prior_mat <- function(m, mod_mat, prior_inclusion) {
  mem <- m[upper.tri(m)]
  source_vec <- prior_inclusion[upper.tri(prior_inclusion)]
  s_in <- source_vec # source inclusion probability
  s_ex <- 1 - source_vec # source exclusion probability
  prod(s_in^(mem) * s_ex^(1 - mem))
}

i_models <- function(hh, models, samp) {
  k <- length(models[1, ])
  if (hh == 1) {
    ii <- seq_len(k)
  } else if (hh == k) {
    ii <- c(k, seq_len(k - 1))
  } else {
    ii <- c(hh, seq_len(hh - 1), (hh + 1):k)
  }
  which(apply(models,
    MARGIN = 1,
    FUN = function(x, t) {
      sum(x == t)
    },
    t = samp[hh, ii]
  ) == k)
}

models_count <- function(samp, models) {
  out <- matrix(0, nrow(models), ncol(models))
  u <- vapply(seq_len(ncol(models)),
    FUN = i_models,
    FUN.VALUE = NA_integer_,
    models = models,
    samp = samp
  )
  for (i in seq_along(u)) {
    out[u[i], i] <- 1
  }
  out
}

#' @importFrom crayon red
mem_post_prob <- function(model, fit) {
  if (model$alternative == "greater") {
    out <-
      vapply(
        seq_len(ncol(fit$samples)),
        FUN = function(j, x, t) {
          return(sum(x[, j] > t[j]) / length(x[, j]))
        },
        FUN.VALUE = NA_real_,
        x = fit$samples,
        t = model$p0
      )
  } else if (model$alternative == "less") {
    out <-
      vapply(
        seq_len(ncol(fit$samples)),
        FUN = function(j, x, t) {
          return(sum(x[, j] < t[j]) / length(x[, j]))
        },
        FUN.VALUE = NA_real_,
        x = fit$samples,
        t = model$p0
      )
  }
  else {
    stop(red("Alternative must be either \"greater\" or \"less\"."))
  }

  names(out) <- model$name
  return(out)
}


#############################################################
### Added Shape and Direction to Posterior Prob #############
#############################################################

#' @importFrom stats pbeta
eval.Post <- function(p0, x, n, omega, w, a, b, alternative = "greater") {
  alpha <- a + omega %*% x
  beta <- b + (omega %*% n - omega %*% x)
  if (alternative == "greater") {
    out <- sum((1 - pbeta(p0, alpha, beta)) * w)
  } else {
    out <- sum((pbeta(p0, alpha, beta)) * w)
  }
  out
}

#############################################################
### Fixed shape paramter specification in gen_post() ########
#############################################################

#' @importFrom stats rbeta
gen_post <- function(x, n, omega, a, b) {
  alpha <- a + omega %*% x
  beta <- b + (omega %*% n - omega %*% x)
  rbeta(1, alpha, beta)
}

#' @importFrom stats rmultinom
samp_post <- function(x, n, omega, w, a, b) {
  return(gen_post(x, n, omega[which(rmultinom(1, 1, w) == 1), ], a, b))
}


sample_posterior_model <- function(model, num_samples = 100000) {
  ret <- replicate(
    num_samples,
    samp_post(
      model$responses,
      model$size,
      model$models,
      model$pweights[[1]],
      model$shape1[1],
      model$shape2[1]
    )
  )
  k <- length(model$responses)
  if (k > 2) {
    ret <- rbind(
      ret,
      foreach(j = seq_len(k - 1)[-1], .combine = rbind) %do% {
        ii <- c(j, seq_len(j - 1), (j + 1):k)
        replicate(
          num_samples,
          samp_post(
            model$responses[ii],
            model$size[ii],
            model$models,
            model$pweights[[j]],
            model$shape1[j],
            model$shape2[j]
          )
        )
      }
    )
  }
  j <- k
  ii <- c(j, 1:(j - 1))
  ret <- rbind(
    ret,
    replicate(
      num_samples,
      samp_post(
        model$responses[ii],
        model$size[ii],
        model$models,
        model$pweights[[j]],
        model$shape1[j],
        model$shape2[j]
      )
    )
  )
  ret <- t(ret)
  dimnames(ret) <- list(NULL, model$name)
  ret
}

#' @importFrom stats rbeta
samp_one_group <- function(x, n, a, b, num_samples = 100000) {
  return(rbeta(num_samples, x + a, n - x + b))
}

#' @importFrom stats pbeta
eval_post_one_group <- function(p0, x, n, a, b, alternative = "greater") {
  alpha <- a + x
  beta <- b + n - x
  if (alternative == "greater") {
    out <- 1 - pbeta(p0, alpha, beta)
  } else {
    out <- pbeta(p0, alpha, beta)
  }
  out
}

#' Cluster Baskets Based on the Posterior Exchangeabilities
#'
#' This is the default function used to cluster cohorts in the
#' \code{basket}, \code{mem_mcmc}, and \code{mem_exact} functions.
#' The approach creates a graph where each vertex is a cohort and the
#' weight between two cohorts is determined by their posterior exchangeability
#' probability. The graph is then clustered using \pkg{igraph}'s
#' \code{louvain} function, which determines the number of clusters and
#' the cluster memberships, and has been shown to perform well with
#' real clinical data.
#' @param m the adjacency matrix.
#' @return A factor variable with cluster memberships for each cohort in
#' the study.
#' @importFrom igraph graph_from_adjacency_matrix E cluster_louvain
#' @seealso basket mem_mcmc mem_exact
#' @export
cluster_membership <- function(m) {
  graph <-
    graph_from_adjacency_matrix(m,
      mode = "undirected",
      weighted = TRUE,
      diag = FALSE
    )
  factor(cluster_louvain(graph, weights = E(graph)$weight)$membership)
}

cluster_comp <- function(basket, cluster_function) {
  n_vec <- length(basket$size)
  map <- basket$map
  name <- basket$name
  all_samples <- basket$samples
  result <- cluster_function(map)
  num_clusters <- length(levels(result))

  cluster_sample <- list()
  cluster_name <- c()
  cluster_element <- list()
  for (k in 1:num_clusters) {
    rank <- which(result == levels(result)[k])

    cluster_element[[k]] <- name[rank]
    cluster_name <- c(cluster_name, paste0("Cluster ", k))

    if (n_vec == 1) {
      sample_vector <- as.vector(all_samples)
    } else {
      sample_vector <- as.vector(all_samples[, rank])
    }
    cluster_sample[[k]] <- sample_vector
  }
  names(cluster_sample) <- cluster_name

  ret <- basket
  ret$p0 <- unique(ret$p0)
  ret$name <- cluster_name
  ret$cluster <- cluster_element
  ret$samples <- cluster_sample

  p0_test <- unique(ret$p0)
  all_cdf <- matrix(0, 0, num_clusters)
  for (kk in seq_along(p0_test)) {
    if (ret$alternative == "greater") {
      res <- unlist(lapply(
        seq_len(num_clusters),
        FUN = function(j, x, t) {
          return(sum(x[[j]] > t) / length(x[[j]]))
        },
        x = cluster_sample,
        t = p0_test[kk]
      ))
    } else {
      res <- unlist(lapply(
        seq_len(num_clusters),
        FUN = function(j, x, t) {
          return(sum(x[[j]] > t) / length(x[[j]]))
        },
        x = cluster_sample,
        t = p0_test[kk]
      ))
    }
    all_cdf <- rbind(all_cdf, res)
  }
  colnames(all_cdf) <- cluster_name
  rownames(all_cdf) <- p0_test
  ret$post_prob <- all_cdf

  ret$mean_est <- unlist(lapply(ret$samples, mean))
  ret$median_est <- unlist(lapply(ret$samples, median))
  ret$hpd <- matrix(unlist(
    lapply(
      ret$samples,
      FUN = boa_hpd,
      alpha = ret$alpha
    )
  ), nrow = 2, byrow = FALSE)
  colnames(ret$hpd) <- cluster_name
  row.names(ret$hpd) <- c("Lower Bound", "Upper Bound")
  ret$ess <- ess_from_hpd(fit = ret, alpha = ret$alpha)
  names(ret$ess) <- cluster_name
  class(ret) <- c("mem_exact", "exchangeability_model")
  ret
}


update_mh <- function(m_old,
                      m,
                      xvec,
                      nvec,
                      avec,
                      bvec,
                      mod_mat,
                      prior_ep,
                      beta_vec, old_dens, prod_vec) {
  k <- max(m, na.rm = TRUE) + 1
  v <- sample(
    seq_len(k - 1),
    seq_len(k - 1)[which(rmultinom(
      1, 1,
      (rev(seq_len(k - 1))^3) / sum(rev(seq_len(k - 1))^3)
    ) == 1)]
  )

  m_prop <- flip_mem(v[1], m_old, m)
  if (length(v) > 1) {
    for (ii in seq_along(v)[-1]) {
      m_prop <- flip_mem(v[ii], m_prop, m)
    }
  }

  if (is.na(old_dens)) {
    old_dens <-
      log_marg_dens_mcmc(m_old, mod_mat, xvec, nvec, avec, bvec, beta_vec,
                         prod_vec) +
      log(mem_prior_mat(m_old, mod_mat, prior_ep))
  }
  new_dens <-
    log_marg_dens_mcmc(m_prop, mod_mat, xvec, nvec, avec, bvec, beta_vec, 
                       prod_vec) +
    log(mem_prior_mat(m_prop, mod_mat, prior_ep))

  rho <- exp(new_dens - old_dens)

  if (rho >= 1) {
    old_dens <- new_dens
    out <- m_prop
  } else {
    if (rbinom(1, 1, rho) == 1) {
      old_dens <- new_dens
      out <- m_prop
    } else {
      out <- m_old
    }
  }
  list(out, old_dens)
}

log_marg_dens_mcmc <- function(m, mod_mat, xvec, nvec, avec, bvec, beta_vec,
                               prod_vec) {
  marg_vec <- rep(NA, dim(m)[1])

  # calculate the product portion of integrated marginal likelihood
  for (i in seq_len(dim(m)[1])) {
    p_vec <- prod(prod_vec^(1 - m[i, ]))
    marg_vec[i] <-
      (beta(avec[i] + m[i, ] %*% xvec, bvec[i] + m[i, ] %*% (nvec - xvec)) /
        beta_vec[i]) * p_vec
  }
  sum(log(marg_vec))
}
