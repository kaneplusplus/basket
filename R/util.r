
#' @importFrom foreach getDoParWorkers
num_workers <- function(factor = 1) {
  getDoParWorkers() * factor
}

boa.hpd <- function(x, alpha) {
  n <- length(x)
  m <- max(1, ceiling(alpha * n))
  y <- sort(x)
  a <- y[seq_len(m)]
  b <- y[(n - m + 1):n]
  i <- order(b - a)[1]
  structure(c(a[i], b[i]), names = c("Lower Bound", "Upper Bound"))
}

MEM.mat <- function(Indices, mod.mat, H) {
  M <- matrix(NA, H, H)
  diag(M) <- rep(1, dim(M)[1])
  for (i in seq_len(length(Indices) - 1)) {
    M[(i + 1):dim(M)[2], i] <- M[i, (i + 1):dim(M)[2]] <- 
      mod.mat[[i]][Indices[i], ]
  }
  M[dim(M)[2], i + 1] <- M[i + 1, dim(M)[2]] <- 
    c(0, 1)[ Indices[length(Indices)] ]
  M
}

logMarg.Dens <- function(I, mod.mat, xvec, nvec, avec, bvec) {
  M <- MEM.mat(I, mod.mat, length(xvec))
  marg.vec <- rep(NA, dim(M)[1])
  # calculate the product portion of integrated marginal likelihood
  prod.vec <- beta(xvec + avec, nvec + bvec - xvec) / beta(avec, bvec) 
  for (i in 1:dim(M)[1]) {
    p.vec <- prod(prod.vec^(1 - M[i, ]))
    marg.vec[i] <- (beta(avec[i] + M[i, ] %*% xvec, bvec[i] + M[i, ] %*% 
      (nvec - xvec)) / beta(avec[i], bvec[i])) * 
      p.vec
  } 
  sum(log(marg.vec))
}

ESS <- function(X, N, Omega, a, b) {
  # a <- b <- 0.5
  alph <- a + Omega %*% X
  beta <- b + (Omega %*% N - Omega %*% X)
  alph + beta
}

mem.Prior <- function(I, mod.mat, pr.Inclus) {
  M <- MEM.mat(I, mod.mat, nrow(pr.Inclus))
  mem <- M[ upper.tri(M) ]
  source.vec <- pr.Inclus[upper.tri(pr.Inclus)]
  s.in <- source.vec # source inclusion probability
  s.ex <- 1 - source.vec # source exclusion probability
  prod(s.in^(mem) * s.ex^(1 - mem))
}

mem.ij <- function(I, mod.mat, pr.Inclus, i, j) {
  MEM.mat(I, mod.mat, nrow(pr.Inclus))[i, j]
}

ij.pep <- function(j, Mod.I, mod.mat, pr.Inclus, i, log.Marg, PRIOR) {
  I.1 <- apply(Mod.I, MARGIN = 1, FUN = mem.ij, mod.mat, pr.Inclus, i, j)
  ((exp(log.Marg) * PRIOR) %*% I.1) / sum(exp(log.Marg) * PRIOR)
}

compPEP <- function(i, J, Mod.I, mod.mat, pr.Inclus, log.Marg, PRIOR) {
  vapply((i + 1):J, FUN = ij.pep, FUN.VALUE = NA_real_,
         Mod.I, mod.mat, pr.Inclus, i, log.Marg, PRIOR)
}

mem.j <- function(I, mod.mat, pr.Inclus, j, m) {
  M <- MEM.mat(I, mod.mat, nrow(pr.Inclus))
  as.numeric(sum(as.numeric(M[j, ] == m)) == length(m))
}

j.weight.1 <- function(i, Mod.I, mod.mat, pr.Inclus, j, log.Marg, PRIOR) {
  I.m <- apply(Mod.I,
    MARGIN = 1, FUN = mem.j, mod.mat, pr.Inclus, j,
    c(1, mod.mat[[1]][i, ])
  )
  ((exp(log.Marg) * PRIOR) %*% I.m) / sum(exp(log.Marg) * PRIOR)
}

j.weight.j <- function(i, Mod.I, mod.mat, pr.Inclus, j, log.Marg, PRIOR) {
  u <- mod.mat[[1]][i, ]
  I.m <- apply(Mod.I,
    MARGIN = 1, FUN = mem.j, mod.mat, pr.Inclus, j,
    c(u[1:(j - 1)], 1, u[j:length(u)])
  )
  ((exp(log.Marg) * PRIOR) %*% I.m) / sum(exp(log.Marg) * PRIOR)
}

j.weight.J <- function(i, Mod.I, mod.mat, pr.Inclus, j, log.Marg, PRIOR) {
  I.m <- apply(Mod.I,
    MARGIN = 1, FUN = mem.j, mod.mat, pr.Inclus, j,
    c(mod.mat[[1]][i, ], 1)
  )
  ((exp(log.Marg) * PRIOR) %*% I.m) / sum(exp(log.Marg) * PRIOR)
}

post.Weights <- function(j, J, Mod.I, mod.mat, pr.Inclus, log.Marg, PRIOR) {
  if (j == 1) {
    out <- vapply(seq_len(nrow(mod.mat[[1]])),
      FUN = j.weight.1, FUN.VALUE = NA_real_,
      Mod.I, mod.mat, pr.Inclus, j, log.Marg, PRIOR
    )
  } else if (j == J) {
    out <- vapply(seq_len(nrow(mod.mat[[1]])),
      FUN = j.weight.J, FUN.VALUE = NA_real_,
      Mod.I, mod.mat, pr.Inclus, j, log.Marg, PRIOR
    )
  } else {
    out <- vapply(seq_len(nrow(mod.mat[[1]])),
      FUN = j.weight.j, FUN.VALUE = NA_real_,
      Mod.I, mod.mat, pr.Inclus, j, log.Marg, PRIOR
    )
  }
  out
}

euc.dist <- function(x1, x2, w = c(1, 1)) {
  if (sum(is.na(x1)) > 1) {
    Inf
  } else {
    sqrt(sum(w * ((x1 - x2)^2)))
  }
}

#' @importFrom stats qbeta
dist.beta.HPD <- function(ess, fit, alpha, jj) {
  al <- fit$mean_est[jj] * ess
  al <- max(1e-2, al)
  be <- ess - al
  be <- max(1e-2, be)
  euc.dist(fit$HPD[, jj], qbeta(
    c(alpha / 2, 1 - alpha / 2),
    al, be
  ))
}

#' @importFrom GenSA GenSA
ESS.from.HPD.i <- function(jj, fit, alpha) {
  # library(GenSA)
  opt <-
    GenSA(
      par = 1,
      fn = dist.beta.HPD,
      lower = 0,
      upper = 10000000,
      # control=list(maxit=pars$DTW.maxit),
      fit = fit,
      alpha = alpha,
      jj = jj
    )
  opt$par
}

#' @importFrom foreach %dopar%
calc.ESS.from.HPD <- function(fit, alpha) {
  ## fit is list with median vec and HPD vec ##
  i <- NULL
  foreach(i = seq_along(fit$mean_est), .combine = c) %dopar% {
    ESS.from.HPD.i(i, fit, alpha)
  }
}

####################################################################
########## MCMC for Bayes with Metropolis-Hastings #################
####################################################################

flip.MEM <- function(v, MOld, M) {
  out <- MOld
  if (MOld[ which(M == v)[1] ] == 1) {
    out[ which(M == v) ] <- 0
  } else {
    out[ which(M == v) ] <- 1
  }
  out
}

mem.Prior.Mat <- function(M, mod.mat, pr.Inclus) {
  # M <- MEM.mat(I, mod.mat, nrow(pr.Inclus));
  mem <- M[upper.tri(M)]
  source.vec <- pr.Inclus[upper.tri(pr.Inclus)]
  s.in <- source.vec # source inclusion probability
  s.ex <- 1 - source.vec # source exclusion probability
  prod(s.in^(mem) * s.ex^(1 - mem))
}

I.models <- function(hh, models, Samp) {
  K <- length(models[1, ])
  if (hh == 1) {
    Ii <- 1:K
  } else if (hh == K) {
    Ii <- c(K, 1:(K - 1))
  } else {
    Ii <- c(hh, 1:(hh - 1), (hh + 1):K)
  }
  which(apply(models, MARGIN = 1,
              FUN = function(x, t) {
                sum(x == t)
              },
              t = Samp[hh, Ii]) == K)
}

models.Count <- function(Samp, models) {
  out <- matrix(0, nrow(models), ncol(models))
  u <- vapply(seq_len(ncol(models)),
              FUN = I.models,
              FUN.VALUE = NA_integer_,
              models = models,
              Samp = Samp)
  for (i in seq_along(u)) {
    out[u[i], i] <- 1
  }
  return(out)
}

#' @importFrom crayon red
mem.PostProb <- function(model, fit) {
  
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
eval.Post <- function(p0, X, N, Omega, w, a, b, alternative = "greater") {
  alph <- a + Omega %*% X
  beta <- b + (Omega %*% N - Omega %*% X)
  if (alternative == "greater") {
    out <- sum((1 - pbeta(p0, alph, beta)) * w)
  } else {
    out <- sum((pbeta(p0, alph, beta)) * w)
  }
  return(out)
}

#############################################################
### Fixed shape paramter specification in gen.Post() ########
#############################################################

#' @importFrom stats rbeta
gen.Post <- function(X, N, Omega, a, b) {
  alph <- a + Omega %*% X
  beta <- b + (Omega %*% N - Omega %*% X)
  return(rbeta(1, alph, beta))
}

#' @importFrom stats rmultinom
samp.Post <- function(X, N, Omega, w, a, b) {
  return(gen.Post(X, N, Omega[which(rmultinom(1, 1, w) == 1), ], a, b))
}

sample_posterior_model <- function(model, num_samples = 100000){   
  ret <- replicate(
    num_samples,
    samp.Post(
      model$responses,
      model$size,
      model$models,
      model$pweights[[1]],
      model$shape1[1],
      model$shape2[1]
    )
  )
  K <- length(model$responses)
  if (K > 2) {
    ret <- rbind(ret,
                 foreach(j = seq_len(K - 1)[-1], .combine = rbind) %do% {
                   Ii <- c(j, seq_len(j-1), (j + 1):K)
                   replicate(
                     num_samples,
                     samp.Post(
                       model$responses[Ii],
                       model$size[Ii],
                       model$models,
                       model$pweights[[j]],
                       model$shape1[j],
                       model$shape2[j]
                     )
                   )
                 })
  }
  j <- K
  Ii <- c(j, 1:(j - 1))
  ret <- rbind(
    ret,
    replicate(
      num_samples,
      samp.Post(
        model$responses[Ii],
        model$size[Ii],
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

#' Cluster Baskets Base on their PEP's
#'
#' This is the default function used to cluster cohorts in the 
#' \code{basket}, \code{mem_mcmc}, and \code{mem_exact} functions. 
#' The approach creates a graph where each vertex is a cohort and the
#' weight between two cohorts is determined by their posterior exchangeability
#' probability. The graph is then clustered using \pkg{igraph}'s 
#' \code{spinglass} function, which determines the number of clusters and
#' the cluster memberships, and has been shown to perform well with 
#' real clinical data.
#' @param PEP the posterior probability matrix.
#' @return A factor variable with cluster memberships for each cohort in 
#' the study.
#' @importFrom igraph graph_from_adjacency_matrix E cluster_spinglass
#' @seealso basket mem_mcmc mem_exact
#' @export
cluster_pep_membership <- function(PEP) {
  graph <-
    graph_from_adjacency_matrix(PEP,
      mode = "undirected",
      weighted = TRUE,
      diag = FALSE
    )
  factor(cluster_spinglass(graph, weights = E(graph)$weight)$membership)
}

clusterComp <- function(basketRet, cluster_function) {
  PEP <- basketRet$PEP
  name <- basketRet$name
  p0 <- unique(basketRet$p0)
  allSamp <- basketRet$samples
  result <- cluster_function(PEP)
  numClusters <- length(levels(result))

  sampleC <- list()
  cName <- c()
  clusterElement <- list()
  for (k in 1:numClusters) {
    rank <- which(result == levels(result)[k])

    cBasket <- name[rank]
    clusterElement[[k]] <- cBasket
    cName <- c(cName, paste0("Cluster ", k))

    sampV <- as.vector(allSamp[, rank])
    sampleC[[k]] <- sampV
  }
  names(sampleC) <- cName

  ret <- basketRet
  ret$p0<-unique(ret$p0)
  ret$name <- cName
  ret$cluster <- clusterElement
  ret$samples <- sampleC

  p0Test <- unique(ret$p0)
  allCDF <- matrix(0, 0, numClusters)
  for (kk in seq_along(p0Test)) {
    if (ret$alternative == "greater") {
      res <- unlist(lapply(
        1:numClusters,
        FUN = function(j, x, t) {
          return(sum(x[[j]] > t) / length(x[[j]]))
        },
        x = sampleC,
        t = p0Test[kk]
      ))
    } else {
      res <- unlist(lapply(
        1:numClusters,
        FUN = function(j, x, t) {
          return(sum(x[[j]] > t) / length(x[[j]]))
        },
        x = sampleC,
        t = p0Test[kk]
      ))
    }
    allCDF <- rbind(allCDF, res)
  }
  colnames(allCDF) <- cName
  rownames(allCDF) <- p0Test
  ret$post.prob <- allCDF

  # ret$accept.rate <- (n.chg) / niter.MCMC
  ret$mean_est <- unlist(lapply(ret$samples, mean))
  ret$median_est <- unlist(lapply(ret$samples, median))
  # ret$PEP <- PEP
  ret$HPD <- matrix(unlist(
    lapply(
      ret$samples,
      # MARGIN = 2,
      FUN = boa.hpd,
      alpha = ret$alpha
    )
  ), nrow = 2, byrow = FALSE)
  colnames(ret$HPD) <- cName
  row.names(ret$HPD) <- c("Lower Bound", "Upper Bound")
  ret$ESS <-
    calc.ESS.from.HPD(fit = ret, alpha = ret$alpha)
  names(ret$ESS) <- cName
  class(ret) <- c("mem_exact", "exchangeability_model")
  ret
}


update.MH <- function(MOld,
                      M,
                      xvec,
                      nvec,
                      avec,
                      bvec,
                      mod.mat,
                      Prior.EP,
                      betaV, oldDens, prod.vec) {
  K <- max(M, na.rm = TRUE) + 1
  v <- sample(seq_len(K-1), 
    seq_len(K-1)[which(rmultinom(1, 1, 
      (rev(seq_len(K-1)) ^ 3) / sum(rev(seq_len(K-1)) ^ 3)) == 1)])
  
  MProp <- flip.MEM(v[1], MOld, M)
  if (length(v) > 1) {
    for (ii in seq_along(v)[-1]) {
      MProp <- flip.MEM(v[ii], MProp, M)
    }
  }
  
  if (is.na(oldDens)) {
    oldDens <- 
      logMarg.DensMCMC(MOld, mod.mat, xvec, nvec, avec, bvec, betaV, prod.vec)+ 
      log(mem.Prior.Mat(MOld, mod.mat, Prior.EP))
  }
  newDens <-
    logMarg.DensMCMC(MProp, mod.mat, xvec, nvec, avec, bvec, betaV, prod.vec) + 
    log(mem.Prior.Mat(MProp, mod.mat, Prior.EP))
  
  rho <- exp(newDens - oldDens)
  
  # cat(rho, rho1, xvec,nvec, avec, bvec)
  if (rho >= 1) {
    oldDens <- newDens
    out <- MProp
  } else {
    if (rbinom(1, 1, rho) == 1) {
      oldDens <- newDens
      out <- MProp
    } else {
      out <- MOld
    }
  }
  list(out, oldDens)
}

logMarg.DensMCMC <- function(M, mod.mat, xvec, nvec, avec, bvec, betaV, 
                             prod.vec) {

  marg.vec <- rep(NA, dim(M)[1])
  
  # calculate the product portion of integrated marginal likelihood
  for (i in seq_len(dim(M)[1])) {
    p.vec <- prod(prod.vec^(1 - M[i, ]))
    marg.vec[i] <-
      (beta(avec[i] + M[i, ] %*% xvec, bvec[i] + M[i, ] %*% (nvec - xvec)) /
         betaV[i]) * p.vec
  }
  sum(log(marg.vec))
}

