


#' @title Update Full Bayes results with different p0 values
#'
#' @description Update Full Bayes results with different p0 values and alternative
#' @param res TODO: WHAT IS THIS?
#' @param p0 the null response rate for the poster probability calculation
#' (default 0.15).
#' @param alternative the alternative case defination (default greater)

#' @examples
#' MHResult1New <- updateResult(MHResult1, 0.25)
#' @importFrom stats median

#' @export
#'

update_result <- function(res, p0 = 0.15, alternative = "greater") {
  if (length(p0) == 1) {
    p0 <- rep(p0, length(res$basketwise$responses))
  }
  ret <- res
  # basket
  ret$basketwise$p0 <- p0
  ret$basketwise$alternative <- alternative


  if (grepl("exact", res$call[1])) {
    xvec <- res$basketwise$responses
    nvec <- res$basketwise$size
    name <- res$basketwise$name
    CDF <- rep(NA, length(xvec))
    models <- res$basketwise$models
    pweights <- res$basketwise$pweights
    shape1 <- res$basketwise$shape1
    shape2 <- res$basketwise$shape2
    names(CDF) <- name
    CDF[1] <- eval.Post(p0[1], xvec, nvec, models, pweights[[1]], shape1[1], shape2[1], alternative)

    K <- length(xvec)
    for (j in 2:(K - 1)) {
      Ii <- c(j, 1:(j - 1), (j + 1):K)
      CDF[j] <- eval.Post(p0[j], xvec[Ii], nvec[Ii], models, pweights[[j]], shape1[j], shape2[j], alternative)
    }
    j <- j + 1
    Ii <- c(j, 1:(j - 1))
    CDF[j] <- eval.Post(p0[j], xvec[Ii], nvec[Ii], models, pweights[[j]], shape1[j], shape2[j], alternative)
    ret$basketwise$post.prob <- CDF
  } else {
    MODEL <- ret$basketwise
    ret$basketwise$post.prob <- mem.PostProb(MODEL, fit = ret$basketwise)
  }
  # cluster
  retB <- ret$basketwise
  sampleC <- ret$clusterwise$samples
  numClusters <- length(ret$clusterwise$name)
  p0Test <- unique(retB$p0)
  allCDF <- matrix(0, 0, 2)
  for (kk in 1:length(p0Test))
  {
    if (retB$alternative == "greater") {
      res1 <- unlist(lapply(
        1:numClusters,
        FUN = function(j, x, t) {
          return(sum(x[[j]] > t) / length(x[[j]]))
        },
        x = sampleC,
        t = p0Test[kk]
      ))
    } else {
      res1 <- unlist(lapply(
        1:numClusters,
        FUN = function(j, x, t) {
          return(sum(x[[j]] > t) / length(x[[j]]))
        },
        x = sampleC,
        t = p0Test[kk]
      ))
    }
    allCDF <- rbind(allCDF, res1)
  }
  colnames(allCDF) <- ret$clusterwise$name
  rownames(allCDF) <- p0Test
  ret$clusterwise$post.prob <- allCDF

  return(ret)
}
