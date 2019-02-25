
#' @export
summary.full_bayes <- function(object, ...) {
  # object <- MHResult1$clusterwise
  if (class(object$post.prob) == "matrix") {
    cdfS <- c()
    for (i in 1:dim(object$post.prob)[1]) {
      ss <- paste("Post.Prob ", rownames(object$post.prob)[i])
      cdfS <- c(cdfS, ss)
    }
  } else {
    cdfS <- "Post.Prob"
  }
  oct <- rbind(object$post.prob, object$HPD, object$ESS, object$mean_est, object$median_est)
  rownames(oct) <- c(cdfS, "HPD LB", "HPD HB", "ESS", "Mean", "Median")
  oct <- round(oct, 3)
  return(oct)
}
