#' @title Update Full Bayes results with different p0 values
#'
#' @description After running either `mem_mcmc` or `mem_exact`, the
#' test can be updated without rerunning the entire analysis. This function
#' provides updating of both the null response rate along with the
#' alternative rerunning relevant test.
#' @param res the result of an mem analysis.
#' @param p0 the null response rate for the poster probability calculation
#' (default 0.15).
#' @param alternative the alternative case definition (default greater)
#' @examples
#' \dontrun{
#' # Create an MEM analysis of the Vemurafenib trial data.
#' data(vemu_wide)
#'
#' mem_analysis <- mem_exact(
#'   vemu_wide$responders,
#'   vemu_wide$evaluable,
#'   vemu_wide$baskets
#' )
#'
#' # Update the null from p0 = 0.15 the default, to p = 0.25.
#' update_p0(mem_analysis, 0.20)
#' }
#' @importFrom stats median
#' @importFrom crayon red
#' @export
update_p0 <- function(res, p0 = 0.15, alternative = "greater") {
  if (length(p0) == 1) {
    p0 <- rep(p0, length(res$basket$responses))
  }
  ret <- res
  # basket
  ret$basket$p0 <- p0
  ret$basket$alternative <- alternative

  if (inherits(res, "mem_exact")) {
    xvec <- res$basket$responses
    nvec <- res$basket$size
    name <- res$basket$name
    cdf <- rep(NA, length(xvec))
    models <- res$basket$models
    pweights <- res$basket$pweights
    shape1 <- res$basket$shape1
    shape2 <- res$basket$shape2
    names(cdf) <- name
    # This pattern appears in at least one other place. Should it be a function?
    cdf[1] <- eval.Post(
      p0[1], xvec, nvec, models, pweights[[1]],
      shape1[1], shape2[1], alternative
    )

    k <- length(xvec)
    for (j in 2:(k - 1)) {
      ii <- c(j, 1:(j - 1), (j + 1):k)
      cdf[j] <- eval.Post(
        p0[j], xvec[ii], nvec[ii], models, pweights[[j]],
        shape1[j], shape2[j], alternative
      )
    }
    j <- j + 1
    ii <- c(j, 1:(j - 1))
    cdf[j] <- eval.Post(
      p0[j], xvec[ii], nvec[ii], models, pweights[[j]],
      shape1[j], shape2[j], alternative
    )
    ret$basket$post_prob <- cdf
  } else if (inherits(res, "mem_mcmc")) {
    model <- ret$basket
    ret$basket$post_prob <- mem_post_prob(model, fit = ret$basket)
  } else {
    stop(
      red("Argument `res` should be of either type `mem_mcmc` or `mem_exact`.")
    )
  }

  ret
}
