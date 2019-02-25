
# project into the support of the beta
to_beta <- function(x, a, b) {
  if (any(b <= a)) {
    stop("The b parameter must be bigger than the a paramter.")
  }
  (x - a) / (b - a)
}

# project into the support of the kbeta
to_kbeta <- function(x, a, b) {
  if (any(b <= a)) {
    stop("The b parameter must be bigger than the a paramter.")
  }
  x * (b - a) + a
}

#' @title The K-Beta Disribution
#'
#' @description Density, distribution function, quantile function, and
#' random geration for the K-Beta distribution with shape parameters
#' \code{shape1}, \code{shape2}, \code{a}, and \code{b}. The shape
#' parameters are the same as those of the beta distribution. The
#' \code{a} parameter defines the smallest value in the support of the
#' distribution and the \code{b} parameter defines the largest.
#' @param x vector of values in the support.
#' @param q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If 