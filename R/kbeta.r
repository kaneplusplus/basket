
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
#' @param n number of observations. If ‘length(n) > 1’, the length is
#'        taken to be the number required.
#' @param shape1 the first non-negative parameters of the Beta distribution.
#' @param shape2 the second non-negative parameters of the Beta distribution.
#' @param a the lower bound of the support of the distribution (default -1).
#' @param b the upper bound of the support of the distribution (default 1).
#' @param log should the distribution value be shown on the log scale?
#' (default FALSE)
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x],
#' otherwise, P[X > x].
#' @aliases dkbeta pkbeta qkbeta rkbeta kbeta
#' @usage dkbeta(x, shape1, shape2, a = -1, b = 1, log = FALSE)
#' pkbeta(q, shape1, shape2, a = -1, b = 1, lower.tail = TRUE, log.p = FALSE)
#' qkbeta(p, shape1, shape2, a = -1, b = 1, lower.tail = TRUE, log.p = FALSE)
#' rkbeta(n, shape1, shape2, a = -1, b = 1)
#' @rdname kbeta
#' @name kbeta
#' @importFrom stats dbeta pbeta qbeta rbeta
#' @examples
#' library(ggplot2)
#' a <- 1
#' b <- 3
#' shape1 <- 0.5
#' shape2 <- 0.5
#' x <- seq(a, b, by = 0.01)
#' d <- data.frame(x = x, y = dkbeta(x, 0.5, 0.5, a, b))
#' ggplot(d, aes(x = x, y = y)) +
#'   geom_area(fill = "black", alpha = 0.7) +
#'   xlim(0, 3)
#' @export
dkbeta <- function(x, shape1, shape2, a = -1, b = 1, log = FALSE) {
  dbeta(to_beta(x, a, b), shape1, shape2, log = log)
}

#' @export
pkbeta <- function(q, shape1, shape2, a = -1, b = 1, lower.tail = TRUE,
                   log.p = FALSE) {
  pbeta(to_beta(q, a, b), shape1, shape2,
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @export
qkbeta <- function(p, shape1, shape2, a = -1, b = 1, lower.tail = TRUE,
                   log.p = FALSE) {
  to_kbeta(
    qbeta(p, shape1, shape2, lower.tail = lower.tail, log.p = log.p),
    a, b
  )
}

#' @export
rkbeta <- function(n, shape1, shape2, a = -1, b = 1) {
  to_kbeta(rbeta(n, shape1, shape2), a, b)
}
