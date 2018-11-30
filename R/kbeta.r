
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
#' @aliases dkbeta pkbeta qkbeta rkbeta
#' @usage dkbeta(x, shape1, shape2, a = 0, b = 1, log = FALSE)
#' pkbeta(q, shape1, shape2, a = 0, b = 1, lower.tail = TRUE, log.p = FALSE)
#' qkbeta(p, shape1, shape2, a = 0, b = 0, lower.tail = TRUE, log.p = FALSE)
#' rkbeta(n, shape1, shape2, a = 0, b = 0)
#' @rdname kbeta
#' @name kbeta
#' @examples
#' a <- 1
#' b <- 3
#' shape1 <- 0.5
#' shape2 <- 0.5
#' x <- seq(a, b, by = 0.01)
#' d <- tibble(x = x, y = dkbeta(x, 0.5, 0.5, a, b))
#' ggplot(d, aes(x = x, y = y)) + 
#'   geom_area(fill = "black", alpha = 0.7) +
#'   xlim(0, 3)
#' @export
dkbeta <- function(x, shape1, shape2, a = 0, b = 1, log = FALSE) {
  dbeta(to_beta(x, a, b), shape1, shape2, log = log)  
}

#' @export
pkbeta <- function(q, shape1, shape2, a = 0, b = 1, lower.tail = TRUE, 
                   log.p = FALSE) {
  pbeta(to_beta(q, a, b), shape1, shape2, lower.tail = lower.tail, 
        log.p = log.p)
}

#' @export
qkbeta <- function(p, shape1, shape2, a = 0, b = 0, lower.tail = TRUE, 
                   log.p = FALSE) {
  to_kbeta(qbeta(p, shape1, shape2, lower.tail = lower.tail, log.p = log.p),
           a, b)
}

#' @export
rkbeta <- function(n, shape1, shape2, a = 0, b = 0) {
  to_kbeta(rbeta(n, shape1, shape2), a, b)
}
