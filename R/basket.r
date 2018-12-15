
#' @export
sample_posterior <- function(model, num_samples = 10000) {
  UseMethod("sample_posterior", model)
}

#' @export
basket <- function(formula, data, p0, method = "mem_empirical_bayes", ...) {
  stop("Not implemneted.")
}

#' @export
basket_name <- function(model) {
  UseMethod("basket_name", model)
}

#' @export
basket_name.default <- function(model) {
  stop(paste("Don't know how to get basket names for model of type",
             class(model)))
}

#' @export
basket_name.exchangeability_model <- function(model) {
  ret <- NULL
  if ("name" %in% names(model)) {
    ret <- model$name
  }
  ret
}

#' @export 
cluster_mean <- function(x, method = "maximizer") {
  UseMethod("cluster_mean", x)
}

#' @export
cluster_mean.exchangeability_model <- function(x, method = "maximizer") {
}

#' @export
mean.exchangeability_model <- function(x, ...) {
  x$mean_est
}

#' @export
median.exchangeability_model <- function(x, na.rm = FALSE, ...) {
  x$median_est
}

#' @export
quantile.exchangeability_model <- function(x, ...) {
  apply(x$samples, 2, quantile, ...)
}

#' @export
print.exchangeability_model <- function(x, 
  digits = max(3L, getOption("digits") - 3L), ...) {

  if (!is.null(x$call)) {
    cat("\nCall:\n")
    print(x$call)
  }

  cat("\nMedian estimate of basket response rates:\n")
  med <- x$median_est
  names(med) <- x$name
  print(zapsmall(med, digits + 1L), digits = digits)
  cat("\n")
  invisible(x)
}
