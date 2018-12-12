
#' @export
sample_posterior <- function(model, num_samples = 10000) {
  UseMethod("sample_posterior", model)
}

#' @export
basket <- function(formula, data, p0, method = "mem_empirical_bayes", ...) {
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
  if ("basket_name" %in% names(model)) {
    ret <- model$basket_names
  }
  ret
}

#' @export
mean.exchangeability_model <- function(model) {
  ret$mean_est
}

#' @export
median.exchangeability_model <- function(x, na.rm = FALSE, ...) {
  ret$median_est
}

#' @export
quantile.exchangeability_model <- function(x, ...) {
  dots <- list(...)
  num_samples <- ifelse(is.null(dots$num_samples), 10000, dots$num_samples)
  apply(sample_posterior, 2, quantile, ...)
}

