
#' @title Sample Posterior Samples from a Basket Trial
#' 
#' @description TODO FINISH THIS
#' @param model the exchangeability model
#' @param num_samples the number of samples to draw. Default 10000
#' @examples
#' # WRITE THIS
#' @export
sample_posterior <- function(model, num_samples = 10000) {
  UseMethod("sample_posterior", model)
}

#' @title Fit Basket Trial Data
#' 
#' @description TODO FINISH THIS
#' @param formula the number of responders, enrollees, and possibly basket
#' names for a basket trial.
#' @param data the basket trial data
#' @param p0 the null (or failure probability) for the basket trial. 
#' @param method the method for fitting the basket trial. This could be
#' mem_empirical_bayes, mem_full_bayes, sem, or a custom fitting function.
#' @param ... other parameters to be passed to the supplied 
#' \code{method} argument.
#' @examples
#' # WRITE THIS
#' @export
basket <- function(formula, data, p0, method = mem_empirical_bayes, ...) {
  stop("Not implemneted.")
}

#' @title The Names of the Baskets
#'
#' @description Retrieve the basket names in an exchangeability model.
#' @param model the model to retrieve the basket names of
#' @examples
#' #' # 5 baskets, each with enrollement size 5
#' trial_sizes <- rep(5, 5)
#'
#' # The response rates for the baskets.
#' resp_rate <- 0.15
#'
#' # The trials: a column of the number of responses and a column of the
#' # the size of each trial.
#' trials <- data.frame(responses=rbinom(trial_sizes, trial_sizes, resp_rate),
#'                      size = trial_sizes, 
#'                      name = paste("Basket", seq_len(5)))
#'
#' mem_empirical_bayes(trials$responses, trials$size, trials$basket)
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

#' @title Calculate the Exchangeability Cluster Means
#'
#' @description TODO FINISH THIS
#' @param x the exchangeability model
#' @param method the exchangeability model to get the clusters for. This
#' can be either "maximzer" or "pep" if a model includes a posterior 
#' exchangeability probability.
#' @examples
#' # WRITE THIS
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
