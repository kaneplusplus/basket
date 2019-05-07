
#' @title Sample Posterior Samples from a Basket Trial
#'
#' @description Sample Posterior Samples from a Basket Trial
#' @param model the exchangeability model
#' @param num_samples the number of samples to draw. Default 10000
#' @examples
#' # 3 baskets, each with enrollement size 5
#' trial_sizes <- rep(5, 3)
#' 
#' # The response rates for the baskets.
#' resp_rate <- 0.15
#' 
#' # The trials: a column of the number of responses and a column of the
#' # the size of each trial.
#' trials <- data.frame(
#'   responses = rbinom(trial_sizes, trial_sizes, resp_rate),
#'   size = trial_sizes,
#'   name = paste("Basket", seq_len(3))
#' )
#' @export
sample_posterior <- function(model, num_samples = 10000) {
  UseMethod("sample_posterior", model)
}

#' @importFrom crayon red
sample_posterior.default <- function(model, num_samples = 10000) {
  stop(red(
    "Don't know how to sample posterior from an object of type",
    paste(class(model), collapse = ", "), "."))
}

sample_posterior.mem <- function(model, num_samples = 10000) {
  sample_posterior_model(model)
}

#' @title The Names of the Baskets
#'
#' @description Retrieve the basket names in an exchangeability model.
#' @param model the model to retrieve the basket names of
#' @importFrom stats rbinom
#' @examples
#' # 3 baskets, each with enrollement size 5
#' trial_sizes <- rep(5, 3)
#' 
#' # The response rates for the baskets.
#' resp_rate <- 0.15
#' 
#' # The trials: a column of the number of responses and a column of the
#' # the size of each trial.
#' trials <- data.frame(
#'   responses = rbinom(trial_sizes, trial_sizes, resp_rate),
#'   size = trial_sizes,
#'   name = paste("Basket", seq_len(3))
#' )
#' 
#' basket_name(mem_mcmc(trials$responses, trials$size, trials$basket))
#' @export
basket_name <- function(model) {
  UseMethod("basket_name", model)
}

#' @importFrom crayon red
#' @export
basket_name.default <- function(model) {
  stop(red(
    "Don't know how to get basket names for model of type",
    paste(class(model), collapse = ", "), "."
  ))
}

#' @export
basket_name.exchangeability_model <- function(model) {
  ret <- NULL
  if ("name" %in% names(model$basket)) {
    ret <- model$basket$name
  }
  ret
}
