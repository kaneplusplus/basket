
#' @title Sample Posterior Samples from a Basket Trial
#'
#' @description Sample Posterior Samples from a Basket Trial 
#' @param model the exchangeability model
#' @param num_samples the number of samples to draw. Default 10000
#' @examples
#' # WRITE THIS
#' @export
sample_posterior <- function(model, num_samples = 10000) {
  UseMethod("sample_posterior", model)
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
#' trials <- data.frame(
#'   responses = rbinom(trial_sizes, trial_sizes, resp_rate),
#'   size = trial_sizes,
#'   name = paste("Basket", seq_len(5))
#' )
#' 
#' mem_empirical_bayes(trials$responses, trials$size, trials$basket)
#' @export
basket_name <- function(model) {
  UseMethod("basket_name", model)
}

#' @export
basket_name.default <- function(model) {
  stop(paste(
    "Don't know how to get basket names for model of type",
    class(model)
  ))
}

#' @export
basket_name.exchangeability_model <- function(model) {
  ret <- NULL
  if ("name" %in% names(model)) {
    ret <- model$name
  }
  ret
}
