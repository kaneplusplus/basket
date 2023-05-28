
#' @title Get the Current Random Seed
#'
#' @description Get the value of the current random seed. If one has not been
#' initialized, then initialize it and return the new seed value.
#' @export
get_seed <- function() {
  if (!exists(".Random.seed")) {
    set.seed(as.integer(Sys.time()))
  }
  .Random.seed
}

#' @title Create a Basket Trial Analysis
#'
#' @description This function creates an analysis modeling the exchangeability
#' and distribution characteristics of cohorts in a basket trial, where
#' a single therapy is given to multiple cohorts. The model is based on the
#' multi-source exchangeability model. This is a generalization of the
#' Bayesian hierarchical model and it accomodates multiple sets of latent
#' factors shared combinations of cohorts.
#' @param responses the number of responses in each basket.
#' @param size the size of each basket.
#' @param name the name of each basket.
#' @param p0 the null response rate for the poster probability calculation.
#' (default 0.15)
#' @param shape1 the first shape parameter(s) for the prior of each basket.
#' (default 0.5)
#' @param shape2 the second shape parameter(s) for the prior of each basket.
#' (default 0.5)
#' @param prior the matrix giving the prior inclusion probability
#' for each pair of baskets. The default is on on the main diagonal and 0.5
#' elsewhere.
#' @param hpd_alpha the highest posterior density trial significance.
#' (default 0.05)
#' @param alternative the alternative case definition (default "greater").
#' @param call the call of the function. (default NULL)
#' @param cluster_function a function to cluster baskets.
#' @param method "mcmc", "lmem" (local MEM),  or "exact". See details for an 
#' explanation.
#' (default "mcmc")
#' @param mcmc_iter if the method is "mcmc" then this spcifies the number of
#' mcmc iterations. Otherwise, it is ignored. (default 200000)
#' @param mcmc_burnin if the method is "mcmc" then this specifies the number of
#' burn-in iterations. (default 50000)
#' @param initial_mem if the method is "mcmc" then this spcifies the initial
#' MEM matrix. Otherwise, it is ignored.
#' @param cluster_analysis if the cluster analysis is conducted.
#' @param seed the random seed for the mcmc calculations. By default this is
#' .Random.seed. If this value is not initialized, then it is first initialized
#' with Sys.time() and then returned.
#' @details The model may be fit using an exact calculation via
#' mcmc, or a local version of MEM. The first conducts posterior inference 
#' through the entire set of
#' exchangeability relationships in the sample domain. This approach is
#' computationally feasible only when the number of cohorts is relatively
#' small. As a rule of thumb this option should be used with a maximum of
#' 20 cohorts.  By default, the second(mcmc) is used and it is based on
#' the Metropolis algorithm and it extends the model's implementation to
#' larger collections of subpopulations. The algorithm initiates with a
#' burn-in period (see mcmc_burnin), which are discarded from the analysis.
#' The third is based on a recent paper by Wei et al. The calculation
#' is computationally efficient and tends to borrow to a lesser extent
#' than other implemented methods.
#' @examples
#' \donttest{
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
#'   name = letters[1:3]
#' )
#'
#' summary(basket(trials$responses, trials$size, trials$name))
#' }
#' @importFrom crayon red
#' @export
basket <- function(responses,
                   size,
                   name,
                   p0 = 0.15,
                   shape1 = 0.5,
                   shape2 = 0.5,
                   prior = diag(length(responses)) / 2 +
                     matrix(0.5,
                       nrow = length(responses),
                       ncol = length(responses)
                     ),
                   hpd_alpha = 0.05,
                   alternative = "greater",
                   call = NULL,
                   cluster_function = cluster_membership,
                   method = c("mcmc", "lmem", "exact"),
                   mcmc_iter = 200000,
                   mcmc_burnin = 50000,
                   initial_mem = round(prior - 0.001),
                   cluster_analysis = FALSE,
                   seed = get_seed()) {
  if (method[1] == "exact") {
    mem_exact(responses, size, name,
      p0 = p0, shape1 = shape1,
      shape2 = shape2, prior = prior,
      hpd_alpha = hpd_alpha, alternative = alternative, seed = seed,
      cluster_analysis = cluster_analysis,
      call = call, cluster_function = cluster_function
    )
  } else if(method[1] == "mcmc") {
    mem_mcmc(responses, size, name,
      p0 = p0, shape1 = shape1, shape2 = shape2,
      prior = prior, hpd_alpha = hpd_alpha,
      alternative = alternative, mcmc_iter = mcmc_iter,
      mcmc_burnin = mcmc_burnin, initial_mem = initial_mem,
      seed = seed,
      cluster_analysis = cluster_analysis,
      call = call,
      cluster_function = cluster_function
    )
  } else if (method[1] == "lmem") {
    if (is.null(call)) {
      call <- match.call()
    }
    partitions = basket_part(responses, size, shape1, shape2)
    clusters = basket_cluster(partitions, shape1, shape2, responses, size)
    post_summary = summary_basket(clusters, 
                                  q = 1 - hpd_alpha, 
                                  p0 = p0, 
                                  basket_name = name,
                                  basket_member = clusters$member)

    post_prob = partitions[which.max(partitions$pp_grp), 1:(ncol(partitions)-1)]
    names(post_prob) = name
    mean_est = as.vector(clusters$alpha / (clusters$alpha + clusters$beta))
    names(mean_est) = name
    median_est = post_summary[,'Median']
    names(median_est) = name
    hpd = t(post_summary[,c(3, 4)])
    rownames(hpd) = c("Lower Bound", "Upper Bound")
    colnames(hpd) = name
    ess = clusters$alpha + clusters$beta
    names(ess) = name
    list(
      call = call,
      basket = list(
        responses = responses,
        size = size,
        name = name,
        p0 = p0,
        alternative = alternative,
        shape1 = shape1,
        shape2 = shape2,
        prior = prior,
        hpd_alpha = hpd_alpha,
        initial = prior,
        post_prob  = post_prob,
        median_est = median_est,
        pep = matrix(clusters$Smat, 
                     nrow = nrow(clusters$Smat), 
                     ncol = ncol(clusters$Smat), 
                     dimnames = list(name, name)),
        hpd = hpd,
        ess = ess
      ),
      cluster = clusters$member
    )
  } else {
    stop(red("Unsupported method."))
  }
}

#' @title Sample Posterior Samples from a Basket Trial
#'
#' @description Sample Posterior Samples from a Basket Trial
#' @param model the exchangeability model
#' @param num_samples the number of samples to draw. Default 10000
#' @examples
#' \donttest{
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
#' }
#' @export
sample_posterior <- function(model, num_samples = 10000) {
  UseMethod("sample_posterior", model)
}

#' @importFrom crayon red
sample_posterior.default <- function(model, num_samples = 10000) {
  stop(red(
    "Don't know how to sample posterior from an object of type",
    paste(class(model), collapse = ", "), "."
  ))
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
#' \donttest{
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
#' }
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
