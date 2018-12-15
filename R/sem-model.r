
#' @title Single-Source Model
#'
#' @description Fit the single source model.
#' @references Berry SM, Broglio KR, Groshen S, Berry DA. Bayesian 
#' hierarchical modeling of patient subpopulations: efficient designs of 
#' phase II oncology clinical trials. Clinical Trials. 2013 Oct;10(5):720-34.
#' @importFrom stats update
#' @importFrom rjags jags.model
#' @export
sem <- function(responses, size, name, p0 = 0.25, mu = 0, 
  M = 0, inv_sig = 1, num_chains = 3, quiet = TRUE, num_iterations = 3000, 
  progress_bar = interactive()) {

  if (length(responses) != length(size)) {
    stop(paste("The length of the responses and size parameters",
               "must be equal."))
  }
  if (length(mu) == 1) {
    mu <- rep(mu, length(size))
  }
  inits <- list(mu = mu, M = M, inv.sig2 = inv_sig)
  data <- list(J = length(responses), N = size, Y = responses)
  jags_model <- jags.model(system.file("inst", "jags-models", 
                                       "single-source-model.jags"),
                           data, inits, n.chains = num_chains, quiet = quiet)
  update(jags_model, num_iterations, progress.bar = progress_bar)

  if (missing(name)) {
    name <- as.character(NULL)
  } else {
    if (!is.character(name) || 
        length(name) != length(size)) {

      stop(paste("The basket name argument must be a character vector with",
                 "one name per basket."))
    }
  }

  ret <- list(jags_model = jags_model, responses = responses, 
              size = size, p0 = p0, mu = mu, M = M, 
              inv_sig = inv_sig, num_chains = num_chains, quiet = quiet, 
              num_iterations = num_iterations)
  class(ret) <- c("sem", "exchangeability_model")
  ret
}

#' @export
summary.sem <- function(object, ...) {
  # TODO: fix this
#S <- coda.samples(m, "pr", n.iter=10000, progress.bar="none"); dput(data,"data.txt")
#  PPROB <- c( sum( c(S[[1]][,1],S[[2]][,1],S[[3]][,1]) > p0 )/length( c(S[[1]][,1],S[[2]][,1],S[[3]][,1]) ),
#              sum( c(S[[1]][,2],S[[2]][,2],S[[3]][,2]) > p0 )/length( c(S[[1]][,2],S[[2]][,2],S[[3]][,2]) ),
#              sum( c(S[[1]][,3],S[[2]][,3],S[[3]][,3]) > p0 )/length( c(S[[1]][,3],S[[2]][,3],S[[3]][,3]) ),
#              sum( c(S[[1]][,4],S[[2]][,4],S[[3]][,4]) > p0 )/length( c(S[[1]][,4],S[[2]][,4],S[[3]][,4]) ),
#              sum( c(S[[1]][,5],S[[2]][,5],S[[3]][,5]) > p0 )/length( c(S[[1]][,5],S[[2]][,5],S[[3]][,5]) ),
#              sum( c(S[[1]][,6],S[[2]][,6],S[[3]][,6]) > p0 )/length( c(S[[1]][,6],S[[2]][,6],S[[3]][,6]) ) )
# names(PPROB) <-  c("NSCLC ","CRC.v ","CRC.vc","  BD  ","ED.LH "," ATC  ")
}
