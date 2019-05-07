
#' @title Basket Trial Analysis
#' @name basket-package
#' @aliases basket-package basket
#' @docType package
#' @references \url{https://github.com/kaneplusplus/basket}
#' @description The R basket package provides for the designs and analysis of 
#' basket trials for multi-source exchangeability models (MEM) 
#' (\url{https://academic.oup.com/biostatistics/article/19/2/169/3930935}) 
#' allowing arms to "share" power with similar arms in a trial. The package 
#' is intended to perform the exact or MCMC computation of the operating 
#' characteristics of different scenarios. Calculations derived from these 
#' analyses include the posterior probabilities, HPD boundaries, effective 
#' sample sizes (ESS), mean and median estimations can be calculated with 
#' this package using the MEM method. Along with providing "basketwise"
#' analyses, the package includes similar calculations for "clusterwise" 
#' analyses where a cluster a set of similar baskets. In addition plotting 
#' tools are provided to visualize basket and cluster density as well as their 
#' exchangeability. The package includes the following main functions:
#' \itemize{
#'   \item{[basket_name()] }{Get the basket names in an analysis}
#'   \item{[basket_pep()] }{Get the Posterior Exchangeability Probability 
#'                          (PEP) matrix for the arms in the trial.}
#'   \item{[basket_map()] }{Get the Maximum A Posteriori (MAP) matrix for the 
#'                          arms in a trial.}
#'   \item{[cluster_pep()] }{Get the Posterior Exchangeability Probability 
#'                           (PEP) matrix for the arms in the trial.}
#'   \item{[cluster_map()] }{Get the Maximum A Posteriori (MAP) matrix for 
#'                           the arms in a trial.}
#'   \item{[sample_posterior()] }{Sample from the posterior distribution of 
#'                                the arms in the trial.}
#'   \item{[mem_exact()] }{Create a basket analysis using the exact method.}
#'   \item{[mem_mcmc()] }{Create a basket analysis using the exact method.}
#'   \item{[summary()] }{Summarize the results of an analysis.}
#'   \item{[update_p0()] }{Update the null that a basket response rate is 
#'                         above a specified value.}
#'   \item{[plot_density()] }{Plot the estimated densities of arms or clusters.}
#'   \item{[plot_pep()] }{Show the exchangeogram of the PEP matrix.}
#'   \item{[plot_mem()] }{Plot the arm prior, MAP, and PEP of a basket trial.}
#' }
NULL
