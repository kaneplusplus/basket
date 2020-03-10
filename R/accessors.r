
#' @title Get the Cluster Members of MEM Models
#'
#' @description Object returned by the `mem_mcmc()` and `mem_exact()` include
#' information about the arms in the trials and the cluster
#' composed of sets of similar arms. This function returns the name
#' of each arm in a cluster.
#' @param x either an exchangeability model or basket object.
#' @return A named list is returned where the name is the cluster name and
#' each element of the list is comprise do a character vector of the
#' baskets in each cluster.
#' @examples
#' \donttest{
#' # Create an MEM analysis of the Vemurafenib trial data.
#' data(vemu_wide)
#'
#' mem_analysis <- mem_exact(
#'   vemu_wide$responders,
#'   vemu_wide$evaluable,
#'   vemu_wide$baskets
#' )
#'
#' # Get the baskets in the clusters.
#' cluster_baskets(mem_analysis)
#' }
#' @export
cluster_baskets <- function(x) {
  UseMethod("cluster_baskets", x)
}

#' @export
cluster_baskets.default <- function(x) {
  stop(red(
    "Don't know how to arms from an object of type",
    paste(class(x), collapse = ", "), "."
  ))
}

#' @export
cluster_baskets.mem_cluster <- function(x) {
  x$cluster
}

#' @export
cluster_baskets.exchangeability_model <- function(x) {
  cluster_baskets(x$cluster)
}

#' @title The Basketwise Posterior Exchangeability Probability Matrix
#'
#' @description MEM analyses include the posterior exchangeability
#' probability (PEP) of included arms giving the probability
#' that any two arms are exchangeable. This function returns the matrix
#' of those probabilities.
#' @param x either an exchangeability model or basket object.
#' @examples
#' \donttest{
#' # Create an MEM analysis of the Vemurafenib trial data.
#' data(vemu_wide)
#'
#' mem_analysis <- mem_exact(
#'   vemu_wide$responders,
#'   vemu_wide$evaluable,
#'   vemu_wide$baskets
#' )
#'
#' # Get the PEP for baskets.
#' basket_pep(mem_analysis)
#' }
#' @export
basket_pep <- function(x) {
  UseMethod("basket_pep", x)
}

#' @importFrom crayon red
#' @export
basket_pep.default <- function(x) {
  stop(red(
    "Don't know how to extract posterior probability matrix from",
    "an object of type",
    paste(class(x), collapse = ", "), "."
  ))
}

#' @export
basket_pep.mem_basket <- function(x) {
  x$pep
}

#' @export
basket_pep.exchangeability_model <- function(x) {
  basket_pep(x$basket)
}

#' @title Get the Basketwise Maximum A Posteriori Probability Matrix
#'
#' @description MEM analyses include the maximum a posterior exchangeability
#' probability (MAP) of included arms indicating whether
#' two arms in the trial are exchangeable. This function returns the matrix
#' of those relationships.
#' @param x either an exchangeability model or basket object.
#' @examples
#' \donttest{
#' # Create an MEM analysis of the Vemurafenib trial data.
#' data(vemu_wide)
#'
#' mem_analysis <- mem_exact(
#'   vemu_wide$responders,
#'   vemu_wide$evaluable,
#'   vemu_wide$baskets
#' )
#'
#' # Get basket MAPs.
#' basket_map(mem_analysis)
#' }
#' @export
basket_map <- function(x) {
  UseMethod("basket_map", x)
}

#' @importFrom crayon red
#' @export
basket_map.default <- function(x) {
  stop(red(
    "Don't know how to extract maximum a posteriori probability",
    "matrix from an object of type", paste(class(x), collapse = ", "), "."
  ))
}

#' @export
basket_map.mem_basket <- function(x) {
  x$map
}

#' @export
basket_map.exchangeability_model <- function(x) {
  basket_map(x$basket)
}

#' @title Get the Clusterwise Posterior Exchangeability Matrix
#'
#' @description MEM analyses include the posterior exchangeability
#' probability (PEP) of clusters of arms giving the probability
#' that any two arms are exchangeable. This function returns the matrix
#' of those probabilities.
#' @param x either an exchangeability model or basket object.
#' @examples
#' \donttest{
#' # Create an MEM analysis of the Vemurafenib trial data.
#' data(vemu_wide)
#'
#' mem_analysis <- mem_exact(
#'   vemu_wide$responders,
#'   vemu_wide$evaluable,
#'   vemu_wide$baskets
#' )
#'
#' # Get cluster PEPs.
#' basket_pep(mem_analysis)
#' }
#' @export
cluster_pep <- function(x) {
  UseMethod("cluster_pep", x)
}

#' @importFrom crayon red
#' @export
cluster_pep.default <- function(x) {
  stop(red(
    "Don't know how to extract posterior probability matrix from",
    "an object of type",
    paste(class(x), collapse = ", "), "."
  ))
}

#' @export
cluster_pep.mem_basket <- function(x) {
  x$pep
}

#' @export
cluster_pep.exchangeability_model <- function(x) {
  cluster_pep(x$basket)
}

#' @title Get the Clusterwise Maximum A Posteriori Probability Matrix
#'
#' @description MEM analyses include the maximum a posterior exchangeability
#' probability (MAP) of included arms indicating whether
#' two arms in the trial are exchangeable. This function returns the matrix
#' of those relationships.
#' @param x either an exchangeability model or basket object.
#' @examples
#' \donttest{
#' # Create an MEM analysis of the Vemurafenib trial data.
#' data(vemu_wide)
#'
#' mem_analysis <- mem_exact(
#'   vemu_wide$responders,
#'   vemu_wide$evaluable,
#'   vemu_wide$baskets
#' )
#'
#' # Get the cluster MAPs.
#' cluster_map(mem_analysis)
#' }
#' @export
cluster_map <- function(x) {
  UseMethod("cluster_map", x)
}

#' @importFrom crayon red
#' @export
cluster_map.default <- function(x) {
  stop(red(
    "Don't know how to extract maximum a posteriori probability",
    "matrix from an object of type",
    paste(class(x), collapse = ", "), "."
  ))
}

#' @export
cluster_map.mem_basket <- function(x) {
  x$map
}

#' @export
cluster_map.exchangeability_model <- function(x) {
  basket_map(x$basket)
}
